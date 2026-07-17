# ruff: noqa: E402
import json
import sqlite3
import shlex
import subprocess
import sys
import tempfile
import unittest
from argparse import Namespace
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / 'src'
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from branchmanager.cli import _resolve_novelty_baseline_datasets, build_parser, cmd_filing_cabinet
from branchmanager.db.interface import Database
from branchmanager.it_desk import run_it_desk_checks, write_it_desk_report
from branchmanager.onboarding import validate_submission, write_onboarding_outputs
from branchmanager.marker_provenance import load_marker_provenance, marker_qc_flag
from branchmanager.personnel import load_exit_requests
from branchmanager.pipeline.classify import _parse_vsearch_match
from branchmanager.pipeline.tree import _orient_tree_input_fasta
from branchmanager.pipeline.workflow_helpers import build_selection_decision
from branchmanager.project_state import import_genome_results
from branchmanager.reporting import write_annual_report
from branchmanager.run_manifest import file_record
from branchmanager.utils.subprocess import run_cmd


class OperationalWorkflowTests(unittest.TestCase):
    def test_onboarding_writes_normalised_per_read_map(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            (root / 'well_A.ab1').write_bytes(b'ABIF-placeholder')
            mapping = root / 'mapping.tsv'
            mapping.write_text(
                'sequence_id\t27F\n'
                'ISO1\twell_A.ab1\n'
            )
            metadata = root / 'project_metadata.tsv'
            metadata.write_text(
                'sequence_id\tpartner_id\talready_sequenced\n'
                'ISO1\tQUB\tno\n'
            )
            result = validate_submission(
                mapping,
                partner_metadata=metadata,
                primers=['27F'],
            )
            self.assertEqual(result['status'], 'PASS')
            outputs = write_onboarding_outputs(root / 'out', result)
            self.assertEqual(Path(outputs['report']).name, 'onboarding_report.tsv')
            self.assertEqual(Path(outputs['summary']).name, 'onboarding_summary.json')
            read_map = Path(outputs['read_map']).read_text()
            self.assertIn('ISO1\t' + str((root / 'well_A.ab1').resolve()) + '\t27F\tforward', read_map)
            submission_header = Path(outputs['normalised']).read_text().splitlines()[0].split('\t')
            self.assertNotIn('source_fasta', submission_header)
            self.assertNotIn('sample_name', submission_header)
            self.assertEqual(
                Path(outputs['report']).read_text().splitlines(),
                ['severity\tline\tsequence_id\tcode\tdetail'],
            )

    def test_onboarding_groups_row_per_read_map_and_resolves_trace_prefixes(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            reads = root / 'reads'
            reads.mkdir()
            forward = reads / 'KKX994_70439994.ab1'
            reverse = reads / 'KKY011_70440011.ab1'
            forward.write_bytes(b'ABIF-forward')
            reverse.write_bytes(b'ABIF-reverse')
            mapping = root / 'mapping.csv'
            mapping.write_text(
                'Sequencing ID,Isolate Number,Read\n'
                'KKX994,SW_0016,Forward\n'
                'KKY011,SW_0016,Reverse\n'
            )
            metadata = root / 'partner.tsv'
            metadata.write_text(
                'sequence_id\tpartner_id\talready_sequenced\n'
                'SW_0016\tQUB\tno\n'
            )

            result = validate_submission(
                mapping,
                partner_metadata=metadata,
                read_dir=reads,
                expected_partner_id='QUB',
                dataset='QUB_01',
            )

            self.assertEqual(result['status'], 'PASS')
            self.assertEqual(result['isolates'], 1)
            self.assertEqual(result['read_files'], 2)
            self.assertEqual(
                {row['direction'] for row in result['normalised_reads']},
                {'forward', 'reverse'},
            )
            self.assertEqual(
                {row['processing_mode'] for row in result['normalised_reads']},
                {'assemble'},
            )
            self.assertEqual(result['normalised'][0]['dataset'], 'QUB_01')

    def test_onboarding_accepts_partner_fasta_with_cumulative_metadata(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            fasta = root / 'partner.fasta.gz'
            import gzip
            with gzip.open(fasta, 'wt') as handle:
                handle.write('>ISO_FASTA_1\nACGTRYN\n')
            metadata = root / 'project_metadata.tsv'
            metadata.write_text(
                'sequence_id\tpartner_id\tselected_for_genome_sequencing\talready_sequenced\n'
                'ISO_FASTA_1\tUoG\tyes\tno\n'
                'OLDER_ISOLATE\tQUB\tno\tyes\n'
            )

            result = validate_submission(fasta=fasta, partner_metadata=metadata)
            outputs = write_onboarding_outputs(root / 'out', result)

            self.assertEqual(result['status'], 'PASS')
            self.assertEqual(result['input_type'], 'fasta')
            self.assertEqual(result['isolates'], 1)
            self.assertEqual(result['read_files'], 0)
            self.assertIn('>ISO_FASTA_1\nACGTRYN', Path(outputs['fasta']).read_text())
            self.assertEqual(
                Path(outputs['report']).read_text().splitlines(),
                ['severity\tline\tsequence_id\tcode\tdetail'],
            )
            submission_header = Path(outputs['normalised']).read_text().splitlines()[0].split('\t')
            self.assertIn('source_fasta', submission_header)
            self.assertNotIn('read_files', submission_header)
            self.assertNotIn('sample_name', submission_header)

    def test_onboarding_rejects_duplicate_ledger_ids_and_partner_mismatch(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            read = root / 'read.ab1'
            read.write_bytes(b'ABIF-placeholder')
            mapping = root / 'mapping.tsv'
            mapping.write_text('sequence_id\tfile\nISO1\tread.ab1\n')
            metadata = root / 'metadata.tsv'
            metadata.write_text(
                'sequence_id\tpartner_id\talready_sequenced\n'
                'ISO1\tUoG\tno\n'
                'ISO1\tUoG\tno\n'
            )

            result = validate_submission(
                mapping,
                partner_metadata=metadata,
                expected_partner_id='QUB',
                dataset='QUB_01',
            )

            codes = {row['code'] for row in result['problems']}
            self.assertEqual(result['status'], 'FAIL')
            self.assertIn('DUPLICATE_METADATA_ID', codes)
            self.assertIn('PARTNER_ID_MISMATCH', codes)
            outputs = write_onboarding_outputs(root / 'out', result)
            report = Path(outputs['report']).read_text()
            self.assertIn('DUPLICATE_METADATA_ID', report)
            self.assertIn('PARTNER_ID_MISMATCH', report)

    def test_direct_fasta_is_review_required_without_provenance(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta = Path(tmpdir) / 'markers.fasta'
            fasta.write_text('>ISO1\nACGT\n')
            rows, source = load_marker_provenance(fasta)
            self.assertIsNone(source)
            self.assertEqual(rows[0]['marker_qc_class'], 'QUALITY_UNVERIFIED')
            self.assertEqual(marker_qc_flag(rows[0]), 'MARKER_QC_REVIEW_REQUIRED')
            decision = build_selection_decision({
                'marker_qc_flag': 'MARKER_QC_REVIEW_REQUIRED',
                'placement_flags': 'MARKER_QC_REVIEW_REQUIRED',
                'classification_identity': 99.0,
                'taxonomy': 'd__Bacteria;s__Example species',
            })
            self.assertEqual(decision['decision'], 'REVIEW BEFORE SELECTION')

    def test_genome_metrics_derive_qc_pass_and_update_collection(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            db = Database(str(root / 'project.sqlite'))
            db.initialise()
            db.insert_sequences([('ISO1', 'ACGT')], dataset='Batch1')
            table = root / 'genomes.tsv'
            table.write_text(
                'sequence_id\tgenome_id\tgenome_status\tcompleteness\tcontamination\tgtdb_taxonomy\tani_cluster\n'
                'ISO1\tG1\tSEQUENCED\t95.2\t1.1\td__Bacteria;s__Example species\tANI_001\n'
            )
            results = import_genome_results(db, table)
            self.assertEqual(results[0]['detail'], 'GENOME_QC_PASSED')
            genome = db.get_genome_records()[0]
            self.assertTrue(genome['genome_qc_pass'])
            self.assertEqual(genome['ani_cluster'], 'ANI_001')
            self.assertTrue(db.get_sequencing_metadata_for_ids(['ISO1'])['ISO1']['selected_for_wgs'])
            self.assertEqual(db.get_isolate_statuses(['ISO1'])['ISO1']['status'], 'GENOME_QC_PASSED')

    def test_conflicting_explicit_genome_qc_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            db = Database(str(root / 'project.sqlite'))
            db.initialise()
            db.insert_sequences([('ISO1', 'ACGT')], dataset='Batch1')
            table = root / 'genomes.tsv'
            table.write_text(
                'sequence_id\tgenome_id\tgenome_status\tgenome_qc_pass\tcompleteness\tcontamination\n'
                'ISO1\tG1\tGENOME_QC_PASSED\tyes\t70\t12\n'
            )
            results = import_genome_results(db, table)
            self.assertEqual(results[0]['result'], 'REJECTED')
            self.assertEqual(db.get_genome_records(), [])

    def test_annual_report_contains_project_ledgers_and_run_history(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            db = Database(str(root / 'project.sqlite'))
            db.initialise()
            db.insert_sequences([('ISO1', 'ACGT'), ('ISO2', 'ACGA')], dataset='Batch1')
            db.upsert_dataset_role('Batch1', 'candidate')
            db.update_isolate_status('ISO1', 'MARKER_QC_PASSED')
            db.update_isolate_status('ISO2', 'PROPOSED')
            db.record_project_run('performance-review:1', 'performance_review', 'COMPLETE', dataset='Batch1')
            outputs = write_annual_report(db, root / 'annual_report')
            self.assertEqual(Path(outputs['html']).name, 'annual_report.html')
            html = Path(outputs['html']).read_text()
            self.assertIn('Cumulative Project Overview', html)
            self.assertIn('Proposed sequences', html)
            self.assertIn('<td>Batch1</td><td>candidate</td><td></td><td>0</td><td>1</td><td>2</td>', html)
            self.assertIn('ISO1', Path(outputs['isolates']).read_text())
            self.assertTrue(Path(outputs['candidates']).is_file())
            self.assertTrue(Path(outputs['removals']).is_file())
            self.assertTrue(Path(outputs['decision_changes']).is_file())

    def test_annual_report_discovers_failed_manifest_and_staged_dataset(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            db = Database(str(root / 'project.sqlite'))
            db.initialise()
            db.insert_sequences([('HUN1', 'ACGT')], dataset='Hungate')
            db.upsert_dataset_role('Hungate', 'baseline', genomes_available=True)
            review_dir = root / 'UoG_01' / '03_performance_review_hiring_panel'
            review_dir.mkdir(parents=True)
            staged = Database(str(review_dir / '.performance_review_project.sqlite'))
            staged.initialise()
            staged.insert_sequences([('UOG1', 'ACGT'), ('UOG2', 'ACGA')], dataset='UoG_01')
            staged.upsert_dataset_role('UoG_01', 'candidate')
            staged.update_isolate_status('UOG1', 'MARKER_QC_PASSED')
            staged.update_isolate_status('UOG2', 'TRACE_REVIEW')
            (review_dir / 'run_manifest.json').write_text(json.dumps({
                'workflow': 'performance_review',
                'status': 'FAILED',
                'started_at': '2026-07-17T07:37:43Z',
                'completed_at': '2026-07-17T07:39:32Z',
                'command': ['branchmanager', 'assistant', '--dataset', 'UoG_01'],
                'error': 'local tree context resolved 49/50 assessed sequences',
            }))

            outputs = write_annual_report(db, root / 'annual_report')
            html = Path(outputs['html']).read_text()
            self.assertIn('<td>UoG_01</td><td>candidate</td><td></td><td>0</td><td>0</td><td>2</td>', html)
            self.assertIn('FAILED', html)
            self.assertIn('UOG1', Path(outputs['isolates']).read_text())
            self.assertIn('NOT YET ASSESSED', Path(outputs['candidates']).read_text())

    def test_exit_interview_removes_active_records_and_retains_audit(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            db = Database(str(root / 'project.sqlite'))
            db.initialise()
            db.insert_sequences([('ISO1', 'ACGT'), ('ISO2', 'ACGA')], dataset='Batch1')
            db.upsert_dataset_role('Batch1', 'candidate')
            db.insert_taxonomy([('ISO1', 'd__Bacteria;s__Example one', 0.99, 'Batch1')])
            db.insert_distances([('ISO2', 'Batch1', 'ISO1', 99.0)])
            db.upsert_classification_evidence([{
                'sequence_id': 'ISO1', 'ref_db': 'GTDB', 'best_hit': 'REF1', 'identity': 99.0,
            }])
            db.upsert_sequencing_metadata([{
                'id': 'ISO1', 'partner_id': 'QUB', 'dataset': 'Batch1',
                'selected_for_sequencing': True, 'selected_for_wgs': False,
            }])
            db.save_assessment_snapshot('review:1', [{'id': 'ISO1', 'taxonomy': 'example'}], dataset='Batch1')
            db.save_selection_round('round:1', [{'sequence_id': 'ISO1', 'role': 'PRIMARY'}])
            db.update_isolate_status('ISO1', 'PROPOSED', detail='test')

            requests = load_exit_requests(['ISO1'], default_reason='supplier withdrew isolate')
            planned = db.plan_sequence_removals(requests)
            self.assertEqual(planned[0]['status'], 'READY')
            db.apply_sequence_removals(planned)

            self.assertNotIn('ISO1', db.get_all_ids())
            self.assertIn('ISO2', db.get_all_ids())
            with db.connect() as conn:
                for table, column in (
                    ('taxonomy', 'id'), ('classification_evidence', 'sequence_id'),
                    ('sequencing_metadata', 'id'), ('assessment_snapshots', 'sequence_id'),
                    ('selection_round_members', 'sequence_id'), ('isolate_status', 'sequence_id'),
                    ('seq_aliases', 'canonical_id'),
                ):
                    self.assertEqual(
                        conn.execute(f'SELECT COUNT(*) FROM {table} WHERE {column} = ?', ('ISO1',)).fetchone()[0],
                        0,
                    )
                self.assertEqual(conn.execute('SELECT COUNT(*) FROM distances').fetchone()[0], 0)
                removal = conn.execute(
                    'SELECT sequence_id, original_dataset, reason FROM sequence_removals'
                ).fetchone()
            self.assertEqual(removal, ('ISO1', 'Batch1', 'supplier withdrew isolate'))

    def test_exit_interview_protects_baselines_and_genome_backed_isolates(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            db = Database(str(Path(tmpdir) / 'project.sqlite'))
            db.initialise()
            db.insert_sequences([('BASE1', 'ACGT')], dataset='Hungate')
            db.insert_sequences([('ISO1', 'ACGA')], dataset='Batch1')
            db.upsert_dataset_role('Hungate', 'baseline', genomes_available=True)
            db.upsert_dataset_role('Batch1', 'candidate')
            db.upsert_genome_records([{
                'genome_id': 'G1', 'sequence_id': 'ISO1', 'genome_status': 'SEQUENCED',
            }])
            requests = load_exit_requests(
                ['BASE1', 'ISO1'], default_reason='test protection',
            )
            planned = db.plan_sequence_removals(requests)
            self.assertEqual(
                [row['status'] for row in planned],
                ['BLOCKED_BASELINE', 'BLOCKED_GENOME_RECORD'],
            )
            allowed = db.plan_sequence_removals(
                requests, allow_baseline=True, allow_genome_records=True,
            )
            self.assertEqual([row['status'] for row in allowed], ['READY', 'READY'])

    def test_run_cmd_never_invokes_a_shell(self):
        completed = subprocess.CompletedProcess(['tool'], 0, stdout=b'', stderr=b'')
        with mock.patch('branchmanager.utils.subprocess.subprocess.run', return_value=completed) as called:
            run_cmd('tool --input "path with spaces.fasta"')
        args, kwargs = called.call_args
        self.assertEqual(args[0], ['tool', '--input', 'path with spaces.fasta'])
        self.assertNotIn('shell', kwargs)

    def test_rootstock_failure_does_not_publish_partial_database(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            fasta = root / 'baseline.fasta'
            fasta.write_text('>BASE1\nACGT\n')
            project = root / 'project.sqlite'
            original = Database(str(project))
            original.initialise()
            original.insert_sequences([('EXISTING', 'AAAA')], dataset='Existing')
            args = Namespace(
                fasta=str(fasta), db=str(project), out=str(root / 'rootstock'),
                dataset='Hungate', taxa=None, taxa_assignments=None, ref=None,
                colours=None, anchors=None, build_tree=False,
            )

            def fail_after_staged_write(staged_args):
                staged = Database(staged_args.db)
                staged.initialise()
                staged.insert_sequences([('PARTIAL', 'CCCC')], dataset='Hungate')
                raise RuntimeError('simulated Filing Cabinet failure')

            with mock.patch('branchmanager.cli._cmd_filing_cabinet_impl', side_effect=fail_after_staged_write):
                with self.assertRaisesRegex(RuntimeError, 'simulated Filing Cabinet failure'):
                    cmd_filing_cabinet(args)

            with original.connect() as conn:
                sequence_ids = {row[0] for row in conn.execute('SELECT id FROM sequences')}
            self.assertEqual(sequence_ids, {'EXISTING'})

    def test_tree_orientation_search_failure_is_fatal(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            query = root / 'query.fasta'
            reference = root / 'reference.fasta'
            query.write_text('>ISO1\nACGT\n')
            reference.write_text('>REF1\nACGT\n')
            with mock.patch('branchmanager.pipeline.tree.run_cmd', side_effect=RuntimeError('vsearch failed')):
                with self.assertRaisesRegex(RuntimeError, 'orientation search failed'):
                    _orient_tree_input_fasta(
                        str(query), ref_fasta=str(reference), anchor_fasta=None,
                        outdir=root, label='queries',
                    )

    def test_classification_match_preserves_alignment_evidence(self):
        hit = _parse_vsearch_match([
            'ISO1', 'GTDB1', '98.5', '1200', '12', '6',
            '1', '1200', '50', '1249', '1400', '1500',
        ])
        self.assertEqual(hit['mismatches'], 12)
        self.assertEqual(hit['gaps'], 6)
        self.assertAlmostEqual(hit['query_coverage'], 85.7142857)
        self.assertEqual(hit['target_coverage'], 80.0)

    def test_older_project_database_is_migrated_deterministically(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / 'old.sqlite'
            with sqlite3.connect(path) as conn:
                conn.executescript(
                    'CREATE TABLE sequences (id TEXT PRIMARY KEY, sequence TEXT, length INTEGER);'
                    'CREATE TABLE taxonomy (id TEXT, taxonomy TEXT, confidence REAL);'
                    'CREATE TABLE taxonomy_alt '
                    '(id TEXT, ref_db TEXT, taxonomy TEXT, confidence REAL, best_hit TEXT, identity REAL);'
                    'CREATE TABLE distances (id TEXT, nearest TEXT, identity REAL);'
                    'CREATE TABLE colors (id TEXT PRIMARY KEY, color TEXT, source TEXT);'
                    'CREATE TABLE sequencing_metadata '
                    '(id TEXT PRIMARY KEY, partner_id TEXT, dataset TEXT, selected_for_wgs INTEGER, '
                    'source_id TEXT, source_file TEXT, raw_selected_value TEXT);'
                )
                conn.execute("INSERT INTO sequences VALUES ('ISO1', 'ACGT', 4)")
                conn.execute("INSERT INTO taxonomy VALUES ('ISO1', 'old', 0.8)")
                conn.execute("INSERT INTO taxonomy VALUES ('ISO1', 'new', 0.9)")
                conn.execute("INSERT INTO colors VALUES ('ISO1', '#123456', 'user')")

            db = Database(str(path))
            db.initialise()
            with db.connect() as conn:
                sequence_columns = {row[1] for row in conn.execute('PRAGMA table_info(sequences)')}
                metadata_columns = {row[1] for row in conn.execute('PRAGMA table_info(sequencing_metadata)')}
                taxonomy_rows = conn.execute(
                    'SELECT taxonomy FROM taxonomy WHERE id = ?', ('ISO1',)
                ).fetchall()
                versions = {row[0] for row in conn.execute('SELECT version FROM schema_migrations')}
                colour_row = conn.execute(
                    'SELECT id, colour, source, dataset FROM colours WHERE id = ?', ('ISO1',)
                ).fetchone()
                operational_tables = {
                    row[0] for row in conn.execute(
                        "SELECT name FROM sqlite_master WHERE type='table'"
                    )
                }
            self.assertIn('dataset', sequence_columns)
            self.assertIn('operational_status', metadata_columns)
            self.assertEqual(taxonomy_rows, [('new',)])
            self.assertEqual(colour_row, ('ISO1', '#123456', 'user', ''))
            self.assertNotIn('colors', operational_tables)
            self.assertTrue({1, 2, 3}.issubset(versions))
            self.assertTrue({'project_runs', 'sequence_provenance', 'genome_records'}.issubset(operational_tables))

    def test_manifest_hashes_raw_read_directories(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir) / 'raw reads'
            root.mkdir()
            (root / 'B.ab1').write_bytes(b'B')
            (root / 'A.ab1').write_bytes(b'A')
            first = file_record(root, role='raw_read_directory')
            second = file_record(root, role='raw_read_directory')
            self.assertEqual(first['kind'], 'directory')
            self.assertEqual(first['file_count'], 2)
            self.assertEqual(first['sha256'], second['sha256'])
            (root / 'A.ab1').write_bytes(b'changed')
            self.assertNotEqual(first['sha256'], file_record(root)['sha256'])

    def test_it_desk_rejects_non_fasta_reference(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            reference = root / 'not_fasta.txt'
            reference.write_text('this is not a fasta file\n')
            checks = run_it_desk_checks(references=[reference], output_dir=root)
            reference_check = next(row for row in checks if row['check'] == 'reference:1')
            self.assertEqual(reference_check['status'], 'FAIL')
            outputs = write_it_desk_report(root / 'it_desk', checks)
            self.assertEqual(Path(outputs['tsv']).name, 'it_desk_report.tsv')
            self.assertEqual(Path(outputs['json']).name, 'it_desk_summary.json')

    def test_orientation_command_preserves_paths_with_spaces(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir) / 'path with spaces'
            root.mkdir()
            query = root / 'query reads.fasta'
            reference = root / 'reference db.fasta'
            query.write_text('>ISO1\nACGT\n')
            reference.write_text('>REF1\nACGT\n')

            def fake_run(command):
                tokens = shlex.split(command)
                self.assertIn(str(query), tokens)
                self.assertIn(str(reference), tokens)
                output = Path(tokens[tokens.index('--blast6out') + 1])
                output.write_text('ISO1\tREF1\t99\t4\t0\t0\t1\t4\t1\t4\t1\t0\n')

            with mock.patch('branchmanager.pipeline.tree.run_cmd', side_effect=fake_run):
                oriented, rows = _orient_tree_input_fasta(
                    str(query), ref_fasta=str(reference), anchor_fasta=None,
                    outdir=root, label='queries',
                )
            self.assertEqual(oriented, str(query))
            self.assertEqual(rows[0]['status'], 'forward_kept')

    def test_registered_filing_cabinet_is_automatic_novelty_baseline(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            db = Database(str(Path(tmpdir) / 'project.sqlite'))
            db.initialise()
            db.upsert_dataset_role('Hungate', 'baseline', genomes_available=True)
            args = Namespace(
                baseline_dataset='Baseline', baseline_fasta=None,
                novelty_baseline_datasets=[],
            )
            self.assertEqual(_resolve_novelty_baseline_datasets(db, args), ['Hungate'])

    def test_office_workflow_commands_are_primary_cli_entries(self):
        parser = build_parser()
        cases = (
            (['mailroom', '--read-dir', 'reads', '--metadata', 'supplier.csv', '--dataset', 'QUB_01', '-o', 'out'], 'mailroom'),
            (['interview', '--mailroom', 'mailroom_out', '-o', 'out'], 'interview'),
            (['filing-cabinet', '--fasta', 'b.fa', '--db', 'p.db', '--dataset', 'Hungate'], 'filing-cabinet'),
            (['performance-review', '--input', 'q.fa', '--db', 'p.db', '--dataset', 'QUB', '-o', 'out'], 'performance-review'),
            (['quarterly-review', '--db', 'p.db', '--genome-budget', '3', '-o', 'out'], 'quarterly-review'),
            (['paper-trail', '--input', 'reads', '-o', 'out'], 'paper-trail'),
            (['onboarding', '--sample-map', 'map.tsv', '--partner-metadata', 'meta.tsv', '-o', 'out'], 'onboarding'),
            (['onboarding', '--fasta', 'markers.fasta', '--partner-metadata', 'meta.tsv', '-o', 'out'], 'onboarding'),
            (['status-meeting', '--db', 'p.db', '--input', 'status.tsv', '-o', 'out'], 'status-meeting'),
            (['records-update', '--db', 'p.db', '--input', 'genomes.tsv', '-o', 'out'], 'records-update'),
            (['exit-interview', '--db', 'p.db', '--sequence-id', 'ISO1', '--reason', 'withdrawn', '-o', 'out'], 'exit-interview'),
            (['annual-report', '--db', 'p.db', '-o', 'out'], 'annual-report'),
            (['it-desk'], 'it-desk'),
            (['assistant', '--sample-map', 'map.tsv', '--partner-metadata', 'meta.tsv', '--db', 'p.db', '--dataset', 'QUB', '--ref', 'gtdb.fa', '-o', 'out'], 'assistant'),
            (['assistant', '--fasta', 'markers.fasta', '--partner-metadata', 'meta.tsv', '--accept-unverified-marker-qc', '--db', 'p.db', '--dataset', 'QUB', '--ref', 'gtdb.fa', '-o', 'out'], 'assistant'),
            (['background-check', '--dataset', 'Hungate=b.fa', '--ref', 'gtdb.fa', '-o', 'out'], 'background-check'),
        )
        for argv, expected in cases:
            with self.subTest(command=expected):
                self.assertEqual(parser.parse_args(argv).command, expected)

    def test_retired_command_aliases_are_not_accepted(self):
        parser = build_parser()
        retired = (
            'preload', 'run', 'evaluate', 'eval', 'portfolio', 'sanger', 'ab1',
            'ab1-to-fasta', 'trace', 'intake', 'status-update', 'ingest-genomes',
            'doctor', 'preclassify', 'assistant-to-the-branch-manager',
            'subtree', 'regen-itol',
        )
        for command in retired:
            with self.subTest(command=command), self.assertRaises(SystemExit):
                parser.parse_args([command])

    def test_duplicate_read_metadata_input_is_not_exposed(self):
        parser = build_parser()
        with self.assertRaises(SystemExit):
            parser.parse_args([
                'paper-trail', '--read-metadata', 'reads.tsv', '-o', 'out',
            ])

    def test_onboarding_requires_cumulative_partner_metadata(self):
        parser = build_parser()
        with self.assertRaises(SystemExit):
            parser.parse_args(['onboarding', '--sample-map', 'batch.tsv', '-o', 'out'])

    def test_colour_options_use_british_spelling(self):
        parser = build_parser()
        filing_args = parser.parse_args([
            'filing-cabinet', '--fasta', 'b.fa', '--db', 'p.db',
            '--dataset', 'Hungate', '--colours', 'baseline.tsv',
        ])
        review_args = parser.parse_args([
            'performance-review', '--input', 'q.fa', '--db', 'p.db',
            '--dataset', 'QUB', '-o', 'out', '--baseline-colours', 'baseline.tsv',
            '--user-colours', 'partner.tsv',
        ])
        self.assertEqual(filing_args.colours, 'baseline.tsv')
        self.assertEqual(review_args.baseline_colours, 'baseline.tsv')
        self.assertEqual(review_args.user_colours, 'partner.tsv')
        paper_trail = parser.parse_args([
            'paper-trail', '--input', 'reads', '-o', 'out',
            '--max-report-image-height', '900',
        ])
        self.assertEqual(paper_trail.max_report_image_height, 900)


if __name__ == '__main__':
    unittest.main()
