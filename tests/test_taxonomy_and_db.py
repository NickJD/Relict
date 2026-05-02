import os
import re
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

from relict.cli import _find_preferred_id_map, _resolve_reference_inputs, build_parser, cmd_preload, cmd_run
from relict.db.interface import Database
from relict.pipeline import classify as classify_pipeline
from relict.pipeline import novelty as novelty_pipeline
from relict.pipeline.classify import _load_taxa_map_from_reference_fasta, validate_reference_taxonomy_consistency
from relict.pipeline.collapse import collapse_fasta_within_taxa
from relict.pipeline.itol import generate_itol_colors, write_dataset_membership_strip
from relict.pipeline.novelty import build_reference_novelty_metrics
from relict.pipeline.tree import (
    _orient_tree_input_fasta,
    _label_internal_nodes,
    _repair_legacy_internal_node_labels,
    build_combined_fasta,
    collect_tree_build_warnings,
    initialise_or_update_tree,
    summarize_alignment_quality,
)
from relict.utils.fasta import read_fasta, reverse_complement
from relict.pipeline.workflow_helpers import (
    _assignment_source_is_fasta,
    build_placement_warning_rows,
    build_sequence_assessment_rows,
    classification_ids_matching_kingdom,
    iter_assignment_rows,
    load_taxonomy_entries_from_assignments,
    merge_combined_taxonomy_rows,
)
from relict.taxonomy import (
    canonicalize_sequence_id,
    parse_taxon_string,
    taxonomy_matches_kingdom,
)


class TaxonomyUtilsTests(unittest.TestCase):
    def test_parse_taxon_string_supports_rank_prefixes(self):
        parsed = parse_taxon_string('d__Bacteria; p__Firmicutes; g__Bacillus')
        self.assertEqual(parsed['d'], 'Bacteria')
        self.assertEqual(parsed['p'], 'Firmicutes')
        self.assertEqual(parsed['g'], 'Bacillus')

    def test_taxonomy_matches_kingdom_is_rank_aware(self):
        self.assertTrue(taxonomy_matches_kingdom('d__Bacteria; p__Firmicutes', 'bacteria'))
        self.assertFalse(taxonomy_matches_kingdom('d__Archaea; p__Euryarchaeota', 'bacteria'))
        self.assertFalse(taxonomy_matches_kingdom('p__Bacteria_like; g__Foo', 'bacteria'))

    def test_canonicalize_sequence_id_preserves_asv_suffixes(self):
        self.assertEqual(canonicalize_sequence_id('ASV_1 some description'), 'ASV_1')
        self.assertEqual(canonicalize_sequence_id('read42:100-250(+)' ), 'read42')
        self.assertEqual(canonicalize_sequence_id('abc|ASV_2#fragment extra'), 'ASV_2')


class DatabaseBehaviorTests(unittest.TestCase):
    def test_taxonomy_replace_is_unique_per_dataset(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            db = Database(os.path.join(tmpdir, 'test.sqlite'))
            db.initialise()
            db.insert_sequences([('ASV_1', 'ACGT')], dataset='run1')
            db.insert_taxonomy([('ASV_1', 'd__Bacteria; p__Firmicutes', 0.9, 'run1')])
            db.insert_taxonomy([('ASV_1', 'd__Bacteria; p__Bacillota', 0.95, 'run1')])
            with db.connect() as conn:
                cur = conn.cursor()
                cur.execute("SELECT taxonomy, confidence, COUNT(*) OVER() FROM taxonomy WHERE id = ? AND dataset = ?", ('ASV_1', 'run1'))
                row = cur.fetchone()
            self.assertIsNotNone(row)
            self.assertEqual(row[0], 'd__Bacteria; p__Bacillota')
            self.assertEqual(row[1], 0.95)
            self.assertEqual(row[2], 1)

    def test_distance_replace_is_unique_per_dataset(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            db = Database(os.path.join(tmpdir, 'test.sqlite'))
            db.initialise()
            db.insert_sequences([('ASV_1', 'ACGT')], dataset='run1')
            db.insert_distances([('ASV_1', 'run1', 'REF1', 97.1)])
            db.insert_distances([('ASV_1', 'run1', 'REF2', 98.2)])
            with db.connect() as conn:
                cur = conn.cursor()
                cur.execute("SELECT nearest, identity, COUNT(*) OVER() FROM distances WHERE id = ? AND dataset = ?", ('ASV_1', 'run1'))
                row = cur.fetchone()
            self.assertIsNotNone(row)
            self.assertEqual(row[0], 'REF2')
            self.assertEqual(row[1], 98.2)
            self.assertEqual(row[2], 1)

    def test_classification_ids_matching_kingdom(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            class_tsv = Path(tmpdir) / 'taxonomy.tsv'
            class_tsv.write_text(
                'ID\tBestHit\tIdentity\tTaxon\tConfidence\n'
                'Q1\tR1\t99.0\td__Bacteria; p__Firmicutes\t0.9\n'
                'Q2\tR2\t95.0\td__Archaea; p__Euryarchaeota\t0.8\n'
                'Q3\tR3\t0.0\tNA\tNA\n'
            )
            matched = classification_ids_matching_kingdom(str(class_tsv), 'bacteria')
            self.assertEqual(matched, {'Q1'})

    def test_merge_combined_taxonomy_rows_overrides_existing_with_run_results(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            db = Database(os.path.join(tmpdir, 'test.sqlite'))
            db.initialise()
            class_tsv = Path(tmpdir) / 'taxonomy.tsv'
            class_tsv.write_text(
                'ID\tBestHit\tIdentity\tTaxon\tConfidence\n'
                'origA\tR1\t99.0\td__Bacteria; p__Bacillota\t0.97\n'
            )
            merged = merge_combined_taxonomy_rows(
                [('S01', 'd__Bacteria; p__Firmicutes', 0.90)],
                str(class_tsv),
                {'origA': 'S01', 'S01': 'S01'},
                db,
            )
            self.assertEqual(merged, [('S01', 'd__Bacteria; p__Bacillota', 0.97)])

    def test_validate_reference_taxonomy_consistency_reports_mismatch(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            ref_fasta = tmp / 'ref.fasta'
            taxa_tsv = tmp / 'taxa.tsv'
            ref_fasta.write_text('>REF_A\nACGT\n>REF_B\nACGT\n')
            taxa_tsv.write_text('FeatureID\tTaxon\tConfidence\nREF_A\td__Bacteria; p__Firmicutes\t0.9\nREF_X\td__Bacteria; p__Bacteroidota\t0.8\n')
            warnings = validate_reference_taxonomy_consistency(str(ref_fasta), str(taxa_tsv))
            cats = {w['category'] for w in warnings}
            self.assertIn('LOW_REFERENCE_TAXONOMY_OVERLAP', cats)
            self.assertIn('REFERENCE_ID_MISSING_TAXONOMY', cats)
            self.assertIn('TAXONOMY_ID_MISSING_REFERENCE', cats)

    def test_load_taxa_map_from_gtdb_reference_fasta_headers(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            ref_fasta = Path(tmpdir) / 'gtdb.fasta'
            ref_fasta.write_text(
                '>RS_GCF_031457235.1 d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria [locus_tag=abc]\nACGT\n'
                '>BAD_HEADER\nACGT\n'
            )
            taxa_map, warnings = _load_taxa_map_from_reference_fasta(str(ref_fasta))
            self.assertIn('RS_GCF_031457235.1', taxa_map)
            self.assertEqual(taxa_map['RS_GCF_031457235.1'][0], 'd__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria')
            self.assertTrue(any(w['category'] == 'REFERENCE_HEADER_WITHOUT_TAXONOMY' for w in warnings))

    def test_iter_assignment_rows_parses_taxonomy_from_same_fasta(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta = Path(tmpdir) / 'gtdb.fasta'
            fasta.write_text(
                '>RS_GCF_031457235.1 d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria [locus_tag=abc]\nACGT\n'
            )
            rows = list(iter_assignment_rows(str(fasta), source_fasta_path=str(fasta)))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]['qid'], 'RS_GCF_031457235.1 d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria [locus_tag=abc]')
            self.assertEqual(rows[0]['tax'], 'd__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria')

    def test_iter_assignment_rows_parses_taxonomy_from_external_fasta(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta = Path(tmpdir) / 'gtdb_assignments.fasta'
            fasta.write_text(
                '>RS_GCF_031457235.1 d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria [locus_tag=abc]\nACGT\n'
            )
            rows = list(iter_assignment_rows(str(fasta), source_fasta_path=str(Path(tmpdir) / 'different_input.fasta')))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]['tax'], 'd__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria')

    def test_load_taxonomy_entries_from_assignments_accepts_same_fasta(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            db = Database(os.path.join(tmpdir, 'test.sqlite'))
            db.initialise()
            fasta = Path(tmpdir) / 'gtdb.fasta'
            header = 'RS_GCF_031457235.1 d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria [locus_tag=abc]'
            fasta.write_text(f'>{header}\nACGT\n')
            entries = load_taxonomy_entries_from_assignments(
                str(fasta),
                {header: 'S01'},
                db,
                'run1',
                source_fasta_path=str(fasta),
            )
            self.assertEqual(entries, [('S01', 'd__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria', None, 'run1')])

    def test_load_taxonomy_entries_from_assignments_accepts_external_fasta(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            db = Database(os.path.join(tmpdir, 'test.sqlite'))
            db.initialise()
            fasta = Path(tmpdir) / 'gtdb_assignments.fasta'
            header = 'RS_GCF_031457235.1 d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria [locus_tag=abc]'
            fasta.write_text(f'>{header}\nACGT\n')
            entries = load_taxonomy_entries_from_assignments(
                str(fasta),
                {header: 'S01'},
                db,
                'run1',
                source_fasta_path=str(Path(tmpdir) / 'different_input.fasta'),
            )
            self.assertEqual(entries, [('S01', 'd__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria', None, 'run1')])

    def test_assignment_source_is_fasta_supports_gzipped_reference(self):
        import gzip

        with tempfile.TemporaryDirectory() as tmpdir:
            fasta = Path(tmpdir) / 'gtdb.fasta.gz'
            with gzip.open(fasta, 'wt') as fh:
                fh.write('>RS_GCF_031457235.1 d__Bacteria;p__Pseudomonadota\nACGT\n')
            self.assertTrue(_assignment_source_is_fasta(str(fasta), source_fasta_path=str(Path(tmpdir) / 'query.fasta')))

    def test_resolve_reference_inputs_treats_external_fasta_assignments_as_reference(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            query = Path(tmpdir) / 'query.fasta'
            query.write_text('>Q1\nACGT\n')
            gtdb = Path(tmpdir) / 'gtdb.fasta'
            gtdb.write_text('>REF1 d__Bacteria;p__Firmicutes\nACGT\n')
            taxa_tsv = Path(tmpdir) / 'taxa.tsv'
            taxa_tsv.write_text('FeatureID\tTaxon\nREF1\td__Bacteria; p__Firmicutes\n')

            ref, taxa, assignments = _resolve_reference_inputs(
                None,
                str(taxa_tsv),
                str(gtdb),
                source_fasta_path=str(query),
                log_prefix='[TEST]',
            )

            self.assertEqual(ref, str(gtdb))
            self.assertEqual(taxa, str(taxa_tsv))
            self.assertIsNone(assignments)

    def test_resolve_reference_inputs_keeps_same_fasta_as_direct_assignments(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            query = Path(tmpdir) / 'query.fasta'
            query.write_text('>Q1 d__Bacteria;p__Firmicutes\nACGT\n')

            ref, taxa, assignments = _resolve_reference_inputs(
                None,
                None,
                str(query),
                source_fasta_path=str(query),
                log_prefix='[TEST]',
            )

            self.assertIsNone(ref)
            self.assertIsNone(taxa)
            self.assertEqual(assignments, str(query))

    def test_build_parser_accepts_taxa_aasignments_alias(self):
        parser = build_parser()
        args = parser.parse_args([
            'preload',
            '--fasta', 'input.fasta',
            '--db', 'db.sqlite',
            '--dataset', 'Hungate',
            '--taxa-aasignments', 'gtdb.fna',
        ])
        self.assertEqual(args.taxa_assignments, 'gtdb.fna')

    def test_build_parser_accepts_no_shorten_ids(self):
        parser = build_parser()
        args = parser.parse_args([
            'run',
            '--input', 'input.fasta',
            '--db', 'db.sqlite',
            '--out', 'outdir',
            '--dataset', 'Run1',
            '--ref', 'ref.fasta',
            '--no-shorten-ids',
        ])
        self.assertFalse(args.shorten_ids)

    def test_find_preferred_id_map_prefers_preload_specific_mapping(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            (tmp / 'user_id_map.tsv').write_text('short_id\toriginal_header\nUSR01\torig_user\n')
            (tmp / 'preload_id_map.tsv').write_text('short_id\toriginal_header\nQUE06\torig_preload\n')
            self.assertEqual(_find_preferred_id_map(tmp), tmp / 'preload_id_map.tsv')

    def test_cmd_preload_writes_preload_id_map(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            fasta = tmp / 'preload.fasta'
            fasta.write_text('>query_alpha full original header\nACGT\n>query_beta another header\nTGCA\n')
            db_path = tmp / 'test.sqlite'
            outdir = tmp / 'preload_out'
            args = Namespace(
                fasta=str(fasta),
                db=str(db_path),
                out=str(outdir),
                dataset='Hungate',
                collapse=False,
                collapse_threshold=99.9,
                taxa=None,
                taxa_assignments=None,
                colors=None,
                classify=False,
                kingdom=None,
                build_tree=False,
                ref=None,
                threads=1,
                anchor_file=None,
            )

            cmd_preload(args)

            map_path = outdir / 'preload_id_map.tsv'
            self.assertTrue(map_path.exists())
            text = map_path.read_text()
            self.assertIn('short_id\toriginal_header', text)
            self.assertIn('query_alpha full original header', text)
            self.assertIn('query_beta another header', text)

    def test_preload_can_keep_canonical_ids_when_shortening_disabled(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            fasta = tmp / 'preload.fasta'
            fasta.write_text('>Alpha description one\nACGT\n>Beta description two\nTGCA\n')
            db = Database(os.path.join(tmpdir, 'test.sqlite'))
            db.initialise()

            alias_entries, mapped_fasta = db.preload_from_files(
                str(fasta),
                dataset='Hungate',
                outdir=str(tmp),
                shorten_ids=False,
            )

            self.assertEqual(alias_entries, [('Alpha', 'Alpha description one'), ('Beta', 'Beta description two')])
            headers = [h for h, _ in read_fasta(mapped_fasta)]
            self.assertEqual(headers, ['Alpha', 'Beta'])

    def test_preload_no_shorten_ids_raises_on_canonical_collision(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            fasta = tmp / 'preload.fasta'
            fasta.write_text('>Alpha description one\nACGT\n>Alpha description two\nTGCA\n')
            db = Database(os.path.join(tmpdir, 'test.sqlite'))
            db.initialise()

            with self.assertRaises(ValueError):
                db.preload_from_files(str(fasta), dataset='Hungate', outdir=str(tmp), shorten_ids=False)

    def test_cmd_run_uses_external_fasta_taxa_assignments_as_reference(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            input_fasta = tmp / 'input.fasta'
            input_fasta.write_text('>origA\nACGT\n')
            ref_fasta = tmp / 'gtdb.fasta'
            ref_fasta.write_text('>REF1 d__Bacteria;p__Firmicutes\nACGT\n')
            outdir = tmp / 'out'
            db_path = tmp / 'test.sqlite'

            def fake_run_classification(input_path, out_path, ref_fasta=None, taxa_tsv=None, threads=None):
                headers = [h for h, _ in read_fasta(input_path)]
                class_path = Path(out_path) / 'taxonomy.tsv'
                class_path.write_text(
                    'ID\tBestHit\tIdentity\tTaxon\tConfidence\n'
                    f'{headers[0]}\tREF1\t99.0\td__Bacteria; p__Firmicutes\t0.95\n'
                )
                return str(class_path)

            def fake_run_novelty(input_path, ref_path, out_path, id_threshold=0.97, threads=None):
                headers = [h for h, _ in read_fasta(input_path)]
                novelty_path = Path(out_path) / 'novelty.tsv'
                novelty_path.write_text(
                    'ID\tNearestIdentity\tNearestHit\tNovel\n'
                    f'{headers[0]}\t99.0\tREF1\tFalse\n'
                )
                return str(novelty_path)

            def fake_build_reference_novelty_metrics(input_path, ref_path, novelty_out, out_path, threads=None):
                headers = [h for h, _ in read_fasta(input_path)]
                metrics_path = Path(out_path) / 'novelty_metrics.tsv'
                metrics_path.write_text(
                    'ID\tNearestIdentity\tNearestHit\tNovel\tMatchesGE99\tMatchesGE97\tMatchesGE95\tNoveltyScore\tCrowding\tSequencingPriority\n'
                    f'{headers[0]}\t99.00\tREF1\tFalse\t1\t1\t1\t1.00\tdense\tLOW\n'
                )
                return str(metrics_path)

            args = Namespace(
                input=str(input_fasta),
                db=str(db_path),
                out=str(outdir),
                ref=None,
                taxa=None,
                dataset='run1',
                min_len=1,
                max_n=10,
                threads=1,
                preload_dir=None,
                kingdom=None,
                collapse=False,
                collapse_threshold=99.8,
                taxa_assignments=str(ref_fasta),
                user_colors=None,
                anchor_file=None,
            )

            with mock.patch('relict.cli.qc.run_qc', side_effect=lambda *a, **k: str(input_fasta)), \
                 mock.patch('relict.cli.derep.run_derep', side_effect=lambda *a, **k: str(input_fasta)), \
                 mock.patch('relict.cli.classify.run_classification', side_effect=fake_run_classification) as classify_mock, \
                 mock.patch('relict.cli.novelty.run_novelty', side_effect=fake_run_novelty) as novelty_mock, \
                 mock.patch('relict.cli.novelty.build_reference_novelty_metrics', side_effect=fake_build_reference_novelty_metrics), \
                 mock.patch('relict.cli.tree.initialise_or_update_tree'), \
                 mock.patch('relict.cli.tree.collect_tree_build_warnings', return_value=[]), \
                 mock.patch('relict.cli.tree.summarize_alignment_quality', return_value=[]), \
                 mock.patch('relict.cli.itol.generate_itol_colors'), \
                 mock.patch('relict.cli.novelty.build_run_novelty_itol'):
                cmd_run(args)

            self.assertEqual(classify_mock.call_args.kwargs['ref_fasta'], str(ref_fasta))
            self.assertEqual(novelty_mock.call_args.args[1], str(ref_fasta))

            db = Database(str(db_path))
            db.initialise()
            with db.connect() as conn:
                cur = conn.cursor()
                cur.execute('SELECT taxonomy FROM taxonomy WHERE dataset = ?', ('run1',))
                rows = cur.fetchall()
            self.assertEqual(rows, [('d__Bacteria; p__Firmicutes',)])

    def test_cmd_run_can_keep_canonical_ids_when_shortening_disabled(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            input_fasta = tmp / 'input.fasta'
            input_fasta.write_text('>origA extra words\nACGT\n')
            ref_fasta = tmp / 'ref.fasta'
            ref_fasta.write_text('>REF1 d__Bacteria;p__Firmicutes\nACGT\n')
            outdir = tmp / 'out'
            db_path = tmp / 'test.sqlite'

            def fake_run_classification(input_path, out_path, ref_fasta=None, taxa_tsv=None, threads=None):
                headers = [h for h, _ in read_fasta(input_path)]
                class_path = Path(out_path) / 'taxonomy.tsv'
                class_path.write_text(
                    'ID\tBestHit\tIdentity\tTaxon\tConfidence\n'
                    f'{headers[0]}\tREF1\t99.0\td__Bacteria; p__Firmicutes\t0.95\n'
                )
                return str(class_path)

            def fake_run_novelty(input_path, ref_path, out_path, id_threshold=0.97, threads=None):
                headers = [h for h, _ in read_fasta(input_path)]
                novelty_path = Path(out_path) / 'novelty.tsv'
                novelty_path.write_text(
                    'ID\tNearestIdentity\tNearestHit\tNovel\n'
                    f'{headers[0]}\t99.0\tREF1\tFalse\n'
                )
                return str(novelty_path)

            def fake_build_reference_novelty_metrics(input_path, ref_path, novelty_out, out_path, threads=None):
                headers = [h for h, _ in read_fasta(input_path)]
                metrics_path = Path(out_path) / 'novelty_metrics.tsv'
                metrics_path.write_text(
                    'ID\tNearestIdentity\tNearestHit\tNovel\tMatchesGE99\tMatchesGE97\tMatchesGE95\tNoveltyScore\tCrowding\tSequencingPriority\n'
                    f'{headers[0]}\t99.00\tREF1\tFalse\t1\t1\t1\t1.00\tdense\tLOW\n'
                )
                return str(metrics_path)

            args = Namespace(
                input=str(input_fasta),
                db=str(db_path),
                out=str(outdir),
                ref=str(ref_fasta),
                taxa=None,
                dataset='run1',
                min_len=1,
                max_n=10,
                threads=1,
                preload_dir=None,
                kingdom=None,
                collapse=False,
                collapse_threshold=99.8,
                taxa_assignments=None,
                user_colors=None,
                anchor_file=None,
                shorten_ids=False,
            )

            with mock.patch('relict.cli.qc.run_qc', side_effect=lambda *a, **k: str(input_fasta)), \
                 mock.patch('relict.cli.derep.run_derep', side_effect=lambda *a, **k: str(input_fasta)), \
                 mock.patch('relict.cli.classify.run_classification', side_effect=fake_run_classification), \
                 mock.patch('relict.cli.novelty.run_novelty', side_effect=fake_run_novelty), \
                 mock.patch('relict.cli.novelty.build_reference_novelty_metrics', side_effect=fake_build_reference_novelty_metrics), \
                 mock.patch('relict.cli.tree.initialise_or_update_tree'), \
                 mock.patch('relict.cli.tree.collect_tree_build_warnings', return_value=[]), \
                 mock.patch('relict.cli.tree.summarize_alignment_quality', return_value=[]), \
                 mock.patch('relict.cli.itol.generate_itol_colors'), \
                 mock.patch('relict.cli.novelty.build_run_novelty_itol'):
                cmd_run(args)

            headers = [h for h, _ in read_fasta(outdir / 'derep_short.fasta')]
            self.assertEqual(headers, ['origA'])
            mapping = (outdir / 'user_id_map.tsv').read_text()
            self.assertIn('origA\torigA extra words', mapping)

    def test_cmd_run_passes_preload_dir_to_tree_builder(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            input_fasta = tmp / 'input.fasta'
            input_fasta.write_text('>origA\nACGT\n')
            ref_fasta = tmp / 'ref.fasta'
            ref_fasta.write_text('>REF1 d__Bacteria;p__Firmicutes\nACGT\n')
            outdir = tmp / 'out'
            preload_dir = tmp / 'preload'
            preload_dir.mkdir()
            db_path = tmp / 'test.sqlite'

            def fake_run_classification(input_path, out_path, ref_fasta=None, taxa_tsv=None, threads=None):
                headers = [h for h, _ in read_fasta(input_path)]
                class_path = Path(out_path) / 'taxonomy.tsv'
                class_path.write_text(
                    'ID\tBestHit\tIdentity\tTaxon\tConfidence\n'
                    f'{headers[0]}\tREF1\t99.0\td__Bacteria; p__Firmicutes\t0.95\n'
                )
                return str(class_path)

            def fake_run_novelty(input_path, ref_path, out_path, id_threshold=0.97, threads=None):
                headers = [h for h, _ in read_fasta(input_path)]
                novelty_path = Path(out_path) / 'novelty.tsv'
                novelty_path.write_text(
                    'ID\tNearestIdentity\tNearestHit\tNovel\n'
                    f'{headers[0]}\t99.0\tREF1\tFalse\n'
                )
                return str(novelty_path)

            def fake_build_reference_novelty_metrics(input_path, ref_path, novelty_out, out_path, threads=None):
                headers = [h for h, _ in read_fasta(input_path)]
                metrics_path = Path(out_path) / 'novelty_metrics.tsv'
                metrics_path.write_text(
                    'ID\tNearestIdentity\tNearestHit\tNovel\tMatchesGE99\tMatchesGE97\tMatchesGE95\tNoveltyScore\tCrowding\tSequencingPriority\n'
                    f'{headers[0]}\t99.00\tREF1\tFalse\t1\t1\t1\t1.00\tdense\tLOW\n'
                )
                return str(metrics_path)

            args = Namespace(
                input=str(input_fasta),
                db=str(db_path),
                out=str(outdir),
                ref=str(ref_fasta),
                taxa=None,
                dataset='run1',
                min_len=1,
                max_n=10,
                threads=1,
                preload_dir=str(preload_dir),
                kingdom=None,
                collapse=False,
                collapse_threshold=99.8,
                taxa_assignments=None,
                user_colors=None,
                anchor_file=None,
                shorten_ids=True,
            )

            with mock.patch('relict.cli.qc.run_qc', side_effect=lambda *a, **k: str(input_fasta)), \
                 mock.patch('relict.cli.derep.run_derep', side_effect=lambda *a, **k: str(input_fasta)), \
                 mock.patch('relict.cli.classify.run_classification', side_effect=fake_run_classification), \
                 mock.patch('relict.cli.novelty.run_novelty', side_effect=fake_run_novelty), \
                 mock.patch('relict.cli.novelty.build_reference_novelty_metrics', side_effect=fake_build_reference_novelty_metrics), \
                 mock.patch('relict.cli.tree.initialise_or_update_tree') as tree_mock, \
                 mock.patch('relict.cli.tree.collect_tree_build_warnings', return_value=[]), \
                 mock.patch('relict.cli.tree.summarize_alignment_quality', return_value=[]), \
                 mock.patch('relict.cli.itol.generate_itol_colors'), \
                 mock.patch('relict.cli.novelty.build_run_novelty_itol'):
                cmd_run(args)

            self.assertEqual(tree_mock.call_args.kwargs['preload_dir'], str(preload_dir))


class OutputHelperTests(unittest.TestCase):
    def test_orient_tree_input_fasta_reverse_complements_minus_strand_hits(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            ref = tmp / 'ref.fasta'
            query = tmp / 'query.fasta'
            ref.write_text('>REF1\nAGTC\n')
            query.write_text('>Q1\nGACT\n')

            def fake_run_cmd(cmd):
                m = re.search(r'--blast6out\s+(\S+)', cmd)
                self.assertIsNotNone(m)
                Path(m.group(1)).write_text('Q1\tREF1\t99.0\t4\t0\t0\t4\t1\t1\t4\t-1\t0\n')

            with mock.patch('relict.pipeline.tree.run_cmd', side_effect=fake_run_cmd):
                oriented, rows = _orient_tree_input_fasta(
                    str(query),
                    ref_fasta=str(ref),
                    anchor_fasta=None,
                    outdir=tmp,
                    label='query',
                    threads=1,
                )

            records = list(read_fasta(oriented))
            self.assertEqual(records, [('Q1', 'AGTC')])
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]['status'], 'reverse_flipped')

    def test_build_combined_fasta_orients_db_and_user_sequences(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            ref = tmp / 'ref.fasta'
            user = tmp / 'user.fasta'
            ref.write_text('>REF1\nAAAA\n')
            user.write_text('>U1\nGACT\n')

            class FakeDb:
                def get_sequences_fasta(self, outpath, dataset=None):
                    Path(outpath).write_text('>D1\nTTTT\n')
                    return True

            def fake_run_cmd(cmd):
                m = re.search(r'--blast6out\s+(\S+)', cmd)
                self.assertIsNotNone(m)
                matches_path = Path(m.group(1))
                input_path = Path(cmd.split()[2])
                hits = []
                for header, _seq in read_fasta(str(input_path)):
                    if header == 'D1':
                        hits.append('D1\tREF1\t99.0\t4\t0\t0\t4\t1\t1\t4\t-1\t0')
                    elif header == 'U1':
                        hits.append('U1\tREF1\t99.0\t4\t0\t0\t4\t1\t1\t4\t-1\t0')
                matches_path.write_text('\n'.join(hits) + ('\n' if hits else ''))

            with mock.patch('relict.pipeline.tree.run_cmd', side_effect=fake_run_cmd):
                combined = build_combined_fasta(
                    str(user),
                    tmp,
                    anchor_file=None,
                    db=FakeDb(),
                    db_dataset=None,
                    ref_fasta=str(ref),
                    threads=1,
                    orientation_summary_path=tmp / 'tree_orientation_summary.tsv',
                    build_mode='initial',
                )

            combined_records = dict(read_fasta(str(combined)))
            self.assertEqual(combined_records['D1'], 'AAAA')
            self.assertEqual(combined_records['U1'], reverse_complement('GACT'))
            summary = (tmp / 'tree_orientation_summary.tsv').read_text()
            self.assertIn('BuildMode\tSourceGroup\tSequenceID', summary)
            self.assertIn('initial\tdb_sequences\tD1\treverse\tTrue', summary)
            self.assertIn('initial\tuser_sequences\tU1\treverse\tTrue', summary)

    def test_initialise_or_update_tree_orients_new_sequences_before_addfragments(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            outdir = tmp / 'tree_out'
            outdir.mkdir()
            current_aln = outdir / 'current_alignment.fasta'
            current_aln.write_text('>BASE\nAAAA\n')
            user = tmp / 'user.fasta'
            user.write_text('>Q1\nGACT\n')
            ref = tmp / 'ref.fasta'
            ref.write_text('>REF1\nAGTC\n')

            def fake_run_cmd(cmd):
                m = re.search(r'--blast6out\s+(\S+)', cmd)
                self.assertIsNotNone(m)
                Path(m.group(1)).write_text('Q1\tREF1\t99.0\t4\t0\t0\t4\t1\t1\t4\t-1\t0\n')

            def fake_addfragments(new_fasta, backbone_aln, output_fasta, threads=4):
                records = list(read_fasta(str(new_fasta)))
                self.assertEqual(records, [('Q1', 'AGTC')])
                Path(output_fasta).write_text('>BASE\nAAAA\n>Q1\nAGTC\n')
                return True

            with mock.patch('relict.pipeline.tree.run_cmd', side_effect=fake_run_cmd), \
                 mock.patch('relict.pipeline.tree._run_mafft_addfragments', side_effect=fake_addfragments), \
                 mock.patch('relict.pipeline.tree._run_fasttree', return_value=True):
                initialise_or_update_tree(
                    ref_fasta=str(ref),
                    user_fasta=str(user),
                    outdir=str(outdir),
                    db=None,
                    db_dataset=None,
                    threads=1,
                )

            summary = (outdir / 'tree_orientation_summary.tsv').read_text()
            self.assertIn('incremental\tnew_sequences\tQ1\treverse\tTrue', summary)

    def test_initialise_or_update_tree_seeds_preload_alignment_and_uses_mafft_add_for_full_length_sequences(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            preload = tmp / 'preload'
            preload.mkdir()
            (preload / 'current_alignment.fasta').write_text('>BASE\n' + ('A' * 1500) + '\n')
            (preload / 'current_tree.nwk').write_text('(BASE:0.1);\n')
            outdir = tmp / 'run_out'
            user = tmp / 'user.fasta'
            user.write_text('>Q1\n' + ('A' * 1400) + '\n')
            ref = tmp / 'ref.fasta'
            ref.write_text('>REF1\n' + ('A' * 1500) + '\n')

            def fake_run_cmd(cmd):
                m = re.search(r'--blast6out\s+(\S+)', cmd)
                self.assertIsNotNone(m)
                Path(m.group(1)).write_text('Q1\tREF1\t99.0\t1400\t0\t0\t1\t1400\t1\t1500\t1\t0\n')

            def fake_add(new_fasta, backbone_aln, output_fasta, threads=4):
                self.assertTrue(Path(backbone_aln).exists())
                self.assertEqual(Path(backbone_aln).read_text(), (preload / 'current_alignment.fasta').read_text())
                Path(output_fasta).write_text('>BASE\n' + ('A' * 1500) + '\n>Q1\n' + ('A' * 1400) + '\n')
                return True

            with mock.patch('relict.pipeline.tree.run_cmd', side_effect=fake_run_cmd), \
                 mock.patch('relict.pipeline.tree._run_mafft_full', return_value=False) as full_mock, \
                 mock.patch('relict.pipeline.tree._run_mafft_add', side_effect=fake_add) as add_mock, \
                 mock.patch('relict.pipeline.tree._run_mafft_addfragments', return_value=False) as addfrag_mock, \
                 mock.patch('relict.pipeline.tree._run_fasttree', return_value=True):
                initialise_or_update_tree(
                    ref_fasta=str(ref),
                    user_fasta=str(user),
                    outdir=str(outdir),
                    db=None,
                    db_dataset=None,
                    threads=1,
                    preload_dir=str(preload),
                )

            self.assertFalse(full_mock.called)
            self.assertTrue(add_mock.called)
            self.assertFalse(addfrag_mock.called)
            self.assertTrue((outdir / 'current_alignment.fasta').exists())

    def test_run_classification_searches_both_strands(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            query = tmp / 'query.fasta'
            ref = tmp / 'ref.fasta'
            query.write_text('>Q1\nACGT\n')
            ref.write_text('>REF1 d__Bacteria;p__Firmicutes\nACGT\n')

            seen = {}

            def fake_run_cmd(cmd):
                seen['cmd'] = cmd
                (tmp / 'matches.tsv').write_text('')

            with mock.patch('relict.pipeline.classify.run_cmd', side_effect=fake_run_cmd):
                classify_pipeline.run_classification(str(query), str(tmp), ref_fasta=str(ref))

            self.assertIn('--strand both', seen['cmd'])

    def test_run_novelty_searches_both_strands(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            query = tmp / 'query.fasta'
            ref = tmp / 'ref.fasta'
            query.write_text('>Q1\nACGT\n')
            ref.write_text('>REF1 d__Bacteria;p__Firmicutes\nACGT\n')

            seen = {}

            def fake_run_cmd(cmd):
                seen['cmd'] = cmd
                (tmp / 'novelty_matches.tsv').write_text('')

            with mock.patch('relict.pipeline.novelty.run_cmd', side_effect=fake_run_cmd):
                novelty_pipeline.run_novelty(str(query), str(ref), str(tmp))

            self.assertIn('--strand both', seen['cmd'])

    def test_write_dataset_membership_strip_writes_expected_rows(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir) / 'membership.itol'
            write_dataset_membership_strip(
                str(out),
                ['A1', 'A2', 'A3'],
                {'A1': 'preload', 'A2': 'run1', 'A3': ''},
            )
            text = out.read_text()
            self.assertIn('DATASET_LABEL,Dataset membership', text)
            self.assertIn('A1,', text)
            self.assertIn('A2,', text)
            self.assertIn('A3,#cccccc', text)

    def test_collapse_helper_falls_back_when_vsearch_missing(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            taxa_groups = {
                'tax1': [('S1', 'ACGT'), ('S2', 'ACGT')],
                'tax2': [('S3', 'TGCA')],
            }
            with mock.patch('relict.pipeline.collapse.shutil.which', return_value=None):
                artifacts = collapse_fasta_within_taxa(
                    taxa_groups,
                    tmpdir,
                    'collapsed.fasta',
                    'collapsed.tsv',
                    'members.tsv',
                    threshold=99.0,
                    threads=1,
                    log_prefix='[TEST COLLAPSE]',
                )
            self.assertTrue(artifacts.collapsed_path.exists())
            fasta_text = artifacts.collapsed_path.read_text()
            self.assertIn('>S1', fasta_text)
            self.assertIn('>S2', fasta_text)
            self.assertIn('>S3', fasta_text)

    def test_internal_node_labeling_does_not_introduce_double_colons(self):
        newick = '((A:0.1,B:0.2)0.95:0.3,C:0.4);'
        labeled = _label_internal_nodes(newick)
        self.assertIn('NODE', labeled)
        self.assertNotIn('::', labeled)

    def test_repair_legacy_internal_node_labels(self):
        broken = '((A:0.1,B:0.2)NODE0001::0.3,C:0.4);'
        repaired = _repair_legacy_internal_node_labels(broken)
        self.assertEqual(repaired, '((A:0.1,B:0.2)NODE0001:0.3,C:0.4);')

    def test_generate_itol_colors_repairs_legacy_tree_labels(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            taxonomy_tsv = tmp / 'combined_taxonomy.tsv'
            taxonomy_tsv.write_text(
                'ID\tTaxon\tConfidence\n'
                'A\td__Bacteria; p__Firmicutes\t0.9\n'
                'B\td__Bacteria; p__Firmicutes\t0.9\n'
            )
            tree_path = tmp / 'current_tree.nwk'
            tree_path.write_text('((A:0.1,B:0.2)NODE0001::0.3,C:0.4);')
            generate_itol_colors(str(taxonomy_tsv), str(tmp), tree_file=str(tree_path))
            repaired = tree_path.read_text()
            self.assertIn('NODE0001:0.3', repaired)
            self.assertNotIn('NODE0001::0.3', repaired)

    def test_tree_warning_helpers_report_partial_sequences_and_alignment_issues(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            user_fasta = tmp / 'user.fasta'
            user_fasta.write_text('>A\nACGT\n>B\nACGTNNNN\n')
            with mock.patch('relict.pipeline.tree.load_anchor_sequences', return_value=[]):
                warnings = collect_tree_build_warnings(str(user_fasta), anchor_file='missing_anchor_file.fasta')
            cats = {w['category'] for w in warnings}
            self.assertIn('LOW_SEQUENCE_COUNT', cats)
            self.assertIn('PARTIAL_16S_SEQUENCES', cats)
            self.assertIn('MISSING_REFERENCE_ANCHORS', cats)

            aln = tmp / 'alignment.fasta'
            aln.write_text('>A\n----\n>B\nAC--\n')
            aln_warnings = summarize_alignment_quality(str(aln))
            aln_cats = {w['category'] for w in aln_warnings}
            self.assertIn('ALL_GAP_ALIGNMENT_ROWS', aln_cats)

    def test_build_placement_warning_rows_flags_low_support(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            db = Database(os.path.join(tmpdir, 'test.sqlite'))
            db.initialise()
            db.insert_sequences([('S01', 'ACGT')], dataset='run1')
            class_tsv = tmp / 'taxonomy.tsv'
            class_tsv.write_text(
                'ID\tBestHit\tIdentity\tTaxon\tConfidence\n'
                'origA\tREF1\t94.0\td__Bacteria; p__Firmicutes\t0.5\n'
            )
            novelty_tsv = tmp / 'novelty.tsv'
            novelty_tsv.write_text(
                'ID\tNearestIdentity\tNearestHit\tNovel\n'
                'origA\t96.0\tREF1\tTrue\n'
            )
            rows = build_placement_warning_rows(str(class_tsv), str(novelty_tsv), {'origA': 'S01'}, db)
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]['id'], 'S01')
            self.assertIn('LOW_CLASSIFICATION_IDENTITY', rows[0]['flags'])
            self.assertIn('LOW_CONFIDENCE', rows[0]['flags'])
            self.assertIn('LOW_NEAREST_IDENTITY', rows[0]['flags'])

    def test_build_reference_novelty_metrics_reports_density_and_priority(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            input_fasta = tmp / 'query.fasta'
            ref_fasta = tmp / 'ref.fasta'
            novelty_tsv = tmp / 'novelty.tsv'
            input_fasta.write_text('>Q1\nACGT\n')
            ref_fasta.write_text('>R1 d__Bacteria;p__Firmicutes\nACGT\n')
            novelty_tsv.write_text('ID\tNearestIdentity\tNearestHit\tNovel\nQ1\t96.5\tR1\tTrue\n')

            def fake_run_cmd(cmd):
                density = tmp / 'novelty_density_matches.tsv'
                density.write_text('Q1\tR1\t99.2\nQ1\tR2\t97.4\nQ1\tR3\t95.1\n')

            with mock.patch('relict.pipeline.novelty.run_cmd', side_effect=fake_run_cmd):
                out = build_reference_novelty_metrics(str(input_fasta), str(ref_fasta), str(novelty_tsv), str(tmp))

            text = Path(out).read_text()
            self.assertIn('NoveltyScore', text)
            self.assertIn('MatchesGE99', text)
            self.assertIn('Q1\t96.50\tR1\tTrue\t1\t2\t3', text)

    def test_build_sequence_assessment_rows_combines_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            db = Database(os.path.join(tmpdir, 'test.sqlite'))
            db.initialise()
            class_tsv = tmp / 'taxonomy.tsv'
            class_tsv.write_text('ID\tBestHit\tIdentity\tTaxon\tConfidence\norigA\tREF1\t98.0\td__Bacteria; p__Firmicutes\t0.9\n')
            novelty_metrics = tmp / 'novelty_metrics.tsv'
            novelty_metrics.write_text(
                'ID\tNearestIdentity\tNearestHit\tNovel\tMatchesGE99\tMatchesGE97\tMatchesGE95\tNoveltyScore\tCrowding\tSequencingPriority\n'
                'S01\t96.50\tREF1\tTrue\t1\t2\t3\t25.00\tsparse\tHIGH\n'
            )
            warning_rows = [{'id': 'S01', 'flags': 'LOW_NEAREST_IDENTITY'}]
            rows = build_sequence_assessment_rows(['S01'], str(class_tsv), str(novelty_metrics), warning_rows, {'origA': 'S01'}, db)
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]['id'], 'S01')
            self.assertEqual(rows[0]['best_hit'], 'REF1')
            self.assertEqual(rows[0]['sequencing_priority'], 'HIGH')
            self.assertEqual(rows[0]['placement_flags'], 'LOW_NEAREST_IDENTITY')


if __name__ == '__main__':
    unittest.main()

