import gzip
import csv
import os
import re
import struct
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

from branchmanager.cli import _find_preferred_id_map, _load_partner_metadata_for_run, _load_performance_review_baseline, _resolve_reference_inputs, _write_output_explanations, build_parser, cmd_filing_cabinet, cmd_performance_review
from branchmanager.db.interface import Database
from branchmanager import mailroom as mailroom_module
from branchmanager.pipeline import classify as classify_pipeline
from branchmanager.pipeline import novelty as novelty_pipeline
from branchmanager.pipeline import neighbourhood as neighbourhood_pipeline
from branchmanager.pipeline import cluster_report as cluster_report_pipeline
from branchmanager.pipeline import background_check as background_check_pipeline
from branchmanager.pipeline import mwl as mwl_pipeline
from branchmanager.pipeline import selection_sets as selection_sets_pipeline
from branchmanager.pipeline import quarterly_review as quarterly_review_pipeline
from branchmanager.pipeline import qc as qc_pipeline
from branchmanager.pipeline import paper_trail as paper_trail_pipeline
from branchmanager.pipeline.classify import _load_taxa_map_from_reference_fasta, validate_reference_taxonomy_consistency
from branchmanager.pipeline.collapse import collapse_fasta_within_taxa
from branchmanager.pipeline.itol import generate_itol_colours, write_dataset_membership_strip
from branchmanager.pipeline.novelty import build_reference_novelty_metrics
from branchmanager.pipeline.tree import (
    _new_sequences_only,
    _orient_tree_input_fasta,
    _label_internal_nodes,
    _repair_internal_node_label_delimiters,
    build_combined_fasta,
    collect_tree_build_warnings,
    initialise_or_update_tree,
    summarise_alignment_quality,
)
from branchmanager.utils.fasta import read_fasta, reverse_complement
from branchmanager.pipeline.workflow_helpers import (
    _assignment_source_is_fasta,
    build_placement_warning_rows,
    build_selection_decision,
    build_sequence_assessment_rows,
    classification_ids_matching_kingdom,
    iter_assignment_rows,
    load_taxonomy_entries_from_assignments,
    merge_combined_taxonomy_rows,
    write_selection_summary_tsv,
)
from branchmanager.taxonomy import (
    canonicalise_sequence_id,
    parse_taxon_string,
    taxonomy_matches_kingdom,
)


def _write_minimal_ab1(path: Path, sequence: str, qualities):
    entries = []
    payload = bytearray()
    header_size = 4 + 2 + 28
    directory_offset = header_size
    directory_size = 2 * 28
    payload_offset = directory_offset + directory_size

    def add_entry(tag, tag_num, elem_type, elem_size, data_bytes):
        data_offset = payload_offset + len(payload)
        payload.extend(data_bytes)
        entries.append((
            tag.encode('ascii'),
            tag_num,
            elem_type,
            elem_size,
            len(data_bytes) // elem_size,
            len(data_bytes),
            data_offset,
            0,
        ))

    add_entry('PBAS', 2, 2, 1, sequence.encode('ascii'))
    add_entry('PCON', 2, 1, 1, bytes(qualities))

    root = struct.pack(
        '>4sIHHIIII',
        b'tdir',
        1,
        1023,
        28,
        len(entries),
        directory_size,
        directory_offset,
        0,
    )
    directory = b''.join(struct.pack('>4sIHHIIII', *entry) for entry in entries)
    path.write_bytes(b'ABIF' + struct.pack('>H', 1) + root + directory + payload)


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
        self.assertTrue(taxonomy_matches_kingdom('d__Eukaryota; k__Fungi; p__Ascomycota', 'fungi'))
        self.assertTrue(taxonomy_matches_kingdom('k__Fungi; p__Basidiomycota', 'fungal'))

    def test_canonicalise_sequence_id_preserves_asv_suffixes(self):
        self.assertEqual(canonicalise_sequence_id('ASV_1 some description'), 'ASV_1')
        self.assertEqual(canonicalise_sequence_id('read42:100-250(+)' ), 'read42')
        self.assertEqual(canonicalise_sequence_id('abc|ASV_2#fragment extra'), 'ASV_2')


class DatabaseBehaviourTests(unittest.TestCase):
    def test_database_initialise_creates_parent_directory(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            db_path = Path(tmpdir) / 'missing' / 'nested' / 'project.sqlite'
            db = Database(str(db_path))
            db.initialise()
            self.assertTrue(db_path.exists())

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

    def test_taxa_map_accepts_gzipped_csv_assignment_table(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            taxa_csv_gz = Path(tmpdir) / 'taxa.csv.gz'
            with gzip.open(taxa_csv_gz, 'wt') as fh:
                fh.write(
                    'sequence_id,taxonomy,confidence\n'
                    'TABLE_REF,d__Bacteria; p__Bacillota; g__Blautia,0.91\n'
                )

            taxa_map = classify_pipeline._load_taxa_map(str(taxa_csv_gz))
            tax, conf = classify_pipeline._lookup_tax('TABLE_REF', taxa_map)

            self.assertEqual(tax, 'd__Bacteria; p__Bacillota; g__Blautia')
            self.assertEqual(conf, 0.91)

    def test_load_taxa_map_from_gzipped_reference_fasta_headers(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            ref_fasta_gz = Path(tmpdir) / 'gtdb.fasta.gz'
            with gzip.open(ref_fasta_gz, 'wt') as fh:
                fh.write(
                    '>REF_FASTA d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria [locus_tag=abc]\nACGT\n'
                )

            taxa_map, warnings = _load_taxa_map_from_reference_fasta(str(ref_fasta_gz))
            tax, conf = classify_pipeline._lookup_tax('REF_FASTA', taxa_map)

            self.assertEqual(tax, 'd__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria')
            self.assertIsNone(conf)
            self.assertFalse(warnings)

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

    def test_iter_assignment_rows_accepts_gzipped_csv(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            assignments = Path(tmpdir) / 'assignments.csv.gz'
            with gzip.open(assignments, 'wt') as fh:
                fh.write(
                    'id,taxon,confidence\n'
                    'S01,d__Bacteria; p__Bacillota,0.88\n'
                )

            rows = list(iter_assignment_rows(str(assignments), source_fasta_path=str(Path(tmpdir) / 'input.fasta')))

            self.assertEqual(rows, [{'qid': 'S01', 'tax': 'd__Bacteria; p__Bacillota', 'confidence': 0.88}])

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

    def test_build_parser_accepts_no_shorten_ids(self):
        parser = build_parser()
        args = parser.parse_args([
            'performance-review',
            '--input', 'input.fasta',
            '--db', 'db.sqlite',
            '--out', 'outdir',
            '--dataset', 'Run1',
            '--ref', 'ref.fasta',
            '--no-shorten-ids',
        ])
        self.assertFalse(args.shorten_ids)

    def test_build_parser_accepts_paper_trail_subcommand(self):
        parser = build_parser()
        args = parser.parse_args([
            'paper-trail',
            '--input', 'reads',
            '-o', 'paper_trail_out',
            '--min-quality', '25',
            '--min-overlap', '30',
            '--primer', '27F',
            '--primer', '907R',
            '--sample-map', 'reads.tsv',
        ])
        self.assertEqual(args.command, 'paper-trail')
        self.assertEqual(args.input, ['reads'])
        self.assertEqual(args.out, 'paper_trail_out')
        self.assertEqual(args.min_quality, 25)
        self.assertEqual(args.min_overlap, 30)
        self.assertEqual(args.min_length, 800)
        self.assertIsNone(args.min_read_length)
        self.assertEqual(args.primers, ['27F', '907R'])
        self.assertEqual(args.sample_map, 'reads.tsv')

    def test_build_parser_accepts_performance_review_mwl_flags(self):
        parser = build_parser()
        args = parser.parse_args([
            'performance-review',
            '--input', 'input.fasta',
            '--db', 'project.sqlite',
            '--out', 'outdir',
            '--dataset', 'PartnerA',
            '--ref', 'gtdb.fna',
            '--mwl', 'MWL.xlsx',
            '--mwl-sheet', 'MWL_V1',
            '--mwl-min-rank', 'order',
            '--partner-metadata', 'partner_metadata.tsv',
            '--baseline-fasta', 'hungate.fna',
            '--baseline-dataset', 'Hungate',
            '--novelty-baseline-dataset', 'CulturedSetB',
            '--baseline-taxa-assignments', 'hungate_taxonomy.tsv',
            '--baseline-skip-classify',
            '--no-baseline-shorten-ids',
            '--sequence-domain', 'archaea',
        ])
        self.assertEqual(args.command, 'performance-review')
        self.assertEqual(args.mwl, 'MWL.xlsx')
        self.assertEqual(args.mwl_sheet, 'MWL_V1')
        self.assertEqual(args.mwl_min_rank, 'order')
        self.assertEqual(args.partner_metadata, 'partner_metadata.tsv')
        self.assertEqual(args.baseline_fasta, 'hungate.fna')
        self.assertEqual(args.baseline_dataset, 'Hungate')
        self.assertEqual(args.novelty_baseline_datasets, ['CulturedSetB'])
        self.assertEqual(args.baseline_taxa_assignments, 'hungate_taxonomy.tsv')
        self.assertTrue(args.baseline_skip_classify)
        self.assertFalse(args.baseline_shorten_ids)
        self.assertEqual(args.sequence_domain, 'archaea')
        self.assertEqual(args.neighbourhood_format, 'png')
        self.assertEqual(args.pangenome_target, 3)
        self.assertEqual(args.candidate_set_size, 4)

        with self.assertRaises(SystemExit):
            parser.parse_args([
                'performance-review', '--input', 'input.fasta', '--db', 'project.sqlite',
                '--out', 'outdir', '--dataset', 'PartnerA',
                '--neighbourhood-format', 'pdf',
            ])

    def test_build_parser_accepts_quarterly_review(self):
        args = build_parser().parse_args([
            'quarterly-review', '--db', 'project.sqlite', '--out', 'quarterly_review_out',
            '--genome-budget', '24', '--backups-per-primary', '2',
            '--tree', 'current_tree.nwk', '--alignment', 'current_alignment.fasta',
            '--assessment', 'sequence_assessment.tsv',
        ])
        self.assertEqual(args.command, 'quarterly-review')
        self.assertEqual(args.genome_budget, 24)
        self.assertEqual(args.backups_per_primary, 2)
        self.assertEqual(args.pangenome_target, 3)
        self.assertEqual(args.alignment, 'current_alignment.fasta')
        self.assertEqual(args.assessment, ['sequence_assessment.tsv'])

    def test_find_preferred_id_map_prefers_filing_cabinet_mapping(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            (tmp / 'user_id_map.tsv').write_text('short_id\toriginal_header\nUSR01\torig_user\n')
            (tmp / 'filing_cabinet_id_map.tsv').write_text('short_id\toriginal_header\nQUE06\torig_baseline\n')
            self.assertEqual(_find_preferred_id_map(tmp), tmp / 'filing_cabinet_id_map.tsv')

    def test_filing_cabinet_writes_office_named_id_map(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            fasta = tmp / 'baseline.fasta'
            fasta.write_text('>query_alpha full original header\nACGT\n>query_beta another header\nTGCA\n')
            db_path = tmp / 'test.sqlite'
            outdir = tmp / 'filing_cabinet_out'
            args = Namespace(
                fasta=str(fasta),
                db=str(db_path),
                out=str(outdir),
                dataset='Hungate',
                collapse=False,
                collapse_threshold=99.9,
                taxa=None,
                taxa_assignments=None,
                colours=None,
                classify=False,
                kingdom=None,
                build_tree=False,
                ref=None,
                threads=1,
                anchor_file=None,
            )

            cmd_filing_cabinet(args)

            map_path = outdir / 'filing_cabinet_id_map.tsv'
            self.assertTrue(map_path.exists())
            text = map_path.read_text()
            self.assertIn('short_id\toriginal_header', text)
            self.assertIn('query_alpha full original header', text)
            self.assertIn('query_beta another header', text)

    def test_filing_cabinet_preserves_supplied_ids_when_shortening_disabled(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            fasta = tmp / 'baseline.fasta'
            fasta.write_text('>Alpha description one\nACGT\n>Beta description two\nTGCA\n')
            db = Database(os.path.join(tmpdir, 'test.sqlite'))
            db.initialise()

            alias_entries, mapped_fasta = db.register_filing_cabinet(
                str(fasta),
                dataset='Hungate',
                outdir=str(tmp),
                shorten_ids=False,
            )

            self.assertEqual(alias_entries, [('Alpha description one', 'Alpha description one'), ('Beta description two', 'Beta description two')])
            headers = [h for h, _ in read_fasta(mapped_fasta)]
            self.assertEqual(headers, ['Alpha description one', 'Beta description two'])

    def test_filing_cabinet_no_shorten_ids_raises_on_exact_id_collision(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            fasta = tmp / 'baseline.fasta'
            fasta.write_text('>Alpha description one\nACGT\n>Alpha description one\nTGCA\n')
            db = Database(os.path.join(tmpdir, 'test.sqlite'))
            db.initialise()

            with self.assertRaises(ValueError):
                db.register_filing_cabinet(str(fasta), dataset='Hungate', outdir=str(tmp), shorten_ids=False)

    def test_performance_review_uses_external_fasta_taxa_assignments_as_reference(self):
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

            def fake_run_novelty(input_path, _ref_path, out_path, id_threshold=0.97, threads=None, **kwargs):
                headers = [h for h, _ in read_fasta(input_path)]
                novelty_path = Path(out_path) / 'novelty.tsv'
                novelty_path.write_text(
                    'ID\tNearestIdentity\tNearestHit\tNovel\n'
                    f'{headers[0]}\t99.0\tREF1\tFalse\n'
                )
                return str(novelty_path)

            def fake_build_reference_novelty_metrics(input_path, _ref_path, out_path, threads=None, **kwargs):
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
                previous_review=None,
                kingdom=None,
                collapse=False,
                collapse_threshold=99.8,
                taxa_assignments=str(ref_fasta),
                user_colours=None,
                anchor_file=None,
            )

            with mock.patch('branchmanager.cli.qc.run_qc', side_effect=lambda *a, **k: str(input_fasta)), \
                 mock.patch('branchmanager.cli.derep.run_derep', side_effect=lambda *a, **k: str(input_fasta)), \
                 mock.patch('branchmanager.cli.classify.run_classification', side_effect=fake_run_classification) as classify_mock, \
                 mock.patch('branchmanager.cli.novelty.run_novelty', side_effect=fake_run_novelty) as novelty_mock, \
                 mock.patch('branchmanager.cli.novelty.build_reference_novelty_metrics', side_effect=fake_build_reference_novelty_metrics), \
                 mock.patch('branchmanager.cli.tree.initialise_or_update_tree'), \
                 mock.patch('branchmanager.cli.tree.collect_tree_build_warnings', return_value=[]), \
                 mock.patch('branchmanager.cli.tree.summarise_alignment_quality', return_value=[]), \
                 mock.patch('branchmanager.cli.itol.generate_itol_colours'), \
                 mock.patch('branchmanager.cli.novelty.build_run_novelty_itol'):
                cmd_performance_review(args)

            self.assertEqual(classify_mock.call_args.kwargs['ref_fasta'], str(ref_fasta))
            self.assertEqual(novelty_mock.call_args.args[1], str(ref_fasta))

            db = Database(str(db_path))
            db.initialise()
            with db.connect() as conn:
                cur = conn.cursor()
                cur.execute('SELECT taxonomy FROM taxonomy WHERE dataset = ?', ('run1',))
                rows = cur.fetchall()
            self.assertEqual(rows, [('d__Bacteria; p__Firmicutes',)])

    def test_performance_review_can_keep_canonical_ids_when_shortening_disabled(self):
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

            def fake_run_novelty(input_path, _ref_path, out_path, id_threshold=0.97, threads=None, **kwargs):
                headers = [h for h, _ in read_fasta(input_path)]
                novelty_path = Path(out_path) / 'novelty.tsv'
                novelty_path.write_text(
                    'ID\tNearestIdentity\tNearestHit\tNovel\n'
                    f'{headers[0]}\t99.0\tREF1\tFalse\n'
                )
                return str(novelty_path)

            def fake_build_reference_novelty_metrics(input_path, _ref_path, out_path, threads=None, **kwargs):
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
                previous_review=None,
                kingdom=None,
                collapse=False,
                collapse_threshold=99.8,
                taxa_assignments=None,
                user_colours=None,
                anchor_file=None,
                shorten_ids=False,
            )

            with mock.patch('branchmanager.cli.qc.run_qc', side_effect=lambda *a, **k: str(input_fasta)), \
                 mock.patch('branchmanager.cli.derep.run_derep', side_effect=lambda *a, **k: str(input_fasta)), \
                 mock.patch('branchmanager.cli.classify.run_classification', side_effect=fake_run_classification), \
                 mock.patch('branchmanager.cli.novelty.run_novelty', side_effect=fake_run_novelty), \
                 mock.patch('branchmanager.cli.novelty.build_reference_novelty_metrics', side_effect=fake_build_reference_novelty_metrics), \
                 mock.patch('branchmanager.cli.tree.initialise_or_update_tree'), \
                 mock.patch('branchmanager.cli.tree.collect_tree_build_warnings', return_value=[]), \
                 mock.patch('branchmanager.cli.tree.summarise_alignment_quality', return_value=[]), \
                 mock.patch('branchmanager.cli.itol.generate_itol_colours'), \
                 mock.patch('branchmanager.cli.novelty.build_run_novelty_itol'):
                cmd_performance_review(args)

            headers = [h for h, _ in read_fasta(outdir / 'intermediate' / 'derep_short.fasta')]
            self.assertEqual(headers, ['origA extra words'])
            mapping = (outdir / 'ids' / 'user_id_map.tsv').read_text()
            self.assertIn('origA extra words\torigA extra words', mapping)
            self.assertIn('origA extra words\torigA', mapping)
            self.assertTrue((outdir / 'assessment' / 'sequence_assessment.tsv').exists())
            self.assertTrue((outdir / 'baseline' / 'baseline_hits.tsv').exists())
            self.assertTrue((outdir / 'taxonomy' / 'ref.tsv').exists())

    def test_performance_review_defaults_to_bacterial_sequence_domain_filter(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            input_fasta = tmp / 'input.fasta'
            input_fasta.write_text('>B1\nACGT\n>A1\nTGCA\n')
            ref_fasta = tmp / 'ref.fasta'
            ref_fasta.write_text('>REF_B d__Bacteria;p__Bacillota\nACGT\n>REF_A d__Archaea;p__Euryarchaeota\nTGCA\n')
            outdir = tmp / 'out'
            db_path = tmp / 'test.sqlite'

            def fake_run_classification(input_path, out_path, ref_fasta=None, taxa_tsv=None, threads=None):
                class_path = Path(out_path) / 'taxonomy.tsv'
                with open(class_path, 'w') as handle:
                    handle.write('ID\tBestHit\tIdentity\tTaxon\tConfidence\n')
                    for h, _s in read_fasta(input_path):
                        if h == 'A1':
                            handle.write(f'{h}\tREF_A\t99.0\td__Archaea; p__Euryarchaeota\t0.95\n')
                        else:
                            handle.write(f'{h}\tREF_B\t99.0\td__Bacteria; p__Bacillota\t0.95\n')
                return str(class_path)

            def fake_run_novelty(input_path, _ref_path, out_path, id_threshold=0.97, threads=None, **kwargs):
                headers = [h for h, _ in read_fasta(input_path)]
                novelty_path = Path(out_path) / 'novelty.tsv'
                novelty_path.write_text(
                    'ID\tNearestIdentity\tNearestHit\tNovel\n' +
                    ''.join(f'{h}\t99.0\tREF_B\tFalse\n' for h in headers)
                )
                return str(novelty_path)

            def fake_build_reference_novelty_metrics(input_path, _ref_path, out_path, threads=None, **kwargs):
                headers = [h for h, _ in read_fasta(input_path)]
                metrics_path = Path(out_path) / 'novelty_metrics.tsv'
                metrics_path.write_text(
                    'ID\tNearestIdentity\tNearestHit\tNovel\tMatchesGE99\tMatchesGE97\tMatchesGE95\tNoveltyScore\tCrowding\tSequencingPriority\n' +
                    ''.join(f'{h}\t99.00\tREF_B\tFalse\t1\t1\t1\t1.00\tdense\tLOW\n' for h in headers)
                )
                return str(metrics_path)

            args = Namespace(
                command='run',
                input=str(input_fasta),
                db=str(db_path),
                out=str(outdir),
                ref=str(ref_fasta),
                taxa=None,
                dataset='run1',
                min_len=1,
                max_n=10,
                threads=1,
                previous_review=None,
                kingdom=None,
                sequence_domain=None,
                collapse=False,
                collapse_threshold=99.8,
                taxa_assignments=None,
                user_colours=None,
                anchor_file=None,
                shorten_ids=False,
            )

            with mock.patch('branchmanager.cli.qc.run_qc', side_effect=lambda *a, **k: str(input_fasta)), \
                 mock.patch('branchmanager.cli.derep.run_derep', side_effect=lambda *a, **k: str(input_fasta)), \
                 mock.patch('branchmanager.cli.classify.run_classification', side_effect=fake_run_classification), \
                 mock.patch('branchmanager.cli.novelty.run_novelty', side_effect=fake_run_novelty), \
                 mock.patch('branchmanager.cli.novelty.build_reference_novelty_metrics', side_effect=fake_build_reference_novelty_metrics), \
                 mock.patch('branchmanager.cli.tree.initialise_or_update_tree'), \
                 mock.patch('branchmanager.cli.tree.collect_tree_build_warnings', return_value=[]), \
                 mock.patch('branchmanager.cli.tree.summarise_alignment_quality', return_value=[]), \
                 mock.patch('branchmanager.cli.itol.generate_itol_colours'), \
                 mock.patch('branchmanager.cli.novelty.build_run_novelty_itol'):
                cmd_performance_review(args)

            db = Database(str(db_path))
            db.initialise()
            with db.connect() as conn:
                cur = conn.cursor()
                cur.execute('SELECT id FROM sequences WHERE dataset = ? ORDER BY id', ('run1',))
                rows = [r[0] for r in cur.fetchall()]
            self.assertEqual(rows, ['B1'])

    def test_performance_review_passes_previous_review_to_tree_builder(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            input_fasta = tmp / 'input.fasta'
            input_fasta.write_text('>origA\nACGT\n')
            ref_fasta = tmp / 'ref.fasta'
            ref_fasta.write_text('>REF1 d__Bacteria;p__Firmicutes\nACGT\n')
            outdir = tmp / 'out'
            previous_review = tmp / 'previous_review'
            previous_review.mkdir()
            db_path = tmp / 'test.sqlite'

            def fake_run_classification(input_path, out_path, ref_fasta=None, taxa_tsv=None, threads=None):
                headers = [h for h, _ in read_fasta(input_path)]
                class_path = Path(out_path) / 'taxonomy.tsv'
                class_path.write_text(
                    'ID\tBestHit\tIdentity\tTaxon\tConfidence\n'
                    f'{headers[0]}\tREF1\t99.0\td__Bacteria; p__Firmicutes\t0.95\n'
                )
                return str(class_path)

            def fake_run_novelty(input_path, _ref_path, out_path, id_threshold=0.97, threads=None, **kwargs):
                headers = [h for h, _ in read_fasta(input_path)]
                novelty_path = Path(out_path) / 'novelty.tsv'
                novelty_path.write_text(
                    'ID\tNearestIdentity\tNearestHit\tNovel\n'
                    f'{headers[0]}\t99.0\tREF1\tFalse\n'
                )
                return str(novelty_path)

            def fake_build_reference_novelty_metrics(input_path, _ref_path, out_path, threads=None, **kwargs):
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
                previous_review=str(previous_review),
                kingdom=None,
                collapse=False,
                collapse_threshold=99.8,
                taxa_assignments=None,
                user_colours=None,
                anchor_file=None,
                shorten_ids=True,
            )

            with mock.patch('branchmanager.cli.qc.run_qc', side_effect=lambda *a, **k: str(input_fasta)), \
                 mock.patch('branchmanager.cli.derep.run_derep', side_effect=lambda *a, **k: str(input_fasta)), \
                 mock.patch('branchmanager.cli.classify.run_classification', side_effect=fake_run_classification), \
                 mock.patch('branchmanager.cli.novelty.run_novelty', side_effect=fake_run_novelty), \
                 mock.patch('branchmanager.cli.novelty.build_reference_novelty_metrics', side_effect=fake_build_reference_novelty_metrics), \
                 mock.patch('branchmanager.cli.tree.initialise_or_update_tree') as tree_mock, \
                 mock.patch('branchmanager.cli.tree.collect_tree_build_warnings', return_value=[]), \
                 mock.patch('branchmanager.cli.tree.summarise_alignment_quality', return_value=[]), \
                 mock.patch('branchmanager.cli.itol.generate_itol_colours'), \
                 mock.patch('branchmanager.cli.novelty.build_run_novelty_itol'):
                cmd_performance_review(args)

            self.assertEqual(tree_mock.call_args.kwargs['previous_review'], str(previous_review))

    def test_performance_review_baseline_loads_sequences_and_taxonomy_before_run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            baseline = tmp / 'hungate.fasta'
            baseline.write_text('>HungateA source sequence\nACGTACGT\n')
            ref = tmp / 'gtdb.fasta'
            ref.write_text('>REF1 d__Bacteria;p__Bacillota\nACGTACGT\n')
            outdir = tmp / 'out'
            db_path = tmp / 'test.sqlite'
            db = Database(str(db_path))
            db.initialise()
            args = Namespace(
                baseline_fasta=str(baseline),
                baseline_dataset='Hungate',
                baseline_taxa_assignments=None,
                baseline_skip_classify=False,
                baseline_colours=None,
                baseline_shorten_ids=False,
                dataset='Batch1',
            )

            def fake_run_classification(input_path, out_path, ref_fasta=None, taxa_tsv=None, threads=None):
                headers = [h for h, _ in read_fasta(input_path)]
                class_path = Path(out_path) / 'taxonomy.tsv'
                class_path.write_text(
                    'ID\tBestHit\tIdentity\tTaxon\tConfidence\n'
                    f'{headers[0]}\tREF1\t99.0\td__Bacteria; p__Bacillota\t0.99\n'
                )
                return str(class_path)

            with mock.patch('branchmanager.cli.classify.run_classification', side_effect=fake_run_classification):
                baseline_out = _load_performance_review_baseline(args, db, str(outdir), str(ref), None, 1)

            self.assertEqual(baseline_out, str(outdir / 'filing_cabinet_baseline'))
            self.assertTrue((outdir / 'filing_cabinet_baseline' / 'baseline_id_map.tsv').exists())
            self.assertTrue((outdir / 'filing_cabinet_baseline' / 'baseline_combined_taxonomy.tsv').exists())

            with db.connect() as conn:
                cur = conn.cursor()
                cur.execute(
                    'SELECT s.id, s.dataset, t.taxonomy FROM sequences s '
                    'LEFT JOIN taxonomy t ON s.id = t.id WHERE s.dataset = ?',
                    ('Hungate',),
                )
                rows = cur.fetchall()
            self.assertEqual(rows, [('HungateA source sequence', 'Hungate', 'd__Bacteria; p__Bacillota')])

    def test_partner_metadata_loads_already_sequenced_status_for_run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            db = Database(str(tmp / 'test.sqlite'))
            db.initialise()
            db.insert_sequences([('IsoA supplied id', 'ACGT'), ('IsoB supplied id', 'ACGA')], dataset='PartnerA')
            metadata = tmp / 'partner_metadata.tsv'
            metadata.write_text(
                'sequence_id\tpartner_id\talready_sequenced\n'
                'IsoA supplied id\tQUB\tyes\n'
                'IsoB supplied id\tUoG\tno\n'
                'MissingIso\tQUB\tyes\n'
            )
            args = Namespace(
                partner_metadata=str(metadata),
                dataset='PartnerA',
                command='performance-review',
            )

            loaded = _load_partner_metadata_for_run(
                args,
                db,
                str(tmp),
                {
                    'IsoA supplied id': 'IsoA supplied id',
                    'IsoB supplied id': 'IsoB supplied id',
                },
                ['IsoA supplied id', 'IsoB supplied id'],
            )

            self.assertTrue(loaded['IsoA supplied id']['selected_for_wgs'])
            self.assertFalse(loaded['IsoB supplied id']['selected_for_wgs'])
            self.assertEqual(loaded['IsoA supplied id']['partner_id'], 'QUB')
            warnings = (tmp / 'partner_metadata_warnings.tsv').read_text()
            self.assertIn('MissingIso\tmetadata_id_not_found_in_project_database_or_current_run', warnings)

    def test_partner_metadata_can_update_an_isolate_from_an_earlier_run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            db = Database(str(tmp / 'test.sqlite'))
            db.initialise()
            db.insert_sequences([('OldIso', 'ACGT')], dataset='PartnerA_001')
            db.insert_sequences([('NewIso', 'ACGA')], dataset='PartnerB_001')
            db.upsert_dataset_role('PartnerA_001', 'candidate')
            db.upsert_dataset_role('PartnerB_001', 'candidate')
            metadata = tmp / 'partner_metadata.tsv'
            metadata.write_text(
                'sequence_id\tpartner_id\talready_sequenced\n'
                'OldIso\tQUB\tyes\n'
                'NewIso\tUoG\tno\n'
            )
            args = Namespace(
                partner_metadata=str(metadata),
                dataset='PartnerB_001',
                command='performance-review',
            )

            loaded = _load_partner_metadata_for_run(
                args,
                db,
                str(tmp),
                {'NewIso': 'NewIso'},
                ['NewIso'],
            )

            self.assertFalse(loaded['NewIso']['selected_for_wgs'])
            stored = db.get_sequencing_metadata_for_ids()
            self.assertTrue(stored['OldIso']['selected_for_wgs'])
            self.assertEqual(stored['OldIso']['dataset'], 'PartnerA_001')

    def test_partner_metadata_persists_pending_selection_without_claiming_a_genome(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            db = Database(str(tmp / 'test.sqlite'))
            db.initialise()
            db.insert_sequences([('PendingIso', 'ACGT')], dataset='QUB_01')
            db.upsert_dataset_role('QUB_01', 'candidate')
            metadata = tmp / 'project_metadata.tsv'
            metadata.write_text(
                'sequence_id\tpartner_id\tselected_for_genome_sequencing\talready_sequenced\n'
                'PendingIso\tQUB\tyes\tno\n'
            )
            args = Namespace(
                partner_metadata=str(metadata), dataset='QUB_01', command='performance-review',
            )

            _load_partner_metadata_for_run(
                args, db, str(tmp), {'PendingIso': 'PendingIso'}, ['PendingIso'],
            )
            stored = db.get_sequencing_metadata_for_ids(['PendingIso'])['PendingIso']

            self.assertTrue(stored['selected_for_sequencing'])
            self.assertFalse(stored['genome_available'])

    def test_database_rejects_reused_id_with_changed_sequence(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            db = Database(str(Path(tmpdir) / 'test.sqlite'))
            db.initialise()
            db.insert_sequences([('StableIso', 'ac gt')], dataset='PartnerA')

            self.assertEqual(db.assert_sequence_compatible('StableIso', 'ACGT'), 'StableIso')
            with self.assertRaisesRegex(ValueError, 'Rolling sequence IDs are immutable'):
                db.assert_sequence_compatible('StableIso', 'ACGA')

    def test_partner_metadata_requires_partner_acronym_column(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            metadata = tmp / 'partner_metadata.tsv'
            metadata.write_text(
                'sequence_id\talready_sequenced\n'
                'IsoA\tyes\n'
            )
            from branchmanager.partner_metadata import load_partner_sequencing_metadata

            with self.assertRaisesRegex(ValueError, 'partner acronym'):
                load_partner_sequencing_metadata(str(metadata))

    def test_partner_metadata_separates_selection_commitment_from_available_genome(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            metadata = Path(tmpdir) / 'partner_metadata.tsv'
            metadata.write_text(
                'sequence_id\tpartner_id\tselected_for_genome_sequencing\talready_sequenced\n'
                'PendingIso\tQUB\tyes\tno\n'
                'GenomeIso\tUoG\tyes\tyes\n'
            )
            from branchmanager.partner_metadata import load_partner_sequencing_metadata

            rows = {row['source_id']: row for row in load_partner_sequencing_metadata(metadata)}

            self.assertTrue(rows['PendingIso']['selected_for_sequencing'])
            self.assertFalse(rows['PendingIso']['genome_available'])
            self.assertFalse(rows['PendingIso']['selected_for_wgs'])
            self.assertTrue(rows['GenomeIso']['selected_for_sequencing'])
            self.assertTrue(rows['GenomeIso']['genome_available'])

    def test_partner_metadata_rejects_duplicate_isolate_rows(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            metadata = Path(tmpdir) / 'partner_metadata.csv'
            metadata.write_text(
                'sequence_id,partner_id,selected_for_genome_sequencing,already_sequenced\n'
                'IsoA,QUB,no,no\n'
                'IsoA,QUB,no,no\n'
            )
            from branchmanager.partner_metadata import load_partner_sequencing_metadata

            with self.assertRaisesRegex(ValueError, 'exactly one row per isolate'):
                load_partner_sequencing_metadata(metadata)

    def test_background_check_reports_high_quality_taxonomy_disagreements(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            query = tmp / 'query.fasta'
            ref = tmp / 'ref.fasta'
            taxa = tmp / 'taxa.tsv'
            outdir = tmp / 'out'

            query.write_text('>Q1\nACGT\n>Q2\nACGT\n')
            ref.write_text('>R1\nACGT\n>R2\nACGT\n>R3\nACGT\n')
            taxa.write_text(
                'FeatureID\tTaxon\tConfidence\n'
                'R1\td__Bacteria; p__Firmicutes; g__Alpha\t0.99\n'
                'R2\td__Bacteria; p__Bacteroidota; g__Beta\t0.98\n'
                'R3\td__Bacteria; p__Firmicutes; g__Alpha\t0.97\n'
            )

            def fake_run_vsearch_pass(**kwargs):
                Path(kwargs['matches_path']).write_text(
                    'Q1\tR1\t99.8\n'
                    'Q1\tR2\t99.7\n'
                    'Q1\tR3\t99.6\n'
                    'Q2\tR1\t99.9\n'
                    'Q2\tR3\t99.8\n'
                )
                return {
                    'Q1': ('R1', 99.8),
                    'Q2': ('R1', 99.9),
                }

            with mock.patch('branchmanager.pipeline.background_check._run_vsearch_pass', side_effect=fake_run_vsearch_pass):
                res = background_check_pipeline.classify_fasta(
                    fasta_path=str(query),
                    ref_fasta=str(ref),
                    outdir=str(outdir),
                    dataset_name='testset',
                    taxa_tsv=str(taxa),
                    threads=1,
                    low_confidence_threshold=0.97,
                    max_hits=3,
                )

            report = Path(res['taxonomic_disagreement_tsv'])
            text = report.read_text()
            self.assertEqual(res['n_taxonomic_disagreements'], 1)
            self.assertIn('Q1\tR1\t99.8', text)
            self.assertIn('p__Bacteroidota', text)
            self.assertNotIn('\nQ2\t', text)

    def test_mwl_annotation_scores_deepest_matching_rank(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            mwl_tsv = Path(tmpdir) / 'mwl.tsv'
            mwl_tsv.write_text(
                'MWL\tHierarchy (Domain : Phylum : Class : Order : Family : Genus)\tFunctional guild / Metabolic role in rumen fermentation\n'
                'MWL6\td__Bacteria; p__Bacillota_A; c__Clostridia; o__Lachnospirales (including f__Lachnospiraceae; g__Blautia)\tAlternative hydrogen sink\n'
                'MWL13\td__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Bacteroidales (including g__Prevotella)\tSuccinate producers\n'
            )
            entries = mwl_pipeline.load_mwl_entries(str(mwl_tsv))
            row = {
                'id': 'S01',
                'taxonomy': 'd__Bacteria; p__Bacillota_A; c__Clostridia; o__Lachnospirales; f__Lachnospiraceae; g__Blautia',
                'classification_identity': '99.0',
                'novelty_score': '42.0',
                'investigation_score': '80.0',
            }

            mwl_pipeline.annotate_assessment_rows([row], entries, min_rank='p')

            self.assertEqual(row['mwl_match'], 'Yes')
            self.assertEqual(row['mwl_id'], 'MWL6')
            self.assertEqual(row['mwl_matched_rank'], 'genus')
            self.assertEqual(row['mwl_score'], '89.10')
            self.assertEqual(row['evaluation_score'], '82.28')


class OutputHelperTests(unittest.TestCase):
    def test_mailroom_builds_paired_ab1_map_with_primer_provenance(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            reads = tmp / 'All_AB1'
            reads.mkdir()
            _write_minimal_ab1(reads / 'KKX994_123_123.ab1', 'ACGT' * 12, [35] * 48)
            _write_minimal_ab1(reads / 'KKY011_456_456.ab1', 'TGCA' * 12, [35] * 48)
            metadata = tmp / 'supplier.csv'
            metadata.write_text(
                'Sequencing ID,Isolate Number,Read\n'
                'KKX994,SW_0016,Forward\n'
                'KKY011,SW_0016,Reverse\n'
            )

            result = mailroom_module.prepare_ab1_map(
                reads,
                metadata,
                tmp / 'mailroom',
                dataset='UoG_01',
                forward_primer='63F',
                reverse_primer='1492R',
            )

            self.assertEqual(result['status'], 'PASS')
            self.assertEqual(result['mapped_reads'], 2)
            with open(result['ab1_map']) as handle:
                rows = list(csv.DictReader(handle, delimiter='\t'))
            self.assertEqual({row['primer'] for row in rows}, {'63F', '1492R'})
            self.assertEqual({row['direction'] for row in rows}, {'forward', 'reverse'})
            self.assertEqual({row['processing_mode'] for row in rows}, {'assemble'})
            self.assertEqual({row['primer_assignment'] for row in rows}, {'configured_for_batch'})
            self.assertTrue(all(row['dataset'] == 'UoG_01' for row in rows))

    def test_mailroom_requires_review_when_primer_is_not_supplied(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            reads = tmp / 'All_AB1'
            reads.mkdir()
            _write_minimal_ab1(reads / 'READ001_123.ab1', 'ACGT' * 12, [35] * 48)
            metadata = tmp / 'supplier.tsv'
            metadata.write_text(
                'sequencing_id\tisolate_number\tread\n'
                'READ001\tISO1\tForward\n'
            )

            result = mailroom_module.prepare_ab1_map(
                reads, metadata, tmp / 'mailroom', dataset='Batch1',
            )

            self.assertEqual(result['status'], 'REVIEW_REQUIRED')
            self.assertEqual(result['unresolved_primers'], 1)
            self.assertIn('UNRESOLVED_PRIMER', Path(result['report']).read_text())
            self.assertIn('\tunknown\tforward\tbest_read\tunresolved', Path(result['ab1_map']).read_text())

    def test_mailroom_reports_missing_and_unmapped_ab1_files(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            reads = tmp / 'All_AB1'
            reads.mkdir()
            _write_minimal_ab1(reads / 'EXTRA001_123.ab1', 'ACGT' * 12, [35] * 48)
            metadata = tmp / 'supplier.csv'
            metadata.write_text(
                'Sequencing ID,Isolate Number,Read,Primer\n'
                'MISSING001,ISO1,Forward,63F\n'
            )

            result = mailroom_module.prepare_ab1_map(
                reads, metadata, tmp / 'mailroom', dataset='Batch1',
            )

            self.assertEqual(result['status'], 'FAIL')
            report = Path(result['report']).read_text()
            self.assertIn('READ_FILE_NOT_FOUND', report)
            self.assertIn('UNMAPPED_AB1_FILE', report)

    def test_read_ab1_extracts_sequence_and_quality(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            ab1 = Path(tmpdir) / 'Iso001_27F.ab1'
            _write_minimal_ab1(ab1, 'NNACGTACGTNN', [5, 5, 35, 36, 37, 38, 39, 40, 30, 31, 5, 5])

            seq, qual = paper_trail_pipeline.read_ab1(ab1)

            self.assertEqual(seq, 'NNACGTACGTNN')
            self.assertEqual(qual[:4], [5, 5, 35, 36])

    def test_paper_trail_trims_orients_and_assembles_primer_reads(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            reads = tmp / 'reads'
            reads.mkdir()
            # Same isolate, forward and reverse-primer reads with an overlap.
            (reads / 'Iso001_27F.fasta').write_text('>Iso001_27F\nNNNNAAAACCCCGGGGTTTTAAAANNNN\n')
            reverse_read = reverse_complement('CCCCGGGGTTTTAAAACCCC')
            (reads / 'Iso001_907R.fasta').write_text(f'>Iso001_907R\nNNNN{reverse_read}NNNN\n')

            out = paper_trail_pipeline.run_paper_trail(
                [str(reads)],
                tmp / 'out',
                min_quality=20,
                window=4,
                min_length=8,
                min_overlap=8,
                min_overlap_identity=0.90,
            )

            assembled = dict(read_fasta(out['assembled_fasta']))
            self.assertIn('Iso001', assembled)
            self.assertEqual(assembled['Iso001'], 'AAAACCCCGGGGTTTTAAAACCCC')
            report = Path(out['assembly_tsv']).read_text()
            self.assertIn('Iso001\tassembled\t2\t2', report)
            placements = Path(out['assembly_placements_tsv']).read_text()
            self.assertIn('Iso001\tIso001_27F\t27F\tforward\tkept\tyes\t1\t20', placements)
            self.assertIn('Iso001\tIso001_907R\t907R\treverse\tkept\tyes\t5\t24', placements)
            self.assertEqual(Path(out['qc_policy_tsv']).name, 'paper_trail_qc_policy.tsv')
            self.assertEqual(Path(out['summary']).name, 'paper_trail_summary.txt')

    def test_paper_trail_gapped_fallback_assembles_indel_shifted_reads(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            reads = tmp / 'reads'
            reads.mkdir()
            (reads / 'IsoGap_27F.fasta').write_text(
                '>IsoGap_27F\nAAAAGTACCGTTAACGATCGTACGATCGTTAA\n'
            )
            reverse_oriented = 'GTACCGTTAAGCGATCGTACGATCGTTAACCCC'
            (reads / 'IsoGap_907R.fasta').write_text(
                f'>IsoGap_907R\n{reverse_complement(reverse_oriented)}\n'
            )

            out = paper_trail_pipeline.run_paper_trail(
                [str(reads)],
                tmp / 'out',
                min_quality=20,
                window=4,
                min_length=8,
                min_overlap=18,
                min_overlap_identity=0.85,
            )

            report = Path(out['assembly_tsv']).read_text()
            self.assertIn('IsoGap\tassembled\t2\t2', report)
            self.assertIn('biopython_pairwise_overlap', report)
            self.assertIn('ambiguous_overlap_conflicts_1', report)

    def test_paper_trail_uses_batch_map_when_filenames_do_not_encode_primer(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            fwd = tmp / 'well_A01.ab1'
            rev = tmp / 'well_A02.ab1'
            _write_minimal_ab1(fwd, 'AAAACCCCGGGGTTTTAAAA', [35] * 20)
            _write_minimal_ab1(rev, reverse_complement('CCCCGGGGTTTTAAAACCCC'), [35] * 20)
            meta = tmp / 'reads.tsv'
            meta.write_text(
                'file\tsequence_id\tprimer\tdirection\n'
                'well_A01.ab1\tIsoA\t27F\tforward\n'
                'well_A02.ab1\tIsoA\t907R\treverse\n'
            )

            out = paper_trail_pipeline.run_paper_trail(
                [str(tmp)],
                tmp / 'out',
                sample_map=str(meta),
                min_length=8,
                min_overlap=8,
                min_overlap_identity=0.90,
            )

            assembled = dict(read_fasta(out['assembled_fasta']))
            self.assertEqual(assembled['IsoA'], 'AAAACCCCGGGGTTTTAAAACCCC')

    def test_paper_trail_accepts_sample_map_and_writes_error_visuals(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            fwd = tmp / 'well_A01.ab1'
            rev = tmp / 'well_A02.ab1'
            forward_seq = 'NNNNAAAACCCCGGGGTTTTAAAANNNN'
            reverse_seq = 'NNNN' + reverse_complement('CCCCGGGGTTTTAAAACCCC') + 'NNNN'
            _write_minimal_ab1(fwd, forward_seq, [5] * 4 + [35] * 20 + [5] * 4)
            _write_minimal_ab1(rev, reverse_seq, [5] * 4 + [35] * 20 + [5] * 4)
            sample_map = tmp / 'sample_reads.tsv'
            sample_map.write_text(
                'isolate_id\t27F\t907R\n'
                'IsoMap\twell_A01.ab1\twell_A02.ab1\n'
            )

            out = paper_trail_pipeline.run_paper_trail(
                [],
                tmp / 'out',
                sample_map=str(sample_map),
                min_quality=20,
                window=4,
                min_length=8,
                min_overlap=8,
                min_overlap_identity=0.90,
            )

            assembled = dict(read_fasta(out['assembled_fasta']))
            self.assertEqual(assembled['IsoMap'], 'AAAACCCCGGGGTTTTAAAACCCC')
            read_qc = Path(out['read_qc_tsv']).read_text()
            self.assertIn('MeanRawErrorProbability', read_qc)
            per_base = Path(out['per_base_error_tsv']).read_text()
            self.assertIn('ErrorProbability', per_base)
            self.assertIn('left_trimmed', per_base)
            self.assertIn('right_trimmed', per_base)
            from PIL import Image
            for key in ('read_error_pngs', 'chromatogram_pngs', 'assembly_pngs'):
                self.assertEqual(len(out[key]), 1)
                image_path = Path(out[key][0])
                self.assertEqual(image_path.read_bytes()[:8], b'\x89PNG\r\n\x1a\n')
                with Image.open(image_path) as image:
                    self.assertEqual(image.format, 'PNG')
                    self.assertGreater(image.width, 1000)
                    self.assertGreater(image.height, 100)
                    if key == 'assembly_pngs':
                        self.assertGreater(image.height, 250)
            self.assertTrue(Path(out['visual_manifest_tsv']).exists())

    def test_paper_trail_paginates_all_tall_visual_reports(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            fasta = tmp / 'reads.fasta'
            fasta.write_text(''.join(
                f'>Iso{index:02d}_27F\n{"ACGT" * 10}\n'
                for index in range(12)
            ))
            output_dir = tmp / 'out'
            stale_dir = output_dir / 'visual_reports' / 'read_error_profiles'
            stale_dir.mkdir(parents=True)
            stale_page = stale_dir / 'read_error_profiles_page_999.png'
            stale_page.write_bytes(b'stale')
            obsolete_single_image = output_dir / 'read_error_profiles.png'
            obsolete_single_image.write_bytes(b'stale')

            out = paper_trail_pipeline.run_paper_trail(
                [fasta],
                output_dir,
                min_length=8,
                min_read_length=8,
                max_report_image_height=600,
            )

            self.assertEqual(len(out['read_error_pngs']), 4)
            self.assertEqual(len(out['chromatogram_pngs']), 6)
            self.assertEqual(len(out['assembly_pngs']), 4)
            self.assertFalse(stale_page.exists())
            self.assertFalse(obsolete_single_image.exists())
            from PIL import Image
            all_pages = (
                out['read_error_pngs']
                + out['chromatogram_pngs']
                + out['assembly_pngs']
            )
            for page in all_pages:
                with Image.open(page) as image:
                    self.assertEqual(image.format, 'PNG')
                    self.assertLessEqual(image.height, 600)
            with open(out['visual_manifest_tsv']) as handle:
                rows = list(csv.DictReader(handle, delimiter='\t'))
            self.assertEqual(len(rows), 14)
            self.assertTrue(all(int(row['HeightPixels']) <= 600 for row in rows))
            self.assertEqual({row['MaxHeightPixels'] for row in rows}, {'600'})

    def test_paper_trail_sample_map_best_read_selects_highest_quality_read(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            low = tmp / 'low_quality.ab1'
            high = tmp / 'high_quality.ab1'
            _write_minimal_ab1(low, 'AAAAAAAAAAAA', [18] * 12)
            _write_minimal_ab1(high, 'CCCCCCCCCC', [40] * 10)
            sample_map = tmp / 'sample_reads.tsv'
            sample_map.write_text(
                'isolate_id\tab1_files\tprocessing_mode\n'
                'IsoBest\tlow_quality.ab1;high_quality.ab1\tbest_read\n'
            )

            out = paper_trail_pipeline.run_paper_trail(
                [],
                tmp / 'out',
                sample_map=str(sample_map),
                min_quality=10,
                window=4,
                min_length=8,
                min_overlap=8,
                min_overlap_identity=0.90,
            )

            assembled = dict(read_fasta(out['assembled_fasta']))
            self.assertEqual(assembled['IsoBest'], 'CCCCCCCCCC')
            report = Path(out['assembly_tsv']).read_text()
            self.assertIn('IsoBest\tbest_read\t2\t1\t10', report)
            self.assertIn('\tbest_read\thigh_quality\t', report)

    def test_paper_trail_filters_final_outputs_below_minimum_length(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            read = tmp / 'IsoShort_27F.ab1'
            _write_minimal_ab1(read, 'ACGTACGTACGT', [35] * 12)

            out = paper_trail_pipeline.run_paper_trail(
                [str(read)],
                tmp / 'out',
                min_quality=20,
                window=4,
                min_length=20,
                min_read_length=8,
            )

            self.assertEqual(dict(read_fasta(out['assembled_fasta'])), {})
            self.assertEqual(dict(read_fasta(out['failed_final_fasta'])), {'IsoShort': 'ACGTACGTACGT'})
            self.assertEqual(dict(read_fasta(out['failed_read_fasta'])), {})
            failed_manifest = Path(out['failed_manifest_tsv']).read_text()
            self.assertIn('final\tIsoShort\tIsoShort\t\tFAIL_QC\tRESEQUENCE\tfiltered_output_length_lt_20', failed_manifest)
            report = Path(out['assembly_tsv']).read_text()
            self.assertIn('IsoShort\tfiltered_output_length_lt_20\t1\t1\t12', report)
            self.assertIn('\tno\toutput_length_lt_20;read_level_warnings\t', report)

    def test_paper_trail_masks_internal_low_quality_and_recommends_resequence(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            read = tmp / 'IsoBad_27F.ab1'
            seq = 'A' * 200 + 'C' * 12 + 'G' * 200
            qual = [35] * 200 + [10] * 12 + [35] * 200
            _write_minimal_ab1(read, seq, qual)

            out = paper_trail_pipeline.run_paper_trail(
                [str(read)],
                tmp / 'out',
                min_quality=20,
                mask_quality=20,
                min_length=300,
                min_read_length=300,
                max_internal_low_quality_run=5,
            )

            self.assertEqual(dict(read_fasta(out['assembled_fasta'])), {})
            self.assertIn('IsoBad|read=IsoBad_27F', '\n'.join(dict(read_fasta(out['failed_read_fasta'])).keys()))
            failed_manifest = Path(out['failed_manifest_tsv']).read_text()
            self.assertIn('read\tIsoBad\tIsoBad_27F', failed_manifest)
            self.assertIn('internal_low_quality_run_gt_5', failed_manifest)
            read_qc = Path(out['read_qc_tsv']).read_text()
            self.assertIn('internal_low_quality_run_gt_5', read_qc)
            self.assertIn('MaskedBases', read_qc)
            recommendations = Path(out['recommendations_tsv']).read_text()
            self.assertIn('IsoBad\tRESEQUENCE\tFAIL_QC', recommendations)
            self.assertIn('failed_no_reads', recommendations)

    def test_paper_trail_accepts_sequencing_id_mapping_with_prefixed_ab1_names(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            reads = tmp / 'All_AB1'
            reads.mkdir()
            fwd = reads / 'KKX994_70439947_70439947.ab1'
            rev = reads / 'KKY011_70440110_70440110.ab1'
            _write_minimal_ab1(fwd, 'AAAACCCCGGGGTTTTAAAA', [35] * 20)
            _write_minimal_ab1(rev, reverse_complement('CCCCGGGGTTTTAAAACCCC'), [35] * 20)
            mapping = tmp / 'ab1_mapping.csv'
            mapping.write_text(
                '\ufeffSequencing ID,Isolate Number,Read\n'
                'KKX994-F ,SW_0016,Forward\n'
                'KKY011-R,SW_0016,Reverse\n'
            )

            out = paper_trail_pipeline.run_paper_trail(
                [str(reads)],
                tmp / 'out',
                sample_map=str(mapping),
                min_length=8,
                min_overlap=8,
                min_overlap_identity=0.90,
            )

            assembled = dict(read_fasta(out['assembled_fasta']))
            self.assertEqual(assembled['SW_0016'], 'AAAACCCCGGGGTTTTAAAACCCC')

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

            with mock.patch('branchmanager.pipeline.tree.run_cmd', side_effect=fake_run_cmd):
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

            with mock.patch('branchmanager.pipeline.tree.run_cmd', side_effect=fake_run_cmd):
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

            with mock.patch('branchmanager.pipeline.tree.run_cmd', side_effect=fake_run_cmd), \
                 mock.patch('branchmanager.pipeline.tree._run_mafft_addfragments', side_effect=fake_addfragments), \
                 mock.patch('branchmanager.pipeline.tree._run_fasttree', return_value=True):
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

    def test_incremental_tree_matching_preserves_complete_ids(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            alignment = tmp / 'current_alignment.fasta'
            alignment.write_text('>Partner isolate_1\nACGT\n')
            user = tmp / 'user.fasta'
            user.write_text(
                '>Partner isolate_1\nACGT\n'
                '>Partner isolate_2\nACGA\n'
                '>prefix|Partner isolate_1\nACGG\n'
            )

            new_records = _new_sequences_only(user, alignment, None, tmp)

            self.assertEqual(
                [header for header, _sequence in new_records],
                ['Partner isolate_2', 'prefix|Partner isolate_1'],
            )

    def test_initialise_or_update_tree_seeds_previous_review_alignment_and_uses_mafft_add_for_full_length_sequences(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            previous_review = tmp / 'previous_review'
            previous_review.mkdir()
            (previous_review / 'current_alignment.fasta').write_text('>BASE\n' + ('A' * 1500) + '\n')
            (previous_review / 'current_tree.nwk').write_text('(BASE:0.1);\n')
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
                self.assertEqual(Path(backbone_aln).read_text(), (previous_review / 'current_alignment.fasta').read_text())
                Path(output_fasta).write_text('>BASE\n' + ('A' * 1500) + '\n>Q1\n' + ('A' * 1400) + '\n')
                return True

            with mock.patch('branchmanager.pipeline.tree.run_cmd', side_effect=fake_run_cmd), \
                 mock.patch('branchmanager.pipeline.tree._run_mafft_full', return_value=False) as full_mock, \
                 mock.patch('branchmanager.pipeline.tree._run_mafft_add', side_effect=fake_add) as add_mock, \
                 mock.patch('branchmanager.pipeline.tree._run_mafft_addfragments', return_value=False) as addfrag_mock, \
                 mock.patch('branchmanager.pipeline.tree._run_fasttree', return_value=True):
                initialise_or_update_tree(
                    ref_fasta=str(ref),
                    user_fasta=str(user),
                    outdir=str(outdir),
                    db=None,
                    db_dataset=None,
                    threads=1,
                    previous_review=str(previous_review),
                )

            self.assertFalse(full_mock.called)
            self.assertTrue(add_mock.called)
            self.assertFalse(addfrag_mock.called)
            self.assertTrue((outdir / 'current_alignment.fasta').exists())

    def test_initialise_or_update_tree_reuses_organised_backbone_in_same_outdir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            outdir = tmp / 'run_out'
            organised_tree = outdir / 'tree'
            organised_tree.mkdir(parents=True)
            prior_alignment = '>BASE\n' + ('A' * 1500) + '\n'
            (organised_tree / 'current_alignment.fasta').write_text(prior_alignment)
            (organised_tree / 'current_tree.nwk').write_text('(BASE:0.1);\n')
            user = tmp / 'user.fasta'
            user.write_text('>Q1\n' + ('A' * 1400) + '\n')
            ref = tmp / 'ref.fasta'
            ref.write_text('>REF1\n' + ('A' * 1500) + '\n')

            def fake_run_cmd(cmd):
                match = re.search(r'--blast6out\s+(\S+)', cmd)
                self.assertIsNotNone(match)
                Path(match.group(1)).write_text(
                    'Q1\tREF1\t99.0\t1400\t0\t0\t1\t1400\t1\t1500\t1\t0\n'
                )

            def fake_add(new_fasta, backbone_aln, output_fasta, threads=4):
                self.assertEqual(Path(backbone_aln).read_text(), prior_alignment)
                Path(output_fasta).write_text(
                    prior_alignment + '>Q1\n' + ('A' * 1400) + '\n'
                )
                return True

            with mock.patch('branchmanager.pipeline.tree.run_cmd', side_effect=fake_run_cmd), \
                 mock.patch('branchmanager.pipeline.tree._run_mafft_full', return_value=False) as full_mock, \
                 mock.patch('branchmanager.pipeline.tree._run_mafft_add', side_effect=fake_add) as add_mock, \
                 mock.patch('branchmanager.pipeline.tree._run_mafft_addfragments', return_value=False), \
                 mock.patch('branchmanager.pipeline.tree._run_fasttree', return_value=True):
                initialise_or_update_tree(
                    ref_fasta=str(ref),
                    user_fasta=str(user),
                    outdir=str(outdir),
                    db=None,
                    db_dataset=None,
                    threads=1,
                )

            self.assertFalse(full_mock.called)
            self.assertTrue(add_mock.called)
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

            with mock.patch('branchmanager.pipeline.classify.run_cmd', side_effect=fake_run_cmd):
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

            with mock.patch('branchmanager.pipeline.novelty.run_cmd', side_effect=fake_run_cmd):
                novelty_pipeline.run_novelty(str(query), str(ref), str(tmp), target_fasta=str(ref))

            self.assertIn('--strand both', seen['cmd'])

    def test_write_dataset_membership_strip_writes_expected_rows(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir) / 'membership.itol'
            write_dataset_membership_strip(
                str(out),
                ['A1', 'A2', 'A3'],
                {'A1': 'previous_review', 'A2': 'run1', 'A3': ''},
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
            with mock.patch('branchmanager.pipeline.collapse.shutil.which', return_value=None):
                artefacts = collapse_fasta_within_taxa(
                    taxa_groups,
                    tmpdir,
                    'collapsed.fasta',
                    'collapsed.tsv',
                    'members.tsv',
                    threshold=99.0,
                    threads=1,
                    log_prefix='[TEST COLLAPSE]',
                )
            self.assertTrue(artefacts.collapsed_path.exists())
            fasta_text = artefacts.collapsed_path.read_text()
            self.assertIn('>S1', fasta_text)
            self.assertIn('>S2', fasta_text)
            self.assertIn('>S3', fasta_text)

    def test_internal_node_labelling_does_not_introduce_double_colons(self):
        newick = '((A:0.1,B:0.2)0.95:0.3,C:0.4);'
        labelled = _label_internal_nodes(newick)
        self.assertIn('NODE', labelled)
        self.assertNotIn('::', labelled)

    def test_repair_internal_node_label_delimiters(self):
        broken = '((A:0.1,B:0.2)NODE0001::0.3,C:0.4);'
        repaired = _repair_internal_node_label_delimiters(broken)
        self.assertEqual(repaired, '((A:0.1,B:0.2)NODE0001:0.3,C:0.4);')

    def test_generate_itol_colours_repairs_malformed_tree_labels(self):
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
            generate_itol_colours(str(taxonomy_tsv), str(tmp), tree_file=str(tree_path))
            repaired = tree_path.read_text()
            self.assertIn('NODE0001:0.3', repaired)
            self.assertNotIn('NODE0001::0.3', repaired)
            self.assertTrue((tmp / 'itol_phylum_colours.itol').exists())
            self.assertTrue((tmp / 'itol_family_colours.itol').exists())
            self.assertTrue((tmp / 'itol_genus_colours.itol').exists())
            self.assertFalse((tmp / 'itol_phylum_colors.itol').exists())

    def test_tree_warning_helpers_report_partial_sequences_and_alignment_issues(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            user_fasta = tmp / 'user.fasta'
            user_fasta.write_text('>A\nACGT\n>B\nACGTNNNN\n')
            with mock.patch('branchmanager.pipeline.tree.load_anchor_sequences', return_value=[]):
                warnings = collect_tree_build_warnings(str(user_fasta), anchor_file='missing_anchor_file.fasta')
            cats = {w['category'] for w in warnings}
            self.assertIn('LOW_SEQUENCE_COUNT', cats)
            self.assertIn('PARTIAL_16S_SEQUENCES', cats)
            self.assertIn('MISSING_REFERENCE_ANCHORS', cats)

            aln = tmp / 'alignment.fasta'
            aln.write_text('>A\n----\n>B\nAC--\n')
            aln_warnings = summarise_alignment_quality(str(aln))
            aln_cats = {w['category'] for w in aln_warnings}
            self.assertIn('ALL_GAP_ALIGNMENT_ROWS', aln_cats)

    def test_run_qc_writes_per_sequence_rejection_reasons(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            fasta = tmp / 'input.fasta'
            fasta.write_text(
                '>keep\nACGTACGT\n'
                '>short\nACG\n'
                '>ambiguous\nACGTNNN\n'
                '>both\nNN\n'
            )

            qc_pipeline.run_qc(str(fasta), str(tmp), min_len=4, max_n=1)

            stats = (tmp / 'qc.stats').read_text()
            rejections = (tmp / 'qc_rejections.tsv').read_text()
            self.assertIn('rejected_total\t3', stats)
            self.assertIn('rejection_details\t', stats)
            self.assertIn('short\t3\t0\ttoo_short\t4\t1', rejections)
            self.assertIn('ambiguous\t7\t3\ttoo_many_n\t4\t1', rejections)
            self.assertIn('both\t2\t2\ttoo_short;too_many_n\t4\t1', rejections)

    def test_write_output_explanations_writes_detailed_output_guide(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            assessment = tmp / 'assessment'
            assessment.mkdir()
            (assessment / 'sequence_assessment.tsv').write_text('ID\tNoveltyScore\nS01\t42.0\n')
            (assessment / 'qc.stats').write_text(
                'total_input\t4\n'
                'kept\t1\n'
                'rejected_total\t3\n'
                'rejected_too_short\t2\n'
                'rejected_too_many_n\t2\n'
                'min_len\t4\n'
                'max_n\t1\n'
            )
            (assessment / 'qc_rejections.tsv').write_text(
                'ID\tLength\tNCount\tReasons\tMinLength\tMaxN\nshort\t3\t0\ttoo_short\t4\t1\n'
            )

            _write_output_explanations(str(tmp))

            guide = (tmp / 'OUTPUT_GUIDE.md').read_text()
            self.assertIn('QC Filtering', guide)
            self.assertIn('Rejected as too short: `2`', guide)
            self.assertIn('NoveltyScore', guide)
            self.assertIn('assessment/qc_rejections.tsv', guide)

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
            input_fasta.write_text('>Q1\nACGT\n')
            ref_fasta.write_text('>R1 d__Bacteria;p__Firmicutes\nACGT\n')

            def fake_run_cmd(cmd):
                parts = cmd.split()
                density = Path(parts[parts.index('--blast6out') + 1])
                density.write_text('Q1\tR1\t99.2\nQ1\tR2\t97.4\nQ1\tR3\t95.1\n')

            with mock.patch('branchmanager.pipeline.novelty.run_cmd', side_effect=fake_run_cmd):
                out = build_reference_novelty_metrics(str(input_fasta), str(ref_fasta), str(tmp))

            text = Path(out).read_text()
            self.assertIn('NoveltyScore', text)
            self.assertIn('MatchesGE99', text)
            self.assertIn('Q1\t99.20\tR1\tFalse\t1\t2\t3', text)

    def test_build_reference_novelty_metrics_reports_baseline_and_project_scores(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            input_fasta = tmp / 'query.fasta'
            ref_fasta = tmp / 'ref.fasta'
            input_fasta.write_text('>Q1\nACGT\n')
            ref_fasta.write_text('>R1 d__Bacteria;p__Firmicutes\nACGT\n')

            db = Database(str(tmp / 'test.sqlite'))
            db.initialise()
            db.insert_sequences([('HungateHit', 'ACGT'), ('PartnerHit', 'ACGT')], dataset='Hungate')
            with db.connect() as conn:
                conn.execute(
                    'UPDATE sequences SET dataset = ? WHERE id = ?',
                    ('PartnerB', 'PartnerHit'),
                )

            def fake_run_cmd(cmd):
                parts = cmd.split()
                out = Path(parts[parts.index('--blast6out') + 1])
                if 'baseline_nearest' in out.name:
                    out.write_text('Q1\tHungateHit\t96.2\n')
                elif 'baseline_density' in out.name:
                    out.write_text('Q1\tHungateHit\t96.2\n')
                elif 'project_nearest' in out.name:
                    out.write_text('Q1\tPartnerHit\t99.4\n')
                elif 'project_density' in out.name:
                    out.write_text('Q1\tPartnerHit\t99.4\n')
                elif 'reference_nearest' in out.name:
                    out.write_text('Q1\tGTDB_REF\t98.8\n')
                elif 'reference_density' in out.name:
                    out.write_text('Q1\tGTDB_REF\t98.8\n')
                else:
                    out.write_text('')

            with mock.patch('branchmanager.pipeline.novelty.run_cmd', side_effect=fake_run_cmd):
                out = build_reference_novelty_metrics(
                    str(input_fasta),
                    str(ref_fasta),
                    str(tmp),
                    db=db,
                    run_dataset='CurrentBatch',
                    baseline_datasets=['Hungate'],
                )

            text = Path(out).read_text()
            self.assertIn('BaselineNoveltyScore', text)
            self.assertIn('ProjectNoveltyScore', text)
            self.assertIn('ReferenceNoveltyScore', text)
            self.assertIn('Q1\t96.20\tHungateHit\tTrue', text)
            self.assertIn('\t99.40\tPartnerHit\tFalse\t1\t1\t1\t', text)
            self.assertIn('\t98.80\tGTDB_REF\tFalse\t0\t1\t1\t', text)

    def test_build_reference_novelty_metrics_reports_selected_wgs_clade_context(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            input_fasta = tmp / 'query.fasta'
            ref_fasta = tmp / 'ref.fasta'
            input_fasta.write_text('>Q1\nACGT\n')
            ref_fasta.write_text('>R1 d__Bacteria;p__Firmicutes\nACGT\n')

            db = Database(str(tmp / 'test.sqlite'))
            db.initialise()
            db.insert_sequences([('Q1', 'ACGT')], dataset='CurrentBatch')
            db.insert_sequences([('SelectedPrev', 'ACGT')], dataset='OlderBatch')
            db.upsert_sequencing_metadata([
                {
                    'id': 'Q1',
                    'partner_id': 'PartnerQ1',
                    'dataset': 'CurrentBatch',
                    'selected_for_wgs': False,
                },
                {
                    'id': 'SelectedPrev',
                    'partner_id': 'PartnerPrev',
                    'dataset': 'OlderBatch',
                    'selected_for_wgs': True,
                },
            ])

            def fake_run_cmd(cmd):
                parts = cmd.split()
                out = Path(parts[parts.index('--blast6out') + 1])
                if 'selected_for_wgs' in out.name:
                    out.write_text('Q1\tSelectedPrev\t99.5\n')
                elif 'project_density' in out.name:
                    out.write_text('Q1\tSelectedPrev\t99.5\n')
                elif 'reference_nearest' in out.name:
                    out.write_text('Q1\tR1\t99.0\n')
                elif 'reference_density' in out.name:
                    out.write_text('Q1\tR1\t99.0\n')
                else:
                    out.write_text('')

            with mock.patch('branchmanager.pipeline.novelty.run_cmd', side_effect=fake_run_cmd):
                out = build_reference_novelty_metrics(
                    str(input_fasta),
                    str(ref_fasta),
                    str(tmp),
                    db=db,
                    run_dataset='CurrentBatch',
                )

            text = Path(out).read_text()
            self.assertIn('GenomeAlreadySequenced', text)
            self.assertIn('RelatedGenomeCladeGE97', text)
            self.assertIn('PartnerQ1\tFalse\tFalse\tSelectedPrev\t99.50\t1\t1\t1\tTrue', text)
            self.assertIn('\t3\t3\tbaseline_genomes_and_project_ledger', text)

    def test_build_sequence_assessment_rows_combines_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            db = Database(os.path.join(tmpdir, 'test.sqlite'))
            db.initialise()
            db.insert_sequences([('REF1', 'ACGT')], dataset='Hungate')
            with db.connect() as conn:
                conn.execute(
                    'INSERT INTO taxonomy (id, dataset, taxonomy, confidence) VALUES (?, ?, ?, ?)',
                    ('REF1', 'Hungate', 'd__Bacteria; p__Firmicutes', 1.0),
                )
            class_tsv = tmp / 'taxonomy.tsv'
            class_tsv.write_text('ID\tBestHit\tIdentity\tTaxon\tConfidence\norigA\tREF1\t98.0\td__Bacteria; p__Firmicutes\t0.9\n')
            novelty_metrics = tmp / 'novelty_metrics.tsv'
            novelty_metrics.write_text(
                'ID\tBaselineNearestIdentity\tBaselineNearestHit\tBaselineNovel\tBaselineMatchesGE99\tBaselineMatchesGE97\tBaselineMatchesGE95\tBaselineNoveltyScore\tBaselineCrowding\tBaselineSequencingPriority\tBaselineDensitySource\n'
                'S01\t96.50\tREF1\tTrue\t1\t2\t3\t25.00\tsparse\tHIGH\tbaseline:Hungate\n'
            )
            warning_rows = [{'id': 'S01', 'flags': 'LOW_NEAREST_IDENTITY'}]
            rows = build_sequence_assessment_rows(['S01'], str(class_tsv), str(novelty_metrics), warning_rows, {'origA': 'S01'}, db)
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]['id'], 'S01')
            self.assertEqual(rows[0]['classification_hit'], 'REF1')
            self.assertEqual(rows[0]['nearest_hit_dataset'], 'Hungate')
            self.assertEqual(rows[0]['nearest_hit_taxonomy'], 'd__Bacteria; p__Firmicutes')
            self.assertEqual(rows[0]['sequencing_priority'], 'HIGH')
            self.assertEqual(rows[0]['placement_flags'], 'LOW_NEAREST_IDENTITY')

    def test_build_sequence_assessment_rows_reads_project_novelty_columns(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            db = Database(os.path.join(tmpdir, 'test.sqlite'))
            db.initialise()
            db.insert_sequences([('H1', 'ACGT')], dataset='Hungate')
            db.insert_sequences([('P1', 'ACGT')], dataset='PartnerB')
            with db.connect() as conn:
                conn.execute(
                    'INSERT INTO taxonomy (id, dataset, taxonomy, confidence) VALUES (?, ?, ?, ?)',
                    ('H1', 'Hungate', 'd__Bacteria; p__Bacillota; g__Hungate', 1.0),
                )
                conn.execute(
                    'INSERT INTO taxonomy (id, dataset, taxonomy, confidence) VALUES (?, ?, ?, ?)',
                    ('P1', 'PartnerB', 'd__Bacteria; p__Bacillota; g__Partner', 1.0),
                )
            class_tsv = tmp / 'taxonomy.tsv'
            class_tsv.write_text('ID\tBestHit\tIdentity\tTaxon\tConfidence\nS01\tREF1\t98.0\td__Bacteria; p__Bacillota\t0.9\n')
            novelty_metrics = tmp / 'novelty_metrics.tsv'
            novelty_metrics.write_text(
                'ID\tNearestIdentity\tNearestHit\tNovel\tMatchesGE99\tMatchesGE97\tMatchesGE95\tNoveltyScore\tCrowding\tSequencingPriority\tDensitySource\t'
                'BaselineNearestIdentity\tBaselineNearestHit\tBaselineNovel\tBaselineMatchesGE99\tBaselineMatchesGE97\tBaselineMatchesGE95\tBaselineNoveltyScore\tBaselineCrowding\tBaselineSequencingPriority\tBaselineDensitySource\t'
                'ProjectNearestIdentity\tProjectNearestHit\tProjectNovel\tProjectMatchesGE99\tProjectMatchesGE97\tProjectMatchesGE95\tProjectNoveltyScore\tProjectCrowding\tProjectSequencingPriority\tProjectDensitySource\t'
                'ReferenceNearestIdentity\tReferenceNearestHit\tReferenceNovel\tReferenceMatchesGE99\tReferenceMatchesGE97\tReferenceMatchesGE95\tReferenceNoveltyScore\tReferenceCrowding\tReferenceSequencingPriority\tReferenceDensitySource\n'
                'S01\t96.20\tH1\tTrue\t0\t0\t1\t55.00\tisolated\tHIGH\tbaseline:Hungate\t'
                '96.20\tH1\tTrue\t0\t0\t1\t55.00\tisolated\tHIGH\tbaseline:Hungate\t'
                '99.40\tP1\tFalse\t1\t1\t2\t26.80\tisolated\tLOW\tproject_collection\t'
                '98.80\tREF1\tFalse\t0\t1\t1\t33.60\tisolated\tLOW\treference_fasta\n'
            )

            rows = build_sequence_assessment_rows(['S01'], str(class_tsv), str(novelty_metrics), [], {}, db)

            self.assertEqual(rows[0]['nearest_hit'], 'H1')
            self.assertEqual(rows[0]['nearest_hit_dataset'], 'Hungate')
            self.assertEqual(rows[0]['project_nearest_hit'], 'P1')
            self.assertEqual(rows[0]['project_nearest_hit_dataset'], 'PartnerB')
            self.assertEqual(rows[0]['project_novelty_score'], '26.80')
            self.assertEqual(rows[0]['reference_nearest_hit'], 'REF1')
            self.assertEqual(rows[0]['reference_nearest_hit_taxonomy'], 'd__Bacteria; p__Bacillota')
            self.assertEqual(rows[0]['reference_novelty_score'], '33.60')

    def test_write_selection_summary_tsv_keeps_board_facing_fields(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir) / 'selection_summary.tsv'
            rows = [
                {
                    'id': 'Iso001',
                    'partner_id': 'QUB',
                    'taxonomy': 'd__Bacteria; p__Bacillota; g__Novel',
                    'classification_identity': '98.4',
                    'nearest_identity': '95.50',
                    'nearest_hit': 'HungateA',
                    'nearest_hit_taxonomy': 'd__Bacteria; p__Bacillota; g__Known',
                    'matches_ge_97': '0',
                    'density_source': 'baseline:Hungate',
                    'project_nearest_identity': '96.10',
                    'project_matches_ge_97': '0',
                    'project_density_source': 'project_collection',
                    'project_novelty_score': '55.0',
                    'reference_nearest_identity': '98.80',
                    'reference_nearest_hit': 'GTDB_REF',
                    'reference_density_source': 'reference_fasta',
                    'sequencing_priority': 'HIGH',
                    'selected_for_genome_sequencing': 'False',
                    'genome_committed_count_same_species': '0',
                    'genome_available_count_same_species': '0',
                    'genome_selected_count_same_species': '0',
                    'pangenome_target': '3',
                    'pangenome_gap': '3',
                    'sequencing_set_id': 'SET_A',
                    'sequencing_set_role': 'PRIMARY',
                    'sequencing_set_rank': '1',
                    'mwl_match': 'Yes',
                    'mwl_matched_rank': 'genus',
                    'mwl_matched_taxon': 'g__Novel',
                    'mwl_score': '88.00',
                    'in_tree': 'Yes',
                    'cluster_representative': 'self',
                    'cluster_size': '1',
                    'placement_flags': '',
                    'local_neighbourhood_figure': 'neighbourhoods/clade_001.png',
                },
                {
                    'id': 'Iso002',
                    'partner_id': 'UoG',
                    'taxonomy': 'd__Bacteria; p__Bacillota',
                    'classification_identity': '99.0',
                    'nearest_identity': '96.80',
                    'nearest_hit': 'HungateB',
                    'nearest_hit_taxonomy': 'd__Bacteria; p__Bacillota',
                    'matches_ge_97': '1',
                    'density_source': 'baseline:Hungate',
                    'project_nearest_identity': '99.20',
                    'project_matches_ge_97': '5',
                    'project_density_source': 'project_collection',
                    'project_novelty_score': '10.0',
                    'reference_nearest_identity': '99.50',
                    'reference_nearest_hit': 'GTDB_REF2',
                    'reference_density_source': 'reference_fasta',
                    'sequencing_priority': 'HIGH',
                    'selected_for_genome_sequencing': 'False',
                    'nearest_genome_hit': 'IsoPrev',
                    'nearest_genome_identity': '99.20',
                    'genome_committed_count_same_species': '3',
                    'genome_available_count_same_species': '2',
                    'genome_selected_count_same_species': '1',
                    'pangenome_target': '3',
                    'pangenome_gap': '0',
                    'sequencing_set_id': 'SET_B',
                    'sequencing_set_role': 'TARGET_MET',
                    'mwl_match': 'No',
                    'in_tree': 'Yes',
                    'cluster_representative': 'self',
                    'cluster_size': '1',
                    'placement_flags': '',
                },
            ]

            path = write_selection_summary_tsv(out, rows)
            text = Path(path).read_text()

            self.assertIn(
                'SequenceID\tPartnerID\tSelectedForGenomeSequencing\t'
                'GenomeAlreadySequenced\tRecommendation\tSequencingSetID',
                text,
            )
            self.assertIn('Iso001\tQUB\tFalse\tNA\tPRIORITISE - SET PRIMARY\tSET_A\tPRIMARY\t1\tHIGH', text)
            self.assertIn('MWL target match (g__Novel)', text)
            self.assertIn('neighbourhoods/clade_001.png', text)
            self.assertIn('Iso002\tUoG\tFalse\tNA\tLOWER PRIORITY - TARGET MET\tSET_B\tTARGET_MET', text)
            self.assertIn('assessment-species pangenome target already met', text)

    def test_selection_decision_does_not_treat_related_available_genome_as_duplicate(self):
        decision = build_selection_decision({
            'taxonomy': 'd__Bacteria; p__Bacillota',
            'classification_identity': '98.0',
            'nearest_identity': '96.0',
            'density_source': 'baseline:Hungate',
            'project_nearest_identity': '96.5',
            'project_matches_ge_97': '0',
            'project_density_source': 'project_collection',
            'reference_nearest_identity': '97.5',
            'reference_density_source': 'reference_fasta',
            'selected_for_genome_sequencing': 'False',
            'nearest_genome_identity': '97.8',
            'mwl_match': 'No',
            'placement_flags': '',
        })
        self.assertEqual(decision['decision'], 'STRONG CANDIDATE')
        self.assertEqual(decision['genome_coverage'], 'RELATED AVAILABLE GENOME - 97.80%')

    def test_sequencing_sets_propose_three_primaries_and_one_backup(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            tree_path = tmp / 'tree.nwk'
            tree_path.write_text(
                '((Q1:0.01,Q2:0.02):0.01,(Q3:0.03,(Q4:0.04,Q5:0.05):0.01):0.02);'
            )
            rows = []
            for index in range(1, 6):
                rows.append({
                    'id': f'Q{index}',
                    'partner_id': 'QUB',
                    'taxonomy': 'd__Bacteria; p__Bacillota; g__Novel; s__Novel species',
                    'classification_identity': '99.0',
                    'classification_confidence': '0.99',
                    'nearest_identity': str(94.0 + index / 10),
                    'density_source': 'baseline:Hungate',
                    'project_nearest_identity': str(96.0 + index / 10),
                    'project_density_source': 'project_collection',
                    'reference_nearest_identity': '97.0',
                    'reference_density_source': 'reference_fasta',
                    'selected_for_genome_sequencing': 'False',
                    'genome_available_count_same_species': '0',
                    'genome_selected_count_same_species': '0',
                    'genome_committed_count_same_species': '0',
                    'pangenome_target': '3',
                    'pangenome_gap': '3',
                    'in_tree': 'Yes',
                    'cluster_representative': 'self',
                    'placement_flags': '',
                })

            output = selection_sets_pipeline.build_sequencing_sets(
                rows,
                tmp / 'sequencing_sets.tsv',
                tree_path=tree_path,
                pangenome_target=3,
                candidate_set_size=4,
            )

            roles = [row['sequencing_set_role'] for row in rows]
            self.assertEqual(roles.count('PRIMARY'), 3)
            self.assertEqual(roles.count('BACKUP'), 1)
            self.assertEqual(roles.count('ALTERNATE'), 1)
            text = Path(output).read_text()
            self.assertIn('CommittedGenomeCount\tCommittedGenomeIDs\tPangenomeTarget\tPangenomeGap', text)
            self.assertIn('\tBACKUP\t4\tHIGH\t', text)

    def test_sequencing_sets_count_every_baseline_as_genome_available(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            db = Database(str(tmp / 'project.sqlite'))
            db.initialise()
            db.upsert_dataset_role('Hungate', 'baseline', genomes_available=True)
            db.insert_sequences([('H1', 'ACGT'), ('H2', 'ACGT')], dataset='Hungate')
            db.insert_sequences([('P1', 'ACGT'), ('Q1', 'ACGT')], dataset='Partners')
            db.upsert_dataset_role('Partners', 'candidate', genomes_available=False)
            taxonomy = 'd__Bacteria; p__Bacillota; g__Covered; s__Covered species'
            db.insert_taxonomy([(sid, taxonomy, 1.0, dataset) for sid, dataset in (
                ('H1', 'Hungate'), ('H2', 'Hungate'), ('P1', 'Partners'), ('Q1', 'Partners'),
            )])
            db.upsert_sequencing_metadata([{
                'id': 'P1', 'partner_id': 'QUB', 'dataset': 'Partners', 'selected_for_wgs': True,
            }])
            rows = [{
                'id': 'Q1', 'partner_id': 'QUB', 'taxonomy': taxonomy,
                'classification_identity': '99.0', 'selected_for_genome_sequencing': 'False',
                'nearest_identity': '99.0', 'density_source': 'baseline:Hungate',
                'project_density_source': 'project_collection', 'reference_density_source': 'reference_fasta',
                'genome_available_count_same_species': '2',
                'genome_selected_count_same_species': '1',
                'genome_committed_count_same_species': '3',
                'pangenome_target': '3', 'pangenome_gap': '0',
                'placement_flags': '',
            }]

            output = selection_sets_pipeline.build_sequencing_sets(
                rows, tmp / 'sequencing_sets.tsv', db=db,
            )

            self.assertEqual(rows[0]['sequencing_set_role'], 'TARGET_MET')
            text = Path(output).read_text()
            self.assertIn('H1;H2;P1', text)
            self.assertIn('\t2\t1\t0\t3\t', text)

    def test_already_sequenced_is_a_state_not_a_new_recommendation(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            row = {
                'id': 'P1', 'partner_id': 'QUB',
                'taxonomy': 'd__Bacteria; g__Covered; s__Covered species',
                'classification_identity': '99.0',
                'selected_for_genome_sequencing': 'True',
                'already_sequenced': 'True',
                'genome_available_count_same_species': '0',
                'genome_selected_count_same_species': '1',
                'genome_committed_count_same_species': '1',
                'pangenome_target': '3', 'pangenome_gap': '2',
                'density_source': 'baseline:Hungate',
                'project_density_source': 'project_collection',
                'reference_density_source': 'reference_fasta',
                'placement_flags': '',
            }

            output = selection_sets_pipeline.build_sequencing_sets(
                [row], Path(tmpdir) / 'sequencing_sets.tsv',
            )

            self.assertEqual(row['sequencing_set_role'], 'SEQUENCED')
            self.assertEqual(row['sequencing_set_rank'], 'NA')
            self.assertIn('\tSEQUENCED\tNA\t', Path(output).read_text())
            self.assertEqual(build_selection_decision(row)['decision'], 'ALREADY SEQUENCED')

    def test_selected_pending_is_committed_but_not_already_sequenced(self):
        row = {
            'id': 'Pending', 'partner_id': 'QUB',
            'taxonomy': 'd__Bacteria; g__Covered; s__Covered species',
            'classification_identity': '99.0', 'classification_confidence': '0.99',
            'selected_for_genome_sequencing': 'True', 'already_sequenced': 'False',
            'genome_available_count_same_species': '0',
            'genome_selected_count_same_species': '0',
            'genome_pending_count_same_species': '1',
            'genome_committed_count_same_species': '1',
            'pangenome_target': '3', 'pangenome_gap': '2',
            'density_source': 'baseline:Hungate',
            'project_density_source': 'project_collection',
            'reference_density_source': 'reference_fasta',
            'placement_flags': '',
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            selection_sets_pipeline.build_sequencing_sets([row], Path(tmpdir) / 'pending_sets.tsv')
            decision = build_selection_decision(row)

        self.assertEqual(row['sequencing_set_role'], 'COMMITTED')
        self.assertEqual(decision['decision'], 'ALREADY SELECTED - GENOME PENDING')
        self.assertIn('GENOME PENDING', decision['genome_coverage'])

    def test_database_keeps_latest_assessment_and_selection_round_audit(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            db = Database(str(Path(tmpdir) / 'project.sqlite'))
            db.initialise()
            db.save_assessment_snapshot('round_1', [{'id': 'A', 'nearest_identity': '99.0'}])
            db.save_assessment_snapshot('round_2', [{'id': 'A', 'nearest_identity': '97.0'}])

            latest = db.get_latest_assessment_rows()
            self.assertEqual(latest['A']['nearest_identity'], '97.0')
            self.assertEqual(latest['A']['_snapshot_id'], 'round_2')

            recommendations = [{
                'sequence_id': 'A', 'role': 'PRIMARY', 'round_rank': 1,
                'source_snapshot': 'round_2',
            }]
            self.assertEqual(
                db.save_selection_round('quarterly_review_1', recommendations, parameters={'GenomeBudget': 1}),
                1,
            )
            with db.connect() as conn:
                stored = conn.execute(
                    'SELECT role, round_rank FROM selection_round_members WHERE round_id = ?',
                    ('quarterly_review_1',),
                ).fetchone()
            self.assertEqual(stored, ('PRIMARY', 1))

    def test_quarterly_review_fills_gaps_then_balances_post_target_diversity(self):
        def row(sequence_id, species, available, nearest_genome, *, sequenced=False):
            return {
                'id': sequence_id,
                'partner_id': 'QUB',
                'taxonomy': f'd__Bacteria; g__Test; s__{species}',
                'classification_identity': '99.5',
                'classification_confidence': '0.99',
                'nearest_identity': '98.0',
                'project_nearest_identity': '98.0',
                'reference_nearest_identity': '98.0',
                'nearest_genome_identity': str(nearest_genome),
                'selected_for_genome_sequencing': str(sequenced),
                'already_sequenced': str(sequenced),
                'genome_committed_count_same_species': str(available),
                'pangenome_target': '3',
                'placement_flags': '',
            }

        rows = [
            row('BetaA', 'Beta species', 2, 99.5),
            row('BetaB', 'Beta species', 2, 99.6),
            row('AlphaA', 'Alpha species', 3, 97.0),
            row('AlphaB', 'Alpha species', 3, 99.8),
            row('GammaGenome', 'Gamma species', 3, 100.0, sequenced=True),
        ]
        recommendations = quarterly_review_pipeline.build_quarterly_review(
            rows,
            genome_budget=2,
            backups_per_primary=1,
            pangenome_target=3,
        )
        by_id = {item['sequence_id']: item for item in recommendations}

        self.assertEqual(by_id['BetaA']['role'], 'PRIMARY')
        self.assertEqual(by_id['BetaA']['round_rank'], 1)
        self.assertEqual(by_id['BetaA']['priority_tier'], 'COVERAGE_GAP')
        self.assertEqual(by_id['AlphaA']['role'], 'PRIMARY')
        self.assertEqual(by_id['AlphaA']['round_rank'], 2)
        self.assertEqual(by_id['AlphaA']['priority_tier'], 'NOVEL_GENOME_NEIGHBOURHOOD')
        self.assertEqual(by_id['BetaB']['role'], 'BACKUP')
        self.assertEqual(by_id['AlphaB']['role'], 'BACKUP')
        self.assertEqual(by_id['GammaGenome']['role'], 'ALREADY_SEQUENCED')

        with tempfile.TemporaryDirectory() as tmpdir:
            outputs = quarterly_review_pipeline.write_quarterly_review_reports(
                tmpdir, 'quarterly_review_test', recommendations, {'GenomeBudget': 2},
            )
            self.assertEqual(Path(outputs['summary']).name, 'quarterly_review_summary.tsv')
            self.assertEqual(Path(outputs['manifest']).name, 'quarterly_review_manifest.tsv')
            self.assertIn('BetaA', Path(outputs['selected']).read_text())
            self.assertIn('GammaGenome', Path(outputs['summary']).read_text())
            self.assertIn('ScientificScope', Path(outputs['manifest']).read_text())

    def test_quarterly_review_refreshes_available_genome_identity_from_current_alignment(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            alignment = Path(tmpdir) / 'current_alignment.fasta'
            alignment.write_text('>Candidate\nACGT\n>Genome\nACGA\n')
            common = {
                'partner_id': 'QUB',
                'taxonomy': 'd__Bacteria; g__Test; s__Test species',
                'classification_identity': '99.0',
                'classification_confidence': '0.99',
                'genome_committed_count_same_species': '1',
                'placement_flags': '',
                'in_tree': 'Yes',
            }
            rows = [
                {
                    **common, 'id': 'Candidate',
                    'selected_for_genome_sequencing': 'False',
                    'nearest_genome_hit': 'OldGenome',
                    'nearest_genome_identity': '99.9',
                },
                {
                    **common, 'id': 'Genome',
                    'selected_for_genome_sequencing': 'True',
                    'already_sequenced': 'True',
                    'nearest_genome_identity': '100.0',
                },
            ]

            recommendations = quarterly_review_pipeline.build_quarterly_review(
                rows, genome_budget=1, backups_per_primary=0,
                alignment_path=alignment,
            )
            candidate = next(row for row in recommendations if row['sequence_id'] == 'Candidate')
            self.assertEqual(candidate['nearest_available_genome'], 'Genome')
            self.assertEqual(candidate['nearest_available_identity'], '75.00')
            self.assertIn('current_alignment:4_comparable_acgt_columns', candidate['nearest_available_identity_source'])

    def test_local_neighbourhood_groups_assessed_sequences_in_one_png(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            tree_path = tmp / 'current_tree.nwk'
            tree_path.write_text(
                '(((Q1:0.01,Q2:0.01)NODE0001:0.01,H1:0.02)NODE0002:0.01,'
                '(P1:0.02,P2:0.02)NODE0003:0.02)NODE0000;'
            )
            db = Database(str(tmp / 'project.sqlite'))
            db.initialise()
            db.insert_sequences([('Q1', 'ACGT'), ('Q2', 'ACGT')], dataset='Current')
            db.insert_sequences([('H1', 'ACGT')], dataset='Hungate')
            db.insert_sequences([('P1', 'ACGT'), ('P2', 'ACGT')], dataset='Prior')
            db.upsert_sequencing_metadata([{
                'id': 'P1',
                'partner_id': 'QUB',
                'dataset': 'Prior',
                'selected_for_wgs': True,
            }])
            rows = [
                {
                    'id': 'Q1', 'cluster_representative': 'self',
                    'taxonomy': 'd__Bacteria; g__Alpha', 'partner_id': 'QUB',
                    'density_source': 'baseline:Hungate',
                    'selected_for_genome_sequencing': 'False',
                },
                {
                    'id': 'Q2', 'cluster_representative': 'self',
                    'taxonomy': 'd__Bacteria; g__Alpha', 'partner_id': 'UoG',
                    'density_source': 'baseline:Hungate',
                    'selected_for_genome_sequencing': 'False',
                },
            ]

            result = neighbourhood_pipeline.generate_local_neighbourhood_visuals(
                tree_path,
                rows,
                db,
                tmp / 'neighbourhoods',
                min_context_leaves=3,
                max_context_leaves=5,
                image_format='png',
            )

            self.assertEqual(len(result['figures']), 1)
            self.assertEqual(rows[0]['local_neighbourhood_figure'], rows[1]['local_neighbourhood_figure'])
            png_path = Path(result['figures'][0])
            self.assertEqual(png_path.suffix, '.png')
            self.assertEqual(png_path.read_bytes()[:8], b'\x89PNG\r\n\x1a\n')
            manifest = Path(result['manifest']).read_text()
            self.assertIn('Q1;Q2', manifest)
            with open(result['manifest']) as handle:
                manifest_rows = list(csv.DictReader(handle, delimiter='\t'))
            self.assertEqual(manifest_rows[0]['TreeLeavesShown'], '3')
            self.assertEqual(manifest_rows[0]['BaselineLeavesShown'], '1')
            self.assertEqual(manifest_rows[0]['SequencedGenomeLeavesShown'], '0')

    def test_local_neighbourhood_forces_nearest_baseline_and_labels_rank_and_pident(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            tree_path = tmp / 'current_tree.nwk'
            tree_path.write_text(
                '((Q1:0.01,Q2:0.01)NODE0001:0.01,'
                '(H1:0.01,(P1:0.01,P2:0.01)NODE0002:0.01)NODE0003:0.01)NODE0000;'
            )
            alignment = tmp / 'current_alignment.fasta'
            alignment.write_text(
                '>Q1\nACGT\n>Q2\nACGA\n>H1\nAGGT\n>P1\nAAAA\n>P2\nTTTT\n'
            )
            db = Database(str(tmp / 'project.sqlite'))
            db.initialise()
            db.insert_sequences([('Q1', 'ACGT'), ('Q2', 'ACGA')], dataset='Current')
            db.insert_sequences([('H1', 'AGGT')], dataset='Hungate')
            db.insert_sequences([('P1', 'AAAA'), ('P2', 'TTTT')], dataset='Prior')
            rows = [
                {
                    'id': 'Q1', 'cluster_representative': 'self',
                    'taxonomy': 'd__Bacteria; g__Alpha', 'partner_id': 'QUB',
                    'density_source': 'baseline:Hungate',
                    'nearest_hit': 'H1', 'nearest_identity': '98.50',
                    'sequencing_set_id': 'BMSET_TEST',
                    'sequencing_set_role': 'PRIMARY', 'sequencing_set_rank': '1',
                    'selected_for_genome_sequencing': 'False',
                },
                {
                    'id': 'Q2', 'cluster_representative': 'self',
                    'taxonomy': 'd__Bacteria; g__Alpha', 'partner_id': 'UoG',
                    'density_source': 'baseline:Hungate',
                    'nearest_hit': 'H1', 'nearest_identity': '97.50',
                    'sequencing_set_id': 'BMSET_TEST',
                    'sequencing_set_role': 'BACKUP', 'sequencing_set_rank': '2',
                    'selected_for_genome_sequencing': 'False',
                },
            ]

            result = neighbourhood_pipeline.generate_local_neighbourhood_visuals(
                tree_path,
                rows,
                db,
                tmp / 'neighbourhoods',
                alignment_path=alignment,
                min_context_leaves=2,
                max_context_leaves=2,
                image_format='png',
            )

            png_path = Path(result['figures'][0])
            self.assertEqual(png_path.read_bytes()[:8], b'\x89PNG\r\n\x1a\n')
            root, _, _ = neighbourhood_pipeline._parse_tree(tree_path)
            targets = [
                {'row': row, 'tree_leaf_id': row['id']}
                for row in rows
            ]
            metadata = neighbourhood_pipeline._load_metadata(
                db, ['Q1', 'Q2', 'H1', 'P1', 'P2'], rows,
            )
            geometry = neighbourhood_pipeline._figure_geometry(
                root,
                targets,
                metadata,
                {'Hungate'},
                aligned_sequences=dict(read_fasta(alignment)),
                identity_anchor='Q1',
                baseline_context={'H1': {'identity': 98.50, 'queries': ['Q1', 'Q2']}},
            )
            labels = '\n'.join(geometry['display_labels'])
            self.assertIn('[P1] Q1', labels)
            self.assertIn('[B2] Q2', labels)
            self.assertIn('H1', labels)
            self.assertIn('nearest baseline context; max 98.50% vsearch', labels)
            self.assertIn('MSA pident 75.00%', labels)
            pairwise = tmp / 'neighbourhoods' / 'clade_001_pairwise_pident.tsv'
            self.assertTrue(pairwise.exists())
            self.assertIn('Q1\tQ2\t75.00\t4', pairwise.read_text())
            with open(result['manifest']) as handle:
                manifest_rows = list(csv.DictReader(handle, delimiter='\t'))
            self.assertEqual(manifest_rows[0]['NearestBaselineHitsShown'], 'H1')
            self.assertEqual(manifest_rows[0]['IdentityAnchor'], 'Q1')
            self.assertEqual(manifest_rows[0]['BaselineLeavesShown'], '1')
            self.assertEqual(manifest_rows[0]['TreeLeavesShown'], '3')
            self.assertEqual(
                rows[0]['local_pairwise_pident_table'],
                'neighbourhoods/clade_001_pairwise_pident.tsv',
            )

    def test_msa_pident_excludes_terminal_gaps_and_ambiguous_columns(self):
        identity, compared = neighbourhood_pipeline._msa_pident('ACGT--N', 'ACG---A')
        self.assertEqual(compared, 3)
        self.assertAlmostEqual(identity, 100.0)

    def test_cluster_reports_write_one_consolidated_csv(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            rows = [
                {
                    'id': 'A', 'cluster_representative': 'self', 'cluster_size': '2',
                    'clustered_members': 'B', 'taxonomy': 'd__Bacteria; g__Alpha',
                    'classification_identity': '99.0', 'nearest_identity': '98.0',
                    'matches_ge_99': '0', 'matches_ge_97': '2', 'matches_ge_95': '3',
                    'novelty_score': '40.0', 'crowding': 'sparse',
                    'sequencing_priority': 'MEDIUM', 'in_tree': 'Yes', 'placement_flags': '',
                },
                {
                    'id': 'B', 'cluster_representative': 'A', 'cluster_size': '2',
                    'clustered_members': '', 'taxonomy': 'd__Bacteria; g__Alpha',
                    'classification_identity': '98.8', 'nearest_identity': '97.8',
                    'matches_ge_99': '0', 'matches_ge_97': '2', 'matches_ge_95': '3',
                    'novelty_score': '42.0', 'crowding': 'sparse',
                    'sequencing_priority': 'MEDIUM', 'in_tree': 'No', 'placement_flags': '',
                },
                {
                    'id': 'C', 'cluster_representative': 'self', 'cluster_size': '1',
                    'clustered_members': '', 'taxonomy': 'd__Bacteria; g__Beta',
                    'classification_identity': '99.2', 'nearest_identity': '99.1',
                    'matches_ge_99': '1', 'matches_ge_97': '4', 'matches_ge_95': '5',
                    'novelty_score': '20.0', 'crowding': 'moderate',
                    'sequencing_priority': 'LOW', 'in_tree': 'Yes', 'placement_flags': '',
                },
            ]

            summary, clusters_csv, backups = cluster_report_pipeline.generate_cluster_reports(
                tmpdir,
                rows,
                tree_path=None,
            )

            self.assertTrue(Path(summary).exists())
            self.assertTrue(Path(backups).exists())
            self.assertEqual(Path(clusters_csv).name, 'clusters.csv')
            self.assertFalse((Path(tmpdir) / 'clusters').exists())
            cluster_text = Path(clusters_csv).read_text()
            self.assertIn('ClusterID,ID,IsRepresentative,BackupRank', cluster_text)
            self.assertIn('A,A,True', cluster_text)
            self.assertIn('A,B,False', cluster_text)
            self.assertIn('C,C,True', cluster_text)


if __name__ == '__main__':
    unittest.main()
