"""BranchManager CLI — lightweight entrypoint for the branchmanager package.

Implements the `preload`, `run`/`evaluate`, and `regen-itol` commands for the src/ layout.
"""
import argparse
import logging
import os
import shutil
import sys
from pathlib import Path

# When invoked directly (python src/branchmanager/PhenGO-Predict.py) the package root (src)
# may not be on sys.path. Ensure the parent of this `branchmanager` package is
# available so absolute imports like `branchmanager.db.interface` work.
try:
    here = Path(__file__).resolve().parent
    src_root = str(here.parent)
    if src_root not in sys.path:
        sys.path.insert(0, src_root)
except Exception:
    pass

from branchmanager.db.interface import Database
from branchmanager.pipeline import classify, tree, itol, qc, derep, novelty
from branchmanager.pipeline import cluster_report as _cluster_report
from branchmanager.pipeline import mwl as _mwl
from branchmanager.pipeline.classify import _derive_db_name as _classify_derive_db_name
from branchmanager.pipeline.collapse import collapse_fasta_within_taxa
from branchmanager.pipeline.workflow_helpers import (
    _assignment_source_is_fasta,
    build_orig_to_short_map as _build_orig_to_short_map_helper,
    build_placement_warning_rows,
    build_sequence_assessment_rows,
    classification_ids_matching_kingdom as _classification_ids_matching_kingdom_helper,
    collect_db_taxonomy_rows,
    iter_assignment_rows,
    load_classification_results_for_dataset,
    load_taxonomy_entries_from_assignments,
    merge_combined_taxonomy_rows,
    prune_dataset_by_kingdom as _prune_dataset_by_kingdom_helper,
    read_combined_taxonomy_ids,
    write_combined_taxonomy_tsv,
    write_baseline_hits_tsv,
    write_placement_warning_tsv,
    write_selection_summary_tsv,
    write_sequence_assessment_tsv,
)
from branchmanager.taxonomy import canonicalize_sequence_id, normalize_domain_query, taxonomy_matches_kingdom
from branchmanager.partner_metadata import load_partner_sequencing_metadata


SEQUENCE_DOMAIN_CHOICES = (
    'bacteria', 'bacterial',
    'archaea', 'archaeal',
    'fungi', 'fungal',
    'eukaryota', 'eukarya',
    'mixed', 'all', 'none',
)

def _find_tree_file_in_dir(d: str):
    """Return path to a tree file in directory d if present, preferring current_tree.nwk."""
    p = Path(d)
    for cand in (p / 'current_tree.nwk', p / 'tree' / 'current_tree.nwk'):
        if cand.exists():
            return str(cand)
    # otherwise search for any .nwk or .tree file in the root or organised tree dir
    for base in (p, p / 'tree'):
        for ext in ('*.nwk', '*.tree', '*.tre'):
            found = next(base.glob(ext), None)
            if found:
                return str(found)
    return None


def _configure_logging(outdir: str):
    """Configure root logger: console INFO + optional file DEBUG in outdir.

    The root logger level must be set to DEBUG so that INFO/DEBUG records
    actually reach the handlers — handler-level filtering alone is not enough
    because Python's logging framework gates records at the logger level first.
    """
    logger = logging.getLogger()
    # remove existing handlers to avoid duplicate messages on repeated calls
    for h in list(logger.handlers):
        logger.removeHandler(h)
    # Set root logger to DEBUG so all records flow through to handlers;
    # each handler then applies its own level filter.
    logger.setLevel(logging.DEBUG)
    fmt = logging.Formatter("%(asctime)s %(levelname)s: %(message)s")
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(fmt)
    logger.addHandler(ch)
    try:
        fh = logging.FileHandler(os.path.join(outdir, 'branchmanager.log'))
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(fmt)
        logger.addHandler(fh)
    except Exception:
        logger.warning("Could not write log file to %s", outdir)


def _build_orig_to_short(alias_entries):
    return _build_orig_to_short_map_helper(alias_entries, Database(':memory:'))


def _write_id_map_tsv(path: str | Path, entries, *, short_header: str = 'short_id', original_header: str = 'original_header'):
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, 'w') as handle:
        handle.write(f'{short_header}\t{original_header}\n')
        for short, original in entries:
            handle.write(f'{short}\t{original}\n')
    return str(p)


def _dataset_sequence_ids(db: Database, dataset: str) -> set[str]:
    with db.connect() as conn:
        cur = conn.cursor()
        cur.execute("SELECT id FROM sequences WHERE dataset = ?", (dataset,))
        return {str(iid) for (iid,) in cur.fetchall()}


def _filter_fasta_to_ids(src: str | Path, dst: str | Path, allowed_ids: set[str]) -> int:
    from branchmanager.utils.fasta import read_fasta, write_fasta
    records = [(h, s) for h, s in read_fasta(str(src)) if str(h) in allowed_ids]
    write_fasta(records, str(dst))
    return len(records)


def _write_partner_metadata_warnings(path: str | Path, warning_rows):
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, 'w') as handle:
        handle.write('SourceID\tReason\n')
        for source_id, reason in warning_rows:
            handle.write(f'{source_id}\t{reason}\n')
    return str(p)


def _load_partner_metadata_for_run(args, db: Database, outdir: str, orig_to_short: dict, run_ids):
    metadata_path = getattr(args, 'partner_metadata', None)
    command = getattr(args, 'command', None)
    if not metadata_path:
        if command in ('evaluate', 'eval'):
            raise SystemExit(
                '[RUN] evaluate requires --partner-metadata / --sequencing-metadata: '
                'a CSV/TSV sidecar table with sequence IDs, partner IDs, and WGS-selected status.'
            )
        return {}

    log = logging.getLogger(__name__)
    try:
        metadata_rows = load_partner_sequencing_metadata(metadata_path)
    except Exception as e:
        raise SystemExit(f'[RUN] Failed to read partner metadata {metadata_path}: {e}')

    run_id_set = {str(x) for x in run_ids}
    resolved_rows = []
    warnings = []
    matched_sources = set()
    for row in metadata_rows:
        source_id = str(row.get('source_id') or '').strip()
        if not source_id:
            continue
        mapped = orig_to_short.get(source_id)
        if not mapped:
            try:
                cid = canonicalize_sequence_id(source_id)
            except Exception:
                cid = None
            if cid:
                mapped = orig_to_short.get(cid)
        if not mapped and source_id in run_id_set:
            mapped = source_id
        if not mapped:
            existing = db.resolve_sequence_id(source_id)
            if existing in run_id_set:
                mapped = existing
        if not mapped:
            warnings.append((source_id, 'metadata_id_not_found_in_current_run'))
            continue

        matched_sources.add(str(mapped))
        resolved_rows.append({
            'id': mapped,
            'partner_id': row.get('partner_id') or source_id,
            'dataset': getattr(args, 'dataset', ''),
            'selected_for_wgs': bool(row.get('selected_for_wgs')),
            'source_id': source_id,
            'source_file': str(metadata_path),
            'raw_selected_value': row.get('raw_selected_value', ''),
        })

    for run_id in sorted(run_id_set - matched_sources):
        warnings.append((run_id, 'run_sequence_missing_from_partner_metadata'))

    inserted = db.upsert_sequencing_metadata(resolved_rows)
    log.info(
        '[RUN] Loaded partner sequencing metadata from %s: %d matched rows, %d warning(s)',
        metadata_path,
        inserted,
        len(warnings),
    )
    if warnings:
        warn_path = _write_partner_metadata_warnings(
            Path(outdir) / 'partner_metadata_warnings.tsv',
            warnings,
        )
        log.warning('[RUN] Partner metadata warnings written to %s', warn_path)

    return db.get_sequencing_metadata_for_ids(run_ids)


def _load_id_map_from_tsv(path: str | Path, db: Database | None = None):
    mapping = {}
    with open(path) as handle:
        next(handle, None)
        for line in handle:
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            short, orig = parts[0], parts[1]
            if not short or not orig:
                continue
            mapping[orig] = short
            mapping[short] = short
            if db is not None:
                try:
                    cid = db._canonical_from_header(orig)
                    if cid:
                        mapping[cid] = short
                except Exception:
                    pass
    return mapping


def _find_preferred_id_map(directory: str | Path):
    p = Path(directory)
    preferred = [
        p / 'preload_id_map.tsv',
        p / 'user_id_map.tsv',
        p / 'user_id_map.csv',
        p / 'ids' / 'preload_id_map.tsv',
        p / 'ids' / 'user_id_map.tsv',
        p / 'ids' / 'user_id_map.csv',
        p / 'ids' / 'baseline_id_map.tsv',
    ]
    for cand in preferred:
        if cand.exists():
            return cand
    try:
        for base in (p, p / 'ids'):
            for cand in base.glob('*_id_map.tsv'):
                if cand.name != 'id_map.tsv':
                    return cand
    except Exception:
        pass
    return None


def _combined_taxonomy_candidates(directory: str | Path):
    p = Path(directory)
    return [
        p / 'preload_combined_taxonomy.tsv',
        p / 'combined_taxonomy.tsv',
        p / 'taxonomy' / 'tree_taxonomy.tsv',
        p / 'baseline' / 'loaded_baseline_taxonomy.tsv',
    ]


def _read_qc_stats(outdir: Path) -> dict:
    stats = {}
    for cand in (outdir / 'assessment' / 'qc.stats', outdir / 'qc.stats'):
        if not cand.exists():
            continue
        try:
            with open(cand) as handle:
                for line in handle:
                    line = line.rstrip('\n')
                    if not line or '\t' not in line:
                        continue
                    key, value = line.split('\t', 1)
                    stats[key] = value
            stats['_path'] = str(cand)
            return stats
        except Exception:
            return {}
    return {}


def _rel_output_path(root: Path, path: str | Path) -> str:
    try:
        return str(Path(path).relative_to(root))
    except Exception:
        return str(path)


def _write_detailed_output_guide(outdir: Path, rows):
    """Write a detailed Markdown guide for files and core metrics."""
    guide = outdir / 'OUTPUT_GUIDE.md'
    qc_stats = _read_qc_stats(outdir)
    qc_rejections = outdir / 'assessment' / 'qc_rejections.tsv'
    if not qc_rejections.exists():
        qc_rejections = outdir / 'qc_rejections.tsv'

    def stat(name: str, default: str = 'NA') -> str:
        return qc_stats.get(name, default)

    file_rows = []
    for name, path, expl in rows:
        rel = _rel_output_path(outdir, path)
        if rel in ('OUTPUT_GUIDE.md',):
            continue
        safe_expl = str(expl).replace('|', '/').replace('\t', ' ')
        file_rows.append(f"| `{rel}` | {safe_expl} |")

    qc_section = [
        "## QC Filtering",
        "",
        "QC happens before dereplication, classification, novelty scoring, and tree building. "
        "Sequences that fail QC are not included in downstream reports because they are too short "
        "for reliable 16S placement or contain too many ambiguous bases.",
        "",
        "Filtering rules:",
        "",
        f"- `too_short`: sequence length is less than `min_len` (`{stat('min_len')}` bp for this run).",
        f"- `too_many_n`: sequence has more than `max_n` ambiguous `N` bases (`{stat('max_n')}` for this run).",
        "- A sequence can have more than one reason in `qc_rejections.tsv`.",
        "",
    ]
    if qc_stats:
        qc_section.extend([
            "This run:",
            "",
            f"- Input sequences: `{stat('total_input')}`",
            f"- Kept after QC: `{stat('kept')}`",
            f"- Rejected total: `{stat('rejected_total', str(max(0, int(stat('total_input', '0')) - int(stat('kept', '0'))) if stat('total_input', '0').isdigit() and stat('kept', '0').isdigit() else 'NA'))}`",
            f"- Rejected as too short: `{stat('rejected_too_short')}`",
            f"- Rejected for too many `N` bases: `{stat('rejected_too_many_n')}`",
            "",
        ])
        if qc_rejections.exists():
            qc_section.append(f"See `{_rel_output_path(outdir, qc_rejections)}` for every rejected sequence, its length, N count, and exact reason.")
        else:
            qc_section.append("Older runs may only include the summary file `qc.stats`; rerun to generate per-sequence `qc_rejections.tsv`.")
        qc_section.append("")
    else:
        qc_section.extend([
            "No `qc.stats` file was found in this output directory. If this was a preload-only or partially completed run, QC may not have been executed.",
            "",
        ])

    lines = [
        "# BranchManager Output Guide",
        "",
        "This guide explains the output files and the main metrics produced by `branchmanager run` / `branchmanager evaluate`.",
        "",
        "## Recommended Reading Order",
        "",
        "1. `assessment/sequence_assessment.tsv` - primary per-sequence decision table.",
        "2. `baseline/baseline_hits.tsv` - closest cultured/baseline hits such as Hungate.",
        "3. `taxonomy/` - taxonomy assignments for each configured reference database.",
        "4. `tree/current_tree.nwk` plus `tree/current_alignment.fasta` and `itol/*.itol` - upload to iTOL for visual inspection.",
        "5. `assessment/novelty_metrics.tsv` - detailed novelty calculations behind the assessment table.",
        "",
        "## Directory Overview",
        "",
        "- `assessment/`: primary reports, novelty metrics, cluster reports, warnings, and raw all-known nearest-hit table.",
        "- `baseline/`: nearest-hit reports and loaded baseline taxonomy/sequences.",
        "- `taxonomy/`: per-database classification outputs and combined taxonomy used for tree metadata.",
        "- `tree/`: MSA, Newick tree, and tree/alignment warning files.",
        "- `itol/`: one iTOL metadata dataset per metadata type, usually colorstrip files.",
        "- `ids/`: short ID to original FASTA header maps.",
        "- `intermediate/`: QC/dereplication/collapse/debug FASTA files and scratch pools.",
        "- `logs/`: pipeline log.",
        "",
        *qc_section,
        "## Main Assessment Metrics",
        "",
        "`assessment/sequence_assessment.tsv` is the main table. Important column groups:",
        "",
        "- `Taxonomy`, `ClassificationHit`, `ClassificationIdentity`, `ClassificationConfidence`: primary reference-database assignment, usually GTDB when `--main-ref GTDB` is used.",
        "- `Taxonomy_<DB>`, `ClassificationHit_<DB>`, `Identity_<DB>`, `Confidence_<DB>`: assignment from additional databases such as GG2, SILVA, or NCBI.",
        "- `NearestHit`, `NearestHitDataset`, `NearestHitTaxonomy`, `NearestIdentity`: closest sequence in the baseline/cultured pool, for example Hungate.",
        "- `AllKnownNearestHit`, `AllKnownNearestIdentity`, `AllKnownNoveltyScore`: same idea, but against all non-current datasets in the project DB.",
        "- `ReferenceNearestHit`, `ReferenceNearestIdentity`, `ReferenceNoveltyScore`: nearest hit and novelty score against the selected external reference FASTA, usually GTDB.",
        "- `PartnerID`, `SelectedForGenomeSequencing`, `CladeAlreadySelectedForGenomeSequencing`, `GenomeSequencingAdjustedPriority`: rolling partner/WGS-selection context from `--partner-metadata`.",
        "- `InTree`, `ClusterRepresentative`, `ClusterSize`, `ClusteredMembers`: whether the sequence itself entered the tree or was represented by another clustered sequence.",
        "- `PlacementFlags`: warnings such as low classification identity, low nearest identity, or novelty/classification disagreement.",
        "",
        "## Novelty Metrics",
        "",
        "`assessment/novelty_metrics.tsv` contains cultured-baseline novelty, all-known novelty, and external reference novelty.",
        "",
        "- Leading `Nearest*`, `NoveltyScore`, `Crowding`, and `SequencingPriority` columns mirror the baseline/cultured comparison when a baseline pool exists; otherwise they mirror all-known novelty.",
        "- `Baseline*` columns compare only against explicit baseline datasets such as Hungate or datasets supplied with `--novelty-baseline-dataset`.",
        "- `AllKnown*` columns compare against every DB dataset except the current run. This includes Hungate plus non-Hungate/prior partner datasets.",
        "- `Reference*` columns compare against the chosen external reference FASTA supplied with `--ref`, usually GTDB. These are separate from the baseline/project novelty scores.",
        "- `SelectedGenome*` and `GenomeSequencing*` columns use the rolling sequencing metadata table in the SQLite DB. A selected neighbour at >=97 percent identity marks the local 16S clade as already represented for genome sequencing.",
        "- `NearestIdentity`: vsearch global-alignment percent identity to the nearest sequence in that pool.",
        "- `Novel`: `True` when nearest identity is below 97 percent.",
        "- `MatchesGE99`, `MatchesGE97`, `MatchesGE95`: number of pool sequences at or above 99, 97, and 95 percent identity. These describe how busy the local neighbourhood is.",
        "- `Crowding`: `isolated` when there is at most one hit at both 99 and 97 percent; `sparse` when <=3 hits at 97 percent; `moderate` when <=10 hits at 97 percent; otherwise `crowded`.",
        "- `NoveltyScore`: 0-100 score where higher means more novel and less crowded. It combines distance from the nearest hit with density bonuses for sparse neighbourhoods.",
        "- `SequencingPriority`: `HIGH` for <97 percent identity with few close neighbours, `MEDIUM` for moderately novel/sparse cases, otherwise `LOW`.",
        "- `DensitySource`: names the pool used, for example `baseline:Hungate`, `all_known`, `target_fasta`, or `reference_fasta` fallback.",
        "",
        "Interpretation: a sequence far from Hungate but close to non-Hungate isolates is likely novel relative to cultured rumen isolate collections, but not necessarily novel relative to everything already supplied to the project.",
        "",
        "## Taxonomy Metrics",
        "",
        "- `ClassificationIdentity` is the vsearch percent identity to the best reference hit used for taxonomy assignment.",
        "- `ClassificationConfidence` is derived from the assignment parser/classifier output; higher means stronger taxonomy support.",
        "- `taxonomy/all_databases.tsv` shows assignments across all configured databases in one place.",
        "- `taxonomy/input_warnings.tsv` flags mismatches between reference FASTA IDs and supplied taxonomy tables.",
        "",
        "## MWL Metrics",
        "",
        "When `--mwl` is supplied, MWL columns are added to `sequence_assessment.tsv` and `assessment/mwl_matches.tsv` is written.",
        "",
        "- `MWLMatch`: whether the GTDB-derived taxonomy matched a Most Wanted List taxon.",
        "- `MWLMatchedRank` and `MWLMatchedTaxon`: deepest rank/taxon that matched the MWL entry.",
        "- `MWLTaxonomicScore`, `MWLIdentity`, `MWLScore`: taxonomic and identity contributions to MWL priority.",
        "- `EvaluationScore`: combined candidate score after MWL contribution is considered.",
        "",
        "## Tree And iTOL Files",
        "",
        "- `tree/current_alignment.fasta`: MSA used to build the tree.",
        "- `tree/current_tree.nwk`: tree to upload to iTOL.",
        "- `itol/phylum.itol`, `itol/family.itol`, `itol/genus.itol`: taxonomy colorstrips.",
        "- `itol/dataset_membership.itol`: dataset-of-origin colorstrip.",
        "- `itol/novelty.itol`: novelty/nearest-identity colorstrip.",
        "",
        "BranchManager keeps iTOL `DATASET_COLORSTRIP` files only. Older branch `TREE_COLORS` and symbol-strip variants duplicated the same metadata and are removed to keep outputs readable.",
        "",
        "## File Catalogue",
        "",
        "| File | Explanation |",
        "|---|---|",
        *file_rows,
        "",
    ]
    try:
        guide.write_text('\n'.join(lines))
    except Exception:
        pass


def _write_output_explanations(outdir: str):
    """Write short human-readable explanation files for common outputs under outdir.

    For each file discovered under outdir we create a companion file named
    <filename>.explain.txt containing a brief description of the file and how
    it is produced. This helps downstream users and CI to understand outputs.
    """
    p = Path(outdir)
    if not p.exists():
        return
    # mapping of filename substrings (lowercase) to explanation text
    patterns = [
        ('tree_taxonomy.tsv', 'Combined taxonomy TSV used for tree/iTOL metadata. Columns: ID\tTaxon\tConfidence.'),
        ('combined_taxonomy.tsv', 'Combined taxonomy TSV. Columns: ID\tTaxon\tConfidence. Used to generate iTOL color/legend files.'),
        ('loaded_baseline_taxonomy.tsv', 'Combined taxonomy TSV for baseline/provided datasets loaded before evaluate. Columns: ID\tTaxon\tConfidence.'),
        ('preload_combined_taxonomy.tsv', 'Combined taxonomy TSV for preload dataset. Columns: ID\tTaxon\tConfidence.'),
        ('all_databases.tsv', 'Per-sequence taxonomic assignments from every configured reference database.'),
        ('baseline_hits.tsv', 'Nearest-hit report against baseline/provided datasets such as Hungate.'),
        ('nearest_all_known_hits_raw.tsv', 'Raw nearest-hit novelty table against all non-current DB datasets; baseline-specific hits are summarised in baseline_hits.tsv and novelty_metrics.tsv.'),
        ('qc_rejections.tsv', 'Per-sequence QC rejection table. Columns: ID, Length, NCount, Reasons, MinLength, MaxN. Reasons explain exactly why each sequence was filtered before downstream analysis.'),
        ('partner_metadata_warnings.tsv', 'Warnings from --partner-metadata mapping. Rows indicate metadata IDs not found in the current run, or run sequences missing from the partner metadata table.'),
        ('qc.stats', 'QC summary with total input, kept count, rejection counts, min_len, max_n, and pointer to qc_rejections.tsv when available.'),
        ('user_id_map.tsv', 'Mapping of runtime sequence IDs to original headers produced when inserting user sequences into the DB. When shortening is disabled these usually match.'),
        ('preload_id_map.tsv', 'Mapping of preload runtime IDs back to original FASTA headers. Use this to trace tree labels back to source records.'),
        ('*_id_map.tsv', 'ID map mapping original headers to runtime DB ids. Useful for iTOL and metadata tracing.'),
        ('collapsed_map.tsv', 'Cluster map for collapsed sequences: rep_id\ttaxonomy\tcount.'),
        ('preload_collapsed_map.tsv', 'Cluster map for preload collapsed sequences: rep_id\ttaxonomy\tcount.'),
        ('collapsed_members.tsv', 'Member->representative mapping (member\trep) for collapsed clusters.'),
        ('preload_collapsed_members.tsv', 'Member->representative mapping for preload collapsed clusters.'),
        ('derep_short.fasta', 'Dereplicated FASTA where sequence headers are the runtime IDs used by the DB. These are preserved source IDs unless --shorten-ids was requested.'),
        ('derep_short_collapsed.fasta', 'Dereplicated FASTA after collapse; representatives for clusters kept with runtime IDs.'),
        ('preload_short_collapsed.fasta', 'Collapsed preload FASTA; representatives retained for tree building.'),
        ('novelty_matches.tsv', 'vsearch BLAST-like output used to compute nearest-neighbour novelty identities.'),
        ('novelty_metrics.tsv', (
            'Per-sequence novelty metrics. The leading Nearest*/NoveltyScore columns mirror the '
            'baseline/cultured comparison when available; explicit Baseline* columns compare against '
            'datasets such as Hungate, and AllKnown* columns compare against every non-current dataset '
            'stored in the DB. Reference* columns compare against the external reference FASTA, usually '
            'GTDB. PartnerID/SelectedForGenomeSequencing and SelectedGenome*/GenomeSequencing* columns '
            'add rolling WGS-selection context. DensitySource columns name the comparison pool.'
        )),
        ('sequence_assessment.tsv', (
            'Unified per-sequence assessment. '
            'COLUMN GROUPS: '
            '(1) TAXONOMY/CLASSIFICATION — Taxonomy, ClassificationHit, ClassificationIdentity, '
            'ClassificationConfidence: derived from the primary reference database (GTDB/SILVA). '
            'ClassificationHit is the reference accession vsearch matched. '
            'Repeated as Taxonomy_<DB>, ClassificationHit_<DB>, Identity_<DB>, Confidence_<DB> '
            'for each additional --alt-ref database. '
            '(2) NOVELTY — NearestHit, NearestIdentity, MatchesGE*, NoveltyScore, Crowding, '
            'SequencingPriority: baseline/cultured novelty when a baseline pool exists. '
            'AllKnown* columns repeat the same metrics against all non-current DB datasets. '
            'Reference* columns repeat the same metrics against the external taxonomy reference '
            'FASTA, usually GTDB. PartnerID/SelectedForGenomeSequencing and GenomeSequencing* '
            'columns report whether this isolate or a nearby 16S clade has already been selected '
            'for WGS and provide a WGS-aware adjusted priority. '
            '(3) TREE/CLUSTER — InTree, ClusterRepresentative, ClusterSize, ClusteredMembers: '
            'records whether the sequence entered the phylogenetic tree directly or was '
            'represented by a cluster representative after --collapse. '
            '(4) MWL — when --mwl is supplied, EvaluationScore and MWL* columns describe '
            'GTDB-based Most Wanted List matches and MWL priority contribution.'
        )),
        ('selection_summary.tsv', (
            'Concise scientific-advisory-board selection table. Contains one row per evaluated '
            'sequence with partner acronym, recommendation, WGS-adjusted priority, key novelty '
            'identities, taxonomy/MWL evidence, selected-clade status, and a short rationale. '
            'Use sequence_assessment.tsv for the full audit trail.'
        )),
        ('mwl_matches.tsv', 'Most Wanted List match report. Contains sequences whose GTDB taxonomy matched an MWL taxon, with matched rank, MWL score, evaluation score, and functional role.'),
        ('taxonomy_input_warnings.tsv', 'Warnings about inconsistencies between the classifier reference FASTA and supplied taxonomy table.'),
        ('tree_build_warnings.tsv', 'Warnings about weak phylogenetic signal, missing anchors, or poor alignment quality.'),
        ('tree_orientation_summary.tsv', 'Sequence-level audit of tree-input orientation checks. Reports which sequences were kept forward, reverse-complemented, or lacked orientation evidence before alignment.'),
        ('placement_warnings.tsv', 'Warnings about low-support placements, low identity matches, or potentially artefactual novelty assignments.'),
        ('taxa_assignments_classout.tsv', 'Synthetic classification-like TSV created when --taxa-assignments provided. Columns: id\tbest\tidentity\ttaxon\tconfidence.'),
        ('dataset_membership.itol', 'iTOL DATASET_COLORSTRIP mapping sequence IDs to dataset colors (membership).'),
        ('novelty.itol', 'iTOL colorstrip showing novelty (nearest identity) for run sequences.'),
        ('preload_dataset.itol', 'iTOL colorstrip for the preload dataset; maps preload ids to the dataset color.'),
        ('itol_dataset_membership.itol', 'iTOL DATASET_COLORSTRIP mapping sequence IDs to dataset colors (membership).'),
        ('itol_novelty.itol', 'iTOL colorstrip showing novelty (nearest identity) for run sequences.'),
        ('itol_dataset_preload.itol', 'iTOL colorstrip for the preload dataset; maps preload ids to the dataset color.'),
        ('.nwk', 'Newick tree file (phylogenetic tree). Commonly named current_tree.nwk.'),
        ('.itol', 'iTOL dataset file (text format) describing colors/strips/legends for visualization in iTOL.'),
        ('.log', 'Log file produced by the pipeline (branchmanager.log) containing debug/info messages')
    ]

    rows = []
    # iterate files in outdir (non-recursive and recursive) and collect explanations
    for fp in sorted(p.rglob('*')):
        if fp.is_dir():
            continue
        if fp.name.endswith('.explain.txt'):
            try:
                fp.unlink()
            except Exception:
                pass
            continue
        name = fp.name
        lname = name.lower()
        expl = None
        for patt, text in patterns:
            # treat wildcard at start/end
            if patt.startswith('*') and patt.endswith('*'):
                key = patt.strip('*').lower()
                if key in lname:
                    expl = text
                    break
            elif patt.startswith('*'):
                key = patt.lstrip('*').lower()
                if lname.endswith(key):
                    expl = text
                    break
            elif patt.endswith('*'):
                key = patt.rstrip('*').lower()
                if lname.startswith(key):
                    expl = text
                    break
            elif patt.startswith('.'):
                if lname.endswith(patt):
                    expl = text
                    break
            else:
                if patt in lname:
                    expl = text
                    break
        if expl is None:
            # generic descriptions based on extension
            if lname.endswith('.fasta') or lname.endswith('.fa'):
                expl = 'FASTA file containing sequences. May be dereplicated, collapsed, or exported from the DB.'
            elif lname.endswith('.tsv') or lname.endswith('.csv'):
                expl = 'Tab/CSV-delimited table used by the pipeline.'
            elif lname.endswith('.uc'):
                expl = 'vsearch UC membership file describing clusters.'
            else:
                expl = 'Pipeline output file.'

        rows.append((fp.name, str(fp), expl))

    try:
        manifest = p / 'OUTPUT_EXPLANATIONS.tsv'
        with open(manifest, 'w') as ef:
            ef.write('File\tPath\tDescription\n')
            for name, path, expl in rows:
                ef.write(f"{name}\t{path}\t{expl}\n")
    except Exception:
        pass
    _write_detailed_output_guide(p, rows)


def _replace_path(src: Path, dst: Path):
    if not src.exists():
        return None
    try:
        if src.resolve() == dst.resolve():
            return str(dst)
    except Exception:
        pass
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        if dst.is_dir():
            shutil.rmtree(dst)
        else:
            dst.unlink()
    shutil.move(str(src), str(dst))
    return str(dst)


def _move_if_exists(outdir: Path, name: str, target_dir: str, target_name: str | None = None):
    return _replace_path(outdir / name, outdir / target_dir / (target_name or name))


def _move_glob(outdir: Path, pattern: str, target_dir: str, rename=None):
    moved = []
    for src in sorted(outdir.glob(pattern)):
        if not src.exists() or src.is_dir():
            continue
        target_name = rename(src) if rename else src.name
        dst = outdir / target_dir / target_name
        got = _replace_path(src, dst)
        if got:
            moved.append(got)
    return moved


def _unlink_glob(outdir: Path, pattern: str):
    for src in sorted(outdir.glob(pattern)):
        try:
            if src.is_dir():
                shutil.rmtree(src)
            else:
                src.unlink()
        except Exception:
            pass


def _organize_run_outputs(outdir: str, *, primary_db_name: str | None = None):
    """Clean and organise evaluate/run outputs into high-level report folders."""
    out = Path(outdir)
    if not out.exists():
        return

    # Remove redundant iTOL representations and bulky/raw classifier caches.
    for patt in (
        'itol_*_symbols.itol',
        'itol_*_tree_colors.txt',
        'tree_colors_with_clades.txt',
        'itol_combined_colors.csv',
        'matches*.tsv',
        'novelty_matches.tsv',
        'novelty_density_matches.tsv',
        'novelty_*_matches.tsv',
        'ref_uncompressed*.fasta',
    ):
        _unlink_glob(out, patt)

    # Overall assessment and prioritisation reports.
    for name in (
        'sequence_assessment.tsv',
        'cluster_summary.tsv',
        'backup_candidates.tsv',
        'placement_warnings.tsv',
        'novelty_metrics.tsv',
        'selection_summary.tsv',
        'rumen_functions_draft.tsv',
        'qc.stats',
        'qc_rejections.tsv',
        'partner_metadata_warnings.tsv',
    ):
        _move_if_exists(out, name, 'assessment')
    _replace_path(out / 'clusters', out / 'assessment' / 'clusters')
    _move_if_exists(out, 'mwl_matches.tsv', 'assessment')

    # Direct baseline/provided-dataset nearest-hit reports.
    _move_if_exists(out, 'baseline_hits.tsv', 'baseline')
    _move_if_exists(out, 'novelty.tsv', 'assessment', 'nearest_all_known_hits_raw.tsv')

    # Taxonomy assignment reports.
    primary_name = primary_db_name or 'primary'
    safe_primary = ''.join(c if c.isalnum() or c in ('_', '-') else '_' for c in str(primary_name)).strip('_') or 'primary'
    _move_if_exists(out, 'taxonomy.tsv', 'taxonomy', f'{safe_primary}.tsv')
    _move_if_exists(out, 'taxonomy_all_dbs.tsv', 'taxonomy', 'all_databases.tsv')
    _move_if_exists(out, 'combined_taxonomy.tsv', 'taxonomy', 'tree_taxonomy.tsv')
    _move_if_exists(out, 'taxonomy_input_warnings.tsv', 'taxonomy', 'input_warnings.tsv')
    _move_if_exists(out, 'taxa_assignments_classout.tsv', 'taxonomy', 'input_assignments.tsv')

    def _rename_taxonomy(src: Path):
        stem = src.stem
        if stem.startswith('taxonomy_'):
            stem = stem[len('taxonomy_'):]
        return f'{stem}.tsv'

    _move_glob(out, 'taxonomy_*.tsv', 'taxonomy', rename=_rename_taxonomy)

    # Tree/MSA deliverables for iTOL upload.
    for name in (
        'current_tree.nwk',
        'current_tree_labeled.nwk',
        'current_alignment.fasta',
        'tree_build_warnings.tsv',
        'tree_orientation_summary.tsv',
    ):
        _move_if_exists(out, name, 'tree')
    _move_glob(out, 'tree_sequences_phylum_*.fasta', 'tree')

    # One iTOL metadata dataset per metadata type.
    itol_renames = {
        'itol_phylum_colors.itol': 'phylum.itol',
        'itol_family_colors.itol': 'family.itol',
        'itol_genus_colors.itol': 'genus.itol',
        'itol_dataset_membership.itol': 'dataset_membership.itol',
        'itol_novelty.itol': 'novelty.itol',
        'itol_user_colors.itol': 'user_colors.itol',
        'itol_dataset_preload.itol': 'preload_dataset.itol',
    }
    for src_name, dst_name in itol_renames.items():
        _move_if_exists(out, src_name, 'itol', dst_name)

    def _rename_function_itol(src: Path):
        stem = src.stem
        if stem.startswith('itol_func_'):
            stem = 'function_' + stem[len('itol_func_'):]
        return stem + src.suffix

    _move_glob(out, 'itol_func_*.itol', 'itol', rename=_rename_function_itol)

    # Stable ID maps and run-sequence processing artefacts.
    for name in ('user_id_map.tsv', 'preload_id_map.tsv'):
        _move_if_exists(out, name, 'ids')

    # Baseline preload internals: keep useful reports, drop raw matches/ref cache.
    baseline_preload = out / 'baseline_preload'
    if baseline_preload.exists():
        _replace_path(baseline_preload / 'baseline_id_map.tsv', out / 'ids' / 'baseline_id_map.tsv')
        _replace_path(
            baseline_preload / 'baseline_combined_taxonomy.tsv',
            out / 'baseline' / 'loaded_baseline_taxonomy.tsv',
        )
        _replace_path(
            baseline_preload / 'taxonomy.tsv',
            out / 'baseline' / 'loaded_baseline_classification.tsv',
        )
        for src in sorted(baseline_preload.glob('preload_*_seqs.fasta')):
            _replace_path(src, out / 'baseline' / src.name)
        try:
            shutil.rmtree(baseline_preload)
        except Exception:
            pass

    for name in (
        'qc.fasta',
        'derep.fasta',
        'derep_short.fasta',
        'derep_short_collapsed.fasta',
        'collapsed_map.tsv',
        'collapsed_members.tsv',
        'submitted_sequences.fasta',
        'db_preload_seqs.fasta',
        'density_query_db.fasta',
        'db_sequences.fasta',
        'combined_input.fasta',
        'new_sequences.fasta',
        'combined_aln.fasta',
        'tree_orientation_ref.fasta',
        'id_map.tsv',
    ):
        _move_if_exists(out, name, 'intermediate')
    _move_glob(out, '*_oriented.fasta', 'intermediate')
    _move_glob(out, '*_unique.fasta', 'intermediate')
    _move_glob(out, 'novelty_*_pool.fasta', 'intermediate')
    _move_glob(out, '*.uc', 'intermediate')

    _move_if_exists(out, 'branchmanager.log', 'logs')


def _classification_ids_matching_kingdom(classification_tsv: str, kingdom: str):
    return _classification_ids_matching_kingdom_helper(classification_tsv, kingdom)


def _normalise_sequence_domain(value: str | None) -> str | None:
    if value is None:
        return None
    got = normalize_domain_query(str(value))
    if got in ('', 'none', 'all', 'mixed'):
        return None
    return got


def _sequence_domain_filter(args, *, default: str | None = 'bacteria') -> str | None:
    """Return the domain filter implied by CLI flags."""
    profile = getattr(args, 'sequence_domain', None)
    if profile:
        return _normalise_sequence_domain(profile)
    return _normalise_sequence_domain(default)


def _sequence_domain_label(domain_filter: str | None) -> str:
    return domain_filter if domain_filter else 'mixed'


def _prune_dataset_by_kingdom(db: Database, dataset: str, kingdom: str | None, log_prefix: str):
    deleted = _prune_dataset_by_kingdom_helper(db, dataset, kingdom)
    if deleted:
        logging.getLogger(__name__).info(
            "%s Removed %d sequences from dataset %s not matching kingdom %s",
            log_prefix,
            deleted,
            dataset,
            kingdom,
        )
    return deleted


def _load_evaluate_baseline(
    args,
    db: Database,
    outdir: str,
    effective_ref: str | None,
    effective_taxa_tsv: str | None,
    threads: int,
):
    """Load an optional baseline FASTA before evaluating the current run."""
    baseline_fasta = getattr(args, 'baseline_fasta', None)
    if not baseline_fasta:
        return None

    baseline_dataset = getattr(args, 'baseline_dataset', None) or 'Baseline'
    run_dataset = getattr(args, 'dataset', None)
    if run_dataset and str(baseline_dataset) == str(run_dataset):
        raise SystemExit(
            "[BASELINE] --baseline-dataset must differ from --dataset so novelty can exclude the current run correctly."
        )

    log = logging.getLogger(__name__)
    baseline_out = Path(outdir) / 'baseline_preload'
    baseline_out.mkdir(parents=True, exist_ok=True)

    log.info(
        "[BASELINE] Loading baseline FASTA %s into dataset=%s before evaluating %s",
        baseline_fasta,
        baseline_dataset,
        run_dataset or '(current run)',
    )
    alias_entries, mapped_fasta = db.preload_from_files(
        baseline_fasta,
        taxa_tsv=None,
        color_csv=getattr(args, 'baseline_colors', None),
        source='baseline',
        dataset=baseline_dataset,
        outdir=str(baseline_out),
        shorten_ids=bool(getattr(args, 'baseline_shorten_ids', False)),
    )

    try:
        if alias_entries:
            map_path = _write_id_map_tsv(baseline_out / 'baseline_id_map.tsv', alias_entries)
            log.info("[BASELINE] Wrote baseline id mapping to %s", map_path)
    except Exception as e:
        log.warning("[BASELINE] Could not write baseline id mapping file: %s", e)

    orig_to_short = _build_orig_to_short(alias_entries)
    for short, _orig in alias_entries or []:
        orig_to_short[short] = short

    baseline_ref, baseline_taxa_tsv, baseline_assignment_tsv = _resolve_reference_inputs(
        effective_ref,
        effective_taxa_tsv,
        getattr(args, 'baseline_taxa_assignments', None),
        source_fasta_path=baseline_fasta,
        log_prefix='[BASELINE]',
    )

    if baseline_assignment_tsv:
        log.info("[BASELINE] Using baseline taxonomy assignments from %s", baseline_assignment_tsv)
        try:
            tax_entries = load_taxonomy_entries_from_assignments(
                baseline_assignment_tsv,
                orig_to_short,
                db,
                baseline_dataset,
                source_fasta_path=baseline_fasta,
            )
        except Exception as e:
            log.warning("[BASELINE] Failed to read baseline taxonomy assignments %s: %s", baseline_assignment_tsv, e)
            tax_entries = []
        if tax_entries:
            db.insert_taxonomy(tax_entries)
            log.info("[BASELINE] Inserted/updated taxonomy for %d baseline ids", len(tax_entries))
    elif not bool(getattr(args, 'baseline_skip_classify', False)):
        if not baseline_ref:
            log.warning("[BASELINE] No reference FASTA available; baseline loaded without taxonomy classification")
        else:
            input_for_classify = str(mapped_fasta) if mapped_fasta else baseline_fasta
            log.info("[BASELINE] Classifying baseline %s against %s", input_for_classify, baseline_ref)
            try:
                class_out = classify.run_classification(
                    input_for_classify,
                    str(baseline_out),
                    ref_fasta=baseline_ref,
                    taxa_tsv=baseline_taxa_tsv,
                    threads=threads,
                )
                tax_entries, dist_entries = load_classification_results_for_dataset(
                    class_out,
                    orig_to_short,
                    db,
                    baseline_dataset,
                )
                if tax_entries:
                    db.insert_taxonomy(tax_entries)
                    log.info("[BASELINE] Inserted/updated taxonomy for %d baseline ids", len(tax_entries))
                if dist_entries:
                    db.insert_distances(dist_entries)
                    log.info("[BASELINE] Inserted/updated nearest-reference distances for %d baseline ids", len(dist_entries))
            except Exception as e:
                log.warning("[BASELINE] Baseline classification failed: %s", e)
    else:
        log.info("[BASELINE] --baseline-skip-classify set; baseline loaded without taxonomy classification")

    baseline_domain = _sequence_domain_filter(args)
    if baseline_domain:
        try:
            _prune_dataset_by_kingdom(db, baseline_dataset, baseline_domain, '[BASELINE]')
        except Exception as e:
            log.warning("[BASELINE] Domain-based pruning failed: %s", e)

    try:
        combined_tax = baseline_out / 'baseline_combined_taxonomy.tsv'
        with db.connect() as conn:
            cur = conn.cursor()
            cur.execute(
                "SELECT s.id, t.taxonomy, t.confidence "
                "FROM sequences s LEFT JOIN taxonomy t ON s.id = t.id "
                "WHERE s.dataset = ?",
                (baseline_dataset,),
            )
            rows = cur.fetchall()
        write_combined_taxonomy_tsv(combined_tax, rows)
        log.info("[BASELINE] Wrote baseline combined taxonomy for %d ids to %s", len(rows), combined_tax)
    except Exception as e:
        log.warning("[BASELINE] Failed to write baseline combined taxonomy: %s", e)

    return str(baseline_out)


def _same_path(path_a: str | None, path_b: str | None) -> bool:
    if not path_a or not path_b:
        return False
    try:
        return os.path.samefile(path_a, path_b)
    except Exception:
        try:
            return Path(path_a).resolve() == Path(path_b).resolve()
        except Exception:
            return str(path_a) == str(path_b)


def _build_alt_databases(args) -> list:
    """Return a list of (ref_fasta, taxa_tsv_or_None, db_name) for alt references."""
    alt_refs = getattr(args, 'alt_ref', None) or []
    alt_taxa = getattr(args, 'alt_taxa', None) or []
    alt_names = getattr(args, 'alt_ref_name', None) or []
    result = []
    for i, ref in enumerate(alt_refs):
        taxa = alt_taxa[i] if i < len(alt_taxa) else None
        raw_name = alt_names[i] if i < len(alt_names) else None
        name = raw_name if raw_name else _classify_derive_db_name(ref)
        result.append((ref, taxa, name))
    return result


def _store_alt_taxonomy_in_db(db: Database, all_results: dict, main_db_name: str):
    """Persist alt-db classification results to taxonomy_alt table.

    Iterates over *all_results* ``{db_name: {qid: (hit, pct, tax, conf)}}``,
    skipping the primary database (its results go into the standard taxonomy table).
    """
    log = logging.getLogger(__name__)
    for db_name, results in all_results.items():
        if db_name == main_db_name:
            continue
        alt_entries = []
        for qid, row in results.items():
            if not isinstance(row, (tuple, list)) or len(row) < 4:
                continue
            hit, pct, tax, conf = row[:4]
            if tax == 'NA' and hit == 'NA':
                continue  # skip unclassified rows to keep table lean
            alt_entries.append((qid, db_name, tax, conf, hit, pct))
        if alt_entries:
            try:
                db.insert_taxonomy_alt(alt_entries)
                log.info("[DB] Stored %d alt-db taxonomy entries for ref_db=%s", len(alt_entries), db_name)
            except Exception as e:
                log.warning("[DB] Failed to store alt-db taxonomy for %s: %s", db_name, e)


def _resolve_reference_inputs(
    ref_fasta: str | None,
    taxa_tsv: str | None,
    taxa_assignments: str | None,
    *,
    source_fasta_path: str | None,
    log_prefix: str,
):
    """Resolve whether --taxa-assignments is a TSV assignment file or a reference FASTA.

    External FASTA/FASTA.gz files provided via --taxa-assignments are treated as
    the effective classifier reference database. TSV inputs (and the special case
    where the assignments file is the same file as the source FASTA) are treated
    as direct input-sequence taxonomy assignments.
    """
    effective_ref = ref_fasta
    effective_taxa = taxa_tsv
    assignments_tsv = None

    if not taxa_assignments:
        return effective_ref, effective_taxa, assignments_tsv

    use_as_reference = False
    try:
        use_as_reference = (
            not _same_path(taxa_assignments, source_fasta_path)
            and _assignment_source_is_fasta(taxa_assignments, source_fasta_path=source_fasta_path)
        )
    except Exception:
        use_as_reference = False

    if use_as_reference:
        log = logging.getLogger(__name__)
        if effective_ref and not _same_path(effective_ref, taxa_assignments):
            log.info(
                "%s Treating --taxa-assignments=%s as the GTDB/reference FASTA and using it instead of --ref=%s (--ref is the preferred flag for reference databases)",
                log_prefix,
                taxa_assignments,
                effective_ref,
            )
        else:
            log.info(
                "%s Treating --taxa-assignments=%s as the GTDB/reference FASTA for classification (--ref is the preferred flag for reference databases)",
                log_prefix,
                taxa_assignments,
            )
        if not effective_taxa:
            log.info(
                "%s No --taxa TSV/CSV provided; taxonomy will be parsed directly from reference FASTA headers",
                log_prefix,
            )
        effective_ref = taxa_assignments
        return effective_ref, effective_taxa, None

    assignments_tsv = taxa_assignments
    return effective_ref, effective_taxa, assignments_tsv


def cmd_preload(args):
    db = Database(args.db)
    db.initialise()
    outdir = args.out or '.'
    threads = int(getattr(args, 'threads', 4) or 4)
    requested_kingdom = _sequence_domain_filter(args)
    os.makedirs(outdir, exist_ok=True)
    _configure_logging(outdir)
    logging.getLogger(__name__).info("[PRELOAD] Starting preload into %s", args.db)
    logging.getLogger(__name__).info(
        "[PRELOAD] Sequence-domain profile: %s",
        _sequence_domain_label(requested_kingdom),
    )

    alias_entries, mapped_fasta = db.preload_from_files(
        args.fasta,
        taxa_tsv=getattr(args, 'taxa', None),
        color_csv=getattr(args, 'colors', None),
        source='preload',
        dataset=getattr(args, 'dataset', 'preload'),
        outdir=outdir,
        shorten_ids=bool(getattr(args, 'shorten_ids', False)),
    )
    try:
        if alias_entries:
            preload_map_path = _write_id_map_tsv(Path(outdir) / 'preload_id_map.tsv', alias_entries)
            logging.getLogger(__name__).info('[PRELOAD] Wrote preload id mapping to %s', preload_map_path)
    except Exception as e:
        logging.getLogger(__name__).warning('[PRELOAD] Could not write preload id mapping file: %s', e)
    effective_ref, effective_taxa_tsv, assignment_tsv = _resolve_reference_inputs(
        getattr(args, 'ref', None),
        getattr(args, 'taxa', None),
        getattr(args, 'taxa_assignments', None),
        source_fasta_path=args.fasta,
        log_prefix='[PRELOAD]',
    )
    classification_requested = bool(getattr(args, 'classify', False) or (getattr(args, 'taxa_assignments', None) and not assignment_tsv))

    # If the user provided a table of predetermined taxa assignments, use it
    # instead of running the classifier. The table should have at least ID and
    # taxonomy columns; TSV, CSV, and .gz variants are supported.
    if assignment_tsv:
        taxa_file = assignment_tsv
        logging.getLogger(__name__).info("[PRELOAD] Using taxa assignments from %s (skipping classifier)", taxa_file)
        # build mapping orig->short from alias_entries
        orig_to_short = _build_orig_to_short(alias_entries)
        try:
            tax_entries = load_taxonomy_entries_from_assignments(
                taxa_file,
                orig_to_short,
                db,
                getattr(args, 'dataset', 'preload'),
                source_fasta_path=args.fasta,
            )
        except Exception as e:
            logging.getLogger(__name__).warning("[PRELOAD] Failed to read taxa assignments file %s: %s", taxa_file, e)
            tax_entries = []

        if tax_entries:
            db.insert_taxonomy(tax_entries)
            logging.getLogger(__name__).info("[PRELOAD] Inserted/updated taxonomy for %d preloaded ids from taxa_assignments", len(tax_entries))

    # If classification requested, run classifier on the mapped fasta (short ids)
    # unless the user supplied a taxa assignments table, in which case use that
    # instead and skip running the external classifier.
    if classification_requested and not assignment_tsv:
        if not effective_ref:
            logging.getLogger(__name__).info("[PRELOAD] Classification requested but no reference FASTA was provided via --ref or --taxa-assignments; skipping classification")
        else:
            input_for_classify = str(mapped_fasta) if mapped_fasta else args.fasta
            logging.getLogger(__name__).info("[PRELOAD] Classifying preloaded fasta %s against %s", input_for_classify, effective_ref)

            alt_databases = _build_alt_databases(args)
            ref_name = getattr(args, 'ref_name', None) or _classify_derive_db_name(effective_ref)
            main_db = getattr(args, 'main_ref', None) or ref_name
            all_results: dict = {}

            if alt_databases:
                logging.getLogger(__name__).info(
                    "[PRELOAD] Multi-database classification: primary=%s, alt=%s, main=%s",
                    ref_name, [n for _, _, n in alt_databases], main_db,
                )
                class_out, all_results = classify.run_all_classifications(
                    input_for_classify, outdir,
                    primary_ref=effective_ref,
                    primary_taxa=effective_taxa_tsv,
                    primary_name=ref_name,
                    alt_refs=alt_databases,
                    threads=threads,
                    main_db=main_db,
                )
                _store_alt_taxonomy_in_db(db, all_results, main_db)
            else:
                class_out = classify.run_classification(input_for_classify, outdir, ref_fasta=effective_ref, taxa_tsv=effective_taxa_tsv, threads=threads)

            # parse classification output and persist taxonomy/distances only for
            # the preloaded sequence ids (short ids present in the DB)
            orig_to_short = _build_orig_to_short(alias_entries)
            try:
                tax_entries, dist_entries = load_classification_results_for_dataset(
                    class_out,
                    orig_to_short,
                    db,
                    getattr(args, 'dataset', 'preload'),
                )

                if tax_entries:
                    db.insert_taxonomy(tax_entries)
                    logging.getLogger(__name__).info("[PRELOAD] Inserted/updated taxonomy for %d preloaded ids", len(tax_entries))
                if dist_entries:
                    db.insert_distances(dist_entries)
                    logging.getLogger(__name__).info("[PRELOAD] Inserted/updated distances for %d preloaded ids", len(dist_entries))
            except Exception as e:
                logging.getLogger(__name__).warning("[PRELOAD] Failed to parse classification output: %s", e)

    # Apply the active sequence-domain profile before tree building so the
    # backbone tree is domain-specific.
    try:
        _prune_dataset_by_kingdom(
            db,
            getattr(args, 'dataset', 'preload'),
            requested_kingdom,
            '[PRELOAD]',
        )
    except Exception as e:
        logging.getLogger(__name__).warning("[PRELOAD] Domain-based pruning failed: %s", e)

    # Optionally build a baseline tree/alignment from the preloaded sequences
    if getattr(args, 'build_tree', False):
        try:
            # prefer mapped fasta (short ids) if present
            user_fasta = str(mapped_fasta) if mapped_fasta else args.fasta
            if requested_kingdom:
                try:
                    allowed_ids = _dataset_sequence_ids(db, getattr(args, 'dataset', 'preload'))
                    domain_fasta = Path(outdir) / f'preload_{requested_kingdom}_seqs.fasta'
                    kept = _filter_fasta_to_ids(user_fasta, domain_fasta, allowed_ids)
                    if kept:
                        logging.getLogger(__name__).info(
                            "[PRELOAD] Filtered preload tree FASTA to %d %s sequence(s): %s",
                            kept,
                            requested_kingdom,
                            domain_fasta,
                        )
                        user_fasta = str(domain_fasta)
                    else:
                        logging.getLogger(__name__).warning(
                            "[PRELOAD] Domain filter %s left no FASTA records for tree building; using unfiltered FASTA",
                            requested_kingdom,
                        )
                except Exception as e:
                    logging.getLogger(__name__).warning("[PRELOAD] Failed to filter tree FASTA by domain: %s", e)

            # optionally collapse preloaded sequences before building tree
            if getattr(args, 'collapse', False):
                # require classification to be run for safe taxon-based collapsing
                if not (classification_requested or assignment_tsv):
                    logging.getLogger(__name__).warning("[PRELOAD COLLAPSE] --collapse requires --classify to be set for safe taxon grouping; skipping collapse")
                else:
                    try:
                        threshold = float(getattr(args, 'collapse_threshold', 99.9))
                    except Exception:
                        threshold = 99.9

                    # build qid -> tax mapping from DB taxonomy for this dataset
                    qid_to_tax = {}
                    try:
                        with db.connect() as conn:
                            cur = conn.cursor()
                            cur.execute("SELECT s.id, t.taxonomy FROM sequences s LEFT JOIN taxonomy t ON s.id = t.id WHERE s.dataset = ?", (getattr(args, 'dataset', 'preload'),))
                            for rid, tax in cur.fetchall():
                                qid_to_tax[rid] = tax
                    except Exception:
                        qid_to_tax = {}

                    # group records by tax and cluster within groups
                    from branchmanager.utils.fasta import read_fasta
                    taxa_groups = {}
                    try:
                        for h, s in read_fasta(user_fasta):
                            tax = qid_to_tax.get(h)
                            taxa_groups.setdefault(tax, []).append((h, s))
                    except Exception:
                        taxa_groups = {None: []}
                    try:
                        artifacts = collapse_fasta_within_taxa(
                            taxa_groups,
                            outdir,
                            'preload_short_collapsed.fasta',
                            'preload_collapsed_map.tsv',
                            'preload_collapsed_members.tsv',
                            threshold=threshold,
                            threads=threads,
                            log_prefix='[PRELOAD COLLAPSE]',
                        )
                        user_fasta = artifacts.collapsed_path
                        logging.getLogger(__name__).info("[PRELOAD COLLAPSE] Wrote collapsed preload fasta %s (reps=%d)", user_fasta, len(artifacts.collapsed_records))
                    except Exception as e:
                        logging.getLogger(__name__).warning("[PRELOAD COLLAPSE] Failed to write collapsed preload fasta: %s", e)

            logging.getLogger(__name__).info("[PRELOAD] Building baseline tree/alignment in %s from preloaded sequences", outdir)
            tree.initialise_or_update_tree(
                ref_fasta=effective_ref,
                user_fasta=user_fasta,
                outdir=outdir,
                db=None,
                threads=threads,
                anchor_file=getattr(args, 'anchors', None),
                tree_method=getattr(args, 'tree_method', 'fasttree'),
            )
            logging.getLogger(__name__).info("[PRELOAD] Baseline tree/alignment written to %s", outdir)
            try:
                warning_rows = tree.collect_tree_build_warnings(user_fasta=str(user_fasta), anchor_file=getattr(args, 'anchors', None), db=None)
                warning_rows.extend(tree.summarize_alignment_quality(str(Path(outdir) / 'current_alignment.fasta')))
                if warning_rows:
                    warn_path = tree.write_tree_warning_tsv(outdir, warning_rows)
                    logging.getLogger(__name__).warning("[PRELOAD] Tree/alignment warnings written to %s", warn_path)
            except Exception as e:
                logging.getLogger(__name__).warning("[PRELOAD] Failed to summarise tree/alignment quality: %s", e)
        except Exception as e:
            logging.getLogger(__name__).warning("[PRELOAD] Failed to build baseline tree: %s", e)

    # If kingdom filter requested for preload, remove any preloaded sequences
    # whose assigned taxonomy indicates they are not the requested kingdom.
    try:
        _prune_dataset_by_kingdom(
            db,
            getattr(args, 'dataset', 'preload'),
            requested_kingdom,
            '[PRELOAD]',
        )
    except Exception as e:
        logging.getLogger(__name__).warning("[PRELOAD] Kingdom-based pruning failed: %s", e)

    # Build combined taxonomy and generate iTOL color files for the preload dataset.
    try:
        out_p = Path(outdir)
        combined_tax = out_p / 'preload_combined_taxonomy.tsv'
        if not combined_tax.exists():
            try:
                with db.connect() as conn:
                    cur = conn.cursor()
                    cur.execute("SELECT s.id, t.taxonomy, t.confidence FROM sequences s LEFT JOIN taxonomy t ON s.id = t.id WHERE s.dataset = ?", (getattr(args, 'dataset', 'preload'),))
                    rows = cur.fetchall()
                write_combined_taxonomy_tsv(combined_tax, rows)
                logging.getLogger(__name__).info("[PRELOAD] Wrote combined taxonomy for %d ids to %s", len(rows), combined_tax)
            except Exception as e:
                logging.getLogger(__name__).warning("[PRELOAD] Failed to build combined taxonomy: %s", e)

        # Build id_map from alias_entries so the generator can emit short ids
        id_map = {}
        try:
            if alias_entries:
                for short, orig in alias_entries:
                    id_map[orig] = short
                    id_map[short] = short
                    try:
                        cid = db._canonical_from_header(orig)
                        if cid:
                            id_map[cid] = short
                    except Exception:
                        pass
        except Exception:
            id_map = {}

        # Call generator (fallback if id_map not accepted)
        try:
            tree_path = Path(outdir) / 'current_tree.nwk'
            tfile = str(tree_path) if tree_path.exists() else _find_tree_file_in_dir(outdir)
            itol.generate_itol_colors(str(combined_tax), outdir, user_color_csv=getattr(args, 'colors', None), id_map=id_map, tree_file=tfile, phylum_groups=getattr(args, 'group_phyla', None))
            logging.getLogger(__name__).info("[PRELOAD] Generated iTOL color files in %s", outdir)
        except TypeError:
            try:
                itol.generate_itol_colors(str(combined_tax), outdir, user_color_csv=getattr(args, 'colors', None), phylum_groups=getattr(args, 'group_phyla', None))
                logging.getLogger(__name__).info("[PRELOAD] Generated iTOL color files in %s (no id_map)", outdir)
            except Exception as e:
                logging.getLogger(__name__).warning("[PRELOAD] Failed to generate iTOL colors: %s", e)
        except Exception as e:
            logging.getLogger(__name__).warning("[PRELOAD] Failed to generate iTOL colors: %s", e)

        # ── Functional annotations (optional) ─────────────────────────────────
        try:
            func_tsv = getattr(args, 'functional', None)
            if func_tsv:
                try:
                    written = itol.write_functional_annotations(str(func_tsv), outdir, id_map=id_map)
                    logging.getLogger(__name__).info("[PRELOAD] Wrote functional annotation iTOL files: %s", ','.join(written) if written else '(none)')
                except Exception as e:
                    logging.getLogger(__name__).warning("[PRELOAD] Functional annotations generation failed: %s", e)
        except Exception:
            pass

        # ── Draft rumen functional groups (optional) ──────────────────────────
        if getattr(args, 'draft_rumen_functions', False):
            try:
                combined_tax_path = str(out_p / 'combined_taxonomy.tsv')
                tsv_out, itol_out = itol.generate_rumen_function_draft(
                    combined_tax_path, outdir, id_map=id_map
                )
                if tsv_out:
                    logging.getLogger(__name__).info(
                        "[PRELOAD] Draft rumen functional annotation: %s", tsv_out
                    )
                if itol_out:
                    logging.getLogger(__name__).info(
                        "[PRELOAD] Rumen functional iTOL file: %s", itol_out
                    )
            except Exception as e:
                logging.getLogger(__name__).warning(
                    "[PRELOAD] Draft rumen functions generation failed: %s", e
                )

        # Write a simple dataset colorstrip mapping preloaded ids to dataset color
        try:
            ds_color = getattr(args, 'dataset_color', None) if hasattr(args, 'dataset_color') else None
            if not ds_color:
                # deterministic dataset-level color (use distinct palette)
                ds_color = itol._name_to_dataset_color(getattr(args, 'dataset', 'preload'))
            itol_path = out_p / 'itol_dataset_preload.itol'
            dataset_label = getattr(args, 'dataset', 'preload')
            with db.connect() as conn:
                cur = conn.cursor()
                cur.execute("SELECT id FROM sequences WHERE dataset = ?", (getattr(args, 'dataset', 'preload'),))
                id_to_color = {iid: ds_color for (iid,) in cur.fetchall()}
            itol.write_dataset_colorstrip(str(itol_path), dataset_label, id_to_color, legend_title=f"{dataset_label} legend")
            logging.getLogger(__name__).info("[PRELOAD] Wrote dataset ITOL colorstrip to %s", itol_path)
        except Exception as e:
            logging.getLogger(__name__).warning("[PRELOAD] Failed to write dataset iTOL colorstrip: %s", e)
    except Exception:
        pass

    # write brief explanations for files produced by preload
    try:
        _write_output_explanations(outdir)
    except Exception:
        pass


def cmd_run(args):
    db = Database(args.db)
    db.initialise()
    outdir = args.out
    threads = int(getattr(args, 'threads', 4) or 4)
    os.makedirs(outdir, exist_ok=True)
    _configure_logging(outdir)
    logging.getLogger(__name__).info("[RUN] Starting run pipeline (input=%s)", args.input)
    domain_filter = _sequence_domain_filter(args)
    logging.getLogger(__name__).info(
        "[RUN] Sequence-domain profile: %s",
        _sequence_domain_label(domain_filter),
    )
    effective_ref, effective_taxa_tsv, assignment_tsv = _resolve_reference_inputs(
        getattr(args, 'ref', None),
        getattr(args, 'taxa', None),
        getattr(args, 'taxa_assignments', None),
        source_fasta_path=args.input,
        log_prefix='[RUN]',
    )
    if not effective_ref:
        raise SystemExit("[RUN] A reference FASTA is required via --ref, or `--taxa-assignments` must point to a GTDB/reference FASTA rather than a taxonomy assignment table.")
    if getattr(args, 'command', None) in ('evaluate', 'eval') and not getattr(args, 'partner_metadata', None):
        raise SystemExit(
            '[RUN] evaluate requires --partner-metadata / --sequencing-metadata: '
            'a CSV/TSV sidecar table with sequence IDs, partner IDs, and WGS-selected status.'
        )

    _load_evaluate_baseline(args, db, outdir, effective_ref, effective_taxa_tsv, threads)

    # QC
    qc_out = qc.run_qc(args.input, outdir, min_len=getattr(args, 'min_len', 800), max_n=getattr(args, 'max_n', 5))

    # derep
    derep_out = derep.run_derep(qc_out, outdir)

    # map user-provided dereplicated IDs to runtime IDs and insert into DB
    from branchmanager.utils.fasta import read_fasta, write_fasta
    mapped_derep = Path(outdir) / 'derep_short.fasta'
    used_ids = set(db.get_all_ids())
    orig_to_short = {}
    mapped_records = []
    skipped_existing = 0

    # If the user requested a kingdom filter run classification early on the
    # dereplicated fasta so we can keep only sequences assigned to the chosen
    # kingdom. This avoids inserting unwanted sequences into the DB.
    early_class_out = None
    allowed_qids = None
    kingdom = domain_filter
    if kingdom:
        kingdom_text = str(kingdom)
        if not effective_ref:
            logging.getLogger(__name__).warning("[RUN] --sequence-domain specified but no reference FASTA was available; cannot classify to filter; proceeding without domain filtering")
            allowed_qids = None
        else:
            try:
                logging.getLogger(__name__).info("[RUN] Running pre-insert classification on dereplicated fasta to filter by domain=%s", kingdom)
                early_class_out = classify.run_classification(str(derep_out), outdir, ref_fasta=effective_ref, taxa_tsv=effective_taxa_tsv, threads=threads)
                allowed_qids = _classification_ids_matching_kingdom(early_class_out, kingdom_text)
                logging.getLogger(__name__).info("[RUN] Kingdom filter: %d dereplicated sequences match %s", len(allowed_qids), kingdom)
            except Exception as e:
                logging.getLogger(__name__).warning("[RUN] Failed to run pre-insert classification for kingdom filtering: %s", e)
                allowed_qids = None

    for h, s in read_fasta(derep_out):
        # if kingdom filtering is active, skip sequences that did not match
        if allowed_qids is not None:
            # try a few candidate header forms for matching
            candidates = [h, h.split('|')[-1], h.split()[-1] if h.split() else h]
            hit = False
            for c in candidates:
                if c in allowed_qids:
                    hit = True
                    break
            if not hit:
                continue

        # If DB already contains this sequence ID (or another alias), we still
        # want to include it in the tree for visualization, but we use the
        # existing DB ID so downstream taxonomy/metadata resolve consistently.
        # The "INSERT OR IGNORE" during insert_sequences will prevent DB duplicates.
        existing_id = db.resolve_sequence_id(h)
        if existing_id and existing_id in used_ids:
            skipped_existing += 1
            short = existing_id
            orig_to_short[h] = short
            try:
                cid = canonicalize_sequence_id(h)
                if cid:
                    orig_to_short[cid] = short
            except Exception:
                pass
            # Still add to mapped_records so it appears in the tree
            mapped_records.append((short, s))
        else:
            try:
                short = db.choose_effective_sequence_id(
                    h,
                    used_ids,
                    shorten_ids=bool(getattr(args, 'shorten_ids', False)),
                )
            except ValueError as e:
                raise SystemExit(f"[RUN] {e}")
            mapped_records.append((short, s))
            orig_to_short[h] = short
            try:
                cid = canonicalize_sequence_id(h)
                if cid:
                    orig_to_short[cid] = short
            except Exception:
                pass
    if skipped_existing:
        logging.getLogger(__name__).info("[DB] Found %d sequences already in DB; keeping them for tree inclusion", skipped_existing)
    write_fasta(mapped_records, str(mapped_derep))
    logging.getLogger(__name__).info("[DB] Mapped %d user sequence IDs to runtime IDs and wrote %s", len(mapped_records), mapped_derep)
    if not mapped_records:
        logging.getLogger(__name__).warning("[DB] No sequences were mapped — check if input sequences were filtered")
    # insert mapped records into DB (dataset provided by user)
    run_dataset = getattr(args, 'dataset', 'user')
    db.insert_sequences(mapped_records, dataset=run_dataset)
    # write mapping file
    try:
        map_path = Path(outdir) / 'user_id_map.tsv'
        _write_id_map_tsv(map_path, ((short, orig) for orig, short in orig_to_short.items()))
        logging.getLogger(__name__).info("[DB] Wrote user id mapping to %s", map_path)
    except Exception as e:
        logging.getLogger(__name__).warning("[DB] Could not write user id mapping file: %s", e)

    for short, _seq in mapped_records:
        orig_to_short[short] = short

    run_metadata = _load_partner_metadata_for_run(
        args,
        db,
        outdir,
        orig_to_short,
        [short for short, _seq in mapped_records],
    )

    # classify (or use external taxa assignments if provided)
    class_out = None
    all_class_results: dict = {}
    run_main_db_name: str = (
        getattr(args, 'main_ref', None)
        or getattr(args, 'ref_name', None)
        or _classify_derive_db_name(effective_ref)
        or 'main'
    )
    if assignment_tsv:
        taxa_file = assignment_tsv
        logging.getLogger(__name__).info("[RUN] Using taxa assignments from %s (skipping classifier)", taxa_file)
        # parse provided taxonomy assignment table and insert taxonomy for mapped runtime ids
        try:
            tax_entries_local = load_taxonomy_entries_from_assignments(
                taxa_file,
                orig_to_short,
                db,
                run_dataset,
                source_fasta_path=args.input,
            )
        except Exception as e:
            logging.getLogger(__name__).warning("[RUN] Failed to read taxa assignments file %s: %s", taxa_file, e)
            tax_entries_local = []

        # persist taxonomy for mapped ids
        if tax_entries_local:
            try:
                db.insert_taxonomy(tax_entries_local)
                logging.getLogger(__name__).info("[DB] Inserted/updated taxonomy for %d ids from taxa_assignments", len(tax_entries_local))
            except Exception as e:
                logging.getLogger(__name__).warning("[DB] Failed to insert taxonomy from taxa_assignments: %s", e)

        # create a synthetic classification-like file so downstream code that
        # expects `class_out` can run unchanged. Format: qid\tbest\tidentity\ttaxon\tconfidence
        try:
            class_out_path = Path(outdir) / 'taxa_assignments_classout.tsv'
            with open(class_out_path, 'w') as cf:
                cf.write('id\tbest\tidentity\ttaxon\tconfidence\n')
                try:
                    for row in iter_assignment_rows(taxa_file, source_fasta_path=args.input):
                        qid = row.get('qid', 'NA')
                        tax = row.get('tax') if row.get('tax') is not None else 'NA'
                        conf = row.get('confidence') if row.get('confidence') is not None else 'NA'
                        cf.write(f"{qid}\t\tNA\t{tax}\t{conf}\n")
                except Exception:
                    pass
            class_out = str(class_out_path)
        except Exception:
            class_out = None
    else:
        # Only run classification here if we did not already run an early
        # pre-insert classification for kingdom filtering and the user did
        # not supply taxa_assignments.
        if 'early_class_out' in locals() and early_class_out:
            class_out = early_class_out
        else:
            alt_databases = _build_alt_databases(args)
            ref_name = getattr(args, 'ref_name', None) or _classify_derive_db_name(effective_ref)
            run_main_db_name = getattr(args, 'main_ref', None) or ref_name

            if alt_databases:
                logging.getLogger(__name__).info(
                    "[RUN] Multi-database classification: primary=%s, alt=%s, main=%s",
                    ref_name, [n for _, _, n in alt_databases], run_main_db_name,
                )
                class_out, all_class_results = classify.run_all_classifications(
                    str(mapped_derep), outdir,
                    primary_ref=effective_ref,
                    primary_taxa=effective_taxa_tsv,
                    primary_name=ref_name,
                    alt_refs=alt_databases,
                    threads=threads,
                    main_db=run_main_db_name,
                )
                _store_alt_taxonomy_in_db(db, all_class_results, run_main_db_name)
            else:
                class_out = classify.run_classification(str(mapped_derep), outdir, ref_fasta=effective_ref, taxa_tsv=effective_taxa_tsv, threads=threads)

    # novelty
    target_fasta = getattr(args, 'target', None)
    novelty_out = novelty.run_novelty(str(mapped_derep), effective_ref, outdir, db=db, run_dataset=run_dataset, threads=threads, target_fasta=target_fasta)
    try:
        novelty_baseline_datasets = []
        if getattr(args, 'baseline_dataset', None):
            novelty_baseline_datasets.append(getattr(args, 'baseline_dataset'))
        novelty_baseline_datasets.extend(getattr(args, 'novelty_baseline_datasets', None) or [])
        novelty_metrics_out = novelty.build_reference_novelty_metrics(
            str(mapped_derep),
            effective_ref,
            novelty_out,
            outdir,
            threads=threads,
            db=db,
            run_dataset=run_dataset,
            target_fasta=target_fasta,
            baseline_datasets=novelty_baseline_datasets,
        )
        logging.getLogger(__name__).info("[NOVELTY] Wrote novelty metrics to %s", novelty_metrics_out)
    except Exception as e:
        logging.getLogger(__name__).warning("[NOVELTY] Failed to build novelty metrics: %s", e)

    # persist taxonomy & distances for user sequences only
    try:
        tax_entries, dist_entries = load_classification_results_for_dataset(
            class_out or '',
            orig_to_short,
            db,
            run_dataset,
        )
    except Exception as e:
        logging.getLogger(__name__).warning("[RUN] Failed to parse classification output: %s", e)
        tax_entries, dist_entries = [], []

    if tax_entries:
        db.insert_taxonomy(tax_entries)
        logging.getLogger(__name__).info("[DB] Inserted/updated taxonomy for %d ids", len(tax_entries))
    if dist_entries:
        db.insert_distances(dist_entries)
        logging.getLogger(__name__).info("[DB] Inserted/updated distances for %d ids", len(dist_entries))
    # Safety-net: remove any run sequences with explicit non-matching kingdom assignments.
    try:
        _prune_dataset_by_kingdom(db, run_dataset, str(kingdom) if kingdom else None, '[RUN]')
    except Exception as e:
        logging.getLogger(__name__).warning("[RUN] Kingdom-based pruning failed: %s", e)

    # Initialise cluster-tracking variables used by both the tree section and
    # the assessment section below.  They will be populated if --collapse is
    # active; otherwise they stay empty and the assessment omits cluster columns.
    tree_fasta = mapped_derep
    run_member_to_rep: dict = {}
    run_rep_to_members: dict = {}

    # update tree/alignment
    try:
        # report DB sequence count for diagnostics so users can verify that
        # preloaded sequences exist and will be used to seed the backbone
        try:
            db_ids = db.get_all_ids()
            logging.getLogger(__name__).info("[TREE] DB contains %d sequences; these will be considered for backbone construction", len(db_ids))
        except Exception:
            logging.getLogger(__name__).info("[TREE] Could not determine DB sequence count before tree build")

        # Optionally collapse highly similar run sequences (for tree readability)
        if getattr(args, 'collapse', False):
            try:
                threshold = float(getattr(args, 'collapse_threshold', 99.9))
            except Exception:
                threshold = 99.9
            # build qid -> tax mapping from classification output so we only
            # cluster sequences that share the same taxonomic assignment
            qid_to_tax = {}
            try:
                if class_out:
                    with open(class_out) as cf:
                        next(cf, None)
                        for line in cf:
                            parts = line.strip().split('\t')
                            if not parts:
                                continue
                            q = parts[0]
                            tax = parts[3] if len(parts) > 3 else None
                            qid_to_tax[q] = tax
            except Exception:
                qid_to_tax = {}

            # group short ids by tax
            taxa_groups = {}
            try:
                for h, s in read_fasta(str(mapped_derep)):
                    tax = qid_to_tax.get(h)
                    taxa_groups.setdefault(tax, []).append((h, s))
            except Exception:
                taxa_groups = {None: []}
            try:
                artifacts = collapse_fasta_within_taxa(
                    taxa_groups,
                    outdir,
                    'derep_short_collapsed.fasta',
                    'collapsed_map.tsv',
                    'collapsed_members.tsv',
                    threshold=threshold,
                    threads=threads,
                    log_prefix='[COLLAPSE]',
                )
                tree_fasta = artifacts.collapsed_path
                run_member_to_rep = artifacts.member_to_rep or {}
                # Build rep -> members mapping
                for mem, rep in run_member_to_rep.items():
                    run_rep_to_members.setdefault(rep, []).append(mem)
                logging.getLogger(__name__).info("[COLLAPSE] Wrote collapsed fasta %s (reps=%d)", tree_fasta, len(artifacts.collapsed_records))
            except Exception as e:
                logging.getLogger(__name__).warning("[COLLAPSE] Failed to write collapsed fasta: %s", e)

        # Optionally filter by phylum before tree building
        phylum_filter = getattr(args, 'phylum', None)
        if phylum_filter:
            try:
                # Parse classification to extract phylum-level assignments
                phylum_ids = set()
                if class_out:
                    with open(class_out) as cf:
                        next(cf, None)  # skip header
                        for line in cf:
                            parts = line.strip().split('\t')
                            if len(parts) >= 4:
                                qid = parts[0]
                                tax = parts[3]  # Full taxonomy string
                                # Extract phylum from taxonomy (first level after kingdom)
                                # Format is usually: k__Kingdom; p__Phylum; c__Class; ...
                                try:
                                    tax_parts = [t.strip() for t in tax.split(';')]
                                    for part in tax_parts:
                                        if part.startswith('p__'):
                                            phylum_name = part[3:].strip()
                                            # Normalize for comparison (case-insensitive, strip underscores)
                                            if phylum_name.lower().replace('_', '') == phylum_filter.lower().replace('_', ''):
                                                phylum_ids.add(qid)
                                            break
                                except Exception:
                                    pass
                
                # Filter tree_fasta to only include sequences in phylum_ids
                if phylum_ids:
                    filtered_fasta = Path(outdir) / f'tree_sequences_phylum_{phylum_filter}.fasta'
                    filtered_records = []
                    total_checked = 0
                    try:
                        for h, s in read_fasta(str(tree_fasta)):
                            total_checked += 1
                            if h in phylum_ids:
                                filtered_records.append((h, s))
                    except Exception:
                        pass

                    if filtered_records:
                        write_fasta(filtered_records, str(filtered_fasta))
                        logging.getLogger(__name__).info(
                            "[PHYLUM] Filtered tree sequences from %d to %d for phylum: %s",
                            total_checked, len(filtered_records), phylum_filter,
                        )
                        tree_fasta = filtered_fasta
                    else:
                        logging.getLogger(__name__).warning("[PHYLUM] No sequences found for phylum '%s' in classification results", phylum_filter)
                else:
                    logging.getLogger(__name__).warning("[PHYLUM] Could not extract phylum information from classification output")
            except Exception as e:
                logging.getLogger(__name__).warning("[PHYLUM] Phylum filtering failed: %s", e)

        # pass the Database object so the tree builder can export existing
        # preloaded sequences from the DB to form the backbone alignment
        _force_rebuild = bool(getattr(args, 'force_rebuild', False))
        try:
            tree_seq_count = sum(1 for line in open(str(tree_fasta)) if line.startswith(">"))
            logging.getLogger(__name__).info("[TREE] About to build/update tree with %d sequences from: %s", tree_seq_count, tree_fasta)
        except Exception as e:
            logging.getLogger(__name__).warning("[TREE] Could not count sequences in tree_fasta: %s", e)
        if _force_rebuild:
            logging.getLogger(__name__).info("[TREE] Force rebuild requested (--force-rebuild/--rebuild-tree); tree will be rebuilt from scratch including all datasets")
        tree.initialise_or_update_tree(
            ref_fasta=effective_ref,
            user_fasta=str(tree_fasta),
            outdir=outdir,
            db=db,
            db_dataset=None,
            threads=threads,
            preload_dir=None if _force_rebuild else getattr(args, 'preload_dir', None),
            force_rebuild=_force_rebuild,
            anchor_file=getattr(args, 'anchors', None),
            tree_method=getattr(args, 'tree_method', 'fasttree'),
        )
        try:
            warning_rows = tree.collect_tree_build_warnings(user_fasta=str(tree_fasta), anchor_file=getattr(args, 'anchors', None), db=db, db_dataset=None)
            warning_rows.extend(tree.summarize_alignment_quality(str(Path(outdir) / 'current_alignment.fasta')))
            if warning_rows:
                warn_path = tree.write_tree_warning_tsv(outdir, warning_rows)
                logging.getLogger(__name__).warning("[TREE] Tree/alignment warnings written to %s", warn_path)
        except Exception as e:
            logging.getLogger(__name__).warning("[TREE] Failed to summarise tree/alignment quality: %s", e)
    except Exception as e:
        logging.getLogger(__name__).warning("[TREE] Tree update failed: %s", e)

    # Build combined taxonomy for iTOL from DB sequences and current run results.
    combined_path = Path(outdir) / 'combined_taxonomy.tsv'
    merged = {}
    order = []

    try:
        preload_dir = getattr(args, 'preload_dir', None)
        preload_ids = None
        if preload_dir:
            try:
                cand = next(
                    (path for path in _combined_taxonomy_candidates(preload_dir) if path.exists()),
                    None,
                )
                if cand is not None:
                    preload_ids = read_combined_taxonomy_ids(cand)
            except Exception:
                preload_ids = None

        base_rows = collect_db_taxonomy_rows(db, preload_ids if preload_ids else None)
        merged_rows = merge_combined_taxonomy_rows(base_rows, class_out or '', orig_to_short, db)
        for rid, tax, conf in merged_rows:
            merged[rid] = (tax if tax is not None else 'NA', conf if conf is not None else 'NA')
            order.append(rid)
    except Exception:
        # fallback: leave merged empty and continue
        pass

    # write merged combined taxonomy preserving order (DB-order first, then new)
    try:
        write_combined_taxonomy_tsv(
            combined_path,
            [(iid, merged.get(iid, ('NA', 'NA'))[0], merged.get(iid, ('NA', 'NA'))[1]) for iid in order],
        )
        logging.getLogger(__name__).info("[ITOL] Wrote combined taxonomy to %s", combined_path)
    except Exception as e:
        logging.getLogger(__name__).warning("[ITOL] Failed to write combined taxonomy: %s", e)

    # Load preload id_map for iTOL when --preload-dir is provided.
    id_map_for_itol = None
    preload_dir = getattr(args, 'preload_dir', None)
    if preload_dir:
        try:
            id_map = {}
            cand_map = _find_preferred_id_map(str(preload_dir))
            if cand_map:
                id_map = _load_id_map_from_tsv(cand_map, db=db)
            if id_map:
                id_map_for_itol = id_map
        except Exception:
            id_map_for_itol = None

        try:
            tree_path = Path(outdir) / 'current_tree.nwk'
            tfile = str(tree_path) if tree_path.exists() else _find_tree_file_in_dir(outdir)
            itol.generate_itol_colors(str(combined_path), outdir, user_color_csv=getattr(args, 'user_colors', None), id_map=id_map_for_itol, tree_file=tfile, phylum_groups=getattr(args, 'group_phyla', None))
            logging.getLogger(__name__).info("[ITOL] Generated iTOL color files in %s", outdir)
        except Exception as e:
            logging.getLogger(__name__).warning("[ITOL] Failed to generate iTOL files: %s", e)
    # Generate iTOL color files for the run output directory.
    try:
        tree_path = Path(outdir) / 'current_tree.nwk'
        tfile = str(tree_path) if tree_path.exists() else _find_tree_file_in_dir(outdir)
        itol.generate_itol_colors(str(combined_path), outdir, user_color_csv=getattr(args, 'user_colors', None), id_map=id_map_for_itol, tree_file=tfile, phylum_groups=getattr(args, 'group_phyla', None))
        logging.getLogger(__name__).info("[ITOL] Generated iTOL color files in %s", outdir)
    except Exception as e:
        logging.getLogger(__name__).warning("[ITOL] Failed to generate iTOL files: %s", e)
    # Optional: write functional annotation iTOL datasets when provided
    try:
        func_tsv = getattr(args, 'functional', None)
        if func_tsv:
            try:
                written = itol.write_functional_annotations(str(func_tsv), outdir, id_map=id_map_for_itol)
                logging.getLogger(__name__).info("[ITOL] Wrote functional annotation iTOL files: %s", ','.join(written) if written else '(none)')
            except Exception as e:
                logging.getLogger(__name__).warning("[ITOL] Functional annotations generation failed: %s", e)
    except Exception:
        pass
    # Draft rumen functional groups (auto-generated from output taxonomy)
    if getattr(args, 'draft_rumen_functions', False):
        try:
            tsv_out, itol_out = itol.generate_rumen_function_draft(
                str(combined_path), outdir, id_map=id_map_for_itol
            )
            if tsv_out:
                logging.getLogger(__name__).info("[ITOL] Draft rumen functional annotation: %s", tsv_out)
            if itol_out:
                logging.getLogger(__name__).info("[ITOL] Rumen functional iTOL file: %s", itol_out)
        except Exception as e:
            logging.getLogger(__name__).warning("[ITOL] Draft rumen functions generation failed: %s", e)
    # produce dataset membership band (preload vs run)
    try:
        combined_tax = combined_path
        ids_in_order = []
        if combined_tax.exists():
            with open(combined_tax) as ct:
                next(ct, None)
                for l in ct:
                    iid = l.strip().split('\t')[0]
                    if iid:
                        ids_in_order.append(iid)
        # Query DB for exact dataset membership per id to support arbitrary dataset names
        ds_map = {}
        try:
            with db.connect() as conn:
                cur = conn.cursor()
                placeholders = ','.join('?' for _ in ids_in_order) if ids_in_order else ''
                if placeholders:
                    cur.execute(f"SELECT id, dataset FROM sequences WHERE id IN ({placeholders})", tuple(ids_in_order))
                    for iid_row, ds in cur.fetchall():
                        ds_map[iid_row] = ds or ''
        except Exception:
            ds_map = {}

        membership_path = Path(outdir) / 'itol_dataset_membership.itol'
        itol.write_dataset_membership_strip(str(membership_path), ids_in_order, ds_map)
        logging.getLogger(__name__).info("[ITOL] Wrote dataset membership ITOL to %s", membership_path)
    except Exception as e:
        logging.getLogger(__name__).warning("[ITOL] Failed to write dataset membership ITOL: %s", e)

    # Also produce a novelty gradient colorstrip for NEW (run) sequences
    try:
        try:
            run_ids = [r[0] for r in mapped_records]
        except Exception:
            run_ids = []
        if not run_ids:
            try:
                with db.connect() as conn:
                    cur = conn.cursor()
                    cur.execute("SELECT id FROM sequences WHERE dataset = ?", (run_dataset,))
                    run_ids = [r[0] for r in cur.fetchall()]
            except Exception:
                run_ids = []

        if run_ids:
            novelty.build_run_novelty_itol(
                outdir,
                run_ids,
                str(mapped_derep),
                db,
                run_dataset,
                orig_to_short,
                threads=getattr(args, 'threads', None),
            )
    except Exception as e:
        logging.getLogger(__name__).warning("[ITOL] Failed to write novelty ITOL: %s", e)

    # placement warning summary combining classification + novelty support
    try:
        if class_out:
            warning_rows = build_placement_warning_rows(class_out, novelty_out, orig_to_short, db)
            if warning_rows:
                warn_path = write_placement_warning_tsv(Path(outdir) / 'placement_warnings.tsv', warning_rows)
                logging.getLogger(__name__).warning("[RUN] Placement warnings written to %s", warn_path)
            try:
                run_ids_for_assessment = [r[0] for r in mapped_records]
                # Determine which sequences are actually in the tree
                try:
                    _tree_ids = set()
                    for _h, _s in read_fasta(str(tree_fasta)):
                        _tree_ids.add(_h)
                except Exception:
                    _tree_ids = None

                # Collect alt-db taxonomy from DB for the run sequences
                alt_ref_dbs: list = []
                alt_taxonomies: dict = {}
                try:
                    alt_ref_dbs = db.get_alt_ref_dbs()
                    if alt_ref_dbs:
                        alt_taxonomies = db.get_taxonomy_alt_for_ids(run_ids_for_assessment)
                        logging.getLogger(__name__).info(
                            "[RUN] Found alt-db taxonomy for ref_dbs: %s", alt_ref_dbs
                        )
                except Exception as e:
                    logging.getLogger(__name__).warning("[RUN] Could not load alt-db taxonomy: %s", e)

                assessment_rows = build_sequence_assessment_rows(
                    run_ids_for_assessment,
                    class_out,
                    str(Path(outdir) / 'novelty_metrics.tsv'),
                    warning_rows,
                    orig_to_short,
                    db,
                    member_to_rep=run_member_to_rep if run_member_to_rep else None,
                    rep_to_members=run_rep_to_members if run_rep_to_members else None,
                    tree_ids=_tree_ids,
                    alt_taxonomies=alt_taxonomies,
                    alt_ref_dbs=alt_ref_dbs,
                )
                mwl_path = getattr(args, 'mwl', None)
                mwl_entries = []
                if mwl_path:
                    try:
                        logging.getLogger(__name__).info(
                            "[MWL] Matching the primary Taxonomy column against MWL. For multi-db runs, set --main-ref to your GTDB reference name."
                        )
                        mwl_entries = _mwl.load_mwl_entries(
                            str(mwl_path),
                            sheet_name=getattr(args, 'mwl_sheet', 'MWL_V1') or 'MWL_V1',
                        )
                        _mwl.annotate_assessment_rows(
                            assessment_rows,
                            mwl_entries,
                            min_rank=getattr(args, 'mwl_min_rank', 'p') or 'p',
                        )
                        mwl_report = _mwl.write_mwl_matches_tsv(
                            Path(outdir) / 'mwl_matches.tsv',
                            assessment_rows,
                        )
                        logging.getLogger(__name__).info(
                            "[MWL] Annotated %d assessment rows with %d MWL entries → %s",
                            len(assessment_rows), len(mwl_entries), mwl_report,
                        )
                    except Exception as e:
                        logging.getLogger(__name__).warning("[MWL] Failed to annotate assessment rows: %s", e)
                assess_path = write_sequence_assessment_tsv(Path(outdir) / 'sequence_assessment.tsv', assessment_rows)
                selection_summary_path = write_selection_summary_tsv(Path(outdir) / 'selection_summary.tsv', assessment_rows)
                baseline_hits_path = write_baseline_hits_tsv(Path(outdir) / 'baseline_hits.tsv', assessment_rows)
                logging.getLogger(__name__).info("[RUN] Wrote sequence assessment to %s", assess_path)
                logging.getLogger(__name__).info("[RUN] Wrote SAB selection summary to %s", selection_summary_path)
                logging.getLogger(__name__).info("[RUN] Wrote nearest baseline hit report to %s", baseline_hits_path)

                # ── Cluster-level reports + phylogenetic isolation ───────────
                try:
                    _tree_nwk = str(Path(outdir) / 'current_tree.nwk')
                    _cluster_summary, _cluster_csvs, _backup_tsv = _cluster_report.generate_cluster_reports(
                        outdir=outdir,
                        assessment_rows=assessment_rows,
                        tree_path=_tree_nwk if Path(_tree_nwk).exists() else None,
                    )
                    # Re-write sequence_assessment.tsv now that phylo_isolation /
                    # investigation_score have been filled in by generate_cluster_reports
                    if mwl_entries:
                        try:
                            _mwl.add_evaluation_scores(assessment_rows)
                            _mwl.write_mwl_matches_tsv(Path(outdir) / 'mwl_matches.tsv', assessment_rows)
                        except Exception as _mwe:
                            logging.getLogger(__name__).warning("[MWL] Failed to refresh MWL evaluation scores: %s", _mwe)
                    write_sequence_assessment_tsv(Path(outdir) / 'sequence_assessment.tsv', assessment_rows)
                    write_selection_summary_tsv(Path(outdir) / 'selection_summary.tsv', assessment_rows)
                    write_baseline_hits_tsv(Path(outdir) / 'baseline_hits.tsv', assessment_rows)
                    if _cluster_summary:
                        logging.getLogger(__name__).info(
                            "[CLUSTER] Wrote cluster summary → %s  (%d per-cluster CSVs in %s/clusters/)",
                            _cluster_summary, len(_cluster_csvs), outdir,
                        )
                    if _backup_tsv:
                        logging.getLogger(__name__).info(
                            "[CLUSTER] Wrote backup candidates table → %s", _backup_tsv,
                        )
                except Exception as _ce:
                    logging.getLogger(__name__).warning("[CLUSTER] Cluster report generation failed: %s", _ce)
                # Emit a user-friendly summary of HIGH priority candidates
                try:
                    def _effective_priority(row):
                        return row.get('genome_sequencing_adjusted_priority') or row.get('sequencing_priority')

                    high_priority = [r for r in assessment_rows if _effective_priority(r) == 'HIGH']
                    medium_priority = [r for r in assessment_rows if _effective_priority(r) == 'MEDIUM']
                    selected_clades = [
                        r for r in assessment_rows
                        if r.get('clade_already_selected_for_genome_sequencing') == 'True'
                    ]
                    collapsed_away = [r for r in assessment_rows if r.get('in_tree') == 'No']
                    logging.getLogger(__name__).info(
                        "[ASSESSMENT SUMMARY] %d sequences assessed: %d HIGH priority for sequencing, "
                        "%d MEDIUM priority. %d have a >=97%% selected-genome neighbour. "
                        "%d sequences were clustered and excluded from tree "
                        "(their representatives are in the tree). "
                        "See assessment/sequence_assessment.tsv for full details.",
                        len(assessment_rows),
                        len(high_priority),
                        len(medium_priority),
                        len(selected_clades),
                        len(collapsed_away),
                    )
                except Exception:
                    pass
            except Exception as e:
                logging.getLogger(__name__).warning("[RUN] Failed to build sequence assessment: %s", e)
    except Exception as e:
        logging.getLogger(__name__).warning("[RUN] Failed to build placement warnings: %s", e)

    # At end of run, keep the externally useful deliverables and organise them
    # into high-level folders.
    try:
        _organize_run_outputs(outdir, primary_db_name=run_main_db_name)
    except Exception as e:
        logging.getLogger(__name__).warning("[RUN] Failed to organise output files: %s", e)

    # Write a manifest after files have been organised.
    try:
        _write_output_explanations(outdir)
    except Exception:
        pass


def _detect_taxon_rank(taxon_query: str, rank_arg: str) -> str:
    """Return the single-letter rank key to filter on.

    Priority:
      1. Explicit ``--rank`` argument (if not 'auto')
      2. GTDB-style prefix embedded in the query (``p__``, ``f__``, ``g__`` …)
      3. Known domain keywords: ``archaea`` / ``bacteria`` → 'd'
      4. Fallback: 'p' (phylum)
    """
    from branchmanager.taxonomy import RANK_ALIASES
    if rank_arg and rank_arg.lower() != 'auto':
        return RANK_ALIASES.get(rank_arg.lower(), rank_arg.lower()[:1])
    if '__' in taxon_query:
        prefix = taxon_query.split('__')[0].lower().strip()
        return RANK_ALIASES.get(prefix, prefix[:1] if prefix else 'p')
    # Bare domain-level keywords — map to 'd' so we filter at domain rank
    if taxon_query.strip().lower() in ('archaea', 'bacteria', 'eukarya', 'eukaryota'):
        return 'd'
    return 'p'


def _strip_rank_prefix(taxon: str) -> str:
    """Strip a GTDB-style rank prefix: 'p__Bacillota' → 'Bacillota'."""
    if '__' in taxon:
        return taxon.split('__', 1)[1].strip()
    return taxon.strip()


def _taxon_name_matches(val: str, query: str) -> bool:
    """Case-insensitive, underscore/space-flexible exact match."""
    import re as _re
    def _n(s): return _re.sub(r'[\s_]+', '_', str(s).strip().lower())
    return bool(val) and _n(val) == _n(query)


def cmd_subtree(args):
    """Build a focused subtree for a specific taxon from an existing DB.

    Fast path: if a ``current_alignment.fasta`` already exists in
    ``--from-dir`` (or the output directory), the matching sequences are
    extracted from that *pre-built* alignment and FastTree is run on the
    subset — no re-alignment needed.

    Slow path (fallback): if no existing alignment is found, the matching
    sequences are exported from the DB and a full MAFFT + FastTree build
    is performed (same as a normal run).
    """
    from branchmanager.pipeline import itol as itol_mod
    from branchmanager.pipeline import tree as tree_mod  # noqa: F401 (used in _build_subtree)
    from branchmanager.pipeline.workflow_helpers import write_combined_taxonomy_tsv
    from branchmanager.taxonomy import parse_taxon_string
    from branchmanager.utils.fasta import read_fasta, write_fasta  # noqa: F401

    outdir = args.out
    os.makedirs(outdir, exist_ok=True)
    _configure_logging(outdir)
    log = logging.getLogger(__name__)

    taxon_query = args.taxon
    rank_key = _detect_taxon_rank(taxon_query, getattr(args, 'rank', 'auto'))
    taxon_clean = _strip_rank_prefix(taxon_query)
    threads = getattr(args, 'threads', 4)
    min_seqs = getattr(args, 'min_seqs', 3)
    ref_fasta = getattr(args, 'ref', None)
    anchor_file = getattr(args, 'anchors', None)
    from_dir = getattr(args, 'from_dir', None) or outdir
    no_tree = getattr(args, 'no_tree', False)

    log.info(
        "[SUBTREE] Taxon query: '%s'  rank: '%s'  clean name: '%s'  from-dir: %s",
        taxon_query, rank_key, taxon_clean, from_dir,
    )

    # ── Query DB ──────────────────────────────────────────────────────────────
    db = Database(args.db)
    db.initialise()

    with db.connect() as conn:
        cur = conn.cursor()
        cur.execute(
            "SELECT s.id, s.sequence, t.taxonomy, t.confidence, s.dataset "
            "FROM sequences s LEFT JOIN taxonomy t ON s.id = t.id"
        )
        all_rows = cur.fetchall()

    # ── Filter to matching taxon ──────────────────────────────────────────────
    matched = []
    for sid, seq, tax, conf, dataset in all_rows:
        if not tax:
            continue
        parsed = parse_taxon_string(tax)
        val = parsed.get(rank_key, '')
        if _taxon_name_matches(val, taxon_clean):
            matched.append((sid, seq or '', tax, conf or 'NA', dataset or ''))

    log.info("[SUBTREE] DB has %d sequences total; %d match taxon '%s' at rank '%s'",
             len(all_rows), len(matched), taxon_clean, rank_key)

    if len(matched) < min_seqs:
        log.warning(
            "[SUBTREE] Only %d sequences matched (minimum required: %d). "
            "Check taxon spelling and rank. Available phyla in DB: %s",
            len(matched), min_seqs,
            sorted({parse_taxon_string(r[2]).get('p', 'unknown')
                    for r in all_rows if r[2]}),
        )
        return

    # Collect matched IDs as a set for fast lookup
    matched_ids = {sid for sid, *_ in matched}

    # ── Write taxonomy TSV ────────────────────────────────────────────────────
    combined_tax_path = Path(outdir) / 'subtree_combined_taxonomy.tsv'
    write_combined_taxonomy_tsv(
        combined_tax_path,
        [(sid, tax, conf) for sid, _, tax, conf, _ in matched],
    )
    log.info("[SUBTREE] Wrote taxonomy for %d sequences → %s", len(matched), combined_tax_path)

    # ── Write summary TSV ─────────────────────────────────────────────────────
    summary_path = Path(outdir) / 'subtree_sequence_list.tsv'
    with open(summary_path, 'w') as sf:
        sf.write('ID\tTaxonomy\tConfidence\tDataset\n')
        for sid, _, tax, conf, dataset in matched:
            sf.write(f"{sid}\t{tax}\t{conf}\t{dataset}\n")
    log.info("[SUBTREE] Wrote sequence list → %s", summary_path)

    # ── Tree building ─────────────────────────────────────────────────────────
    tree_path = Path(outdir) / 'subtree_tree.nwk'

    if not no_tree:
        _build_subtree(
            matched=matched,
            matched_ids=matched_ids,
            from_dir=from_dir,
            outdir=outdir,
            tree_path=tree_path,
            ref_fasta=ref_fasta,
            anchor_file=anchor_file,
            threads=threads,
            taxon_clean=taxon_clean,
            log=log,
        )
    else:
        log.info("[SUBTREE] Tree build skipped (--no-tree).")

    # ── iTOL files ────────────────────────────────────────────────────────────
    try:
        tfile = str(tree_path) if tree_path.exists() else None
        itol_mod.generate_itol_colors(
            str(combined_tax_path),
            outdir,
            tree_file=tfile,
            phylum_groups=getattr(args, 'group_phyla', None),
        )
        log.info("[SUBTREE] Generated iTOL colour files in %s", outdir)
    except Exception as e:
        log.warning("[SUBTREE] iTOL colour generation failed: %s", e)

    # Optional: write functional annotation datasets when provided
    try:
        func_tsv = getattr(args, 'functional', None)
        if func_tsv:
            try:
                written = itol_mod.write_functional_annotations(str(func_tsv), outdir, id_map=None)
                log.info("[SUBTREE] Wrote functional annotation iTOL files: %s", ','.join(written) if written else '(none)')
            except Exception as e:
                log.warning("[SUBTREE] Functional annotations generation failed: %s", e)
    except Exception:
        pass
    # Draft rumen functional groups
    if getattr(args, 'draft_rumen_functions', False):
        try:
            tsv_out, itol_out = itol_mod.generate_rumen_function_draft(
                str(combined_tax_path), outdir, id_map=None
            )
            if tsv_out:
                log.info("[SUBTREE] Draft rumen functional annotation: %s", tsv_out)
            if itol_out:
                log.info("[SUBTREE] Rumen functional iTOL file: %s", itol_out)
        except Exception as e:
            log.warning("[SUBTREE] Draft rumen functions generation failed: %s", e)

    # ── Dataset membership strip ──────────────────────────────────────────────
    # One colour per dataset label so users can see which sequences came from
    # which dataset (preload, run batch, etc.) in the same tree view.
    try:
        ids_in_order = [sid for sid, *_ in matched]
        ds_map = {sid: (dataset or 'unknown') for sid, _, _tax, _conf, dataset in matched}
        membership_path = Path(outdir) / 'itol_dataset_membership.itol'
        itol_mod.write_dataset_membership_strip(
            str(membership_path), ids_in_order, ds_map,
            dataset_label='Dataset membership',
        )
        log.info("[SUBTREE] Wrote dataset membership strip → %s", membership_path)
    except Exception as e:
        log.warning("[SUBTREE] Dataset membership strip failed: %s", e)

    # ── Done ──────────────────────────────────────────────────────────────────
    log.info(
        "[SUBTREE] Complete. %d sequences, taxon='%s', output=%s",
        len(matched), taxon_query, outdir,
    )
    print(
        f"[subtree] Done: {len(matched)} sequences for '{taxon_query}' "
        f"→ {outdir}"
    )


def _build_subtree(
    matched, matched_ids, from_dir, outdir, tree_path,
    ref_fasta, anchor_file, threads, taxon_clean, log,
):
    """Build a FastTree for the given matched sequences.

    Fast path: filter an existing ``current_alignment.fasta`` from
    ``from_dir`` to the matched IDs + any reference anchors already
    embedded in that alignment, then run FastTree directly — no
    re-alignment required.

    Slow path: write the raw sequences as a FASTA and run the full
    ``initialise_or_update_tree`` pipeline (MAFFT → FastTree).
    """
    from branchmanager.pipeline.tree import (
        _run_fasttree, _make_unique_fasta,
        is_ref_anchor, get_anchor_file,
    )
    from branchmanager.utils.fasta import read_fasta, write_fasta
    from branchmanager.pipeline import tree as tree_mod
    import re

    out = Path(outdir)
    resolved_anchor = get_anchor_file(anchor_file)

    # ── Fast path ─────────────────────────────────────────────────────────────
    # Look for an existing alignment in from_dir (or outdir)
    aln_candidates = [
        Path(from_dir) / 'current_alignment.fasta',
        Path(from_dir) / 'tree' / 'current_alignment.fasta',
        Path(outdir) / 'current_alignment.fasta',
        Path(outdir) / 'tree' / 'current_alignment.fasta',
    ]
    existing_aln = next((p for p in aln_candidates if p.exists()), None)

    if existing_aln:
        log.info("[SUBTREE] Fast path: filtering existing alignment %s", existing_aln)
        # Build a normalised ID lookup for the matched IDs
        def _norm_id(x):
            x = str(x).split()[0]
            if '|' in x:
                x = x.split('|')[-1]
            x = re.sub(r'_[0-9]+$', '', x)
            return x.lower()

        matched_norm = {_norm_id(sid) for sid in matched_ids}

        kept = []
        anchors_included = 0
        for header, seq in read_fasta(str(existing_aln)):
            hid = header.split()[0]
            if is_ref_anchor(hid):
                kept.append((header, seq))
                anchors_included += 1
            elif hid in matched_ids or _norm_id(hid) in matched_norm:
                kept.append((header, seq))

        n_data = len(kept) - anchors_included
        log.info(
            "[SUBTREE] Fast path: %d data sequences + %d anchors extracted from alignment",
            n_data, anchors_included,
        )

        if n_data < 3:
            log.warning(
                "[SUBTREE] Only %d data sequences found in existing alignment "
                "(IDs may not match). Falling back to slow path.", n_data,
            )
            existing_aln = None   # trigger slow path below
        else:
            subtree_aln = out / 'subtree_alignment.fasta'
            write_fasta(kept, str(subtree_aln))
            fasta_for_tree, id_map = _make_unique_fasta(str(subtree_aln), out)
            if _run_fasttree(Path(fasta_for_tree), tree_path):
                _finalise_tree_subtree(out, id_map, tree_path, log)
            else:
                log.warning("[SUBTREE] FastTree failed on fast path")
            return

    # ── Slow path ─────────────────────────────────────────────────────────────
    log.info("[SUBTREE] Slow path: full MAFFT + FastTree build")
    seqs_fasta = out / 'subtree_input_sequences.fasta'
    write_fasta([(sid, seq) for sid, seq, *_ in matched], str(seqs_fasta))

    # Empty FASTA as user_fasta; all sequences go in via the FASTA directly
    # by passing them as user_fasta and disabling DB pull (db=None)
    try:
        tree_mod.initialise_or_update_tree(
            ref_fasta=ref_fasta or '',
            user_fasta=str(seqs_fasta),
            outdir=outdir,
            db=None,
            threads=threads,
            anchor_file=anchor_file,
            force_rebuild=True,
        )
        # The standard pipeline writes to current_tree.nwk; copy/rename
        default_tree = out / 'current_tree.nwk'
        if default_tree.exists():
            import shutil
            shutil.copy2(str(default_tree), str(tree_path))
            log.info("[SUBTREE] Slow path tree written → %s", tree_path)
    except Exception as e:
        log.warning("[SUBTREE] Slow path tree build failed: %s", e)


def _finalise_tree_subtree(out: Path, id_map: dict, tree_path: Path, log) -> None:
    """Remap IDs, prune anchors, and label internal nodes for subtree output."""
    from branchmanager.pipeline.tree import (
        _repair_internal_node_label_delimiters, _label_internal_nodes,
        _prune_anchor_leaves, REF_ANCHOR_PREFIX,
    )
    if not tree_path.exists():
        return
    newick = tree_path.read_text()
    if id_map:
        for new_id, orig_id in id_map.items():
            newick = newick.replace(new_id, orig_id)
    before = newick.count(REF_ANCHOR_PREFIX)
    newick = _prune_anchor_leaves(newick)
    after = newick.count(REF_ANCHOR_PREFIX)
    if before:
        log.info("[SUBTREE] Pruned %d anchor leaves from subtree newick", before - after)
    newick = _repair_internal_node_label_delimiters(newick)
    newick = _label_internal_nodes(newick)
    tree_path.write_text(newick)
    log.info("[SUBTREE] Subtree finalised → %s", tree_path)


def cmd_regen_itol(args):
    db = Database(args.db)
    db.initialise()
    outdir = args.out
    os.makedirs(outdir, exist_ok=True)
    _configure_logging(outdir)
    log = logging.getLogger(__name__)

    try:
        with db.connect() as conn:
            cur = conn.cursor()
            # If the outdir contains a preload combined taxonomy file, prefer
            # to use only those IDs so the regenerated iTOL matches preload
            # statistics. This mirrors the behaviour used during `run` when a
            # --preload-dir is supplied.
            try:
                p = Path(outdir)
                # Prefer an existing combined taxonomy file in the outdir, but
                # only if it contains data (header + >=1 data row). Fall back
                # to the DB-wide query otherwise to avoid regenerating empty
                # iTOL outputs when a stub file exists.
                preload_file = None
                for cand in _combined_taxonomy_candidates(p):
                    if not cand.exists():
                        continue
                    try:
                        # count lines cheaply
                        with open(cand) as pf:
                            cnt = sum(1 for _ in pf)
                        if cnt > 1:
                            preload_file = cand
                            break
                        else:
                            # file exists but only header or empty -> ignore
                            continue
                    except Exception:
                        continue
                preload_ids = None
                if preload_file is not None:
                    preload_ids = read_combined_taxonomy_ids(preload_file)
                else:
                    preload_ids = None
            except Exception:
                preload_ids = None

            preload_id_list = preload_ids if isinstance(preload_ids, list) else []

            # fetch id, taxonomy, confidence and dataset for all sequences (or filter)
            if getattr(args, 'include_datasets', None):
                ds_list = [d.strip() for d in args.include_datasets.split(',') if d.strip()]
                placeholders = ','.join('?' for _ in ds_list)
                cur.execute(f"SELECT s.id, t.taxonomy, t.confidence, s.dataset FROM sequences s LEFT JOIN taxonomy t ON s.id = t.id WHERE s.dataset IN ({placeholders})", tuple(ds_list))
            elif preload_id_list:
                placeholders = ','.join('?' for _ in preload_id_list)
                cur.execute(f"SELECT s.id, t.taxonomy, t.confidence, s.dataset FROM sequences s LEFT JOIN taxonomy t ON s.id = t.id WHERE s.id IN ({placeholders})", tuple(preload_id_list))
            else:
                cur.execute("SELECT s.id, t.taxonomy, t.confidence, s.dataset FROM sequences s LEFT JOIN taxonomy t ON s.id = t.id")
            rows = cur.fetchall()
    except Exception as e:
        log.warning("[REGEN-ITOL] Failed to query DB: %s", e)
        return

    try:
        kingdom = _sequence_domain_filter(args, default=None)
        if kingdom:
            kingdom_text = str(kingdom)
            rows = [row for row in rows if row[1] and taxonomy_matches_kingdom(str(row[1]), kingdom_text)]
            log.info("[REGEN-ITOL] Retained %d rows after kingdom filter=%s", len(rows), kingdom)
    except Exception as e:
        log.warning("[REGEN-ITOL] Kingdom filter failed: %s", e)

    combined_path = Path(outdir) / 'combined_taxonomy.tsv'
    try:
        write_combined_taxonomy_tsv(combined_path, [(rid, tax, conf) for rid, tax, conf, ds in rows])
        log.info("[REGEN-ITOL] Wrote combined taxonomy for %d ids to %s", len(rows), combined_path)
    except Exception as e:
        log.warning("[REGEN-ITOL] Failed to write combined taxonomy: %s", e)
        return

    # call itol generator
    try:
        tree_path = Path(outdir) / 'current_tree.nwk'
        tfile = str(tree_path) if tree_path.exists() else _find_tree_file_in_dir(outdir)
        itol.generate_itol_colors(str(combined_path), outdir, tree_file=tfile, phylum_groups=getattr(args, 'group_phyla', None))
        log.info("[REGEN-ITOL] Generated iTOL color files in %s", outdir)
    except Exception as e:
        log.warning("[REGEN-ITOL] Color generation failed: %s", e)

    # Optional: write functional annotation datasets when provided
    try:
        func_tsv = getattr(args, 'functional', None)
        if func_tsv:
            try:
                written = itol.write_functional_annotations(str(func_tsv), outdir, id_map=None)
                log.info("[REGEN-ITOL] Wrote functional annotation iTOL files: %s", ','.join(written) if written else '(none)')
            except Exception as e:
                log.warning("[REGEN-ITOL] Functional annotations generation failed: %s", e)
    except Exception:
        pass
    # Draft rumen functional groups
    if getattr(args, 'draft_rumen_functions', False):
        try:
            tsv_out, itol_out = itol.generate_rumen_function_draft(str(combined_path), outdir, id_map=None)
            if tsv_out:
                log.info("[REGEN-ITOL] Draft rumen functional annotation: %s", tsv_out)
            if itol_out:
                log.info("[REGEN-ITOL] Rumen functional iTOL file: %s", itol_out)
        except Exception as e:
            log.warning("[REGEN-ITOL] Draft rumen functions generation failed: %s", e)

    # build dataset membership strip
    try:
        ids_in_order = [r[0] for r in rows]
        ds_map = {r[0]: (r[3] or '') for r in rows}
        membership_path = Path(outdir) / 'itol_dataset_membership.itol'
        itol.write_dataset_membership_strip(str(membership_path), ids_in_order, ds_map)
        log.info("[REGEN-ITOL] Wrote dataset membership ITOL to %s", membership_path)
    except Exception as e:
        log.warning("[REGEN-ITOL] Failed to build/write dataset membership ITOL: %s", e)

    # write explanations for regenerated outputs
    try:
        _write_output_explanations(outdir)
    except Exception:
        pass


def build_parser():
    parser = argparse.ArgumentParser(
        prog='branchmanager',
        description=(
            'BranchManager — marker-gene QC, taxonomy, novelty scoring, and isolate prioritisation toolkit.\n\n'
            'Subcommands:\n'
            '  preclassify Pre-classify reference FASTA collections (Hungate, SILVA …) once and reuse.\n'
            '  preload     Load a baseline dataset (e.g. Hungate) and build the backbone tree.\n'
            '  evaluate    Core partner-sequence evaluation workflow (alias: run/eval).\n'
            '  run         Process new sequences against the baseline; score novelty and update the tree.\n'
            '  subtree     Extract a focused tree and iTOL files for a specific taxon from an existing DB.\n'
            '  regen-itol  Regenerate iTOL colour files from an existing DB without re-running analysis.\n\n'
            'Typical workflow:\n'
            '  0. branchmanager preclassify --dataset hungate16s=hungate.fasta --ref gtdb.fna --taxa gtdb_tax.tsv -o preclassify_out\n'
            '  1. branchmanager preload  --fasta baseline.fasta --db project.db --dataset Hungate --taxa-assignments preclassify_out/pipeline_taxonomy.tsv --build-tree -o preload_out\n'
            '  2. branchmanager evaluate --input new_seqs.fasta --partner-metadata new_seqs_metadata.tsv --db project.db --dataset Batch1  --ref gtdb.fna --baseline-fasta hungate.fasta --baseline-dataset Hungate --mwl MWL.xlsx -o eval_out\n'
            '  3. branchmanager subtree  --db project.db --taxon archaea --from-dir preload_out -o archaea_out\n'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest='command')

    # ── preload ───────────────────────────────────────────────────────────────
    preload = sub.add_parser(
        'preload',
        help='Load a baseline dataset and build the reference tree.',
        description=(
            'Load a baseline FASTA dataset (e.g. Hungate 16S) into the DB, optionally classify\n'
            'sequences against a reference (GTDB/SILVA), collapse near-identical sequences,\n'
            'and build the backbone phylogenetic tree.\n\n'
            'This is always the first step. All subsequent `run` calls measure novelty against\n'
            'sequences stored by this command.'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    preload.add_argument('--fasta', required=True,
        help='Input FASTA file containing the baseline sequences to load.')
    preload.add_argument('--db', required=True,
        help='Path to the BranchManager SQLite database (created if it does not exist).')
    preload.add_argument('-o', '--out', required=False, default='.',
        help='Output directory for tree, iTOL files, and reports (default: current directory).')
    preload.add_argument('--dataset', required=True,
        help='Label for this dataset stored in the DB (e.g. Hungate). Used to colour iTOL strips.')
    preload.add_argument('--shorten-ids', dest='shorten_ids',
        action=argparse.BooleanOptionalAction, default=False,
        help='Replace input headers with compact IDs (e.g. HUN001). Default is to preserve the IDs exactly as supplied.')
    preload.add_argument('--classify', action='store_true',
        help='Classify sequences against --ref and store taxonomy in the DB. Requires --ref.')
    preload.add_argument('--build-tree', action='store_true',
        help='Build the backbone MAFFT + FastTree phylogenetic tree after loading.')
    preload.add_argument('--ref', required=False,
        help='Reference FASTA (GTDB/SILVA reps) for classification and tree orientation. Preferred over --taxa-assignments for externally classified inputs.')
    preload.add_argument('--taxa', required=False,
        help='Taxonomy table matching IDs in --ref (TSV/CSV, optionally .gz: id<TAB>lineage or id,lineage). Optional when --ref FASTA headers already contain lineages.')
    preload.add_argument('--ref-name', dest='ref_name', required=False, default=None,
        help='Display name for the primary reference database (default: derived from --ref filename). Used to label taxonomy columns.')
    preload.add_argument('--alt-ref', dest='alt_ref', action='append', default=None, metavar='FASTA',
        help='Additional reference FASTA to classify against (repeatable). Produces extra taxonomy columns in output files.')
    preload.add_argument('--alt-taxa', dest='alt_taxa', action='append', default=None, metavar='TABLE',
        help='Taxonomy TSV/CSV for the corresponding --alt-ref (positionally paired; repeatable; .gz accepted).')
    preload.add_argument('--alt-ref-name', dest='alt_ref_name', action='append', default=None, metavar='NAME',
        help='Display name for the corresponding --alt-ref (positionally paired; repeatable). Default: derived from filename.')
    preload.add_argument('--main-ref', dest='main_ref', required=False, default=None,
        help='Name of the reference database to use as the primary taxonomy source (default: primary --ref). Must match one of the --ref-name / --alt-ref-name values.')
    preload.add_argument('--taxa-assignments', '--taxa-aasignments',
        dest='taxa_assignments', required=False,
        help='Pre-computed taxonomy assignments for the INPUT sequences (TSV/CSV, optionally .gz: query_id + lineage, or a FASTA with embedded lineages). Use this instead of --classify when you already have taxonomy.')
    preload.add_argument('--collapse', action='store_true',
        help='Collapse sequences that share ≥ --collapse-threshold identity AND the same taxonomy into a single representative for the tree. Saves time and reduces visual clutter.')
    preload.add_argument('--collapse-threshold', type=float, default=99.8,
        help='Identity threshold (percent) for collapsing duplicate-like sequences (default: 99.8).')
    preload.add_argument('--sequence-domain', '--organism-domain', dest='sequence_domain',
        choices=SEQUENCE_DOMAIN_CHOICES, default=None,
        help=(
            'Sequence/domain profile for this preload. Default behavior is bacteria. '
            'Use archaea for archaeal 16S, fungi for fungal/eukaryotic runs with suitable refs/anchors, '
            'or mixed/all/none to disable domain filtering.'
        ))
    preload.add_argument('--anchors', required=False, default=None,
        help='Custom reference anchor FASTA for tree topology scaffolding. Defaults to the 26-sequence bundled anchor set (src/branchmanager/data/reference_anchors.fasta).')
    preload.add_argument('--threads', type=int, required=False, default=4,
        help='Number of CPU threads for MAFFT and VSEARCH (default: 4).')
    preload.add_argument('--tree-method', dest='tree_method',
        choices=['fasttree', 'iqtree', 'iqtree-fast'], default='fasttree',
        help=(
            'Phylogenetic tree-building backend (default: fasttree). '
            'fasttree: approximate ML, GTR+CAT — fast. '
            'iqtree: full ML, GTR+G+I — more accurate and stable topology (slower; '
            'recommended for publication-quality trees or when FastTree produces unstable clades). '
            'iqtree-fast: IQ-TREE 2 with -fast flag — good compromise for exploratory runs. '
            'Requires iqtree2 in PATH when using iqtree/iqtree-fast.'
        ),
    )
    preload.add_argument('--colors', required=False,
        help='CSV file mapping sequence IDs to custom hex colours for iTOL (columns: id, color).')
    preload.add_argument(
        '--group-phyla', dest='group_phyla', action='append', default=None, metavar='SPEC',
        help=(
            'Collapse multiple phyla into a single colour in iTOL legends. Repeatable. '
            'Formats: "archaea" (all archaeal phyla), "bacteria" (all bacterial phyla), '
            '"Bacillota,Bacillota_I" (explicit list; label = first name), '
            '"Firmicutes:Bacillota,Bacillota_I" (named group).'
        ),
    )
    preload.add_argument('--functional', dest='functional', required=False, default=None,
        help=(
            'TSV file mapping sequence IDs to functional attributes (pathways, functions, '
            'traits, scores, etc.). Header row required; first column = sequence ID; '
            'subsequent columns = one functional attribute each. '
            'One iTOL file is generated per column: binary (0/1/yes/no) → DATASET_BINARY, '
            'numeric → DATASET_SIMPLEBAR, categorical → DATASET_COLORSTRIP.'
        ),
    )
    preload.add_argument('--draft-rumen-functions', dest='draft_rumen_functions',
        action='store_true', default=False,
        help=(
            'Auto-generate a draft rumen functional-group annotation from the output taxonomy. '
            'Maps each sequence to a broad ruminant microbiome functional category '
            '(e.g. Cellulolytic/Fibrolytic, Methanogenic Archaea, Butyrate Producers) '
            'and writes rumen_functions_draft.tsv + itol_func_Rumen_Functional_Group.itol. '
            'The draft TSV can be edited and re-supplied via --functional in future runs.'
        ),
    )

    # ── run ───────────────────────────────────────────────────────────────────
    run = sub.add_parser(
        'run',
        aliases=['evaluate', 'eval'],
        help='Process new sequences against the baseline; score novelty and update the tree.',
        description=(
            'Evaluate new partner 16S isolate sequences against the project baseline.\n\n'
            'The workflow classifies against GTDB (primary), optionally cross-checks NCBI/GG2/SILVA\n'
            'as --alt-ref databases, scores novelty and neighbourhood density against prior partner\n'
            'and preload/baseline sequences, updates the tree, and optionally matches GTDB taxonomy against\n'
            'the Most Wanted List via --mwl.\n\n'
            'Provide --baseline-fasta for context datasets such as Hungate when they have not already\n'
            'been loaded with `branchmanager preload`. Novelty is always relative to YOUR submitted data,\n'
            'not the full external reference. Each successive run extends the baseline, so scores\n'
            'become increasingly precise.'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    run.add_argument('--input', required=True,
        help='FASTA file of new sequences to analyse.')
    run.add_argument('--db', required=True,
        help='Path to the BranchManager SQLite database (must have been initialised with `branchmanager preload`).')
    run.add_argument('-o', '--out', required=True,
        help='Output directory for this run (sequence_assessment.tsv, novelty_metrics.tsv, tree, iTOL files, etc.).')
    run.add_argument('--dataset', required=True,
        help='Label for this batch of sequences stored in the DB (e.g. Batch1). Used in iTOL dataset-membership strip.')
    run.add_argument('--ref', required=False,
        help='Reference FASTA (GTDB/SILVA reps) used for classification and tree orientation. Same file used in preload.')
    run.add_argument('--taxa', required=False,
        help='Taxonomy table matching IDs in --ref (TSV/CSV, optionally .gz: id<TAB>lineage or id,lineage).')
    run.add_argument('--ref-name', dest='ref_name', required=False, default=None,
        help='Display name for the primary reference database (default: derived from --ref filename).')
    run.add_argument('--alt-ref', dest='alt_ref', action='append', default=None, metavar='FASTA',
        help='Additional reference FASTA to classify against (repeatable). Adds extra taxonomy columns to sequence_assessment.tsv and taxonomy_all_dbs.tsv.')
    run.add_argument('--alt-taxa', dest='alt_taxa', action='append', default=None, metavar='TABLE',
        help='Taxonomy TSV/CSV for the corresponding --alt-ref (positionally paired; repeatable; .gz accepted).')
    run.add_argument('--alt-ref-name', dest='alt_ref_name', action='append', default=None, metavar='NAME',
        help='Display name for the corresponding --alt-ref (positionally paired; repeatable).')
    run.add_argument('--main-ref', dest='main_ref', required=False, default=None,
        help='Name of the reference database to treat as primary (drives the main Taxonomy column). Default: primary --ref.')
    run.add_argument('--mwl', dest='mwl', required=False, default=None,
        help='Most Wanted List workbook/TSV/CSV. GTDB taxonomy is matched against this list and MWL columns are added to sequence_assessment.tsv.')
    run.add_argument('--mwl-sheet', dest='mwl_sheet', required=False, default='MWL_V1',
        help='Sheet name to read from an MWL .xlsx workbook (default: MWL_V1).')
    run.add_argument('--mwl-min-rank', dest='mwl_min_rank', required=False, default='p',
        choices=['domain', 'd', 'phylum', 'p', 'class', 'c', 'order', 'o', 'family', 'f', 'genus', 'g', 'species', 's'],
        help='Minimum matched rank required for an MWL hit (default: phylum). Domain-only MWL entries still match at domain.')
    run.add_argument('--partner-metadata', '--sequencing-metadata', dest='partner_metadata', required=False, default=None,
        help='CSV/TSV sidecar table for this run with sequence IDs, partner IDs, and whether each isolate was selected for full-genome sequencing. .gz is accepted. Required when using the evaluate/eval alias.')
    run.add_argument('--baseline-fasta', dest='baseline_fasta', required=False, default=None,
        help='Optional baseline/context FASTA to load before evaluating the new sequences (e.g. Hungate 16S).')
    run.add_argument('--baseline-dataset', dest='baseline_dataset', required=False, default='Baseline',
        help='Dataset label for --baseline-fasta in the DB and default cultured-baseline novelty pool (default: Baseline). Must differ from --dataset.')
    run.add_argument('--novelty-baseline-dataset', dest='novelty_baseline_datasets', action='append', default=[],
        help='Existing DB dataset label to include in the baseline/cultured novelty pool (repeatable; useful for Hungate plus other cultured isolate sets).')
    run.add_argument('--baseline-taxa-assignments', dest='baseline_taxa_assignments', required=False, default=None,
        help='Pre-computed taxonomy for --baseline-fasta (TSV/CSV, optionally .gz: sequence_id + lineage + optional confidence, or embedded-lineage FASTA). Skips baseline classification.')
    run.add_argument('--baseline-skip-classify', dest='baseline_skip_classify', action='store_true', default=False,
        help='Load --baseline-fasta into the DB without classifying it. Novelty still uses the baseline sequences, but taxonomy/iTOL context may be sparse.')
    run.add_argument('--baseline-colors', dest='baseline_colors', required=False, default=None,
        help='Optional CSV with baseline sequence colors, same format as preload --colors.')
    run.add_argument('--baseline-shorten-ids', dest='baseline_shorten_ids',
        action=argparse.BooleanOptionalAction, default=False,
        help='Replace baseline FASTA headers with compact IDs. Default is to preserve the IDs exactly as supplied.')
    run.add_argument('--taxa-assignments', '--taxa-aasignments',
        dest='taxa_assignments', required=False,
        help='Pre-computed taxonomy for the INPUT sequences (TSV/CSV, optionally .gz: query_id + lineage, or embedded-lineage FASTA).')
    run.add_argument('--preload-dir', dest='preload_dir', required=False,
        help='Path to the preload output directory. Used to seed the tree backbone alignment so only new sequences need aligning.')
    run.add_argument('--shorten-ids', dest='shorten_ids',
        action=argparse.BooleanOptionalAction, default=False,
        help='Replace input headers with compact IDs. Default is to preserve the IDs exactly as supplied.')
    run.add_argument('--min-len', dest='min_len', type=int, default=800,
        help='Minimum sequence length to retain (bp, default: 800). Shorter sequences are filtered out.')
    run.add_argument('--max-n', dest='max_n', type=int, default=5,
        help='Maximum number of ambiguous (N) bases allowed (default: 5).')
    run.add_argument('--collapse', action='store_true',
        help='Collapse near-identical same-taxonomy sequences into representatives for the tree.')
    run.add_argument('--collapse-threshold', type=float, default=99.8,
        help='Identity threshold (percent) for collapsing (default: 99.8).')
    run.add_argument('--sequence-domain', '--organism-domain', dest='sequence_domain',
        choices=SEQUENCE_DOMAIN_CHOICES, default=None,
        help=(
            'Sequence/domain profile for this evaluate/run. Omitted means bacteria. '
            'Use archaea for archaeal runs, fungi for fungal/eukaryotic runs with suitable references, '
            'or mixed/all/none to disable domain filtering. Provide domain-specific --ref/--alt-ref, '
            '--baseline-fasta, --preload-dir, and --anchors as needed.'
        ))
    run.add_argument('--phylum', required=False,
        help='Filter iTOL output to sequences assigned to this phylum (e.g. Bacillota). Does not affect novelty scoring.')
    run.add_argument('--target', required=False, default=None,
        help=(
            'FASTA of sequences to measure novelty against instead of the DB. '
            'Leave unset to use all sequences previously stored in the DB (recommended).'
        ),
    )
    run.add_argument('--force-rebuild', '--rebuild-tree', dest='force_rebuild', action='store_true', default=False,
        help='Rebuild the entire tree from scratch even when an existing alignment is present. '
             'When combined with --preload-dir, ignores the preload backbone and jointly estimates '
             'tree topology across all datasets. (--rebuild-tree is an alias for this flag.)')
    run.add_argument('--anchors', required=False, default=None,
        help='Custom reference anchor FASTA for tree scaffolding. Defaults to bundled anchors.')
    run.add_argument('--threads', dest='threads', type=int, default=4,
        help='CPU threads for MAFFT and VSEARCH (default: 4).')
    run.add_argument('--tree-method', dest='tree_method',
        choices=['fasttree', 'iqtree', 'iqtree-fast'], default='fasttree',
        help=(
            'Phylogenetic tree-building backend (default: fasttree). '
            'fasttree: approximate ML, GTR+CAT. '
            'iqtree: full ML, GTR+G+I (recommended for production/publication runs). '
            'iqtree-fast: IQ-TREE 2 with -fast flag (good for exploratory incremental runs).'
        ),
    )
    run.add_argument('--user-colors', dest='user_colors', required=False,
        help='CSV file mapping sequence IDs to custom hex colours for iTOL (columns: id, color).')
    run.add_argument(
        '--group-phyla', dest='group_phyla', action='append', default=None, metavar='SPEC',
        help=(
            'Collapse multiple phyla into one colour in iTOL legends. Repeatable. '
            'Formats: "archaea", "bacteria", "Bacillota,Bacillota_I", "Firmicutes:Bacillota,Bacillota_I".'
        ),
    )
    run.add_argument('--functional', dest='functional', required=False, default=None,
        help=(
            'TSV file mapping sequence IDs to functional attributes. '
            'Header row required; first column = sequence ID; subsequent columns = functional attributes. '
            'Generates one iTOL file per column (DATASET_BINARY / DATASET_SIMPLEBAR / DATASET_COLORSTRIP).'
        ),
    )
    run.add_argument('--draft-rumen-functions', dest='draft_rumen_functions',
        action='store_true', default=False,
        help=(
            'Auto-generate a draft rumen functional-group annotation from the output taxonomy. '
            'Writes rumen_functions_draft.tsv and itol_func_Rumen_Functional_Group.itol. '
            'The draft TSV can be edited and re-supplied via --functional in later runs.'
        ),
    )

    # ── regen-itol ────────────────────────────────────────────────────────────
    regen = sub.add_parser(
        'regen-itol',
        help='Regenerate iTOL colour files from an existing DB without re-running analysis.',
        description=(
            'Re-generate all iTOL colour strips (phylum, family, genus, dataset membership)\n'
            'from the taxonomy already stored in the DB. Useful after changing --group-phyla\n'
            'options or after manually editing the database.'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    regen.add_argument('--db', required=True,
        help='Path to the BranchManager SQLite database.')
    regen.add_argument('-o', '--out', required=True,
        help='Output directory where iTOL files will be written (should be the preload or run output dir).')
    regen.add_argument('--include-datasets', required=False,
        help='Comma-separated list of dataset names to include (default: all datasets in the DB).')
    regen.add_argument('--sequence-domain', '--organism-domain', dest='sequence_domain',
        choices=SEQUENCE_DOMAIN_CHOICES, default=None,
        help='Optional domain profile filter for regenerated outputs: bacteria, archaea, fungi, or mixed/all/none.')
    regen.add_argument(
        '--group-phyla', dest='group_phyla', action='append', default=None, metavar='SPEC',
        help=(
            'Collapse multiple phyla into one colour in iTOL legends. Repeatable. '
            'Formats: "archaea", "bacteria", "Bacillota,Bacillota_I", "Firmicutes:Bacillota,Bacillota_I".'
        ),
    )
    regen.add_argument('--functional', dest='functional', required=False, default=None,
        help=(
            'TSV file mapping sequence IDs to functional attributes. '
            'Generates one iTOL file per column (DATASET_BINARY / DATASET_SIMPLEBAR / DATASET_COLORSTRIP).'
        ),
    )
    regen.add_argument('--draft-rumen-functions', dest='draft_rumen_functions',
        action='store_true', default=False,
        help='Auto-generate rumen functional-group iTOL annotation from stored taxonomy.',
    )

    # ── subtree ───────────────────────────────────────────────────────────────
    subtree = sub.add_parser(
        'subtree',
        help='Build a focused tree and iTOL files for a specific taxon from an existing DB.',
        description=(
            'Extract all sequences matching a given taxon from the DB and build a focused\n'
            'phylogenetic tree for that group only.\n\n'
            'Fast path: if --from-dir points to a directory containing current_alignment.fasta\n'
            '(a preload or run output), sequences are sliced from the pre-built alignment and\n'
            'FastTree is run directly — no MAFFT re-alignment needed (~seconds for hundreds of seqs).\n\n'
            'Slow path: if no existing alignment is found, a full MAFFT + FastTree build is run.\n\n'
            'Taxon formats accepted:\n'
            '  archaea, bacteria          → all sequences at domain level\n'
            '  Bacillota, Bacteroidota    → phylum name (auto-detected)\n'
            '  p__Bacillota               → GTDB-prefixed phylum\n'
            '  f__Lachnospiraceae         → GTDB-prefixed family\n'
            '  g__Ruminococcus            → GTDB-prefixed genus\n'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    subtree.add_argument('--db', required=True,
        help='Path to the BranchManager SQLite database.')
    subtree.add_argument('-o', '--out', required=True,
        help='Output directory for the subtree results.')
    subtree.add_argument('--taxon', required=True,
        help='Taxon to extract (see description above for accepted formats).')
    subtree.add_argument(
        '--rank', required=False, default='auto',
        choices=['auto', 'domain', 'd', 'phylum', 'p', 'class', 'c',
                 'order', 'o', 'family', 'f', 'genus', 'g', 'species', 's'],
        help='Taxonomic rank to filter on. Default "auto" detects the rank from the taxon name or prefix.',
    )
    subtree.add_argument('--from-dir', dest='from_dir', required=False, default=None,
        help=(
            'Existing preload or run output directory containing current_alignment.fasta. '
            'Enables the fast path (sequences extracted from the existing MSA; no re-alignment).'
        ),
    )
    subtree.add_argument('--ref', required=False,
        help='Reference FASTA for orientation correction (slow-path full build only).')
    subtree.add_argument('--anchors', required=False, default=None,
        help='Custom reference anchor FASTA. Defaults to bundled anchors (26 NCBI RefSeq sequences).')
    subtree.add_argument('--threads', type=int, default=4,
        help='CPU threads for FastTree / MAFFT (default: 4).')
    subtree.add_argument('--min-seqs', dest='min_seqs', type=int, default=3,
        help='Minimum sequences required to proceed with tree building (default: 3).')
    subtree.add_argument('--no-tree', dest='no_tree', action='store_true', default=False,
        help='Skip tree building; only write taxonomy TSV, sequence list, and iTOL colour files.')
    subtree.add_argument(
        '--group-phyla', dest='group_phyla', action='append', default=None, metavar='SPEC',
        help=(
            'Collapse multiple phyla into one colour in iTOL legends. Repeatable. '
            'Formats: "archaea", "bacteria", "Bacillota,Bacillota_I", "Firmicutes:Bacillota,Bacillota_I".'
        ),
    )
    subtree.add_argument('--functional', dest='functional', required=False, default=None,
        help=(
            'TSV file mapping sequence IDs to functional attributes. Header row required; first column = sequence ID; subsequent columns = functional attributes. '
            'Generates one iTOL file per column (DATASET_BINARY / DATASET_SIMPLEBAR / DATASET_COLORSTRIP).'
        ),
    )
    subtree.add_argument('--draft-rumen-functions', dest='draft_rumen_functions',
        action='store_true', default=False,
        help='Auto-generate rumen functional-group iTOL annotation from stored taxonomy.',
    )

    # ── sanger / AB1 processing ───────────────────────────────────────────────
    sanger = sub.add_parser(
        'sanger',
        aliases=['ab1', 'ab1-to-fasta'],
        help='Convert Sanger AB1/sequence reads to trimmed FASTA and assemble primer reads per isolate.',
        description=(
            'Process Sanger chromatogram reads before evaluate.\n\n'
            'Inputs may be AB1/ABI files, FASTA, or FASTQ. AB1 files are base-called from '
            'PBAS/PCON tags, quality-trimmed, oriented by primer direction, and optionally '
            'assembled into one consensus 16S sequence per isolate. For example, 27F reads '
            'are kept forward and 907R reads are reverse-complemented before overlap assembly.\n\n'
            'If --sample-map/--read-metadata are omitted, sequence IDs and primer names are '
            'inferred from filenames such as Iso001_27F.ab1 and Iso001_907R.ab1.'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sanger.add_argument(
        '--input', nargs='+', required=False, default=[],
        help='AB1/ABI/FASTA/FASTQ files or directories to process. Directories are searched recursively by default. Optional when --sample-map lists the read files.',
    )
    sanger.add_argument('-o', '--out', required=True,
        help='Output directory for assembled.fasta, read_qc.tsv, visual reports, and assembly_report.tsv.')
    sanger.add_argument(
        '--sample-map', required=False, default=None,
        help=(
            'Optional CSV/TSV one row per isolate/sample. Use sequence_id/isolate_id/sample_id plus '
            'an ab1_files/read_files column containing ; separated files, or separate primer columns '
            'such as 27F and 907R. Add processing_mode/tags=assemble or best_read to override '
            'per-isolate handling. Relative paths are resolved next to the mapping file.'
        ),
    )
    sanger.add_argument(
        '--read-metadata', required=False, default=None,
        help=(
            'Optional CSV/TSV mapping read files to sequence_id, primer, and direction. '
            'Columns: file/read_file, sequence_id, primer, direction. Also accepts the same sample-level format as --sample-map.'
        ),
    )
    sanger.add_argument(
        '--primer', dest='primers', action='append', default=None,
        help='Primer name to recognise in filenames (repeatable). Defaults include common 16S primers such as 27F, 907R, and 1492R.',
    )
    sanger.add_argument('--min-quality', dest='min_quality', type=int, default=20,
        help='Phred cutoff for Mott-style end trimming (default: 20).')
    sanger.add_argument('--window', type=int, default=20,
        help='Window size retained in reports for trimming context; rigorous trimming uses the Phred cutoff directly.')
    sanger.add_argument('--min-length', dest='min_length', type=int, default=800,
        help='Minimum final sequence length to write to assembled.fasta (default: 800 bp).')
    sanger.add_argument('--min-read-length', dest='min_read_length', type=int, default=None,
        help='Minimum trimmed read length to retain before assembly/best-read selection. Defaults to --min-length.')
    sanger.add_argument('--min-mean-quality', dest='min_mean_quality', type=float, default=25.0,
        help='Minimum mean Phred score after trimming/masking for read and final QC (default: 25).')
    sanger.add_argument('--mask-quality', dest='mask_quality', type=int, default=20,
        help='Mask internal bases below this Phred score to N before assembly (default: 20).')
    sanger.add_argument('--max-read-expected-errors', dest='max_read_expected_errors', type=float, default=8.0,
        help='Maximum expected base-call errors allowed per retained read (default: 8).')
    sanger.add_argument('--max-output-expected-errors', dest='max_output_expected_errors', type=float, default=5.0,
        help='Maximum expected base-call errors allowed in the final isolate sequence (default: 5).')
    sanger.add_argument('--max-n-percent', dest='max_n_percent', type=float, default=1.0,
        help='Maximum percent N allowed after masking for reads and final output (default: 1.0).')
    sanger.add_argument('--max-internal-lowq-run', dest='max_internal_low_quality_run', type=int, default=20,
        help='Maximum internal low-quality/ambiguous run length before read failure (default: 20 bp).')
    sanger.add_argument('--max-conflict-density', dest='max_conflict_density', type=float, default=1.0,
        help='Maximum overlap conflicts per 100 final bases before final QC failure (default: 1.0).')
    sanger.add_argument('--quality-difference', dest='quality_difference', type=int, default=10,
        help='Minimum Phred difference required to choose one conflicting overlap base over another (default: 10).')
    sanger.add_argument('--allow-missing-quality', dest='allow_missing_quality',
        action='store_true', default=False,
        help='Allow AB1 reads missing PCON quality scores to pass with warnings. By default they fail QC.')
    sanger.add_argument('--min-overlap', dest='min_overlap', type=int, default=40,
        help='Minimum overlap length for assembling multiple primer reads (default: 40 bp).')
    sanger.add_argument('--min-overlap-identity', dest='min_overlap_identity', type=float, default=0.85,
        help='Minimum overlap identity for assembly, 0-1 (default: 0.85).')
    sanger.add_argument('--assemble', dest='assemble',
        action=argparse.BooleanOptionalAction, default=True,
        help='Assemble multiple reads per sequence_id when possible (default). Use --no-assemble to keep the best read.')
    sanger.add_argument('--recursive', dest='recursive',
        action=argparse.BooleanOptionalAction, default=True,
        help='Search input directories recursively (default).')

    # ── preclassify ───────────────────────────────────────────────────────────
    preclassify = sub.add_parser(
        'preclassify',
        help='Pre-classify reference FASTA collections (Hungate, SILVA, RDP …) so on-the-fly classification is not needed at run time.',
        description=(
            'Classify one or more reference FASTA collections against a reference database\n'
            'using vsearch and save the results so the main pipeline can reuse them without\n'
            're-classifying on every run.\n\n'
            'Known dataset names (use with --dataset NAME=FASTA):\n'
            '  hungate16s / hungate   Hungate1000 curated rumen 16S sequences\n'
            '  silva                  SILVA ribosomal RNA database\n'
            '  rdp                    Ribosomal Database Project\n'
            '  homd                   Human Oral Microbiome Database\n'
            '  greengenes2 / gg2      GreenGenes2 16S reference\n'
            '  ncbi16s                NCBI 16S rRNA RefSeq collection\n'
            '  gtdb                   GTDB 16S representative sequences\n'
            '  (any other name)       Treated as a custom dataset\n\n'
            'Outputs written to --out:\n'
            '  {name}_classification.tsv   Full classification (human-readable)\n'
            '  {name}_taxonomy.tsv         Condensed per-dataset taxonomy\n'
            '  {name}_taxonomic_disagreement.tsv  High-quality hits with conflicting taxa\n'
            '  combined_taxonomy.tsv       All datasets merged (with Dataset column)\n'
            '  pipeline_taxonomy.tsv       All datasets merged; pass to --taxa-assignments\n'
            '  preclassify_summary.txt     Plain-text summary with usage examples\n\n'
            'Example:\n'
            '  branchmanager preclassify \\\n'
            '    --dataset hungate16s=/data/hungate.fasta \\\n'
            '    --dataset silva=/data/silva_16s.fasta \\\n'
            '    --ref /data/gtdb_ssu_reps.fna \\\n'
            '    --taxa /data/gtdb_taxonomy.tsv.gz \\\n'
            '    --threads 8 -o preclassify_out/\n\n'
            'Then use the output in a preload:\n'
            '  branchmanager preload --fasta hungate.fasta \\\n'
            '    --taxa-assignments preclassify_out/pipeline_taxonomy.tsv \\\n'
            '    --db project.db --dataset Hungate -o preload_out/'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    preclassify.add_argument(
        '--dataset', dest='datasets', nargs='+', metavar='NAME FASTA',
        required=True,
        help=(
            'One or more NAME FASTA pairs, e.g.: '
            '--dataset Hungate /data/hun.fasta QUB /data/qub.fasta. '
            'Alternatively use NAME=FASTA pairs (repeatable): '
            '--dataset Hungate=/data/hun.fasta --dataset QUB=/data/qub.fasta. '
            'Known names: hungate16s, silva, rdp, homd, greengenes2, ncbi16s, gtdb. '
            'Any other name is treated as a custom dataset.'
        ),
    )
    preclassify.add_argument(
        '--ref', required=True,
        help='Reference FASTA (e.g. GTDB/SILVA reps) used by vsearch for classification.',
    )
    preclassify.add_argument(
        '--taxa', required=False, default=None,
        help='Taxonomy table matching IDs in --ref (TSV/CSV, optionally .gz: id<TAB>lineage or id,lineage). '
             'If omitted, taxonomy is parsed directly from reference FASTA headers.',
    )
    preclassify.add_argument(
        '-o', '--out', required=True,
        help='Output directory where all classification files will be written.',
    )
    preclassify.add_argument(
        '--threads', type=int, default=4,
        help='CPU threads for vsearch (default: 4).',
    )
    preclassify.add_argument(
        '--min-identity', dest='min_identity', type=float, default=0.80,
        help='Minimum vsearch alignment identity threshold 0–1 (default: 0.80 = 80%%).',
    )
    preclassify.add_argument(
        '--max-hits', dest='max_hits', type=int, default=10,
        help=(
            'Number of candidate hits vsearch collects per query '
            '(--maxaccepts / --maxhits).  The best hit by %% identity is '
            'selected after collection.  Higher values are more thorough but '
            'slower (default: 10).'
        ),
    )
    preclassify.add_argument(
        '--max-rejects', dest='max_rejects', type=int, default=256,
        help=(
            'vsearch --maxrejects: maximum number of non-matching candidate '
            'sequences examined before giving up on a query.  Raising this '
            'above vsearch\'s default of 32 helps classify ambiguous sequences '
            '(e.g. those with many N\'s) (default: 256).'
        ),
    )
    preclassify.add_argument(
        '--low-confidence-threshold', dest='low_confidence_threshold',
        type=float, default=0.97,
        help=(
            'Identity threshold below which a classified hit is flagged as low-confidence '
            'and written to *_low_confidence.tsv for manual review. '
            '0–1 (default: 0.97 = 97%%, the traditional species-level cutoff).'
        ),
    )

    return parser


def cmd_preclassify(args):
    """Handler for the ``preclassify`` subcommand."""
    from branchmanager.pipeline import preclassify as _preclassify_mod

    outdir = args.out
    os.makedirs(outdir, exist_ok=True)
    _configure_logging(outdir)
    log = logging.getLogger(__name__)

    # Parse datasets — supports two input styles:
    #   Style A (flat pairs):  --dataset Name1 /path1 Name2 /path2 ...
    #   Style B (NAME=FASTA):  --dataset Name1=/path1 --dataset Name2=/path2
    #   Style C (mixed):       --dataset Name1=/path1 Name2 /path2 ...
    raw = args.datasets or []
    # Flatten in case we ever switch back to action='append'
    flat: list = []
    for item in raw:
        if isinstance(item, list):
            flat.extend(item)
        else:
            flat.append(item)

    datasets = []
    i = 0
    while i < len(flat):
        token = flat[i]
        if '=' in token:
            # NAME=FASTA form
            name, _, fasta = token.partition('=')
            datasets.append((name.strip(), fasta.strip()))
            i += 1
        elif i + 1 < len(flat) and not flat[i + 1].startswith('-'):
            # NAME FASTA pair (next token is not a flag)
            datasets.append((token.strip(), flat[i + 1].strip()))
            i += 2
        else:
            raise SystemExit(
                f"[PRECLASSIFY] Cannot parse dataset spec at position {i}: '{token}'. "
                "Expected 'NAME FASTA' pairs or 'NAME=FASTA' pairs."
            )

    if not datasets:
        raise SystemExit("[PRECLASSIFY] At least one dataset is required. "
                         "Use: --dataset NAME /path/to/file.fasta")

    # Validate FASTA paths exist
    for name, fasta in datasets:
        if not os.path.exists(fasta):
            raise SystemExit(
                f"[PRECLASSIFY] FASTA file not found for dataset '{name}': {fasta}"
            )

    log.info("[PRECLASSIFY] Starting pre-classification for %d dataset(s)", len(datasets))

    pipeline_tsv = _preclassify_mod.run_preclassify(
        datasets=datasets,
        ref_fasta=args.ref,
        outdir=outdir,
        taxa_tsv=getattr(args, 'taxa', None),
        threads=int(getattr(args, 'threads', 4) or 4),
        min_identity=float(getattr(args, 'min_identity', 0.80) or 0.80),
        low_confidence_threshold=float(
            getattr(args, 'low_confidence_threshold', 0.97) or 0.97
        ),
        max_hits=int(getattr(args, 'max_hits', 10) or 10),
        max_rejects=int(getattr(args, 'max_rejects', 256) or 256),
    )

    log.info("[PRECLASSIFY] Done. Pipeline-ready taxonomy → %s", pipeline_tsv)
    print(
        f"[preclassify] Done.\n"
        f"  Pipeline taxonomy : {pipeline_tsv}\n"
        f"  Summary           : {os.path.join(outdir, 'preclassify_summary.txt')}\n\n"
        f"Use in preload:\n"
        f"  branchmanager preload --fasta <fasta> --taxa-assignments {pipeline_tsv} "
        f"--db project.db --dataset <name> -o preload_out/"
    )


def cmd_sanger(args):
    """Handler for the ``sanger`` / ``ab1`` subcommand."""
    from branchmanager.pipeline import sanger as _sanger_mod

    outdir = args.out
    os.makedirs(outdir, exist_ok=True)
    _configure_logging(outdir)
    primers = getattr(args, 'primers', None) or _sanger_mod.DEFAULT_PRIMERS
    inputs = getattr(args, 'input', None) or []
    sample_map = getattr(args, 'sample_map', None)
    read_metadata = getattr(args, 'read_metadata', None)
    if not inputs and not sample_map and not read_metadata:
        raise SystemExit('[sanger] Provide --input and/or --sample-map/--read-metadata.')
    outputs = _sanger_mod.run_sanger(
        inputs,
        outdir,
        read_metadata=read_metadata,
        sample_map=sample_map,
        primers=primers,
        min_quality=int(getattr(args, 'min_quality', 20) or 20),
        window=int(getattr(args, 'window', 20) or 20),
        min_length=int(getattr(args, 'min_length', 800) or 800),
        min_read_length=getattr(args, 'min_read_length', None),
        min_mean_quality=float(getattr(args, 'min_mean_quality', 25.0) or 25.0),
        mask_quality=int(getattr(args, 'mask_quality', 20) or 20),
        max_read_expected_errors=float(getattr(args, 'max_read_expected_errors', 8.0) or 8.0),
        max_output_expected_errors=float(getattr(args, 'max_output_expected_errors', 5.0) or 5.0),
        max_n_percent=float(getattr(args, 'max_n_percent', 1.0) or 1.0),
        max_internal_low_quality_run=int(getattr(args, 'max_internal_low_quality_run', 20) or 20),
        max_conflict_density=float(getattr(args, 'max_conflict_density', 1.0) or 1.0),
        quality_difference=int(getattr(args, 'quality_difference', 10) or 10),
        allow_missing_quality=bool(getattr(args, 'allow_missing_quality', False)),
        min_overlap=int(getattr(args, 'min_overlap', 40) or 40),
        min_overlap_identity=float(getattr(args, 'min_overlap_identity', 0.85) or 0.85),
        assemble=bool(getattr(args, 'assemble', True)),
        recursive=bool(getattr(args, 'recursive', True)),
    )
    logging.getLogger(__name__).info("[SANGER] Final assembled FASTA: %s", outputs['assembled_fasta'])
    print(
        "[sanger] Done.\n"
        f"  Assembled FASTA : {outputs['assembled_fasta']}\n"
        f"  Trimmed reads   : {outputs['trimmed_fasta']}\n"
        f"  Read QC         : {outputs['read_qc_tsv']}\n"
        f"  Per-base errors : {outputs['per_base_error_tsv']}\n"
        f"  Assembly report : {outputs['assembly_tsv']}\n\n"
        f"  Resequence list : {outputs['recommendations_tsv']}\n"
        f"  QC policy       : {outputs['qc_policy_tsv']}\n\n"
        f"  Read visual     : {outputs['read_error_svg']}\n"
        f"  Assembly visual : {outputs['assembly_svg']}\n\n"
        "Use in evaluate:\n"
        f"  branchmanager evaluate --input {outputs['assembled_fasta']} --partner-metadata <metadata.tsv> ..."
    )


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == 'preload':
        cmd_preload(args)
    elif args.command in ('run', 'evaluate', 'eval'):
        cmd_run(args)
    elif args.command == 'regen-itol':
        cmd_regen_itol(args)
    elif args.command == 'subtree':
        cmd_subtree(args)
    elif args.command == 'preclassify':
        cmd_preclassify(args)
    elif args.command in ('sanger', 'ab1', 'ab1-to-fasta'):
        cmd_sanger(args)
    else:
        parser.print_help()


if __name__ == '__main__':
    main()
