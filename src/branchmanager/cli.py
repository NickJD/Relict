"""BranchManager CLI — lightweight entrypoint for the branchmanager package.

Implements the BranchManager office workflow and technical reporting utilities.
"""
import argparse
import csv
import json
import logging
import os
import shutil
import sys
from pathlib import Path

# When this file is invoked directly, the package root (src)
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
from branchmanager.pipeline import (
    classify,
    itol,
    neighbourhood,
    novelty,
    paper_trail as _paper_trail_module,
    qc,
    quarterly_review,
    tree,
)
from branchmanager.pipeline import cluster_report as _cluster_report
from branchmanager.pipeline import mwl as _mwl
from branchmanager.pipeline import selection_sets as _selection_sets
from branchmanager.pipeline.classify import _derive_db_name as _classify_derive_db_name
from branchmanager.pipeline.collapse import collapse_fasta_within_taxa
from branchmanager.pipeline.workflow_helpers import (
    _assignment_source_is_fasta,
    build_orig_to_short_map as _build_orig_to_short_map_helper,
    build_placement_warning_rows,
    build_selection_decision,
    build_sequence_assessment_rows,
    classification_ids_matching_kingdom as _classification_ids_matching_kingdom_helper,
    collect_db_taxonomy_rows,
    iter_assignment_rows,
    iter_classification_rows,
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
from branchmanager.taxonomy import canonicalise_sequence_id, normalise_domain_query, taxonomy_matches_kingdom
from branchmanager.partner_metadata import load_partner_sequencing_metadata
from branchmanager.run_manifest import RunManifest, utc_now


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
        if command == 'performance-review':
            raise SystemExit(
                '[PERFORMANCE REVIEW] --partner-metadata / --sequencing-metadata is required: '
                'a cumulative CSV/TSV ledger with sequence IDs, partner acronyms, optional selection commitments, and already-sequenced status.'
            )
        return {}

    log = logging.getLogger(__name__)
    try:
        metadata_rows = load_partner_sequencing_metadata(metadata_path)
    except Exception as e:
        raise SystemExit(f'[PERFORMANCE REVIEW] Failed to read partner metadata {metadata_path}: {e}')

    run_id_set = {str(x) for x in run_ids}
    with db.connect() as conn:
        sequence_datasets = {
            str(sequence_id): str(dataset or '')
            for sequence_id, dataset in conn.execute('SELECT id, dataset FROM sequences')
        }
    dataset_roles = db.get_dataset_roles()
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
                cid = canonicalise_sequence_id(source_id)
            except Exception:
                cid = None
            if cid:
                mapped = orig_to_short.get(cid)
        if not mapped and source_id in run_id_set:
            mapped = source_id
        if not mapped:
            existing = db.resolve_sequence_id(source_id)
            if existing:
                mapped = existing
        if not mapped:
            warnings.append((source_id, 'metadata_id_not_found_in_project_database_or_current_run'))
            continue

        mapped_dataset = sequence_datasets.get(str(mapped), '')
        if dataset_roles.get(mapped_dataset, {}).get('role') == 'baseline':
            warnings.append((source_id, 'metadata_id_belongs_to_baseline_dataset'))
            continue

        matched_sources.add(str(mapped))
        resolved_rows.append({
            'id': mapped,
            'partner_id': row.get('partner_id') or source_id,
            # Preserve the dataset in which an earlier isolate entered the
            # project instead of relabelling it as part of the current batch.
            'dataset': mapped_dataset or getattr(args, 'dataset', ''),
            'selected_for_sequencing': bool(row.get('selected_for_sequencing')),
            'selected_for_wgs': bool(row.get('selected_for_wgs')),
            'source_id': source_id,
            'source_file': str(metadata_path),
            'raw_selected_value': row.get('raw_selected_value', ''),
            'raw_commitment_value': row.get('raw_commitment_value', ''),
        })

    for run_id in sorted(run_id_set - matched_sources):
        warnings.append((run_id, 'run_sequence_missing_from_partner_metadata'))

    inserted = db.upsert_sequencing_metadata(resolved_rows)
    log.info(
        '[PERFORMANCE REVIEW] Loaded partner sequencing metadata from %s: %d matched rows, %d warning(s)',
        metadata_path,
        inserted,
        len(warnings),
    )
    if warnings:
        warn_path = _write_partner_metadata_warnings(
            Path(outdir) / 'partner_metadata_warnings.tsv',
            warnings,
        )
        log.warning('[PERFORMANCE REVIEW] Partner metadata warnings written to %s', warn_path)

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
        p / 'filing_cabinet_id_map.tsv',
        p / 'user_id_map.tsv',
        p / 'user_id_map.csv',
        p / 'ids' / 'filing_cabinet_id_map.tsv',
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
        p / 'filing_cabinet_combined_taxonomy.tsv',
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
        "QC happens before classification, novelty scoring, and tree building. "
        "Sequences that fail QC are not included in downstream reports because they are too short "
        "for reliable 16S placement or contain too many ambiguous bases.",
        "",
        "Filtering rules:",
        "",
        f"- `too_short`: sequence length is less than `min_len` (`{stat('min_len')}` bp for this run).",
        f"- `too_many_n`: ambiguous `N` bases exceed `max_n_percent` (`{stat('max_n_percent')}`% for this run).",
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
            "No `qc.stats` file was found in this output directory. If this was Filing Cabinet only or a partially completed Performance Review, QC may not have been executed.",
            "",
        ])

    lines = [
        "# BranchManager Output Guide",
        "",
        "This guide explains the output files and metrics produced by the BranchManager Performance Review.",
        "",
        "## Run Scope",
        "",
        "Assessment, selection, per-database taxonomy, novelty, baseline-hit, and neighbourhood target rows are limited to the dataset supplied to this run. Prior partner datasets and cultured baselines remain available as comparison context. Neighbourhood context leaves, the combined tree taxonomy, MSA, tree, and iTOL membership files are cumulative project-context outputs by design.",
        "",
        "## Recommended Reading Order",
        "",
        "1. `assessment/sequencing_sets.tsv` - proposed primary, backup, and alternate isolates by GTDB species/local clade.",
        "2. `assessment/selection_summary.tsv` - decision-facing per-isolate recommendations and evidence.",
        "3. `assessment/neighbourhoods/` - labelled local-clade figures linked from both reports.",
        "4. `assessment/sequence_assessment.tsv` - full per-sequence audit table.",
        "5. `baseline/baseline_hits.tsv`, `taxonomy/`, `tree/`, and `itol/` - supporting assignments and phylogeny.",
        "6. `assessment/novelty_metrics.tsv` - detailed component calculations.",
        "7. `performance_review_dashboard.html` - compact linked view of the current Hiring Panel recommendations.",
        "",
        "## Directory Overview",
        "",
        "- `assessment/`: selection/audit reports, local-clade figures, novelty metrics, cluster reports, and warnings.",
        "- `baseline/`: nearest-hit reports and loaded baseline taxonomy/sequences.",
        "- `taxonomy/`: per-database classification outputs and combined taxonomy used for tree metadata.",
        "- `quality/`: marker provenance and reference-UCHIME results carried into selection evidence.",
        "- `tree/`: MSA, Newick tree, and tree/alignment warning files.",
        "- `itol/`: one iTOL metadata dataset per metadata type, usually colour-strip files.",
        "- `ids/`: short ID to original FASTA header maps.",
        "- `intermediate/`: QC, optional collapse, debug FASTA files, and scratch pools.",
        "- `logs/`: pipeline log.",
        "",
        *qc_section,
        "## Main Assessment Metrics",
        "",
        "`assessment/selection_summary.tsv` is the decision-facing table. It separates scientific value from redundancy and marker-evidence quality:",
        "",
        "- `Recommendation`: primary, backup, secondary, review, target-met, or already-sequenced action; the reason is explicit in `RecommendationReason`.",
        "- `SequencingSetID`, `SelectionGroupType`, `SequencingSetRole`, `SequencingSetRank`: connect each isolate to its nine-member diversity panel. `BMEXT_*` identifies a baseline-pangenome extension; `BMSET_*` identifies a candidate-only group. `PRIMARY` rows fill the target, `BACKUP` rows protect against DNA-extraction failure, and redundancy/boundary rows are retained without a rank.",
        "- `BaselineExtensionStatus`: explains whether exact GTDB species agreement, >=98.65% cultured-baseline identity, and >=95% query coverage support membership in a baseline-pangenome extension.",
        "- `EvidenceQuality`: whether the marker classification and placement warnings are strong enough to support selection.",
        "- `MarkerQC`, `MarkerReview`: Paper Trail/Merge Meeting evidence and any explicit manual decision. Unreviewed warning/unverified markers cannot become PRIMARY/BACKUP.",
        "- `CulturedGap`: distance from the cultured baseline: large below 97 percent, moderate from 97 to below 98.65 percent, small at or above 98.65 percent.",
        "- `ProjectCoverage`: number of other partner candidates at or above 97 percent identity, including the current rolling collection but excluding the query itself.",
        "- `ReferenceContext`: divergence from the external reference using 94.5 and 98.65 percent full-length 16S heuristic boundaries.",
        "- `CommittedGenomesSameAssessmentSpecies`, `SelectedPendingGenomesSameAssessmentSpecies`, `PangenomeTarget`, `PangenomeGap`: exact GTDB-species coverage and commitments. Baselines and completed partner genomes are available; selected pending isolates reserve planned slots. `SpeciesContext` adds a soft near-species or cluster-level note without relaxing the hard species gate.",
        "- `LocalTreeFigure`: relative path to the grouped local-clade PNG for visual inspection.",
        "",
        "These boundaries are decision-support heuristics, not declarations of a new species or genus. Genome analysis remains necessary for formal novelty claims.",
        "",
        "`assessment/sequence_assessment.tsv` is the main table. Important column groups:",
        "",
        "- `GTDBTaxonomy`, `GTDBClassificationHit`, `GTDBClassificationIdentity`, `GTDBClassificationConfidence`, `GTDBQueryCoverage`: authoritative GTDB assignment and supporting alignment extent.",
        "- `Taxonomy_<DB>`, `ClassificationHit_<DB>`, `Identity_<DB>`, `Confidence_<DB>`: assignment from additional databases such as GG2, SILVA, or NCBI.",
        "- `BaselineNearestHit`, `BaselineNearestHitDataset`, `BaselineNearestHitTaxonomy`, `BaselineNearestIdentity`: closest cultured isolate, for example Hungate.",
        "- `ProjectNearestHit`, `ProjectNearestIdentity`, `ProjectNoveltyScore`: comparison against all partner candidates in the rolling project collection, excluding self.",
        "- `GTDBReferenceNearestHit`, `GTDBReferenceNearestIdentity`, `GTDBReferenceNoveltyScore`: nearest hit and context against the GTDB reference FASTA.",
        "- `PartnerID`, `SelectedForGenomeSequencing`, `GenomeAlreadySequenced`, `BaselineGenomesSameAssessmentSpecies`, `SequencedPartnerGenomesSameAssessmentSpecies`, `SelectedPendingGenomesSameAssessmentSpecies`, `PangenomeGap`: rolling genome-collection context. Assessment species is GTDB for bacterial/archaeal Performance Reviews.",
        "- `SequencingSet*`: proposed primary/backup membership based on pangenome gap, evidence quality, and phylogenetic spread.",
        "- `InTree`, `ClusterRepresentative`, `ClusterSize`, `ClusteredMembers`: whether the sequence itself entered the tree or was represented by another clustered sequence.",
        "- `PlacementFlags`: warnings such as low classification identity, low nearest identity, or novelty/classification disagreement.",
        "",
        "## Novelty Metrics",
        "",
        "`assessment/novelty_metrics.tsv` contains cultured-baseline novelty, rolling project novelty, GTDB-reference context, and genome-collection coverage.",
        "",
        "- Leading `Nearest*`, `NoveltyScore`, `Crowding`, and `SequencingPriority` columns mirror the cultured-rumen comparison when a baseline pool exists; otherwise they mirror project novelty.",
        "- `Hungate*`, `SecondaryBaseline*`, and `CulturedRumen*` columns keep priority, secondary, and combined cultured-rumen evidence separate.",
        "- `Project*` columns compare against all partner candidate datasets, including the current rolling collection and excluding each query's self-hit. Baseline datasets are not mixed into this pool.",
        "- `Reference*` columns compare against the chosen external reference FASTA supplied with `--ref`, usually GTDB. These are separate from the baseline/project novelty scores.",
        "- `GenomeCollection*` columns compare against every baseline genome and every partner isolate with a genome already available. Exact same-GTDB-species counts drive the nine-genome target, while baseline identity and query coverage prevent near-identical cultured isolates from consuming diversity-panel ranks.",
        "- `NearestIdentity`: vsearch global-alignment percent identity to the nearest sequence in that pool.",
        "- `NearestQueryCoverage`, `NearestAlignmentLength`: how much marker evidence supports that identity; low-coverage identity must not be treated as a full-length equivalent.",
        "- `Novel`: `True` when nearest identity is below 97 percent.",
        "- `MatchesGE99`, `MatchesGE97`, `MatchesGE95`: number of pool sequences at or above 99, 97, and 95 percent identity. These describe how busy the local neighbourhood is.",
        "- `Crowding`: `isolated` when there is at most one hit at both 99 and 97 percent; `sparse` when <=3 hits at 97 percent; `moderate` when <=10 hits at 97 percent; otherwise `crowded`.",
        "- `NoveltyScore`: 0-100 score where higher means more novel and less crowded. It combines distance from the nearest hit with density bonuses for sparse neighbourhoods.",
        "- `SequencingPriority`: `HIGH` for <97 percent identity with few close neighbours, `MEDIUM` for moderately novel/sparse cases, otherwise `LOW`.",
        "- `DensitySource`: names the pool used, for example `baseline:Hungate`, `project_collection`, `target_fasta`, or `reference_fasta` fallback.",
        "",
        "Interpretation: a sequence far from Hungate but close to non-Hungate isolates is likely novel relative to cultured rumen isolate collections, but not necessarily novel relative to everything already supplied to the project.",
        "",
        "## Taxonomy Metrics",
        "",
        "- `ClassificationIdentity` is the vsearch percent identity to the best reference hit used for taxonomy assignment.",
        "- `ClassificationConfidence` is derived from the assignment parser/classifier output; higher means stronger taxonomy support.",
        "- `QueryCoverage`, `TargetCoverage`, `AlignmentLength`, `Mismatches`, and `Gaps` expose the alignment behind each best hit.",
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
        "- `itol/phylum.itol`, `itol/family.itol`, `itol/genus.itol`: taxonomy colour strips.",
        "- `itol/dataset_membership.itol`: dataset-of-origin colour strip.",
        "- `itol/novelty.itol`: novelty/nearest-identity colour strip.",
        "- `assessment/neighbourhoods/clade_*.png`: compact labelled subtrees with selection-set role/rank and MSA pident to the P1 isolate. Nearest baseline hits are forced into the displayed context.",
        "- `assessment/neighbourhoods/clade_*_pairwise_pident.tsv`: all displayed leaf-to-leaf MSA percent identities and compared-column counts.",
        "- `assessment/neighbourhoods/neighbourhood_manifest.tsv`: maps each assessed ID to its figure/pident table, P1 identity anchor, forced baseline hits, and displayed context.",
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
        ('combined_taxonomy.tsv', 'Combined taxonomy TSV. Columns: ID\tTaxon\tConfidence. Used to generate iTOL colour/legend files.'),
        ('loaded_baseline_taxonomy.tsv', 'Combined taxonomy TSV for baseline/provided datasets loaded before Performance Review. Columns: ID\tTaxon\tConfidence.'),
        ('filing_cabinet_combined_taxonomy.tsv', 'Combined taxonomy TSV for the Filing Cabinet baseline. Columns: ID\tTaxon\tConfidence.'),
        ('all_databases.tsv', 'Per-sequence taxonomic assignments from every configured reference database.'),
        ('baseline_hits.tsv', 'Nearest-hit report against baseline/provided datasets such as Hungate.'),
        ('nearest_project_hits_raw.tsv', 'Raw nearest-hit table produced for the newest submission; rolling project and baseline comparisons are summarised in novelty_metrics.tsv.'),
        ('qc_rejections.tsv', 'Per-sequence QC rejection table. Columns: ID, Length, NCount, NPercent, Reasons, MinLength, MaxNPercent. Reasons explain exactly why each sequence was filtered before downstream analysis.'),
        ('partner_metadata_warnings.tsv', 'Warnings from --partner-metadata mapping. Rows indicate metadata IDs not found in the project database/current run, baseline IDs supplied as partner metadata, or current-run sequences missing from the metadata table.'),
        ('marker_qc_provenance.tsv', 'Per-isolate bridge from Paper Trail/Merge Meeting raw reads and marker QC into Performance Review, including source checksums, manual review, and chimera call.'),
        ('chimera_screen.tsv', 'Reference-UCHIME marker screen. CHIMERA, INDETERMINATE, SKIPPED, or NOT_RUN calls force selection evidence to review.'),
        ('decision_changes.tsv', 'Difference between the two latest stored assessment snapshots, including recommendation, set role, pangenome gap, and tracked evidence changes.'),
        ('run_manifest.json', 'Machine-readable workflow manifest with input/output SHA256 checksums, software/tool versions, stage status, timestamps, warnings, and failure details.'),
        ('performance_review_dashboard.html', 'Compact linked decision dashboard for the current Performance Review/Hiring Panel.'),
        ('qc.stats', 'QC summary with total input, kept count, rejection counts, min_len, max_n_percent, and pointer to qc_rejections.tsv when available.'),
        ('user_id_map.tsv', 'Mapping of runtime sequence IDs to original headers produced when inserting user sequences into the DB. When shortening is disabled these usually match.'),
        ('filing_cabinet_id_map.tsv', 'Mapping of Filing Cabinet runtime IDs back to original FASTA headers. Use this to trace tree labels back to source records.'),
        ('*_id_map.tsv', 'ID map mapping original headers to runtime DB ids. Useful for iTOL and metadata tracing.'),
        ('collapsed_map.tsv', 'Cluster map for collapsed sequences: rep_id\ttaxonomy\tcount.'),
        ('filing_cabinet_collapsed_map.tsv', 'Cluster map for collapsed Filing Cabinet sequences: rep_id\ttaxonomy\tcount.'),
        ('collapsed_members.tsv', 'Member->representative mapping (member\trep) for collapsed clusters.'),
        ('filing_cabinet_collapsed_members.tsv', 'Member-to-representative mapping for collapsed Filing Cabinet clusters.'),
        ('current_dataset_sequences.fasta', 'QC-passing current-dataset FASTA with the runtime IDs used by the database. Exact marker duplicates retain separate isolate IDs.'),
        ('current_dataset_collapsed.fasta', 'Optional collapsed tree-view FASTA; assessment and database rows still retain every current-dataset isolate ID.'),
        ('filing_cabinet_collapsed.fasta', 'Collapsed Filing Cabinet FASTA; representatives retained for tree building.'),
        ('novelty_matches.tsv', 'vsearch BLAST-like output used to compute nearest-neighbour novelty identities.'),
        ('novelty_metrics.tsv', (
            'Per-sequence novelty metrics. The leading Nearest*/NoveltyScore columns mirror the '
            'cultured-rumen comparison when available; explicit Hungate*, SecondaryBaseline*, and '
            'CulturedRumen* columns separate baseline-tier evidence, and Project* columns compare against all partner candidates '
            'stored in the DB. Reference* columns compare against the external reference FASTA, usually '
            'GTDB. GenomeCollection* and Pangenome* columns add rolling genome-coverage context. '
            'DensitySource columns name the comparison pool.'
        )),
        ('sequence_assessment.tsv', (
            'Unified per-sequence assessment. '
            'COLUMN GROUPS: '
            '(1) TAXONOMY/CLASSIFICATION — GTDBTaxonomy, GTDBClassificationHit, '
            'GTDBClassificationIdentity, GTDBClassificationConfidence: authoritative GTDB assignment. '
            'Repeated as Taxonomy_<DB>, ClassificationHit_<DB>, Identity_<DB>, Confidence_<DB> '
            'for each additional --alt-ref database. '
            '(2) NOVELTY — NearestHit, NearestIdentity, MatchesGE*, NoveltyScore, Crowding, '
            'SequencingPriority: baseline/cultured novelty when a baseline pool exists. '
            'Project* columns repeat the same metrics against the rolling partner collection. '
            'Reference* columns repeat the same metrics against the external taxonomy reference '
            'FASTA, usually GTDB. GenomeCollection* and Pangenome* columns report available '
            'baseline genomes, already-sequenced partner genomes, and remaining same-species coverage. '
            '(3) TREE/CLUSTER — InTree, ClusterRepresentative, ClusterSize, ClusteredMembers: '
            'records whether the sequence entered the phylogenetic tree directly or was '
            'represented by a cluster representative after --collapse. '
            '(4) MWL — when --mwl is supplied, EvaluationScore and MWL* columns describe '
            'GTDB-based Most Wanted List matches and MWL priority contribution.'
        )),
        ('selection_summary.tsv', (
            'Concise scientific-advisory-board selection table. Contains one row per assessed '
            'sequence with a transparent decision, marker evidence quality, cultured gap, project '
            'coverage, external-reference context, taxonomy/MWL evidence, available-genome coverage, '
            'local-tree figure, and a short rationale. '
            'Use sequence_assessment.tsv for the full audit trail.'
        )),
        ('sequencing_sets.tsv', (
            'Rolling clade-level nine-member genome-sequencing diversity panel. Every candidate is '
            'placed in a BMEXT baseline-pangenome extension or a BMSET candidate-only group. A baseline '
            'extension requires exact GTDB species agreement, close marker identity, and adequate query '
            'coverage. PRIMARY rows '
            'fill the pangenome target; BACKUP and DIVERSITY_CANDIDATE rows preserve extraction-failure '
            'resilience and additional within-group spread. BASELINE_REDUNDANT rows remain auditable but '
            'are excluded using the reported identity and query-coverage thresholds; '
            'PANGENOME_BOUNDARY_REVIEW rows require lineage review before grouping.'
        )),
        ('neighbourhood_manifest.tsv', 'Maps each assessed sequence to its grouped local-clade image and pairwise-pident table, including the P1 identity anchor, forced nearest-baseline hits, displayed leaf count, assessed peers, baseline leaves, and already-sequenced leaves.'),
        ('visual_report_manifest.tsv', 'Page index for Paper Trail PNG reports. Records each page file, read/isolate range, dimensions, configured height ceiling, and any split-isolate continuation.'),
        ('read_error_profiles_page_', 'Height-bounded page of per-read Phred quality, trim windows, and internally masked low-quality regions.'),
        ('trace_chromatograms_page_', 'Height-bounded page of AB1 dye-channel traces, retained windows, and mixed-peak markers.'),
        ('assembly_overview_page_', 'Height-bounded page aligning every read to consensus coordinates, with read contribution, assembly status, and QC decision.'),
        ('assembly_read_placements.tsv', 'Auditable per-read consensus coordinates and contribution status used by the assembly overview figures.'),
        ('pairwise_pident.tsv', 'Long-form pairwise MSA percent identity for every displayed leaf pair. Identity is identical A/C/G/T bases divided by jointly unambiguous A/C/G/T MSA columns; ComparableACGTColumns reports the overlap denominator. Terminal gaps and ambiguous bases are excluded.'),
        ('clade_', 'Local phylogenetic-neighbourhood image with full sequence IDs, recommendation role/rank labels, MSA pident to P1, dataset/taxonomy context, forced nearest-baseline markers, already-sequenced markers, and primary/backup recommendation stars.'),
        ('clusters.csv', 'Consolidated cluster-membership table containing every cluster and member isolate, representative status, backup rank, taxonomy, novelty, crowding, phylogenetic isolation, and placement fields.'),
        ('cluster_summary.tsv', 'One row per cluster with consensus taxonomy, novelty/crowding summaries, phylogenetic isolation, investigation score, and backup availability.'),
        ('mwl_matches.tsv', 'Most Wanted List match report. Contains sequences whose GTDB taxonomy matched an MWL taxon, with matched rank, MWL score, evaluation score, and functional role.'),
        ('taxonomy_input_warnings.tsv', 'Warnings about inconsistencies between the classifier reference FASTA and supplied taxonomy table.'),
        ('tree_build_warnings.tsv', 'Warnings about weak phylogenetic signal, missing anchors, or poor alignment quality.'),
        ('tree_orientation_summary.tsv', 'Sequence-level audit of tree-input orientation checks. Reports which sequences were kept forward, reverse-complemented, or lacked orientation evidence before alignment.'),
        ('placement_warnings.tsv', 'Warnings about low-support placements, low identity matches, or potentially artefactual novelty assignments.'),
        ('taxa_assignments_classout.tsv', 'Synthetic classification-like TSV created when --taxa-assignments provided. Columns: id\tbest\tidentity\ttaxon\tconfidence.'),
        ('dataset_membership.itol', 'iTOL DATASET_COLORSTRIP mapping sequence IDs to dataset colours (membership).'),
        ('novelty.itol', 'iTOL colour strip showing novelty (nearest identity) for run sequences.'),
        ('filing_cabinet_dataset.itol', 'iTOL colour strip for the Filing Cabinet dataset.'),
        ('itol_dataset_membership.itol', 'iTOL DATASET_COLORSTRIP mapping sequence IDs to dataset colours (membership).'),
        ('itol_novelty.itol', 'iTOL colour strip showing novelty (nearest identity) for run sequences.'),
        ('.nwk', 'Newick tree file (phylogenetic tree). Commonly named current_tree.nwk.'),
        ('.itol', 'iTOL dataset file (text format) describing colours/strips/legends for visualisation in iTOL.'),
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


def _organise_run_outputs(outdir: str, *, primary_db_name: str | None = None):
    """Clean and organise Performance Review outputs into report folders."""
    out = Path(outdir)
    if not out.exists():
        return

    # Remove transient classifier outputs that are represented in final reports.
    for patt in (
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
        'clusters.csv',
        'backup_candidates.tsv',
        'placement_warnings.tsv',
        'novelty_metrics.tsv',
        'selection_summary.tsv',
        'sequencing_sets.tsv',
        'rumen_functions_draft.tsv',
        'qc.stats',
        'qc_rejections.tsv',
        'partner_metadata_warnings.tsv',
    ):
        _move_if_exists(out, name, 'assessment')
    for name in (
        'marker_qc_provenance.tsv', 'chimera_screen.tsv', 'chimera_uchime.tsv',
        'chimera_flagged.fasta', 'chimera_passed.fasta',
    ):
        _move_if_exists(out, name, 'quality')
    _replace_path(out / 'neighbourhoods', out / 'assessment' / 'neighbourhoods')
    _move_if_exists(out, 'mwl_matches.tsv', 'assessment')

    # Direct baseline/provided-dataset nearest-hit reports.
    _move_if_exists(out, 'baseline_hits.tsv', 'baseline')
    _move_if_exists(out, 'novelty.tsv', 'assessment', 'nearest_project_hits_raw.tsv')

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
        'current_tree_labelled.nwk',
        'current_alignment.fasta',
        'tree_build_warnings.tsv',
        'tree_orientation_summary.tsv',
    ):
        _move_if_exists(out, name, 'tree')
    _move_glob(out, 'tree_sequences_phylum_*.fasta', 'tree')

    # One iTOL metadata dataset per metadata type.
    itol_renames = {
        'itol_phylum_colours.itol': 'phylum.itol',
        'itol_family_colours.itol': 'family.itol',
        'itol_genus_colours.itol': 'genus.itol',
        'itol_dataset_membership.itol': 'dataset_membership.itol',
        'itol_baseline_tier.itol': 'baseline_tier.itol',
        'itol_novelty.itol': 'novelty.itol',
        'itol_user_colours.itol': 'user_colours.itol',
        'filing_cabinet_dataset.itol': 'filing_cabinet_dataset.itol',
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
    for name in ('user_id_map.tsv', 'filing_cabinet_id_map.tsv'):
        _move_if_exists(out, name, 'ids')

    # Inline Filing Cabinet internals: keep useful reports, drop raw matches/ref cache.
    baseline_filing_cabinet = out / 'filing_cabinet_baseline'
    if baseline_filing_cabinet.exists():
        _replace_path(baseline_filing_cabinet / 'baseline_id_map.tsv', out / 'ids' / 'baseline_id_map.tsv')
        _replace_path(
            baseline_filing_cabinet / 'baseline_combined_taxonomy.tsv',
            out / 'baseline' / 'loaded_baseline_taxonomy.tsv',
        )
        _replace_path(
            baseline_filing_cabinet / 'taxonomy.tsv',
            out / 'baseline' / 'loaded_baseline_classification.tsv',
        )
        for src in sorted(baseline_filing_cabinet.glob('filing_cabinet_*_sequences.fasta')):
            _replace_path(src, out / 'baseline' / src.name)
        try:
            shutil.rmtree(baseline_filing_cabinet)
        except Exception:
            pass

    for name in (
        'qc.fasta',
        'current_dataset_sequences.fasta',
        'current_dataset_collapsed.fasta',
        'collapsed_map.tsv',
        'collapsed_members.tsv',
        'submitted_sequences.fasta',
        'project_collection_reference.fasta',
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
    got = normalise_domain_query(str(value))
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


def _load_performance_review_baseline(
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
    baseline_tier = _normalise_cli_baseline_tier(
        baseline_dataset, getattr(args, 'baseline_tier', None),
    )
    db.upsert_dataset_role(
        baseline_dataset, 'baseline', genomes_available=True,
        baseline_tier=baseline_tier,
    )
    baseline_out = Path(outdir) / 'filing_cabinet_baseline'
    baseline_out.mkdir(parents=True, exist_ok=True)

    log.info(
        "[BASELINE] Loading baseline FASTA %s into dataset=%s before evaluating %s",
        baseline_fasta,
        baseline_dataset,
        run_dataset or '(current run)',
    )
    alias_entries, mapped_fasta = db.register_filing_cabinet(
        baseline_fasta,
        taxa_tsv=None,
        colour_csv=getattr(args, 'baseline_colours', None),
        source='baseline',
        dataset=baseline_dataset,
        outdir=str(baseline_out),
        shorten_ids=bool(getattr(args, 'baseline_shorten_ids', False)),
        baseline_tier=baseline_tier,
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
            raise RuntimeError(f'Failed to read baseline taxonomy assignments {baseline_assignment_tsv}: {e}') from e
        if not tax_entries:
            raise RuntimeError(
                f'No baseline taxonomy rows mapped to {baseline_dataset}; check baseline FASTA/table IDs.'
            )
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
                raise RuntimeError(f'Baseline classification failed: {e}') from e
    else:
        log.info("[BASELINE] --baseline-skip-classify set; baseline loaded without taxonomy classification")

    baseline_domain = _sequence_domain_filter(args)
    if baseline_domain:
        try:
            _prune_dataset_by_kingdom(db, baseline_dataset, baseline_domain, '[BASELINE]')
        except Exception as e:
            raise RuntimeError(f'Baseline domain filtering failed: {e}') from e

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
        evidence = []
        for qid, row in results.items():
            if not isinstance(row, (tuple, list)) or len(row) < 11:
                continue
            evidence.append({
                'sequence_id': qid, 'ref_db': db_name, 'best_hit': row[0], 'identity': row[1],
                'query_coverage': row[4], 'target_coverage': row[5], 'alignment_length': row[6],
                'query_length': row[7], 'target_length': row[8], 'mismatches': row[9], 'gaps': row[10],
            })
        db.upsert_classification_evidence(evidence)


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


def _normalise_cli_baseline_tier(dataset: str, explicit: str | None = None) -> str:
    value = str(explicit or '').strip().lower()
    if value in {'primary', 'hungate'}:
        value = 'priority'
    if value in {'priority', 'secondary'}:
        return value
    dataset_name = str(dataset or '').strip().lower()
    return 'priority' if 'hungate' in dataset_name else 'secondary'


def _resolve_novelty_baseline_datasets(db, args) -> list[str]:
    """Return every registered cultured baseline plus explicit current-run additions."""
    roles = db.get_dataset_roles()
    names = [
        name for name, details in roles.items()
        if details.get('role') == 'baseline'
    ]
    default_name = getattr(args, 'baseline_dataset', None)
    if default_name and (
        getattr(args, 'baseline_fasta', None)
        or roles.get(default_name, {}).get('role') == 'baseline'
    ):
        names.append(default_name)
    names.extend(getattr(args, 'novelty_baseline_datasets', None) or [])
    return list(dict.fromkeys(str(name) for name in names if str(name).strip()))


def _cmd_filing_cabinet_impl(args):
    db = Database(args.db)
    db.initialise()
    outdir = args.out or '.'
    threads = int(getattr(args, 'threads', 4) or 4)
    requested_kingdom = _sequence_domain_filter(args)
    os.makedirs(outdir, exist_ok=True)
    _configure_logging(outdir)
    logging.getLogger(__name__).info("[FILING CABINET] Starting baseline registration in %s", args.db)
    logging.getLogger(__name__).info(
        "[FILING CABINET] Sequence-domain profile: %s",
        _sequence_domain_label(requested_kingdom),
    )
    baseline_tier = _normalise_cli_baseline_tier(
        getattr(args, 'dataset', 'Baseline'), getattr(args, 'baseline_tier', None),
    )
    db.upsert_dataset_role(
        getattr(args, 'dataset', 'Baseline'),
        'baseline',
        genomes_available=True,
        baseline_tier=baseline_tier,
    )

    alias_entries, mapped_fasta = db.register_filing_cabinet(
        args.fasta,
        taxa_tsv=getattr(args, 'taxa', None),
        colour_csv=getattr(args, 'colours', None),
        source='filing_cabinet',
        dataset=getattr(args, 'dataset', 'Baseline'),
        outdir=outdir,
        shorten_ids=bool(getattr(args, 'shorten_ids', False)),
        baseline_tier=baseline_tier,
    )
    try:
        if alias_entries:
            id_map_path = _write_id_map_tsv(Path(outdir) / 'filing_cabinet_id_map.tsv', alias_entries)
            logging.getLogger(__name__).info('[FILING CABINET] Wrote baseline ID mapping to %s', id_map_path)
    except Exception as e:
        logging.getLogger(__name__).warning('[FILING CABINET] Could not write baseline ID mapping file: %s', e)
    effective_ref, effective_taxa_tsv, assignment_tsv = _resolve_reference_inputs(
        getattr(args, 'ref', None),
        getattr(args, 'taxa', None),
        getattr(args, 'taxa_assignments', None),
        source_fasta_path=args.fasta,
        log_prefix='[FILING CABINET]',
    )
    classification_requested = bool(getattr(args, 'classify', False) or (getattr(args, 'taxa_assignments', None) and not assignment_tsv))

    # If the user provided a table of predetermined taxa assignments, use it
    # instead of running the classifier. The table should have at least ID and
    # taxonomy columns; TSV, CSV, and .gz variants are supported.
    if assignment_tsv:
        taxa_file = assignment_tsv
        logging.getLogger(__name__).info("[FILING CABINET] Using taxa assignments from %s (skipping classifier)", taxa_file)
        # build mapping orig->short from alias_entries
        orig_to_short = _build_orig_to_short(alias_entries)
        try:
            tax_entries = load_taxonomy_entries_from_assignments(
                taxa_file,
                orig_to_short,
                db,
                getattr(args, 'dataset', 'Baseline'),
                source_fasta_path=args.fasta,
            )
        except Exception as e:
            raise RuntimeError(f'Could not read Filing Cabinet taxonomy assignments {taxa_file}: {e}') from e

        if not tax_entries:
            raise RuntimeError(
                f'Filing Cabinet taxonomy assignments {taxa_file} did not map to any baseline sequence IDs'
            )
        db.insert_taxonomy(tax_entries)
        logging.getLogger(__name__).info("[FILING CABINET] Inserted/updated taxonomy for %d baseline IDs from taxa_assignments", len(tax_entries))

    # If classification requested, run classifier on the mapped fasta (short ids)
    # unless the user supplied a taxa assignments table, in which case use that
    # instead and skip running the external classifier.
    if classification_requested and not assignment_tsv:
        if not effective_ref:
            raise RuntimeError(
                'Filing Cabinet classification was requested but no reference FASTA was provided '
                'via --ref or --taxa-assignments'
            )
        else:
            input_for_classify = str(mapped_fasta) if mapped_fasta else args.fasta
            logging.getLogger(__name__).info("[FILING CABINET] Classifying baseline FASTA %s against %s", input_for_classify, effective_ref)

            alt_databases = _build_alt_databases(args)
            ref_name = getattr(args, 'ref_name', None) or _classify_derive_db_name(effective_ref)
            main_db = getattr(args, 'main_ref', None) or ref_name
            all_results: dict = {}

            if alt_databases:
                logging.getLogger(__name__).info(
                    "[FILING CABINET] Multi-database classification: primary=%s, alt=%s, main=%s",
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
            # the registered baseline sequence IDs (short ids present in the DB)
            orig_to_short = _build_orig_to_short(alias_entries)
            try:
                tax_entries, dist_entries = load_classification_results_for_dataset(
                    class_out,
                    orig_to_short,
                    db,
                    getattr(args, 'dataset', 'Baseline'),
                )

                if tax_entries:
                    db.insert_taxonomy(tax_entries)
                    logging.getLogger(__name__).info("[FILING CABINET] Inserted/updated taxonomy for %d baseline IDs", len(tax_entries))
                if dist_entries:
                    db.insert_distances(dist_entries)
                    logging.getLogger(__name__).info("[FILING CABINET] Inserted/updated distances for %d baseline IDs", len(dist_entries))
            except Exception as e:
                raise RuntimeError(f'Could not persist Filing Cabinet classification output: {e}') from e

    # Apply the active sequence-domain profile before tree building so the
    # backbone tree is domain-specific.
    try:
        _prune_dataset_by_kingdom(
            db,
            getattr(args, 'dataset', 'Baseline'),
            requested_kingdom,
            '[FILING CABINET]',
        )
    except Exception as e:
        raise RuntimeError(f'Filing Cabinet domain filtering failed: {e}') from e

    # Optionally build a baseline tree/alignment from the registered baseline sequences
    if getattr(args, 'build_tree', False):
        try:
            # prefer mapped fasta (short ids) if present
            user_fasta = str(mapped_fasta) if mapped_fasta else args.fasta
            if requested_kingdom:
                try:
                    allowed_ids = _dataset_sequence_ids(db, getattr(args, 'dataset', 'Baseline'))
                    domain_fasta = Path(outdir) / f'filing_cabinet_{requested_kingdom}_sequences.fasta'
                    kept = _filter_fasta_to_ids(user_fasta, domain_fasta, allowed_ids)
                    if kept:
                        logging.getLogger(__name__).info(
                            "[FILING CABINET] Filtered baseline tree FASTA to %d %s sequence(s): %s",
                            kept,
                            requested_kingdom,
                            domain_fasta,
                        )
                        user_fasta = str(domain_fasta)
                    else:
                        raise RuntimeError(
                            f'Filing Cabinet domain filter {requested_kingdom} left no sequences for tree building'
                        )
                except Exception as e:
                    raise RuntimeError(f'Could not prepare domain-filtered Filing Cabinet tree FASTA: {e}') from e

            # optionally collapse registered baseline sequences before building tree
            if getattr(args, 'collapse', False):
                # require classification to be run for safe taxon-based collapsing
                if not (classification_requested or assignment_tsv):
                    raise RuntimeError('--collapse requires Filing Cabinet taxonomy assignments or --classify')
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
                            cur.execute("SELECT s.id, t.taxonomy FROM sequences s LEFT JOIN taxonomy t ON s.id = t.id WHERE s.dataset = ?", (getattr(args, 'dataset', 'Baseline'),))
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
                        artefacts = collapse_fasta_within_taxa(
                            taxa_groups,
                            outdir,
                            'filing_cabinet_collapsed.fasta',
                            'filing_cabinet_collapsed_map.tsv',
                            'filing_cabinet_collapsed_members.tsv',
                            threshold=threshold,
                            threads=threads,
                            log_prefix='[FILING CABINET COLLAPSE]',
                            strict=True,
                        )
                        user_fasta = artefacts.collapsed_path
                        logging.getLogger(__name__).info("[FILING CABINET COLLAPSE] Wrote collapsed baseline FASTA %s (reps=%d)", user_fasta, len(artefacts.collapsed_records))
                    except Exception as e:
                        raise RuntimeError(f'Filing Cabinet collapse failed: {e}') from e

            logging.getLogger(__name__).info("[FILING CABINET] Building baseline tree/alignment in %s from registered baseline sequences", outdir)
            tree.initialise_or_update_tree(
                ref_fasta=effective_ref,
                user_fasta=user_fasta,
                outdir=outdir,
                db=None,
                threads=threads,
                anchor_file=getattr(args, 'anchors', None),
                tree_method=getattr(args, 'tree_method', 'fasttree'),
            )
            logging.getLogger(__name__).info("[FILING CABINET] Baseline tree/alignment written to %s", outdir)
            try:
                warning_rows = tree.collect_tree_build_warnings(user_fasta=str(user_fasta), anchor_file=getattr(args, 'anchors', None), db=None)
                warning_rows.extend(tree.summarise_alignment_quality(str(Path(outdir) / 'current_alignment.fasta')))
                if warning_rows:
                    warn_path = tree.write_tree_warning_tsv(outdir, warning_rows)
                    logging.getLogger(__name__).warning("[FILING CABINET] Tree/alignment warnings written to %s", warn_path)
            except Exception as e:
                logging.getLogger(__name__).warning("[FILING CABINET] Failed to summarise tree/alignment quality: %s", e)
        except Exception as e:
            raise RuntimeError(f'Filing Cabinet tree construction failed: {e}') from e

    # If a sequence-domain filter is requested for the Filing Cabinet, remove any registered baseline sequences
    # whose assigned taxonomy indicates they are not the requested kingdom.
    try:
        _prune_dataset_by_kingdom(
            db,
            getattr(args, 'dataset', 'Baseline'),
            requested_kingdom,
            '[FILING CABINET]',
        )
    except Exception as e:
        raise RuntimeError(f'Final Filing Cabinet domain validation failed: {e}') from e

    # Build combined taxonomy and generate iTOL colour files for the Filing Cabinet dataset.
    try:
        out_p = Path(outdir)
        combined_tax = out_p / 'filing_cabinet_combined_taxonomy.tsv'
        if not combined_tax.exists():
            try:
                with db.connect() as conn:
                    cur = conn.cursor()
                    cur.execute("SELECT s.id, t.taxonomy, t.confidence FROM sequences s LEFT JOIN taxonomy t ON s.id = t.id WHERE s.dataset = ?", (getattr(args, 'dataset', 'Baseline'),))
                    rows = cur.fetchall()
                write_combined_taxonomy_tsv(combined_tax, rows)
                logging.getLogger(__name__).info("[FILING CABINET] Wrote combined taxonomy for %d ids to %s", len(rows), combined_tax)
            except Exception as e:
                logging.getLogger(__name__).warning("[FILING CABINET] Failed to build combined taxonomy: %s", e)

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
            itol.generate_itol_colours(str(combined_tax), outdir, user_colour_csv=getattr(args, 'colours', None), id_map=id_map, tree_file=tfile, phylum_groups=getattr(args, 'group_phyla', None))
            logging.getLogger(__name__).info("[FILING CABINET] Generated iTOL colour files in %s", outdir)
        except TypeError:
            try:
                itol.generate_itol_colours(str(combined_tax), outdir, user_colour_csv=getattr(args, 'colours', None), phylum_groups=getattr(args, 'group_phyla', None))
                logging.getLogger(__name__).info("[FILING CABINET] Generated iTOL colour files in %s (no id_map)", outdir)
            except Exception as e:
                logging.getLogger(__name__).warning("[FILING CABINET] Failed to generate iTOL colours: %s", e)
        except Exception as e:
            logging.getLogger(__name__).warning("[FILING CABINET] Failed to generate iTOL colours: %s", e)

        try:
            func_tsv = getattr(args, 'functional', None)
            if func_tsv:
                try:
                    written = itol.write_functional_annotations(str(func_tsv), outdir, id_map=id_map)
                    logging.getLogger(__name__).info("[FILING CABINET] Wrote functional annotation iTOL files: %s", ','.join(written) if written else '(none)')
                except Exception as e:
                    logging.getLogger(__name__).warning("[FILING CABINET] Functional annotations generation failed: %s", e)
        except Exception:
            pass

        if getattr(args, 'draft_rumen_functions', False):
            try:
                combined_tax_path = str(out_p / 'combined_taxonomy.tsv')
                tsv_out, itol_out = itol.generate_rumen_function_draft(
                    combined_tax_path, outdir, id_map=id_map
                )
                if tsv_out:
                    logging.getLogger(__name__).info(
                        "[FILING CABINET] Draft rumen functional annotation: %s", tsv_out
                    )
                if itol_out:
                    logging.getLogger(__name__).info(
                        "[FILING CABINET] Rumen functional iTOL file: %s", itol_out
                    )
            except Exception as e:
                logging.getLogger(__name__).warning(
                    "[FILING CABINET] Draft rumen functions generation failed: %s", e
                )

        try:
            ds_colour = itol._name_to_dataset_colour(getattr(args, 'dataset', 'Baseline'))
            itol_path = out_p / 'filing_cabinet_dataset.itol'
            dataset_label = getattr(args, 'dataset', 'Baseline')
            with db.connect() as conn:
                cur = conn.cursor()
                cur.execute("SELECT id FROM sequences WHERE dataset = ?", (getattr(args, 'dataset', 'Baseline'),))
                id_to_colour = {iid: ds_colour for (iid,) in cur.fetchall()}
            itol.write_dataset_colourstrip(str(itol_path), dataset_label, id_to_colour, legend_title=f"{dataset_label} legend")
            logging.getLogger(__name__).info("[FILING CABINET] Wrote dataset iTOL colour strip to %s", itol_path)
            tier_path = out_p / 'baseline_tier.itol'
            itol.write_baseline_tier_strip(str(tier_path), list(id_to_colour), {iid: baseline_tier for iid in id_to_colour})
            logging.getLogger(__name__).info("[FILING CABINET] Wrote baseline-tier iTOL colour strip to %s", tier_path)
        except Exception as e:
            logging.getLogger(__name__).warning("[FILING CABINET] Failed to write dataset iTOL colour strip: %s", e)
    except Exception:
        pass

    # write brief explanations for files produced by the Filing Cabinet
    try:
        _write_output_explanations(outdir)
    except Exception:
        pass


def cmd_filing_cabinet(args):
    """Register the baseline Filing Cabinet without publishing a partial database."""
    import fcntl
    import sqlite3

    outdir = Path(args.out or '.').expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    original = Path(args.db).expanduser().resolve()
    original.parent.mkdir(parents=True, exist_ok=True)
    manifest = RunManifest(outdir, 'filing_cabinet')
    manifest.add_input(args.fasta, role='baseline_fasta')
    for role, source in (
        ('taxonomy', getattr(args, 'taxa', None)),
        ('taxonomy_assignments', getattr(args, 'taxa_assignments', None)),
        ('reference', getattr(args, 'ref', None)),
        ('colours', getattr(args, 'colours', None)),
        ('anchors', getattr(args, 'anchors', None)),
        ('functional_annotations', getattr(args, 'functional', None)),
    ):
        if source:
            manifest.add_input(source, role=role)
    for source in getattr(args, 'alt_ref', None) or []:
        manifest.add_input(source, role='alternate_reference')
    for source in getattr(args, 'alt_taxa', None) or []:
        manifest.add_input(source, role='alternate_reference_taxonomy')

    staged_path = outdir / '.filing_cabinet_project.sqlite'
    lock_path = original.with_name(original.name + '.lock')
    with open(lock_path, 'a+') as lock_handle:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        try:
            staged_path.unlink(missing_ok=True)
            if original.is_file():
                manifest.add_input(original, role='project_database_before')
                with sqlite3.connect(str(original)) as source, sqlite3.connect(str(staged_path)) as target:
                    source.backup(target)
            else:
                Database(str(staged_path)).initialise()

            staged_args = argparse.Namespace(**vars(args))
            staged_args.db = str(staged_path)
            manifest.add_stage('filing_cabinet', 'RUNNING')
            _cmd_filing_cabinet_impl(staged_args)

            staged_db = Database(str(staged_path))
            staged_db.initialise()
            with staged_db.connect() as conn:
                integrity = conn.execute('PRAGMA integrity_check').fetchone()[0]
                baseline_count = conn.execute(
                    'SELECT COUNT(*) FROM sequences WHERE dataset = ?',
                    (getattr(args, 'dataset', 'Baseline'),),
                ).fetchone()[0]
            if integrity != 'ok':
                raise RuntimeError(f'staged Filing Cabinet database failed integrity_check: {integrity}')
            if baseline_count < 1:
                raise RuntimeError('Filing Cabinet registration produced no baseline sequences')

            for role, path, required in (
                ('id_map', outdir / 'filing_cabinet_id_map.tsv', False),
                ('combined_taxonomy', outdir / 'filing_cabinet_combined_taxonomy.tsv', True),
                ('alignment', outdir / 'current_alignment.fasta', bool(getattr(args, 'build_tree', False))),
                ('tree', outdir / 'current_tree.nwk', bool(getattr(args, 'build_tree', False))),
            ):
                manifest.add_output(path, role=role, required=required)
            missing = manifest.verify_required_outputs()
            if missing:
                raise RuntimeError('required Filing Cabinet outputs are missing: ' + ', '.join(missing))

            manifest.add_stage('filing_cabinet', 'COMPLETE', detail=f'{baseline_count} baseline sequences')
            manifest.finish('COMPLETE')
            staged_db.record_project_run(
                f'filing-cabinet:{manifest.data["started_at"]}', 'filing_cabinet', 'COMPLETE',
                dataset=getattr(args, 'dataset', 'Baseline'), manifest_path=str(manifest.json_path),
                started_at=manifest.data['started_at'], completed_at=manifest.data['completed_at'],
            )

            publish_tmp = original.with_name(original.name + '.publishing')
            publish_tmp.unlink(missing_ok=True)
            with sqlite3.connect(str(staged_path)) as source, sqlite3.connect(str(publish_tmp)) as target:
                source.backup(target)
                if target.execute('PRAGMA integrity_check').fetchone()[0] != 'ok':
                    raise RuntimeError('published Filing Cabinet database failed integrity_check')
            os.replace(publish_tmp, original)
            staged_path.unlink(missing_ok=True)
            manifest.add_output(original, role='project_database')
        except BaseException as exc:
            manifest.finish('FAILED', error=exc)
            raise
        finally:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)
    print(
        f'[filing-cabinet] Registered {baseline_count} baseline sequence(s).\n'
        f'  Project database: {original}\n'
        f'  Backbone output : {outdir}'
    )


def _cmd_performance_review_impl(args):
    db = Database(args.db)
    db.initialise()
    outdir = args.out
    threads = int(getattr(args, 'threads', 4) or 4)
    os.makedirs(outdir, exist_ok=True)
    _configure_logging(outdir)
    logging.getLogger(__name__).info("[PERFORMANCE REVIEW] Starting run pipeline (input=%s)", args.input)
    domain_filter = _sequence_domain_filter(args)
    logging.getLogger(__name__).info(
        "[PERFORMANCE REVIEW] Sequence-domain profile: %s",
        _sequence_domain_label(domain_filter),
    )
    effective_ref, effective_taxa_tsv, assignment_tsv = _resolve_reference_inputs(
        getattr(args, 'ref', None),
        getattr(args, 'taxa', None),
        getattr(args, 'taxa_assignments', None),
        source_fasta_path=args.input,
        log_prefix='[PERFORMANCE REVIEW]',
    )
    if not effective_ref:
        raise SystemExit("[PERFORMANCE REVIEW] A reference FASTA is required via --ref, or `--taxa-assignments` must point to a GTDB/reference FASTA rather than a taxonomy assignment table.")
    if getattr(args, 'command', None) == 'performance-review' and not getattr(args, 'partner_metadata', None):
        raise SystemExit(
            '[PERFORMANCE REVIEW] --partner-metadata / --sequencing-metadata is required: '
            'a cumulative CSV/TSV ledger with sequence IDs, partner acronyms, optional selection commitments, and already-sequenced status.'
        )

    _load_performance_review_baseline(args, db, outdir, effective_ref, effective_taxa_tsv, threads)
    if (
        hasattr(args, 'skip_chimera_check')
        and getattr(args, 'command', None) == 'performance-review'
    ):
        baseline_datasets = db.get_dataset_names_by_role('baseline')
        if not baseline_datasets:
            raise SystemExit(
                '[PERFORMANCE REVIEW] No cultured baseline is registered. Supply --baseline-fasta for this run '
                'or establish it first with `branchmanager filing-cabinet`.'
            )

    # QC
    qc_out = qc.run_qc(
        args.input,
        outdir,
        min_len=getattr(args, 'min_len', 800),
        max_n_percent=getattr(args, 'max_n_percent', 5.0),
    )

    # Preserve every isolate ID, including exact marker-sequence duplicates.
    # Distinct isolates can still contain different genomes or provide useful
    # DNA-extraction backups; optional tree collapsing handles visual density.
    candidate_input = qc_out

    # Map user-provided IDs to runtime IDs and insert them into the DB.
    from branchmanager.utils.fasta import read_fasta, write_fasta
    current_dataset_fasta = Path(outdir) / 'current_dataset_sequences.fasta'
    used_ids = set(db.get_all_ids())
    orig_to_short = {}
    mapped_records = []
    skipped_existing = 0

    # If the user requested a kingdom filter, classify the QC-passing input
    # before insertion so only sequences assigned to that domain are retained.
    early_class_out = None
    allowed_qids = None
    kingdom = domain_filter
    if kingdom:
        kingdom_text = str(kingdom)
        if not effective_ref:
            logging.getLogger(__name__).warning("[PERFORMANCE REVIEW] --sequence-domain specified but no reference FASTA was available; cannot classify to filter; proceeding without domain filtering")
            allowed_qids = None
        else:
            try:
                logging.getLogger(__name__).info("[PERFORMANCE REVIEW] Running pre-insert classification to filter by domain=%s", kingdom)
                early_class_out = classify.run_classification(str(candidate_input), outdir, ref_fasta=effective_ref, taxa_tsv=effective_taxa_tsv, threads=threads)
                allowed_qids = _classification_ids_matching_kingdom(early_class_out, kingdom_text)
                logging.getLogger(__name__).info("[PERFORMANCE REVIEW] Domain filter: %d QC-passing sequences match %s", len(allowed_qids), kingdom)
            except Exception as e:
                raise SystemExit(f'[PERFORMANCE REVIEW] Domain classification/filtering failed: {e}')

    for h, s in read_fasta(candidate_input):
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
        # want to include it in the tree for visualisation, but we use the
        # existing DB ID so downstream taxonomy/metadata resolve consistently.
        # The "INSERT OR IGNORE" during insert_sequences will prevent DB duplicates.
        existing_id = db.resolve_sequence_id(h)
        if existing_id and existing_id in used_ids:
            try:
                db.assert_sequence_compatible(existing_id, s)
            except ValueError as e:
                raise SystemExit(f'[DB] {e}')
            skipped_existing += 1
            short = existing_id
            orig_to_short[h] = short
            try:
                cid = canonicalise_sequence_id(h)
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
                raise SystemExit(f"[PERFORMANCE REVIEW] {e}")
            mapped_records.append((short, s))
            orig_to_short[h] = short
            try:
                cid = canonicalise_sequence_id(h)
                if cid:
                    orig_to_short[cid] = short
            except Exception:
                pass
    if skipped_existing:
        logging.getLogger(__name__).info("[DB] Found %d sequences already in DB; keeping them for tree inclusion", skipped_existing)
    write_fasta(mapped_records, str(current_dataset_fasta))
    logging.getLogger(__name__).info("[DB] Mapped %d user sequence IDs to runtime IDs and wrote %s", len(mapped_records), current_dataset_fasta)
    if not mapped_records:
        raise SystemExit(
            '[PERFORMANCE REVIEW] No sequences passed sequence/domain QC. Project state was not accepted; '
            'inspect qc.stats, classification, and the requested --sequence-domain.'
        )
    # insert mapped records into DB (dataset provided by user)
    run_dataset = getattr(args, 'dataset', 'user')
    db.upsert_dataset_role(run_dataset, 'candidate', genomes_available=False)
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

    # Carry Paper Trail/Merge Meeting quality and raw-read provenance into the persistent
    # project record. Direct FASTA inputs are intentionally QUALITY_UNVERIFIED
    # unless the user supplies an explicit acceptance or review decision.
    try:
        from branchmanager.marker_provenance import load_marker_provenance, write_marker_qc_bridge
        marker_provenance, marker_qc_source = load_marker_provenance(
            args.input,
            qc_path=getattr(args, 'marker_qc', None),
            review_path=getattr(args, 'marker_review', None),
            id_map=orig_to_short,
            accept_unverified=bool(getattr(args, 'accept_unverified_marker_qc', False)),
        )
        mapped_id_set = {short for short, _sequence in mapped_records}
        marker_provenance = [
            row for row in marker_provenance if row.get('sequence_id') in mapped_id_set
        ]
        db.upsert_sequence_provenance(marker_provenance)
        marker_bridge_path = write_marker_qc_bridge(
            Path(outdir) / 'marker_qc_provenance.tsv', marker_provenance,
        )
        existing_statuses = db.get_isolate_statuses(mapped_id_set)
        terminal_states = {
            'SAB_APPROVED', 'DNA_EXTRACTION_PENDING', 'DNA_EXTRACTION_FAILED',
            'LIBRARY_PENDING', 'SEQUENCED', 'GENOME_QC_FAILED', 'GENOME_QC_PASSED', 'WITHDRAWN',
        }
        for row in marker_provenance:
            sequence_id = row['sequence_id']
            if existing_statuses.get(sequence_id, {}).get('status') in terminal_states:
                continue
            qc_class = str(row.get('marker_qc_class') or '').upper()
            review = str(row.get('manual_review_status') or '').upper()
            status = 'MARKER_QC_PASSED' if qc_class == 'PASS_HIGH_CONFIDENCE' or review == 'APPROVED' else 'TRACE_REVIEW'
            db.update_isolate_status(
                sequence_id, status,
                detail=str(row.get('marker_qc_reasons') or ''),
                source_file=marker_qc_source or str(Path(args.input).resolve()),
            )
        logging.getLogger(__name__).info(
            '[PERFORMANCE REVIEW] Marker provenance for %d sequence(s) written to %s',
            len(marker_provenance), marker_bridge_path,
        )
    except Exception as exc:
        raise SystemExit(f'[PERFORMANCE REVIEW] Marker-QC provenance could not be established: {exc}')

    try:
        if not hasattr(args, 'skip_chimera_check') or getattr(args, 'skip_chimera_check', False):
            chimera_results = {
                sequence_id: {'call': 'SKIPPED', 'score': None}
                for sequence_id, _sequence in mapped_records
            }
            logging.getLogger(__name__).warning(
                '[PERFORMANCE REVIEW] Chimera screening explicitly skipped; affected marker evidence will require review.'
            )
        else:
            from branchmanager.pipeline import chimera as _chimera
            chimera_report, chimera_results = _chimera.run_reference_screen(
                str(current_dataset_fasta),
                getattr(args, 'chimera_ref', None) or effective_ref,
                outdir,
                threads=threads,
            )
            logging.getLogger(__name__).info('[PERFORMANCE REVIEW] Chimera report written to %s', chimera_report)
        for row in marker_provenance:
            result = chimera_results.get(row['sequence_id'], {'call': 'INDETERMINATE', 'score': None})
            row['chimera_call'] = result.get('call', 'INDETERMINATE')
            row['chimera_score'] = result.get('score')
        db.upsert_sequence_provenance(marker_provenance)
        write_marker_qc_bridge(Path(outdir) / 'marker_qc_provenance.tsv', marker_provenance)
    except Exception as exc:
        raise SystemExit(f'[PERFORMANCE REVIEW] Reference-based chimera screening failed: {exc}')

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
        logging.getLogger(__name__).info("[PERFORMANCE REVIEW] Using taxa assignments from %s (skipping classifier)", taxa_file)
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
            raise RuntimeError(f'Failed to read taxonomy assignments {taxa_file}: {e}') from e
        if not tax_entries_local:
            raise RuntimeError('No taxonomy assignment rows mapped to the current Performance Review IDs.')

        # persist taxonomy for mapped ids
        if tax_entries_local:
            try:
                db.insert_taxonomy(tax_entries_local)
                logging.getLogger(__name__).info("[DB] Inserted/updated taxonomy for %d ids from taxa_assignments", len(tax_entries_local))
            except Exception as e:
                raise RuntimeError(f'Failed to store provided taxonomy assignments: {e}') from e

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
                except Exception as exc:
                    raise RuntimeError(f'Failed to normalise provided taxonomy assignments: {exc}') from exc
            class_out = str(class_out_path)
        except Exception as exc:
            raise RuntimeError(f'Failed to create normalised taxonomy evidence: {exc}') from exc
    else:
        # Only run classification here if we did not already run an early
        # pre-insert classification for kingdom filtering and the user did
        # not supply taxa_assignments.
        alt_databases = _build_alt_databases(args)
        ref_name = getattr(args, 'ref_name', None) or _classify_derive_db_name(effective_ref)
        run_main_db_name = getattr(args, 'main_ref', None) or ref_name

        if alt_databases:
            if early_class_out:
                logging.getLogger(__name__).info(
                    "[PERFORMANCE REVIEW] Domain-filter preclassification will not be reused for final output; "
                    "running multi-database classification so alternate taxonomy is persisted."
                )
            logging.getLogger(__name__).info(
                "[PERFORMANCE REVIEW] Multi-database classification: primary=%s, alt=%s, main=%s",
                ref_name, [n for _, _, n in alt_databases], run_main_db_name,
            )
            class_out, all_class_results = classify.run_all_classifications(
                str(current_dataset_fasta), outdir,
                primary_ref=effective_ref,
                primary_taxa=effective_taxa_tsv,
                primary_name=ref_name,
                alt_refs=alt_databases,
                threads=threads,
                main_db=run_main_db_name,
            )
            _store_alt_taxonomy_in_db(db, all_class_results, run_main_db_name)
        elif early_class_out:
            class_out = early_class_out
        else:
            class_out = classify.run_classification(
                str(current_dataset_fasta), outdir, ref_fasta=effective_ref, taxa_tsv=effective_taxa_tsv, threads=threads
            )

    if (
        getattr(args, 'command', None) == 'performance-review'
        and domain_filter in ('bacteria', 'archaea')
        and 'gtdb' not in str(run_main_db_name).lower()
    ):
        raise SystemExit(
            '[PERFORMANCE REVIEW] bacterial/archaeal assessment requires GTDB as the authoritative database. '
            'Name the GTDB reference with --ref-name GTDB or select it with --main-ref GTDB. '
            'Other databases remain available through --alt-ref for reporting only.'
        )

    # Persist authoritative taxonomy before rolling novelty/genome coverage so
    # current-run candidates can be counted by GTDB species immediately.
    try:
        tax_entries, dist_entries = load_classification_results_for_dataset(
            class_out or '',
            orig_to_short,
            db,
            run_dataset,
        )
    except Exception as e:
        raise SystemExit(f'[PERFORMANCE REVIEW] Authoritative classification output could not be parsed: {e}')
    if tax_entries:
        db.insert_taxonomy(tax_entries)
        logging.getLogger(__name__).info("[DB] Inserted/updated taxonomy for %d ids", len(tax_entries))
    if dist_entries:
        db.insert_distances(dist_entries)
        logging.getLogger(__name__).info("[DB] Inserted/updated distances for %d ids", len(dist_entries))
    classification_evidence = []
    for row in iter_classification_rows(class_out or ''):
        query_id = str(row.get('qid') or '')
        mapped = orig_to_short.get(query_id) or (query_id if query_id in {item[0] for item in mapped_records} else None)
        if not mapped:
            continue
        classification_evidence.append({
            'sequence_id': mapped, 'ref_db': run_main_db_name,
            'best_hit': row.get('best'), 'identity': row.get('identity'),
            'query_coverage': row.get('query_coverage'), 'target_coverage': row.get('target_coverage'),
            'alignment_length': row.get('alignment_length'), 'query_length': row.get('query_length'),
            'target_length': row.get('target_length'), 'mismatches': row.get('mismatches'),
            'gaps': row.get('gaps'),
        })
    db.upsert_classification_evidence(classification_evidence)

    # Register explicit cultured baselines, then migrate unlabelled historical
    # run datasets into the rolling candidate collection.
    novelty_baseline_datasets = _resolve_novelty_baseline_datasets(db, args)
    for baseline_name in novelty_baseline_datasets:
        db.upsert_dataset_role(
            baseline_name, 'baseline', genomes_available=True,
            baseline_tier=_normalise_cli_baseline_tier(baseline_name),
        )
    try:
        registered_roles = db.get_dataset_roles()
        with db.connect() as conn:
            datasets_in_db = [row[0] for row in conn.execute(
                "SELECT DISTINCT dataset FROM sequences WHERE dataset IS NOT NULL AND dataset != ''"
            ).fetchall()]
        for dataset_name in datasets_in_db:
            if dataset_name not in registered_roles:
                db.upsert_dataset_role(dataset_name, 'candidate', genomes_available=False)
    except Exception as e:
        raise RuntimeError(f'Could not establish rolling dataset roles: {e}') from e

    # Assess only this submission. Novelty searches still use every registered
    # candidate and baseline in the database as rolling comparison context.
    metrics_query_fasta = str(current_dataset_fasta)

    # novelty
    target_fasta = getattr(args, 'target', None)
    novelty_out = novelty.run_novelty(str(current_dataset_fasta), effective_ref, outdir, db=db, run_dataset=run_dataset, threads=threads, target_fasta=target_fasta)
    try:
        novelty_metrics_out = novelty.build_reference_novelty_metrics(
            metrics_query_fasta,
            effective_ref,
            outdir,
            threads=threads,
            db=db,
            run_dataset=run_dataset,
            target_fasta=target_fasta,
            baseline_datasets=novelty_baseline_datasets,
            pangenome_target=getattr(
                args, 'pangenome_target', _selection_sets.DEFAULT_PANGENOME_TARGET,
            ),
        )
        logging.getLogger(__name__).info("[NOVELTY] Wrote novelty metrics to %s", novelty_metrics_out)
    except Exception as e:
        raise SystemExit(f'[PERFORMANCE REVIEW] Required novelty metrics failed: {e}')

    # Safety-net: remove any run sequences with explicit non-matching kingdom assignments.
    try:
        _prune_dataset_by_kingdom(db, run_dataset, str(kingdom) if kingdom else None, '[PERFORMANCE REVIEW]')
    except Exception as e:
        raise RuntimeError(f'Final Performance Review domain filtering failed: {e}') from e

    # Initialise cluster-tracking variables used by both the tree section and
    # the assessment section below.  They will be populated if --collapse is
    # active; otherwise they stay empty and the assessment omits cluster columns.
    tree_fasta = current_dataset_fasta
    run_member_to_rep: dict = {}
    run_rep_to_members: dict = {}

    # update tree/alignment
    try:
        # report DB sequence count for diagnostics so users can verify that
        # registered baseline sequences exist and will be used to seed the backbone
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
                for h, s in read_fasta(str(current_dataset_fasta)):
                    tax = qid_to_tax.get(h)
                    taxa_groups.setdefault(tax, []).append((h, s))
            except Exception:
                taxa_groups = {None: []}
            try:
                artefacts = collapse_fasta_within_taxa(
                    taxa_groups,
                    outdir,
                    'current_dataset_collapsed.fasta',
                    'collapsed_map.tsv',
                    'collapsed_members.tsv',
                    threshold=threshold,
                    threads=threads,
                    log_prefix='[COLLAPSE]',
                    strict=True,
                )
                tree_fasta = artefacts.collapsed_path
                run_member_to_rep = artefacts.member_to_rep or {}
                # Build rep -> members mapping
                for mem, rep in run_member_to_rep.items():
                    run_rep_to_members.setdefault(rep, []).append(mem)
                logging.getLogger(__name__).info("[COLLAPSE] Wrote collapsed fasta %s (reps=%d)", tree_fasta, len(artefacts.collapsed_records))
            except Exception as e:
                raise RuntimeError(f'Requested Performance Review collapse failed: {e}') from e

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
                                            # Normalise for comparison (case-insensitive, strip underscores)
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
                        raise RuntimeError(
                            f"No tree sequences remained for requested phylum '{phylum_filter}'"
                        )
                else:
                    raise RuntimeError(
                        f"Could not find requested phylum '{phylum_filter}' in classification output"
                    )
            except Exception as e:
                raise RuntimeError(f'Phylum filtering failed: {e}') from e

        # pass the Database object so the tree builder can export existing
        # registered baseline sequences from the DB to form the backbone alignment
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
            previous_review=None if _force_rebuild else getattr(args, 'previous_review', None),
            force_rebuild=_force_rebuild,
            anchor_file=getattr(args, 'anchors', None),
            tree_method=getattr(args, 'tree_method', 'fasttree'),
        )
        try:
            warning_rows = tree.collect_tree_build_warnings(user_fasta=str(tree_fasta), anchor_file=getattr(args, 'anchors', None), db=db, db_dataset=None)
            warning_rows.extend(tree.summarise_alignment_quality(str(Path(outdir) / 'current_alignment.fasta')))
            if warning_rows:
                warn_path = tree.write_tree_warning_tsv(outdir, warning_rows)
                logging.getLogger(__name__).warning("[TREE] Tree/alignment warnings written to %s", warn_path)
        except Exception as e:
            logging.getLogger(__name__).warning("[TREE] Failed to summarise tree/alignment quality: %s", e)
    except Exception as e:
        raise SystemExit(f'[PERFORMANCE REVIEW] Required MSA/tree update failed: {e}')

    # Build combined taxonomy for iTOL from DB sequences and current run results.
    combined_path = Path(outdir) / 'combined_taxonomy.tsv'
    merged = {}
    order = []

    try:
        previous_review = getattr(args, 'previous_review', None)
        previous_review_ids = None
        if previous_review:
            try:
                cand = next(
                    (path for path in _combined_taxonomy_candidates(previous_review) if path.exists()),
                    None,
                )
                if cand is not None:
                    previous_review_ids = read_combined_taxonomy_ids(cand)
            except Exception:
                previous_review_ids = None

        base_rows = collect_db_taxonomy_rows(db, previous_review_ids if previous_review_ids else None)
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

    # Load previous-review ID map for iTOL when --previous-review is provided.
    id_map_for_itol = None
    previous_review = getattr(args, 'previous_review', None)
    if previous_review:
        try:
            id_map = {}
            cand_map = _find_preferred_id_map(str(previous_review))
            if cand_map:
                id_map = _load_id_map_from_tsv(cand_map, db=db)
            if id_map:
                id_map_for_itol = id_map
        except Exception:
            id_map_for_itol = None

        try:
            tree_path = Path(outdir) / 'current_tree.nwk'
            tfile = str(tree_path) if tree_path.exists() else _find_tree_file_in_dir(outdir)
            itol.generate_itol_colours(str(combined_path), outdir, user_colour_csv=getattr(args, 'user_colours', None), id_map=id_map_for_itol, tree_file=tfile, phylum_groups=getattr(args, 'group_phyla', None))
            logging.getLogger(__name__).info("[ITOL] Generated iTOL colour files in %s", outdir)
        except Exception as e:
            logging.getLogger(__name__).warning("[ITOL] Failed to generate iTOL files: %s", e)
    # Generate iTOL colour files for the Performance Review output directory.
    try:
        tree_path = Path(outdir) / 'current_tree.nwk'
        tfile = str(tree_path) if tree_path.exists() else _find_tree_file_in_dir(outdir)
        itol.generate_itol_colours(str(combined_path), outdir, user_colour_csv=getattr(args, 'user_colours', None), id_map=id_map_for_itol, tree_file=tfile, phylum_groups=getattr(args, 'group_phyla', None))
        logging.getLogger(__name__).info("[ITOL] Generated iTOL colour files in %s", outdir)
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
    # produce dataset membership band (Filing Cabinet versus review batch)
    try:
        combined_tax = combined_path
        ids_in_order = []
        if combined_tax.exists():
            with open(combined_tax) as ct:
                next(ct, None)
                for line in ct:
                    iid = line.strip().split('\t')[0]
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
        tier_map = {}
        try:
            with db.connect() as conn:
                cur = conn.cursor()
                placeholders = ','.join('?' for _ in ids_in_order) if ids_in_order else ''
                if placeholders:
                    cur.execute(
                        "SELECT s.id, COALESCE(NULLIF(s.baseline_tier, ''), r.baseline_tier), COALESCE(r.role, '') "
                        "FROM sequences s LEFT JOIN dataset_roles r ON r.dataset = s.dataset "
                        f"WHERE s.id IN ({placeholders})",
                        tuple(ids_in_order),
                    )
                    for iid_row, tier, role in cur.fetchall():
                        if role == 'baseline' and str(tier or '').strip().lower() in {'priority', 'secondary'}:
                            tier_map[iid_row] = str(tier).strip().lower()
        except Exception:
            tier_map = {}
        if tier_map:
            tier_path = Path(outdir) / 'itol_baseline_tier.itol'
            itol.write_baseline_tier_strip(str(tier_path), ids_in_order, tier_map)
            logging.getLogger(__name__).info("[ITOL] Wrote baseline-tier ITOL to %s", tier_path)
    except Exception as e:
        logging.getLogger(__name__).warning("[ITOL] Failed to write dataset membership ITOL: %s", e)

    # Add a novelty-gradient colour strip for new sequences.
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
                str(current_dataset_fasta),
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
                logging.getLogger(__name__).warning("[PERFORMANCE REVIEW] Placement warnings written to %s", warn_path)
            try:
                run_ids_for_assessment = [record[0] for record in mapped_records]
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
                            "[PERFORMANCE REVIEW] Found alt-db taxonomy for ref_dbs: %s", alt_ref_dbs
                        )
                except Exception as e:
                    logging.getLogger(__name__).warning("[PERFORMANCE REVIEW] Could not load alt-db taxonomy: %s", e)

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
                from branchmanager.marker_provenance import marker_qc_flag
                project_provenance = db.get_sequence_provenance(run_ids_for_assessment)
                for assessment_row in assessment_rows:
                    provenance = project_provenance.get(str(assessment_row.get('id'))) or {}
                    qc_flag = marker_qc_flag(provenance)
                    assessment_row.update({
                        'marker_qc_class': provenance.get('marker_qc_class', 'QUALITY_UNVERIFIED'),
                        'marker_qc_recommendation': provenance.get('marker_qc_recommendation', 'MANUAL_REVIEW'),
                        'marker_qc_reasons': provenance.get('marker_qc_reasons', 'no_marker_qc_provenance'),
                        'marker_manual_review_status': provenance.get('manual_review_status', 'NOT_REVIEWED'),
                        'marker_qc_flag': qc_flag or 'MARKER_QC_PASSED',
                        'marker_source_manifest': provenance.get('source_manifest', ''),
                        'chimera_call': provenance.get('chimera_call', 'NOT_RUN'),
                        'chimera_score': provenance.get('chimera_score', 'NA'),
                    })
                    if qc_flag:
                        flags = [flag for flag in str(assessment_row.get('placement_flags') or '').split(';') if flag]
                        if qc_flag not in flags:
                            flags.append(qc_flag)
                        assessment_row['placement_flags'] = ';'.join(flags)
                    chimera_call = str(provenance.get('chimera_call') or 'NOT_RUN').upper()
                    if chimera_call != 'PASS':
                        chimera_flag = {
                            'CHIMERA': 'CHIMERA_DETECTED',
                            'INDETERMINATE': 'CHIMERA_INDETERMINATE',
                            'SKIPPED': 'CHIMERA_CHECK_SKIPPED',
                        }.get(chimera_call, 'CHIMERA_NOT_RUN')
                        flags = [flag for flag in str(assessment_row.get('placement_flags') or '').split(';') if flag]
                        if chimera_flag not in flags:
                            flags.append(chimera_flag)
                        assessment_row['placement_flags'] = ';'.join(flags)
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
                        raise RuntimeError(f'MWL annotation failed: {e}') from e
                _tree_nwk = str(Path(outdir) / 'current_tree.nwk')
                try:
                    _cluster_summary, _clusters_csv, _backup_tsv = _cluster_report.generate_cluster_reports(
                        outdir=outdir,
                        assessment_rows=assessment_rows,
                        tree_path=_tree_nwk if Path(_tree_nwk).exists() else None,
                    )
                    if mwl_entries:
                        try:
                            _mwl.add_evaluation_scores(assessment_rows)
                            _mwl.write_mwl_matches_tsv(Path(outdir) / 'mwl_matches.tsv', assessment_rows)
                        except Exception as _mwe:
                            logging.getLogger(__name__).warning("[MWL] Failed to refresh MWL evaluation scores: %s", _mwe)
                    if _cluster_summary:
                        logging.getLogger(__name__).info(
                            "[CLUSTER] Wrote cluster summary → %s",
                            _cluster_summary,
                        )
                    if _clusters_csv:
                        logging.getLogger(__name__).info(
                            "[CLUSTER] Wrote consolidated cluster membership → %s", _clusters_csv,
                        )
                    if _backup_tsv:
                        logging.getLogger(__name__).info(
                            "[CLUSTER] Wrote backup candidates table → %s", _backup_tsv,
                        )
                except Exception as _ce:
                    raise RuntimeError(f'Cluster report generation failed: {_ce}') from _ce

                # Group nearby assessed leaves into compact, labelled local-clade figures.
                try:
                    if Path(_tree_nwk).exists():
                        _neighbourhood_result = neighbourhood.generate_local_neighbourhood_visuals(
                            tree_path=_tree_nwk,
                            assessment_rows=assessment_rows,
                            db=db,
                            outdir=Path(outdir) / 'neighbourhoods',
                            alignment_path=Path(outdir) / 'current_alignment.fasta',
                            image_format=getattr(args, 'neighbourhood_format', 'png'),
                        )
                        logging.getLogger(__name__).info(
                            "[NEIGHBOURHOOD] %d figure(s), %d/%d assessed sequences resolved",
                            len(_neighbourhood_result['figures']),
                            _neighbourhood_result['resolved'],
                            len(assessment_rows),
                        )
                        if _neighbourhood_result['resolved'] != len(assessment_rows):
                            logging.getLogger(__name__).warning(
                                "[NEIGHBOURHOOD] Local tree context resolved %d/%d assessed "
                                "sequences; unresolved rows remain in neighbourhood_manifest.tsv "
                                "with LocalNeighbourhoodFigure=NA",
                                _neighbourhood_result['resolved'],
                                len(assessment_rows),
                            )
                    else:
                        logging.getLogger(__name__).warning(
                            "[NEIGHBOURHOOD] No current_tree.nwk was produced; local-clade figures were skipped"
                        )
                except Exception as _ne:
                    raise RuntimeError(f'Local-clade figure generation failed: {_ne}') from _ne

                # Propose rolling primary + backup isolates within each GTDB
                # species or, where species is unresolved, each local tree clade.
                try:
                    sequencing_sets_path = _selection_sets.build_sequencing_sets(
                        assessment_rows,
                        Path(outdir) / 'sequencing_sets.tsv',
                        tree_path=_tree_nwk if Path(_tree_nwk).exists() else None,
                        db=db,
                        pangenome_target=getattr(
                            args, 'pangenome_target', _selection_sets.DEFAULT_PANGENOME_TARGET,
                        ),
                        candidate_set_size=getattr(
                            args, 'candidate_set_size', _selection_sets.DEFAULT_CANDIDATE_SET_SIZE,
                        ),
                        baseline_redundancy_identity=getattr(
                            args, 'baseline_redundancy_identity',
                            _selection_sets.DEFAULT_BASELINE_REDUNDANCY_IDENTITY,
                        ),
                        baseline_redundancy_min_query_coverage=getattr(
                            args, 'baseline_redundancy_min_query_coverage',
                            _selection_sets.DEFAULT_BASELINE_REDUNDANCY_QUERY_COVERAGE,
                        ),
                        baseline_extension_min_identity=getattr(
                            args, 'baseline_extension_min_identity',
                            _selection_sets.DEFAULT_BASELINE_EXTENSION_MIN_IDENTITY,
                        ),
                        baseline_extension_min_query_coverage=getattr(
                            args, 'baseline_extension_min_query_coverage',
                            _selection_sets.DEFAULT_BASELINE_EXTENSION_MIN_QUERY_COVERAGE,
                        ),
                    )
                    logging.getLogger(__name__).info(
                        "[SELECTION SETS] Wrote rolling clade candidate sets to %s",
                        sequencing_sets_path,
                    )
                    current_statuses = db.get_isolate_statuses(
                        [row.get('id') for row in assessment_rows if row.get('id')]
                    )
                    protected_states = {
                        'SAB_APPROVED', 'DNA_EXTRACTION_PENDING', 'DNA_EXTRACTION_FAILED',
                        'LIBRARY_PENDING', 'SEQUENCED', 'GENOME_QC_FAILED',
                        'GENOME_QC_PASSED', 'WITHDRAWN',
                    }
                    for row in assessment_rows:
                        if str(row.get('sequencing_set_role') or '').upper() not in {'PRIMARY', 'BACKUP'}:
                            continue
                        sequence_id = str(row.get('id') or '')
                        if current_statuses.get(sequence_id, {}).get('status') in protected_states:
                            continue
                        db.update_isolate_status(
                            sequence_id, 'PROPOSED',
                            detail=str(row.get('sequencing_set_reason') or ''),
                            source_file=str(sequencing_sets_path),
                        )
                except Exception as _se:
                    raise RuntimeError(f'Candidate-set generation failed: {_se}') from _se

                # Candidate-set roles/ranks are assigned after local contexts
                # are known. Re-render the same deterministic figure groups so
                # PRIMARY/BACKUP/ALTERNATE labels appear in the final images.
                try:
                    if Path(_tree_nwk).exists():
                        neighbourhood.generate_local_neighbourhood_visuals(
                            tree_path=_tree_nwk,
                            assessment_rows=assessment_rows,
                            db=db,
                            outdir=Path(outdir) / 'neighbourhoods',
                            alignment_path=Path(outdir) / 'current_alignment.fasta',
                            image_format=getattr(args, 'neighbourhood_format', 'png'),
                        )
                except Exception as _ne:
                    raise RuntimeError(f'Ranked local-clade figure refresh failed: {_ne}') from _ne

                # Final report write after tree-derived fields have been attached.
                assess_path = write_sequence_assessment_tsv(
                    Path(outdir) / 'sequence_assessment.tsv', assessment_rows,
                    assessment_db_name=run_main_db_name,
                )
                selection_summary_path = write_selection_summary_tsv(
                    Path(outdir) / 'selection_summary.tsv', assessment_rows,
                    assessment_db_name=run_main_db_name,
                )
                baseline_hits_path = write_baseline_hits_tsv(Path(outdir) / 'baseline_hits.tsv', assessment_rows)
                try:
                    snapshot_id = f'{run_dataset}:{Path(outdir).resolve()}'
                    snapshot_source = Path(outdir).resolve() / 'assessment' / 'sequence_assessment.tsv'
                    saved = db.save_assessment_snapshot(
                        snapshot_id,
                        assessment_rows,
                        dataset=run_dataset,
                        source_path=str(snapshot_source),
                    )
                    logging.getLogger(__name__).info(
                        '[ASSESSMENT] Saved %d rows to project snapshot %s', saved, snapshot_id,
                    )
                except Exception as e:
                    raise RuntimeError(f'Could not persist project assessment snapshot: {e}') from e
                logging.getLogger(__name__).info("[PERFORMANCE REVIEW] Wrote sequence assessment to %s", assess_path)
                logging.getLogger(__name__).info("[PERFORMANCE REVIEW] Wrote SAB selection summary to %s", selection_summary_path)
                logging.getLogger(__name__).info("[PERFORMANCE REVIEW] Wrote nearest baseline hit report to %s", baseline_hits_path)

                # Emit a concise summary using the same decisions as selection_summary.tsv.
                try:
                    decisions = [build_selection_decision(row)['decision'] for row in assessment_rows]
                    logging.getLogger(__name__).info(
                        "[ASSESSMENT SUMMARY] %d assessed: %d PRIMARY, %d BACKUP, "
                        "%d DIVERSITY/SECONDARY, %d BASELINE REDUNDANT, %d REVIEW, "
                        "%d TARGET MET, %d ALREADY SEQUENCED. "
                        "See assessment/selection_summary.tsv for the decision table.",
                        len(assessment_rows),
                        sum(d == 'PRIORITISE - SET PRIMARY' for d in decisions),
                        sum(d == 'RESERVE - SET BACKUP' for d in decisions),
                        sum(d in ('STRONG CANDIDATE', 'SECONDARY CANDIDATE', 'SECONDARY - STRAIN DIVERSITY') for d in decisions),
                        sum(d == 'EXCLUDE - BASELINE REDUNDANT' for d in decisions),
                        sum(d.startswith('REVIEW') for d in decisions),
                        sum(d == 'LOWER PRIORITY - TARGET MET' for d in decisions),
                        sum(d == 'ALREADY SEQUENCED' for d in decisions),
                    )
                except Exception:
                    pass
            except Exception as e:
                raise RuntimeError(f'Required sequence assessment failed: {e}') from e
    except Exception as e:
        raise SystemExit(f'[PERFORMANCE REVIEW] Required assessment/placement reporting failed: {e}')

    # At end of run, keep the externally useful deliverables and organise them
    # into high-level folders.
    try:
        _organise_run_outputs(outdir, primary_db_name=run_main_db_name)
    except Exception as e:
        raise SystemExit(f'[PERFORMANCE REVIEW] Failed to organise required outputs: {e}')

    # Write a manifest after files have been organised.
    try:
        _write_output_explanations(outdir)
    except Exception:
        pass


def cmd_performance_review(args):
    """Run Performance Review transactionally and publish only after validation."""
    import fcntl
    import sqlite3
    import uuid

    # Internal callers predating the CLI command namespace are retained for
    # focused unit tests; every user-facing CLI invocation takes the atomic path.
    if not hasattr(args, 'command') or not hasattr(args, 'skip_chimera_check'):
        return _cmd_performance_review_impl(args)
    original_db = Path(args.db).expanduser().resolve()
    if str(args.db) == ':memory:':
        raise SystemExit('[PERFORMANCE REVIEW] A persistent --db path is required for the rolling workflow.')
    outdir = Path(args.out).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    original_db.parent.mkdir(parents=True, exist_ok=True)
    staged_db = outdir / '.performance_review_project.sqlite'
    lock_path = original_db.with_name(original_db.name + '.lock')
    manifest = RunManifest(outdir, 'performance_review')
    for role, source in (
        ('marker_fasta', args.input), ('partner_metadata', getattr(args, 'partner_metadata', None)),
        ('primary_reference', getattr(args, 'ref', None)), ('reference_taxonomy', getattr(args, 'taxa', None)),
        ('baseline_fasta', getattr(args, 'baseline_fasta', None)), ('mwl', getattr(args, 'mwl', None)),
        ('marker_qc', getattr(args, 'marker_qc', None)), ('marker_review', getattr(args, 'marker_review', None)),
        ('baseline_taxonomy', getattr(args, 'baseline_taxa_assignments', None)),
        ('chimera_reference', getattr(args, 'chimera_ref', None)),
        ('target_collection', getattr(args, 'target', None)),
        ('anchors', getattr(args, 'anchors', None)),
    ):
        if source:
            manifest.add_input(source, role=role, required=True)
    for source in getattr(args, 'alt_ref', None) or []:
        manifest.add_input(source, role='alternate_reference')
    for source in getattr(args, 'alt_taxa', None) or []:
        manifest.add_input(source, role='alternate_reference_taxonomy')
    run_id = f"performance-review:{getattr(args, 'dataset', 'dataset')}:{utc_now()}:{uuid.uuid4().hex[:8]}"

    with open(lock_path, 'a+') as lock_handle:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        try:
            if staged_db.exists():
                staged_db.unlink()
            if original_db.is_file():
                with sqlite3.connect(str(original_db)) as source, sqlite3.connect(str(staged_db)) as target:
                    source.backup(target)
                manifest.add_input(original_db, role='project_database_before')
            staged_args = argparse.Namespace(**vars(args))
            staged_args.db = str(staged_db)
            manifest.add_stage('performance_review', 'RUNNING', detail='taxonomy, novelty, MSA/tree, and assessment')
            _cmd_performance_review_impl(staged_args)

            safe_main = ''.join(
                char if char.isalnum() or char in ('_', '-') else '_'
                for char in str(getattr(args, 'main_ref', None) or getattr(args, 'ref_name', None) or 'GTDB')
            ).strip('_') or 'GTDB'
            required = {
                'full_assessment': outdir / 'assessment' / 'sequence_assessment.tsv',
                'board_summary': outdir / 'assessment' / 'selection_summary.tsv',
                'sequencing_sets': outdir / 'assessment' / 'sequencing_sets.tsv',
                'novelty_metrics': outdir / 'assessment' / 'novelty_metrics.tsv',
                'tree': outdir / 'tree' / 'current_tree.nwk',
                'alignment': outdir / 'tree' / 'current_alignment.fasta',
                'taxonomy': outdir / 'taxonomy' / f'{safe_main}.tsv',
            }
            if not required['taxonomy'].is_file():
                taxonomy_candidates = sorted((outdir / 'taxonomy').glob('*.tsv'))
                if taxonomy_candidates:
                    required['taxonomy'] = taxonomy_candidates[0]
            for role, path in required.items():
                manifest.add_output(path, role=role)
            missing = manifest.verify_required_outputs()
            if missing:
                raise RuntimeError('required Performance Review outputs are missing: ' + ', '.join(missing))

            from branchmanager.reporting import write_decision_changes, write_performance_review_dashboard
            staged = Database(str(staged_db))
            changes = write_decision_changes(staged, outdir / 'assessment' / 'decision_changes.tsv')
            manifest.add_output(changes, role='decision_changes')
            dashboard = write_performance_review_dashboard(
                outdir / 'assessment' / 'selection_summary.tsv', outdir,
            )
            manifest.add_output(dashboard, role='performance_review_dashboard')
            manifest.add_stage('performance_review', 'COMPLETE')
            manifest.add_stage('hiring_panel', 'COMPLETE', detail='candidate sets and board summary')
            manifest.finish('COMPLETE')
            staged.record_project_run(
                run_id, 'performance_review', 'COMPLETE', dataset=getattr(args, 'dataset', ''),
                manifest_path=str(manifest.json_path),
                started_at=manifest.data['started_at'], completed_at=manifest.data['completed_at'],
            )

            publish_tmp = original_db.with_name(original_db.name + '.publishing')
            if publish_tmp.exists():
                publish_tmp.unlink()
            with sqlite3.connect(str(staged_db)) as source, sqlite3.connect(str(publish_tmp)) as target:
                source.backup(target)
                if target.execute('PRAGMA integrity_check').fetchone()[0] != 'ok':
                    raise RuntimeError('staged project database failed SQLite integrity_check')
            os.replace(publish_tmp, original_db)
            staged_db.unlink(missing_ok=True)
        except BaseException as exc:
            manifest.finish('FAILED', error=exc)
            raise
        finally:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)


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


def cmd_org_chart(args):
    """Build a focused phylogenetic tree for a specific taxon from an existing DB.

    Fast path: if a ``current_alignment.fasta`` already exists in
    ``--from-dir`` (or the output directory), the matching sequences are
    extracted from that *pre-built* alignment and FastTree is run on the
    subset — no re-alignment needed.

    Slow path (fallback): if no existing alignment is found, the matching
    sequences are exported from the DB and a full MAFFT + FastTree build
    is performed (same as a normal run).
    """
    from branchmanager.pipeline import itol as itol_mod
    from branchmanager.pipeline import tree as tree_mod  # noqa: F401 (used in _build_org_chart)
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
        "[ORG-CHART] Taxon query: '%s'  rank: '%s'  clean name: '%s'  from-dir: %s",
        taxon_query, rank_key, taxon_clean, from_dir,
    )

    db = Database(args.db)
    db.initialise()

    with db.connect() as conn:
        cur = conn.cursor()
        cur.execute(
            "SELECT s.id, s.sequence, t.taxonomy, t.confidence, s.dataset "
            "FROM sequences s LEFT JOIN taxonomy t ON s.id = t.id"
        )
        all_rows = cur.fetchall()

    matched = []
    for sid, seq, tax, conf, dataset in all_rows:
        if not tax:
            continue
        parsed = parse_taxon_string(tax)
        val = parsed.get(rank_key, '')
        if _taxon_name_matches(val, taxon_clean):
            matched.append((sid, seq or '', tax, conf or 'NA', dataset or ''))

    log.info("[ORG-CHART] DB has %d sequences total; %d match taxon '%s' at rank '%s'",
             len(all_rows), len(matched), taxon_clean, rank_key)

    if len(matched) < min_seqs:
        log.warning(
            "[ORG-CHART] Only %d sequences matched (minimum required: %d). "
            "Check taxon spelling and rank. Available phyla in DB: %s",
            len(matched), min_seqs,
            sorted({parse_taxon_string(r[2]).get('p', 'unknown')
                    for r in all_rows if r[2]}),
        )
        return

    # Collect matched IDs as a set for fast lookup
    matched_ids = {sid for sid, *_ in matched}

    combined_tax_path = Path(outdir) / 'org_chart_combined_taxonomy.tsv'
    write_combined_taxonomy_tsv(
        combined_tax_path,
        [(sid, tax, conf) for sid, _, tax, conf, _ in matched],
    )
    log.info("[ORG-CHART] Wrote taxonomy for %d sequences → %s", len(matched), combined_tax_path)

    summary_path = Path(outdir) / 'org_chart_sequence_list.tsv'
    with open(summary_path, 'w') as sf:
        sf.write('ID\tTaxonomy\tConfidence\tDataset\n')
        for sid, _, tax, conf, dataset in matched:
            sf.write(f"{sid}\t{tax}\t{conf}\t{dataset}\n")
    log.info("[ORG-CHART] Wrote sequence list → %s", summary_path)

    tree_path = Path(outdir) / 'org_chart_tree.nwk'

    if not no_tree:
        _build_org_chart(
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
        log.info("[ORG-CHART] Tree build skipped (--no-tree).")

    try:
        tfile = str(tree_path) if tree_path.exists() else None
        itol_mod.generate_itol_colours(
            str(combined_tax_path),
            outdir,
            tree_file=tfile,
            phylum_groups=getattr(args, 'group_phyla', None),
        )
        log.info("[ORG-CHART] Generated iTOL colour files in %s", outdir)
    except Exception as e:
        log.warning("[ORG-CHART] iTOL colour generation failed: %s", e)

    # Optional: write functional annotation datasets when provided
    try:
        func_tsv = getattr(args, 'functional', None)
        if func_tsv:
            try:
                written = itol_mod.write_functional_annotations(str(func_tsv), outdir, id_map=None)
                log.info("[ORG-CHART] Wrote functional annotation iTOL files: %s", ','.join(written) if written else '(none)')
            except Exception as e:
                log.warning("[ORG-CHART] Functional annotations generation failed: %s", e)
    except Exception:
        pass
    # Draft rumen functional groups
    if getattr(args, 'draft_rumen_functions', False):
        try:
            tsv_out, itol_out = itol_mod.generate_rumen_function_draft(
                str(combined_tax_path), outdir, id_map=None
            )
            if tsv_out:
                log.info("[ORG-CHART] Draft rumen functional annotation: %s", tsv_out)
            if itol_out:
                log.info("[ORG-CHART] Rumen functional iTOL file: %s", itol_out)
        except Exception as e:
            log.warning("[ORG-CHART] Draft rumen functions generation failed: %s", e)

    # One colour per dataset label so users can see which sequences came from
    # which dataset (Filing Cabinet, review batch, etc.) in the same tree view.
    try:
        ids_in_order = [sid for sid, *_ in matched]
        ds_map = {sid: (dataset or 'unknown') for sid, _, _tax, _conf, dataset in matched}
        membership_path = Path(outdir) / 'itol_dataset_membership.itol'
        itol_mod.write_dataset_membership_strip(
            str(membership_path), ids_in_order, ds_map,
            dataset_label='Dataset membership',
        )
        log.info("[ORG-CHART] Wrote dataset membership strip → %s", membership_path)
    except Exception as e:
        log.warning("[ORG-CHART] Dataset membership strip failed: %s", e)

    log.info(
        "[ORG-CHART] Complete. %d sequences, taxon='%s', output=%s",
        len(matched), taxon_query, outdir,
    )
    print(
        f"[org-chart] Done: {len(matched)} sequences for '{taxon_query}' "
        f"→ {outdir}"
    )


def _build_org_chart(
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
        is_ref_anchor,
    )
    from branchmanager.utils.fasta import read_fasta, write_fasta
    from branchmanager.pipeline import tree as tree_mod
    import re

    out = Path(outdir)
    # Reuse an existing alignment when available.
    aln_candidates = [
        Path(from_dir) / 'current_alignment.fasta',
        Path(from_dir) / 'tree' / 'current_alignment.fasta',
        Path(outdir) / 'current_alignment.fasta',
        Path(outdir) / 'tree' / 'current_alignment.fasta',
    ]
    existing_aln = next((p for p in aln_candidates if p.exists()), None)

    if existing_aln:
        log.info("[ORG-CHART] Fast path: filtering existing alignment %s", existing_aln)
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
            "[ORG-CHART] Fast path: %d data sequences + %d anchors extracted from alignment",
            n_data, anchors_included,
        )

        if n_data < 3:
            log.warning(
                "[ORG-CHART] Only %d data sequences found in existing alignment "
                "(IDs may not match). Falling back to slow path.", n_data,
            )
            existing_aln = None   # trigger slow path below
        else:
            org_chart_alignment = out / 'org_chart_alignment.fasta'
            write_fasta(kept, str(org_chart_alignment))
            fasta_for_tree, id_map = _make_unique_fasta(str(org_chart_alignment), out)
            if _run_fasttree(Path(fasta_for_tree), tree_path):
                _finalise_org_chart_tree(out, id_map, tree_path, log)
            else:
                log.warning("[ORG-CHART] FastTree failed on fast path")
            return

    log.info("[ORG-CHART] Slow path: full MAFFT + FastTree build")
    seqs_fasta = out / 'org_chart_input_sequences.fasta'
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
            log.info("[ORG-CHART] Slow path tree written → %s", tree_path)
    except Exception as e:
        log.warning("[ORG-CHART] Slow path tree build failed: %s", e)


def _finalise_org_chart_tree(out: Path, id_map: dict, tree_path: Path, log) -> None:
    """Remap IDs, prune anchors, and label internal nodes for Org Chart output."""
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
        log.info("[ORG-CHART] Pruned %d anchor leaves from focused-tree newick", before - after)
    newick = _repair_internal_node_label_delimiters(newick)
    newick = _label_internal_nodes(newick)
    tree_path.write_text(newick)
    log.info("[ORG-CHART] Focused tree finalised → %s", tree_path)


def cmd_label_maker(args):
    db = Database(args.db)
    db.initialise()
    outdir = args.out
    os.makedirs(outdir, exist_ok=True)
    _configure_logging(outdir)
    log = logging.getLogger(__name__)

    try:
        with db.connect() as conn:
            cur = conn.cursor()
            # If the outdir contains a Filing Cabinet combined taxonomy file, prefer
            # to use only those IDs so the regenerated iTOL matches the previous review
            # statistics. This mirrors the behaviour used during `run` when a
            # --previous-review is supplied.
            try:
                p = Path(outdir)
                # Prefer an existing combined taxonomy file in the outdir, but
                # only if it contains data (header + >=1 data row). Fall back
                # to the DB-wide query otherwise to avoid regenerating empty
                # iTOL outputs when a stub file exists.
                previous_review_file = None
                for cand in _combined_taxonomy_candidates(p):
                    if not cand.exists():
                        continue
                    try:
                        # count lines cheaply
                        with open(cand) as pf:
                            cnt = sum(1 for _ in pf)
                        if cnt > 1:
                            previous_review_file = cand
                            break
                        else:
                            # file exists but only header or empty -> ignore
                            continue
                    except Exception:
                        continue
                previous_review_ids = None
                if previous_review_file is not None:
                    previous_review_ids = read_combined_taxonomy_ids(previous_review_file)
                else:
                    previous_review_ids = None
            except Exception:
                previous_review_ids = None

            previous_review_id_list = previous_review_ids if isinstance(previous_review_ids, list) else []

            # fetch id, taxonomy, confidence and dataset for all sequences (or filter)
            if getattr(args, 'include_datasets', None):
                ds_list = [d.strip() for d in args.include_datasets.split(',') if d.strip()]
                placeholders = ','.join('?' for _ in ds_list)
                cur.execute(f"SELECT s.id, t.taxonomy, t.confidence, s.dataset FROM sequences s LEFT JOIN taxonomy t ON s.id = t.id WHERE s.dataset IN ({placeholders})", tuple(ds_list))
            elif previous_review_id_list:
                placeholders = ','.join('?' for _ in previous_review_id_list)
                cur.execute(f"SELECT s.id, t.taxonomy, t.confidence, s.dataset FROM sequences s LEFT JOIN taxonomy t ON s.id = t.id WHERE s.id IN ({placeholders})", tuple(previous_review_id_list))
            else:
                cur.execute("SELECT s.id, t.taxonomy, t.confidence, s.dataset FROM sequences s LEFT JOIN taxonomy t ON s.id = t.id")
            rows = cur.fetchall()
    except Exception as e:
        log.warning("[LABEL-MAKER] Failed to query DB: %s", e)
        return

    try:
        kingdom = _sequence_domain_filter(args, default=None)
        if kingdom:
            kingdom_text = str(kingdom)
            rows = [row for row in rows if row[1] and taxonomy_matches_kingdom(str(row[1]), kingdom_text)]
            log.info("[LABEL-MAKER] Retained %d rows after kingdom filter=%s", len(rows), kingdom)
    except Exception as e:
        log.warning("[LABEL-MAKER] Kingdom filter failed: %s", e)

    combined_path = Path(outdir) / 'combined_taxonomy.tsv'
    try:
        write_combined_taxonomy_tsv(combined_path, [(rid, tax, conf) for rid, tax, conf, ds in rows])
        log.info("[LABEL-MAKER] Wrote combined taxonomy for %d ids to %s", len(rows), combined_path)
    except Exception as e:
        log.warning("[LABEL-MAKER] Failed to write combined taxonomy: %s", e)
        return

    # call itol generator
    try:
        tree_path = Path(outdir) / 'current_tree.nwk'
        tfile = str(tree_path) if tree_path.exists() else _find_tree_file_in_dir(outdir)
        itol.generate_itol_colours(str(combined_path), outdir, tree_file=tfile, phylum_groups=getattr(args, 'group_phyla', None))
        log.info("[LABEL-MAKER] Generated iTOL colour files in %s", outdir)
    except Exception as e:
        log.warning("[LABEL-MAKER] Colour generation failed: %s", e)

    # Optional: write functional annotation datasets when provided
    try:
        func_tsv = getattr(args, 'functional', None)
        if func_tsv:
            try:
                written = itol.write_functional_annotations(str(func_tsv), outdir, id_map=None)
                log.info("[LABEL-MAKER] Wrote functional annotation iTOL files: %s", ','.join(written) if written else '(none)')
            except Exception as e:
                log.warning("[LABEL-MAKER] Functional annotations generation failed: %s", e)
    except Exception:
        pass
    # Draft rumen functional groups
    if getattr(args, 'draft_rumen_functions', False):
        try:
            tsv_out, itol_out = itol.generate_rumen_function_draft(str(combined_path), outdir, id_map=None)
            if tsv_out:
                log.info("[LABEL-MAKER] Draft rumen functional annotation: %s", tsv_out)
            if itol_out:
                log.info("[LABEL-MAKER] Rumen functional iTOL file: %s", itol_out)
        except Exception as e:
            log.warning("[LABEL-MAKER] Draft rumen functions generation failed: %s", e)

    # build dataset membership strip
    try:
        ids_in_order = [r[0] for r in rows]
        ds_map = {r[0]: (r[3] or '') for r in rows}
        membership_path = Path(outdir) / 'itol_dataset_membership.itol'
        itol.write_dataset_membership_strip(str(membership_path), ids_in_order, ds_map)
        log.info("[LABEL-MAKER] Wrote dataset membership iTOL file to %s", membership_path)
    except Exception as e:
        log.warning("[LABEL-MAKER] Failed to build/write dataset membership iTOL file: %s", e)

    # write explanations for regenerated outputs
    try:
        _write_output_explanations(outdir)
    except Exception:
        pass


def build_parser():
    _Fmt = type('_Fmt', (
        argparse.ArgumentDefaultsHelpFormatter,
        argparse.RawDescriptionHelpFormatter,
    ), {'max_help_position': 40, 'width': 120})
    paper_trail_policy = _paper_trail_module.DEFAULT_QC_POLICY
    selection_target = _selection_sets.DEFAULT_PANGENOME_TARGET
    selection_panel_size = _selection_sets.DEFAULT_CANDIDATE_SET_SIZE
    baseline_redundancy_identity = _selection_sets.DEFAULT_BASELINE_REDUNDANCY_IDENTITY
    baseline_redundancy_coverage = _selection_sets.DEFAULT_BASELINE_REDUNDANCY_QUERY_COVERAGE
    baseline_extension_identity = _selection_sets.DEFAULT_BASELINE_EXTENSION_MIN_IDENTITY
    baseline_extension_coverage = _selection_sets.DEFAULT_BASELINE_EXTENSION_MIN_QUERY_COVERAGE
    parser = argparse.ArgumentParser(
        prog='branchmanager',
        description=(
            'BranchManager — marker-gene QC, taxonomy, novelty scoring, and isolate prioritisation toolkit.\n\n'
            'Subcommands:\n'
            '  mailroom           Inventory AB1 deliveries and build a batch map.\n'
            '  interview          Run standalone AB1 QC from a Mailroom batch.\n'
            '  background-check   Pre-classify reference collections once and reuse the evidence.\n'
            '  filing-cabinet     Register the cultured baseline and backbone tree.\n'
            '  onboarding         Validate a partner submission before project state changes.\n'
            '  paper-trail        Interpret chromatograms and assemble primer reads.\n'
            '  performance-review Assess taxonomy, novelty, phylogeny, and sequencing candidates.\n'
            '  status-meeting     Import factual isolate lifecycle changes.\n'
            '  records-update     Import completed-genome, QC, GTDB, and ANI evidence.\n'
            '  quarterly-review   Nominate a later genome tranche from the complete project.\n'
            '  exit-interview     Remove isolates from active project state with a retained audit.\n'
            '  annual-report      Produce the cumulative end-of-project report.\n'
            '  it-desk            Check dependencies and project inputs before production.\n'
            '  assistant          Run the complete raw-trace-to-Hiring-Panel workflow.\n'
            '  org-chart          Extract a focused tree and iTOL files for a specific taxon.\n'
            '  label-maker        Regenerate iTOL annotation files from stored taxonomy.\n\n'
            'Typical workflow:\n'
            '  P. branchmanager mailroom --read-dir All_AB1 --metadata supplier.csv --dataset QUB_01 --forward-primer 63F -o QUB_01\n'
            '  Q. branchmanager interview --mailroom QUB_01 -o QUB_01_interview\n'
            '  0. branchmanager background-check --dataset hungate16s=hungate.fasta --ref gtdb.fna --taxa gtdb_tax.tsv -o background_check_out\n'
            '  1. branchmanager filing-cabinet --fasta baseline.fasta --db project.db --dataset Hungate --taxa-assignments background_check_out/pipeline_taxonomy.tsv --build-tree -o filing_cabinet_out\n'
            '  2. branchmanager performance-review --input new_seqs.fasta --partner-metadata new_seqs_metadata.tsv --db project.db --dataset Batch1 --ref gtdb.fna --mwl MWL.xlsx -o review_out\n'
            '  3. branchmanager quarterly-review --db project.db --genome-budget 24 --tree review_out/tree/current_tree.nwk -o quarterly_review_01\n'
            '  4. branchmanager org-chart --db project.db --taxon archaea --from-dir filing_cabinet_out -o archaea_out\n'
        ),
        formatter_class=_Fmt,
    )
    sub = parser.add_subparsers(dest='command')

    filing_cabinet_parser = sub.add_parser(
        'filing-cabinet',
        help='Filing Cabinet: register a cultured baseline and build the reference tree.',
        description=(
            'Filing Cabinet registers a baseline FASTA dataset (e.g. Hungate 16S) in the DB, optionally classifies\n'
            'sequences against a reference (GTDB/SILVA), collapse near-identical sequences,\n'
            'and build the backbone phylogenetic tree.\n\n'
            'This is always the first step. All subsequent Performance Reviews measure novelty against\n'
            'sequences stored by this command.'
        ),
        formatter_class=_Fmt,
    )
    filing_cabinet_parser.add_argument('--fasta', required=True,
        help='Input FASTA file containing the baseline sequences to load.')
    filing_cabinet_parser.add_argument('--db', required=True,
        help='Path to the BranchManager SQLite database (created if it does not exist).')
    filing_cabinet_parser.add_argument('-o', '--out', required=False, default='.',
        help='Output directory for tree, iTOL files, and reports (default: current directory).')
    filing_cabinet_parser.add_argument('--dataset', required=True,
        help='Label for this dataset stored in the DB (e.g. Hungate). Used to colour iTOL strips.')
    filing_cabinet_parser.add_argument('--baseline-tier', choices=['priority', 'secondary'], default=None,
        help='Cultured baseline tier for this dataset: priority for Hungate, secondary for other rumen isolate genomes. Inferred from --dataset when omitted.')
    filing_cabinet_parser.add_argument('--shorten-ids', dest='shorten_ids',
        action=argparse.BooleanOptionalAction, default=False,
        help='Replace input headers with compact IDs (e.g. HUN001). Default is to preserve the IDs exactly as supplied.')
    filing_cabinet_parser.add_argument('--classify', action='store_true',
        help='Classify sequences against --ref and store taxonomy in the DB. Requires --ref.')
    filing_cabinet_parser.add_argument('--build-tree', action='store_true',
        help='Build the backbone MAFFT + FastTree phylogenetic tree after loading.')
    filing_cabinet_parser.add_argument('--ref', required=False,
        help='Reference FASTA (GTDB/SILVA reps) for classification and tree orientation. Preferred over --taxa-assignments for externally classified inputs.')
    filing_cabinet_parser.add_argument('--taxa', required=False,
        help='Taxonomy table matching IDs in --ref (TSV/CSV, optionally .gz: id<TAB>lineage or id,lineage). Optional when --ref FASTA headers already contain lineages.')
    filing_cabinet_parser.add_argument('--ref-name', dest='ref_name', required=False, default=None,
        help='Display name for the primary reference database (default: derived from --ref filename). Used to label taxonomy columns.')
    filing_cabinet_parser.add_argument('--alt-ref', dest='alt_ref', action='append', default=None, metavar='FASTA',
        help='Additional reference FASTA to classify against (repeatable). Produces extra taxonomy columns in output files.')
    filing_cabinet_parser.add_argument('--alt-taxa', dest='alt_taxa', action='append', default=None, metavar='TABLE',
        help='Taxonomy TSV/CSV for the corresponding --alt-ref (positionally paired; repeatable; .gz accepted).')
    filing_cabinet_parser.add_argument('--alt-ref-name', dest='alt_ref_name', action='append', default=None, metavar='NAME',
        help='Display name for the corresponding --alt-ref (positionally paired; repeatable). Default: derived from filename.')
    filing_cabinet_parser.add_argument('--main-ref', dest='main_ref', required=False, default=None,
        help='Name of the reference database to use as the primary taxonomy source (default: primary --ref). Must match one of the --ref-name / --alt-ref-name values.')
    filing_cabinet_parser.add_argument('--taxa-assignments',
        dest='taxa_assignments', required=False,
        help='Pre-computed taxonomy assignments for the INPUT sequences (TSV/CSV, optionally .gz: query_id + lineage, or a FASTA with embedded lineages). Use this instead of --classify when you already have taxonomy.')
    filing_cabinet_parser.add_argument('--collapse', action='store_true',
        help='Collapse sequences that share ≥ --collapse-threshold identity AND the same taxonomy into a single representative for the tree. Saves time and reduces visual clutter.')
    filing_cabinet_parser.add_argument('--collapse-threshold', type=float, default=99.8,
        help='Identity threshold (percent) for collapsing duplicate-like sequences (default: 99.8).')
    filing_cabinet_parser.add_argument('--sequence-domain', '--organism-domain', dest='sequence_domain',
        choices=SEQUENCE_DOMAIN_CHOICES, default=None,
        help=(
            'Sequence/domain profile for this Filing Cabinet. Default behaviour is bacteria. '
            'Use archaea for archaeal 16S, fungi for fungal/eukaryotic runs with suitable refs/anchors, '
            'or mixed/all/none to disable domain filtering.'
        ))
    filing_cabinet_parser.add_argument('--anchors', required=False, default=None,
        help='Custom reference anchor FASTA for tree topology scaffolding. Defaults to the 26-sequence bundled anchor set (src/branchmanager/data/reference_anchors.fasta).')
    filing_cabinet_parser.add_argument('--threads', type=int, required=False, default=4,
        help='Number of CPU threads for MAFFT and VSEARCH (default: 4).')
    filing_cabinet_parser.add_argument('--tree-method', dest='tree_method',
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
    filing_cabinet_parser.add_argument('--colours', required=False,
        help='CSV file mapping sequence IDs to custom hex colours for iTOL (columns: id, colour).')
    filing_cabinet_parser.add_argument(
        '--group-phyla', dest='group_phyla', action='append', default=None, metavar='SPEC',
        help=(
            'Collapse multiple phyla into a single colour in iTOL legends. Repeatable. '
            'Formats: "archaea" (all archaeal phyla), "bacteria" (all bacterial phyla), '
            '"Bacillota,Bacillota_I" (explicit list; label = first name), '
            '"Firmicutes:Bacillota,Bacillota_I" (named group).'
        ),
    )
    filing_cabinet_parser.add_argument('--functional', dest='functional', required=False, default=None,
        help=(
            'TSV file mapping sequence IDs to functional attributes (pathways, functions, '
            'traits, scores, etc.). Header row required; first column = sequence ID; '
            'subsequent columns = one functional attribute each. '
            'One iTOL file is generated per column: binary (0/1/yes/no) → DATASET_BINARY, '
            'numeric → DATASET_SIMPLEBAR, categorical → DATASET_COLORSTRIP.'
        ),
    )
    filing_cabinet_parser.add_argument('--draft-rumen-functions', dest='draft_rumen_functions',
        action='store_true', default=False,
        help=(
            'Auto-generate a draft rumen functional-group annotation from the output taxonomy. '
            'Maps each sequence to a broad ruminant microbiome functional category '
            '(e.g. Cellulolytic/Fibrolytic, Methanogenic Archaea, Butyrate Producers) '
            'and writes rumen_functions_draft.tsv + itol_func_Rumen_Functional_Group.itol. '
            'The draft TSV can be edited and re-supplied via --functional in future runs.'
        ),
    )

    performance_review_parser = sub.add_parser(
        'performance-review',
        help='Performance Review: assess taxonomy, novelty, phylogeny, and sequencing candidates.',
        description=(
            'Assess new partner 16S isolate sequences against the project baseline.\n\n'
            'The workflow classifies against GTDB (primary), optionally cross-checks NCBI/GG2/SILVA\n'
            'as --alt-ref databases, scores novelty against separate cultured-baseline and rolling partner\n'
            'collections, updates the tree, and optionally matches GTDB taxonomy against\n'
            'the Most Wanted List via --mwl.\n\n'
            'Provide --baseline-fasta for context datasets such as Hungate when they have not already\n'
            'been registered with `branchmanager filing-cabinet`. Novelty is always relative to YOUR submitted data,\n'
            'not the full external reference. Each successive Performance Review extends the partner collection and\n'
            'updates same-species genome coverage and candidate sets.'
        ),
        formatter_class=_Fmt,
    )
    performance_review_parser.add_argument('--input', required=True,
        help='FASTA file of new sequences to analyse.')
    performance_review_parser.add_argument('--db', required=True,
        help='Path to the BranchManager SQLite database (must have been initialised with `branchmanager filing-cabinet`).')
    performance_review_parser.add_argument('-o', '--out', required=True,
        help='Output directory for this run (sequence_assessment.tsv, novelty_metrics.tsv, tree, iTOL files, etc.).')
    performance_review_parser.add_argument('--dataset', required=True,
        help='Label for this batch of sequences stored in the DB (e.g. Batch1). Used in iTOL dataset-membership strip.')
    performance_review_parser.add_argument('--ref', required=False,
        help='Reference FASTA (GTDB/SILVA reps) used for classification and tree orientation. Use the same file as the Filing Cabinet.')
    performance_review_parser.add_argument('--taxa', required=False,
        help='Taxonomy table matching IDs in --ref (TSV/CSV, optionally .gz: id<TAB>lineage or id,lineage).')
    performance_review_parser.add_argument('--ref-name', dest='ref_name', required=False, default=None,
        help='Display name for the primary reference database (default: derived from --ref filename).')
    performance_review_parser.add_argument('--alt-ref', dest='alt_ref', action='append', default=None, metavar='FASTA',
        help='Additional reference FASTA to classify against (repeatable). Adds extra taxonomy columns to sequence_assessment.tsv and taxonomy_all_dbs.tsv.')
    performance_review_parser.add_argument('--alt-taxa', dest='alt_taxa', action='append', default=None, metavar='TABLE',
        help='Taxonomy TSV/CSV for the corresponding --alt-ref (positionally paired; repeatable; .gz accepted).')
    performance_review_parser.add_argument('--alt-ref-name', dest='alt_ref_name', action='append', default=None, metavar='NAME',
        help='Display name for the corresponding --alt-ref (positionally paired; repeatable).')
    performance_review_parser.add_argument('--main-ref', dest='main_ref', required=False, default=None,
        help='Name of the authoritative assessment database. Bacterial/archaeal Performance Review requires GTDB; other databases are reporting cross-checks.')
    performance_review_parser.add_argument('--mwl', dest='mwl', required=False, default=None,
        help='Most Wanted List workbook/TSV/CSV. GTDB taxonomy is matched against this list and MWL columns are added to sequence_assessment.tsv.')
    performance_review_parser.add_argument('--mwl-sheet', dest='mwl_sheet', required=False, default='MWL_V1',
        help='Sheet name to read from an MWL .xlsx workbook (default: MWL_V1).')
    performance_review_parser.add_argument('--mwl-min-rank', dest='mwl_min_rank', required=False, default='p',
        choices=['domain', 'd', 'phylum', 'p', 'class', 'c', 'order', 'o', 'family', 'f', 'genus', 'g', 'species', 's'],
        help='Minimum matched rank required for an MWL hit (default: phylum). Domain-only MWL entries still match at domain.')
    performance_review_parser.add_argument('--partner-metadata', '--sequencing-metadata', dest='partner_metadata', required=False, default=None,
        help='Cumulative CSV/TSV ledger with sequence IDs, partner acronyms, optional selected_for_genome_sequencing, and required already_sequenced status. Selection is a commitment; already_sequenced means a genome is available. .gz is accepted.')
    performance_review_parser.add_argument('--marker-qc', dest='marker_qc', default=None,
        help='Paper Trail/Merge Meeting assembly_report.tsv sidecar. Auto-discovered beside --input when omitted.')
    performance_review_parser.add_argument('--marker-review', dest='marker_review', default=None,
        help='CSV/TSV manual-review decisions for PASS_WITH_WARNINGS or unverified marker sequences.')
    performance_review_parser.add_argument('--accept-unverified-marker-qc', action='store_true', default=False,
        help='Explicitly accept FASTA inputs lacking Paper Trail provenance. Recorded in the manifest and assessment; use only for independently validated sequences.')
    performance_review_parser.add_argument('--baseline-fasta', dest='baseline_fasta', required=False, default=None,
        help='Optional baseline/context FASTA to load before evaluating the new sequences (e.g. Hungate 16S).')
    performance_review_parser.add_argument('--baseline-dataset', dest='baseline_dataset', required=False, default='Baseline',
        help='Dataset label for --baseline-fasta in the DB and default cultured-baseline novelty pool (default: Baseline). Must differ from --dataset.')
    performance_review_parser.add_argument('--baseline-tier', choices=['priority', 'secondary'], default=None,
        help='Cultured baseline tier for --baseline-fasta: priority for Hungate, secondary for other rumen isolate genomes. Inferred from --baseline-dataset when omitted.')
    performance_review_parser.add_argument('--novelty-baseline-dataset', dest='novelty_baseline_datasets', action='append', default=[],
        help='Existing DB dataset label to include in the baseline/cultured novelty pool (repeatable; useful for Hungate plus other cultured isolate sets).')
    performance_review_parser.add_argument('--baseline-taxa-assignments', dest='baseline_taxa_assignments', required=False, default=None,
        help='Pre-computed taxonomy for --baseline-fasta (TSV/CSV, optionally .gz: sequence_id + lineage + optional confidence, or embedded-lineage FASTA). Skips baseline classification.')
    performance_review_parser.add_argument('--baseline-skip-classify', dest='baseline_skip_classify', action='store_true', default=False,
        help='Load --baseline-fasta into the DB without classifying it. Novelty still uses the baseline sequences, but taxonomy/iTOL context may be sparse.')
    performance_review_parser.add_argument('--baseline-colours', dest='baseline_colours', required=False, default=None,
        help='Optional CSV with baseline sequence colours, in the same format as filing-cabinet --colours.')
    performance_review_parser.add_argument('--baseline-shorten-ids', dest='baseline_shorten_ids',
        action=argparse.BooleanOptionalAction, default=False,
        help='Replace baseline FASTA headers with compact IDs. Default is to preserve the IDs exactly as supplied.')
    performance_review_parser.add_argument('--taxa-assignments',
        dest='taxa_assignments', required=False,
        help='Pre-computed taxonomy for the INPUT sequences (TSV/CSV, optionally .gz: query_id + lineage, or embedded-lineage FASTA).')
    performance_review_parser.add_argument('--previous-review', dest='previous_review', required=False,
        help='Path to the previous Filing Cabinet or Performance Review output directory. Used to seed the tree alignment so only new sequences need aligning.')
    performance_review_parser.add_argument('--shorten-ids', dest='shorten_ids',
        action=argparse.BooleanOptionalAction, default=False,
        help='Replace input headers with compact IDs. Default is to preserve the IDs exactly as supplied.')
    performance_review_parser.add_argument('--min-len', dest='min_len', type=int, default=800,
        help='Minimum sequence length to retain (bp, default: 800). Shorter sequences are filtered out.')
    performance_review_parser.add_argument('--max-n-percent', dest='max_n_percent', type=float, default=5.0,
        help='Maximum percentage of ambiguous (N) bases allowed (default: 5.0).')
    performance_review_parser.add_argument('--chimera-ref', dest='chimera_ref', default=None,
        help='Curated chimera-free marker reference for UCHIME. Defaults to the primary classification reference.')
    performance_review_parser.add_argument('--skip-chimera-check', dest='skip_chimera_check', action='store_true', default=False,
        help='Explicitly skip reference-based chimera screening. The omission is recorded and marker evidence is downgraded.')
    performance_review_parser.add_argument('--collapse', action='store_true',
        help='Collapse near-identical same-taxonomy sequences into representatives for the tree.')
    performance_review_parser.add_argument('--collapse-threshold', type=float, default=99.8,
        help='Identity threshold (percent) for collapsing (default: 99.8).')
    performance_review_parser.add_argument('--sequence-domain', '--organism-domain', dest='sequence_domain',
        choices=SEQUENCE_DOMAIN_CHOICES, default=None,
        help=(
            'Sequence/domain profile for this Performance Review. Omitted means bacteria. '
            'Use archaea for archaeal runs, fungi for fungal/eukaryotic runs with suitable references, '
            'or mixed/all/none to disable domain filtering. Provide domain-specific --ref/--alt-ref, '
            '--baseline-fasta, --previous-review, and --anchors as needed.'
        ))
    performance_review_parser.add_argument('--phylum', required=False,
        help='Filter iTOL output to sequences assigned to this phylum (e.g. Bacillota). Does not affect novelty scoring.')
    performance_review_parser.add_argument('--target', required=False, default=None,
        help=(
            'FASTA of sequences to measure novelty against instead of the DB. '
            'Leave unset to use the full rolling partner-candidate collection, including this batch with self-hits removed (recommended).'
        ),
    )
    performance_review_parser.add_argument('--force-rebuild', '--rebuild-tree', dest='force_rebuild', action='store_true', default=False,
        help='Rebuild the entire tree from scratch even when an existing alignment is present. '
             'When combined with --previous-review, ignores the previous alignment and jointly estimates '
             'tree topology across all datasets. (--rebuild-tree is an alias for this flag.)')
    performance_review_parser.add_argument('--anchors', required=False, default=None,
        help='Custom reference anchor FASTA for tree scaffolding. Defaults to bundled anchors.')
    performance_review_parser.add_argument('--threads', dest='threads', type=int, default=4,
        help='CPU threads for MAFFT and VSEARCH (default: 4).')
    performance_review_parser.add_argument('--tree-method', dest='tree_method',
        choices=['fasttree', 'iqtree', 'iqtree-fast'], default='fasttree',
        help=(
            'Phylogenetic tree-building backend (default: fasttree). '
            'fasttree: approximate ML, GTR+CAT. '
            'iqtree: full ML, GTR+G+I (recommended for production/publication runs). '
            'iqtree-fast: IQ-TREE 2 with -fast flag (good for exploratory incremental runs).'
        ),
    )
    performance_review_parser.add_argument(
        '--neighbourhood-format', dest='neighbourhood_format',
        choices=['png'], default='png',
        help='Image format for local phylogenetic-neighbourhood figures (PNG only).',
    )
    performance_review_parser.add_argument(
        '--pangenome-target', dest='pangenome_target', type=int, default=selection_target,
        help='Target number of committed genomes per GTDB species/local pangenome group (default: 9). Baseline isolates count because their genomes are available.',
    )
    performance_review_parser.add_argument(
        '--candidate-set-size', dest='candidate_set_size', type=int, default=selection_panel_size,
        help='Maximum ranked diversity-panel candidates per GTDB species/local clade (default: 9).',
    )
    performance_review_parser.add_argument(
        '--baseline-redundancy-identity', type=float, default=baseline_redundancy_identity,
        help='Exclude uncommitted candidates at or above this nearest cultured-baseline marker identity when coverage also passes its threshold (default: 99.8).',
    )
    performance_review_parser.add_argument(
        '--baseline-redundancy-min-query-coverage', type=float, default=baseline_redundancy_coverage,
        help='Minimum query coverage required before baseline near-identity exclusion is applied (default: 95).',
    )
    performance_review_parser.add_argument(
        '--baseline-extension-min-identity', type=float, default=baseline_extension_identity,
        help='Minimum cultured-baseline 16S identity for a same-species candidate to enter a baseline-pangenome extension (default: 98.65).',
    )
    performance_review_parser.add_argument(
        '--baseline-extension-min-query-coverage', type=float, default=baseline_extension_coverage,
        help='Minimum query coverage for admission to a baseline-pangenome extension (default: 95).',
    )
    performance_review_parser.add_argument('--user-colours', dest='user_colours', required=False,
        help='CSV file mapping sequence IDs to custom hex colours for iTOL (columns: id, colour).')
    performance_review_parser.add_argument(
        '--group-phyla', dest='group_phyla', action='append', default=None, metavar='SPEC',
        help=(
            'Collapse multiple phyla into one colour in iTOL legends. Repeatable. '
            'Formats: "archaea", "bacteria", "Bacillota,Bacillota_I", "Firmicutes:Bacillota,Bacillota_I".'
        ),
    )
    performance_review_parser.add_argument('--functional', dest='functional', required=False, default=None,
        help=(
            'TSV file mapping sequence IDs to functional attributes. '
            'Header row required; first column = sequence ID; subsequent columns = functional attributes. '
            'Generates one iTOL file per column (DATASET_BINARY / DATASET_SIMPLEBAR / DATASET_COLORSTRIP).'
        ),
    )
    performance_review_parser.add_argument('--draft-rumen-functions', dest='draft_rumen_functions',
        action='store_true', default=False,
        help=(
            'Auto-generate a draft rumen functional-group annotation from the output taxonomy. '
            'Writes rumen_functions_draft.tsv and itol_func_Rumen_Functional_Group.itol. '
            'The draft TSV can be edited and re-supplied via --functional in later runs.'
        ),
    )

    quarterly_review_parser = sub.add_parser(
        'quarterly-review',
        help='Select the next project-wide genome tranche after rolling Performance Reviews.',
        description=(
            'Reconsider every accumulated partner isolate in one auditable selection round.\n\n'
            'Residual nine-genome coverage gaps are filled first. Remaining budget is allocated\n'
            'across species using marginal tree/marker diversity, existing genome representation,\n'
            'cultured-baseline novelty, GTDB-reference context, and MWL evidence. Recommendations\n'
            'never change already_sequenced status. Use genome ANI/phylogenomics to validate\n'
            'strain-level diversity once genome data are available.'
        ),
        formatter_class=_Fmt,
    )
    quarterly_review_parser.add_argument('--db', required=True,
        help='Rolling BranchManager SQLite database containing all partner and baseline sequences.')
    quarterly_review_parser.add_argument('-o', '--out', required=True,
        help='Output directory for quarterly_review_summary.tsv, next_genome_set.tsv, and quarterly_review_manifest.tsv.')
    quarterly_review_parser.add_argument('--genome-budget', type=int, required=True,
        help='Number of new PRIMARY genomes to nominate in this round. This is intentionally explicit.')
    quarterly_review_parser.add_argument('--backups-per-primary', type=int, default=1,
        help='Number of extraction-failure/diversity backups to nominate per primary (default: 1).')
    quarterly_review_parser.add_argument('--pangenome-target', type=int, default=selection_target,
        help='Initial exact-GTDB-species/local-group coverage target before expansion (default: 9).')
    quarterly_review_parser.add_argument('--baseline-redundancy-identity', type=float, default=baseline_redundancy_identity,
        help='Exclude uncommitted candidates at or above this cultured-baseline marker identity (default: 99.8).')
    quarterly_review_parser.add_argument('--baseline-redundancy-min-query-coverage', type=float, default=baseline_redundancy_coverage,
        help='Minimum query coverage required for baseline redundancy exclusion (default: 95).')
    quarterly_review_parser.add_argument('--baseline-extension-min-identity', type=float, default=baseline_extension_identity,
        help='Minimum same-species cultured-baseline marker identity for baseline-pangenome extension membership (default: 98.65).')
    quarterly_review_parser.add_argument('--baseline-extension-min-query-coverage', type=float, default=baseline_extension_coverage,
        help='Minimum query coverage for baseline-pangenome extension membership (default: 95).')
    quarterly_review_parser.add_argument('--assessment', action='append', default=None, metavar='TSV',
        help='Import a full sequence_assessment.tsv before selection (repeatable; later files supersede earlier rows). Future Performance Reviews store snapshots automatically.')
    quarterly_review_parser.add_argument('--partner-metadata', '--sequencing-metadata', dest='partner_metadata', default=None,
        help='Optional cumulative metadata TSV/CSV used to refresh confirmed already_sequenced status before this round.')
    quarterly_review_parser.add_argument('--tree', default=None,
        help='Latest cumulative Newick tree. Enables marginal patristic-diversity ranking.')
    quarterly_review_parser.add_argument('--alignment', default=None,
        help='Latest cumulative MSA. Recomputes nearest currently available genome context after metadata updates.')
    quarterly_review_parser.add_argument('--from-dir', default=None,
        help='Latest Performance Review output directory; used to locate tree/current_tree.nwk when --tree is omitted.')
    quarterly_review_parser.add_argument('--round-id', default=None,
        help='Stable round label. Default: UTC timestamp such as quarterly_review_20260713T120000Z.')
    quarterly_review_parser.add_argument('--neighbourhood-format', choices=['png'], default='png',
        help='Image format for Quarterly Review local-neighbourhood figures (PNG only).')
    quarterly_review_parser.add_argument('--include-moderate-evidence', action='store_true', default=False,
        help='Allow MODERATE marker-evidence rows into the Quarterly Review. By default they are REVIEW only.')

    label_maker_parser = sub.add_parser(
        'label-maker',
        help='Label Maker: regenerate iTOL annotation files from stored taxonomy.',
        description=(
            'Re-generate all iTOL colour strips (phylum, family, genus, dataset membership)\n'
            'from the taxonomy already stored in the DB. Useful after changing --group-phyla\n'
            'options or after manually editing the database.'
        ),
        formatter_class=_Fmt,
    )
    label_maker_parser.add_argument('--db', required=True,
        help='Path to the BranchManager SQLite database.')
    label_maker_parser.add_argument('-o', '--out', required=True,
        help='Output directory where iTOL files will be written (normally a Filing Cabinet or Performance Review output directory).')
    label_maker_parser.add_argument('--include-datasets', required=False,
        help='Comma-separated list of dataset names to include (default: all datasets in the DB).')
    label_maker_parser.add_argument('--sequence-domain', '--organism-domain', dest='sequence_domain',
        choices=SEQUENCE_DOMAIN_CHOICES, default=None,
        help='Optional domain profile filter for regenerated outputs: bacteria, archaea, fungi, or mixed/all/none.')
    label_maker_parser.add_argument(
        '--group-phyla', dest='group_phyla', action='append', default=None, metavar='SPEC',
        help=(
            'Collapse multiple phyla into one colour in iTOL legends. Repeatable. '
            'Formats: "archaea", "bacteria", "Bacillota,Bacillota_I", "Firmicutes:Bacillota,Bacillota_I".'
        ),
    )
    label_maker_parser.add_argument('--functional', dest='functional', required=False, default=None,
        help=(
            'TSV file mapping sequence IDs to functional attributes. '
            'Generates one iTOL file per column (DATASET_BINARY / DATASET_SIMPLEBAR / DATASET_COLORSTRIP).'
        ),
    )
    label_maker_parser.add_argument('--draft-rumen-functions', dest='draft_rumen_functions',
        action='store_true', default=False,
        help='Auto-generate rumen functional-group iTOL annotation from stored taxonomy.',
    )

    org_chart_parser = sub.add_parser(
        'org-chart',
        help='Org Chart: build a focused tree and iTOL files for a specific taxon.',
        description=(
            'Extract all sequences matching a given taxon from the DB and build a focused\n'
            'phylogenetic tree for that group only.\n\n'
            'Fast path: if --from-dir points to a directory containing current_alignment.fasta\n'
            '(a Filing Cabinet or Performance Review output), sequences are sliced from the pre-built alignment and\n'
            'FastTree is run directly — no MAFFT re-alignment needed (~seconds for hundreds of seqs).\n\n'
            'Slow path: if no existing alignment is found, a full MAFFT + FastTree build is run.\n\n'
            'Taxon formats accepted:\n'
            '  archaea, bacteria          → all sequences at domain level\n'
            '  Bacillota, Bacteroidota    → phylum name (auto-detected)\n'
            '  p__Bacillota               → GTDB-prefixed phylum\n'
            '  f__Lachnospiraceae         → GTDB-prefixed family\n'
            '  g__Ruminococcus            → GTDB-prefixed genus\n'
        ),
        formatter_class=_Fmt,
    )
    org_chart_parser.add_argument('--db', required=True,
        help='Path to the BranchManager SQLite database.')
    org_chart_parser.add_argument('-o', '--out', required=True,
        help='Output directory for the focused-tree results.')
    org_chart_parser.add_argument('--taxon', required=True,
        help='Taxon to extract (see description above for accepted formats).')
    org_chart_parser.add_argument(
        '--rank', required=False, default='auto',
        choices=['auto', 'domain', 'd', 'phylum', 'p', 'class', 'c',
                 'order', 'o', 'family', 'f', 'genus', 'g', 'species', 's'],
        help='Taxonomic rank to filter on. Default "auto" detects the rank from the taxon name or prefix.',
    )
    org_chart_parser.add_argument('--from-dir', dest='from_dir', required=False, default=None,
        help=(
            'Existing Filing Cabinet or Performance Review output directory containing current_alignment.fasta. '
            'Enables the fast path (sequences extracted from the existing MSA; no re-alignment).'
        ),
    )
    org_chart_parser.add_argument('--ref', required=False,
        help='Reference FASTA for orientation correction (slow-path full build only).')
    org_chart_parser.add_argument('--anchors', required=False, default=None,
        help='Custom reference anchor FASTA. Defaults to bundled anchors (26 NCBI RefSeq sequences).')
    org_chart_parser.add_argument('--threads', type=int, default=4,
        help='CPU threads for FastTree / MAFFT (default: 4).')
    org_chart_parser.add_argument('--min-seqs', dest='min_seqs', type=int, default=3,
        help='Minimum sequences required to proceed with tree building (default: 3).')
    org_chart_parser.add_argument('--no-tree', dest='no_tree', action='store_true', default=False,
        help='Skip tree building; only write taxonomy TSV, sequence list, and iTOL colour files.')
    org_chart_parser.add_argument(
        '--group-phyla', dest='group_phyla', action='append', default=None, metavar='SPEC',
        help=(
            'Collapse multiple phyla into one colour in iTOL legends. Repeatable. '
            'Formats: "archaea", "bacteria", "Bacillota,Bacillota_I", "Firmicutes:Bacillota,Bacillota_I".'
        ),
    )
    org_chart_parser.add_argument('--functional', dest='functional', required=False, default=None,
        help=(
            'TSV file mapping sequence IDs to functional attributes. Header row required; first column = sequence ID; subsequent columns = functional attributes. '
            'Generates one iTOL file per column (DATASET_BINARY / DATASET_SIMPLEBAR / DATASET_COLORSTRIP).'
        ),
    )
    org_chart_parser.add_argument('--draft-rumen-functions', dest='draft_rumen_functions',
        action='store_true', default=False,
        help='Auto-generate rumen functional-group iTOL annotation from stored taxonomy.',
    )

    paper_trail_parser = sub.add_parser(
        'paper-trail',
        help='Paper Trail and Merge Meeting: convert AB1 reads and assemble primer sequences.',
        description=(
            'Process Sanger chromatogram reads before Performance Review.\n\n'
            'Inputs may be AB1/ABI files, FASTA, or FASTQ. AB1 files are base-called from '
            'PBAS/PCON tags, quality-trimmed, oriented by primer direction, and optionally '
            'assembled into one consensus 16S sequence per isolate. For example, 27F reads '
            'are kept forward and 907R reads are reverse-complemented before overlap assembly.\n\n'
            'If --sample-map is omitted, sequence IDs and primer names are '
            'inferred from filenames such as Iso001_27F.ab1 and Iso001_907R.ab1.'
        ),
        formatter_class=_Fmt,
    )
    paper_trail_parser.add_argument(
        '--input', nargs='+', required=False, default=[],
        help='AB1/ABI/FASTA/FASTQ files or directories to process. Directories are searched recursively by default. Optional when --sample-map lists the read files.',
    )
    paper_trail_parser.add_argument('-o', '--out', required=True,
        help='Output directory for assembled.fasta, read_qc.tsv, visual reports, and assembly_report.tsv.')
    paper_trail_parser.add_argument(
        '--sample-map', required=False, default=None,
        help=(
            'One CSV/TSV per partner batch. Prefer one row per read with sequence_id, read_file, '
            'primer, direction, and processing_mode. Multiple rows may share a sequence_id when '
            'different primer reads belong to one isolate. Relative paths are resolved next to the map.'
        ),
    )
    paper_trail_parser.add_argument(
        '--primer', dest='primers', action='append', default=None,
        help='Primer name to recognise in filenames (repeatable). Defaults include common 16S primers such as 27F, 907R, and 1492R.',
    )
    paper_trail_parser.add_argument(
        '--primer-sequence', dest='primer_sequences', action='append', default=None, metavar='NAME=SEQUENCE',
        help='Primer oligonucleotide used for confident leading-primer removal (repeatable). Built-in common 16S primers are used by default.',
    )
    paper_trail_parser.add_argument('--trim-primers', dest='trim_primers',
        action=argparse.BooleanOptionalAction, default=True,
        help='Remove confidently matched primer sequence from the retained read (default: enabled).')
    paper_trail_parser.add_argument('--secondary-peak-ratio', dest='secondary_peak_ratio', type=float, default=paper_trail_policy.secondary_peak_ratio,
        help='Secondary/called chromatogram peak ratio considered mixed (default: 0.33).')
    paper_trail_parser.add_argument('--max-mixed-peak-percent', dest='max_mixed_peak_percent', type=float, default=paper_trail_policy.max_mixed_peak_percent,
        help='Maximum percent retained high-quality positions with mixed peaks before read QC failure (default: 15).')
    paper_trail_parser.add_argument('--mixed-peak-min-quality', dest='mixed_peak_min_quality', type=int, default=paper_trail_policy.mixed_peak_min_quality,
        help='Minimum Phred score for a mixed-peak position to count (default: 20).')
    paper_trail_parser.add_argument('--screen-ref', dest='screen_ref', default=None,
        help='Optional marker reference FASTA for independent per-primer taxonomy concordance screening; disagreements require manual review and retain the read-level assignments.')
    paper_trail_parser.add_argument('--screen-taxa', dest='screen_taxa', default=None,
        help='Optional taxonomy table corresponding to --screen-ref; FASTA header taxonomy is used when omitted.')
    paper_trail_parser.add_argument('--threads', type=int, default=4,
        help='CPU threads for optional primer-read taxonomy screening (default: 4).')
    paper_trail_parser.add_argument('--min-quality', dest='min_quality', type=int, default=paper_trail_policy.min_quality,
        help='Phred cutoff for Mott-style end trimming (default: 20).')
    paper_trail_parser.add_argument('--min-length', dest='min_length', type=int, default=paper_trail_policy.min_final_length,
        help='Minimum final sequence length to write to assembled.fasta (default: 800 bp).')
    paper_trail_parser.add_argument('--min-read-length', dest='min_read_length', type=int, default=paper_trail_policy.min_read_length,
        help='Minimum trimmed read length allowed to contribute to an assembly or best-read decision (default: 300 bp).')
    paper_trail_parser.add_argument('--min-mean-quality', dest='min_mean_quality', type=float, default=paper_trail_policy.min_mean_quality,
        help='Minimum mean Phred score after trimming/masking for read and final QC (default: 20).')
    paper_trail_parser.add_argument('--mask-quality', dest='mask_quality', type=int, default=paper_trail_policy.mask_quality,
        help='Mask internal bases below this Phred score to N before assembly (default: 20).')
    paper_trail_parser.add_argument('--max-read-expected-errors', dest='max_read_expected_errors', type=float, default=paper_trail_policy.max_read_expected_errors,
        help='Maximum expected base-call errors allowed per retained read (default: 8).')
    paper_trail_parser.add_argument('--max-output-expected-errors', dest='max_output_expected_errors', type=float, default=paper_trail_policy.max_output_expected_errors,
        help='Maximum expected base-call errors allowed in the final isolate sequence (default: 5).')
    paper_trail_parser.add_argument('--warn-n-percent', dest='warn_n_percent', type=float, default=paper_trail_policy.warn_n_percent,
        help='Percent N above which a passing sequence requires manual review (default: 3.0).')
    paper_trail_parser.add_argument('--max-n-percent', dest='max_n_percent', type=float, default=paper_trail_policy.max_n_percent,
        help='Maximum percent N allowed after masking before QC failure (default: 5.0).')
    paper_trail_parser.add_argument('--warn-internal-lowq-run', dest='warn_internal_low_quality_run', type=int, default=paper_trail_policy.warn_internal_low_quality_run,
        help='Internal low-quality run length above which manual review is required (default: 5 bp).')
    paper_trail_parser.add_argument('--max-internal-lowq-run', dest='max_internal_low_quality_run', type=int, default=paper_trail_policy.max_internal_low_quality_run,
        help='Maximum internal low-quality/ambiguous run length before read failure (default: 20 bp).')
    paper_trail_parser.add_argument('--max-conflict-density', dest='max_conflict_density', type=float, default=paper_trail_policy.max_conflict_density,
        help='Maximum overlap conflicts per 100 final bases before final QC failure (default: 1.0).')
    paper_trail_parser.add_argument('--max-ambiguous-overlap-conflicts-without-review', dest='max_ambiguous_overlap_conflicts_without_review', type=int, default=paper_trail_policy.max_ambiguous_overlap_conflicts_without_review,
        help='Unresolved overlap conflicts tolerated as N before manual review is required (default: 2).')
    paper_trail_parser.add_argument('--quality-difference', dest='quality_difference', type=int, default=paper_trail_policy.quality_difference,
        help='Minimum Phred-scaled posterior odds required to resolve a conflicting overlap base (default: 10).')
    paper_trail_parser.add_argument('--allow-missing-quality', dest='allow_missing_quality',
        action='store_true', default=False,
        help='Allow AB1 reads missing PCON quality scores to pass with warnings. By default they fail QC.')
    paper_trail_parser.add_argument('--min-overlap', dest='min_overlap', type=int, default=paper_trail_policy.min_overlap,
        help='Minimum overlap length for assembling multiple primer reads (default: 40 bp).')
    paper_trail_parser.add_argument('--min-overlap-identity', dest='min_overlap_identity', type=float, default=paper_trail_policy.min_overlap_identity,
        help='Minimum overlap identity for assembly, 0-1 (default: 0.95).')
    paper_trail_parser.add_argument('--assemble', dest='assemble',
        action=argparse.BooleanOptionalAction, default=True,
        help='Assemble multiple reads per sequence_id when possible (default). Use --no-assemble to keep the best read.')
    paper_trail_parser.add_argument('--recursive', dest='recursive',
        action=argparse.BooleanOptionalAction, default=True,
        help='Search input directories recursively (default).')
    paper_trail_parser.add_argument('--max-report-image-height', dest='max_report_image_height',
        type=int, default=2400,
        help='Maximum height in pixels for each visual-report PNG; larger reports are split into numbered pages (minimum: 600, default: 2400).')

    mailroom_parser = sub.add_parser(
        'mailroom',
        help='Mailroom: inventory an AB1 delivery and build its batch map.',
        description=(
            'Reconcile a directory of AB1/ABI chromatograms with supplier metadata and write '
            'a validated one-row-per-read ab1_map.tsv. Primer names are taken from supplier '
            'metadata, embedded ABIF fields, or explicit batch-level forward/reverse settings. '
            'Unresolved primers are reported for review and are never silently inferred as fact.'
        ),
        formatter_class=_Fmt,
    )
    mailroom_parser.add_argument('--read-dir', required=True,
        help='Directory containing the AB1/ABI files for one partner batch.')
    mailroom_parser.add_argument('--metadata', required=True,
        help='Supplier CSV/TSV with sequencing/read ID, isolate ID, and forward/reverse direction. Primer and processing-mode columns are optional.')
    mailroom_parser.add_argument('--dataset', required=True,
        help='Stable batch label written to every map row, for example UoG_01.')
    mailroom_parser.add_argument('--forward-primer', default=None,
        help='Confirmed primer name applied to forward rows lacking a supplied/embedded primer.')
    mailroom_parser.add_argument('--reverse-primer', default=None,
        help='Confirmed primer name applied to reverse rows lacking a supplied/embedded primer.')
    mailroom_parser.add_argument('--processing-mode', choices=['auto', 'assemble', 'best_read'], default='auto',
        help='Default read handling. Auto selects assemble for multi-read isolates and best_read for single-read isolates.')
    mailroom_parser.add_argument('--recursive', action=argparse.BooleanOptionalAction, default=True,
        help='Search the AB1 directory recursively (default: enabled).')
    mailroom_parser.add_argument('-o', '--out', required=True,
        help='Output directory for ab1_map.tsv, AB1 inventory, discrepancy report, and summary.')

    interview_parser = sub.add_parser(
        'interview',
        help='Interview: run standalone AB1 conversion, assembly, and QC from Mailroom output.',
        description=(
            'Run the same versioned AB1 QC policy used by Assistant, starting from a completed '
            'Mailroom output directory or its ab1_map.tsv. Interview writes converted and assembled '
            'FASTA files, QC tables, chromatogram and assembly PNGs, manual-review records, and '
            'resequencing recommendations. It does not read or modify a BranchManager project database.'
        ),
        formatter_class=_Fmt,
    )
    interview_parser.add_argument(
        '--mailroom', required=True,
        help='Mailroom output directory, or its ab1_map.tsv file.',
    )
    interview_parser.add_argument(
        '--primer-sequence', dest='primer_sequences', action='append', default=None,
        metavar='NAME=SEQUENCE',
        help='Confirmed primer oligonucleotide for leading-primer removal (repeatable).',
    )
    interview_parser.add_argument(
        '--trim-primers', dest='trim_primers', action=argparse.BooleanOptionalAction,
        default=True, help='Remove confidently matched leading primer sequences (default: enabled).',
    )
    interview_parser.add_argument(
        '--screen-ref', dest='screen_ref', default=None,
        help='Optional marker reference FASTA for independent per-read taxonomy concordance screening; disagreements trigger manual review and report each read assignment.',
    )
    interview_parser.add_argument(
        '--screen-taxa', dest='screen_taxa', default=None,
        help='Optional taxonomy table corresponding to --screen-ref.',
    )
    interview_parser.add_argument(
        '--threads', type=int, default=4,
        help='CPU threads for optional taxonomy screening (default: 4).',
    )
    interview_parser.add_argument(
        '--allow-missing-quality', action='store_true', default=False,
        help='Allow traces missing PCON quality scores to pass with warnings.',
    )
    interview_parser.add_argument(
        '--allow-mailroom-review', action='store_true', default=False,
        help='Process a REVIEW_REQUIRED Mailroom batch after an explicit review; FAIL batches are always rejected.',
    )
    interview_parser.add_argument(
        '--max-report-image-height', type=int, default=2400,
        help='Maximum PNG height before automatic pagination (minimum: 600, default: 2400).',
    )
    interview_parser.add_argument(
        '-o', '--out', required=True,
        help='Output directory for standalone AB1 QC, converted FASTA files, and visual reports.',
    )

    onboarding_parser = sub.add_parser(
        'onboarding',
        help='Onboarding: validate partner IDs, metadata, and either raw-read ownership or a supplied FASTA.',
        formatter_class=_Fmt,
    )
    onboarding_inputs = onboarding_parser.add_mutually_exclusive_group(required=True)
    onboarding_inputs.add_argument('--sample-map',
        help='CSV/TSV mapping isolate IDs to raw read files for an AB1/primer-read submission.')
    onboarding_inputs.add_argument('--fasta',
        help='Partner-supplied marker FASTA to validate without running Paper Trail/Merge Meeting.')
    onboarding_parser.add_argument('--partner-metadata', required=True,
        help='Cumulative project metadata CSV/TSV. Required for every submission and updated as selection/genome statuses change.')
    onboarding_parser.add_argument('--partner-id', default=None,
        help='Expected partner acronym for this submission, used to validate cumulative-ledger ownership.')
    onboarding_parser.add_argument('--dataset', default=None,
        help='Unique partner batch label, for example QUB_01 or UoG_02.')
    onboarding_parser.add_argument('--read-dir', default=None,
        help='Optional base directory for relative read paths in the sample map.')
    onboarding_parser.add_argument('--primer', dest='primers', action='append', default=None,
        help='Primer column/name to recognise (repeatable).')
    onboarding_parser.add_argument('-o', '--out', required=True,
        help='Output directory for the normalised submission and Onboarding report.')

    status_meeting_parser = sub.add_parser(
        'status-meeting',
        help='Status Meeting: import factual isolate lifecycle changes into the project ledger.',
        formatter_class=_Fmt,
    )
    status_meeting_parser.add_argument('--db', required=True)
    status_meeting_parser.add_argument('--input', required=True,
        help='CSV/TSV with sequence_id, status, and optional detail.')
    status_meeting_parser.add_argument('-o', '--out', required=True)

    records_update_parser = sub.add_parser(
        'records-update',
        help='Records Update: import completed genome/QC/GTDB/ANI evidence.',
        formatter_class=_Fmt,
    )
    records_update_parser.add_argument('--db', required=True)
    records_update_parser.add_argument('--input', required=True,
        help='CSV/TSV with sequence_id, genome_id/accession, genome_status, QC, and optional GTDB/ANI fields.')
    records_update_parser.add_argument('--min-completeness', type=float, default=90.0,
        help='Minimum estimated genome completeness for automatic QC pass (default: 90).')
    records_update_parser.add_argument('--max-contamination', type=float, default=5.0,
        help='Maximum estimated contamination for automatic QC pass (default: 5).')
    records_update_parser.add_argument('-o', '--out', required=True)

    exit_interview_parser = sub.add_parser(
        'exit-interview',
        help='Exit Interview: withdraw sequences from active project analyses with an audit trail.',
        description=(
            'Preview or apply removal of one or more sequences from the project database. '
            'Applying an Exit Interview removes active sequence, taxonomy, distance, assessment, '
            'selection, provenance, status, and metadata records while retaining an immutable '
            'tombstone. Baseline and genome-backed records are protected by default.'
        ),
        formatter_class=_Fmt,
    )
    exit_source = exit_interview_parser.add_mutually_exclusive_group(required=True)
    exit_source.add_argument('--sequence-id', action='append',
        help='Sequence/isolate ID to remove; repeat for multiple IDs.')
    exit_source.add_argument('--input',
        help='CSV/TSV, optionally gzipped, with sequence_id and reason columns.')
    exit_interview_parser.add_argument('--db', required=True)
    exit_interview_parser.add_argument('--reason', default='',
        help='Required reason for direct IDs, or fallback reason for input rows.')
    exit_interview_parser.add_argument('--apply', action='store_true', default=False,
        help='Apply the reviewed plan. Without this flag, Exit Interview is preview-only.')
    exit_interview_parser.add_argument('--allow-baseline', action='store_true', default=False,
        help='Explicitly allow removal from a registered cultured baseline.')
    exit_interview_parser.add_argument('--allow-genome-records', action='store_true', default=False,
        help='Explicitly allow removal of an isolate that has genome evidence.')
    exit_interview_parser.add_argument('--backup', action=argparse.BooleanOptionalAction, default=True,
        help='Write project_before_exit_interview.sqlite before applying (default: enabled).')
    exit_interview_parser.add_argument('-o', '--out', required=True)

    annual_report_parser = sub.add_parser(
        'annual-report',
        help='Annual Report: create a point-in-time project overview, ledgers, and decision-change report.',
        formatter_class=_Fmt,
    )
    annual_report_parser.add_argument('--db', required=True)
    annual_report_parser.add_argument('-o', '--out', required=True)

    it_desk_parser = sub.add_parser(
        'it-desk',
        help='IT Desk: check dependencies, references, database integrity, and output location.',
        formatter_class=_Fmt,
    )
    it_desk_parser.add_argument('--db', default=None)
    it_desk_parser.add_argument('--ref', dest='references', action='append', default=[])
    it_desk_parser.add_argument('--tree-method', choices=['fasttree', 'iqtree', 'iqtree-fast'], default='fasttree')
    it_desk_parser.add_argument('-o', '--out', default='branchmanager_it_desk')
    it_desk_parser.add_argument('--strict', action='store_true', default=False,
        help='Return a non-zero exit status when any required check fails.')

    assistant_parser = sub.add_parser(
        'assistant',
        help='Assistant to the Branch Manager: run Onboarding through the Hiring Panel.',
        formatter_class=_Fmt,
    )
    assistant_input = assistant_parser.add_mutually_exclusive_group(required=True)
    assistant_input.add_argument('--sample-map', help='AB1/primer-read sample map; runs Paper Trail and Merge Meeting.')
    assistant_input.add_argument('--fasta', help='Partner-supplied marker FASTA; bypasses Paper Trail and Merge Meeting.')
    assistant_parser.add_argument('--partner-metadata', required=True)
    assistant_parser.add_argument('--partner-id', default=None,
        help='Expected partner acronym for validation, for example QUB or UoG.')
    assistant_parser.add_argument('--read-dir', default=None)
    assistant_parser.add_argument('--db', required=True)
    assistant_parser.add_argument('--dataset', required=True)
    assistant_parser.add_argument('--ref', required=True)
    assistant_parser.add_argument('--taxa', default=None)
    assistant_parser.add_argument('--ref-name', default='GTDB')
    assistant_parser.add_argument('--alt-ref', action='append', default=[])
    assistant_parser.add_argument('--alt-taxa', action='append', default=[])
    assistant_parser.add_argument('--alt-ref-name', action='append', default=[])
    assistant_parser.add_argument('--baseline-fasta', default=None)
    assistant_parser.add_argument('--baseline-dataset', default='Baseline')
    assistant_parser.add_argument('--baseline-tier', choices=['priority', 'secondary'], default=None)
    assistant_parser.add_argument('--baseline-taxa-assignments', default=None)
    assistant_parser.add_argument('--mwl', default=None)
    assistant_parser.add_argument('--sequence-domain', choices=SEQUENCE_DOMAIN_CHOICES, default='bacteria')
    assistant_parser.add_argument('--threads', type=int, default=4,
        help='CPU threads for vsearch and MAFFT (default: 4).')
    assistant_parser.add_argument('--tree-method', choices=['fasttree', 'iqtree', 'iqtree-fast'], default='fasttree')
    assistant_parser.add_argument('--pangenome-target', type=int, default=selection_target)
    assistant_parser.add_argument('--candidate-set-size', type=int, default=selection_panel_size)
    assistant_parser.add_argument('--baseline-redundancy-identity', type=float, default=baseline_redundancy_identity)
    assistant_parser.add_argument('--baseline-redundancy-min-query-coverage', type=float, default=baseline_redundancy_coverage)
    assistant_parser.add_argument('--baseline-extension-min-identity', type=float, default=baseline_extension_identity)
    assistant_parser.add_argument('--baseline-extension-min-query-coverage', type=float, default=baseline_extension_coverage)
    assistant_parser.add_argument('--primer', dest='primers', action='append', default=None)
    assistant_parser.add_argument('--primer-sequence', dest='primer_sequences', action='append', default=None,
        metavar='NAME=SEQUENCE', help='Primer sequence used for IUPAC-aware trimming; repeatable.')
    assistant_parser.add_argument('--trim-primers', dest='trim_primers',
        action=argparse.BooleanOptionalAction, default=True)
    assistant_parser.add_argument('--min-quality', type=int, default=paper_trail_policy.min_quality,
        help='Phred cutoff for Mott-style end trimming (default: 20).')
    assistant_parser.add_argument('--min-marker-length', dest='min_length', type=int, default=paper_trail_policy.min_final_length,
        help='Minimum final assembled sequence length in bp (default: 800).')
    assistant_parser.add_argument('--min-read-length', type=int, default=paper_trail_policy.min_read_length,
        help='Minimum trimmed primer-read length allowed to contribute to an assembly (default: 300 bp).')
    assistant_parser.add_argument('--min-mean-quality', type=float, default=paper_trail_policy.min_mean_quality,
        help='Minimum mean Phred quality after trimming (default: 20.0).')
    assistant_parser.add_argument('--mask-quality', type=int, default=paper_trail_policy.mask_quality,
        help='Mask internal bases below this Phred score to N (default: 20).')
    assistant_parser.add_argument('--max-read-expected-errors', type=float, default=paper_trail_policy.max_read_expected_errors,
        help='Maximum expected errors per retained read (default: 8.0).')
    assistant_parser.add_argument('--max-output-expected-errors', type=float, default=paper_trail_policy.max_output_expected_errors,
        help='Maximum expected errors in the final assembled sequence (default: 5.0).')
    assistant_parser.add_argument('--warn-n-percent', type=float, default=paper_trail_policy.warn_n_percent,
        help='Percent N above which a passing sequence requires manual review (default: 3.0).')
    assistant_parser.add_argument('--max-n-percent', type=float, default=paper_trail_policy.max_n_percent,
        help='Maximum percent N bases allowed after masking before QC failure (default: 5.0).')
    assistant_parser.add_argument('--warn-internal-lowq-run', dest='warn_internal_low_quality_run', type=int, default=paper_trail_policy.warn_internal_low_quality_run,
        help='Internal low-quality run length above which manual review is required (default: 5 bp).')
    assistant_parser.add_argument('--max-internal-lowq-run', dest='max_internal_low_quality_run', type=int, default=paper_trail_policy.max_internal_low_quality_run,
        help='Maximum internal low-quality/ambiguous run before QC failure (default: 20 bp).')
    assistant_parser.add_argument('--secondary-peak-ratio', type=float, default=paper_trail_policy.secondary_peak_ratio,
        help='Secondary/called peak ratio threshold for mixed-base detection (default: 0.33).')
    assistant_parser.add_argument('--max-mixed-peak-percent', type=float, default=paper_trail_policy.max_mixed_peak_percent,
        help='Maximum percent mixed-peak positions before read QC failure (default: 15.0).')
    assistant_parser.add_argument('--mixed-peak-min-quality', dest='mixed_peak_min_quality', type=int, default=paper_trail_policy.mixed_peak_min_quality,
        help='Minimum Phred score for a mixed-peak position to count (default: 20).')
    assistant_parser.add_argument('--max-conflict-density', type=float, default=paper_trail_policy.max_conflict_density,
        help='Maximum overlap conflicts per 100 final bases before QC failure (default: 1.0).')
    assistant_parser.add_argument('--max-ambiguous-overlap-conflicts-without-review', type=int, default=paper_trail_policy.max_ambiguous_overlap_conflicts_without_review,
        help='Unresolved overlap conflicts tolerated as N before manual review is required (default: 2).')
    assistant_parser.add_argument('--quality-difference', type=int, default=paper_trail_policy.quality_difference,
        help='Phred-scaled posterior odds required to resolve an overlap conflict (default: 10).')
    assistant_parser.add_argument('--allow-missing-quality', action='store_true', default=False,
        help='Allow traces missing PCON quality scores to pass with warnings.')
    assistant_parser.add_argument('--max-report-image-height', type=int, default=2400,
        help='Maximum Paper Trail visual-report PNG height before automatic pagination (minimum: 600, default: 2400).')
    assistant_parser.add_argument('--min-overlap', type=int, default=paper_trail_policy.min_overlap)
    assistant_parser.add_argument('--min-overlap-identity', type=float, default=paper_trail_policy.min_overlap_identity)
    assistant_parser.add_argument('--chimera-ref', default=None,
        help='Curated reference FASTA for UCHIME; defaults to the primary reference.')
    assistant_parser.add_argument('--skip-chimera-check', action='store_true', default=False,
        help='Audited override: continue with all markers marked for review.')
    assistant_parser.add_argument('--marker-qc', default=None,
        help='Marker-QC sidecar for --fasta submissions, such as a reviewed assembly_report.tsv.')
    assistant_parser.add_argument('--marker-review', default=None,
        help='CSV/TSV decisions for PASS_WITH_WARNINGS markers, with sequence_id, decision, reviewer, and optional notes.')
    assistant_parser.add_argument('--accept-unverified-marker-qc', action='store_true', default=False,
        help='Audited acceptance of a partner FASTA without BranchManager marker-QC provenance.')
    assistant_parser.add_argument('-o', '--out', required=True)

    background_check_parser = sub.add_parser(
        'background-check',
        help='Background Check: pre-classify reference collections for reuse.',
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
            '  background_check_summary.txt  Plain-text summary with usage examples\n\n'
            'Example:\n'
            '  branchmanager background-check \\\n'
            '    --dataset hungate16s=/data/hungate.fasta \\\n'
            '    --dataset silva=/data/silva_16s.fasta \\\n'
            '    --ref /data/gtdb_ssu_reps.fna \\\n'
            '    --taxa /data/gtdb_taxonomy.tsv.gz \\\n'
            '    --threads 8 -o background_check_out/\n\n'
            'Then use the output in the Filing Cabinet:\n'
            '  branchmanager filing-cabinet --fasta hungate.fasta \\\n'
            '    --taxa-assignments background_check_out/pipeline_taxonomy.tsv \\\n'
            '    --db project.db --dataset Hungate -o filing_cabinet_out/'
        ),
        formatter_class=_Fmt,
    )
    background_check_parser.add_argument(
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
    background_check_parser.add_argument(
        '--ref', required=True,
        help='Reference FASTA (e.g. GTDB/SILVA reps) used by vsearch for classification.',
    )
    background_check_parser.add_argument(
        '--taxa', required=False, default=None,
        help='Taxonomy table matching IDs in --ref (TSV/CSV, optionally .gz: id<TAB>lineage or id,lineage). '
             'If omitted, taxonomy is parsed directly from reference FASTA headers.',
    )
    background_check_parser.add_argument(
        '-o', '--out', required=True,
        help='Output directory where all classification files will be written.',
    )
    background_check_parser.add_argument(
        '--threads', type=int, default=4,
        help='CPU threads for vsearch (default: 4).',
    )
    background_check_parser.add_argument(
        '--min-identity', dest='min_identity', type=float, default=0.80,
        help='Minimum vsearch alignment identity threshold 0–1 (default: 0.80 = 80%%).',
    )
    background_check_parser.add_argument(
        '--max-hits', dest='max_hits', type=int, default=10,
        help=(
            'Number of candidate hits vsearch collects per query '
            '(--maxaccepts / --maxhits).  The best hit by %% identity is '
            'selected after collection.  Higher values are more thorough but '
            'slower (default: 10).'
        ),
    )
    background_check_parser.add_argument(
        '--max-rejects', dest='max_rejects', type=int, default=256,
        help=(
            'vsearch --maxrejects: maximum number of non-matching candidate '
            'sequences examined before giving up on a query.  Raising this '
            'above vsearch\'s default of 32 helps classify ambiguous sequences '
            '(e.g. those with many N\'s) (default: 256).'
        ),
    )
    background_check_parser.add_argument(
        '--low-confidence-threshold', dest='low_confidence_threshold',
        type=float, default=0.97,
        help=(
            'Identity threshold below which a classified hit is flagged as low-confidence '
            'and written to *_low_confidence.tsv for manual review. '
            '0–1 (default: 0.97 = 97%%, the traditional species-level cutoff).'
        ),
    )

    return parser


def cmd_background_check(args):
    """Run the reusable Background Check classification workflow."""
    from branchmanager.pipeline import background_check as _background_check_module

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
                f"[BACKGROUND CHECK] Cannot parse dataset spec at position {i}: '{token}'. "
                "Expected 'NAME FASTA' pairs or 'NAME=FASTA' pairs."
            )

    if not datasets:
        raise SystemExit("[BACKGROUND CHECK] At least one dataset is required. "
                         "Use: --dataset NAME /path/to/file.fasta")

    # Validate FASTA paths exist
    for name, fasta in datasets:
        if not os.path.exists(fasta):
            raise SystemExit(
                f"[BACKGROUND CHECK] FASTA file not found for dataset '{name}': {fasta}"
            )

    log.info("[BACKGROUND CHECK] Starting pre-classification for %d dataset(s)", len(datasets))

    pipeline_tsv = _background_check_module.run_background_check(
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

    log.info("[BACKGROUND CHECK] Done. Pipeline-ready taxonomy → %s", pipeline_tsv)
    print(
        f"[background-check] Done.\n"
        f"  Pipeline taxonomy : {pipeline_tsv}\n"
        f"  Summary           : {os.path.join(outdir, 'background_check_summary.txt')}\n\n"
        f"Use in the Filing Cabinet:\n"
        f"  branchmanager filing-cabinet --fasta <fasta> --taxa-assignments {pipeline_tsv} "
        f"--db project.db --dataset <name> -o filing_cabinet_out/"
    )


def _resolve_mailroom_interview_input(
    source: str | os.PathLike,
    *,
    allow_review: bool = False,
) -> tuple[str, str, str]:
    """Resolve and validate the immutable Mailroom hand-off used by Interview."""
    supplied = Path(source).expanduser().resolve()
    if supplied.is_dir():
        map_path = supplied / 'ab1_map.tsv'
        summary_path = supplied / 'mailroom_summary.json'
    elif supplied.is_file():
        map_path = supplied
        summary_path = supplied.parent / 'mailroom_summary.json'
    else:
        raise SystemExit(f'[interview] Mailroom input does not exist: {supplied}')

    if map_path.name != 'ab1_map.tsv' or not map_path.is_file():
        raise SystemExit(
            f'[interview] Expected a Mailroom ab1_map.tsv, but it was not found at {map_path}'
        )
    if not summary_path.is_file():
        raise SystemExit(
            f'[interview] Expected the accompanying Mailroom summary at {summary_path}'
        )

    try:
        summary = json.loads(summary_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise SystemExit(f'[interview] Could not read {summary_path}: {exc}') from exc
    status = str(summary.get('status', '')).strip().upper()
    if status not in {'PASS', 'REVIEW_REQUIRED', 'FAIL'}:
        raise SystemExit(
            f'[interview] Mailroom summary has an invalid or missing status: {status or "<blank>"}'
        )
    if status == 'FAIL':
        raise SystemExit(
            '[interview] Mailroom status is FAIL. Correct the Mailroom report before AB1 QC.'
        )
    if status == 'REVIEW_REQUIRED' and not allow_review:
        raise SystemExit(
            '[interview] Mailroom status is REVIEW_REQUIRED. Resolve its findings or repeat with '
            '--allow-mailroom-review after documenting the review.'
        )

    with open(map_path, newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fields = {str(field).strip() for field in (reader.fieldnames or [])}
        required = {
            'sequence_id', 'dataset', 'read_file', 'primer', 'direction',
            'processing_mode', 'primer_assignment',
        }
        missing = sorted(required - fields)
        if missing:
            raise SystemExit(
                '[interview] Mailroom ab1_map.tsv is missing required column(s): '
                + ', '.join(missing)
            )
        if next(reader, None) is None:
            raise SystemExit('[interview] Mailroom ab1_map.tsv contains no reads.')

    return str(map_path), str(summary_path), status



def _infer_qc_dataset_from_path(value) -> str:
    raw = str(value or '').strip()
    if not raw:
        return ''
    try:
        parts = [part for part in Path(raw).parts if part not in {'', '/'}]
    except TypeError:
        return ''
    for index, part in enumerate(parts):
        if part == 'mailroom' and index > 0:
            return parts[index - 1]
    for part in reversed(parts):
        if part.startswith('01_interview_') and len(part) > len('01_interview_'):
            return part[len('01_interview_'):]
    for index, part in enumerate(parts):
        if part in {'01_onboarding', '02_paper_trail_merge_meeting'} and index > 0:
            return parts[index - 1]
    return ''

def cmd_interview(args):
    """Run standalone AB1 QC from a validated Mailroom batch."""
    sample_map, summary, status = _resolve_mailroom_interview_input(
        args.mailroom,
        allow_review=bool(getattr(args, 'allow_mailroom_review', False)),
    )
    args.sample_map = sample_map
    args.mailroom_summary = summary
    args.mailroom_status = status
    args.input = []
    cmd_paper_trail(args)


def cmd_paper_trail(args):
    """Handler for standalone or Assistant Paper Trail chromatogram QC."""
    is_interview = getattr(args, 'command', None) == 'interview'
    command_name = 'interview' if is_interview else 'paper-trail'
    workflow_name = 'interview' if is_interview else 'paper_trail'
    outdir = args.out
    os.makedirs(outdir, exist_ok=True)
    _configure_logging(outdir)
    primers = getattr(args, 'primers', None) or _paper_trail_module.DEFAULT_PRIMERS
    inputs = getattr(args, 'input', None) or []
    sample_map = getattr(args, 'sample_map', None)
    if not inputs and not sample_map:
        raise SystemExit(f'[{command_name}] Provide --input and/or --sample-map.')
    primer_sequences = dict(_paper_trail_module.DEFAULT_PRIMER_SEQUENCES)
    for specification in getattr(args, 'primer_sequences', None) or []:
        if '=' not in str(specification):
            raise SystemExit(f'[{command_name}] --primer-sequence must use NAME=SEQUENCE.')
        name, sequence = str(specification).split('=', 1)
        sequence = ''.join(sequence.split()).upper()
        if not name.strip() or len(sequence) < 12:
            raise SystemExit(
                f'[{command_name}] --primer-sequence requires a name and at least 12 bases.'
            )
        primer_sequences[name.strip().upper()] = sequence
    manifest = RunManifest(outdir, workflow_name)
    dataset_hint = (
        _infer_qc_dataset_from_path(sample_map)
        or _infer_qc_dataset_from_path(getattr(args, 'mailroom_summary', None))
        or _infer_qc_dataset_from_path(outdir)
    )
    if dataset_hint:
        manifest.data['dataset'] = dataset_hint
        manifest.write()
    for source, role in (
        (sample_map, 'sample_map'),
        (getattr(args, 'mailroom_summary', None), 'mailroom_summary'),
        (getattr(args, 'screen_ref', None), 'screen_reference'),
        (getattr(args, 'screen_taxa', None), 'screen_taxonomy'),
    ):
        if source:
            manifest.add_input(source, role=role)
    if getattr(args, 'mailroom_status', 'PASS') != 'PASS':
        manifest.warn(
            'Interview proceeded after explicit acceptance of a REVIEW_REQUIRED Mailroom batch.'
        )
    for source in inputs:
        manifest.add_input(
            source,
            role='raw_read' if Path(source).is_file() else 'raw_read_directory',
        )
    manifest.add_stage('paper_trail', 'RUNNING', detail='base calls, chromatogram evidence, and trimming')
    policy = _paper_trail_module.DEFAULT_QC_POLICY

    def setting(name, default):
        value = getattr(args, name, None)
        return default if value is None else value

    try:
        outputs = _paper_trail_module.run_paper_trail(
            inputs,
            outdir,
            sample_map=sample_map,
            primers=primers,
            primer_sequences=primer_sequences,
            trim_primers=bool(getattr(args, 'trim_primers', True)),
            secondary_peak_ratio=float(setting('secondary_peak_ratio', policy.secondary_peak_ratio)),
            max_mixed_peak_percent=float(setting('max_mixed_peak_percent', policy.max_mixed_peak_percent)),
            mixed_peak_min_quality=int(setting('mixed_peak_min_quality', policy.mixed_peak_min_quality)),
            screen_ref=getattr(args, 'screen_ref', None),
            screen_taxa=getattr(args, 'screen_taxa', None),
            threads=int(setting('threads', 4)),
            min_quality=int(setting('min_quality', policy.min_quality)),
            min_length=int(setting('min_length', policy.min_final_length)),
            min_read_length=int(setting('min_read_length', policy.min_read_length)),
            min_mean_quality=float(setting('min_mean_quality', policy.min_mean_quality)),
            mask_quality=int(setting('mask_quality', policy.mask_quality)),
            max_read_expected_errors=float(setting('max_read_expected_errors', policy.max_read_expected_errors)),
            max_output_expected_errors=float(setting('max_output_expected_errors', policy.max_output_expected_errors)),
            warn_n_percent=float(setting('warn_n_percent', policy.warn_n_percent)),
            max_n_percent=float(setting('max_n_percent', policy.max_n_percent)),
            warn_internal_low_quality_run=int(setting('warn_internal_low_quality_run', policy.warn_internal_low_quality_run)),
            max_internal_low_quality_run=int(setting('max_internal_low_quality_run', policy.max_internal_low_quality_run)),
            max_conflict_density=float(setting('max_conflict_density', policy.max_conflict_density)),
            max_ambiguous_overlap_conflicts_without_review=int(setting('max_ambiguous_overlap_conflicts_without_review', policy.max_ambiguous_overlap_conflicts_without_review)),
            quality_difference=int(setting('quality_difference', policy.quality_difference)),
            allow_missing_quality=bool(getattr(args, 'allow_missing_quality', False)),
            min_overlap=int(setting('min_overlap', policy.min_overlap)),
            min_overlap_identity=float(setting('min_overlap_identity', policy.min_overlap_identity)),
            assemble=bool(getattr(args, 'assemble', True)),
            recursive=bool(getattr(args, 'recursive', True)),
            max_report_image_height=int(getattr(args, 'max_report_image_height', 2400) or 2400),
        )
        manifest.add_stage('paper_trail', 'COMPLETE')
        manifest.add_stage('merge_meeting', 'COMPLETE', detail='per-isolate assembly or best-read selection')
        for key in (
            'raw_fasta', 'trimmed_fasta', 'assembled_fasta', 'read_qc_tsv',
            'per_base_error_tsv', 'assembly_tsv', 'assembly_placements_tsv',
            'recommendations_tsv', 'marker_review_template_tsv', 'qc_policy_tsv',
            'failed_final_fasta',
            'failed_read_fasta', 'failed_manifest_tsv', 'failed_read_manifest_tsv',
            'failed_qc_guide', 'visual_manifest_tsv', 'summary',
        ):
            manifest.add_output(outputs[key], role=key)
        for key in ('read_error_pngs', 'chromatogram_pngs', 'assembly_pngs'):
            for path in outputs[key]:
                manifest.add_output(path, role=key)
        manifest.finish('COMPLETE')
    except Exception as exc:
        manifest.finish('FAILED', error=exc)
        raise
    logging.getLogger(__name__).info(
        "[%s] Final assembled FASTA: %s", command_name.upper(), outputs['assembled_fasta'],
    )
    completion = (
        'Standalone AB1 QC complete.' if is_interview
        else 'Paper Trail + Merge Meeting complete.'
    )
    next_step = (
        'No project database was read or changed. After resolving manual-review decisions, '
        'provide assembled.fasta and assembly_report.tsv to Assistant as a reviewed FASTA submission.'
        if is_interview else
        'Use in Performance Review:\n'
        f"  branchmanager performance-review --input {outputs['assembled_fasta']} --partner-metadata <metadata.tsv> ..."
    )
    print(
        f"[{command_name}] {completion}\n"
        f"  Assembled FASTA : {outputs['assembled_fasta']}\n"
        f"  Trimmed reads   : {outputs['trimmed_fasta']}\n"
        f"  Read QC         : {outputs['read_qc_tsv']}\n"
        f"  Per-base errors : {outputs['per_base_error_tsv']}\n"
        f"  Assembly report : {outputs['assembly_tsv']}\n"
        f"  Read placements : {outputs['assembly_placements_tsv']}\n\n"
        f"  Resequence list : {outputs['recommendations_tsv']}\n"
        f"  Review template : {outputs['marker_review_template_tsv']}\n"
        f"  QC policy       : {outputs['qc_policy_tsv']}\n"
        f"  Failed QC seqs  : {outputs['failed_qc_dir']}\n\n"
        f"  Failed isolates : {outputs['failed_manifest_tsv']}\n"
        f"  Failed reads    : {outputs['failed_read_manifest_tsv']}\n\n"
        f"  Failure guide   : {outputs['failed_qc_guide']}\n\n"
        f"  Read visuals    : {len(outputs['read_error_pngs'])} page(s)\n"
        f"  Chromatograms   : {len(outputs['chromatogram_pngs'])} page(s)\n"
        f"  Assembly visuals: {len(outputs['assembly_pngs'])} page(s)\n"
        f"  Visual manifest : {outputs['visual_manifest_tsv']}\n\n"
        f"{next_step}"
    )


def cmd_mailroom(args):
    from branchmanager.mailroom import prepare_ab1_map

    manifest = RunManifest(args.out, 'mailroom')
    manifest.add_input(args.read_dir, role='ab1_directory')
    manifest.add_input(args.metadata, role='supplier_metadata')
    try:
        result = prepare_ab1_map(
            args.read_dir,
            args.metadata,
            args.out,
            dataset=args.dataset,
            forward_primer=args.forward_primer,
            reverse_primer=args.reverse_primer,
            processing_mode=args.processing_mode,
            recursive=args.recursive,
        )
        for role in ('ab1_map', 'inventory', 'report', 'summary'):
            manifest.add_output(result[role], role=role)
        manifest.add_stage(
            'mailroom', result['status'],
            detail=f"{result['mapped_reads']} mapped reads; {result['errors']} errors",
        )
        manifest.finish('COMPLETE' if result['status'] == 'PASS' else 'FAILED')
    except Exception as exc:
        manifest.finish('FAILED', error=exc)
        raise
    print(
        f"[mailroom] {result['status']}: {result['physical_ab1_files']} AB1 file(s), "
        f"{result['mapped_reads']} mapped read(s), {result['isolates']} metadata isolate(s), "
        f"{result['mapped_isolates']} with mapped reads, "
        f"{result['errors']} error(s), {result['unresolved_primers']} unresolved primer(s).\n"
        f"  Batch map : {result['ab1_map']}\n"
        f"  Inventory : {result['inventory']}\n"
        f"  Report    : {result['report']}\n"
        f"  Summary   : {result['summary']}"
    )
    if result['status'] != 'PASS':
        raise SystemExit(2)


def cmd_onboarding(args):
    from branchmanager.onboarding import validate_submission, write_onboarding_outputs
    from branchmanager.pipeline.paper_trail import DEFAULT_PRIMERS

    manifest = RunManifest(args.out, 'onboarding')
    sample_map = getattr(args, 'sample_map', None)
    fasta = getattr(args, 'fasta', None)
    manifest.add_input(sample_map or fasta, role='sample_map' if sample_map else 'marker_fasta')
    if args.partner_metadata:
        manifest.add_input(args.partner_metadata, role='partner_metadata')
    try:
        result = validate_submission(
            sample_map,
            fasta=fasta,
            partner_metadata=args.partner_metadata,
            read_dir=args.read_dir,
            primers=args.primers or DEFAULT_PRIMERS,
            expected_partner_id=getattr(args, 'partner_id', None),
            dataset=getattr(args, 'dataset', None),
        )
        outputs = write_onboarding_outputs(args.out, result)
        for role, path in outputs.items():
            manifest.add_output(path, role=role)
        manifest.add_stage('onboarding', result['status'], detail=f"{result['isolates']} isolates; {result['errors']} errors")
        manifest.finish('COMPLETE' if result['status'] == 'PASS' else 'FAILED')
    except Exception as exc:
        manifest.finish('FAILED', error=exc)
        raise
    print(
        f"[onboarding] {result['status']}: {result['isolates']} isolate(s), "
        f"input={result['input_type']}, {result['read_files']} read file(s), {result['errors']} error(s).\n"
        f"  Normalised map : {outputs['normalised']}\n"
        f"  Validation     : {outputs['report']}"
        + (f"\n  Normalised FASTA: {outputs['fasta']}" if outputs.get('fasta') else '')
    )
    if result['status'] != 'PASS':
        raise SystemExit(2)


def _cmd_project_import(args, *, genomes: bool):
    import fcntl
    import sqlite3
    from branchmanager.project_state import import_genome_results, import_status_updates, write_import_report

    workflow = 'records_update' if genomes else 'status_meeting'
    manifest = RunManifest(args.out, workflow)
    manifest.add_input(args.input, role='genome_results' if genomes else 'status_updates')
    original = Path(args.db).expanduser().resolve()
    if not original.is_file():
        raise SystemExit(f'[records-update] Project database does not exist: {original}')
    manifest.add_input(original, role='project_database_before')
    staged_path = Path(args.out).resolve() / f'.{workflow}_project.sqlite'
    lock_path = original.with_name(original.name + '.lock')
    with open(lock_path, 'a+') as lock_handle:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        try:
            staged_path.unlink(missing_ok=True)
            with sqlite3.connect(str(original)) as source, sqlite3.connect(str(staged_path)) as target:
                source.backup(target)
            db = Database(str(staged_path))
            db.initialise()
            rows = import_genome_results(
                db, args.input,
                min_completeness=float(getattr(args, 'min_completeness', 90.0)),
                max_contamination=float(getattr(args, 'max_contamination', 5.0)),
            ) if genomes else import_status_updates(db, args.input)
            report_name = 'records_update_report.tsv' if genomes else 'status_meeting_report.tsv'
            report = write_import_report(Path(args.out) / report_name, rows)
            manifest.add_output(report, role='import_report')
            rejected = sum(row.get('result') == 'REJECTED' for row in rows)
            manifest.add_stage(workflow, 'COMPLETE' if rejected == 0 else 'COMPLETE_WITH_WARNINGS', detail=f'{len(rows) - rejected} imported; {rejected} rejected')
            manifest.finish('COMPLETE' if rejected == 0 else 'COMPLETE_WITH_WARNINGS')
            db.record_project_run(
                f'{workflow.replace("_", "-")}:{manifest.data["started_at"]}', workflow, manifest.data['status'],
                manifest_path=str(manifest.json_path), started_at=manifest.data['started_at'],
                completed_at=manifest.data['completed_at'],
            )
            publish_tmp = original.with_name(original.name + '.publishing')
            publish_tmp.unlink(missing_ok=True)
            with sqlite3.connect(str(staged_path)) as source, sqlite3.connect(str(publish_tmp)) as target:
                source.backup(target)
                if target.execute('PRAGMA integrity_check').fetchone()[0] != 'ok':
                    raise RuntimeError('staged Records Update database failed integrity_check')
            os.replace(publish_tmp, original)
            staged_path.unlink(missing_ok=True)
        except BaseException as exc:
            manifest.finish('FAILED', error=exc)
            raise
        finally:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)
    command_label = 'records-update' if genomes else 'status-meeting'
    print(f'[{command_label}] {len(rows) - rejected} row(s) imported; {rejected} rejected.\n  Report: {report}')


def cmd_exit_interview(args):
    """Preview or atomically apply an audited sequence withdrawal."""
    import fcntl
    import sqlite3
    from branchmanager.personnel import (
        load_exit_requests, write_departing_fasta, write_exit_interview_report,
    )

    output = Path(args.out).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    original = Path(args.db).expanduser().resolve()
    if not original.is_file():
        raise SystemExit(f'[exit-interview] Project database does not exist: {original}')
    requests = load_exit_requests(
        args.sequence_id, input_path=args.input, default_reason=args.reason,
    )
    manifest = RunManifest(output, 'exit_interview')
    manifest.add_input(original, role='project_database_before')
    if args.input:
        manifest.add_input(args.input, role='exit_interview_requests')

    if not args.apply:
        db = Database(str(original))
        rows = db.plan_sequence_removals(
            requests,
            allow_baseline=bool(args.allow_baseline),
            allow_genome_records=bool(args.allow_genome_records),
        )
        report = write_exit_interview_report(output / 'exit_interview_plan.tsv', rows)
        manifest.add_output(report, role='exit_interview_plan')
        ready = sum(row.get('status') == 'READY' for row in rows)
        blocked = len(rows) - ready
        status = 'COMPLETE' if blocked == 0 else 'COMPLETE_WITH_WARNINGS'
        manifest.add_stage('exit_interview', status, detail=f'{ready} ready; {blocked} blocked')
        manifest.finish(status)
        print(
            f'[exit-interview] Preview only: {ready} sequence(s) ready; {blocked} blocked.\n'
            f'  Plan: {report}\n'
            '  Review the plan, then repeat the command with --apply.'
        )
        return

    staged_path = output / '.exit_interview_project.sqlite'
    lock_path = original.with_name(original.name + '.lock')
    blocked_rows = []
    report = ''
    removed_count = 0
    with open(lock_path, 'a+') as lock_handle:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        try:
            staged_path.unlink(missing_ok=True)
            with sqlite3.connect(str(original)) as source, sqlite3.connect(str(staged_path)) as target:
                source.backup(target)
            db = Database(str(staged_path))
            db.initialise()
            rows = db.plan_sequence_removals(
                requests,
                allow_baseline=bool(args.allow_baseline),
                allow_genome_records=bool(args.allow_genome_records),
            )
            blocked_rows = [row for row in rows if row.get('status') != 'READY']
            if blocked_rows:
                report = write_exit_interview_report(output / 'exit_interview_plan.tsv', rows)
                manifest.add_output(report, role='exit_interview_plan')
                manifest.add_stage(
                    'exit_interview', 'COMPLETE_WITH_WARNINGS',
                    detail=f'0 removed; {len(blocked_rows)} blocked',
                )
                manifest.finish('COMPLETE_WITH_WARNINGS')
            else:
                archive = write_departing_fasta(output / 'departing_sequences.fasta', rows)
                manifest.add_output(archive, role='departing_sequences')
                if args.backup:
                    backup = output / 'project_before_exit_interview.sqlite'
                    backup.unlink(missing_ok=True)
                    with sqlite3.connect(str(original)) as source, sqlite3.connect(str(backup)) as target:
                        source.backup(target)
                    manifest.add_output(backup, role='project_database_backup')
                db.apply_sequence_removals(
                    rows, source_request=str(Path(args.input).resolve()) if args.input else 'command_line',
                )
                report = write_exit_interview_report(output / 'exit_interview_report.tsv', rows)
                manifest.add_output(report, role='exit_interview_report')
                with db.connect() as conn:
                    integrity = conn.execute('PRAGMA integrity_check').fetchone()[0]
                    remaining = conn.execute(
                        'SELECT COUNT(*) FROM sequences WHERE id IN ('
                        + ','.join('?' for _ in rows) + ')',
                        tuple(row['sequence_id'] for row in rows),
                    ).fetchone()[0]
                if integrity != 'ok':
                    raise RuntimeError(f'staged Exit Interview database failed integrity_check: {integrity}')
                if remaining:
                    raise RuntimeError(f'{remaining} planned sequence(s) remained after Exit Interview')
                removed_count = len(rows)
                manifest.add_stage('exit_interview', 'COMPLETE', detail=f'{removed_count} removed')
                manifest.finish('COMPLETE')
                db.record_project_run(
                    f'exit-interview:{manifest.data["started_at"]}', 'exit_interview', 'COMPLETE',
                    manifest_path=str(manifest.json_path), started_at=manifest.data['started_at'],
                    completed_at=manifest.data['completed_at'],
                )
                publish_tmp = original.with_name(original.name + '.publishing')
                publish_tmp.unlink(missing_ok=True)
                with sqlite3.connect(str(staged_path)) as source, sqlite3.connect(str(publish_tmp)) as target:
                    source.backup(target)
                    if target.execute('PRAGMA integrity_check').fetchone()[0] != 'ok':
                        raise RuntimeError('published Exit Interview database failed integrity_check')
                os.replace(publish_tmp, original)
                staged_path.unlink(missing_ok=True)
        except BaseException as exc:
            if not blocked_rows:
                manifest.finish('FAILED', error=exc)
            raise
        finally:
            staged_path.unlink(missing_ok=True)
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)
    if blocked_rows:
        print(
            f'[exit-interview] No changes applied; {len(blocked_rows)} request(s) are blocked.\n'
            f'  Plan: {report}'
        )
        raise SystemExit(2)
    print(
        f'[exit-interview] Removed {removed_count} sequence(s) from active project state.\n'
        f'  Audit report: {report}\n'
        f'  Recovery FASTA: {output / "departing_sequences.fasta"}'
    )


def cmd_annual_report(args):
    from branchmanager.reporting import write_annual_report

    db = Database(args.db)
    db.initialise()
    manifest = RunManifest(args.out, 'annual_report')
    manifest.add_input(args.db, role='project_database')
    try:
        outputs = write_annual_report(db, args.out)
        for role, path in outputs.items():
            manifest.add_output(path, role=role)
        manifest.add_stage('annual_report', 'COMPLETE')
        manifest.finish('COMPLETE')
    except Exception as exc:
        manifest.finish('FAILED', error=exc)
        raise
    print(f"[annual-report] Cumulative project report: {outputs['html']}")


def cmd_it_desk(args):
    from branchmanager.it_desk import run_it_desk_checks, write_it_desk_report

    rows = run_it_desk_checks(
        db_path=args.db, references=args.references, output_dir=args.out,
        tree_method=args.tree_method,
    )
    outputs = write_it_desk_report(args.out, rows)
    for row in rows:
        print(f"[{row['status']}] {row['check']}: {row['detail']}")
    print(f"[it-desk] {outputs['status']}; report: {outputs['tsv']}")
    if args.strict and outputs['status'] != 'PASS':
        raise SystemExit(2)


def cmd_assistant(args):
    """Run the Assistant workflow with explicit scientific stage boundaries."""
    from branchmanager.it_desk import run_it_desk_checks, write_it_desk_report

    root = Path(args.out).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    manifest = RunManifest(root, 'assistant')
    for role, source in (
        ('sample_map', args.sample_map), ('marker_fasta', args.fasta),
        ('partner_metadata', args.partner_metadata), ('marker_qc', args.marker_qc),
        ('marker_review', args.marker_review),
        ('primary_reference', args.ref), ('reference_taxonomy', args.taxa),
        ('baseline_fasta', args.baseline_fasta), ('mwl', args.mwl),
        ('baseline_taxonomy', args.baseline_taxa_assignments),
        ('chimera_reference', args.chimera_ref),
    ):
        if source:
            manifest.add_input(source, role=role)
    for source in args.alt_ref or []:
        manifest.add_input(source, role='alternate_reference')
    for source in args.alt_taxa or []:
        manifest.add_input(source, role='alternate_reference_taxonomy')
    try:
        it_desk_dir = root / '00_it_desk'
        checks = run_it_desk_checks(
            db_path=args.db,
            references=[args.ref, *(args.alt_ref or [])],
            output_dir=root,
            tree_method=args.tree_method,
        )
        it_desk_report = write_it_desk_report(it_desk_dir, checks)
        manifest.add_output(it_desk_report['tsv'], role='it_desk_report')
        if it_desk_report['status'] != 'PASS':
            raise RuntimeError(f'IT Desk checks failed; inspect {it_desk_report["tsv"]}')

        onboarding_dir = root / '01_onboarding'
        onboarding_args = argparse.Namespace(
            sample_map=args.sample_map, fasta=args.fasta, partner_metadata=args.partner_metadata,
            read_dir=args.read_dir, primers=args.primers, partner_id=args.partner_id,
            dataset=args.dataset, out=str(onboarding_dir),
        )
        cmd_onboarding(onboarding_args)
        manifest.add_stage('onboarding', 'COMPLETE')

        paper_trail_dir = root / '02_paper_trail_merge_meeting'
        if args.sample_map:
            trace_args = argparse.Namespace(
                command='paper-trail', input=[], out=str(paper_trail_dir),
                sample_map=str(onboarding_dir / 'normalised_read_map.tsv'), primers=args.primers,
                primer_sequences=args.primer_sequences, trim_primers=args.trim_primers,
                min_quality=args.min_quality, min_length=args.min_length,
                min_read_length=args.min_read_length, min_mean_quality=args.min_mean_quality,
                mask_quality=args.mask_quality,
                max_read_expected_errors=args.max_read_expected_errors,
                max_output_expected_errors=args.max_output_expected_errors,
                warn_n_percent=args.warn_n_percent, max_n_percent=args.max_n_percent,
                warn_internal_low_quality_run=args.warn_internal_low_quality_run,
                max_internal_low_quality_run=args.max_internal_low_quality_run,
                secondary_peak_ratio=args.secondary_peak_ratio,
                max_mixed_peak_percent=args.max_mixed_peak_percent,
                mixed_peak_min_quality=args.mixed_peak_min_quality,
                max_conflict_density=args.max_conflict_density,
                max_ambiguous_overlap_conflicts_without_review=args.max_ambiguous_overlap_conflicts_without_review,
                quality_difference=args.quality_difference,
                allow_missing_quality=args.allow_missing_quality,
                min_overlap=args.min_overlap, min_overlap_identity=args.min_overlap_identity,
                screen_ref=args.ref, screen_taxa=args.taxa,
                threads=args.threads,
                max_report_image_height=args.max_report_image_height,
            )
            cmd_paper_trail(trace_args)
            marker_input = paper_trail_dir / 'assembled.fasta'
            marker_qc = paper_trail_dir / 'assembly_report.tsv'
            manifest.add_stage('paper_trail', 'COMPLETE')
            manifest.add_stage('merge_meeting', 'COMPLETE')
        else:
            if not args.marker_qc and not args.accept_unverified_marker_qc:
                raise RuntimeError(
                    'FASTA Assistant submissions require --marker-qc or the explicit '
                    '--accept-unverified-marker-qc audit flag'
                )
            marker_input = onboarding_dir / 'normalised_input.fasta'
            marker_qc = Path(args.marker_qc) if args.marker_qc else None
            manifest.add_stage('paper_trail', 'SKIPPED', detail='partner supplied FASTA')
            manifest.add_stage('merge_meeting', 'SKIPPED', detail='partner supplied FASTA')

        performance_review_dir = root / '03_performance_review_hiring_panel'
        performance_review_args = argparse.Namespace(
            command='performance-review', input=str(marker_input),
            partner_metadata=args.partner_metadata, db=args.db, dataset=args.dataset,
            out=str(performance_review_dir), ref=args.ref, taxa=args.taxa, ref_name=args.ref_name,
            alt_ref=args.alt_ref, alt_taxa=args.alt_taxa, alt_ref_name=args.alt_ref_name,
            main_ref=args.ref_name, baseline_fasta=args.baseline_fasta,
            baseline_dataset=args.baseline_dataset,
            baseline_tier=args.baseline_tier,
            baseline_taxa_assignments=args.baseline_taxa_assignments,
            baseline_skip_classify=False, baseline_colours=None, baseline_shorten_ids=False,
            novelty_baseline_datasets=[], mwl=args.mwl, mwl_sheet='MWL_V1', mwl_min_rank='p',
            taxa_assignments=None, previous_review=None, shorten_ids=False,
            min_len=800, max_n_percent=args.max_n_percent,
            sequence_domain=args.sequence_domain,
            phylum=None, target=None, force_rebuild=False, anchors=None,
            threads=args.threads, tree_method=args.tree_method,
            neighbourhood_format='png', pangenome_target=args.pangenome_target,
            candidate_set_size=args.candidate_set_size, user_colours=None,
            baseline_redundancy_identity=args.baseline_redundancy_identity,
            baseline_redundancy_min_query_coverage=args.baseline_redundancy_min_query_coverage,
            baseline_extension_min_identity=args.baseline_extension_min_identity,
            baseline_extension_min_query_coverage=args.baseline_extension_min_query_coverage,
            group_phyla=None, functional=None, draft_rumen_functions=False,
            marker_qc=str(marker_qc) if marker_qc else None,
            marker_review=args.marker_review,
            accept_unverified_marker_qc=bool(args.accept_unverified_marker_qc),
            chimera_ref=args.chimera_ref or args.ref,
            skip_chimera_check=args.skip_chimera_check, collapse=False,
        )
        cmd_performance_review(performance_review_args)
        manifest.add_stage('performance_review', 'COMPLETE')
        manifest.add_stage('hiring_panel', 'COMPLETE')
        for role, path in (
            ('marker_fasta', marker_input),
            ('performance_review_dashboard', performance_review_dir / 'performance_review_dashboard.html'),
            ('selection_summary', performance_review_dir / 'assessment' / 'selection_summary.tsv'),
            ('sequencing_sets', performance_review_dir / 'assessment' / 'sequencing_sets.tsv'),
        ):
            manifest.add_output(path, role=role)
        missing = manifest.verify_required_outputs()
        if missing:
            raise RuntimeError('Assistant workflow outputs missing: ' + ', '.join(missing))
        manifest.finish('COMPLETE')
    except BaseException as exc:
        manifest.finish('FAILED', error=exc)
        raise
    print(
        '[assistant] Onboarding through Hiring Panel complete.\n'
        f'  Decision dashboard: {performance_review_dir / "performance_review_dashboard.html"}\n'
        f'  Sequencing sets  : {performance_review_dir / "assessment" / "sequencing_sets.tsv"}'
    )


def _cmd_quarterly_review_impl(args):
    """Run a project-wide post-coverage genome-selection round."""
    from datetime import datetime, timezone

    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    _configure_logging(str(outdir))
    log = logging.getLogger(__name__)
    db = Database(args.db)
    db.initialise()

    # A final cumulative metadata file may be supplied without re-running the
    # expensive classification/tree workflow. It updates factual genome status
    # only; Quarterly Review recommendations never write back to that status.
    if getattr(args, 'partner_metadata', None):
        with db.connect() as conn:
            sequence_rows = conn.execute('SELECT id, dataset FROM sequences').fetchall()
        roles = db.get_dataset_roles()
        candidate_ids = [
            str(sequence_id) for sequence_id, dataset in sequence_rows
            if roles.get(str(dataset or ''), {}).get('role') != 'baseline'
        ]
        setattr(args, 'dataset', 'QuarterlyReview')
        _load_partner_metadata_for_run(args, db, str(outdir), {}, candidate_ids)

    round_id = getattr(args, 'round_id', None) or datetime.now(timezone.utc).strftime(
        'quarterly_review_%Y%m%dT%H%M%SZ'
    )
    assessment_paths = list(getattr(args, 'assessment', None) or [])
    if assessment_paths:
        try:
            imported_rows = quarterly_review.read_assessment_tables(assessment_paths)
        except Exception as exc:
            raise SystemExit(f'[QUARTERLY REVIEW] Failed to import assessment table: {exc}')
        if not imported_rows:
            raise SystemExit('[QUARTERLY REVIEW] The supplied assessment table(s) contained no sequence rows.')
        snapshot_id = f'quarterly-review-import:{round_id}'
        db.save_assessment_snapshot(
            snapshot_id,
            imported_rows,
            dataset='quarterly_review_import',
            source_path=';'.join(str(Path(path).resolve()) for path in assessment_paths),
        )
        log.info('[QUARTERLY REVIEW] Imported %d latest assessment rows into %s', len(imported_rows), snapshot_id)

    latest = db.get_latest_assessment_rows()
    if not latest:
        raise SystemExit(
            '[QUARTERLY REVIEW] No assessment snapshots are stored. Run Performance Review once with the current '
            'BranchManager version or provide --assessment path/to/sequence_assessment.tsv.'
        )

    tree_path = getattr(args, 'tree', None)
    if not tree_path and getattr(args, 'from_dir', None):
        tree_path = _find_tree_file_in_dir(args.from_dir)
    if tree_path and not Path(tree_path).exists():
        raise SystemExit(f'[QUARTERLY REVIEW] Tree file not found: {tree_path}')
    if not tree_path:
        log.warning(
            '[QUARTERLY REVIEW] No cumulative tree supplied; ranking will use marker identities and other '
            'assessment evidence without marginal patristic distance.'
        )

    alignment_path = getattr(args, 'alignment', None)
    if not alignment_path and getattr(args, 'from_dir', None):
        for candidate in (
            Path(args.from_dir) / 'tree' / 'current_alignment.fasta',
            Path(args.from_dir) / 'current_alignment.fasta',
        ):
            if candidate.exists():
                alignment_path = str(candidate)
                break
    if not alignment_path and tree_path:
        candidate = Path(tree_path).parent / 'current_alignment.fasta'
        if candidate.exists():
            alignment_path = str(candidate)
    if alignment_path and not Path(alignment_path).exists():
        raise SystemExit(f'[QUARTERLY REVIEW] Alignment file not found: {alignment_path}')

    parameters = {
        'GenomeBudget': int(args.genome_budget),
        'BackupsPerPrimary': int(args.backups_per_primary),
        'PangenomeTarget': int(args.pangenome_target),
        'BaselineRedundancyIdentity': float(args.baseline_redundancy_identity),
        'BaselineRedundancyMinQueryCoverage': float(args.baseline_redundancy_min_query_coverage),
        'BaselineExtensionMinIdentity': float(args.baseline_extension_min_identity),
        'BaselineExtensionMinQueryCoverage': float(args.baseline_extension_min_query_coverage),
        'IncludeModerateEvidence': bool(args.include_moderate_evidence),
        'Tree': str(Path(tree_path).resolve()) if tree_path else 'None',
        'Alignment': str(Path(alignment_path).resolve()) if alignment_path else 'None',
        'AssessmentRows': len(latest),
    }
    recommendations = quarterly_review.build_quarterly_review(
        list(latest.values()),
        db=db,
        tree_path=tree_path,
        alignment_path=alignment_path,
        genome_budget=args.genome_budget,
        backups_per_primary=args.backups_per_primary,
        pangenome_target=args.pangenome_target,
        baseline_redundancy_identity=args.baseline_redundancy_identity,
        baseline_redundancy_min_query_coverage=args.baseline_redundancy_min_query_coverage,
        baseline_extension_min_identity=args.baseline_extension_min_identity,
        baseline_extension_min_query_coverage=args.baseline_extension_min_query_coverage,
        include_moderate=args.include_moderate_evidence,
    )

    # Quarterly Review recommendations need round-specific figures. Reusing the most
    # recent Performance Review image would display stale PRIMARY/BACKUP stars.
    if tree_path:
        try:
            recommendation_by_id = {row['sequence_id']: row for row in recommendations}
            visual_rows = [quarterly_review.normalise_assessment_row(row) for row in latest.values()]
            for row in visual_rows:
                recommendation = recommendation_by_id.get(row['id'], {})
                role = str(recommendation.get('role') or '')
                row['sequencing_set_role'] = {
                    'PRIMARY': 'PRIMARY',
                    'BACKUP': 'BACKUP',
                    'ALREADY_SELECTED': 'COMMITTED',
                    'ALREADY_SEQUENCED': 'SEQUENCED',
                }.get(role, '')
                row['sequencing_set_rank'] = (
                    recommendation.get('round_rank', 'NA') if role in {'PRIMARY', 'BACKUP'} else 'NA'
                )
                if role == 'PRIMARY':
                    row['sequencing_set_badge'] = f'P{recommendation.get("round_rank", "")}'
                elif role == 'BACKUP':
                    row['sequencing_set_badge'] = (
                        f'B{recommendation.get("round_rank", "")}.{recommendation.get("backup_rank", "")}'
                    )
                row['selected_for_genome_sequencing'] = str(role == 'ALREADY_SELECTED')
                row['already_sequenced'] = str(role == 'ALREADY_SEQUENCED')
            visual_result = neighbourhood.generate_local_neighbourhood_visuals(
                tree_path=tree_path,
                assessment_rows=visual_rows,
                db=db,
                outdir=outdir / 'neighbourhoods',
                alignment_path=alignment_path,
                image_format=args.neighbourhood_format,
            )
            for row in visual_rows:
                recommendation = recommendation_by_id.get(row['id'])
                if recommendation is not None:
                    recommendation['local_tree_figure'] = row.get('local_neighbourhood_figure', 'NA')
            parameters['NeighbourhoodFigures'] = len(visual_result['figures'])
            parameters['NeighbourhoodFormat'] = args.neighbourhood_format
        except Exception as exc:
            log.warning('[QUARTERLY REVIEW] Could not generate round-specific neighbourhood figures: %s', exc)
    outputs = quarterly_review.write_quarterly_review_reports(outdir, round_id, recommendations, parameters)
    db.save_selection_round(round_id, recommendations, mode='quarterly_review', parameters=parameters)

    primary_count = sum(row.get('role') == 'PRIMARY' for row in recommendations)
    backup_count = sum(row.get('role') == 'BACKUP' for row in recommendations)
    log.info(
        '[QUARTERLY REVIEW] Round %s nominated %d primary and %d backup isolate(s)',
        round_id, primary_count, backup_count,
    )
    print(
        '[quarterly-review] Done.\n'
        f'  Round             : {round_id}\n'
        f'  Primary nominees  : {primary_count}\n'
        f'  Backup nominees   : {backup_count}\n'
        f'  Next genome set   : {outputs["selected"]}\n'
        f'  Full audit        : {outputs["summary"]}\n'
        f'  Round manifest    : {outputs["manifest"]}'
    )


def cmd_quarterly_review(args):
    import fcntl
    import sqlite3

    manifest = RunManifest(args.out, 'quarterly_review')
    original = Path(args.db).expanduser().resolve()
    if not original.is_file():
        raise SystemExit(f'[QUARTERLY REVIEW] Project database does not exist: {original}')
    manifest.add_input(original, role='project_database')
    for source in getattr(args, 'assessment', None) or []:
        manifest.add_input(source, role='assessment_import')
    for role, source in (
        ('partner_metadata', getattr(args, 'partner_metadata', None)),
        ('tree', getattr(args, 'tree', None)), ('alignment', getattr(args, 'alignment', None)),
    ):
        if source:
            manifest.add_input(source, role=role)
    staged_path = Path(args.out).resolve() / '.quarterly_review_project.sqlite'
    lock_path = original.with_name(original.name + '.lock')
    with open(lock_path, 'a+') as lock_handle:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        try:
            staged_path.unlink(missing_ok=True)
            with sqlite3.connect(str(original)) as source, sqlite3.connect(str(staged_path)) as target:
                source.backup(target)
            staged_args = argparse.Namespace(**vars(args))
            staged_args.db = str(staged_path)
            manifest.add_stage('quarterly_review', 'RUNNING')
            _cmd_quarterly_review_impl(staged_args)
            for role, name in (
                ('quarterly_review_summary', 'quarterly_review_summary.tsv'),
                ('next_genome_set', 'next_genome_set.tsv'),
                ('quarterly_review_parameters', 'quarterly_review_manifest.tsv'),
            ):
                manifest.add_output(Path(args.out) / name, role=role)
            missing = manifest.verify_required_outputs()
            if missing:
                raise RuntimeError('required Quarterly Review outputs are missing: ' + ', '.join(missing))
            manifest.add_stage('quarterly_review', 'COMPLETE')
            manifest.finish('COMPLETE')
            db = Database(str(staged_path))
            db.initialise()
            db.record_project_run(
                str(getattr(args, 'round_id', None) or f'quarterly-review:{manifest.data["started_at"]}'),
                'quarterly_review', 'COMPLETE', manifest_path=str(manifest.json_path),
                started_at=manifest.data['started_at'], completed_at=manifest.data['completed_at'],
            )
            publish_tmp = original.with_name(original.name + '.publishing')
            publish_tmp.unlink(missing_ok=True)
            with sqlite3.connect(str(staged_path)) as source, sqlite3.connect(str(publish_tmp)) as target:
                source.backup(target)
                if target.execute('PRAGMA integrity_check').fetchone()[0] != 'ok':
                    raise RuntimeError('staged Quarterly Review database failed integrity_check')
            os.replace(publish_tmp, original)
            staged_path.unlink(missing_ok=True)
        except BaseException as exc:
            manifest.finish('FAILED', error=exc)
            raise
        finally:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == 'filing-cabinet':
        cmd_filing_cabinet(args)
    elif args.command == 'performance-review':
        cmd_performance_review(args)
    elif args.command == 'quarterly-review':
        cmd_quarterly_review(args)
    elif args.command == 'label-maker':
        cmd_label_maker(args)
    elif args.command == 'org-chart':
        cmd_org_chart(args)
    elif args.command == 'background-check':
        cmd_background_check(args)
    elif args.command == 'paper-trail':
        cmd_paper_trail(args)
    elif args.command == 'mailroom':
        cmd_mailroom(args)
    elif args.command == 'interview':
        cmd_interview(args)
    elif args.command == 'onboarding':
        cmd_onboarding(args)
    elif args.command == 'status-meeting':
        _cmd_project_import(args, genomes=False)
    elif args.command == 'records-update':
        _cmd_project_import(args, genomes=True)
    elif args.command == 'exit-interview':
        cmd_exit_interview(args)
    elif args.command == 'annual-report':
        cmd_annual_report(args)
    elif args.command == 'it-desk':
        cmd_it_desk(args)
    elif args.command == 'assistant':
        cmd_assistant(args)
    else:
        parser.print_help()


if __name__ == '__main__':
    main()
