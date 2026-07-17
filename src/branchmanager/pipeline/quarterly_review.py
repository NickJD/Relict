"""Project-wide post-coverage Quarterly Review genome-selection rounds."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional

from branchmanager.pipeline.neighbourhood import _load_alignment_sequences, _msa_pident, pairwise_leaf_distances
from branchmanager.pipeline.selection_sets import (
    DEFAULT_BASELINE_EXTENSION_MIN_IDENTITY,
    DEFAULT_BASELINE_EXTENSION_MIN_QUERY_COVERAGE,
    DEFAULT_BASELINE_REDUNDANCY_IDENTITY,
    DEFAULT_BASELINE_REDUNDANCY_QUERY_COVERAGE,
    DEFAULT_PANGENOME_TARGET,
    baseline_extension_evidence,
    baseline_redundancy_evidence,
)
from branchmanager.taxonomy import normalise_taxon_name, parse_taxon_string


MWL_STRENGTH = {
    'species': 5, 's': 5, 'genus': 4, 'g': 4, 'family': 3, 'f': 3,
    'order': 2, 'o': 2, 'class': 1, 'c': 1, 'phylum': 0, 'p': 0,
}


def _value(row: dict, *names: str, default='NA'):
    for name in names:
        value = row.get(name)
        if value not in (None, ''):
            return value
    return default


def _float(value, default: Optional[float] = None) -> Optional[float]:
    try:
        if value in (None, '', 'NA', 'None'):
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def _int(value, default: int = 0) -> int:
    try:
        if value in (None, '', 'NA', 'None'):
            return default
        return int(float(value))
    except (TypeError, ValueError):
        return default


def _true(value) -> bool:
    return str(value or '').strip().lower() in {'1', 'true', 'yes', 'y'}


def _primary_taxonomy(row: dict) -> str:
    direct = _value(row, 'taxonomy', 'GTDBTaxonomy', 'Taxonomy', default='')
    if direct:
        return str(direct)
    excluded = ('Baseline', 'Project', 'Reference', 'NearestHit')
    for key, value in row.items():
        if key.endswith('Taxonomy') and not any(token in key for token in excluded) and value:
            return str(value)
    return ''


def normalise_assessment_row(row: dict, *, source_path: str = '') -> dict:
    """Normalise a sequence_assessment.tsv row or an in-memory assessment row."""
    reference_identity = _value(
        row, 'reference_nearest_identity', 'GTDBReferenceNearestIdentity',
        'ReferenceNearestIdentity', default='NA',
    )
    normalised = {
        'id': str(_value(row, 'id', 'ID', 'SequenceID', default='')).strip(),
        'partner_id': _value(row, 'partner_id', 'PartnerID'),
        'taxonomy': _primary_taxonomy(row),
        'classification_identity': _value(
            row, 'classification_identity', 'GTDBClassificationIdentity', 'ClassificationIdentity',
        ),
        'classification_confidence': _value(
            row, 'classification_confidence', 'GTDBClassificationConfidence', 'ClassificationConfidence',
        ),
        'classification_query_coverage': _value(
            row, 'classification_query_coverage', 'GTDBQueryCoverage', 'QueryCoverage',
        ),
        'baseline_evidence_class': _value(row, 'baseline_evidence_class', 'BaselineEvidenceClass'),
        'hungate_nearest_hit': _value(row, 'hungate_nearest_hit', 'HungateNearestHit'),
        'hungate_nearest_hit_taxonomy': _value(row, 'hungate_nearest_hit_taxonomy', 'HungateNearestHitTaxonomy'),
        'hungate_nearest_identity': _value(row, 'hungate_nearest_identity', 'HungateNearestIdentity'),
        'hungate_nearest_query_coverage': _value(row, 'hungate_nearest_query_coverage', 'HungateNearestQueryCoverage'),
        'secondary_baseline_nearest_hit': _value(row, 'secondary_baseline_nearest_hit', 'SecondaryBaselineNearestHit'),
        'secondary_baseline_nearest_hit_taxonomy': _value(row, 'secondary_baseline_nearest_hit_taxonomy', 'SecondaryBaselineNearestHitTaxonomy'),
        'secondary_baseline_nearest_identity': _value(row, 'secondary_baseline_nearest_identity', 'SecondaryBaselineNearestIdentity'),
        'secondary_baseline_nearest_query_coverage': _value(row, 'secondary_baseline_nearest_query_coverage', 'SecondaryBaselineNearestQueryCoverage'),
        'cultured_rumen_nearest_hit': _value(row, 'cultured_rumen_nearest_hit', 'CulturedRumenNearestHit', 'nearest_hit', 'BaselineNearestHit'),
        'cultured_rumen_nearest_hit_taxonomy': _value(row, 'cultured_rumen_nearest_hit_taxonomy', 'CulturedRumenNearestHitTaxonomy', 'nearest_hit_taxonomy', 'BaselineNearestHitTaxonomy'),
        'cultured_rumen_nearest_identity': _value(row, 'cultured_rumen_nearest_identity', 'CulturedRumenNearestIdentity', 'nearest_identity', 'BaselineNearestIdentity'),
        'cultured_rumen_nearest_query_coverage': _value(row, 'cultured_rumen_nearest_query_coverage', 'CulturedRumenNearestQueryCoverage', 'nearest_query_coverage', 'BaselineNearestQueryCoverage'),
        'nearest_hit': _value(row, 'nearest_hit', 'CulturedRumenNearestHit', 'BaselineNearestHit'),
        'nearest_hit_taxonomy': _value(
            row, 'nearest_hit_taxonomy', 'CulturedRumenNearestHitTaxonomy', 'BaselineNearestHitTaxonomy',
        ),
        'nearest_identity': _value(row, 'nearest_identity', 'CulturedRumenNearestIdentity', 'BaselineNearestIdentity'),
        'nearest_query_coverage': _value(
            row, 'nearest_query_coverage', 'CulturedRumenNearestQueryCoverage', 'BaselineNearestQueryCoverage',
        ),
        'density_source': _value(row, 'density_source', 'CulturedRumenSource', 'BaselineSource', default=''),
        'project_nearest_identity': _value(row, 'project_nearest_identity', 'ProjectNearestIdentity'),
        'project_density_source': _value(row, 'project_density_source', 'ProjectSource', default=''),
        'reference_nearest_identity': reference_identity,
        'reference_density_source': _value(row, 'reference_density_source', 'GTDBReferenceSource', default=''),
        'selected_for_genome_sequencing': _value(
            row, 'selected_for_genome_sequencing', 'SelectedForGenomeSequencing', default='False',
        ),
        'already_sequenced': _value(
            row, 'already_sequenced', 'GenomeAlreadySequenced', default='False',
        ),
        'nearest_genome_hit': _value(row, 'nearest_genome_hit', 'NearestGenomeHit'),
        'nearest_genome_identity': _value(row, 'nearest_genome_identity', 'NearestGenomeIdentity'),
        'nearest_genome_identity_source': _value(
            row, 'nearest_genome_identity_source', 'NearestGenomeIdentitySource',
            default='assessment_snapshot',
        ),
        'hungate_genome_count_same_species': _value(
            row, 'hungate_genome_count_same_species', 'HungateGenomesSameAssessmentSpecies',
            'HungateGenomeCountSameAssessmentSpecies', default='0',
        ),
        'secondary_baseline_genome_count_same_species': _value(
            row, 'secondary_baseline_genome_count_same_species', 'SecondaryBaselineGenomesSameAssessmentSpecies',
            'SecondaryBaselineGenomeCountSameAssessmentSpecies', default='0',
        ),
        'cultured_rumen_genome_count_same_species': _value(
            row, 'cultured_rumen_genome_count_same_species', 'CulturedRumenGenomesSameAssessmentSpecies',
            'CulturedRumenGenomeCountSameAssessmentSpecies', 'BaselineGenomesSameAssessmentSpecies',
            'BaselineGenomesSameSpecies', default='0',
        ),
        'genome_available_count_same_species': _value(
            row, 'genome_available_count_same_species', 'CulturedRumenGenomesSameAssessmentSpecies',
            'CulturedRumenGenomeCountSameAssessmentSpecies', 'BaselineGenomesSameAssessmentSpecies',
            'BaselineGenomesSameSpecies', default='0',
        ),
        'genome_selected_count_same_species': _value(
            row, 'genome_selected_count_same_species', 'SequencedPartnerGenomesSameAssessmentSpecies',
            'SelectedPartnerGenomesSameAssessmentSpecies', 'SequencedPartnerGenomesSameSpecies',
            'SelectedPartnerGenomesSameSpecies', default='0',
        ),
        'genome_pending_count_same_species': _value(
            row, 'genome_pending_count_same_species',
            'SelectedPendingGenomesSameAssessmentSpecies', default='0',
        ),
        'genome_committed_count_same_species': _value(
            row, 'genome_committed_count_same_species', 'CommittedGenomesSameAssessmentSpecies',
            'AvailableGenomesSameAssessmentSpecies', 'CommittedGenomesSameSpecies', 'AvailableGenomesSameSpecies',
            default='0',
        ),
        'pangenome_target': _value(row, 'pangenome_target', 'PangenomeTarget', default='9'),
        'pangenome_gap': _value(row, 'pangenome_gap', 'PangenomeGap', default='0'),
        'mwl_matched_rank': _value(row, 'mwl_matched_rank', 'MWLMatchedRank'),
        'mwl_matched_taxon': _value(row, 'mwl_matched_taxon', 'MWLMatchedTaxon'),
        'mwl_score': _value(row, 'mwl_score', 'MWLScore', default='0'),
        'phylo_isolation': _value(row, 'phylo_isolation', 'PhyloIsolation'),
        'local_neighbourhood_figure': _value(
            row, 'local_neighbourhood_figure', 'LocalNeighbourhoodFigure', 'LocalTreeFigure',
        ),
        'cluster_representative': _value(row, 'cluster_representative', 'ClusterRepresentative', default='self'),
        'in_tree': _value(row, 'in_tree', 'InTree', default='Unknown'),
        'placement_flags': _value(row, 'placement_flags', 'PlacementFlags', default=''),
        'marker_qc_class': _value(row, 'marker_qc_class', 'MarkerQCClass', 'MarkerQC', default='QUALITY_UNVERIFIED'),
        'marker_manual_review_status': _value(row, 'marker_manual_review_status', 'MarkerManualReviewStatus', 'MarkerReview', default='NOT_REVIEWED'),
        'marker_qc_flag': _value(row, 'marker_qc_flag', 'MarkerQCFlag', default=''),
        'evidence_quality': _value(row, 'evidence_quality', 'EvidenceQuality', default=''),
        '_snapshot_id': _value(row, '_snapshot_id', default=''),
        '_snapshot_source_path': _value(row, '_snapshot_source_path', default=source_path),
    }
    return normalised


def read_assessment_tables(paths: Iterable[str | Path]) -> List[dict]:
    """Read one or more full assessment tables; later files replace earlier rows."""
    latest = {}
    for path in paths:
        source = Path(path)
        with open(source, newline='') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            for raw in reader:
                row = normalise_assessment_row(raw, source_path=str(source.resolve()))
                if row['id']:
                    latest[row['id']] = row
    return list(latest.values())


def _species(taxonomy: object) -> str:
    value = str(parse_taxon_string(str(taxonomy or '')).get('s') or '').strip()
    if not value or value.lower() in {'na', 'none', 'unclassified'}:
        return ''
    return value


def _group(row: dict) -> tuple[str, str]:
    species = _species(row.get('taxonomy'))
    if species:
        return f'species:{normalise_taxon_name(species)}', species
    parsed = parse_taxon_string(str(row.get('taxonomy') or ''))
    for rank in ('g', 'f', 'o', 'c', 'p', 'd'):
        taxon = str(parsed.get(rank) or '').strip()
        if taxon:
            figure = str(row.get('local_neighbourhood_figure') or 'no_tree_context')
            return f'local:{figure}|{rank}:{normalise_taxon_name(taxon)}', taxon
    return f'local:{row.get("local_neighbourhood_figure") or "unclassified"}', 'unclassified'


def _tree_leaf(row: dict) -> str:
    sequence_id = str(row.get('id') or '')
    representative = str(row.get('cluster_representative') or '')
    if str(row.get('in_tree') or '').lower() == 'yes':
        return sequence_id
    if representative not in {'', 'self', 'duplicate', 'N/A', 'NA'}:
        return representative
    return sequence_id


def _evidence_quality(row: dict) -> str:
    explicit = str(row.get('evidence_quality') or '').strip().upper()
    if explicit in {'HIGH', 'MODERATE', 'LOW'}:
        return explicit
    flag_set = {f.strip() for f in str(row.get('placement_flags') or '').upper().split(';') if f.strip()}
    identity = _float(row.get('classification_identity'))
    confidence = _float(row.get('classification_confidence'))
    coverage = _float(row.get('classification_query_coverage'))
    taxonomy = str(row.get('taxonomy') or '').strip().lower()
    # Use exact flag membership to avoid CHIMERA_INDETERMINATE matching as CHIMERA.
    if flag_set & {
        'NO_REFERENCE_HIT', 'NO_CLASSIFICATION', 'CHIMERA', 'CHIMERA_CONFIRMED',
        'VERY_SHORT', 'HIGH_N_CONTENT', 'MARKER_QC_FAILED',
        'MARKER_QC_REVIEW_REQUIRED', 'MARKER_QC_UNVERIFIED',
    }:
        return 'LOW'
    # MARKER_QC_REVIEW_APPROVED: fall through to standard checks.
    if identity is not None and identity < 90.0:
        return 'LOW'
    if coverage is not None and coverage < 80.0:
        return 'LOW'
    if identity is None and taxonomy in {'', 'na', 'none'}:
        return 'LOW'
    if flag_set & {'LOW_CLASSIFICATION', 'LOW_CONFIDENCE', 'DISAGREEMENT', 'CONFLICT'}:
        return 'MODERATE'
    if (identity is not None and identity < 95.0) or (confidence is not None and confidence < 0.8) or (coverage is not None and coverage < 90.0):
        return 'MODERATE'
    return 'HIGH'


def _available_ids_by_species(db) -> Dict[str, List[str]]:
    if db is None:
        return {}
    roles = db.get_dataset_roles()
    with db.connect() as conn:
        records = conn.execute(
            "SELECT s.id, s.dataset, t.taxonomy, COALESCE(m.selected_for_wgs, 0) "
            "FROM sequences s LEFT JOIN taxonomy t ON s.id = t.id "
            "LEFT JOIN sequencing_metadata m ON s.id = m.id ORDER BY t.rowid"
        ).fetchall()
    available = defaultdict(set)
    for sequence_id, dataset, taxonomy, already_sequenced in records:
        species = _species(taxonomy)
        dataset_role = roles.get(str(dataset or ''), {})
        is_baseline = dataset_role.get('role') == 'baseline' and dataset_role.get('genomes_available')
        if species and (is_baseline or bool(already_sequenced)):
            available[normalise_taxon_name(species)].add(str(sequence_id))
    return {key: sorted(ids) for key, ids in available.items()}


def _genome_feedback_by_species(db) -> dict:
    """Summarise QC-passed genome taxonomy and ANI clusters by genome species."""
    if db is None:
        return {}
    feedback = defaultdict(lambda: {'genome_ids': set(), 'ani_clusters': set()})
    with db.connect() as conn:
        rows = conn.execute(
            'SELECT g.genome_id, g.sequence_id, g.gtdb_taxonomy, g.ani_cluster, t.taxonomy '
            'FROM genome_records g LEFT JOIN taxonomy t ON t.id=g.sequence_id '
            'WHERE g.genome_qc_pass=1 ORDER BY g.genome_id, t.rowid'
        ).fetchall()
    for genome_id, sequence_id, genome_taxonomy, ani_cluster, marker_taxonomy in rows:
        species = _species(genome_taxonomy or marker_taxonomy)
        if not species:
            continue
        key = normalise_taxon_name(species)
        feedback[key]['genome_ids'].add(str(sequence_id or genome_id))
        if ani_cluster:
            feedback[key]['ani_clusters'].add(str(ani_cluster))
    return {
        key: {
            'genome_ids': sorted(value['genome_ids']),
            'ani_clusters': sorted(value['ani_clusters']),
        }
        for key, value in feedback.items()
    }


def _distance(left: str, right: str, distances: Dict[tuple[str, str], float]) -> Optional[float]:
    return distances.get(tuple(sorted((left, right))))


def _marginal_distance(row: dict, anchors: List[str], distances: Dict[tuple[str, str], float]) -> Optional[float]:
    leaf = _tree_leaf(row)
    values = [
        value for anchor in anchors
        if anchor != leaf and (value := _distance(leaf, anchor, distances)) is not None
    ]
    return min(values) if values else None


def _refresh_nearest_available_from_alignment(
    grouped: Dict[str, List[dict]],
    anchors_by_group: Dict[str, List[str]],
    alignment_path: Optional[str | Path],
) -> int:
    """Refresh nearest genome context after genome-availability metadata changes."""
    aligned = _load_alignment_sequences(alignment_path)
    if not aligned:
        return 0
    refreshed = 0
    for key, members in grouped.items():
        anchors = anchors_by_group.get(key, [])
        for row in members:
            query_id = _tree_leaf(row)
            query = aligned.get(query_id)
            if not query:
                continue
            hits = []
            for anchor in anchors:
                if anchor == query_id or anchor not in aligned:
                    continue
                identity, compared = _msa_pident(query, aligned[anchor])
                if identity is not None and compared > 0:
                    hits.append((identity, compared, anchor))
            if not hits:
                continue
            identity, compared, anchor = max(hits, key=lambda item: (item[0], item[1], item[2]))
            row['nearest_genome_hit'] = anchor
            row['nearest_genome_identity'] = f'{identity:.2f}'
            row['nearest_genome_identity_source'] = f'current_alignment:{compared}_comparable_acgt_columns'
            refreshed += 1
    return refreshed


def _divergence(identity) -> float:
    value = _float(identity)
    return max(0.0, 100.0 - value) if value is not None and value > 0 else 0.0


def _baseline_evidence_rank(row: dict) -> int:
    ranks = {
        'NOVEL TO BOTH BASELINES': 4,
        'HUNGATE GAP / SECONDARY COVERED': 3,
        'HUNGATE COVERED': 2,
        'NO CULTURED BASELINE AVAILABLE': 1,
        'CULTURED BASELINE REDUNDANT': 0,
    }
    value = str(row.get('baseline_evidence_class') or '').strip().upper()
    return ranks.get(value, 1)


def _candidate_key(row: dict, anchors: List[str], distances: Dict[tuple[str, str], float]) -> tuple:
    patristic = _marginal_distance(row, anchors, distances)
    mwl = MWL_STRENGTH.get(str(row.get('mwl_matched_rank') or '').lower(), 0)
    mwl_score = _float(row.get('mwl_score'), 0.0) or 0.0
    return (
        0 if patristic is not None else 1,
        -(patristic or 0.0),
        -_baseline_evidence_rank(row),
        -_divergence(row.get('nearest_genome_identity')),
        -_divergence(row.get('nearest_identity')),
        -_divergence(row.get('project_nearest_identity')),
        -_divergence(row.get('reference_nearest_identity')),
        -(_float(row.get('phylo_isolation'), 0.0) or 0.0),
        -mwl,
        -mwl_score,
        str(row.get('id') or ''),
    )


def _tier(row: dict, gap: int) -> str:
    if gap > 0:
        return 'COVERAGE_GAP'
    nearest = _float(row.get('nearest_genome_identity'))
    if nearest is not None and 0 < nearest < 98.65:
        return 'NOVEL_GENOME_NEIGHBOURHOOD'
    if nearest is not None and nearest >= 99.0:
        return 'REDUNDANT_EXTENSION'
    return 'DIVERSITY_EXPANSION'


def _reason(row: dict, tier: str, available_count: int, target: int, patristic: Optional[float]) -> str:
    parts = []
    if tier == 'COVERAGE_GAP':
        parts.append(f'fills remaining {target}-genome coverage gap')
    elif tier == 'NOVEL_GENOME_NEIGHBOURHOOD':
        parts.append('no close genome representation at the 98.65% marker threshold')
    elif tier == 'DIVERSITY_EXPANSION':
        parts.append(f'post-target diversity expansion from {available_count} available genomes')
    else:
        parts.append('post-target extension despite close existing genome representation')
    nearest = _float(row.get('nearest_genome_identity'))
    if nearest is not None and nearest > 0:
        parts.append(f'nearest available genome {nearest:.2f}% marker identity')
    if patristic is not None:
        parts.append(f'marginal patristic distance {patristic:.5f}')
    baseline = _float(row.get('nearest_identity'))
    if baseline is not None and baseline > 0:
        parts.append(f'nearest cultured baseline {baseline:.2f}%')
    mwl_rank = str(row.get('mwl_matched_rank') or '')
    if mwl_rank not in {'', 'NA', 'None'}:
        parts.append(f'MWL {mwl_rank}-level context')
    return '; '.join(parts)


def build_quarterly_review(
    assessment_rows: List[dict],
    *,
    db=None,
    tree_path: Optional[str | Path] = None,
    alignment_path: Optional[str | Path] = None,
    genome_budget: int,
    backups_per_primary: int = 1,
    pangenome_target: int = DEFAULT_PANGENOME_TARGET,
    baseline_redundancy_identity: float = DEFAULT_BASELINE_REDUNDANCY_IDENTITY,
    baseline_redundancy_min_query_coverage: float = DEFAULT_BASELINE_REDUNDANCY_QUERY_COVERAGE,
    baseline_extension_min_identity: float = DEFAULT_BASELINE_EXTENSION_MIN_IDENTITY,
    baseline_extension_min_query_coverage: float = DEFAULT_BASELINE_EXTENSION_MIN_QUERY_COVERAGE,
    include_moderate: bool = False,
) -> List[dict]:
    """Select a globally balanced next sequencing tranche without changing status."""
    budget = max(1, int(genome_budget))
    backup_count = max(0, int(backups_per_primary))
    target = max(1, int(pangenome_target))
    redundancy_identity = float(baseline_redundancy_identity)
    redundancy_coverage = float(baseline_redundancy_min_query_coverage)
    extension_identity = float(baseline_extension_min_identity)
    extension_coverage = float(baseline_extension_min_query_coverage)
    if not 0.0 <= redundancy_identity <= 100.0:
        raise ValueError('baseline redundancy identity must be between 0 and 100')
    if not 0.0 <= redundancy_coverage <= 100.0:
        raise ValueError('baseline redundancy minimum query coverage must be between 0 and 100')
    if not 0.0 <= extension_identity <= redundancy_identity:
        raise ValueError('baseline extension minimum identity must be between 0 and the redundancy identity threshold')
    if not 0.0 <= extension_coverage <= 100.0:
        raise ValueError('baseline extension minimum query coverage must be between 0 and 100')
    rows = [normalise_assessment_row(row) for row in assessment_rows]
    rows = [row for row in rows if row['id']]

    metadata = db.get_sequencing_metadata_for_ids([row['id'] for row in rows]) if db else {}
    for row in rows:
        if row['id'] in metadata:
            row['selected_for_genome_sequencing'] = str(bool(metadata[row['id']]['selected_for_sequencing']))
            row['already_sequenced'] = str(bool(metadata[row['id']]['selected_for_wgs']))

    grouped = defaultdict(list)
    group_names = {}
    for row in rows:
        key, name = _group(row)
        row['_quarterly_review_group'] = key
        grouped[key].append(row)
        group_names[key] = name

    available_species = _available_ids_by_species(db)
    genome_feedback = _genome_feedback_by_species(db)
    anchors_by_group = {}
    available_anchors_by_group = {}
    counts_by_group = {}
    baseline_counts_by_group = {}
    hungate_counts_by_group = {}
    secondary_counts_by_group = {}
    for key, members in grouped.items():
        species = _species(members[0].get('taxonomy'))
        available_anchors = list(available_species.get(normalise_taxon_name(species), [])) if species else []
        available_anchors.extend(row['id'] for row in members if _true(row.get('already_sequenced')))
        available_anchors_by_group[key] = sorted(set(available_anchors))
        committed_anchors = list(available_anchors)
        committed_anchors.extend(row['id'] for row in members if _true(row.get('selected_for_genome_sequencing')))
        anchors_by_group[key] = sorted(set(committed_anchors))
        reported = max((_int(row.get('genome_committed_count_same_species')) for row in members), default=0)
        counts_by_group[key] = max(reported, len(anchors_by_group[key]))
        hungate_count = max((_int(row.get('hungate_genome_count_same_species')) for row in members), default=0)
        secondary_count = max((_int(row.get('secondary_baseline_genome_count_same_species')) for row in members), default=0)
        prior_baseline_count = max((_int(row.get('genome_available_count_same_species')) for row in members), default=0)
        if hungate_count + secondary_count == 0 and prior_baseline_count:
            hungate_count = prior_baseline_count
        hungate_counts_by_group[key] = hungate_count
        secondary_counts_by_group[key] = secondary_count
        baseline_counts_by_group[key] = hungate_count + secondary_count

    _refresh_nearest_available_from_alignment(grouped, available_anchors_by_group, alignment_path)

    leaves = {_tree_leaf(row) for row in rows}
    for anchors in anchors_by_group.values():
        leaves.update(anchors)
    distances = {}
    if tree_path and Path(tree_path).exists():
        distances = pairwise_leaf_distances(tree_path, leaves)

    recommendations = {}
    candidates = defaultdict(list)
    for key, members in grouped.items():
        species_name = _species(members[0].get('taxonomy'))
        if species_name and hungate_counts_by_group[key] > 0:
            selection_group_type = 'HUNGATE_BASELINE_EXTENSION'
        elif species_name and secondary_counts_by_group[key] > 0:
            selection_group_type = 'SECONDARY_BASELINE_EXTENSION'
        else:
            selection_group_type = 'CANDIDATE_ONLY_GROUP'
        is_baseline_extension = selection_group_type != 'CANDIDATE_ONLY_GROUP'
        for row in members:
            evidence = _evidence_quality(row)
            base = {
                'sequence_id': row['id'],
                'partner_id': row.get('partner_id', 'NA'),
                'assessment_species': group_names[key],
                'quarterly_review_group': key,
                'selection_group_type': selection_group_type,
                'evidence_quality': evidence,
                'available_genomes_before': counts_by_group[key],
                'hungate_baseline_genomes_before': hungate_counts_by_group[key],
                'secondary_baseline_genomes_before': secondary_counts_by_group[key],
                'cultured_rumen_baseline_genomes_before': baseline_counts_by_group[key],
                'pangenome_target': target,
                'coverage_gap_before': max(0, target - counts_by_group[key]),
                'qc_passed_genome_ids': ';'.join(
                    genome_feedback.get(normalise_taxon_name(group_names[key]), {}).get('genome_ids', [])
                ) or 'None',
                'available_ani_clusters': ';'.join(
                    genome_feedback.get(normalise_taxon_name(group_names[key]), {}).get('ani_clusters', [])
                ) or 'None',
                'nearest_available_genome': row.get('nearest_genome_hit', 'NA'),
                'nearest_available_identity': row.get('nearest_genome_identity', 'NA'),
                'nearest_available_identity_source': row.get('nearest_genome_identity_source', 'assessment_snapshot'),
                'baseline_evidence_class': row.get('baseline_evidence_class', 'NA'),
                'hungate_nearest_hit': row.get('hungate_nearest_hit', 'NA'),
                'hungate_nearest_identity': row.get('hungate_nearest_identity', 'NA'),
                'hungate_nearest_query_coverage': row.get('hungate_nearest_query_coverage', 'NA'),
                'secondary_baseline_nearest_hit': row.get('secondary_baseline_nearest_hit', 'NA'),
                'secondary_baseline_nearest_identity': row.get('secondary_baseline_nearest_identity', 'NA'),
                'secondary_baseline_nearest_query_coverage': row.get('secondary_baseline_nearest_query_coverage', 'NA'),
                'cultured_rumen_nearest_hit': row.get('cultured_rumen_nearest_hit', row.get('nearest_hit', 'NA')),
                'cultured_rumen_nearest_identity': row.get('cultured_rumen_nearest_identity', row.get('nearest_identity', 'NA')),
                'cultured_rumen_nearest_query_coverage': row.get('cultured_rumen_nearest_query_coverage', row.get('nearest_query_coverage', 'NA')),
                'baseline_nearest_hit': row.get('nearest_hit', 'NA'),
                'baseline_nearest_identity': row.get('nearest_identity', 'NA'),
                'baseline_nearest_query_coverage': row.get('nearest_query_coverage', 'NA'),
                'baseline_redundancy_identity_threshold': f'{redundancy_identity:.2f}',
                'baseline_redundancy_min_query_coverage': f'{redundancy_coverage:.2f}',
                'baseline_redundancy_status': 'ELIGIBLE',
                'baseline_extension_status': 'NOT_APPLICABLE_CANDIDATE_GROUP',
                'baseline_extension_min_identity': f'{extension_identity:.2f}',
                'baseline_extension_min_query_coverage': f'{extension_coverage:.2f}',
                'project_nearest_identity': row.get('project_nearest_identity', 'NA'),
                'gtdb_reference_nearest_identity': row.get('reference_nearest_identity', 'NA'),
                'mwl_matched_rank': row.get('mwl_matched_rank', 'NA'),
                'mwl_matched_taxon': row.get('mwl_matched_taxon', 'NA'),
                'mwl_score': row.get('mwl_score', 'NA'),
                'local_tree_figure': row.get('local_neighbourhood_figure', 'NA'),
                'source_snapshot': row.get('_snapshot_id') or row.get('_snapshot_source_path') or 'NA',
                'round_rank': 'NA',
                'backup_rank': 'NA',
                'backup_for': 'NA',
                'marginal_patristic_distance': 'NA',
            }
            if _true(row.get('already_sequenced')):
                base.update({
                    'role': 'ALREADY_SEQUENCED', 'priority_tier': 'ALREADY_SEQUENCED',
                    'baseline_redundancy_status': 'NOT_APPLICABLE_COMMITTED',
                    'baseline_extension_status': 'NOT_APPLICABLE_COMMITTED',
                    'recommendation_reason': 'genome already available; excluded from new recommendations',
                })
                recommendations[row['id']] = base
            elif _true(row.get('selected_for_genome_sequencing')):
                base.update({
                    'role': 'ALREADY_SELECTED', 'priority_tier': 'COMMITTED_PENDING',
                    'baseline_redundancy_status': 'NOT_APPLICABLE_COMMITTED',
                    'baseline_extension_status': 'NOT_APPLICABLE_COMMITTED',
                    'recommendation_reason': 'already selected for sequencing; genome is not yet available',
                })
                recommendations[row['id']] = base
            else:
                extension_eligible = False
                extension_status = base['baseline_extension_status']
                extension_reason = 'not a baseline-pangenome extension group'
                if is_baseline_extension:
                    extension_eligible, extension_status, extension_reason = baseline_extension_evidence(
                        row, extension_identity, extension_coverage,
                    )
                    base['baseline_extension_status'] = extension_status

                redundant, identity, coverage = baseline_redundancy_evidence(
                    row, redundancy_identity, redundancy_coverage,
                )
                retained_for_gap = bool(
                    redundant and is_baseline_extension and extension_eligible
                    and base['coverage_gap_before'] > 0
                )
                if retained_for_gap:
                    base['baseline_redundancy_status'] = 'RETAINED_FOR_PANGENOME_GAP'

                review_evidence = evidence == 'LOW' or (evidence == 'MODERATE' and not include_moderate)
                if redundant and not retained_for_gap:
                    base.update({
                        'role': 'BASELINE_REDUNDANT',
                        'priority_tier': 'NEAR_IDENTICAL_CULTURED_BASELINE',
                        'baseline_redundancy_status': 'EXCLUDED_NEAR_IDENTICAL_BASELINE',
                        'baseline_extension_status': (
                            'EXCLUDED_NEAR_IDENTICAL_BASELINE' if is_baseline_extension
                            else base['baseline_extension_status']
                        ),
                        'recommendation_reason': (
                            f'nearest cultured baseline {identity:.2f}% across {coverage:.2f}% of the query; '
                            f'excluded at >={redundancy_identity:.2f}% identity and '
                            f'>={redundancy_coverage:.2f}% query coverage; MWL evidence does not override redundancy'
                        ),
                    })
                    if review_evidence:
                        base['recommendation_reason'] += '; marker evidence also requires review'
                    recommendations[row['id']] = base
                elif review_evidence:
                    base.update({
                        'role': 'REVIEW', 'priority_tier': 'EVIDENCE_REVIEW',
                        'recommendation_reason': f'{evidence.lower()} marker evidence; excluded pending review',
                    })
                    if retained_for_gap:
                        base['baseline_redundancy_status'] = 'RETAINED_FOR_PANGENOME_GAP_PENDING_REVIEW'
                        base['recommendation_reason'] += '; same-species pangenome gap remains'
                    else:
                        base['baseline_redundancy_status'] = 'NOT_EVALUATED_EVIDENCE'
                        base['baseline_extension_status'] = 'NOT_EVALUATED_EVIDENCE'
                    recommendations[row['id']] = base
                elif is_baseline_extension and not extension_eligible:
                    base.update({
                        'role': 'PANGENOME_BOUNDARY_REVIEW',
                        'priority_tier': 'PANGENOME_BOUNDARY_REVIEW',
                        'baseline_redundancy_status': 'NOT_EVALUATED_PANGENOME_BOUNDARY',
                        'recommendation_reason': (
                            f'not admitted to baseline-pangenome extension: {extension_reason}; '
                            'review as a possible separate candidate lineage'
                        ),
                    })
                    recommendations[row['id']] = base
                else:
                    row['_quarterly_review_base'] = base
                    candidates[key].append(row)

    primaries = []
    while len(primaries) < budget:
        group_heads = []
        for key, remaining in candidates.items():
            if not remaining:
                continue
            anchors = anchors_by_group[key]
            winner = min(remaining, key=lambda row: _candidate_key(row, anchors, distances))
            gap = max(0, target - counts_by_group[key])
            nearest = _float(winner.get('nearest_genome_identity'))
            no_close = nearest is not None and 0 < nearest < 98.65
            priority = (
                0 if gap > 0 else 1,
                -gap,
                0 if no_close else 1,
                counts_by_group[key],
                *_candidate_key(winner, anchors, distances),
            )
            group_heads.append((priority, key, winner))
        if not group_heads:
            break
        _, key, winner = min(group_heads, key=lambda item: item[0])
        candidates[key].remove(winner)
        rank = len(primaries) + 1
        available_before = counts_by_group[key]
        gap = max(0, target - available_before)
        patristic = _marginal_distance(winner, anchors_by_group[key], distances)
        tier = _tier(winner, gap)
        result = dict(winner['_quarterly_review_base'])
        result.update({
            'role': 'PRIMARY', 'round_rank': rank, 'priority_tier': tier,
            'available_genomes_before': available_before,
            'coverage_gap_before': gap,
            'marginal_patristic_distance': f'{patristic:.6f}' if patristic is not None else 'NA',
            'recommendation_reason': _reason(winner, tier, available_before, target, patristic),
        })
        recommendations[winner['id']] = result
        primaries.append((key, winner, result))
        anchors_by_group[key].append(_tree_leaf(winner))
        counts_by_group[key] += 1

    used_backups = set()
    for key, primary, primary_result in primaries:
        for backup_rank in range(1, backup_count + 1):
            remaining = [row for row in candidates[key] if row['id'] not in used_backups]
            if not remaining:
                break
            backup = min(remaining, key=lambda row: _candidate_key(row, anchors_by_group[key], distances))
            used_backups.add(backup['id'])
            patristic = _marginal_distance(backup, anchors_by_group[key], distances)
            result = dict(backup['_quarterly_review_base'])
            result.update({
                'role': 'BACKUP', 'round_rank': primary_result['round_rank'],
                'backup_rank': backup_rank, 'backup_for': primary['id'],
                'priority_tier': primary_result['priority_tier'],
                'marginal_patristic_distance': f'{patristic:.6f}' if patristic is not None else 'NA',
                'recommendation_reason': (
                    f'backup for {primary["id"]}; extraction-failure resilience and alternate diversity capture'
                ),
            })
            recommendations[backup['id']] = result

    for key, remaining in candidates.items():
        for row in remaining:
            if row['id'] in recommendations:
                continue
            base = dict(row['_quarterly_review_base'])
            gap = max(0, target - counts_by_group[key])
            base.update({
                'role': 'NOT_SELECTED',
                'priority_tier': _tier(row, gap),
                'recommendation_reason': 'eligible but outside the current genome budget and backup allocation',
            })
            recommendations[row['id']] = base

    role_order = {
        'PRIMARY': 0, 'BACKUP': 1, 'REVIEW': 2, 'NOT_SELECTED': 3,
        'PANGENOME_BOUNDARY_REVIEW': 4, 'BASELINE_REDUNDANT': 5,
        'ALREADY_SELECTED': 6, 'ALREADY_SEQUENCED': 7,
    }
    return sorted(
        recommendations.values(),
        key=lambda row: (
            role_order.get(str(row.get('role')), 99),
            _int(row.get('round_rank'), 10**9),
            _int(row.get('backup_rank'), 10**9),
            str(row.get('sequence_id')),
        ),
    )


QUARTERLY_REVIEW_FIELDS = [
    'round_id', 'sequence_id', 'partner_id', 'role', 'round_rank', 'backup_rank', 'backup_for',
    'priority_tier', 'assessment_species', 'evidence_quality', 'available_genomes_before',
    'hungate_baseline_genomes_before', 'secondary_baseline_genomes_before',
    'cultured_rumen_baseline_genomes_before', 'selection_group_type',
    'qc_passed_genome_ids', 'available_ani_clusters',
    'pangenome_target', 'coverage_gap_before', 'nearest_available_genome',
    'nearest_available_identity', 'nearest_available_identity_source',
    'marginal_patristic_distance', 'baseline_evidence_class',
    'hungate_nearest_hit', 'hungate_nearest_identity', 'hungate_nearest_query_coverage',
    'secondary_baseline_nearest_hit', 'secondary_baseline_nearest_identity',
    'secondary_baseline_nearest_query_coverage',
    'cultured_rumen_nearest_hit', 'cultured_rumen_nearest_identity',
    'cultured_rumen_nearest_query_coverage', 'baseline_nearest_hit',
    'baseline_nearest_identity', 'baseline_nearest_query_coverage',
    'baseline_redundancy_identity_threshold', 'baseline_redundancy_min_query_coverage',
    'baseline_redundancy_status', 'baseline_extension_status',
    'baseline_extension_min_identity', 'baseline_extension_min_query_coverage',
    'project_nearest_identity', 'gtdb_reference_nearest_identity',
    'mwl_matched_rank', 'mwl_matched_taxon', 'mwl_score', 'local_tree_figure', 'source_snapshot',
    'recommendation_reason',
]


def write_quarterly_review_reports(outdir: str | Path, round_id: str, rows: List[dict], parameters: dict) -> dict:
    output = Path(outdir)
    output.mkdir(parents=True, exist_ok=True)
    prepared = []
    for row in rows:
        item = dict(row)
        item['round_id'] = round_id
        prepared.append(item)

    summary = output / 'quarterly_review_summary.tsv'
    with open(summary, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=QUARTERLY_REVIEW_FIELDS, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(prepared)

    selected = output / 'next_genome_set.tsv'
    with open(selected, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=QUARTERLY_REVIEW_FIELDS, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(row for row in prepared if row.get('role') in {'PRIMARY', 'BACKUP'})

    manifest = output / 'quarterly_review_manifest.tsv'
    counts = defaultdict(int)
    for row in prepared:
        counts[str(row.get('role'))] += 1
    with open(manifest, 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['Parameter', 'Value'])
        writer.writerow(['RoundID', round_id])
        for key in sorted(parameters):
            writer.writerow([key, parameters[key]])
        for role in (
            'PRIMARY', 'BACKUP', 'REVIEW', 'NOT_SELECTED', 'PANGENOME_BOUNDARY_REVIEW', 'BASELINE_REDUNDANT',
            'ALREADY_SELECTED', 'ALREADY_SEQUENCED',
        ):
            writer.writerow([f'{role}Count', counts[role]])
        writer.writerow(['ScientificScope', 'Marker-gene diversity expansion; confirm strain diversity with genome-derived ANI/phylogenomics'])

    return {'summary': str(summary), 'selected': str(selected), 'manifest': str(manifest)}
