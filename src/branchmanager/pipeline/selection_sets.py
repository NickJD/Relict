"""Rolling, clade-aware candidate sets for genome-sequencing decisions."""

from __future__ import annotations

import csv
import hashlib
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional

from branchmanager.pipeline.neighbourhood import pairwise_leaf_distances
from branchmanager.taxonomy import normalise_taxon_name, parse_taxon_string


RANKS = [('s', 'species'), ('g', 'genus'), ('f', 'family'), ('o', 'order'),
         ('c', 'class'), ('p', 'phylum'), ('d', 'domain')]
MWL_STRENGTH = {
    'species': 5, 's': 5, 'genus': 4, 'g': 4, 'family': 3, 'f': 3,
    'order': 2, 'o': 2, 'class': 1, 'c': 1, 'phylum': 0, 'p': 0,
    'domain': 0, 'd': 0,
}


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


def _selected(row: dict) -> bool:
    return str(row.get('selected_for_genome_sequencing') or '').strip().lower() == 'true'


def _available(row: dict) -> bool:
    return str(row.get('already_sequenced') or '').strip().lower() == 'true'


def _committed(row: dict) -> bool:
    return _selected(row) or _available(row)


def _evidence_quality(row: dict) -> str:
    flags = str(row.get('placement_flags') or '').upper()
    identity = _float(row.get('classification_identity'))
    confidence = _float(row.get('classification_confidence'))
    coverage = _float(row.get('classification_query_coverage'))
    taxonomy = str(row.get('taxonomy') or '').strip().lower()
    if any(flag in flags for flag in (
        'NO_REFERENCE_HIT', 'NO_CLASSIFICATION', 'CHIMERA', 'VERY_SHORT', 'HIGH_N_CONTENT',
        'MARKER_QC_FAILED', 'MARKER_QC_REVIEW_REQUIRED', 'MARKER_QC_UNVERIFIED',
    )):
        return 'LOW'
    if 'MARKER_QC_REVIEW_APPROVED' in flags:
        return 'MODERATE'
    if identity is not None and identity < 90.0:
        return 'LOW'
    if coverage is not None and coverage < 80.0:
        return 'LOW'
    if identity is None and taxonomy in ('', 'na', 'none'):
        return 'LOW'
    if any(flag in flags for flag in ('LOW_CLASSIFICATION', 'LOW_CONFIDENCE', 'DISAGREEMENT', 'CONFLICT')):
        return 'MODERATE'
    if (identity is not None and identity < 95.0) or (confidence is not None and confidence < 0.8) or (coverage is not None and coverage < 90.0):
        return 'MODERATE'
    return 'HIGH'


def _species(taxonomy: object) -> str:
    value = str(parse_taxon_string(str(taxonomy or '')).get('s') or '').strip()
    return value if value and value.lower() not in ('na', 'none', 'unclassified') else ''


def _group_for_row(row: dict) -> tuple[str, str, str]:
    parsed = parse_taxon_string(str(row.get('taxonomy') or ''))
    species = str(parsed.get('s') or '').strip()
    if species:
        return f'species:{normalise_taxon_name(species)}', 'assessment species', species

    rank_key = 'unclassified'
    taxon = 'unclassified'
    rank_name = 'unclassified'
    for key, name in RANKS[1:]:
        value = str(parsed.get(key) or '').strip()
        if value:
            rank_key = f'{key}:{normalise_taxon_name(value)}'
            taxon = value
            rank_name = f'assessment {name}'
            break
    figure = str(row.get('local_neighbourhood_figure') or 'no_tree_context')
    return f'local:{figure}|{rank_key}', f'local clade + {rank_name}', taxon


def _stable_set_id(group_key: str) -> str:
    digest = hashlib.sha1(group_key.encode('utf-8')).hexdigest()[:8].upper()
    return f'BMSET_{digest}'


def _novel_looking(rows: Iterable[dict]) -> bool:
    for row in rows:
        baseline = _float(row.get('nearest_identity'))
        reference = _float(row.get('reference_nearest_identity'))
        mwl_strength = MWL_STRENGTH.get(str(row.get('mwl_matched_rank') or '').lower(), 0)
        if baseline is not None and baseline < 98.65:
            return True
        if reference is not None and reference < 98.65:
            return True
        if mwl_strength >= 3:
            return True
    return False


def _exceptional_after_target(row: dict) -> bool:
    """Flag unusually divergent or specific-priority evidence after target is met."""
    baseline = _float(row.get('nearest_identity'))
    reference = _float(row.get('reference_nearest_identity'))
    mwl_strength = MWL_STRENGTH.get(str(row.get('mwl_matched_rank') or '').lower(), 0)
    return (
        (baseline is not None and baseline < 97.0)
        or (reference is not None and reference < 94.5)
        or mwl_strength >= 3
    )


def _evidence_tuple(row: dict) -> tuple:
    baseline = _float(row.get('nearest_identity'), 100.0)
    project = _float(row.get('project_nearest_identity'), 100.0)
    reference = _float(row.get('reference_nearest_identity'), 100.0)
    isolation = _float(row.get('phylo_isolation'), 0.0)
    confidence = _float(row.get('classification_confidence'), 0.0)
    mwl = MWL_STRENGTH.get(str(row.get('mwl_matched_rank') or '').lower(), 0)
    return (
        mwl,
        100.0 - float(baseline),
        100.0 - float(project),
        100.0 - float(reference),
        float(isolation),
        float(confidence),
        str(row.get('id') or ''),
    )


def _tree_leaf(row: dict) -> str:
    sid = str(row.get('id') or '')
    representative = str(row.get('cluster_representative') or '')
    if str(row.get('in_tree') or '').lower() == 'yes':
        return sid
    return representative if representative not in ('', 'self', 'duplicate', 'N/A', 'NA') else sid


def _distance(left: str, right: str, distances: Dict[tuple[str, str], float]) -> Optional[float]:
    return distances.get(tuple(sorted((left, right))))


def _farthest_first(rows: List[dict], anchors: List[str], distances: Dict[tuple[str, str], float]) -> List[tuple[dict, Optional[float]]]:
    remaining = list(rows)
    chosen: List[tuple[dict, Optional[float]]] = []
    selected_leaves = [anchor for anchor in anchors if anchor]
    while remaining:
        ranked = []
        for row in remaining:
            leaf = _tree_leaf(row)
            observed = [
                value for anchor in selected_leaves
                if (value := _distance(leaf, anchor, distances)) is not None
            ]
            spread = min(observed) if observed else None
            ranked.append((spread is not None, spread or 0.0, _evidence_tuple(row), row))
        has_spread, spread, _, winner = max(ranked, key=lambda item: (item[0], item[1], item[2]))
        chosen.append((winner, spread if has_spread else None))
        selected_leaves.append(_tree_leaf(winner))
        remaining.remove(winner)
    return chosen


def _committed_ids_by_species(db) -> Dict[str, List[str]]:
    if db is None:
        return {}
    try:
        roles = db.get_dataset_roles()
        with db.connect() as conn:
            rows = conn.execute(
                "SELECT s.id, s.dataset, t.taxonomy, COALESCE(m.selected_for_wgs, 0), "
                "COALESCE(m.selected_for_sequencing, 0) "
                "FROM sequences s LEFT JOIN taxonomy t ON s.id = t.id "
                "LEFT JOIN sequencing_metadata m ON s.id = m.id ORDER BY t.rowid"
            ).fetchall()
    except Exception:
        return {}
    records = {}
    for sid, dataset, taxonomy, available, selected in rows:
        entry = records.setdefault(str(sid), {'dataset': str(dataset or ''), 'taxonomy': '', 'available': False, 'selected': False})
        if taxonomy:
            entry['taxonomy'] = str(taxonomy)
        entry['available'] = entry['available'] or bool(available)
        entry['selected'] = entry['selected'] or bool(selected)
    result: Dict[str, List[str]] = defaultdict(list)
    for sid, entry in records.items():
        role = roles.get(entry['dataset'], {})
        committed = (
            (role.get('role') == 'baseline' and role.get('genomes_available'))
            or entry['available'] or entry['selected']
        )
        species = _species(entry['taxonomy'])
        if committed and species:
            result[normalise_taxon_name(species)].append(sid)
    return {key: sorted(set(values)) for key, values in result.items()}


def build_sequencing_sets(
    assessment_rows: List[dict],
    output_path: str | Path,
    *,
    tree_path: Optional[str | Path] = None,
    db=None,
    pangenome_target: int = 3,
    candidate_set_size: int = 4,
) -> str:
    """Assign rolling primary/backup roles within assessment species or local clades."""
    target = max(1, int(pangenome_target or 3))
    set_size = max(target, int(candidate_set_size or 4))
    grouped: Dict[str, List[dict]] = defaultdict(list)
    group_metadata = {}
    for row in assessment_rows:
        key, basis, taxon = _group_for_row(row)
        grouped[key].append(row)
        group_metadata[key] = (basis, taxon)

    committed_by_species = _committed_ids_by_species(db)
    leaf_ids = [_tree_leaf(row) for row in assessment_rows]
    leaf_ids.extend(
        sid for committed_ids in committed_by_species.values() for sid in committed_ids
    )
    distances = {}
    if tree_path and Path(tree_path).exists():
        try:
            distances = pairwise_leaf_distances(tree_path, leaf_ids)
        except Exception:
            distances = {}
    output_rows = []
    for group_key in sorted(grouped):
        members = grouped[group_key]
        basis, taxon = group_metadata[group_key]
        set_id = _stable_set_id(group_key)
        species = _species(members[0].get('taxonomy'))
        species_key = normalise_taxon_name(species) if species else ''
        committed_ids = committed_by_species.get(species_key, []) if species_key else [
            str(row.get('id')) for row in members if _committed(row)
        ]
        baseline_count = max((_int(row.get('genome_available_count_same_species')) for row in members), default=0)
        selected_count = max((_int(row.get('genome_selected_count_same_species')) for row in members), default=0)
        pending_count = max((_int(row.get('genome_pending_count_same_species')) for row in members), default=0)
        committed_count = max(baseline_count + selected_count + pending_count, len(committed_ids))
        gap = max(0, target - committed_count)
        novel_clade = _novel_looking(members)

        for row in members:
            row.update({
                'sequencing_set_id': set_id,
                'sequencing_set_role': 'ALTERNATE',
                'sequencing_set_rank': 'NA',
                'sequencing_set_reason': 'not selected into the proposed working set',
                'selection_group_basis': basis,
                'selection_group_taxon': taxon,
                'novel_looking_clade': 'Yes' if novel_clade else 'No',
            })

        eligible = []
        for row in members:
            if _available(row):
                row['sequencing_set_role'] = 'SEQUENCED'
                row['sequencing_set_reason'] = 'genome already sequenced and available'
            elif _selected(row):
                row['sequencing_set_role'] = 'COMMITTED'
                row['sequencing_set_reason'] = 'already selected for genome sequencing; genome not yet available'
            elif _evidence_quality(row) != 'HIGH':
                row['sequencing_set_role'] = 'REVIEW_EVIDENCE'
                row['sequencing_set_reason'] = 'marker evidence requires review before selection'
            elif gap == 0 and _exceptional_after_target(row):
                eligible.append(row)
            elif gap == 0:
                row['sequencing_set_role'] = 'TARGET_MET'
                row['sequencing_set_reason'] = f'{target}-genome assessment-species target already met'
            else:
                eligible.append(row)

        ranked = _farthest_first(eligible, committed_ids, distances)
        recommended = ranked[:min(1 if gap == 0 else set_size, len(ranked))]
        for rank, (row, spread) in enumerate(recommended, start=1):
            role = 'DIVERSITY_CANDIDATE' if gap == 0 else ('PRIMARY' if rank <= gap else 'BACKUP')
            row['sequencing_set_role'] = role
            row['sequencing_set_rank'] = str(rank)
            spread_text = f'; patristic spread {spread:.5f}' if spread is not None else ''
            if role == 'PRIMARY':
                reason = f'fills genome {committed_count + rank} of target {target}{spread_text}'
            elif role == 'DIVERSITY_CANDIDATE':
                reason = f'exceptional novelty/MWL evidence after minimum genome target was met{spread_text}'
            else:
                reason = f'backup for DNA-extraction failure and strain-diversity capture{spread_text}'
            row['sequencing_set_reason'] = reason

        for row in members:
            output_rows.append({
                'SetID': set_id,
                'GroupBasis': basis,
                'AssessmentTaxon': taxon,
                'NovelLookingClade': 'Yes' if novel_clade else 'No',
                'BaselineGenomes': baseline_count,
                'SequencedPartnerGenomes': selected_count,
                'SelectedPendingGenomes': pending_count,
                'CommittedGenomeCount': committed_count,
                'CommittedGenomeIDs': ';'.join(committed_ids) or 'None',
                'PangenomeTarget': target,
                'PangenomeGap': gap,
                'CandidateID': row.get('id', 'NA'),
                'PartnerID': row.get('partner_id', 'NA'),
                'SetRole': row.get('sequencing_set_role', 'NA'),
                'SetRank': row.get('sequencing_set_rank', 'NA'),
                'EvidenceQuality': _evidence_quality(row),
                'BaselineNearestIdentity': row.get('nearest_identity', 'NA'),
                'ProjectNearestIdentity': row.get('project_nearest_identity', 'NA'),
                'GTDBReferenceIdentity': row.get('reference_nearest_identity', 'NA'),
                'MWLMatchedRank': row.get('mwl_matched_rank', 'NA'),
                'PhyloIsolation': row.get('phylo_isolation', 'NA'),
                'LocalTreeFigure': row.get('local_neighbourhood_figure', 'NA'),
                'SelectionReason': row.get('sequencing_set_reason', 'NA'),
            })

    role_order = {'PRIMARY': 0, 'BACKUP': 1, 'DIVERSITY_CANDIDATE': 2,
                  'COMMITTED': 3, 'SEQUENCED': 4, 'ALTERNATE': 5,
                  'REVIEW_EVIDENCE': 6, 'TARGET_MET': 7}
    output_rows.sort(key=lambda row: (
        row['SetID'], role_order.get(str(row['SetRole']), 99),
        _int(row['SetRank'], 999), str(row['CandidateID']),
    ))
    fields = list(output_rows[0]) if output_rows else [
        'SetID', 'GroupBasis', 'AssessmentTaxon', 'NovelLookingClade', 'BaselineGenomes',
        'SequencedPartnerGenomes', 'SelectedPendingGenomes', 'CommittedGenomeCount', 'CommittedGenomeIDs',
        'PangenomeTarget', 'PangenomeGap', 'CandidateID', 'PartnerID', 'SetRole',
        'SetRank', 'EvidenceQuality', 'BaselineNearestIdentity', 'ProjectNearestIdentity',
        'GTDBReferenceIdentity', 'MWLMatchedRank', 'PhyloIsolation', 'LocalTreeFigure',
        'SelectionReason',
    ]
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(output_rows)
    return str(path)
