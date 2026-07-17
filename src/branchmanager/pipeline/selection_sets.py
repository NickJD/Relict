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
EVIDENCE_STRENGTH = {'HIGH': 2, 'MODERATE': 1, 'LOW': 0}
DEFAULT_PANGENOME_TARGET = 9
DEFAULT_CANDIDATE_SET_SIZE = 9
DEFAULT_BASELINE_REDUNDANCY_IDENTITY = 99.8
DEFAULT_BASELINE_REDUNDANCY_QUERY_COVERAGE = 95.0
DEFAULT_BASELINE_EXTENSION_MIN_IDENTITY = 98.65
DEFAULT_BASELINE_EXTENSION_MIN_QUERY_COVERAGE = 95.0


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
    flag_set = {f.strip() for f in str(row.get('placement_flags') or '').upper().split(';') if f.strip()}
    identity = _float(row.get('classification_identity'))
    confidence = _float(row.get('classification_confidence'))
    coverage = _float(row.get('classification_query_coverage'))
    taxonomy = str(row.get('taxonomy') or '').strip().lower()
    # Use exact flag matching to avoid CHIMERA_INDETERMINATE being caught by a 'CHIMERA' substring check.
    if flag_set & {
        'NO_REFERENCE_HIT', 'NO_CLASSIFICATION', 'CHIMERA', 'CHIMERA_CONFIRMED',
        'VERY_SHORT', 'HIGH_N_CONTENT', 'MARKER_QC_FAILED',
        'MARKER_QC_REVIEW_REQUIRED', 'MARKER_QC_UNVERIFIED',
    }:
        return 'LOW'
    # MARKER_QC_REVIEW_APPROVED: manual approval has addressed QC concerns; fall through to
    # standard identity/coverage checks so clean assemblies can still reach HIGH.
    if identity is not None and identity < 90.0:
        return 'LOW'
    if coverage is not None and coverage < 80.0:
        return 'LOW'
    if identity is None and taxonomy in ('', 'na', 'none'):
        return 'LOW'
    if flag_set & {'LOW_CLASSIFICATION', 'LOW_CONFIDENCE', 'DISAGREEMENT', 'CONFLICT'}:
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


def _stable_set_id(group_key: str, prefix: str = 'BMSET') -> str:
    digest = hashlib.sha1(group_key.encode('utf-8')).hexdigest()[:8].upper()
    return f'{prefix}_{digest}'


def _novel_looking(rows: Iterable[dict]) -> bool:
    for row in rows:
        evidence_class = str(row.get('baseline_evidence_class') or '').upper()
        if evidence_class in {'NOVEL TO BOTH BASELINES', 'HUNGATE GAP / SECONDARY COVERED'}:
            return True
        baseline = _float(row.get('nearest_identity'))
        reference = _float(row.get('reference_nearest_identity'))
        if baseline is not None and baseline < 98.65:
            return True
        if reference is not None and reference < 98.65:
            return True
    return False


def _first_present(*values):
    for value in values:
        if value not in (None, '', 'NA', 'None'):
            return value
    return None


def _baseline_hit_records(row: dict, *, include_combined: bool = True) -> list[dict]:
    records = [
        {
            'tier': 'priority',
            'label': 'Hungate',
            'hit': row.get('hungate_nearest_hit'),
            'identity': _float(row.get('hungate_nearest_identity')),
            'coverage': _float(row.get('hungate_nearest_query_coverage')),
            'taxonomy': row.get('hungate_nearest_hit_taxonomy'),
            'dataset': row.get('hungate_nearest_hit_dataset'),
        },
        {
            'tier': 'secondary',
            'label': 'SecondaryBaseline',
            'hit': row.get('secondary_baseline_nearest_hit'),
            'identity': _float(row.get('secondary_baseline_nearest_identity')),
            'coverage': _float(row.get('secondary_baseline_nearest_query_coverage')),
            'taxonomy': row.get('secondary_baseline_nearest_hit_taxonomy'),
            'dataset': row.get('secondary_baseline_nearest_hit_dataset'),
        },
    ]
    if include_combined:
        records.append({
            'tier': 'cultured_rumen',
            'label': 'CulturedRumen',
            'hit': _first_present(row.get('cultured_rumen_nearest_hit'), row.get('nearest_hit')),
            'identity': _float(_first_present(row.get('cultured_rumen_nearest_identity'), row.get('nearest_identity'))),
            'coverage': _float(_first_present(row.get('cultured_rumen_nearest_query_coverage'), row.get('nearest_query_coverage'))),
            'taxonomy': _first_present(row.get('cultured_rumen_nearest_hit_taxonomy'), row.get('nearest_hit_taxonomy')),
            'dataset': _first_present(row.get('cultured_rumen_nearest_hit_dataset'), row.get('nearest_hit_dataset')),
        })
    deduped = []
    seen = set()
    for record in records:
        hit = str(record.get('hit') or '').strip()
        if hit in {'', 'NA', 'None'}:
            continue
        key = (record['label'], hit)
        if key in seen:
            continue
        seen.add(key)
        deduped.append(record)
    return deduped


def _baseline_evidence_rank(row: dict) -> int:
    evidence_class = str(row.get('baseline_evidence_class') or '').strip().upper()
    ranks = {
        'NOVEL TO BOTH BASELINES': 4,
        'HUNGATE GAP / SECONDARY COVERED': 3,
        'HUNGATE COVERED': 2,
        'NO CULTURED BASELINE AVAILABLE': 1,
        'CULTURED BASELINE REDUNDANT': 0,
    }
    if evidence_class in ranks:
        return ranks[evidence_class]
    identity = _float(row.get('nearest_identity'))
    coverage = _float(row.get('nearest_query_coverage'))
    if identity is not None and coverage is not None and identity >= DEFAULT_BASELINE_REDUNDANCY_IDENTITY and coverage >= DEFAULT_BASELINE_REDUNDANCY_QUERY_COVERAGE:
        return 0
    if identity is not None and identity < DEFAULT_BASELINE_EXTENSION_MIN_IDENTITY:
        return 4
    return 1


def baseline_redundancy_evidence(
    row: dict,
    identity_threshold: float = DEFAULT_BASELINE_REDUNDANCY_IDENTITY,
    min_query_coverage: float = DEFAULT_BASELINE_REDUNDANCY_QUERY_COVERAGE,
) -> tuple[bool, Optional[float], Optional[float]]:
    """Return whether a candidate is near-identical to either cultured baseline tier."""
    hits = _baseline_hit_records(row, include_combined=False)
    if not hits:
        hits = _baseline_hit_records(row, include_combined=True)
    if not hits:
        identity = _float(_first_present(row.get('cultured_rumen_nearest_identity'), row.get('nearest_identity')))
        coverage = _float(_first_present(row.get('cultured_rumen_nearest_query_coverage'), row.get('nearest_query_coverage')))
        redundant = bool(
            identity is not None and coverage is not None
            and identity >= float(identity_threshold)
            and coverage >= float(min_query_coverage)
        )
        return redundant, identity, coverage
    observed = [record for record in hits if record.get('identity') is not None]
    best = max(observed, key=lambda record: record['identity'], default=None)
    redundant_hits = [
        record for record in hits
        if record.get('identity') is not None
        and record.get('coverage') is not None
        and record['identity'] >= float(identity_threshold)
        and record['coverage'] >= float(min_query_coverage)
    ]
    if redundant_hits:
        winner = max(redundant_hits, key=lambda record: (record['identity'], record['coverage'] or 0.0))
        return True, winner.get('identity'), winner.get('coverage')
    return False, (best or {}).get('identity'), (best or {}).get('coverage')


def baseline_extension_evidence(
    row: dict,
    min_identity: float = DEFAULT_BASELINE_EXTENSION_MIN_IDENTITY,
    min_query_coverage: float = DEFAULT_BASELINE_EXTENSION_MIN_QUERY_COVERAGE,
) -> tuple[bool, str, str]:
    """Test whether a candidate can extend either cultured-rumen baseline tier."""
    candidate_species = _species(row.get('taxonomy'))
    if not candidate_species:
        return False, 'CANDIDATE_SPECIES_UNRESOLVED', 'candidate GTDB species is unresolved'
    hits = _baseline_hit_records(row, include_combined=False)
    if not hits:
        hits = _baseline_hit_records(row, include_combined=True)
    if not hits:
        return False, 'BASELINE_HIT_UNAVAILABLE', 'no cultured-baseline hit is available'

    failures: list[tuple[str, str]] = []
    for record in hits:
        hit = str(record.get('hit') or '').strip()
        baseline_species = _species(record.get('taxonomy'))
        label = str(record.get('label') or 'Baseline')
        status_prefix = {
            'Hungate': 'HUNGATE_BASELINE',
            'SecondaryBaseline': 'SECONDARY_BASELINE',
            'CulturedRumen': 'BASELINE',
        }.get(label, 'BASELINE')
        eligible_status = {
            'Hungate': 'ELIGIBLE_HUNGATE_BASELINE_EXTENSION',
            'SecondaryBaseline': 'ELIGIBLE_SECONDARY_BASELINE_EXTENSION',
            'CulturedRumen': 'ELIGIBLE_BASELINE_PANGENOME_EXTENSION',
        }.get(label, 'ELIGIBLE_BASELINE_PANGENOME_EXTENSION')
        if not baseline_species:
            failures.append((
                f'{status_prefix}_SPECIES_UNRESOLVED',
                f'{label} hit {hit} has unresolved species',
            ))
            continue
        if normalise_taxon_name(candidate_species) != normalise_taxon_name(baseline_species):
            failures.append((
                f'{status_prefix}_SPECIES_MISMATCH',
                f'candidate species {candidate_species} differs from {label} species {baseline_species}',
            ))
            continue
        identity = record.get('identity')
        if identity is None or identity < float(min_identity):
            observed = 'unavailable' if identity is None else f'{identity:.2f}%'
            failures.append((
                f'{status_prefix}_IDENTITY_BELOW_EXTENSION_THRESHOLD',
                f'{label} identity {observed} is below {float(min_identity):.2f}%',
            ))
            continue
        coverage = record.get('coverage')
        if coverage is None or coverage < float(min_query_coverage):
            observed = 'unavailable' if coverage is None else f'{coverage:.2f}%'
            failures.append((
                f'{status_prefix}_COVERAGE_BELOW_EXTENSION_THRESHOLD',
                f'{label} query coverage {observed} is below {float(min_query_coverage):.2f}%',
            ))
            continue
        return (
            True,
            eligible_status,
            f'exact species agreement with {hit}; {identity:.2f}% identity across {coverage:.2f}% of the query',
        )

    if failures:
        return False, failures[0][0], failures[0][1]
    return False, 'BASELINE_EXTENSION_CRITERIA_NOT_MET', 'no cultured-baseline tier met extension criteria'


def _evidence_tuple(row: dict) -> tuple:
    baseline = _float(row.get('nearest_identity'), 100.0)
    project = _float(row.get('project_nearest_identity'), 100.0)
    reference = _float(row.get('reference_nearest_identity'), 100.0)
    isolation = _float(row.get('phylo_isolation'), 0.0)
    confidence = _float(row.get('classification_confidence'), 0.0)
    mwl = MWL_STRENGTH.get(str(row.get('mwl_matched_rank') or '').lower(), 0)
    mwl_score = _float(row.get('mwl_score'), 0.0)
    return (
        _baseline_evidence_rank(row),
        EVIDENCE_STRENGTH.get(_evidence_quality(row), 0),
        100.0 - float(baseline),
        100.0 - float(project),
        100.0 - float(reference),
        float(isolation),
        mwl,
        float(mwl_score),
        float(confidence),
        str(row.get('id') or ''),
    )


def _ranking_reason(row: dict, opening: str, spread: Optional[float]) -> str:
    parts = [opening]
    if spread is not None:
        parts.append(f'marginal patristic distance {spread:.5f}')
    baseline = _float(row.get('nearest_identity'))
    if baseline is not None and baseline > 0:
        parts.append(f'nearest cultured baseline {baseline:.2f}%')
    project = _float(row.get('project_nearest_identity'))
    if project is not None and project > 0:
        parts.append(f'nearest project marker {project:.2f}%')
    mwl_rank = str(row.get('mwl_matched_rank') or '').strip()
    if mwl_rank not in {'', 'NA', 'None'}:
        parts.append(f'MWL {mwl_rank}-level context')
    return '; '.join(parts)


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
    pangenome_target: int = DEFAULT_PANGENOME_TARGET,
    candidate_set_size: int = DEFAULT_CANDIDATE_SET_SIZE,
    baseline_redundancy_identity: float = DEFAULT_BASELINE_REDUNDANCY_IDENTITY,
    baseline_redundancy_min_query_coverage: float = DEFAULT_BASELINE_REDUNDANCY_QUERY_COVERAGE,
    baseline_extension_min_identity: float = DEFAULT_BASELINE_EXTENSION_MIN_IDENTITY,
    baseline_extension_min_query_coverage: float = DEFAULT_BASELINE_EXTENSION_MIN_QUERY_COVERAGE,
) -> str:
    """Assign rolling nine-member diversity panels within species or local clades."""
    target = max(1, int(pangenome_target or DEFAULT_PANGENOME_TARGET))
    set_size = max(target, int(candidate_set_size or DEFAULT_CANDIDATE_SET_SIZE))
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
        species = _species(members[0].get('taxonomy'))
        species_key = normalise_taxon_name(species) if species else ''
        committed_ids = committed_by_species.get(species_key, []) if species_key else [
            str(row.get('id')) for row in members if _committed(row)
        ]
        hungate_count = max((_int(row.get('hungate_genome_count_same_species')) for row in members), default=0)
        secondary_count = max((_int(row.get('secondary_baseline_genome_count_same_species')) for row in members), default=0)
        prior_baseline_count = max((_int(row.get('genome_available_count_same_species')) for row in members), default=0)
        if hungate_count + secondary_count == 0 and prior_baseline_count:
            hungate_count = prior_baseline_count
        baseline_count = hungate_count + secondary_count
        selected_count = max((_int(row.get('genome_selected_count_same_species')) for row in members), default=0)
        pending_count = max((_int(row.get('genome_pending_count_same_species')) for row in members), default=0)
        committed_count = max(baseline_count + selected_count + pending_count, len(committed_ids))
        gap = max(0, target - committed_count)
        novel_clade = _novel_looking(members)
        if species and hungate_count > 0:
            group_type = 'HUNGATE_BASELINE_EXTENSION'
            basis = 'Hungate-anchored pangenome extension: exact species + close full-length marker'
            set_prefix = 'BMEXT'
        elif species and secondary_count > 0:
            group_type = 'SECONDARY_BASELINE_EXTENSION'
            basis = 'secondary cultured-rumen baseline extension: exact species + close full-length marker'
            set_prefix = 'BMEXT'
        else:
            group_type = 'CANDIDATE_ONLY_GROUP'
            set_prefix = 'BMSET'
        is_baseline_extension = group_type != 'CANDIDATE_ONLY_GROUP'
        set_id = _stable_set_id(group_key, set_prefix)

        def same_species_hits(hit_key: str, tax_key: str) -> set[str]:
            return {
                str(row.get(hit_key)) for row in members
                if str(row.get(hit_key) or '') not in {'', 'NA', 'None'}
                and _species(row.get(tax_key))
                and normalise_taxon_name(_species(row.get(tax_key))) == species_key
            }

        hungate_anchor_ids = sorted(same_species_hits('hungate_nearest_hit', 'hungate_nearest_hit_taxonomy'))
        secondary_anchor_ids = sorted(same_species_hits('secondary_baseline_nearest_hit', 'secondary_baseline_nearest_hit_taxonomy'))
        cultured_anchor_ids = sorted(same_species_hits('nearest_hit', 'nearest_hit_taxonomy'))
        baseline_anchor_ids = sorted(set(hungate_anchor_ids + secondary_anchor_ids + cultured_anchor_ids))

        for row in members:
            row.update({
                'sequencing_set_id': set_id,
                'sequencing_set_role': 'ALTERNATE',
                'sequencing_set_rank': 'NA',
                'sequencing_set_reason': 'not selected into the proposed working set',
                'selection_group_basis': basis,
                'selection_group_type': group_type,
                'selection_group_taxon': taxon,
                'novel_looking_clade': 'Yes' if novel_clade else 'No',
                'baseline_redundancy_status': 'ELIGIBLE',
                'baseline_redundancy_identity_threshold': f'{redundancy_identity:.2f}',
                'baseline_redundancy_min_query_coverage': f'{redundancy_coverage:.2f}',
                'baseline_extension_status': 'NOT_APPLICABLE_CANDIDATE_GROUP',
                'baseline_extension_min_identity': f'{extension_identity:.2f}',
                'baseline_extension_min_query_coverage': f'{extension_coverage:.2f}',
                'selection_diversity_distance': 'NA',
            })

        eligible = []
        redundant_count = 0
        for row in members:
            if _available(row):
                row['sequencing_set_role'] = 'SEQUENCED'
                row['sequencing_set_reason'] = 'genome already sequenced and available'
                row['baseline_redundancy_status'] = 'NOT_APPLICABLE_COMMITTED'
                row['baseline_extension_status'] = 'NOT_APPLICABLE_COMMITTED'
                continue
            if _selected(row):
                row['sequencing_set_role'] = 'COMMITTED'
                row['sequencing_set_reason'] = 'already selected for genome sequencing; genome not yet available'
                row['baseline_redundancy_status'] = 'NOT_APPLICABLE_COMMITTED'
                row['baseline_extension_status'] = 'NOT_APPLICABLE_COMMITTED'
                continue

            evidence_quality = _evidence_quality(row)
            extension_eligible = False
            extension_status = row.get('baseline_extension_status', 'NOT_APPLICABLE_CANDIDATE_GROUP')
            extension_reason = 'not a baseline-pangenome extension group'
            if is_baseline_extension:
                extension_eligible, extension_status, extension_reason = baseline_extension_evidence(
                    row, extension_identity, extension_coverage,
                )
                row['baseline_extension_status'] = extension_status

            redundant, identity, coverage = baseline_redundancy_evidence(
                row, redundancy_identity, redundancy_coverage,
            )
            retained_for_gap = bool(
                redundant and is_baseline_extension and extension_eligible and gap > 0
            )
            if retained_for_gap:
                row['baseline_redundancy_status'] = 'RETAINED_FOR_PANGENOME_GAP'

            if redundant and not retained_for_gap:
                redundant_count += 1
                row['sequencing_set_role'] = 'BASELINE_REDUNDANT'
                row['baseline_redundancy_status'] = 'EXCLUDED_NEAR_IDENTICAL_BASELINE'
                if is_baseline_extension:
                    row['baseline_extension_status'] = 'EXCLUDED_NEAR_IDENTICAL_BASELINE'
                row['sequencing_set_reason'] = (
                    f'nearest cultured baseline {identity:.2f}% across {coverage:.2f}% of the query; '
                    f'excluded at >={redundancy_identity:.2f}% identity and '
                    f'>={redundancy_coverage:.2f}% query coverage; MWL evidence does not override redundancy'
                )
                if evidence_quality == 'LOW':
                    row['sequencing_set_reason'] += '; marker evidence also requires review'
            elif evidence_quality == 'LOW':
                row['sequencing_set_role'] = 'REVIEW_EVIDENCE'
                row['sequencing_set_reason'] = 'marker evidence requires review before selection'
                if retained_for_gap:
                    row['baseline_redundancy_status'] = 'RETAINED_FOR_PANGENOME_GAP_PENDING_REVIEW'
                    row['sequencing_set_reason'] += '; same-species pangenome gap remains'
                else:
                    row['baseline_redundancy_status'] = 'NOT_EVALUATED_EVIDENCE'
                    row['baseline_extension_status'] = 'NOT_EVALUATED_EVIDENCE'
            elif is_baseline_extension and not extension_eligible:
                row['sequencing_set_role'] = 'PANGENOME_BOUNDARY_REVIEW'
                row['baseline_redundancy_status'] = 'NOT_EVALUATED_PANGENOME_BOUNDARY'
                row['sequencing_set_reason'] = (
                    f'not admitted to baseline-pangenome extension: {extension_reason}; '
                    'review as a possible separate candidate lineage'
                )
            else:
                eligible.append(row)

        ranked = _farthest_first(eligible, committed_ids, distances)
        recommended = ranked[:min(set_size, len(ranked))]
        for rank, (row, spread) in enumerate(recommended, start=1):
            role = 'DIVERSITY_CANDIDATE' if gap == 0 else ('PRIMARY' if rank <= gap else 'BACKUP')
            row['sequencing_set_role'] = role
            row['sequencing_set_rank'] = str(rank)
            row['selection_diversity_distance'] = f'{spread:.6f}' if spread is not None else 'NA'
            if role == 'PRIMARY':
                opening = f'fills genome {committed_count + rank} of target {target}'
            elif role == 'DIVERSITY_CANDIDATE':
                opening = f'ranked diversity expansion after the {target}-genome count target was met'
            else:
                opening = 'ranked backup for DNA-extraction failure and additional within-group diversity'
            row['sequencing_set_reason'] = _ranking_reason(row, opening, spread)

        ranked_count = len(recommended)
        unfilled_ranks = max(0, set_size - ranked_count)

        for row in members:
            output_rows.append({
                'SetID': set_id,
                'GroupType': group_type,
                'GroupBasis': basis,
                'AssessmentTaxon': taxon,
                'BaselineAnchorIDs': ';'.join(baseline_anchor_ids) or 'None',
                'HungateAnchorIDs': ';'.join(hungate_anchor_ids) or 'None',
                'SecondaryBaselineAnchorIDs': ';'.join(secondary_anchor_ids) or 'None',
                'NovelLookingClade': 'Yes' if novel_clade else 'No',
                'HungateBaselineGenomes': hungate_count,
                'SecondaryBaselineGenomes': secondary_count,
                'CulturedRumenBaselineGenomes': baseline_count,
                'BaselineGenomes': baseline_count,
                'SequencedPartnerGenomes': selected_count,
                'SelectedPendingGenomes': pending_count,
                'CommittedGenomeCount': committed_count,
                'CommittedGenomeIDs': ';'.join(committed_ids) or 'None',
                'PangenomeTarget': target,
                'PangenomeGap': gap,
                'DiversityPanelSize': set_size,
                'EligibleCandidateCount': len(eligible),
                'RankedCandidateCount': ranked_count,
                'BaselineRedundantCount': redundant_count,
                'UnfilledPanelRanks': unfilled_ranks,
                'CandidateID': row.get('id', 'NA'),
                'PartnerID': row.get('partner_id', 'NA'),
                'SetRole': row.get('sequencing_set_role', 'NA'),
                'SetRank': row.get('sequencing_set_rank', 'NA'),
                'EvidenceQuality': _evidence_quality(row),
                'BaselineEvidenceClass': row.get('baseline_evidence_class', 'NA'),
                'HungateNearestIdentity': row.get('hungate_nearest_identity', 'NA'),
                'HungateNearestQueryCoverage': row.get('hungate_nearest_query_coverage', 'NA'),
                'SecondaryBaselineNearestIdentity': row.get('secondary_baseline_nearest_identity', 'NA'),
                'SecondaryBaselineNearestQueryCoverage': row.get('secondary_baseline_nearest_query_coverage', 'NA'),
                'CulturedRumenNearestIdentity': row.get('cultured_rumen_nearest_identity', row.get('nearest_identity', 'NA')),
                'CulturedRumenNearestQueryCoverage': row.get('cultured_rumen_nearest_query_coverage', row.get('nearest_query_coverage', 'NA')),
                'BaselineNearestIdentity': row.get('nearest_identity', 'NA'),
                'BaselineNearestQueryCoverage': row.get('nearest_query_coverage', 'NA'),
                'BaselineRedundancyStatus': row.get('baseline_redundancy_status', 'NA'),
                'BaselineRedundancyIdentityThreshold': f'{redundancy_identity:.2f}',
                'BaselineRedundancyMinQueryCoverage': f'{redundancy_coverage:.2f}',
                'BaselineExtensionStatus': row.get('baseline_extension_status', 'NA'),
                'BaselineExtensionMinIdentity': f'{extension_identity:.2f}',
                'BaselineExtensionMinQueryCoverage': f'{extension_coverage:.2f}',
                'ProjectNearestIdentity': row.get('project_nearest_identity', 'NA'),
                'GTDBReferenceIdentity': row.get('reference_nearest_identity', 'NA'),
                'MWLMatchedRank': row.get('mwl_matched_rank', 'NA'),
                'MWLScore': row.get('mwl_score', 'NA'),
                'PhyloIsolation': row.get('phylo_isolation', 'NA'),
                'SelectionDiversityDistance': row.get('selection_diversity_distance', 'NA'),
                'LocalTreeFigure': row.get('local_neighbourhood_figure', 'NA'),
                'SelectionReason': row.get('sequencing_set_reason', 'NA'),
            })

    role_order = {'PRIMARY': 0, 'BACKUP': 1, 'DIVERSITY_CANDIDATE': 2,
                  'COMMITTED': 3, 'SEQUENCED': 4, 'ALTERNATE': 5,
                  'REVIEW_EVIDENCE': 6, 'PANGENOME_BOUNDARY_REVIEW': 7,
                  'BASELINE_REDUNDANT': 8, 'TARGET_MET': 9}
    output_rows.sort(key=lambda row: (
        row['SetID'], role_order.get(str(row['SetRole']), 99),
        _int(row['SetRank'], 999), str(row['CandidateID']),
    ))
    fields = list(output_rows[0]) if output_rows else [
        'SetID', 'GroupType', 'GroupBasis', 'AssessmentTaxon', 'BaselineAnchorIDs',
        'HungateAnchorIDs', 'SecondaryBaselineAnchorIDs',
        'NovelLookingClade', 'HungateBaselineGenomes', 'SecondaryBaselineGenomes',
        'CulturedRumenBaselineGenomes', 'BaselineGenomes',
        'SequencedPartnerGenomes', 'SelectedPendingGenomes', 'CommittedGenomeCount', 'CommittedGenomeIDs',
        'PangenomeTarget', 'PangenomeGap', 'DiversityPanelSize', 'EligibleCandidateCount',
        'RankedCandidateCount', 'BaselineRedundantCount', 'UnfilledPanelRanks',
        'CandidateID', 'PartnerID', 'SetRole', 'SetRank', 'EvidenceQuality',
        'BaselineEvidenceClass', 'HungateNearestIdentity', 'HungateNearestQueryCoverage',
        'SecondaryBaselineNearestIdentity', 'SecondaryBaselineNearestQueryCoverage',
        'CulturedRumenNearestIdentity', 'CulturedRumenNearestQueryCoverage',
        'BaselineNearestIdentity', 'BaselineNearestQueryCoverage', 'BaselineRedundancyStatus',
        'BaselineRedundancyIdentityThreshold', 'BaselineRedundancyMinQueryCoverage',
        'BaselineExtensionStatus', 'BaselineExtensionMinIdentity',
        'BaselineExtensionMinQueryCoverage',
        'ProjectNearestIdentity', 'GTDBReferenceIdentity', 'MWLMatchedRank', 'MWLScore', 'PhyloIsolation',
        'SelectionDiversityDistance', 'LocalTreeFigure', 'SelectionReason',
    ]
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(output_rows)
    return str(path)
