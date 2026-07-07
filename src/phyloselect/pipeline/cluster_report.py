"""
cluster_report.py — Cluster-level summarisation and per-cluster detail files.

After `phyloselect run --collapse`, multiple sequences are grouped into clusters
(representatives + members).  This module:

  1. Aggregates per-sequence assessment rows into *cluster-level summaries*
     (cluster_summary.tsv) so humans can prioritise whole groups at a glance.

  2. Writes per-cluster detail CSV files (clusters/<rep_id>.csv) containing
     every member's individual metrics, including backup rank.

  3. Computes a *phylogenetic isolation* score from the Newick tree using the
     leaf's own branch length as a proxy for how much unique evolution it
     represents.  Longer branch ⟹ more isolated ⟹ higher investigation
     priority.

  4. Produces an *InvestigationScore* composite that combines:
       • NoveltyScore       (sequence-similarity distance from known DB)
       • PhyloIsolation     (branch-length isolation in the tree)
       • TaxonomyConflict   (alt-ref DBs disagree on assignment)
       • NeighborhoodSparsity (density counts at 97%/99%)

  5. Ranks **backup candidates** within each cluster so that if the primary
     candidate cannot be genome-sequenced (DNA quality failure, PCR dropout,
     chimera confirmed, etc.) the user immediately knows which sequence to
     try next.  Written to backup_candidates.tsv (wide format, one row per
     primary) and embedded in each per-cluster CSV.

Suggested additional analyses (see docstring at bottom of file).
"""

from __future__ import annotations

import csv
import re
import statistics
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# ─────────────────────────────────────────────────────────────────────────────
# Minimal Newick parser — extract (leaf_name, own_branch_length) pairs
# ─────────────────────────────────────────────────────────────────────────────

def _tokenise(nwk: str) -> list:
    """Return a flat list of tokens from a Newick string."""
    tokens = []
    i = 0
    n = len(nwk)
    while i < n:
        c = nwk[i]
        if c in '(),;':
            tokens.append(c)
            i += 1
        elif c in ' \t\r\n':
            i += 1
        else:
            # Collect a label+branch-length token until a delimiter
            j = i
            while j < n and nwk[j] not in '(),;\t\r\n':
                j += 1
            tokens.append(nwk[i:j].strip())
            i = j
    return [t for t in tokens if t]


class _Node:
    __slots__ = ('name', 'branch_length', 'children')

    def __init__(self, name: str = '', branch_length: float = 0.0):
        self.name = name
        self.branch_length = branch_length
        self.children: list[_Node] = []


def _parse_node(tokens: list, pos: int) -> Tuple[_Node, int]:
    """Recursive-descent parse; returns (node, next_pos)."""
    node = _Node()
    if pos < len(tokens) and tokens[pos] == '(':
        pos += 1  # consume '('
        while True:
            child, pos = _parse_node(tokens, pos)
            node.children.append(child)
            if pos >= len(tokens):
                break
            if tokens[pos] == ')':
                pos += 1  # consume ')'
                break
            if tokens[pos] == ',':
                pos += 1  # consume ','
    # Parse optional label and branch length that follow the node body
    if pos < len(tokens) and tokens[pos] not in ('(', ')', ',', ';'):
        token = tokens[pos]
        pos += 1
        # token may be  Label:0.123  or  :0.123  or  Label  or  0.963  (bootstrap)
        if ':' in token:
            parts = token.split(':', 1)
            node.name = parts[0]
            try:
                node.branch_length = float(parts[1])
            except ValueError:
                node.branch_length = 0.0
        else:
            # lone label (internal) or bootstrap value (ignored)
            node.name = token
    return node, pos


def _collect_leaves(node: _Node, leaf_branches: Dict[str, float]) -> None:
    if not node.children:
        # leaf
        name = node.name.strip()
        if name and not re.match(r'^\d*\.?\d+$', name):
            # skip purely-numeric strings (bootstrap values that leaked through)
            leaf_branches[name] = node.branch_length
        return
    for child in node.children:
        _collect_leaves(child, leaf_branches)


def compute_phylogenetic_isolation(tree_path: str) -> Dict[str, float]:
    """Parse *tree_path* (Newick) and return {leaf_id: own_branch_length}.

    The leaf's own branch length to its parent is used as a proxy for
    phylogenetic isolation — a long terminal branch means the sequence has no
    close relative in the tree and represents unique evolutionary history.

    Returns an empty dict if the file is missing or unparseable.
    """
    try:
        nwk = Path(tree_path).read_text()
    except Exception:
        return {}

    # Strip newlines / trailing semicolons
    nwk = nwk.strip().rstrip(';')
    if not nwk:
        return {}

    try:
        tokens = _tokenise(nwk)
        root, _ = _parse_node(tokens, 0)
    except Exception:
        return {}

    leaf_branches: Dict[str, float] = {}
    _collect_leaves(root, leaf_branches)

    # Normalise: scale so the max = 1.0 to make scores comparable across runs
    if not leaf_branches:
        return {}
    max_bl = max(leaf_branches.values()) or 1.0
    return {k: round(v / max_bl, 6) for k, v in leaf_branches.items()}


# ─────────────────────────────────────────────────────────────────────────────
# Composite investigation score
# ─────────────────────────────────────────────────────────────────────────────

def compute_investigation_composite(
    novelty_score: float,
    phylo_isolation: float,           # 0–1 normalised
    matches_ge_97: int,
    matches_ge_99: int,
    taxonomy_conflict: bool = False,  # alt-ref DBs disagree?
    min_nearest_identity: float = 100.0,
) -> float:
    """Compute a composite *InvestigationScore* on a 0–100 scale.

    Weighting rationale
    -------------------
    • NoveltyScore (30 pts max)  — already encodes distance + density
    • PhyloIsolation (25 pts max) — uniqueness in the evolutionary tree
    • NeighborhoodSparsity (25 pts max) — how many near-identical seqs exist
    • TaxonomyConflict (10 pts) — disagreement between ref DBs signals ambiguity
    • DistanceBonus (10 pts) — raw sequence distance from nearest known seq

    Higher score = more worth investigating.
    """
    # Contribution 1: NoveltyScore (already 0–100), scaled to 30 pts
    c_novelty = min(30.0, (novelty_score / 100.0) * 30.0)

    # Contribution 2: Phylogenetic isolation (0–1 → 0–25 pts)
    c_phylo = min(25.0, phylo_isolation * 25.0)

    # Contribution 3: Neighbourhood sparsity (no near-identical neighbours)
    c_sparse = 0.0
    if matches_ge_99 == 0:
        c_sparse += 15.0
    elif matches_ge_99 <= 2:
        c_sparse += 8.0
    if matches_ge_97 <= 3:
        c_sparse += 10.0
    elif matches_ge_97 <= 8:
        c_sparse += 5.0
    c_sparse = min(25.0, c_sparse)

    # Contribution 4: Taxonomy conflict bonus
    c_conflict = 10.0 if taxonomy_conflict else 0.0

    # Contribution 5: Raw distance bonus (below 97% = truly novel lineage)
    c_distance = 0.0
    if min_nearest_identity < 90.0:
        c_distance = 10.0
    elif min_nearest_identity < 95.0:
        c_distance = 7.0
    elif min_nearest_identity < 97.0:
        c_distance = 4.0

    return round(min(100.0, c_novelty + c_phylo + c_sparse + c_conflict + c_distance), 2)


def _has_taxonomy_conflict(row: dict) -> bool:
    """Return True if alt-ref databases disagree with the primary taxonomy."""
    primary_tax = str(row.get('taxonomy') or '').lower().strip()
    if not primary_tax or primary_tax == 'na':
        return False
    for key, val in row.items():
        if key.startswith('alt_tax_') and val not in (None, 'NA', 'na', ''):
            alt_tax = str(val).lower().strip()
            if alt_tax and alt_tax != 'na' and alt_tax != primary_tax:
                return True
    return False


def _safe_float(val, default: float = 0.0) -> float:
    try:
        return float(val)
    except (TypeError, ValueError):
        return default


def _safe_int(val, default: int = 0) -> int:
    try:
        return int(val)
    except (TypeError, ValueError):
        return default


# ─────────────────────────────────────────────────────────────────────────────
# Cluster aggregation
# ─────────────────────────────────────────────────────────────────────────────

def _taxonomy_consensus(taxonomies: List[str]) -> Tuple[str, float, str]:
    """Return (consensus_taxonomy, agreement_fraction, mode).

    mode: 'unanimous' | 'majority' | 'mixed'
    """
    clean = [t for t in taxonomies if t and t.lower() not in ('na', 'none', '')]
    if not clean:
        return 'NA', 0.0, 'mixed'
    counts: Dict[str, int] = {}
    for t in clean:
        counts[t] = counts.get(t, 0) + 1
    best = max(counts, key=counts.__getitem__)
    frac = counts[best] / len(clean)
    if frac == 1.0:
        mode = 'unanimous'
    elif frac >= 0.5:
        mode = 'majority'
    else:
        mode = 'mixed'
    return best, round(frac, 3), mode


def _collect_unique_flags(rows: List[dict]) -> str:
    seen = set()
    for row in rows:
        flags = str(row.get('placement_flags') or '')
        for f in flags.split(';'):
            f = f.strip()
            if f:
                seen.add(f)
    return ';'.join(sorted(seen))


def aggregate_cluster_rows(
    assessment_rows: List[dict],
    phylo_isolation: Optional[Dict[str, float]] = None,
) -> List[dict]:
    """Group *assessment_rows* by cluster and return one summary dict per cluster.

    Singletons (ClusterSize=1, ClusterRepresentative='self') are included as
    single-member clusters so every sequence appears in cluster_summary.tsv.

    Each returned dict has the columns written by ``write_cluster_summary_tsv``.
    """
    # ── Group rows by representative ─────────────────────────────────────────
    iso: Dict[str, float] = phylo_isolation or {}
    clusters: Dict[str, List[dict]] = {}

    for row in assessment_rows:
        rep = row.get('cluster_representative', 'self')
        # Sequences that ARE the representative (cluster_rep='self') use their
        # own id as the representative key.
        if rep in ('self', 'N/A', 'duplicate', None, ''):
            rep = row['id']
        clusters.setdefault(rep, []).append(row)

    # ── Aggregate per cluster ─────────────────────────────────────────────────
    cluster_rows = []
    for rep_id, members in clusters.items():
        # The representative row is the one whose id matches rep_id, or first
        rep_row = next((r for r in members if r['id'] == rep_id), members[0])

        all_ids: List[str] = [str(r['id']) for r in members]
        member_ids_not_rep = [i for i in all_ids if i != rep_id]

        # Novelty aggregation
        novelty_scores = [_safe_float(r.get('novelty_score')) for r in members]
        nearest_ids = [_safe_float(r.get('nearest_identity'), 100.0) for r in members]
        mean_novelty = round(statistics.mean(novelty_scores), 2) if novelty_scores else 0.0
        max_novelty = round(max(novelty_scores), 2) if novelty_scores else 0.0
        min_nearest = round(min(nearest_ids), 2) if nearest_ids else 100.0
        mean_nearest = round(statistics.mean(nearest_ids), 2) if nearest_ids else 100.0

        # Taxonomy consensus
        tax_list = [str(r.get('taxonomy') or '') for r in members]
        consensus_tax, tax_agreement, tax_mode = _taxonomy_consensus(tax_list)

        # Sequencing priority: highest in cluster wins
        _prio_rank = {'HIGH': 3, 'MEDIUM': 2, 'LOW': 1, 'NA': 0}
        priorities = [str(r.get('sequencing_priority', 'NA')) for r in members]
        best_priority = max(priorities, key=lambda p: _prio_rank.get(p, 0))

        # Placement flags (union of all unique flags)
        unified_flags = _collect_unique_flags(members)

        # Phylogenetic isolation scores
        iso_scores = [iso.get(str(r['id']), 0.0) for r in members]
        rep_iso = iso.get(str(rep_id), 0.0)
        mean_iso = round(statistics.mean(iso_scores), 6) if iso_scores else 0.0

        # Taxonomy conflict (any alt-ref disagrees with primary?)
        any_conflict = any(_has_taxonomy_conflict(r) for r in members)

        # Density aggregates
        mean_ge97 = statistics.mean([_safe_int(r.get('matches_ge_97')) for r in members])
        mean_ge99 = statistics.mean([_safe_int(r.get('matches_ge_99')) for r in members])
        min_ge97 = min(_safe_int(r.get('matches_ge_97')) for r in members)
        min_ge99 = min(_safe_int(r.get('matches_ge_99')) for r in members)

        # Composite investigation scores — computed for every member
        scored_members: List[Tuple[float, dict]] = []
        for r in members:
            r_iso = iso.get(str(r['id']), 0.0)
            inv = compute_investigation_composite(
                novelty_score=_safe_float(r.get('novelty_score')),
                phylo_isolation=r_iso,
                matches_ge_97=_safe_int(r.get('matches_ge_97')),
                matches_ge_99=_safe_int(r.get('matches_ge_99')),
                taxonomy_conflict=_has_taxonomy_conflict(r),
                min_nearest_identity=_safe_float(r.get('nearest_identity'), 100.0),
            )
            scored_members.append((inv, r))
        scored_members.sort(key=lambda x: x[0], reverse=True)
        inv_scores = [s for s, _ in scored_members]
        cluster_investigation = round(inv_scores[0], 2) if inv_scores else 0.0

        # Backup IDs: members ordered by InvestigationScore, excluding the top one
        backup_ids = [r['id'] for _, r in scored_members[1:]]

        # InTree status for the representative
        rep_in_tree = rep_row.get('in_tree', 'Unknown')

        cluster_rows.append({
            'ClusterID': f"CLUSTER_{rep_id}",
            'RepresentativeID': rep_id,
            'ClusterSize': len(members),
            'TaxonomyConsensus': consensus_tax,
            'TaxonomyAgreement': tax_agreement,
            'TaxonomyConsensusMode': tax_mode,
            'MultiDBConflict': 'Yes' if any_conflict else ('N/A' if not any(r for r in members if any(k.startswith('alt_tax_') for k in r)) else 'No'),
            'MeanNoveltyScore': mean_novelty,
            'MaxNoveltyScore': max_novelty,
            'MinNearestIdentity': min_nearest,
            'MeanNearestIdentity': mean_nearest,
            'BestSequencingPriority': best_priority,
            'RepPhyloIsolation': rep_iso,
            'MeanPhyloIsolation': mean_iso,
            'MeanMatchesGE97': round(mean_ge97, 1),
            'MinMatchesGE97': min_ge97,
            'MeanMatchesGE99': round(mean_ge99, 1),
            'MinMatchesGE99': min_ge99,
            'PlacementFlags': unified_flags,
            'RepInTree': rep_in_tree,
            'ClusterInvestigationScore': cluster_investigation,
            'NBackupsAvailable': len(backup_ids),
            'Backup1ID': backup_ids[0] if len(backup_ids) > 0 else 'NA',
            'Backup2ID': backup_ids[1] if len(backup_ids) > 1 else 'NA',
            'Backup3ID': backup_ids[2] if len(backup_ids) > 2 else 'NA',
            # member list (not rep for singletons)
            'MemberIDs': ';'.join(sorted(member_ids_not_rep)),
            # stash scored_members for backup_candidates writer (popped later)
            '_scored_members': scored_members,
        })

    # Sort by investigation score descending so the most interesting are first
    cluster_rows.sort(key=lambda r: r['ClusterInvestigationScore'], reverse=True)
    return cluster_rows


# ─────────────────────────────────────────────────────────────────────────────
# Writers
# ─────────────────────────────────────────────────────────────────────────────

_CLUSTER_SUMMARY_COLUMNS = [
    'ClusterID',
    'RepresentativeID',
    'ClusterSize',
    'TaxonomyConsensus',
    'TaxonomyAgreement',
    'TaxonomyConsensusMode',
    'MultiDBConflict',
    'MeanNoveltyScore',
    'MaxNoveltyScore',
    'MinNearestIdentity',
    'MeanNearestIdentity',
    'BestSequencingPriority',
    'RepPhyloIsolation',
    'MeanPhyloIsolation',
    'MeanMatchesGE97',
    'MinMatchesGE97',
    'MeanMatchesGE99',
    'MinMatchesGE99',
    'PlacementFlags',
    'RepInTree',
    'ClusterInvestigationScore',
    'NBackupsAvailable',
    'Backup1ID',
    'Backup2ID',
    'Backup3ID',
    'MemberIDs',
]


def write_cluster_summary_tsv(path: str | Path, cluster_rows: List[dict]) -> str:
    """Write cluster_summary.tsv — one row per cluster, sorted by InvestigationScore."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, 'w', newline='') as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=_CLUSTER_SUMMARY_COLUMNS,
            delimiter='\t',
            extrasaction='ignore',
            restval='NA',
        )
        writer.writeheader()
        writer.writerows(cluster_rows)
    return str(p)


# ─────────────────────────────────────────────────────────────────────────────
# Backup candidate helpers
# ─────────────────────────────────────────────────────────────────────────────

def _backup_rationale(row: dict, rank: int) -> str:
    """Build a short human-readable string explaining why this sequence is a backup.

    Examples:
        "Rank2; InvestigationScore=74.5; NoveltyScore=82.1; InTree=True; Crowding=sparse"
    """
    parts = [
        f"Rank{rank}",
        f"InvestigationScore={row.get('investigation_score', 'NA')}",
        f"NoveltyScore={row.get('novelty_score', 'NA')}",
        f"NearestIdentity={row.get('nearest_identity', 'NA')}",
        f"InTree={row.get('in_tree', 'Unknown')}",
        f"Crowding={row.get('crowding', 'NA')}",
    ]
    flags = str(row.get('placement_flags') or '').strip()
    if flags:
        parts.append(f"Flags={flags}")
    return '; '.join(parts)


def write_backup_candidates_tsv(
    path: str | Path,
    cluster_rows: List[dict],
    phylo_isolation: Optional[Dict[str, float]] = None,
    n_backups: int = 5,
) -> str:
    """Write ``backup_candidates.tsv`` — one row per *primary* sequence with
    ordered backups listed as additional columns.

    The primary sequence is the top-ranked member of each cluster
    (highest InvestigationScore).  For singletons the primary has no
    backups; for multi-member clusters the remaining members are ranked
    Backup_1 … Backup_N.

    Parameters
    ----------
    cluster_rows  : Output of ``aggregate_cluster_rows`` (contains the
                    ``_scored_members`` stash).
    phylo_isolation : {seq_id: normalised_branch_length} mapping.
    n_backups     : Maximum number of backup columns to write (default 5).
    """
    iso = phylo_isolation or {}
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)

    # Build fieldnames
    base_fields = [
        'PrimaryID',
        'ClusterID',
        'ClusterSize',
        'PrimaryInvestigationScore',
        'PrimaryNoveltyScore',
        'PrimaryNearestIdentity',
        'PrimaryPhyloIsolation',
        'PrimaryTaxonomy',
        'PrimaryClassificationHit',
        'PrimarySequencingPriority',
        'PrimaryInTree',
        'PrimaryCrowding',
        'PrimaryPlacementFlags',
        'NBackupsAvailable',
    ]
    for i in range(1, n_backups + 1):
        base_fields += [
            f'Backup{i}_ID',
            f'Backup{i}_InvestigationScore',
            f'Backup{i}_NoveltyScore',
            f'Backup{i}_NearestIdentity',
            f'Backup{i}_PhyloIsolation',
            f'Backup{i}_InTree',
            f'Backup{i}_Crowding',
            f'Backup{i}_Rationale',
        ]

    rows_out = []
    for cr in cluster_rows:
        scored = cr.get('_scored_members', [])
        if not scored:
            continue
        primary_score, primary_row = scored[0]
        backups = scored[1:]

        primary_iso = float(iso.get(str(primary_row.get('id', '')), 0.0))

        out: dict = {
            'PrimaryID': primary_row.get('id', 'NA'),
            'ClusterID': cr['ClusterID'],
            'ClusterSize': cr['ClusterSize'],
            'PrimaryInvestigationScore': round(primary_score, 2),
            'PrimaryNoveltyScore': primary_row.get('novelty_score', 'NA'),
            'PrimaryNearestIdentity': primary_row.get('nearest_identity', 'NA'),
            'PrimaryPhyloIsolation': round(primary_iso, 6),
            'PrimaryTaxonomy': primary_row.get('taxonomy', 'NA'),
            'PrimaryClassificationHit': primary_row.get('classification_hit',
                                          primary_row.get('best_hit', 'NA')),
            'PrimarySequencingPriority': primary_row.get('sequencing_priority', 'NA'),
            'PrimaryInTree': primary_row.get('in_tree', 'Unknown'),
            'PrimaryCrowding': primary_row.get('crowding', 'NA'),
            'PrimaryPlacementFlags': primary_row.get('placement_flags', ''),
            'NBackupsAvailable': len(backups),
        }

        for i, (bscore, brow) in enumerate(backups[:n_backups], start=1):
            b_iso = float(iso.get(str(brow.get('id', '')), 0.0))
            # Attach investigation_score temporarily for rationale
            brow_with_score = dict(brow)
            brow_with_score['investigation_score'] = round(bscore, 2)
            out[f'Backup{i}_ID'] = brow.get('id', 'NA')
            out[f'Backup{i}_InvestigationScore'] = round(bscore, 2)
            out[f'Backup{i}_NoveltyScore'] = brow.get('novelty_score', 'NA')
            out[f'Backup{i}_NearestIdentity'] = brow.get('nearest_identity', 'NA')
            out[f'Backup{i}_PhyloIsolation'] = round(b_iso, 6)
            out[f'Backup{i}_InTree'] = brow.get('in_tree', 'Unknown')
            out[f'Backup{i}_Crowding'] = brow.get('crowding', 'NA')
            out[f'Backup{i}_Rationale'] = _backup_rationale(brow_with_score, i)

        # Fill missing backup slots with NA
        for i in range(len(backups) + 1, n_backups + 1):
            for suffix in ('ID', 'InvestigationScore', 'NoveltyScore',
                           'NearestIdentity', 'PhyloIsolation', 'InTree',
                           'Crowding', 'Rationale'):
                out[f'Backup{i}_{suffix}'] = 'NA'

        rows_out.append(out)

    # Sort by PrimaryInvestigationScore descending so top candidates are first
    rows_out.sort(key=lambda r: _safe_float(r.get('PrimaryInvestigationScore')), reverse=True)

    with open(p, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=base_fields, extrasaction='ignore',
                                restval='NA', delimiter='\t')
        writer.writeheader()
        writer.writerows(rows_out)

    return str(p)


def write_per_cluster_csvs(
    outdir: str | Path,
    assessment_rows: List[dict],
    phylo_isolation: Optional[Dict[str, float]] = None,
    cluster_rows: Optional[List[dict]] = None,
) -> List[str]:
    """Write one CSV per cluster to ``outdir/clusters/<rep_id>.csv``.

    Each CSV contains every sequence that belongs to the cluster (including the
    representative), with all assessment columns plus:
      - ``IsRepresentative``  (True/False)
      - ``PhyloIsolation``    (normalised leaf branch length)
      - ``InvestigationScore``(composite per-sequence score)
      - ``BackupRank``        (0 = primary candidate, 1 = first backup, 2 = …)
      - ``PrimaryID``         (ID of the primary sequence for this cluster)
      - ``BackupRationale``   (human-readable reason for this rank)
    """
    clusters_dir = Path(outdir) / 'clusters'
    clusters_dir.mkdir(parents=True, exist_ok=True)

    iso = phylo_isolation or {}

    # Build a lookup of {rep_id: scored_members} from cluster_rows stash
    scored_lookup: Dict[str, List[Tuple[float, dict]]] = {}
    if cluster_rows:
        for cr in cluster_rows:
            rid = cr.get('RepresentativeID', '')
            sm: Optional[List[Tuple[float, dict]]] = cr.get('_scored_members')  # type: ignore[assignment]
            if rid and sm is not None:
                scored_lookup[rid] = sm

    # Build representative → [member rows] map
    cluster_map: Dict[str, List[dict]] = {}
    for row in assessment_rows:
        rep = row.get('cluster_representative', 'self')
        if rep in ('self', 'N/A', 'duplicate', None, ''):
            rep = row['id']
        cluster_map.setdefault(rep, []).append(row)

    written_paths = []

    # Discover alt-db column keys from the first row
    alt_keys: List[str] = []
    if assessment_rows:
        alt_keys = sorted(k for k in assessment_rows[0] if k.startswith('alt_tax_'))

    # Build fieldnames dynamically
    base_fields = [
        'ID', 'IsRepresentative', 'BackupRank', 'PrimaryID', 'BackupRationale',
        'Taxonomy', 'ClassificationHit', 'ClassificationIdentity', 'ClassificationConfidence',
        'NearestHit', 'NearestIdentity',
        'MatchesGE99', 'MatchesGE97', 'MatchesGE95',
        'NoveltyScore', 'Crowding', 'SequencingPriority',
        'PhyloIsolation', 'InvestigationScore',
        'InTree', 'ClusterRepresentative', 'ClusterSize', 'ClusteredMembers',
        'PlacementFlags',
    ]
    # Add alt-db column groups
    alt_fields: List[str] = []
    for key in alt_keys:
        safe = key[len('alt_tax_'):]
        alt_fields += [
            f'Taxonomy_{safe}',
            f'ClassificationHit_{safe}',
            f'Identity_{safe}',
            f'Confidence_{safe}',
        ]
    fieldnames = base_fields + alt_fields

    for rep_id, member_rows in cluster_map.items():
        # Safe filename: strip non-alphanumeric characters
        safe_name = re.sub(r'[^A-Za-z0-9_.\-]+', '_', rep_id)[:80]
        out_path = clusters_dir / f'{safe_name}.csv'

        # Determine backup ranking for this cluster
        # Use pre-computed scored_members if available, else score on the fly
        sm: List[Tuple[float, dict]] = scored_lookup.get(rep_id, [])
        if not sm:
            for r in member_rows:
                r_iso = float(iso.get(str(r['id']), 0.0))
                inv = compute_investigation_composite(
                    novelty_score=_safe_float(r.get('novelty_score')),
                    phylo_isolation=r_iso,
                    matches_ge_97=_safe_int(r.get('matches_ge_97')),
                    matches_ge_99=_safe_int(r.get('matches_ge_99')),
                    taxonomy_conflict=_has_taxonomy_conflict(r),
                    min_nearest_identity=_safe_float(r.get('nearest_identity'), 100.0),
                )
                sm.append((inv, r))
            sm.sort(key=lambda x: x[0], reverse=True)

        # Build rank lookup {seq_id: rank}  (0 = primary)
        rank_map: Dict[str, int] = {r['id']: i for i, (_, r) in enumerate(sm)}
        primary_id: str = sm[0][1]['id'] if sm else rep_id

        with open(out_path, 'w', newline='') as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction='ignore', restval='NA')
            writer.writeheader()

            # Write rows in backup-rank order (primary first, then backups)
            sorted_rows = sorted(member_rows, key=lambda r: rank_map.get(r['id'], 999))

            for row in sorted_rows:
                seq_id = str(row['id'])
                r_iso: float = float(iso.get(seq_id, 0.0))
                brank = rank_map.get(seq_id, 0)
                inv_score = compute_investigation_composite(
                    novelty_score=_safe_float(row.get('novelty_score')),
                    phylo_isolation=_safe_float(r_iso),
                    matches_ge_97=_safe_int(row.get('matches_ge_97')),
                    matches_ge_99=_safe_int(row.get('matches_ge_99')),
                    taxonomy_conflict=_has_taxonomy_conflict(row),
                    min_nearest_identity=_safe_float(row.get('nearest_identity'), 100.0),
                )
                row_with_score: Dict[str, object] = dict(row)
                row_with_score['investigation_score'] = round(inv_score, 2)
                rationale = (
                    'Primary candidate (highest InvestigationScore in cluster)'
                    if brank == 0
                    else _backup_rationale(row_with_score, brank)
                )

                # Flatten alt-db columns to human-friendly names
                flat_alt: dict = {}
                for key in alt_keys:
                    safe = key[len('alt_tax_'):]
                    flat_alt[f'Taxonomy_{safe}'] = row.get(f'alt_tax_{safe}', 'NA')
                    flat_alt[f'ClassificationHit_{safe}'] = row.get(f'alt_class_hit_{safe}', 'NA')
                    flat_alt[f'Identity_{safe}'] = row.get(f'alt_ident_{safe}', 'NA')
                    flat_alt[f'Confidence_{safe}'] = row.get(f'alt_conf_{safe}', 'NA')

                out_row = {
                    'ID': seq_id,
                    'IsRepresentative': seq_id == rep_id,
                    'BackupRank': brank,
                    'PrimaryID': primary_id,
                    'BackupRationale': rationale,
                    'Taxonomy': row.get('taxonomy', 'NA'),
                    'ClassificationHit': row.get('classification_hit', row.get('best_hit', 'NA')),
                    'ClassificationIdentity': row.get('classification_identity', 'NA'),
                    'ClassificationConfidence': row.get('classification_confidence', 'NA'),
                    'NearestHit': row.get('nearest_hit', 'NA'),
                    'NearestIdentity': row.get('nearest_identity', 'NA'),
                    'MatchesGE99': row.get('matches_ge_99', 'NA'),
                    'MatchesGE97': row.get('matches_ge_97', 'NA'),
                    'MatchesGE95': row.get('matches_ge_95', 'NA'),
                    'NoveltyScore': row.get('novelty_score', 'NA'),
                    'Crowding': row.get('crowding', 'NA'),
                    'SequencingPriority': row.get('sequencing_priority', 'NA'),
                    'PhyloIsolation': round(r_iso, 6),
                    'InvestigationScore': inv_score,
                    'InTree': row.get('in_tree', 'Unknown'),
                    'ClusterRepresentative': row.get('cluster_representative', 'N/A'),
                    'ClusterSize': row.get('cluster_size', '1'),
                    'ClusteredMembers': row.get('clustered_members', ''),
                    'PlacementFlags': row.get('placement_flags', ''),
                    **flat_alt,
                }
                writer.writerow(out_row)

        written_paths.append(str(out_path))

    return written_paths


# ─────────────────────────────────────────────────────────────────────────────
# High-level entry point
# ─────────────────────────────────────────────────────────────────────────────

def generate_cluster_reports(
    outdir: str | Path,
    assessment_rows: List[dict],
    tree_path: Optional[str] = None,
) -> Tuple[Optional[str], List[str], Optional[str]]:
    """Generate all cluster-level output files.

    Parameters
    ----------
    outdir          : Run output directory.
    assessment_rows : Per-sequence rows from ``build_sequence_assessment_rows``.
    tree_path       : Path to ``current_tree.nwk``; used for phylo isolation.

    Returns
    -------
    (cluster_summary_path, [per_cluster_csv_paths], backup_candidates_path)
      Any entry may be None/empty if *assessment_rows* is empty.
    """
    if not assessment_rows:
        return None, [], None

    import logging as _logging
    log = _logging.getLogger(__name__)

    # Phylogenetic isolation scores
    phylo_isolation: Dict[str, float] = {}
    if tree_path and Path(tree_path).exists():
        phylo_isolation = compute_phylogenetic_isolation(tree_path)
        log.info(
            "[CLUSTER] Computed phylogenetic isolation for %d tree leaves from %s",
            len(phylo_isolation), tree_path,
        )
    else:
        log.info(
            "[CLUSTER] No tree file available at %s; phylogenetic isolation will be 0",
            tree_path,
        )

    # Attach phylo isolation and composite score to every assessment row (in-place)
    for row in assessment_rows:
        p_iso = float(phylo_isolation.get(str(row['id']), 0.0))
        row['phylo_isolation'] = round(p_iso, 6)
        row['investigation_score'] = compute_investigation_composite(
            novelty_score=_safe_float(row.get('novelty_score')),
            phylo_isolation=p_iso,
            matches_ge_97=_safe_int(row.get('matches_ge_97')),
            matches_ge_99=_safe_int(row.get('matches_ge_99')),
            taxonomy_conflict=_has_taxonomy_conflict(row),
            min_nearest_identity=_safe_float(row.get('nearest_identity'), 100.0),
        )

    # Cluster aggregation (also builds _scored_members stash per cluster)
    cluster_rows = aggregate_cluster_rows(assessment_rows, phylo_isolation=phylo_isolation)

    summary_path = write_cluster_summary_tsv(
        Path(outdir) / 'cluster_summary.tsv', cluster_rows
    )
    log.info(
        "[CLUSTER] Wrote cluster summary (%d clusters) to %s", len(cluster_rows), summary_path
    )

    # Backup candidates table
    backup_path = write_backup_candidates_tsv(
        Path(outdir) / 'backup_candidates.tsv',
        cluster_rows,
        phylo_isolation=phylo_isolation,
    )
    log.info("[CLUSTER] Wrote backup candidates table to %s", backup_path)

    per_cluster_paths = write_per_cluster_csvs(
        outdir, assessment_rows, phylo_isolation, cluster_rows=cluster_rows
    )
    log.info(
        "[CLUSTER] Wrote %d per-cluster CSV files to %s/clusters/",
        len(per_cluster_paths), outdir,
    )

    return summary_path, per_cluster_paths, backup_path


# ─────────────────────────────────────────────────────────────────────────────
# Suggested additional analyses (comments for future work / explanation)
# ─────────────────────────────────────────────────────────────────────────────
"""
ADDITIONAL ANALYSES THAT COULD BE ADDED
========================================

1. **Sister-clade richness**
   Count the number of sequences in the immediate sister clade to each
   sequence in the tree.  A sequence whose sister clade contains many known
   sequences is well-contextualised; one whose sister clade is empty or sparse
   deserves further investigation.
   Implementation: parse tree, walk up to LCA, count leaves.

2. **Multi-run novelty trajectory**
   Track how NoveltyScore for each sequence changes across successive `phyloselect run`
   calls.  A sequence that remains novel even after adding more data is a stronger
   candidate than one that becomes crowded quickly.
   Implementation: store novelty_score per run in a new DB table `novelty_history`.

3. **Taxonomy confidence trajectory at different rank levels**
   Compare classification confidence at species vs genus vs family.  A sequence
   with high family-level confidence but low species-level confidence is at an
   interesting boundary.
   Implementation: parse GTDB/SILVA taxonomy strings and compute rank-specific
   confidence ladder.

4. **Functional annotation context**
   If a functional annotation file is provided (via --functional), annotate each
   cluster with the dominant predicted function.  Clusters with high novelty AND
   unique functional profiles are the top candidates.
   Implementation: join functional annotation file to cluster members.

5. **16S V-region coverage mapping**
   Map each sequence against V1-V9 region primers to determine which hypervariable
   regions are present.  Sequences covering V3-V5 alone vs V1-V9 have very
   different resolution and should not be directly compared.
   Implementation: vsearch against a short V-region primer/probe database.

6. **Per-cluster MSA and local divergence score**
   For clusters of ≥3 sequences, align members against each other and compute
   the column-wise entropy.  High within-cluster entropy at certain positions may
   indicate sequencing errors vs genuine sequence heterogeneity.
   Implementation: MAFFT on cluster members, Shannon entropy per column.

7. **Redundancy rank within cluster**
   Within each cluster, rank members by how many unique sequence positions they
   contribute relative to the representative.  This guides which members to
   sequence further or which are truly redundant.
   Implementation: pairwise alignment to cluster representative.
"""








