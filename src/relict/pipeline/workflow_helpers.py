from __future__ import annotations

import gzip
import os
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Set, Tuple

from relict.taxonomy import parse_reference_header_taxonomy, taxonomy_matches_kingdom
from relict.utils.fasta import read_fasta


ClassificationRow = Dict[str, object]


def _assignment_source_is_fasta(assignments_path: str, source_fasta_path: Optional[str] = None) -> bool:
    try:
        if source_fasta_path and os.path.samefile(assignments_path, source_fasta_path):
            return True
    except Exception:
        pass

    try:
        if str(assignments_path).endswith('.gz'):
            handle_ctx = gzip.open(assignments_path, 'rt')
        else:
            handle_ctx = open(assignments_path, 'rt')
        with handle_ctx as handle:
            for line in handle:
                stripped = line.strip()
                if not stripped:
                    continue
                return stripped.startswith('>')
    except Exception:
        return False
    return False


def build_orig_to_short_map(alias_entries, db) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    if not alias_entries:
        return mapping
    for short, orig in alias_entries:
        mapping[orig] = short
        mapping[short] = short
        try:
            cid = db._canonical_from_header(orig)
            if cid:
                mapping[cid] = short
        except Exception:
            pass
    return mapping


def resolve_short_id(qid: str, orig_to_short: Dict[str, str], db) -> Optional[str]:
    if qid in orig_to_short:
        return orig_to_short[qid]
    try:
        cid = db._canonical_from_header(qid)
    except Exception:
        cid = None
    if isinstance(cid, str) and cid in orig_to_short:
        return orig_to_short[cid]
    try:
        lk = qid.lower()
        for key, value in orig_to_short.items():
            if isinstance(key, str) and key.lower() == lk:
                return value
    except Exception:
        pass
    return None


def iter_classification_rows(classification_tsv: str) -> Iterator[ClassificationRow]:
    with open(classification_tsv) as handle:
        next(handle, None)
        for line in handle:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 1 or not parts[0]:
                continue
            try:
                identity = float(parts[2]) if len(parts) > 2 and parts[2] not in ('', 'NA') else None
            except Exception:
                identity = None
            try:
                confidence = float(parts[4]) if len(parts) > 4 and parts[4] not in ('', 'NA') else None
            except Exception:
                confidence = None
            yield {
                'qid': parts[0],
                'best': parts[1] if len(parts) > 1 and parts[1] != '' else None,
                'identity': identity,
                'tax': parts[3] if len(parts) > 3 and parts[3] != '' else None,
                'confidence': confidence,
            }


def iter_assignment_rows(assignments_path: str, source_fasta_path: Optional[str] = None) -> Iterator[dict[str, object]]:
    """Yield normalized taxonomy-assignment rows from TSV or FASTA headers.

    If `assignments_path` is the same file as `source_fasta_path`, taxonomy is
    parsed from GTDB-style FASTA headers.
    """
    use_fasta_headers = _assignment_source_is_fasta(assignments_path, source_fasta_path=source_fasta_path)

    if use_fasta_headers:
        for header, _ in read_fasta(assignments_path):
            _, tax = parse_reference_header_taxonomy(header)
            yield {
                'qid': header,
                'tax': tax,
                'confidence': None,
            }
        return

    with open(assignments_path) as tf:
        for line in tf:
            if not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            if not parts:
                continue
            try:
                conf = float(parts[2]) if len(parts) > 2 and parts[2] not in ('', 'NA') else None
            except Exception:
                conf = None
            yield {
                'qid': parts[0],
                'tax': parts[1] if len(parts) > 1 else None,
                'confidence': conf,
            }


def load_taxonomy_entries_from_assignments(taxa_file: str, orig_to_short: Dict[str, str], db, dataset: str, source_fasta_path: Optional[str] = None):
    entries = []
    for row in iter_assignment_rows(taxa_file, source_fasta_path=source_fasta_path):
        mapped = resolve_short_id(str(row['qid']), orig_to_short, db)
        if not mapped:
            continue
        entries.append((mapped, row.get('tax'), row.get('confidence'), dataset))
    return entries


def load_classification_results_for_dataset(classification_tsv: str, orig_to_short: Dict[str, str], db, dataset: str):
    tax_entries = []
    dist_entries = []
    for row in iter_classification_rows(classification_tsv):
        mapped = resolve_short_id(str(row['qid']), orig_to_short, db)
        if not mapped:
            continue
        tax_entries.append((mapped, row['tax'], row['confidence'], dataset))
        dist_entries.append((mapped, dataset, row['best'], row['identity']))
    return tax_entries, dist_entries


def classification_ids_matching_kingdom(classification_tsv: str, kingdom: str):
    matched = set()
    if not classification_tsv or not kingdom:
        return matched
    for row in iter_classification_rows(classification_tsv):
        raw_tax = row.get('tax')
        tax = str(raw_tax) if raw_tax is not None else None
        if taxonomy_matches_kingdom(tax, kingdom):
            matched.add(str(row['qid']))
    return matched


def prune_dataset_by_kingdom(db, dataset: str, kingdom: str | None):
    if not kingdom:
        return 0
    to_delete = []
    with db.connect() as conn:
        cur = conn.cursor()
        cur.execute(
            "SELECT s.id, t.taxonomy FROM sequences s LEFT JOIN taxonomy t ON s.id = t.id WHERE s.dataset = ?",
            (dataset,),
        )
        for iid, tax in cur.fetchall():
            if tax is None:
                continue
            tax_text = str(tax).strip()
            if not tax_text or tax_text.lower() in ('na', 'none'):
                continue
            if not taxonomy_matches_kingdom(tax_text, kingdom):
                to_delete.append(iid)
        for iid in to_delete:
            cur.execute('DELETE FROM distances WHERE id = ?', (iid,))
            cur.execute('DELETE FROM taxonomy WHERE id = ?', (iid,))
            cur.execute('DELETE FROM colors WHERE id = ?', (iid,))
            cur.execute('DELETE FROM sequences WHERE id = ?', (iid,))
        conn.commit()
    return len(to_delete)


def read_combined_taxonomy_ids(path: str | Path):
    ids = []
    p = Path(path)
    if not p.exists():
        return ids
    with open(p) as handle:
        next(handle, None)
        for line in handle:
            iid = line.strip().split('\t')[0]
            if iid:
                ids.append(iid)
    return ids


def collect_db_taxonomy_rows(db, ids: Optional[Sequence[str]] = None):
    with db.connect() as conn:
        cur = conn.cursor()
        if ids:
            placeholders = ','.join('?' for _ in ids)
            cur.execute(
                f"SELECT s.id, t.taxonomy, t.confidence FROM sequences s LEFT JOIN taxonomy t ON s.id = t.id WHERE s.id IN ({placeholders})",
                tuple(ids),
            )
        else:
            cur.execute("SELECT s.id, t.taxonomy, t.confidence FROM sequences s LEFT JOIN taxonomy t ON s.id = t.id")
        return cur.fetchall()


def merge_combined_taxonomy_rows(base_rows, classification_tsv: str, orig_to_short: Dict[str, str], db):
    merged = {}
    order: List[str] = []
    for rid, tax, conf in base_rows:
        merged[rid] = (tax if tax is not None else 'NA', conf if conf is not None else 'NA')
        order.append(rid)
    if classification_tsv:
        for row in iter_classification_rows(classification_tsv):
            mapped = resolve_short_id(str(row['qid']), orig_to_short, db)
            if not mapped:
                continue
            if mapped not in merged:
                order.append(mapped)
            merged[mapped] = (
                row['tax'] if row['tax'] is not None else 'NA',
                row['confidence'] if row['confidence'] is not None else 'NA',
            )
    return [(iid, merged[iid][0], merged[iid][1]) for iid in order]


def write_combined_taxonomy_tsv(path: str | Path, rows: Iterable[Tuple[str, object, object]]):
    p = Path(path)
    with open(p, 'w') as outf:
        outf.write('ID\tTaxon\tConfidence\n')
        for iid, tax, conf in rows:
            outf.write(f"{iid}\t{tax if tax is not None else 'NA'}\t{conf if conf is not None else 'NA'}\n")
    return str(p)


def iter_novelty_rows(novelty_tsv: str):
    with open(novelty_tsv) as handle:
        next(handle, None)
        for line in handle:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4:
                continue
            try:
                nearest_identity = float(parts[1]) if parts[1] not in ('', 'NA') else None
            except Exception:
                nearest_identity = None
            yield {
                'qid': parts[0],
                'nearest_identity': nearest_identity,
                'nearest_hit': parts[2] if len(parts) > 2 else None,
                'novel': parts[3] if len(parts) > 3 else None,
            }


def build_placement_warning_rows(classification_tsv: str, novelty_tsv: Optional[str], orig_to_short: Dict[str, str], db):
    novelty_by_qid = {}
    if novelty_tsv:
        try:
            for row in iter_novelty_rows(novelty_tsv):
                novelty_by_qid[str(row['qid'])] = row
        except Exception:
            novelty_by_qid = {}

    warning_rows = []
    for row in iter_classification_rows(classification_tsv):
        qid = str(row['qid'])
        mapped = resolve_short_id(qid, orig_to_short, db) or qid
        flags = []
        best = row.get('best')
        identity = row.get('identity')
        confidence = row.get('confidence')
        tax = row.get('tax')
        nov = novelty_by_qid.get(qid, {})
        nearest_identity = nov.get('nearest_identity')

        if not best or best == 'NA' or identity in (None, 0, 0.0):
            flags.append('NO_REFERENCE_HIT')
        if isinstance(identity, (int, float)) and float(identity) < 95.0:
            flags.append('LOW_CLASSIFICATION_IDENTITY')
        if isinstance(confidence, (int, float)) and float(confidence) < 0.8:
            flags.append('LOW_CONFIDENCE')
        if isinstance(nearest_identity, (int, float)) and float(nearest_identity) < 97.0:
            flags.append('LOW_NEAREST_IDENTITY')
            if tax not in (None, '', 'NA'):
                flags.append('NOVEL_BUT_ASSIGNED')
        if nov and nearest_identity in (None, 0, 0.0) and tax not in (None, '', 'NA'):
            flags.append('ASSIGNED_WITHOUT_PLACEMENT_SUPPORT')

        if flags:
            warning_rows.append({
                'id': mapped,
                'query_id': qid,
                'taxonomy': tax if tax is not None else 'NA',
                'classification_identity': identity if identity is not None else 'NA',
                'classification_confidence': confidence if confidence is not None else 'NA',
                'nearest_identity': nearest_identity if nearest_identity is not None else 'NA',
                'nearest_hit': nov.get('nearest_hit', 'NA') if nov else 'NA',
                'flags': ';'.join(flags),
            })
    return warning_rows


def write_placement_warning_tsv(path: str | Path, warning_rows):
    p = Path(path)
    with open(p, 'w') as fh:
        fh.write('ID\tQueryID\tTaxonomy\tClassificationIdentity\tClassificationConfidence\tNearestIdentity\tNearestHit\tFlags\n')
        for row in warning_rows:
            fh.write(
                f"{row['id']}\t{row['query_id']}\t{row['taxonomy']}\t{row['classification_identity']}\t"
                f"{row['classification_confidence']}\t{row['nearest_identity']}\t{row['nearest_hit']}\t{row['flags']}\n"
            )
    return str(p)
    return str(p)


def build_sequence_assessment_rows(
    run_ids,
    class_out: str,
    novelty_metrics_tsv: str,
    placement_warning_rows,
    orig_to_short: Dict[str, str],
    db,
    member_to_rep: Optional[Dict[str, str]] = None,
    rep_to_members: Optional[Dict[str, List[str]]] = None,
    tree_ids: Optional[set] = None,
    alt_taxonomies: Optional[Dict[str, Dict[str, tuple]]] = None,
    alt_ref_dbs: Optional[List[str]] = None,
):
    """Build per-sequence assessment rows.

    Optional clustering info (``member_to_rep``, ``rep_to_members``) is
    included when ``--collapse`` was used.  Tree membership (``tree_ids``)
    records whether a sequence entered the phylogenetic tree or was filtered
    out because a cluster representative was used in its place.

    ``alt_taxonomies`` is an optional dict ``{seq_id: {ref_db: (tax, conf, hit, identity)}}``
    produced when multiple reference databases are used for classification.
    ``alt_ref_dbs`` is the ordered list of alternative ref-db names to include.

    Cluster columns added to each row
    -----------------------------------
    InTree              : 'Yes' / 'No' / 'Unknown'
    ClusterRepresentative : seq-id of the representative for this cluster,
                            or 'self' if this sequence IS in the tree as its own entry,
                            or 'duplicate' if removed as an exact duplicate by dereplication,
                            or 'N/A' if no clustering was done.
    ClusterSize         : number of sequences in this cluster (1 if singleton).
    ClusteredMembers    : semi-colon-separated list of member IDs (empty for singletons).
    """
    class_by_id = {}
    for row in iter_classification_rows(class_out):
        mapped = resolve_short_id(str(row['qid']), orig_to_short, db)
        if mapped:
            class_by_id[mapped] = row

    novelty_by_id = {}
    try:
        with open(novelty_metrics_tsv) as fh:
            next(fh, None)
            for line in fh:
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 10:
                    continue
                novelty_by_id[parts[0]] = {
                    'nearest_identity': parts[1],
                    'nearest_hit': parts[2],
                    'novel': parts[3],
                    'matches_ge_99': parts[4],
                    'matches_ge_97': parts[5],
                    'matches_ge_95': parts[6],
                    'novelty_score': parts[7],
                    'crowding': parts[8],
                    'sequencing_priority': parts[9],
                    'density_source': parts[10] if len(parts) > 10 else 'NA',
                }
    except FileNotFoundError:
        import logging as _logging
        _logging.getLogger(__name__).warning(
            "[ASSESSMENT] Novelty metrics file not found: %s — novelty columns will be NA",
            novelty_metrics_tsv,
        )

    warnings_by_id = {row['id']: row for row in placement_warning_rows}

    # Normalise cluster dicts
    _m2r: Dict[str, str] = member_to_rep or {}
    _r2m: Dict[str, List[str]] = rep_to_members or {}
    _tree: Optional[set] = tree_ids  # may be None (collapse not used)

    assessment_rows = []
    for iid in run_ids:
        c = class_by_id.get(iid, {})
        n = novelty_by_id.get(iid, {})
        w = warnings_by_id.get(iid, {})

        # ── Cluster info ────────────────────────────────────────────────────
        if _m2r or _r2m:
            if iid in _r2m:
                # This sequence IS the cluster representative
                members = _r2m[iid]
                cluster_rep = 'self'
                cluster_size = str(len(members) + 1)
                clustered_members = ';'.join(sorted(members))
            elif iid in _m2r:
                # This sequence was collapsed into a representative
                cluster_rep = _m2r[iid]
                # find total size from rep's member list
                members = _r2m.get(cluster_rep, [])
                cluster_size = str(len(members) + 1)
                clustered_members = ''
            else:
                # Singleton — no cluster
                cluster_rep = 'self'
                cluster_size = '1'
                clustered_members = ''
        else:
            cluster_rep = 'N/A'
            cluster_size = '1'
            clustered_members = ''

        # ── Tree membership ─────────────────────────────────────────────────
        if _tree is not None:
            in_tree = 'Yes' if iid in _tree else 'No'
        elif cluster_rep not in ('N/A', 'self'):
            in_tree = 'No'  # was collapsed; rep is in tree
        else:
            in_tree = 'Yes'  # if no collapse, all sequences enter the tree

        # Resolve None values returned by iter_classification_rows to 'NA'.
        _conf = c.get('confidence')
        _ident = c.get('identity')
        _best = c.get('best')
        _tax = c.get('tax')

        # If a sequence is not in the tree and was not collapsed, it was
        # removed as an exact duplicate during the dereplication step.
        if in_tree == 'No' and cluster_rep == 'self':
            cluster_rep = 'duplicate'

        assessment_rows.append({
            'id': iid,
            # ClassificationHit = best vsearch hit in the REFERENCE DB (GTDB/SILVA) → drives Taxonomy
            'taxonomy': _tax if _tax is not None else 'NA',
            'classification_hit': _best if _best is not None else 'NA',
            'classification_identity': _ident if _ident is not None else 'NA',
            'classification_confidence': _conf if _conf is not None else 'NA',
            # NearestHit = best vsearch hit among PRELOAD sequences already in the DB → drives NoveltyScore
            'nearest_hit': n.get('nearest_hit', 'NA'),
            'nearest_identity': n.get('nearest_identity', 'NA'),
            'matches_ge_99': n.get('matches_ge_99', 'NA'),
            'matches_ge_97': n.get('matches_ge_97', 'NA'),
            'matches_ge_95': n.get('matches_ge_95', 'NA'),
            'novelty_score': n.get('novelty_score', 'NA'),
            'crowding': n.get('crowding', 'NA'),
            'sequencing_priority': n.get('sequencing_priority', 'NA'),
            'placement_flags': w.get('flags', ''),
            'in_tree': in_tree,
            'cluster_representative': cluster_rep,
            'cluster_size': cluster_size,
            'clustered_members': clustered_members,
            # Phylogenetic isolation + composite investigation score
            # (populated by cluster_report.generate_cluster_reports; default 'NA' here)
            'phylo_isolation': 'NA',
            'investigation_score': 'NA',
            # Alt-db taxonomy (keys: 'alt_tax_<refdb>', 'alt_class_hit_<refdb>', etc.)
            **_build_alt_tax_cols(iid, alt_taxonomies or {}, alt_ref_dbs or []),
        })
    return assessment_rows


def _build_alt_tax_cols(
    seq_id: str,
    alt_taxonomies: Dict[str, Dict[str, tuple]],
    alt_ref_dbs: List[str],
) -> dict:
    """Return a flat dict of alt-db taxonomy columns for one sequence.

    Keys use the prefix ``alt_`` so ``write_sequence_assessment_tsv`` can
    auto-discover them:
      alt_tax_<safe>       — taxonomy string from that reference DB
      alt_class_hit_<safe> — best vsearch hit (accession) in that reference DB
      alt_ident_<safe>     — vsearch identity to that hit
      alt_conf_<safe>      — classification confidence score
    """
    import re as _re
    cols: dict = {}
    per_seq = alt_taxonomies.get(seq_id, {})
    for ref_db in alt_ref_dbs:
        safe = _re.sub(r'[^A-Za-z0-9]+', '_', ref_db).strip('_')
        entry = per_seq.get(ref_db, (None, None, None, None))
        tax, conf, hit, ident = entry if len(entry) == 4 else (None, None, None, None)
        cols[f'alt_tax_{safe}'] = tax if tax is not None else 'NA'
        cols[f'alt_class_hit_{safe}'] = hit if hit is not None else 'NA'
        cols[f'alt_ident_{safe}'] = f"{float(ident):.2f}" if ident is not None else 'NA'
        cols[f'alt_conf_{safe}'] = f"{float(conf):.4f}" if conf is not None else 'NA'
    return cols


def write_sequence_assessment_tsv(path: str | Path, rows):
    """Write sequence_assessment.tsv.

    Column groups
    -------------
    Taxonomy / ClassificationHit / ClassificationIdentity / ClassificationConfidence
        — from the primary REFERENCE DATABASE (GTDB, SILVA, etc.).
          ClassificationHit is the specific reference accession vsearch matched.

    NearestHit / NearestIdentity / MatchesGE* / NoveltyScore / Crowding / SequencingPriority
        — from the PRELOAD sequences already stored in the DB.
          NearestHit is the closest sequence YOU have previously submitted.
          These columns drive the novelty assessment.

    Taxonomy_<DB> / ClassificationHit_<DB> / Identity_<DB> / Confidence_<DB>  (repeated per alt-ref DB)
        — same information from each additional reference database supplied via
          --alt-ref / --alt-taxa.  The <DB> suffix is the database name derived
          from the filename or --alt-ref-name.
    """
    import re as _re
    p = Path(path)

    # Discover alt-db column groups from the first row that contains them
    alt_ref_safes: List[str] = []
    for row in rows:
        alt_ref_safes = sorted(
            k[len('alt_tax_'):]
            for k in row
            if k.startswith('alt_tax_')
        )
        break

    with open(p, 'w') as fh:
        # ── Fixed columns ────────────────────────────────────────────────────
        base_header = (
            'ID'
            '\tTaxonomy'                    # primary ref DB  ─┐ classification
            '\tClassificationHit'           # ref accession    │ (reference DB)
            '\tClassificationIdentity'      # vsearch %id      │
            '\tClassificationConfidence'    # score 0–1       ─┘
            '\tNearestHit'                  # preload seq ─┐ novelty
            '\tNearestIdentity'             # %id           │ (preload DB)
            '\tMatchesGE99'                 # density        │
            '\tMatchesGE97'                 #               │
            '\tMatchesGE95'                 #              ─┘
            '\tNoveltyScore'
            '\tCrowding'
            '\tSequencingPriority'
            '\tPhyloIsolation'              # normalised leaf branch length (0–1); higher = more isolated in tree
            '\tInvestigationScore'          # composite score (0–100) combining novelty + phylo isolation + density + taxonomy conflict
            '\tInTree'
            '\tClusterRepresentative'
            '\tClusterSize'
            '\tClusteredMembers'
            '\tPlacementFlags'
        )
        # ── Alt-db columns (one group per additional reference database) ─────
        alt_header = ''
        for safe in alt_ref_safes:
            alt_header += (
                f'\tTaxonomy_{safe}'
                f'\tClassificationHit_{safe}'
                f'\tIdentity_{safe}'
                f'\tConfidence_{safe}'
            )
        fh.write(base_header + alt_header + '\n')

        for row in rows:
            base = (
                f"{row['id']}"
                f"\t{row['taxonomy']}"
                f"\t{row.get('classification_hit', row.get('best_hit', 'NA'))}"
                f"\t{row['classification_identity']}"
                f"\t{row['classification_confidence']}"
                f"\t{row['nearest_hit']}"
                f"\t{row['nearest_identity']}"
                f"\t{row['matches_ge_99']}"
                f"\t{row['matches_ge_97']}"
                f"\t{row['matches_ge_95']}"
                f"\t{row['novelty_score']}"
                f"\t{row['crowding']}"
                f"\t{row['sequencing_priority']}"
                f"\t{row.get('phylo_isolation', 'NA')}"
                f"\t{row.get('investigation_score', 'NA')}"
                f"\t{row.get('in_tree', 'Unknown')}"
                f"\t{row.get('cluster_representative', 'N/A')}"
                f"\t{row.get('cluster_size', '1')}"
                f"\t{row.get('clustered_members', '')}"
                f"\t{row['placement_flags']}"
            )
            alt_part = ''
            for safe in alt_ref_safes:
                alt_part += (
                    f"\t{row.get(f'alt_tax_{safe}', 'NA')}"
                    f"\t{row.get(f'alt_class_hit_{safe}', 'NA')}"
                    f"\t{row.get(f'alt_ident_{safe}', 'NA')}"
                    f"\t{row.get(f'alt_conf_{safe}', 'NA')}"
                )
            fh.write(base + alt_part + '\n')
    return str(p)

