from __future__ import annotations

import gzip
import os
import csv
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Set, Tuple

from branchmanager.taxonomy_io import iter_taxonomy_assignment_rows
from branchmanager.taxonomy import parse_reference_header_taxonomy, taxonomy_matches_kingdom
from branchmanager.utils.fasta import read_fasta


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
    """Yield normalized taxonomy-assignment rows from table files or FASTA headers.

    If `assignments_path` is the same file as `source_fasta_path`, taxonomy is
    parsed from GTDB-style FASTA headers. Otherwise TSV/CSV and .gz variants are
    parsed as ID-to-taxonomy assignment tables.
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

    for row in iter_taxonomy_assignment_rows(assignments_path):
        yield {
            'qid': row.get('id'),
            'tax': row.get('taxonomy'),
            'confidence': row.get('confidence'),
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
    run_id_set = {str(iid) for iid in run_ids}
    class_by_id = {}
    for row in iter_classification_rows(class_out):
        qid = str(row['qid'])
        mapped = resolve_short_id(qid, orig_to_short, db)
        if not mapped and qid in run_id_set:
            mapped = qid
        if mapped:
            class_by_id[mapped] = row

    novelty_by_id = {}
    try:
        with open(novelty_metrics_tsv) as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            for row in reader:
                qid = row.get('ID') or row.get('id')
                if not qid:
                    continue

                def val(*names, default='NA'):
                    for name in names:
                        got = row.get(name)
                        if got not in (None, ''):
                            return got
                    return default

                baseline_source = row.get('BaselineDensitySource') or ''
                baseline_available = baseline_source not in ('', 'none', 'NA')

                def primary_val(baseline_name, fallback_name, default='NA'):
                    if baseline_available:
                        got = row.get(baseline_name)
                        if got not in (None, ''):
                            return got
                    return val(fallback_name, default=default)

                novelty_by_id[qid] = {
                    'nearest_identity': primary_val('BaselineNearestIdentity', 'NearestIdentity'),
                    'nearest_hit': primary_val('BaselineNearestHit', 'NearestHit'),
                    'novel': primary_val('BaselineNovel', 'Novel'),
                    'matches_ge_99': primary_val('BaselineMatchesGE99', 'MatchesGE99'),
                    'matches_ge_97': primary_val('BaselineMatchesGE97', 'MatchesGE97'),
                    'matches_ge_95': primary_val('BaselineMatchesGE95', 'MatchesGE95'),
                    'novelty_score': primary_val('BaselineNoveltyScore', 'NoveltyScore'),
                    'crowding': primary_val('BaselineCrowding', 'Crowding'),
                    'sequencing_priority': primary_val('BaselineSequencingPriority', 'SequencingPriority'),
                    'density_source': primary_val('BaselineDensitySource', 'DensitySource'),
                    'all_known_nearest_identity': val('AllKnownNearestIdentity', 'NearestIdentity'),
                    'all_known_nearest_hit': val('AllKnownNearestHit', 'NearestHit'),
                    'all_known_novel': val('AllKnownNovel', 'Novel'),
                    'all_known_matches_ge_99': val('AllKnownMatchesGE99', 'MatchesGE99'),
                    'all_known_matches_ge_97': val('AllKnownMatchesGE97', 'MatchesGE97'),
                    'all_known_matches_ge_95': val('AllKnownMatchesGE95', 'MatchesGE95'),
                    'all_known_novelty_score': val('AllKnownNoveltyScore', 'NoveltyScore'),
                    'all_known_crowding': val('AllKnownCrowding', 'Crowding'),
                    'all_known_sequencing_priority': val('AllKnownSequencingPriority', 'SequencingPriority'),
                    'all_known_density_source': val('AllKnownDensitySource', 'DensitySource'),
                    'reference_nearest_identity': val('ReferenceNearestIdentity'),
                    'reference_nearest_hit': val('ReferenceNearestHit'),
                    'reference_novel': val('ReferenceNovel'),
                    'reference_matches_ge_99': val('ReferenceMatchesGE99'),
                    'reference_matches_ge_97': val('ReferenceMatchesGE97'),
                    'reference_matches_ge_95': val('ReferenceMatchesGE95'),
                    'reference_novelty_score': val('ReferenceNoveltyScore'),
                    'reference_crowding': val('ReferenceCrowding'),
                    'reference_sequencing_priority': val('ReferenceSequencingPriority'),
                    'reference_density_source': val('ReferenceDensitySource'),
                    'partner_id': val('PartnerID'),
                    'selected_for_genome_sequencing': val('SelectedForGenomeSequencing'),
                    'nearest_selected_genome_hit': val('NearestSelectedGenomeHit'),
                    'nearest_selected_genome_identity': val('NearestSelectedGenomeIdentity'),
                    'selected_genome_matches_ge_99': val('SelectedGenomeMatchesGE99'),
                    'selected_genome_matches_ge_97': val('SelectedGenomeMatchesGE97'),
                    'selected_genome_matches_ge_95': val('SelectedGenomeMatchesGE95'),
                    'clade_already_selected_for_genome_sequencing': val('CladeAlreadySelectedForGenomeSequencing'),
                    'genome_sequencing_adjusted_novelty_score': val('GenomeSequencingAdjustedNoveltyScore'),
                    'genome_sequencing_adjusted_priority': val('GenomeSequencingAdjustedPriority'),
                    'genome_sequencing_metadata_source': val('GenomeSequencingMetadataSource'),
                }
    except FileNotFoundError:
        import logging as _logging
        _logging.getLogger(__name__).warning(
            "[ASSESSMENT] Novelty metrics file not found: %s — novelty columns will be NA",
            novelty_metrics_tsv,
        )

    nearest_meta = {}
    nearest_ids = sorted({
        str(row.get('nearest_hit'))
        for row in novelty_by_id.values()
        if row.get('nearest_hit') not in (None, '', 'NA')
    } | {
        str(row.get('all_known_nearest_hit'))
        for row in novelty_by_id.values()
        if row.get('all_known_nearest_hit') not in (None, '', 'NA')
    } | {
        str(row.get('reference_nearest_hit'))
        for row in novelty_by_id.values()
        if row.get('reference_nearest_hit') not in (None, '', 'NA')
    })
    if nearest_ids:
        try:
            with db.connect() as conn:
                cur = conn.cursor()
                placeholders = ','.join('?' for _ in nearest_ids)
                cur.execute(
                    "SELECT s.id, s.dataset, t.taxonomy "
                    "FROM sequences s LEFT JOIN taxonomy t ON s.id = t.id "
                    f"WHERE s.id IN ({placeholders})",
                    tuple(nearest_ids),
                )
                for hit_id, dataset, taxonomy in cur.fetchall():
                    nearest_meta[str(hit_id)] = {
                        'dataset': dataset or 'NA',
                        'taxonomy': taxonomy or 'NA',
                    }
        except Exception:
            nearest_meta = {}

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
        reference_hit = n.get('reference_nearest_hit', 'NA')
        reference_hit_taxonomy = nearest_meta.get(str(reference_hit), {}).get('taxonomy', 'NA')
        if reference_hit not in (None, '', 'NA', 'None') and _best is not None and str(reference_hit) == str(_best):
            reference_hit_taxonomy = _tax if _tax is not None else 'NA'

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
            'nearest_hit_dataset': nearest_meta.get(str(n.get('nearest_hit')), {}).get('dataset', 'NA'),
            'nearest_hit_taxonomy': nearest_meta.get(str(n.get('nearest_hit')), {}).get('taxonomy', 'NA'),
            'nearest_identity': n.get('nearest_identity', 'NA'),
            'matches_ge_99': n.get('matches_ge_99', 'NA'),
            'matches_ge_97': n.get('matches_ge_97', 'NA'),
            'matches_ge_95': n.get('matches_ge_95', 'NA'),
            'novelty_score': n.get('novelty_score', 'NA'),
            'crowding': n.get('crowding', 'NA'),
            'sequencing_priority': n.get('sequencing_priority', 'NA'),
            'density_source': n.get('density_source', 'NA'),
            'all_known_nearest_hit': n.get('all_known_nearest_hit', 'NA'),
            'all_known_nearest_hit_dataset': nearest_meta.get(str(n.get('all_known_nearest_hit')), {}).get('dataset', 'NA'),
            'all_known_nearest_hit_taxonomy': nearest_meta.get(str(n.get('all_known_nearest_hit')), {}).get('taxonomy', 'NA'),
            'all_known_nearest_identity': n.get('all_known_nearest_identity', 'NA'),
            'all_known_matches_ge_99': n.get('all_known_matches_ge_99', 'NA'),
            'all_known_matches_ge_97': n.get('all_known_matches_ge_97', 'NA'),
            'all_known_matches_ge_95': n.get('all_known_matches_ge_95', 'NA'),
            'all_known_novelty_score': n.get('all_known_novelty_score', 'NA'),
            'all_known_crowding': n.get('all_known_crowding', 'NA'),
            'all_known_sequencing_priority': n.get('all_known_sequencing_priority', 'NA'),
            'all_known_density_source': n.get('all_known_density_source', 'NA'),
            'reference_nearest_hit': reference_hit,
            'reference_nearest_hit_taxonomy': reference_hit_taxonomy,
            'reference_nearest_identity': n.get('reference_nearest_identity', 'NA'),
            'reference_matches_ge_99': n.get('reference_matches_ge_99', 'NA'),
            'reference_matches_ge_97': n.get('reference_matches_ge_97', 'NA'),
            'reference_matches_ge_95': n.get('reference_matches_ge_95', 'NA'),
            'reference_novelty_score': n.get('reference_novelty_score', 'NA'),
            'reference_crowding': n.get('reference_crowding', 'NA'),
            'reference_sequencing_priority': n.get('reference_sequencing_priority', 'NA'),
            'reference_density_source': n.get('reference_density_source', 'NA'),
            'partner_id': n.get('partner_id', 'NA'),
            'selected_for_genome_sequencing': n.get('selected_for_genome_sequencing', 'NA'),
            'nearest_selected_genome_hit': n.get('nearest_selected_genome_hit', 'NA'),
            'nearest_selected_genome_identity': n.get('nearest_selected_genome_identity', 'NA'),
            'selected_genome_matches_ge_99': n.get('selected_genome_matches_ge_99', 'NA'),
            'selected_genome_matches_ge_97': n.get('selected_genome_matches_ge_97', 'NA'),
            'selected_genome_matches_ge_95': n.get('selected_genome_matches_ge_95', 'NA'),
            'clade_already_selected_for_genome_sequencing': n.get('clade_already_selected_for_genome_sequencing', 'NA'),
            'genome_sequencing_adjusted_novelty_score': n.get('genome_sequencing_adjusted_novelty_score', 'NA'),
            'genome_sequencing_adjusted_priority': n.get('genome_sequencing_adjusted_priority', 'NA'),
            'genome_sequencing_metadata_source': n.get('genome_sequencing_metadata_source', 'NA'),
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

    NearestHit / NearestHitDataset / NearestHitTaxonomy / NearestIdentity / MatchesGE* / NoveltyScore / Crowding / SequencingPriority
        — from the baseline/cultured sequences already stored in the DB
          when a baseline pool is configured, otherwise from the all-known pool.
          NearestHit is the closest sequence YOU have previously submitted.
          These columns drive the novelty assessment.

    AllKnownNearestHit / AllKnownNearestIdentity / AllKnownNoveltyScore / ...
        — same novelty metrics against all non-current datasets stored in
          the DB, so cultured-baseline novelty can be compared to wider
          project-level novelty.

    ReferenceNearestHit / ReferenceNearestIdentity / ReferenceNoveltyScore / ...
        — same nearest-hit and density metrics against the selected external
          reference FASTA, usually the main GTDB reference supplied with --ref.
          This separates "closest cultured/project isolate" from "closest
          taxonomic reference sequence".

    PartnerID / SelectedForGenomeSequencing / NearestSelectedGenomeHit / ...
        — rolling partner metadata and selected-for-WGS neighbourhood metrics.
          If a nearby sequence at >=97% identity has already been selected for
          WGS, the clade is marked as already represented.

    Taxonomy_<DB> / ClassificationHit_<DB> / Identity_<DB> / Confidence_<DB>  (repeated per alt-ref DB)
        — same information from each additional reference database supplied via
          --alt-ref / --alt-taxa.  The <DB> suffix is the database name derived
          from the filename or --alt-ref-name.
    """
    import re as _re
    p = Path(path)

    # Discover alt-db column groups from the first row that contains them
    alt_ref_safes: List[str] = []
    include_mwl = False
    for row in rows:
        alt_ref_safes = sorted(
            k[len('alt_tax_'):]
            for k in row
            if k.startswith('alt_tax_')
        )
        include_mwl = any(k.startswith('mwl_') or k == 'evaluation_score' for k in row)
        break

    with open(p, 'w') as fh:
        # ── Fixed columns ────────────────────────────────────────────────────
        mwl_header = ''
        if include_mwl:
            mwl_header = (
                '\tEvaluationScore'          # optional MWL-aware score; InvestigationScore blended with MWLScore
                '\tMWLMatch'
                '\tMWLID'
                '\tMWLMatchedRank'
                '\tMWLMatchedTaxon'
                '\tMWLTaxonomicScore'
                '\tMWLIdentity'
                '\tMWLScore'
                '\tMWLRole'
            )
        base_header = (
            'ID'
            '\tTaxonomy'                    # primary ref DB  ─┐ classification
            '\tClassificationHit'           # ref accession    │ (reference DB)
            '\tClassificationIdentity'      # vsearch %id      │
            '\tClassificationConfidence'    # score 0–1       ─┘
            '\tNearestHit'                  # preload seq ─┐ novelty
            '\tNearestHitDataset'           # dataset label │
            '\tNearestHitTaxonomy'          # taxonomy      │
            '\tNearestIdentity'             # %id           │ (preload DB)
            '\tMatchesGE99'                 # density        │
            '\tMatchesGE97'                 #               │
            '\tMatchesGE95'                 #              ─┘
            '\tNoveltyScore'
            '\tCrowding'
            '\tSequencingPriority'
            '\tDensitySource'
            '\tAllKnownNearestHit'
            '\tAllKnownNearestHitDataset'
            '\tAllKnownNearestHitTaxonomy'
            '\tAllKnownNearestIdentity'
            '\tAllKnownMatchesGE99'
            '\tAllKnownMatchesGE97'
            '\tAllKnownMatchesGE95'
            '\tAllKnownNoveltyScore'
            '\tAllKnownCrowding'
            '\tAllKnownSequencingPriority'
            '\tAllKnownDensitySource'
            '\tReferenceNearestHit'
            '\tReferenceNearestHitTaxonomy'
            '\tReferenceNearestIdentity'
            '\tReferenceMatchesGE99'
            '\tReferenceMatchesGE97'
            '\tReferenceMatchesGE95'
            '\tReferenceNoveltyScore'
            '\tReferenceCrowding'
            '\tReferenceSequencingPriority'
            '\tReferenceDensitySource'
            '\tPartnerID'
            '\tSelectedForGenomeSequencing'
            '\tNearestSelectedGenomeHit'
            '\tNearestSelectedGenomeIdentity'
            '\tSelectedGenomeMatchesGE99'
            '\tSelectedGenomeMatchesGE97'
            '\tSelectedGenomeMatchesGE95'
            '\tCladeAlreadySelectedForGenomeSequencing'
            '\tGenomeSequencingAdjustedNoveltyScore'
            '\tGenomeSequencingAdjustedPriority'
            '\tGenomeSequencingMetadataSource'
            '\tPhyloIsolation'              # normalised leaf branch length (0–1); higher = more isolated in tree
            '\tInvestigationScore'          # composite score (0–100) combining novelty + phylo isolation + density + taxonomy conflict
            + mwl_header +
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
                f"\t{row.get('nearest_hit_dataset', 'NA')}"
                f"\t{row.get('nearest_hit_taxonomy', 'NA')}"
                f"\t{row['nearest_identity']}"
                f"\t{row['matches_ge_99']}"
                f"\t{row['matches_ge_97']}"
                f"\t{row['matches_ge_95']}"
                f"\t{row['novelty_score']}"
                f"\t{row['crowding']}"
                f"\t{row['sequencing_priority']}"
                f"\t{row.get('density_source', 'NA')}"
                f"\t{row.get('all_known_nearest_hit', 'NA')}"
                f"\t{row.get('all_known_nearest_hit_dataset', 'NA')}"
                f"\t{row.get('all_known_nearest_hit_taxonomy', 'NA')}"
                f"\t{row.get('all_known_nearest_identity', 'NA')}"
                f"\t{row.get('all_known_matches_ge_99', 'NA')}"
                f"\t{row.get('all_known_matches_ge_97', 'NA')}"
                f"\t{row.get('all_known_matches_ge_95', 'NA')}"
                f"\t{row.get('all_known_novelty_score', 'NA')}"
                f"\t{row.get('all_known_crowding', 'NA')}"
                f"\t{row.get('all_known_sequencing_priority', 'NA')}"
                f"\t{row.get('all_known_density_source', 'NA')}"
                f"\t{row.get('reference_nearest_hit', 'NA')}"
                f"\t{row.get('reference_nearest_hit_taxonomy', 'NA')}"
                f"\t{row.get('reference_nearest_identity', 'NA')}"
                f"\t{row.get('reference_matches_ge_99', 'NA')}"
                f"\t{row.get('reference_matches_ge_97', 'NA')}"
                f"\t{row.get('reference_matches_ge_95', 'NA')}"
                f"\t{row.get('reference_novelty_score', 'NA')}"
                f"\t{row.get('reference_crowding', 'NA')}"
                f"\t{row.get('reference_sequencing_priority', 'NA')}"
                f"\t{row.get('reference_density_source', 'NA')}"
                f"\t{row.get('partner_id', 'NA')}"
                f"\t{row.get('selected_for_genome_sequencing', 'NA')}"
                f"\t{row.get('nearest_selected_genome_hit', 'NA')}"
                f"\t{row.get('nearest_selected_genome_identity', 'NA')}"
                f"\t{row.get('selected_genome_matches_ge_99', 'NA')}"
                f"\t{row.get('selected_genome_matches_ge_97', 'NA')}"
                f"\t{row.get('selected_genome_matches_ge_95', 'NA')}"
                f"\t{row.get('clade_already_selected_for_genome_sequencing', 'NA')}"
                f"\t{row.get('genome_sequencing_adjusted_novelty_score', 'NA')}"
                f"\t{row.get('genome_sequencing_adjusted_priority', 'NA')}"
                f"\t{row.get('genome_sequencing_metadata_source', 'NA')}"
                f"\t{row.get('phylo_isolation', 'NA')}"
                f"\t{row.get('investigation_score', 'NA')}"
            )
            if include_mwl:
                base += (
                    f"\t{row.get('evaluation_score', 'NA')}"
                    f"\t{row.get('mwl_match', 'NA')}"
                    f"\t{row.get('mwl_id', 'NA')}"
                    f"\t{row.get('mwl_matched_rank', 'NA')}"
                    f"\t{row.get('mwl_matched_taxon', 'NA')}"
                    f"\t{row.get('mwl_taxonomic_score', 'NA')}"
                    f"\t{row.get('mwl_identity', 'NA')}"
                    f"\t{row.get('mwl_score', 'NA')}"
                    f"\t{row.get('mwl_role', 'NA')}"
                )
            base += (
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


def _as_float(value, default=None):
    try:
        if value in (None, '', 'NA', 'None'):
            return default
        return float(value)
    except Exception:
        return default


def _board_recommendation(row: dict) -> tuple[str, str]:
    """Return a concise WGS-selection recommendation and short rationale."""
    flags = []
    placement_flags = str(row.get('placement_flags', '') or '')
    adjusted_priority = row.get('genome_sequencing_adjusted_priority')
    raw_priority = row.get('sequencing_priority', 'NA')
    priority = adjusted_priority if adjusted_priority not in (None, '', 'NA') else raw_priority
    selected = str(row.get('selected_for_genome_sequencing', 'NA')).lower() == 'true'
    clade_selected = str(row.get('clade_already_selected_for_genome_sequencing', 'NA')).lower() == 'true'
    mwl_match = str(row.get('mwl_match', 'No') or 'No')
    baseline_identity = _as_float(row.get('nearest_identity'))
    reference_identity = _as_float(row.get('reference_nearest_identity'))
    classification_identity = _as_float(row.get('classification_identity'))

    if selected:
        recommendation = 'Already selected'
        flags.append('current isolate marked selected')
    elif clade_selected:
        recommendation = 'Deprioritise - clade represented'
        flags.append('nearby selected genome >=97% 16S identity')
    elif placement_flags:
        recommendation = 'Review before selection'
        flags.append(f'placement flags: {placement_flags}')
    elif priority == 'HIGH':
        recommendation = 'Strong candidate'
    elif priority == 'MEDIUM':
        recommendation = 'Secondary candidate'
    elif mwl_match == 'Yes':
        recommendation = 'Review MWL candidate'
    else:
        recommendation = 'Lower priority'

    if mwl_match == 'Yes':
        matched = row.get('mwl_matched_taxon', 'NA')
        rank = row.get('mwl_matched_rank', 'NA')
        flags.append(f'MWL match {rank}:{matched}')
    if baseline_identity is not None:
        flags.append(f'baseline nearest {baseline_identity:.2f}%')
    if reference_identity is not None:
        flags.append(f'reference nearest {reference_identity:.2f}%')
    if classification_identity is not None and classification_identity < 95.0:
        flags.append(f'low taxonomy identity {classification_identity:.2f}%')

    return recommendation, '; '.join(flags) if flags else 'No major flags'


def write_selection_summary_tsv(path: str | Path, rows):
    """Write a concise SAB-facing WGS selection summary.

    This is intentionally much smaller than sequence_assessment.tsv. It keeps
    the evidence needed for a selection discussion while leaving the full audit
    trail in the main assessment table.
    """
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    headers = [
        'SequenceID',
        'PartnerID',
        'Recommendation',
        'AdjustedPriority',
        'AdjustedNoveltyScore',
        'RawNoveltyPriority',
        'BaselineNearestIdentity',
        'BaselineNearestHit',
        'BaselineNearestHitTaxonomy',
        'AllKnownNearestIdentity',
        'ReferenceNearestIdentity',
        'ReferenceNearestHit',
        'Taxonomy',
        'ClassificationIdentity',
        'MWLMatch',
        'MWLMatchedRank',
        'MWLMatchedTaxon',
        'MWLScore',
        'SelectedForGenomeSequencing',
        'CladeAlreadySelectedForGenomeSequencing',
        'NearestSelectedGenomeHit',
        'NearestSelectedGenomeIdentity',
        'InTree',
        'ClusterRepresentative',
        'ClusterSize',
        'KeyRationale',
    ]
    with open(p, 'w') as fh:
        fh.write('\t'.join(headers) + '\n')
        for row in rows:
            recommendation, rationale = _board_recommendation(row)
            adjusted_priority = row.get('genome_sequencing_adjusted_priority')
            if adjusted_priority in (None, '', 'NA'):
                adjusted_priority = row.get('sequencing_priority', 'NA')
            adjusted_score = row.get('genome_sequencing_adjusted_novelty_score')
            if adjusted_score in (None, '', 'NA'):
                adjusted_score = row.get('novelty_score', 'NA')
            values = [
                row.get('id', 'NA'),
                row.get('partner_id', 'NA'),
                recommendation,
                adjusted_priority,
                adjusted_score,
                row.get('sequencing_priority', 'NA'),
                row.get('nearest_identity', 'NA'),
                row.get('nearest_hit', 'NA'),
                row.get('nearest_hit_taxonomy', 'NA'),
                row.get('all_known_nearest_identity', 'NA'),
                row.get('reference_nearest_identity', 'NA'),
                row.get('reference_nearest_hit', 'NA'),
                row.get('taxonomy', 'NA'),
                row.get('classification_identity', 'NA'),
                row.get('mwl_match', 'NA'),
                row.get('mwl_matched_rank', 'NA'),
                row.get('mwl_matched_taxon', 'NA'),
                row.get('mwl_score', 'NA'),
                row.get('selected_for_genome_sequencing', 'NA'),
                row.get('clade_already_selected_for_genome_sequencing', 'NA'),
                row.get('nearest_selected_genome_hit', 'NA'),
                row.get('nearest_selected_genome_identity', 'NA'),
                row.get('in_tree', 'Unknown'),
                row.get('cluster_representative', 'N/A'),
                row.get('cluster_size', '1'),
                rationale,
            ]
            safe = [str(v if v is not None else 'NA').replace('\t', ' ').replace('\n', ' ') for v in values]
            fh.write('\t'.join(safe) + '\n')
    return str(p)


def write_baseline_hits_tsv(path: str | Path, rows):
    """Write a concise nearest-baseline hit report for evaluated sequences."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, 'w') as fh:
        fh.write(
            'ID\tNearestBaselineHit\tBaselineDataset\tNearestIdentity\t'
            'NearestBaselineTaxonomy\tNoveltyScore\tCrowding\tSequencingPriority\n'
        )
        for row in rows:
            hit = row.get('nearest_hit', 'NA')
            if hit in (None, '', 'NA'):
                continue
            source = row.get('density_source', 'NA')
            if source not in (None, '', 'NA', 'none') and not str(source).startswith('baseline'):
                continue
            fh.write(
                f"{row.get('id', 'NA')}\t"
                f"{hit}\t"
                f"{row.get('nearest_hit_dataset', 'NA')}\t"
                f"{row.get('nearest_identity', 'NA')}\t"
                f"{row.get('nearest_hit_taxonomy', 'NA')}\t"
                f"{row.get('novelty_score', 'NA')}\t"
                f"{row.get('crowding', 'NA')}\t"
                f"{row.get('sequencing_priority', 'NA')}\n"
            )
    return str(p)
