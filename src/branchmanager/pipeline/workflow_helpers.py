from __future__ import annotations

import gzip
import os
import csv
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

from branchmanager.taxonomy_io import iter_taxonomy_assignment_rows
from branchmanager.taxonomy import (
    normalise_taxon_name,
    parse_reference_header_taxonomy,
    parse_taxon_string,
    taxonomy_matches_kingdom,
)
from branchmanager.utils.fasta import read_fasta


ClassificationRow = Dict[str, object]


def _optional_numeric(parts: Sequence[str], index: int, *, integer: bool = False):
    if index >= len(parts) or parts[index] in ('', 'NA', 'None'):
        return None
    try:
        value = float(parts[index])
        return int(value) if integer else value
    except (TypeError, ValueError):
        return None


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
                'query_coverage': _optional_numeric(parts, 5),
                'target_coverage': _optional_numeric(parts, 6),
                'alignment_length': _optional_numeric(parts, 7, integer=True),
                'query_length': _optional_numeric(parts, 8, integer=True),
                'target_length': _optional_numeric(parts, 9, integer=True),
                'mismatches': _optional_numeric(parts, 10, integer=True),
                'gaps': _optional_numeric(parts, 11, integer=True),
            }


def iter_assignment_rows(assignments_path: str, source_fasta_path: Optional[str] = None) -> Iterator[dict[str, object]]:
    """Yield normalised taxonomy-assignment rows from table files or FASTA headers.

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
            cur.execute('DELETE FROM colours WHERE id = ?', (iid,))
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
        query_coverage = row.get('query_coverage')
        tax = row.get('tax')
        nov = novelty_by_qid.get(qid, {})
        nearest_identity = nov.get('nearest_identity')

        if not best or best == 'NA' or identity in (None, 0, 0.0):
            flags.append('NO_REFERENCE_HIT')
        if isinstance(identity, (int, float)) and float(identity) < 95.0:
            flags.append('LOW_CLASSIFICATION_IDENTITY')
        if isinstance(confidence, (int, float)) and float(confidence) < 0.8:
            flags.append('LOW_CONFIDENCE')
        if isinstance(query_coverage, (int, float)) and float(query_coverage) < 90.0:
            flags.append('LOW_QUERY_COVERAGE')
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


def _taxonomy_agreement(query_taxonomy, comparison_taxonomy) -> tuple[str, str]:
    """Return deepest shared rank and first explicit conflict for two lineages."""
    ranks = [('d', 'domain'), ('p', 'phylum'), ('c', 'class'), ('o', 'order'),
             ('f', 'family'), ('g', 'genus'), ('s', 'species')]
    query = parse_taxon_string(str(query_taxonomy or ''))
    comparison = parse_taxon_string(str(comparison_taxonomy or ''))
    if not query or not comparison:
        return 'NA', 'NA'
    deepest = 'none'
    for key, name in ranks:
        left = normalise_taxon_name(query.get(key, ''))
        right = normalise_taxon_name(comparison.get(key, ''))
        if not left or not right:
            continue
        if left != right:
            return deepest, f'{name}:{query.get(key)}!={comparison.get(key)}'
        deepest = name
    return deepest, 'none'


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
                            or 'duplicate' if omitted from the final tree as an exact duplicate,
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

    # Reload stored evidence when a submitted ID was already present in the
    # project database or was not represented in the current classifier output.
    if run_id_set:
        try:
            placeholders = ','.join('?' for _ in run_id_set)
            with db.connect() as conn:
                taxonomy_rows = conn.execute(
                    f"SELECT id, taxonomy, confidence FROM taxonomy "
                    f"WHERE id IN ({placeholders}) ORDER BY rowid",
                    tuple(sorted(run_id_set)),
                ).fetchall()
                distance_rows = conn.execute(
                    f"SELECT id, nearest, identity FROM distances "
                    f"WHERE id IN ({placeholders}) ORDER BY rowid",
                    tuple(sorted(run_id_set)),
                ).fetchall()
            db_classification = {}
            for sid, taxonomy, confidence in taxonomy_rows:
                entry = db_classification.setdefault(str(sid), {'qid': str(sid)})
                if taxonomy not in (None, ''):
                    entry['tax'] = taxonomy
                if confidence is not None:
                    entry['confidence'] = confidence
            for sid, nearest, identity in distance_rows:
                entry = db_classification.setdefault(str(sid), {'qid': str(sid)})
                if nearest not in (None, ''):
                    entry['best'] = nearest
                if identity is not None:
                    entry['identity'] = identity
            for sid, entry in db_classification.items():
                class_by_id.setdefault(sid, entry)
            evidence_rows = db.get_classification_evidence(run_id_set)
            evidence_by_id = {}
            for (sid, _ref_db), evidence in evidence_rows.items():
                evidence_by_id.setdefault(sid, evidence)
            for sid, evidence in evidence_by_id.items():
                class_by_id.setdefault(sid, {'qid': sid}).update({
                    key: evidence.get(key) for key in (
                        'query_coverage', 'target_coverage', 'alignment_length', 'query_length',
                        'target_length', 'mismatches', 'gaps',
                    )
                })
        except Exception:
            pass

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

                def baseline_val(name, default='NA'):
                    if baseline_available:
                        got = row.get(name)
                        if got not in (None, ''):
                            return got
                    return default

                novelty_by_id[qid] = {
                    'nearest_identity': baseline_val('BaselineNearestIdentity'),
                    'nearest_hit': baseline_val('BaselineNearestHit'),
                    'nearest_query_coverage': baseline_val('BaselineNearestQueryCoverage'),
                    'nearest_alignment_length': baseline_val('BaselineNearestAlignmentLength'),
                    'novel': baseline_val('BaselineNovel'),
                    'matches_ge_99': baseline_val('BaselineMatchesGE99'),
                    'matches_ge_97': baseline_val('BaselineMatchesGE97'),
                    'matches_ge_95': baseline_val('BaselineMatchesGE95'),
                    'novelty_score': baseline_val('BaselineNoveltyScore'),
                    'crowding': baseline_val('BaselineCrowding'),
                    'sequencing_priority': baseline_val('BaselineSequencingPriority'),
                    'density_source': baseline_val('BaselineDensitySource'),
                    'project_nearest_identity': val('ProjectNearestIdentity'),
                    'project_nearest_hit': val('ProjectNearestHit'),
                    'project_nearest_query_coverage': val('ProjectNearestQueryCoverage'),
                    'project_nearest_alignment_length': val('ProjectNearestAlignmentLength'),
                    'project_novel': val('ProjectNovel'),
                    'project_matches_ge_99': val('ProjectMatchesGE99'),
                    'project_matches_ge_97': val('ProjectMatchesGE97'),
                    'project_matches_ge_95': val('ProjectMatchesGE95'),
                    'project_novelty_score': val('ProjectNoveltyScore'),
                    'project_crowding': val('ProjectCrowding'),
                    'project_sequencing_priority': val('ProjectSequencingPriority'),
                    'project_density_source': val('ProjectDensitySource'),
                    'reference_nearest_identity': val('ReferenceNearestIdentity'),
                    'reference_nearest_hit': val('ReferenceNearestHit'),
                    'reference_nearest_query_coverage': val('ReferenceNearestQueryCoverage'),
                    'reference_nearest_alignment_length': val('ReferenceNearestAlignmentLength'),
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
                    'already_sequenced': val('GenomeAlreadySequenced'),
                    'nearest_genome_hit': val('NearestGenomeHit'),
                    'nearest_genome_identity': val('NearestGenomeIdentity'),
                    'genome_collection_matches_ge_99': val('GenomeCollectionMatchesGE99'),
                    'genome_collection_matches_ge_97': val('GenomeCollectionMatchesGE97'),
                    'genome_collection_matches_ge_95': val('GenomeCollectionMatchesGE95'),
                    'related_genome_clade_ge_97': val('RelatedGenomeCladeGE97'),
                    'genome_committed_count_same_species': val(
                        'CommittedGenomeCountSameAssessmentSpecies',
                        default=val('AvailableGenomeCountSameAssessmentSpecies',
                                    default=val('GenomeCommittedCountSameAssessmentSpecies')),
                    ),
                    'genome_available_count_same_species': val(
                        'BaselineGenomeCountSameAssessmentSpecies',
                        default=val('GenomeAvailableCountSameAssessmentSpecies'),
                    ),
                    'genome_selected_count_same_species': val(
                        'SequencedPartnerGenomeCountSameAssessmentSpecies',
                        default=val('GenomeSelectedCountSameAssessmentSpecies'),
                    ),
                    'genome_pending_count_same_species': val(
                        'SelectedPendingGenomeCountSameAssessmentSpecies',
                    ),
                    'pangenome_target': val('PangenomeTarget', default='9'),
                    'pangenome_gap': val('PangenomeGap'),
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
        str(row.get('project_nearest_hit'))
        for row in novelty_by_id.values()
        if row.get('project_nearest_hit') not in (None, '', 'NA')
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

        # Some tree builders can omit an exact duplicate even though every
        # isolate remains in the assessment and project database.
        if in_tree == 'No' and cluster_rep == 'self':
            cluster_rep = 'duplicate'

        baseline_hit_taxonomy = nearest_meta.get(str(n.get('nearest_hit')), {}).get('taxonomy', 'NA')
        baseline_agreement_rank, baseline_taxonomy_conflict = _taxonomy_agreement(
            _tax,
            baseline_hit_taxonomy,
        )

        assessment_rows.append({
            'id': iid,
            # ClassificationHit = best vsearch hit in the REFERENCE DB (GTDB/SILVA) → drives Taxonomy
            'taxonomy': _tax if _tax is not None else 'NA',
            'classification_hit': _best if _best is not None else 'NA',
            'classification_identity': _ident if _ident is not None else 'NA',
            'classification_confidence': _conf if _conf is not None else 'NA',
            'classification_query_coverage': c.get('query_coverage', 'NA'),
            'classification_target_coverage': c.get('target_coverage', 'NA'),
            'classification_alignment_length': c.get('alignment_length', 'NA'),
            'classification_query_length': c.get('query_length', 'NA'),
            'classification_target_length': c.get('target_length', 'NA'),
            'classification_mismatches': c.get('mismatches', 'NA'),
            'classification_gaps': c.get('gaps', 'NA'),
            # NearestHit = best vsearch hit among registered baseline sequences already in the DB → drives NoveltyScore
            'nearest_hit': n.get('nearest_hit', 'NA'),
            'nearest_hit_dataset': nearest_meta.get(str(n.get('nearest_hit')), {}).get('dataset', 'NA'),
            'nearest_hit_taxonomy': baseline_hit_taxonomy,
            'baseline_taxonomy_agreement_rank': baseline_agreement_rank,
            'baseline_taxonomy_conflict': baseline_taxonomy_conflict,
            'nearest_identity': n.get('nearest_identity', 'NA'),
            'nearest_query_coverage': n.get('nearest_query_coverage', 'NA'),
            'nearest_alignment_length': n.get('nearest_alignment_length', 'NA'),
            'novel': n.get('novel', 'NA'),
            'matches_ge_99': n.get('matches_ge_99', 'NA'),
            'matches_ge_97': n.get('matches_ge_97', 'NA'),
            'matches_ge_95': n.get('matches_ge_95', 'NA'),
            'novelty_score': n.get('novelty_score', 'NA'),
            'crowding': n.get('crowding', 'NA'),
            'sequencing_priority': n.get('sequencing_priority', 'NA'),
            'density_source': n.get('density_source', 'NA'),
            'project_nearest_hit': n.get('project_nearest_hit', 'NA'),
            'project_nearest_hit_dataset': nearest_meta.get(str(n.get('project_nearest_hit')), {}).get('dataset', 'NA'),
            'project_nearest_hit_taxonomy': nearest_meta.get(str(n.get('project_nearest_hit')), {}).get('taxonomy', 'NA'),
            'project_nearest_identity': n.get('project_nearest_identity', 'NA'),
            'project_nearest_query_coverage': n.get('project_nearest_query_coverage', 'NA'),
            'project_nearest_alignment_length': n.get('project_nearest_alignment_length', 'NA'),
            'project_novel': n.get('project_novel', 'NA'),
            'project_matches_ge_99': n.get('project_matches_ge_99', 'NA'),
            'project_matches_ge_97': n.get('project_matches_ge_97', 'NA'),
            'project_matches_ge_95': n.get('project_matches_ge_95', 'NA'),
            'project_novelty_score': n.get('project_novelty_score', 'NA'),
            'project_crowding': n.get('project_crowding', 'NA'),
            'project_sequencing_priority': n.get('project_sequencing_priority', 'NA'),
            'project_density_source': n.get('project_density_source', 'NA'),
            'reference_nearest_hit': reference_hit,
            'reference_nearest_hit_taxonomy': reference_hit_taxonomy,
            'reference_nearest_identity': n.get('reference_nearest_identity', 'NA'),
            'reference_nearest_query_coverage': n.get('reference_nearest_query_coverage', 'NA'),
            'reference_nearest_alignment_length': n.get('reference_nearest_alignment_length', 'NA'),
            'reference_novel': n.get('reference_novel', 'NA'),
            'reference_matches_ge_99': n.get('reference_matches_ge_99', 'NA'),
            'reference_matches_ge_97': n.get('reference_matches_ge_97', 'NA'),
            'reference_matches_ge_95': n.get('reference_matches_ge_95', 'NA'),
            'reference_novelty_score': n.get('reference_novelty_score', 'NA'),
            'reference_crowding': n.get('reference_crowding', 'NA'),
            'reference_sequencing_priority': n.get('reference_sequencing_priority', 'NA'),
            'reference_density_source': n.get('reference_density_source', 'NA'),
            'partner_id': n.get('partner_id', 'NA'),
            'selected_for_genome_sequencing': n.get('selected_for_genome_sequencing', 'NA'),
            'already_sequenced': n.get('already_sequenced', 'NA'),
            'nearest_genome_hit': n.get('nearest_genome_hit', 'NA'),
            'nearest_genome_identity': n.get('nearest_genome_identity', 'NA'),
            'genome_collection_matches_ge_99': n.get('genome_collection_matches_ge_99', 'NA'),
            'genome_collection_matches_ge_97': n.get('genome_collection_matches_ge_97', 'NA'),
            'genome_collection_matches_ge_95': n.get('genome_collection_matches_ge_95', 'NA'),
            'related_genome_clade_ge_97': n.get('related_genome_clade_ge_97', 'NA'),
            'genome_committed_count_same_species': n.get('genome_committed_count_same_species', 'NA'),
            'genome_available_count_same_species': n.get('genome_available_count_same_species', 'NA'),
            'genome_selected_count_same_species': n.get('genome_selected_count_same_species', 'NA'),
            'genome_pending_count_same_species': n.get('genome_pending_count_same_species', 'NA'),
            'pangenome_target': n.get('pangenome_target', '9'),
            'pangenome_gap': n.get('pangenome_gap', 'NA'),
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


def write_sequence_assessment_tsv(path: str | Path, rows, assessment_db_name: str = 'GTDB'):
    """Write the complete, component-based per-isolate assessment audit."""
    p = Path(path)
    rows = list(rows)
    label = 'GTDB' if 'gtdb' in str(assessment_db_name).lower() else ''.join(
        char if char.isalnum() else '_' for char in str(assessment_db_name or 'Assessment')
    ).strip('_') or 'Assessment'
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
    fields = [
        ('ID', 'id'), ('PartnerID', 'partner_id'),
        (f'{label}Taxonomy', 'taxonomy'), (f'{label}ClassificationHit', 'classification_hit'),
        (f'{label}ClassificationIdentity', 'classification_identity'),
        (f'{label}ClassificationConfidence', 'classification_confidence'),
        (f'{label}QueryCoverage', 'classification_query_coverage'),
        (f'{label}TargetCoverage', 'classification_target_coverage'),
        (f'{label}AlignmentLength', 'classification_alignment_length'),
        (f'{label}QueryLength', 'classification_query_length'),
        (f'{label}TargetLength', 'classification_target_length'),
        (f'{label}Mismatches', 'classification_mismatches'),
        (f'{label}Gaps', 'classification_gaps'),
        ('BaselineNearestHit', 'nearest_hit'), ('BaselineNearestHitDataset', 'nearest_hit_dataset'),
        ('BaselineNearestHitTaxonomy', 'nearest_hit_taxonomy'),
        ('BaselineTaxonomyAgreementRank', 'baseline_taxonomy_agreement_rank'),
        ('BaselineTaxonomyConflict', 'baseline_taxonomy_conflict'),
        ('BaselineNearestIdentity', 'nearest_identity'), ('BaselineNovel', 'novel'),
        ('BaselineNearestQueryCoverage', 'nearest_query_coverage'),
        ('BaselineNearestAlignmentLength', 'nearest_alignment_length'),
        ('BaselineMatchesGE99', 'matches_ge_99'), ('BaselineMatchesGE97', 'matches_ge_97'),
        ('BaselineMatchesGE95', 'matches_ge_95'), ('BaselineNoveltyScore', 'novelty_score'),
        ('BaselineCrowding', 'crowding'), ('BaselinePriority', 'sequencing_priority'),
        ('BaselineSource', 'density_source'),
        ('ProjectNearestHit', 'project_nearest_hit'),
        ('ProjectNearestHitDataset', 'project_nearest_hit_dataset'),
        ('ProjectNearestHitTaxonomy', 'project_nearest_hit_taxonomy'),
        ('ProjectNearestIdentity', 'project_nearest_identity'), ('ProjectNovel', 'project_novel'),
        ('ProjectNearestQueryCoverage', 'project_nearest_query_coverage'),
        ('ProjectNearestAlignmentLength', 'project_nearest_alignment_length'),
        ('ProjectMatchesGE99', 'project_matches_ge_99'), ('ProjectMatchesGE97', 'project_matches_ge_97'),
        ('ProjectMatchesGE95', 'project_matches_ge_95'), ('ProjectNoveltyScore', 'project_novelty_score'),
        ('ProjectCrowding', 'project_crowding'), ('ProjectPriority', 'project_sequencing_priority'),
        ('ProjectSource', 'project_density_source'),
        (f'{label}ReferenceNearestHit', 'reference_nearest_hit'),
        (f'{label}ReferenceNearestHitTaxonomy', 'reference_nearest_hit_taxonomy'),
        (f'{label}ReferenceNearestIdentity', 'reference_nearest_identity'),
        (f'{label}ReferenceNearestQueryCoverage', 'reference_nearest_query_coverage'),
        (f'{label}ReferenceNearestAlignmentLength', 'reference_nearest_alignment_length'),
        (f'{label}ReferenceNovel', 'reference_novel'),
        (f'{label}ReferenceMatchesGE99', 'reference_matches_ge_99'),
        (f'{label}ReferenceMatchesGE97', 'reference_matches_ge_97'),
        (f'{label}ReferenceMatchesGE95', 'reference_matches_ge_95'),
        (f'{label}ReferenceNoveltyScore', 'reference_novelty_score'),
        (f'{label}ReferenceCrowding', 'reference_crowding'),
        (f'{label}ReferenceSource', 'reference_density_source'),
        ('SelectedForGenomeSequencing', 'selected_for_genome_sequencing'),
        ('GenomeAlreadySequenced', 'already_sequenced'),
        ('NearestGenomeHit', 'nearest_genome_hit'), ('NearestGenomeIdentity', 'nearest_genome_identity'),
        ('GenomeCollectionMatchesGE99', 'genome_collection_matches_ge_99'),
        ('GenomeCollectionMatchesGE97', 'genome_collection_matches_ge_97'),
        ('GenomeCollectionMatchesGE95', 'genome_collection_matches_ge_95'),
        ('RelatedGenomeCladeGE97', 'related_genome_clade_ge_97'),
        ('BaselineGenomesSameAssessmentSpecies', 'genome_available_count_same_species'),
        ('SequencedPartnerGenomesSameAssessmentSpecies', 'genome_selected_count_same_species'),
        ('SelectedPendingGenomesSameAssessmentSpecies', 'genome_pending_count_same_species'),
        ('CommittedGenomesSameAssessmentSpecies', 'genome_committed_count_same_species'),
        ('PangenomeTarget', 'pangenome_target'), ('PangenomeGap', 'pangenome_gap'),
        ('GenomeCollectionSource', 'genome_sequencing_metadata_source'),
        ('MarkerQCClass', 'marker_qc_class'),
        ('MarkerQCRecommendation', 'marker_qc_recommendation'),
        ('MarkerManualReviewStatus', 'marker_manual_review_status'),
        ('MarkerQCFlag', 'marker_qc_flag'),
        ('MarkerQCReasons', 'marker_qc_reasons'),
        ('MarkerSourceManifest', 'marker_source_manifest'),
        ('ChimeraCall', 'chimera_call'), ('UCHIMEScore', 'chimera_score'),
    ]
    if include_mwl:
        fields.extend([
            ('MWLMatch', 'mwl_match'), ('MWLID', 'mwl_id'),
            ('MWLMatchedRank', 'mwl_matched_rank'), ('MWLMatchedTaxon', 'mwl_matched_taxon'),
            ('MWLTaxonomicScore', 'mwl_taxonomic_score'), ('MWLIdentity', 'mwl_identity'),
            ('MWLScore', 'mwl_score'), ('MWLRole', 'mwl_role'),
        ])
    fields.extend([
        ('SequencingSetID', 'sequencing_set_id'), ('SequencingSetRole', 'sequencing_set_role'),
        ('SequencingSetRank', 'sequencing_set_rank'),
        ('SelectionDiversityDistance', 'selection_diversity_distance'),
        ('SelectionGroupType', 'selection_group_type'),
        ('BaselineRedundancyStatus', 'baseline_redundancy_status'),
        ('BaselineRedundancyIdentityThreshold', 'baseline_redundancy_identity_threshold'),
        ('BaselineRedundancyMinQueryCoverage', 'baseline_redundancy_min_query_coverage'),
        ('BaselineExtensionStatus', 'baseline_extension_status'),
        ('BaselineExtensionMinIdentity', 'baseline_extension_min_identity'),
        ('BaselineExtensionMinQueryCoverage', 'baseline_extension_min_query_coverage'),
        ('SequencingSetReason', 'sequencing_set_reason'),
        ('SelectionGroupBasis', 'selection_group_basis'), ('SelectionGroupTaxon', 'selection_group_taxon'),
        ('NovelLookingClade', 'novel_looking_clade'),
        ('PhyloIsolation', 'phylo_isolation'), ('LocalNeighbourhoodFigure', 'local_neighbourhood_figure'),
        ('LocalPairwisePidentTable', 'local_pairwise_pident_table'),
        ('TreeContextLeafCount', 'tree_context_leaf_count'),
        ('AssessedSequencesInTreeContext', 'assessed_sequences_in_tree_context'),
        ('InTree', 'in_tree'), ('ClusterRepresentative', 'cluster_representative'),
        ('ClusterSize', 'cluster_size'), ('ClusteredMembers', 'clustered_members'),
        ('PlacementFlags', 'placement_flags'),
    ])
    for safe in alt_ref_safes:
        fields.extend([
            (f'Taxonomy_{safe}', f'alt_tax_{safe}'),
            (f'ClassificationHit_{safe}', f'alt_class_hit_{safe}'),
            (f'Identity_{safe}', f'alt_ident_{safe}'),
            (f'Confidence_{safe}', f'alt_conf_{safe}'),
        ])

    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, 'w', newline='') as fh:
        writer = csv.writer(fh, delimiter='\t', lineterminator='\n')
        writer.writerow([header for header, _ in fields])
        for row in rows:
            writer.writerow([row.get(key, 'NA') for _, key in fields])
    return str(p)


def _as_float(value, default=None):
    try:
        if value in (None, '', 'NA', 'None'):
            return default
        return float(value)
    except Exception:
        return default


def _as_int(value, default=None):
    try:
        if value in (None, '', 'NA', 'None'):
            return default
        return int(float(value))
    except Exception:
        return default


def _pool_available(row: dict, source_key: str) -> bool:
    return str(row.get(source_key) or '').lower() not in ('', 'na', 'none')


def _cultured_gap(row: dict) -> str:
    identity = _as_float(row.get('nearest_identity'))
    if not _pool_available(row, 'density_source') or identity is None:
        return 'UNKNOWN - no cultured baseline'
    if identity < 97.0:
        return 'LARGE - no match >=97%'
    if identity < 98.65:
        return 'MODERATE - nearest 97-98.65%'
    return 'SMALL - close match >=98.65%'


def _project_coverage(row: dict) -> str:
    if not _pool_available(row, 'project_density_source'):
        return 'UNKNOWN - no project collection'
    count = _as_int(row.get('project_matches_ge_97'), 0)
    if count == 0:
        return 'UNCOVERED - 0 neighbours >=97%'
    if count <= 3:
        return f'SPARSE - {count} neighbour(s) >=97%'
    if count <= 10:
        return f'MODERATE - {count} neighbours >=97%'
    return f'DENSE - {count} neighbours >=97%'


def _reference_context(row: dict) -> str:
    identity = _as_float(row.get('reference_nearest_identity'))
    if not _pool_available(row, 'reference_density_source') or identity is None:
        return 'UNKNOWN - no external reference comparison'
    if identity < 94.5:
        return 'HIGH DIVERGENCE - <94.5%'
    if identity < 98.65:
        return 'MODERATE DIVERGENCE - 94.5-98.65%'
    return 'CLOSE REFERENCE - >=98.65%'


def _genome_coverage(row: dict) -> str:
    if str(row.get('already_sequenced') or '').lower() == 'true':
        return 'CURRENT ISOLATE ALREADY SEQUENCED'
    if str(row.get('selected_for_genome_sequencing') or '').lower() == 'true':
        return 'CURRENT ISOLATE SELECTED - GENOME PENDING'
    committed = _as_int(row.get('genome_committed_count_same_species'))
    target = _as_int(row.get('pangenome_target'), 9)
    gap = _as_int(row.get('pangenome_gap'))
    if committed is not None and gap is not None:
        if gap <= 0:
            return f'TARGET MET - {committed}/{target} assessment-species genomes'
        return f'PANGENOME GAP - {committed}/{target} genomes committed; {gap} required'
    identity = _as_float(row.get('nearest_genome_identity'))
    if identity is None or identity <= 0:
        return 'NOT REPRESENTED'
    if identity >= 98.65:
        return f'NEAR-DUPLICATE AVAILABLE GENOME - {identity:.2f}%'
    if identity >= 97.0:
        return f'RELATED AVAILABLE GENOME - {identity:.2f}%'
    return f'NOT CLOSELY REPRESENTED - nearest {identity:.2f}%'


def _evidence_quality(row: dict) -> str:
    flags = {
        flag.strip().upper()
        for flag in str(row.get('placement_flags') or '').split(';')
        if flag.strip()
    }
    classification_identity = _as_float(row.get('classification_identity'))
    classification_confidence = _as_float(row.get('classification_confidence'))
    query_coverage = _as_float(row.get('classification_query_coverage'))
    taxonomy = str(row.get('taxonomy') or '').strip().lower()
    marker_flag = str(row.get('marker_qc_flag') or '').upper()
    if marker_flag in {'MARKER_QC_FAILED', 'MARKER_QC_REVIEW_REQUIRED', 'MARKER_QC_UNVERIFIED'}:
        return 'LOW'
    # MARKER_QC_REVIEW_APPROVED: manual approval has addressed QC concerns; fall through to
    # standard identity/coverage checks so clean assemblies can still reach HIGH.
    # Use exact flag membership to avoid CHIMERA_INDETERMINATE matching as CHIMERA.
    severe = bool(flags & {'NO_CLASSIFICATION', 'CHIMERA', 'CHIMERA_CONFIRMED', 'VERY_SHORT', 'HIGH_N_CONTENT'})
    if severe or (classification_identity is not None and classification_identity < 90.0):
        return 'LOW'
    if query_coverage is not None and query_coverage < 80.0:
        return 'LOW'
    if classification_identity is None and taxonomy in ('', 'na', 'none'):
        return 'LOW'
    moderate = bool(flags & {'LOW_CLASSIFICATION', 'LOW_CONFIDENCE', 'DISAGREEMENT', 'CONFLICT'})
    if (
        moderate
        or (classification_identity is not None and classification_identity < 95.0)
        or (classification_confidence is not None and classification_confidence < 0.8)
        or (query_coverage is not None and query_coverage < 90.0)
    ):
        return 'MODERATE'
    return 'HIGH'


def build_selection_decision(row: dict) -> dict:
    """Build transparent, categorical genome-sequencing recommendation evidence."""
    cultured_gap = _cultured_gap(row)
    project_coverage = _project_coverage(row)
    reference_context = _reference_context(row)
    genome_coverage = _genome_coverage(row)
    evidence_quality = _evidence_quality(row)
    selected = genome_coverage == 'CURRENT ISOLATE ALREADY SEQUENCED'
    pending = genome_coverage == 'CURRENT ISOLATE SELECTED - GENOME PENDING'
    pangenome_gap = _as_int(row.get('pangenome_gap'))
    target_met = pangenome_gap == 0
    set_role = str(row.get('sequencing_set_role') or '').upper()
    baseline_redundant = set_role == 'BASELINE_REDUNDANT'
    boundary_review = set_role == 'PANGENOME_BOUNDARY_REVIEW'
    mwl_match = str(row.get('mwl_match') or '').lower() == 'yes'
    mwl_rank = str(row.get('mwl_matched_rank') or '').lower()
    strong_mwl = mwl_match and mwl_rank in ('species', 's', 'genus', 'g', 'family', 'f')
    moderate_mwl = mwl_match and mwl_rank in ('order', 'o')
    baseline_identity = _as_float(row.get('nearest_identity'))
    project_identity = _as_float(row.get('project_nearest_identity'))
    reference_identity = _as_float(row.get('reference_nearest_identity'))

    positive = []
    caution = []
    if selected:
        caution.append('current isolate already has a genome available')
    elif pending:
        caution.append('current isolate is already selected; genome pending')
    if strong_mwl:
        matched = row.get('mwl_matched_taxon', 'NA')
        positive.append(f'MWL target match ({matched})')
    elif moderate_mwl:
        positive.append(f'MWL order-level context ({row.get("mwl_matched_taxon", "NA")})')
    if baseline_identity is not None and _pool_available(row, 'density_source') and baseline_identity < 97.0:
        positive.append(f'no close cultured isolate ({baseline_identity:.2f}%)')
    if project_identity is not None and _pool_available(row, 'project_density_source') and project_identity < 97.0:
        positive.append(f'no close project isolate ({project_identity:.2f}%)')
    if reference_identity is not None and _pool_available(row, 'reference_density_source') and reference_identity < 98.65:
        positive.append(f'divergent from external reference ({reference_identity:.2f}%)')
    if evidence_quality != 'HIGH':
        caution.append(f'{evidence_quality.lower()} marker evidence')
    if baseline_redundant:
        caution.append('near-identical to a cultured baseline marker at high query coverage')
    if boundary_review:
        status = str(row.get('baseline_extension_status') or 'extension criteria not met')
        caution.append(f'not admitted to the baseline-pangenome extension ({status})')
    if target_met:
        caution.append('assessment-species pangenome target already met')
    elif pangenome_gap is not None:
        positive.append(f'{pangenome_gap} genome(s) still required for pangenome target')
    if project_coverage.startswith('DENSE'):
        caution.append('dense prior project neighbourhood')
    if cultured_gap.startswith('SMALL'):
        caution.append('close cultured-baseline match')
    if cultured_gap.startswith('UNKNOWN') and project_coverage.startswith('UNKNOWN'):
        caution.append('baseline and prior-project comparisons unavailable')

    if selected:
        decision = 'ALREADY SEQUENCED'
    elif pending:
        decision = 'ALREADY SELECTED - GENOME PENDING'
    elif evidence_quality == 'LOW':
        decision = 'REVIEW BEFORE SELECTION'
    elif baseline_redundant:
        decision = 'EXCLUDE - BASELINE REDUNDANT'
    elif boundary_review:
        decision = 'REVIEW - PANGENOME BOUNDARY'
    elif set_role == 'PRIMARY':
        decision = 'PRIORITISE - SET PRIMARY'
    elif set_role == 'BACKUP':
        decision = 'RESERVE - SET BACKUP'
    elif set_role == 'DIVERSITY_CANDIDATE':
        decision = 'SECONDARY - STRAIN DIVERSITY'
    elif target_met:
        decision = 'LOWER PRIORITY - TARGET MET'
    elif set_role == 'ALTERNATE':
        decision = 'SECONDARY CANDIDATE'
    elif strong_mwl or len(positive) >= 2:
        decision = 'STRONG CANDIDATE'
    elif positive:
        decision = 'SECONDARY CANDIDATE'
    elif cultured_gap.startswith('SMALL') and project_coverage.startswith(('MODERATE', 'DENSE')):
        decision = 'LOWER PRIORITY - LIMITED ADDED VALUE'
    elif cultured_gap.startswith('UNKNOWN') and project_coverage.startswith('UNKNOWN'):
        decision = 'REVIEW BEFORE SELECTION'
    else:
        decision = 'SECONDARY CANDIDATE'

    reasons = positive[:3] + caution[:2]
    if not reasons:
        reasons.append('no strong novelty or redundancy signal')
    return {
        'decision': decision,
        'evidence_quality': evidence_quality,
        'cultured_gap': cultured_gap,
        'project_coverage': project_coverage,
        'reference_context': reference_context,
        'genome_coverage': genome_coverage,
        'decision_reason': '; '.join(reasons),
    }


def _board_recommendation(row: dict) -> tuple[str, str]:
    """Compatibility helper returning the transparent decision and rationale."""
    support = build_selection_decision(row)
    return support['decision'], support['decision_reason']


def write_selection_summary_tsv(path: str | Path, rows, assessment_db_name: str = 'GTDB'):
    """Write the decision-facing table while retaining auditable evidence."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    label = 'GTDB' if 'gtdb' in str(assessment_db_name).lower() else ''.join(
        char if char.isalnum() else '_' for char in str(assessment_db_name or 'Assessment')
    ).strip('_') or 'Assessment'
    headers = [
        'SequenceID',
        'PartnerID',
        'SelectedForGenomeSequencing',
        'GenomeAlreadySequenced',
        'Recommendation',
        'SequencingSetID',
        'SequencingSetRole',
        'SequencingSetRank',
        'SelectionDiversityDistance',
        'SelectionGroupType',
        'EvidenceQuality',
        'MarkerQC',
        'MarkerReview',
        f'{label}Taxonomy',
        f'{label}ClassificationIdentity',
        f'{label}QueryCoverage',
        'CulturedGap',
        'BaselineNearestHit',
        'BaselineNearestIdentity',
        'BaselineNearestQueryCoverage',
        'BaselineRedundancyStatus',
        'BaselineExtensionStatus',
        'BaselineTaxonomyAgreementRank',
        'BaselineTaxonomyConflict',
        'BaselineNoveltyScore',
        'ProjectCoverage',
        'ProjectNearestIdentity',
        'ProjectNearestQueryCoverage',
        'ProjectCloseNeighboursGE97',
        'ProjectNoveltyScore',
        f'{label}ReferenceContext',
        f'{label}ReferenceNearestIdentity',
        f'{label}ReferenceNearestQueryCoverage',
        'MWLMatchedRank',
        'MWLMatchedTaxon',
        'MWLScore',
        'BaselineGenomesSameSpecies',
        'SequencedPartnerGenomesSameSpecies',
        'SelectedPendingGenomesSameSpecies',
        'CommittedGenomesSameSpecies',
        'PangenomeTarget',
        'PangenomeGap',
        'GenomeCoverage',
        'LocalTreeFigure',
        'RecommendationReason',
    ]
    with open(p, 'w') as fh:
        fh.write('\t'.join(headers) + '\n')
        prepared = [(build_selection_decision(row), row) for row in rows]
        decision_order = {
            'PRIORITISE - SET PRIMARY': 0,
            'RESERVE - SET BACKUP': 1,
            'STRONG CANDIDATE': 2,
            'SECONDARY - STRAIN DIVERSITY': 3,
            'SECONDARY CANDIDATE': 4,
            'REVIEW BEFORE SELECTION': 5,
            'REVIEW - PANGENOME BOUNDARY': 6,
            'EXCLUDE - BASELINE REDUNDANT': 7,
            'LOWER PRIORITY - LIMITED ADDED VALUE': 8,
            'LOWER PRIORITY - TARGET MET': 9,
            'ALREADY SELECTED - GENOME PENDING': 10,
            'ALREADY SEQUENCED': 11,
        }
        prepared.sort(key=lambda item: (
            decision_order.get(item[0]['decision'], 99),
            str(item[1].get('partner_id', '')),
            str(item[1].get('id', '')),
        ))
        for support, row in prepared:
            values = [
                row.get('id', 'NA'),
                row.get('partner_id', 'NA'),
                row.get('selected_for_genome_sequencing', 'NA'),
                row.get('already_sequenced', 'NA'),
                support['decision'],
                row.get('sequencing_set_id', 'NA'),
                row.get('sequencing_set_role', 'NA'),
                row.get('sequencing_set_rank', 'NA'),
                row.get('selection_diversity_distance', 'NA'),
                row.get('selection_group_type', 'NA'),
                support['evidence_quality'],
                row.get('marker_qc_class', 'QUALITY_UNVERIFIED'),
                row.get('marker_manual_review_status', 'NOT_REVIEWED'),
                row.get('taxonomy', 'NA'),
                row.get('classification_identity', 'NA'),
                row.get('classification_query_coverage', 'NA'),
                support['cultured_gap'],
                row.get('nearest_hit', 'NA'),
                row.get('nearest_identity', 'NA'),
                row.get('nearest_query_coverage', 'NA'),
                row.get('baseline_redundancy_status', 'NA'),
                row.get('baseline_extension_status', 'NA'),
                row.get('baseline_taxonomy_agreement_rank', 'NA'),
                row.get('baseline_taxonomy_conflict', 'NA'),
                row.get('novelty_score', 'NA'),
                support['project_coverage'],
                row.get('project_nearest_identity', 'NA'),
                row.get('project_nearest_query_coverage', 'NA'),
                row.get('project_matches_ge_97', 'NA'),
                row.get('project_novelty_score', 'NA'),
                support['reference_context'],
                row.get('reference_nearest_identity', 'NA'),
                row.get('reference_nearest_query_coverage', 'NA'),
                row.get('mwl_matched_rank', 'NA'),
                row.get('mwl_matched_taxon', 'NA'),
                row.get('mwl_score', 'NA'),
                row.get('genome_available_count_same_species', 'NA'),
                row.get('genome_selected_count_same_species', 'NA'),
                row.get('genome_pending_count_same_species', 'NA'),
                row.get('genome_committed_count_same_species', 'NA'),
                row.get('pangenome_target', 'NA'),
                row.get('pangenome_gap', 'NA'),
                support['genome_coverage'],
                row.get('local_neighbourhood_figure', 'NA'),
                support['decision_reason'],
            ]
            safe = [str(v if v is not None else 'NA').replace('\t', ' ').replace('\n', ' ') for v in values]
            fh.write('\t'.join(safe) + '\n')
    return str(p)


def write_baseline_hits_tsv(path: str | Path, rows):
    """Write a concise nearest-baseline hit report for assessed sequences."""
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
