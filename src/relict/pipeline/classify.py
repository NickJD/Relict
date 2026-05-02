"""
classify.py — Taxonomic classification for Relict.

Uses VSEARCH --usearch_global against a reference database (GreenGenes2,
SILVA, or custom) to assign taxonomy to query sequences.

Key design decisions
--------------------
- Reference anchors (RELICT_REF_* headers) are excluded from output.
- The taxonomy lookup map is built once at load time for performance.
- vsearch is run at --id 0.8 minimum identity — this is for taxonomy only.
  Novelty detection (novelty.py) runs its own separate vsearch at a lower
  threshold to compute true nearest-neighbour distances.
- Output columns: ID  BestHit  Identity  Taxon  Confidence
  Confidence comes from the taxa TSV when provided (e.g. QIIME2 classifier score,
  0–1 scale).  When taxonomy is read directly from reference FASTA headers (no
  --taxa flag), Confidence is derived from the VSEARCH alignment identity
  (identity / 100), so a 98.5% identity hit yields 0.9850.
"""

from __future__ import annotations

import gzip
import logging
import os
import re
from pathlib import Path
from typing import Optional

from relict.taxonomy import parse_reference_header_taxonomy, reference_lookup_keys
from relict.utils.fasta import read_fasta, write_fasta
from relict.utils.subprocess import run_cmd
from relict.pipeline.tree import is_ref_anchor

logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# ID normalisation helpers
# ─────────────────────────────────────────────────────────────────────────────

def _norm_id(x: str) -> str:
    """
    Light normalisation: take first whitespace token, prefer last pipe field,
    strip coordinate suffixes, collapse non-alnum to underscore.
    """
    if not x:
        return ""
    x = str(x).split()[0]
    if "|" in x:
        x = x.split("|")[-1]
    if "#" in x:
        x = x.split("#")[0]
    if ":" in x:
        x = x.split(":")[0]
    x = x.strip("()[]{}").strip()
    x = re.sub(r"[^0-9A-Za-z]+", "_", x).strip("_")
    return x


def _canon_id(x: str) -> str:
    """
    Aggressive canonicalisation: remove node fragments, coordinates,
    prefer last pipe field, collapse to underscore.
    """
    if not x:
        return ""
    y = str(x).split()[0]
    if "#" in y:
        y = y.split("#")[0]
    y = re.sub(r":\d+-\d+\(.*\)$", "", y)
    y = re.sub(r":\d+-\d+$", "", y)
    if "|" in y:
        y = y.split("|")[-1]
    y = re.sub(r"[^0-9A-Za-z]+", "_", y).strip("_")
    return y


# ─────────────────────────────────────────────────────────────────────────────
# Taxonomy map loading
# ─────────────────────────────────────────────────────────────────────────────

def _load_taxa_map(taxa_tsv: str) -> dict[str, tuple[str, Optional[float]]]:
    """
    Load a taxonomy TSV into a lookup dict.

    Accepts two formats:
      FeatureID  Taxon  Confidence       (QIIME2 / GreenGenes2 style)
      ID  Taxon  Confidence              (plain two-column)

    Builds THREE index keys per entry for resilient lookup:
      - raw ID
      - _norm_id(ID)
      - _canon_id(ID)

    Returns {key: (taxon_string, confidence_float_or_None)}
    """
    taxa_map: dict[str, tuple[str, Optional[float]]] = {}

    try:
        if str(taxa_tsv).endswith('.gz'):
            fh_ctx = gzip.open(taxa_tsv, 'rt')
        else:
            fh_ctx = open(taxa_tsv, 'rt')
        with fh_ctx as fh:
            first_line = fh.readline().strip()
            # Detect and skip header row
            is_header = any(
                kw in first_line.lower()
                for kw in ("feature", "taxon", "taxonomy", "id\t")
            )
            if not is_header:
                # First line is data — process it
                _add_taxa_row(first_line, taxa_map)

            for line in fh:
                line = line.strip()
                if line:
                    _add_taxa_row(line, taxa_map)

        logger.info("[CLASSIFY] Loaded %d unique taxa IDs from %s",
                    len(taxa_map) // 3, taxa_tsv)
    except Exception as e:
        logger.warning("[CLASSIFY] Failed to load taxa file %s: %s", taxa_tsv, e)

    return taxa_map


def _load_taxa_map_from_reference_fasta(ref_fasta: str) -> tuple[dict[str, tuple[str, Optional[float]]], list[dict[str, str]]]:
    """Build taxonomy lookups directly from reference FASTA headers.

    Supports GTDB 16S headers where the accession is followed by a lineage.
    Returns (taxa_map, warning_rows).
    """
    taxa_map: dict[str, tuple[str, Optional[float]]] = {}
    warning_rows = []
    total = 0
    parsed = 0
    for header, _ in read_fasta(ref_fasta):
        total += 1
        ref_id, taxonomy = parse_reference_header_taxonomy(header)
        if not ref_id or not taxonomy:
            warning_rows.append({
                'category': 'REFERENCE_HEADER_WITHOUT_TAXONOMY',
                'subject_id': ref_id or str(header).split()[0],
                'detail': 'Reference FASTA header did not contain a parseable lineage string.',
            })
            continue
        parsed += 1
        for key in reference_lookup_keys(ref_id):
            if key and key not in taxa_map:
                taxa_map[key] = (taxonomy, None)
    if total and parsed / float(total) < 0.8:
        warning_rows.insert(0, {
            'category': 'LOW_HEADER_TAXONOMY_PARSE_RATE',
            'subject_id': f'{parsed}/{total}',
            'detail': 'Fewer than 80% of reference FASTA headers contained parseable taxonomy. Check header format or provide --taxa explicitly.',
        })
    return taxa_map, warning_rows


def _add_taxa_row(line: str, taxa_map: dict) -> None:
    parts = line.split("\t")
    if len(parts) < 2:
        return
    fid = parts[0].strip()
    tax = parts[1].strip()
    try:
        conf: Optional[float] = float(parts[2]) if len(parts) > 2 else None
    except (ValueError, IndexError):
        conf = None

    for key in list(reference_lookup_keys(fid)) + [_norm_id(fid), _canon_id(fid)]:
        if key and key not in taxa_map:
            taxa_map[key] = (tax, conf)


def _lookup_tax(sid: str, taxa_map: dict) -> tuple[Optional[str], Optional[float]]:
    """
    Look up a VSEARCH hit ID in the taxa_map.
    Tries: raw → norm → canon → token-level splitting.
    Returns (taxon, confidence) or (None, None).
    """
    if not taxa_map:
        return None, None

    for key in list(reference_lookup_keys(sid)) + [_norm_id(sid), _canon_id(sid)]:
        if key and key in taxa_map:
            return taxa_map[key]

    # Token-level splitting for complex IDs like RS_GCF_000001405.40
    tokens = re.split(r"[|_:\-.]+", sid)
    for tok in tokens:
        if not tok:
            continue
        for key in (tok, _norm_id(tok), _canon_id(tok)):
            if key and key in taxa_map:
                return taxa_map[key]

    return None, None


# ─────────────────────────────────────────────────────────────────────────────
# Reference decompression
# ─────────────────────────────────────────────────────────────────────────────

def _ensure_uncompressed(ref_fasta: str, outdir: str) -> str:
    if not str(ref_fasta).endswith(".gz"):
        return ref_fasta
    ref_unc = os.path.join(outdir, "ref_uncompressed.fasta")
    if not os.path.exists(ref_unc):
        logger.info("[CLASSIFY] Decompressing %s → %s", ref_fasta, ref_unc)
        records = list(read_fasta(ref_fasta))
        write_fasta(records, ref_unc)
    return ref_unc


def _iter_taxa_file_rows(taxa_tsv: str):
    if str(taxa_tsv).endswith('.gz'):
        fh_ctx = gzip.open(taxa_tsv, 'rt')
    else:
        fh_ctx = open(taxa_tsv, 'rt')
    with fh_ctx as fh:
        first_line = fh.readline().strip()
        is_header = any(kw in first_line.lower() for kw in ('feature', 'taxon', 'taxonomy', 'id\t'))
        if not is_header and first_line:
            yield first_line.split('\t')
        for line in fh:
            line = line.strip()
            if line:
                yield line.split('\t')


def validate_reference_taxonomy_consistency(ref_fasta: Optional[str], taxa_tsv: Optional[str]):
    """Return warning rows about mismatches between reference FASTA and taxonomy TSV."""
    warnings = []
    if not ref_fasta or not taxa_tsv:
        return warnings

    ref_ids = []
    ref_key_to_id = {}
    try:
        for header, _ in read_fasta(ref_fasta):
            ref_ids.append(header)
            for key in reference_lookup_keys(header):
                ref_key_to_id.setdefault(key, header)
    except Exception as e:
        return [{'category': 'REFERENCE_READ_FAILED', 'subject_id': str(ref_fasta), 'detail': str(e)}]

    taxa_ids = []
    taxa_key_to_id = {}
    try:
        for parts in _iter_taxa_file_rows(taxa_tsv):
            if len(parts) < 2:
                continue
            tid = parts[0].strip()
            if not tid:
                continue
            taxa_ids.append(tid)
            for key in reference_lookup_keys(tid):
                taxa_key_to_id.setdefault(key, tid)
    except Exception as e:
        return [{'category': 'TAXONOMY_READ_FAILED', 'subject_id': str(taxa_tsv), 'detail': str(e)}]

    matched_ref = [rid for rid in ref_ids if any(key in taxa_key_to_id for key in reference_lookup_keys(rid))]
    unmatched_ref = [rid for rid in ref_ids if rid not in matched_ref]
    unmatched_taxa = [tid for tid in taxa_ids if not any(key in ref_key_to_id for key in reference_lookup_keys(tid))]

    if ref_ids:
        overlap = len(matched_ref) / float(len(ref_ids))
        if overlap < 0.8:
            warnings.append({
                'category': 'LOW_REFERENCE_TAXONOMY_OVERLAP',
                'subject_id': f'{len(matched_ref)}/{len(ref_ids)}',
                'detail': f'Only {overlap:.1%} of reference FASTA IDs matched taxonomy TSV IDs; taxonomy source may be inconsistent with the classifier reference.',
            })
    for rid in unmatched_ref[:25]:
        warnings.append({'category': 'REFERENCE_ID_MISSING_TAXONOMY', 'subject_id': rid, 'detail': 'Reference FASTA entry had no matching taxonomy row.'})
    for tid in unmatched_taxa[:25]:
        warnings.append({'category': 'TAXONOMY_ID_MISSING_REFERENCE', 'subject_id': tid, 'detail': 'Taxonomy TSV entry had no matching reference FASTA record.'})
    return warnings


def write_reference_taxonomy_warnings(outdir: str, warning_rows):
    path = Path(outdir) / 'taxonomy_input_warnings.tsv'
    with open(path, 'w') as fh:
        fh.write('Category\tSubjectID\tDetail\n')
        for row in warning_rows:
            fh.write(f"{row.get('category', 'WARNING')}\t{row.get('subject_id', '')}\t{row.get('detail', '')}\n")
    return str(path)


# ─────────────────────────────────────────────────────────────────────────────
# Main classification function
# ─────────────────────────────────────────────────────────────────────────────

def run_classification(
    input_fasta: str,
    outdir: str,
    ref_fasta: Optional[str] = None,
    taxa_tsv: Optional[str] = None,
    threads: Optional[int] = None,
) -> str:
    """
    Classify sequences in input_fasta against ref_fasta using VSEARCH.

    Parameters
    ----------
    input_fasta  : Query FASTA (QC-filtered, dereplicated)
    outdir       : Output directory
    ref_fasta    : Reference FASTA database for VSEARCH search
    taxa_tsv     : Optional TSV mapping reference IDs to taxonomy strings.
                   If omitted, only the best-hit accession is reported.
    threads      : VSEARCH thread count

    Output
    ------
    outdir/taxonomy.tsv  — TSV with columns:
        ID  BestHit  Identity  Taxon  Confidence

    Reference anchor sequences (RELICT_REF_*) are excluded from output.

    Returns
    -------
    Path to taxonomy.tsv
    """
    output = os.path.join(outdir, "taxonomy.tsv")

    if ref_fasta is None:
        raise ValueError("ref_fasta must be provided to run_classification")

    ref_to_use = _ensure_uncompressed(ref_fasta, outdir)

    # ── Build taxa lookup map once ────────────────────────────────────────────
    taxa_map: dict = {}
    if taxa_tsv:
        try:
            warning_rows = validate_reference_taxonomy_consistency(ref_to_use, taxa_tsv)
            if warning_rows:
                warn_path = write_reference_taxonomy_warnings(outdir, warning_rows)
                logger.warning('[CLASSIFY] Reference/taxonomy consistency warnings written to %s', warn_path)
        except Exception as e:
            logger.warning('[CLASSIFY] Failed to validate reference/taxonomy consistency: %s', e)
        taxa_map = _load_taxa_map(taxa_tsv)
    else:
        try:
            taxa_map, warning_rows = _load_taxa_map_from_reference_fasta(ref_to_use)
            if warning_rows:
                warn_path = write_reference_taxonomy_warnings(outdir, warning_rows)
                logger.warning('[CLASSIFY] Reference-header taxonomy warnings written to %s', warn_path)
            if taxa_map:
                logger.info('[CLASSIFY] Using taxonomy parsed directly from reference FASTA headers')
        except Exception as e:
            logger.warning('[CLASSIFY] Failed to parse taxonomy from reference FASTA headers: %s', e)
            taxa_map = {}

    # ── Run VSEARCH ───────────────────────────────────────────────────────────
    matches_path = os.path.join(outdir, "matches.tsv")
    thread_flag = f" --threads {int(threads)}" if threads and int(threads) > 0 else ""

    cmd = (
        f"vsearch --usearch_global {input_fasta}"
        f" --db {ref_to_use}"
        f" --id 0.8"
        f" --strand both"
        f" --blast6out {matches_path}"
        f" --maxaccepts 1 --maxhits 1"
        f" --query_cov 0.7"   # require 70% query coverage
        f"{thread_flag}"
    )
    logger.info("[CLASSIFY] Running vsearch classification")
    run_cmd(cmd)
    logger.info("[CLASSIFY] vsearch done → %s", matches_path)

    # ── Parse matches and write taxonomy table ────────────────────────────────
    # Build a dict of all query sequences present in the input FASTA so we
    # can emit a row for every sequence (including those with no hit).
    all_query_ids = [
        h for h, _ in read_fasta(input_fasta)
        if not is_ref_anchor(h)
    ]

    # Parse best hits from matches file
    best_hits: dict[str, tuple[str, float]] = {}  # query_id -> (hit_id, pct_identity)
    try:
        with open(matches_path) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                qid, sid = parts[0], parts[1]
                try:
                    pct = float(parts[2])
                except ValueError:
                    continue
                # Keep highest identity hit (should only be one since maxhits=1)
                if qid not in best_hits or pct > best_hits[qid][1]:
                    best_hits[qid] = (sid, pct)
    except FileNotFoundError:
        logger.warning("[CLASSIFY] matches.tsv not found — classification may have failed")

    # Write output
    written = 0
    no_hit = 0
    with open(output, "w") as out_fh:
        out_fh.write("ID\tBestHit\tIdentity\tTaxon\tConfidence\n")
        for qid in all_query_ids:
            if qid not in best_hits:
                out_fh.write(f"{qid}\tNA\t0.0\tNA\tNA\n")
                no_hit += 1
                continue
            sid, pct = best_hits[qid]
            tax, conf = _lookup_tax(sid, taxa_map)
            # Use alignment identity as confidence proxy when the taxa source
            # has no confidence column (e.g. taxonomy parsed from FASTA headers).
            if conf is None:
                conf = round(pct / 100.0, 4)
            out_fh.write(
                f"{qid}\t{sid}\t{pct:.2f}"
                f"\t{tax if tax is not None else 'NA'}"
                f"\t{conf:.4f}\n"
            )
            written += 1

    logger.info(
        "[CLASSIFY] %d/%d sequences classified, %d with no hit → %s",
        written, len(all_query_ids), no_hit, output,
    )
    return output