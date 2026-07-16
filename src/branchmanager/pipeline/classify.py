"""
classify.py — Taxonomic classification for BranchManager.

Uses VSEARCH --usearch_global against a reference database (GreenGenes2,
SILVA, or custom) to assign taxonomy to query sequences.

Key design decisions
--------------------
- Reference anchors (BRANCHMANAGER_REF_* headers) are excluded from output.
- The taxonomy lookup map is built once at load time for performance.
- vsearch is run at --id 0.8 minimum identity — this is for taxonomy only.
  Novelty detection (novelty.py) runs its own separate vsearch at a lower
  threshold to compute true nearest-neighbour distances.
- Output columns: ID  BestHit  Identity  Taxon  Confidence
  Confidence comes from the taxa TSV/CSV when provided (e.g. QIIME2 classifier
  score, 0–1 scale). When taxonomy is read directly from reference FASTA headers
  (no --taxa flag), Confidence is derived from the VSEARCH alignment identity
  (identity / 100), so a 98.5% identity hit yields 0.9850.
"""

from __future__ import annotations

import logging
import os
import re
import shlex
from pathlib import Path
from typing import Optional

from branchmanager.taxonomy import parse_reference_header_taxonomy, reference_lookup_keys
from branchmanager.taxonomy_io import iter_taxonomy_assignment_rows
from branchmanager.utils.fasta import read_fasta
from branchmanager.utils.subprocess import run_cmd
from branchmanager.pipeline.tree import is_ref_anchor

logger = logging.getLogger(__name__)


VSEARCH_USERFIELDS = 'query+target+id+alnlen+mism+gaps+qlo+qhi+tlo+thi+ql+tl'


def _parse_vsearch_match(parts: list[str]) -> Optional[dict]:
    if len(parts) < 3:
        return None
    try:
        identity = float(parts[2])
    except (TypeError, ValueError):
        return None
    result = {
        'query': parts[0], 'target': parts[1], 'identity': identity,
        'alignment_length': None, 'mismatches': None, 'gaps': None,
        'query_length': None, 'target_length': None,
        'query_coverage': None, 'target_coverage': None,
    }
    if len(parts) >= 12:
        try:
            result.update({
                'alignment_length': int(float(parts[3])),
                'mismatches': int(float(parts[4])),
                'gaps': int(float(parts[5])),
                'query_length': int(float(parts[10])),
                'target_length': int(float(parts[11])),
            })
            if result['query_length'] > 0:
                result['query_coverage'] = 100.0 * result['alignment_length'] / result['query_length']
            if result['target_length'] > 0:
                result['target_coverage'] = 100.0 * result['alignment_length'] / result['target_length']
        except (TypeError, ValueError, ZeroDivisionError):
            pass
    return result


def _alignment_columns(hit: Optional[dict]) -> list[str]:
    if not hit:
        return ['NA'] * 7
    def number(key, digits=2):
        value = hit.get(key)
        return 'NA' if value is None else f'{float(value):.{digits}f}'
    def integer(key):
        value = hit.get(key)
        return 'NA' if value is None else str(int(value))
    return [
        number('query_coverage'), number('target_coverage'), integer('alignment_length'),
        integer('query_length'), integer('target_length'), integer('mismatches'), integer('gaps'),
    ]


# ID normalisation helpers

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


# Taxonomy map loading

def _load_taxa_map(taxa_tsv: str) -> dict[str, tuple[str, Optional[float]]]:
    """
    Load a taxonomy TSV/CSV into a lookup dict.

    Accepts two formats:
      FeatureID  Taxon  Confidence       (QIIME2 / GreenGenes2 style)
      ID  Taxon  Confidence              (plain two-column)
      and equivalent comma-delimited CSV files, including .gz variants.

    Builds THREE index keys per entry for resilient lookup:
      - raw ID
      - _norm_id(ID)
      - _canon_id(ID)

    Returns {key: (taxon_string, confidence_float_or_None)}
    """
    taxa_map: dict[str, tuple[str, Optional[float]]] = {}

    try:
        # Log start of taxonomy loading
        try:
            size_mb = os.path.getsize(taxa_tsv) / (1024 * 1024)
            logger.info("[CLASSIFY] Loading taxonomy map from %s (%.1f MB)...", taxa_tsv, size_mb)
        except Exception:
            logger.info("[CLASSIFY] Loading taxonomy map from %s...", taxa_tsv)

        rows_loaded = 0
        for row in iter_taxonomy_assignment_rows(taxa_tsv):
            fid = str(row.get('id', '')).strip()
            tax = str(row.get('taxonomy', '')).strip()
            conf = row.get('confidence')
            _add_taxa_entry(fid, tax, conf, taxa_map)
            rows_loaded += 1

        logger.info("[CLASSIFY] Loaded %d taxonomy assignment rows from %s",
                    rows_loaded, taxa_tsv)
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
        for key in list(reference_lookup_keys(header)) + list(reference_lookup_keys(ref_id)):
            if key and key not in taxa_map:
                taxa_map[key] = (taxonomy, None)
    if total and parsed / float(total) < 0.8:
        warning_rows.insert(0, {
            'category': 'LOW_HEADER_TAXONOMY_PARSE_RATE',
            'subject_id': f'{parsed}/{total}',
            'detail': 'Fewer than 80% of reference FASTA headers contained parseable taxonomy. Check header format or provide --taxa explicitly.',
        })
    return taxa_map, warning_rows


def _add_taxa_entry(fid: str, tax: str, conf: object, taxa_map: dict) -> None:
    fid = str(fid or '').strip()
    tax = str(tax or '').strip()
    if not fid or not tax:
        return
    try:
        conf_value: Optional[float] = float(conf) if conf is not None else None
    except (ValueError, TypeError):
        conf_value = None

    for key in list(reference_lookup_keys(fid)) + [_norm_id(fid), _canon_id(fid)]:
        if key and key not in taxa_map:
            taxa_map[key] = (tax, conf_value)


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


# Reference decompression

def _ensure_uncompressed(ref_fasta: str, outdir: str, out_name: str = "ref_uncompressed.fasta") -> str:
    """Decompress *ref_fasta* to *outdir/<out_name>* if it is gzipped.

    ``out_name`` is exposed so callers can give each database its own cached
    filename and avoid different databases overwriting each other's decompressed
    copy.
    """
    if not str(ref_fasta).endswith(".gz"):
        return ref_fasta
    ref_unc = os.path.join(outdir, out_name)
    if os.path.exists(ref_unc):
        # File already exists - check if it's non-empty
        try:
            size = os.path.getsize(ref_unc)
            if size > 0:
                logger.info("[CLASSIFY] Using existing decompressed reference: %s (%.1f MB)", ref_unc, size / (1024*1024))
                return ref_unc
            else:
                logger.warning("[CLASSIFY] Existing decompressed file is empty; re-decompressing")
                os.remove(ref_unc)
        except Exception as e:
            logger.warning("[CLASSIFY] Could not check existing file: %s; re-decompressing", e)

    # Get compressed file size for progress estimation
    import time
    try:
        compressed_size = os.path.getsize(ref_fasta)
        logger.info("[CLASSIFY] Decompressing %s (%.1f MB compressed) → %s",
                   ref_fasta, compressed_size / (1024*1024), ref_unc)
        if compressed_size > 100 * 1024 * 1024:  # > 100 MB
            logger.warning("[CLASSIFY] Large reference file detected (%.1f MB compressed). "
                          "Decompression may take several minutes depending on disk I/O speed. "
                          "Please be patient...", compressed_size / (1024*1024))
    except Exception:
        logger.info("[CLASSIFY] Decompressing %s → %s", ref_fasta, ref_unc)

    import gzip as _gzip
    import shutil as _shutil
    # Decompress using fast byte-level streaming
    start_time = time.time()

    try:
        with _gzip.open(ref_fasta, 'rb') as _src, open(ref_unc, 'wb') as _dst:
            _shutil.copyfileobj(_src, _dst, length=1 << 20)  # 1 MiB buffer

        elapsed = time.time() - start_time
        final_size = os.path.getsize(ref_unc)
        final_mb = final_size / (1024 * 1024)
        logger.info("[CLASSIFY] Decompression complete: %.1f MB written in %.1f seconds (%.1f MB/s)",
                    final_mb, elapsed, final_mb / elapsed if elapsed > 0 else 0)
    except Exception as e:
        logger.error("[CLASSIFY] Decompression failed: %s", e)
        # Clean up partial file
        try:
            if os.path.exists(ref_unc):
                os.remove(ref_unc)
        except Exception:
            pass
        raise

    return ref_unc


def _iter_taxa_file_rows(taxa_tsv: str):
    for row in iter_taxonomy_assignment_rows(taxa_tsv):
        values = [str(row.get('id', '')), str(row.get('taxonomy', ''))]
        conf = row.get('confidence')
        if conf is not None:
            values.append(str(conf))
        yield values


def validate_reference_taxonomy_consistency(ref_fasta: Optional[str], taxa_tsv: Optional[str]):
    """Return warning rows about mismatches between reference FASTA and taxonomy table.

    Skips validation for large reference files (>500 MB) as it's too expensive.
    """
    warnings = []
    if not ref_fasta or not taxa_tsv:
        return warnings

    # Skip validation for large reference files - too expensive to parse
    try:
        ref_size = os.path.getsize(ref_fasta)
        if ref_size > 500 * 1024 * 1024:  # 500 MB threshold
            logger.info(
                '[CLASSIFY] Skipping reference/taxonomy consistency validation for large reference '
                'file (%.1f MB). This check is too expensive for production databases.',
                ref_size / (1024 * 1024)
            )
            return warnings
    except Exception:
        pass  # If we can't get file size, proceed with validation

    ref_ids = []
    ref_key_to_id = {}
    try:
        logger.info('[CLASSIFY] Validating reference/taxonomy consistency (this may take a minute for large files)...')
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
                'detail': f'Only {overlap:.1%} of reference FASTA IDs matched taxonomy table IDs; taxonomy source may be inconsistent with the classifier reference.',
            })
    for rid in unmatched_ref[:25]:
        warnings.append({'category': 'REFERENCE_ID_MISSING_TAXONOMY', 'subject_id': rid, 'detail': 'Reference FASTA entry had no matching taxonomy row.'})
    for tid in unmatched_taxa[:25]:
        warnings.append({'category': 'TAXONOMY_ID_MISSING_REFERENCE', 'subject_id': tid, 'detail': 'Taxonomy table entry had no matching reference FASTA record.'})
    return warnings


def write_reference_taxonomy_warnings(outdir: str, warning_rows):
    path = Path(outdir) / 'taxonomy_input_warnings.tsv'
    with open(path, 'w') as fh:
        fh.write('Category\tSubjectID\tDetail\n')
        for row in warning_rows:
            fh.write(f"{row.get('category', 'WARNING')}\t{row.get('subject_id', '')}\t{row.get('detail', '')}\n")
    return str(path)


# Main classification function

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
    taxa_tsv     : Optional TSV/CSV mapping reference IDs to taxonomy strings.
                   If omitted, only the best-hit accession is reported.
    threads      : VSEARCH thread count

    Output
    ------
    outdir/taxonomy.tsv  — TSV with columns:
        ID  BestHit  Identity  Taxon  Confidence

    Reference anchor sequences (BRANCHMANAGER_REF_*) are excluded from output.

    Returns
    -------
    Path to taxonomy.tsv
    """
    output = os.path.join(outdir, "taxonomy.tsv")

    if ref_fasta is None:
        raise ValueError("ref_fasta must be provided to run_classification")

    ref_to_use = _ensure_uncompressed(ref_fasta, outdir)

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

    matches_path = os.path.join(outdir, "matches.tsv")
    thread_flag = f" --threads {int(threads)}" if threads and int(threads) > 0 else ""

    cmd = (
        f"vsearch --usearch_global {shlex.quote(str(input_fasta))}"
        f" --db {shlex.quote(str(ref_to_use))}"
        f" --id 0.8"
        f" --strand both"
        f" --userout {shlex.quote(str(matches_path))}"
        f" --userfields {VSEARCH_USERFIELDS}"
        f" --maxaccepts 1 --maxhits 1"
        f" --query_cov 0.7"   # require 70% query coverage
        f"{thread_flag}"
    )
    logger.info("[CLASSIFY] Running vsearch classification")
    run_cmd(cmd)
    logger.info("[CLASSIFY] vsearch done → %s", matches_path)

    # Build a dict of all query sequences present in the input FASTA so we
    # can emit a row for every sequence (including those with no hit).
    all_query_ids = [
        h for h, _ in read_fasta(input_fasta)
        if not is_ref_anchor(h)
    ]

    # Parse best hits from matches file
    best_hits: dict[str, dict] = {}
    try:
        with open(matches_path) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                hit = _parse_vsearch_match(parts)
                if not hit:
                    continue
                # Keep highest identity hit (should only be one since maxhits=1)
                qid = hit['query']
                if qid not in best_hits or hit['identity'] > best_hits[qid]['identity']:
                    best_hits[qid] = hit
    except FileNotFoundError:
        logger.warning("[CLASSIFY] matches.tsv not found — classification may have failed")

    # Write output
    written = 0
    no_hit = 0
    with open(output, "w") as out_fh:
        out_fh.write("ID\tBestHit\tIdentity\tTaxon\tConfidence\tQueryCoverage\tTargetCoverage\tAlignmentLength\tQueryLength\tTargetLength\tMismatches\tGaps\n")
        for qid in all_query_ids:
            if qid not in best_hits:
                out_fh.write(f"{qid}\tNA\t0.0\tNA\tNA\t" + '\t'.join(_alignment_columns(None)) + "\n")
                no_hit += 1
                continue
            hit = best_hits[qid]
            sid, pct = hit['target'], hit['identity']
            tax, conf = _lookup_tax(sid, taxa_map)
            # Use alignment identity as confidence proxy when the taxa source
            # has no confidence column (e.g. taxonomy parsed from FASTA headers).
            if conf is None:
                conf = round(pct / 100.0, 4)
            out_fh.write(
                f"{qid}\t{sid}\t{pct:.2f}"
                f"\t{tax if tax is not None else 'NA'}"
                f"\t{conf:.4f}\t" + '\t'.join(_alignment_columns(hit)) + "\n"
            )
            written += 1

    logger.info(
        "[CLASSIFY] %d/%d sequences classified, %d with no hit → %s",
        written, len(all_query_ids), no_hit, output,
    )
    return output


# Multi-database classification

def _derive_db_name(fasta_path: str) -> str:
    """Derive a short, filesystem-safe name from a reference FASTA path.

    Strips common suffixes (_reps, _ssu, _16s, _rna, _seqs, _NR<digits>) before
    sanitising and truncating to 24 characters.
    """
    import re as _re
    base = Path(fasta_path).stem
    base = _re.sub(r'_(reps?|ssu_?reps?|16[sS]|rna|nr\d*|seqs?)$', '', base, flags=_re.IGNORECASE)
    safe = _re.sub(r'[^A-Za-z0-9]+', '_', base).strip('_')[:24]
    return safe if safe else Path(fasta_path).stem[:16]


def run_classification_single(
    input_fasta: str,
    outdir: str,
    ref_fasta: str,
    taxa_tsv: Optional[str] = None,
    threads: Optional[int] = None,
    db_name: str = 'main',
) -> tuple[str, dict]:
    """Classify *input_fasta* against a single reference database.

    Parameters
    ----------
    db_name : Label for this database.  When ``'main'``, outputs are written
              to ``taxonomy.tsv`` / ``matches.tsv`` (the canonical filenames
              expected by downstream code).  Any other name uses
              ``taxonomy_<db_name>.tsv`` / ``matches_<db_name>.tsv``.

    Returns
    -------
    (tsv_path, results_dict) where
    results_dict = {qid: (hit_id, pct_identity, taxon, confidence)}
    """
    import re as _re

    if db_name == 'main':
        tsv_name = 'taxonomy.tsv'
        match_name = 'matches.tsv'
    else:
        safe = _re.sub(r'[^A-Za-z0-9]+', '_', db_name).strip('_')
        tsv_name = f'taxonomy_{safe}.tsv'
        match_name = f'matches_{safe}.tsv'

    output = os.path.join(outdir, tsv_name)
    matches_path = os.path.join(outdir, match_name)

    # Use a db-specific decompressed filename so multiple gzipped references
    # used in the same run do not overwrite each other's cached copy.
    if db_name == 'main':
        decomp_name = 'ref_uncompressed.fasta'
    else:
        safe_decomp = _re.sub(r'[^A-Za-z0-9]+', '_', db_name).strip('_')
        decomp_name = f'ref_uncompressed_{safe_decomp}.fasta'
    ref_to_use = _ensure_uncompressed(ref_fasta, outdir, out_name=decomp_name)

    # Build taxa lookup map.
    # When a --alt-taxa TSV/CSV is supplied it is used directly; otherwise taxonomy
    # is parsed from the reference FASTA headers (GTDB-style lineage strings).
    taxa_map: dict = {}
    if taxa_tsv:
        logger.info("[CLASSIFY][%s] Loading taxa map from table: %s", db_name, taxa_tsv)
        taxa_map = _load_taxa_map(taxa_tsv)
    else:
        try:
            taxa_map, warning_rows = _load_taxa_map_from_reference_fasta(ref_to_use)
            if taxa_map:
                logger.info(
                    "[CLASSIFY][%s] Parsed taxonomy from %d reference FASTA headers",
                    db_name, len(taxa_map) // 3,
                )
            if warning_rows:
                warn_path = write_reference_taxonomy_warnings(outdir, warning_rows)
                logger.warning("[CLASSIFY][%s] Reference-header taxonomy warnings → %s", db_name, warn_path)
        except Exception as e:
            logger.warning("[CLASSIFY][%s] Could not parse taxa from FASTA headers: %s", db_name, e)

    # Run vsearch
    thread_flag = f" --threads {int(threads)}" if threads and int(threads) > 0 else ""
    cmd = (
        f"vsearch --usearch_global {shlex.quote(str(input_fasta))}"
        f" --db {shlex.quote(str(ref_to_use))}"
        f" --id 0.8"
        f" --strand both"
        f" --userout {shlex.quote(str(matches_path))}"
        f" --userfields {VSEARCH_USERFIELDS}"
        f" --maxaccepts 1 --maxhits 1"
        f" --query_cov 0.7"
        f"{thread_flag}"
    )
    logger.info("[CLASSIFY] Running vsearch for db=%s", db_name)
    run_cmd(cmd)
    logger.info("[CLASSIFY] vsearch done for db=%s → %s", db_name, matches_path)

    all_query_ids = [h for h, _ in read_fasta(input_fasta) if not is_ref_anchor(h)]

    # Parse best hits from matches file
    best_hits: dict[str, dict] = {}
    try:
        with open(matches_path) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                hit = _parse_vsearch_match(parts)
                if not hit:
                    continue
                qid = hit['query']
                if qid not in best_hits or hit['identity'] > best_hits[qid]['identity']:
                    best_hits[qid] = hit
    except FileNotFoundError:
        logger.warning("[CLASSIFY][%s] matches file not found", db_name)

    results: dict = {}
    written = 0
    no_hit = 0
    with open(output, "w") as out_fh:
        out_fh.write("ID\tBestHit\tIdentity\tTaxon\tConfidence\tQueryCoverage\tTargetCoverage\tAlignmentLength\tQueryLength\tTargetLength\tMismatches\tGaps\n")
        for qid in all_query_ids:
            if qid not in best_hits:
                out_fh.write(f"{qid}\tNA\t0.0\tNA\tNA\t" + '\t'.join(_alignment_columns(None)) + "\n")
                results[qid] = ('NA', 0.0, 'NA', 'NA', None, None, None, None, None, None, None)
                no_hit += 1
                continue
            hit = best_hits[qid]
            sid, pct = hit['target'], hit['identity']
            tax, conf = _lookup_tax(sid, taxa_map)
            if conf is None:
                conf = round(pct / 100.0, 4)
            tax_str = tax if tax is not None else 'NA'
            out_fh.write(
                f"{qid}\t{sid}\t{pct:.2f}\t{tax_str}\t{conf:.4f}\t"
                + '\t'.join(_alignment_columns(hit)) + "\n"
            )
            results[qid] = (
                sid, pct, tax_str, conf, hit.get('query_coverage'), hit.get('target_coverage'),
                hit.get('alignment_length'), hit.get('query_length'), hit.get('target_length'),
                hit.get('mismatches'), hit.get('gaps'),
            )
            written += 1

    logger.info(
        "[CLASSIFY][%s] %d/%d sequences classified, %d with no hit → %s",
        db_name, written, len(all_query_ids), no_hit, output,
    )
    return output, results


def run_all_classifications(
    input_fasta: str,
    outdir: str,
    primary_ref: str,
    primary_taxa: Optional[str] = None,
    primary_name: str = 'main',
    alt_refs: Optional[list] = None,
    threads: Optional[int] = None,
    main_db: Optional[str] = None,
) -> tuple[str, dict]:
    """Classify *input_fasta* against the primary and all alternative databases.

    Parameters
    ----------
    primary_ref   : Path to the primary reference FASTA.
    primary_taxa  : Optional taxa TSV/CSV for the primary reference.
    primary_name  : Display name for the primary database (default ``'main'``).
    alt_refs      : List of ``(ref_fasta, taxa_table_or_None, db_name)`` tuples.
    main_db       : Which database (by name) should be written to the canonical
                    ``taxonomy.tsv``.  Defaults to *primary_name*.

    Returns
    -------
    (primary_taxonomy_tsv, all_results)
    where ``all_results`` = ``{db_name: {qid: (hit, pct, tax, conf)}}``.

    Side effects
    ------------
    * ``taxonomy.tsv``                   - main DB
    * ``taxonomy_<name>.tsv``            – one file per additional DB
    * ``taxonomy_all_dbs.tsv``           – wide merged table (multi-DB only)
    """
    effective_main = main_db or primary_name

    all_dbs = [(primary_ref, primary_taxa, primary_name)] + (alt_refs or [])
    all_results: dict = {}

    for ref, taxa, name in all_dbs:
        # The database whose name matches effective_main writes to taxonomy.tsv
        db_label = 'main' if name == effective_main else name
        _, results = run_classification_single(input_fasta, outdir, ref, taxa, threads, db_label)
        all_results[name] = results

    # Write wide merged taxonomy when more than one database was used
    if len(all_dbs) > 1:
        all_query_ids = [h for h, _ in read_fasta(input_fasta) if not is_ref_anchor(h)]
        _write_merged_taxonomy(outdir, all_results, effective_main, all_query_ids)

    primary_tsv = os.path.join(outdir, 'taxonomy.tsv')
    return primary_tsv, all_results


def _write_merged_taxonomy(
    outdir: str,
    all_results: dict,
    primary_name: str,
    all_query_ids: list,
) -> str:
    """Write ``taxonomy_all_dbs.tsv`` — a wide table with one set of columns per DB.

    The primary database's columns appear first (``BestHit``, ``Identity``,
    ``Taxon``, ``Confidence``).  Each additional database adds a suffixed set
    (``BestHit_<DB>``, ``Identity_<DB>``, ``Taxon_<DB>``, ``Confidence_<DB>``).
    """
    import re as _re
    output = os.path.join(outdir, 'taxonomy_all_dbs.tsv')

    db_names = [primary_name] + [n for n in all_results if n != primary_name]

    with open(output, 'w') as fh:
        metric_names = ['QueryCoverage', 'TargetCoverage', 'AlignmentLength', 'QueryLength', 'TargetLength', 'Mismatches', 'Gaps']
        header = ['ID', 'BestHit', 'Identity', 'Taxon', 'Confidence', *metric_names]
        for name in db_names[1:]:
            safe = _re.sub(r'[^A-Za-z0-9]+', '_', name).strip('_')
            header += [f'BestHit_{safe}', f'Identity_{safe}', f'Taxon_{safe}', f'Confidence_{safe}']
            header += [f'{metric}_{safe}' for metric in metric_names]
        fh.write('\t'.join(header) + '\n')

        primary_res = all_results.get(primary_name, {})
        for qid in all_query_ids:
            result = primary_res.get(qid, ('NA', 0.0, 'NA', 'NA', None, None, None, None, None, None, None))
            hit, pct, tax, conf = result[:4]
            row = [
                qid,
                str(hit),
                f"{pct:.2f}" if isinstance(pct, float) else str(pct),
                str(tax),
                f"{conf:.4f}" if isinstance(conf, float) else str(conf),
            ]
            row += ['NA' if value is None else (f'{value:.2f}' if isinstance(value, float) else str(value)) for value in result[4:11]]
            for name in db_names[1:]:
                alt_result = all_results.get(name, {}).get(qid, ('NA', 0.0, 'NA', 'NA', None, None, None, None, None, None, None))
                a_hit, a_pct, a_tax, a_conf = alt_result[:4]
                row += [
                    str(a_hit),
                    f"{a_pct:.2f}" if isinstance(a_pct, float) else str(a_pct),
                    str(a_tax),
                    f"{a_conf:.4f}" if isinstance(a_conf, float) else str(a_conf),
                ]
                row += ['NA' if value is None else (f'{value:.2f}' if isinstance(value, float) else str(value)) for value in alt_result[4:11]]
            fh.write('\t'.join(row) + '\n')

    logger.info("[CLASSIFY] Wrote merged multi-db taxonomy → %s", output)
    return output
