"""
novelty.py — Novelty detection for BranchManager.

Novelty is expressed as the identity percentage of the best hit against
the reference database. This gives a gradient from 0 (no hit at all, or
very distant) to 100 (identical to a known sequence).

The distinction from classify.py is important:
  - classify.py runs vsearch at --id 0.8 to find the best taxonomic hit
  - novelty.py runs vsearch at --id 0.0 (no threshold) to find the NEAREST
    neighbour regardless of distance, so we always get a real distance metric

The two vsearch runs serve different purposes and MUST NOT share matches.tsv.
classify.py  → outdir/matches.tsv         (taxonomy, min 80% identity)
novelty.py   → outdir/novelty_matches.tsv (distances, no threshold)

Output columns: ID\tNearestIdentity\tNearestHit\tNovel
  NearestIdentity : float 0–100, identity to best reference hit
                    0.0 means no hit found at any identity level
  NearestHit      : accession of best reference sequence
  Novel           : True if NearestIdentity < id_threshold * 100
"""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Optional

from branchmanager.pipeline import tree as tree_pipeline
from branchmanager.utils.fasta import read_fasta, write_fasta
from branchmanager.utils.subprocess import run_cmd

logger = logging.getLogger(__name__)


def run_novelty(
    input_fasta: str,
    ref_fasta: str,
    outdir: str,
    db=None,
    run_dataset: Optional[str] = None,
    id_threshold: float = 0.97,
    threads: Optional[int] = None,
    target_fasta: Optional[str] = None,
) -> str:
    """
    Compute novelty scores for all sequences in input_fasta.

    Novelty is measured against previously submitted sequences (preload + other
    runs already in the DB).  This means a sequence is "novel" if it is more
    than ``1 - id_threshold`` different from everything the user has already
    submitted to BranchManager, not from a global reference.

    If ``target_fasta`` is given, novelty is measured against those sequences
    instead of (or as a supplement to) the DB-derived sequences.

    Parameters
    ----------
    input_fasta    : Query sequences (dereplicated, QC-passed)
    ref_fasta      : Reference database FASTA (kept for API compatibility; not
                     used directly — novelty is against submitted sequences)
    outdir         : Output directory
    db             : Database instance to pull other submitted sequences from
    run_dataset    : Current run dataset name (used to exclude self from comparison)
    id_threshold   : Sequences with identity below this are marked Novel=True
                     (default 0.97 = 97% identity, classic species boundary)
    threads        : VSEARCH thread count
    target_fasta   : Explicit FASTA of target sequences; overrides DB lookup

    Returns
    -------
    Path to novelty.tsv
    """
    output = f"{outdir}/novelty.tsv"
    matches_path = f"{outdir}/novelty_matches.tsv"

    # ── Build comparison database ─────────────────────────────────────────────
    # Priority: explicit target_fasta > DB-derived submitted sequences > nothing
    db_fasta = Path(outdir) / "submitted_sequences.fasta"
    submission_db_exists = False

    if target_fasta and Path(target_fasta).exists():
        # Explicit target sequences override DB lookup
        db_fasta = Path(target_fasta)
        submission_db_exists = True
        logger.info("[NOVELTY] Using explicit --target FASTA for novelty comparison: %s", target_fasta)

    elif db is not None:
        try:
            with db.connect() as conn:
                cur = conn.cursor()
                # Get sequences from ALL datasets except the current run.
                # This intentionally includes every preload dataset so that
                # novelty is measured against the full user-submitted baseline.
                if run_dataset:
                    cur.execute(
                        "SELECT id, sequence FROM sequences WHERE dataset IS NOT ? OR dataset IS NULL",
                        (run_dataset,)
                    )
                else:
                    cur.execute("SELECT id, sequence FROM sequences")
                rows = cur.fetchall()

            # Log the distinct datasets contributing to the comparison pool
            try:
                with db.connect() as conn:
                    cur = conn.cursor()
                    if run_dataset:
                        cur.execute(
                            "SELECT DISTINCT dataset FROM sequences "
                            "WHERE dataset IS NOT ? OR dataset IS NULL",
                            (run_dataset,),
                        )
                    else:
                        cur.execute("SELECT DISTINCT dataset FROM sequences")
                    pool_datasets = [r[0] for r in cur.fetchall()]
                logger.info(
                    "[NOVELTY] Comparison pool contains %d distinct dataset(s): %s",
                    len(pool_datasets),
                    pool_datasets,
                )
            except Exception:
                pass

            records = [(r[0], r[1]) for r in rows if r[1]]
            if records:
                write_fasta(records, str(db_fasta))
                submission_db_exists = True
                logger.info(
                    "[NOVELTY] Built submission comparison DB from %d previously submitted sequences "
                    "(preload + prior run datasets). Novelty is expressed relative to YOUR submitted "
                    "sequences, not against the global reference.",
                    len(records),
                )
        except Exception as e:
            logger.warning("[NOVELTY] Failed to build submission database: %s", e)

    # Run vsearch against submission database if available
    if submission_db_exists and db_fasta.exists():
        thread_flag = f" --threads {int(threads)}" if threads and int(threads) > 0 else ""
        cmd = (
            f"vsearch --usearch_global {input_fasta} "
            f"--db {str(db_fasta)} "
            f"--id 0.5 "                          # low floor — we want real distances
            f"--strand both "
            f"--blast6out {matches_path} "
            f"--maxaccepts 1 --maxhits 1 "
            f"--query_cov 0.7"                    # require 70% query coverage to avoid spurious hits
            f"{thread_flag}"
        )
        logger.info("[NOVELTY] Running vsearch against previously submitted sequences")
        try:
            run_cmd(cmd)
            logger.info("[NOVELTY] vsearch done → %s", matches_path)
        except Exception as e:
            logger.warning("[NOVELTY] vsearch against submission DB failed: %s", e)
    else:
        logger.info(
            "[NOVELTY] No previously submitted sequences available and no --target FASTA provided; "
            "all input sequences will be treated as novel. Consider running `branchmanager preload` first "
            "to build a baseline database."
        )

    # Parse best hits: query_id -> (pct_identity, hit_id)
    best_hits: dict[str, tuple[float, str]] = _parse_best_hits(matches_path)

    novel_threshold_pct = id_threshold * 100.0

    with open(output, "w") as fh:
        fh.write("ID\tNearestIdentity\tNearestHit\tNovel\n")
        for header, _ in read_fasta(input_fasta):
            # Skip reference anchors — they are scaffolding, not results
            if tree_pipeline.is_ref_anchor(header):
                continue
            pct_id, hit_id = best_hits.get(header, (0.0, "None"))
            is_novel = pct_id < novel_threshold_pct
            fh.write(f"{header}\t{pct_id:.2f}\t{hit_id}\t{is_novel}\n")

    n_total = sum(1 for _ in read_fasta(input_fasta) if not tree_pipeline.is_ref_anchor(_[0]))
    n_novel = sum(1 for v in best_hits.values() if v[0] < novel_threshold_pct)
    n_unmatched = n_total - len(best_hits)

    logger.info(
        "[NOVELTY] %d sequences: %d novel (<%d%%), %d unmatched (treated as novel), %d known (by submission history)",
        n_total,
        n_novel + n_unmatched,
        int(novel_threshold_pct),
        n_unmatched,
        n_total - n_novel - n_unmatched,
    )
    return output


def _parse_best_hits(matches_path: str) -> dict[str, tuple[float, str]]:
    """
    Parse VSEARCH blast6 output.
    Returns {query_id: (pct_identity_float, hit_id_str)}.
    Only the best (highest identity) hit per query is kept.
    """
    best: dict[str, tuple[float, str]] = {}
    try:
        with open(matches_path) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                qid = parts[0]
                hit_id = parts[1]
                try:
                    pct_id = float(parts[2])
                except ValueError:
                    continue
                if qid not in best or pct_id > best[qid][0]:
                    best[qid] = (pct_id, hit_id)
    except FileNotFoundError:
        logger.warning("[NOVELTY] matches file not found: %s", matches_path)
    return best


def _ensure_uncompressed(ref_fasta: str, outdir: str) -> str:
    """Decompress gzipped reference fasta to outdir if needed. Returns path to use."""
    if not str(ref_fasta).endswith(".gz"):
        return ref_fasta
    ref_unc = os.path.join(outdir, "ref_uncompressed.fasta")
    if os.path.exists(ref_unc):
        # File already exists - check if it's non-empty
        try:
            size = os.path.getsize(ref_unc)
            if size > 0:
                logger.info("[NOVELTY] Using existing decompressed reference: %s (%.1f MB)", ref_unc, size / (1024*1024))
                return ref_unc
            else:
                logger.warning("[NOVELTY] Existing decompressed file is empty; re-decompressing")
                os.remove(ref_unc)
        except Exception as e:
            logger.warning("[NOVELTY] Could not check existing file: %s; re-decompressing", e)

    # Get compressed file size for progress estimation
    import time
    try:
        compressed_size = os.path.getsize(ref_fasta)
        logger.info("[NOVELTY] Decompressing %s (%.1f MB compressed) → %s",
                   ref_fasta, compressed_size / (1024*1024), ref_unc)
        if compressed_size > 100 * 1024 * 1024:  # > 100 MB
            logger.warning("[NOVELTY] Large reference file detected (%.1f MB compressed). "
                          "Decompression may take several minutes depending on disk I/O speed. "
                          "Please be patient...", compressed_size / (1024*1024))
    except Exception:
        logger.info("[NOVELTY] Decompressing %s → %s", ref_fasta, ref_unc)

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
        logger.info("[NOVELTY] Decompression complete: %.1f MB written in %.1f seconds (%.1f MB/s)",
                    final_mb, elapsed, final_mb / elapsed if elapsed > 0 else 0)
    except Exception as e:
        logger.error("[NOVELTY] Decompression failed: %s", e)
        # Clean up partial file
        try:
            if os.path.exists(ref_unc):
                os.remove(ref_unc)
        except Exception:
            pass
        raise

    return ref_unc


def build_novelty_itol(
    novelty_tsv: str,
    outdir: str,
) -> str:
    """
    Generate an iTOL DATASET_COLORSTRIP file showing novelty as a colour
    gradient from red (<90% identity) through orange/yellow to green (100%).

    This replaces the binary True/False novelty flag with a proper gradient
    that matches the Hungate-style figure.
    """
    output = os.path.join(outdir, "itol_novelty.itol")

    # Identity bins → colour (matches the legend in your existing figures)
    def _identity_color(pct: float) -> str:
        if pct == 0.0:
            return "#8B0000"   # dark red — no hit
        elif pct < 90:
            return "#FF0000"   # red
        elif pct < 91:
            return "#FF2200"
        elif pct < 92:
            return "#FF5500"
        elif pct < 93:
            return "#FF7700"
        elif pct < 94:
            return "#FF9900"
        elif pct < 95:
            return "#FFBB00"
        elif pct < 96:
            return "#DDDD00"
        elif pct < 97:
            return "#AADD00"
        elif pct < 98:
            return "#77DD00"
        elif pct < 99:
            return "#44CC00"
        elif pct < 100:
            return "#22BB00"
        else:
            return "#00AA00"   # dark green — 100% match

    legend_bins = [
        ("<90%",  "#FF0000"),
        ("90%",   "#FF2200"),
        ("91%",   "#FF5500"),
        ("92%",   "#FF7700"),
        ("93%",   "#FF9900"),
        ("94%",   "#FFBB00"),
        ("95%",   "#DDDD00"),
        ("96%",   "#AADD00"),
        ("97%",   "#77DD00"),
        ("98%",   "#44CC00"),
        ("99%",   "#22BB00"),
        ("100%",  "#00AA00"),
    ]

    id_color_pairs = []
    try:
        with open(novelty_tsv) as fh:
            next(fh)  # skip header
            for line in fh:
                parts = line.strip().split("\t")
                if len(parts) < 2:
                    continue
                qid = parts[0]
                try:
                    pct = float(parts[1])
                except ValueError:
                    pct = 0.0
                id_color_pairs.append((qid, _identity_color(pct)))
    except FileNotFoundError:
        logger.warning("[NOVELTY] novelty.tsv not found: %s", novelty_tsv)
        return output

    with open(output, "w") as fh:
        fh.write("DATASET_COLORSTRIP\n")
        fh.write("SEPARATOR COMMA\n")
        fh.write("DATASET_LABEL,Novelty (identity)\n")
        fh.write("COLOR,#FF0000\n")
        fh.write("MARGIN,5\n")
        fh.write("SHOW_INTERNAL,0\n")
        fh.write("LEGEND_TITLE,Novelty (nearest identity)\n")
        fh.write("LEGEND_SHAPES," + ",".join(["1"] * len(legend_bins)) + "\n")
        fh.write("LEGEND_COLORS," + ",".join(c for _, c in legend_bins) + "\n")
        fh.write("LEGEND_LABELS," + ",".join(l for l, _ in legend_bins) + "\n")
        fh.write("DATA\n")
        for qid, color in id_color_pairs:
            fh.write(f"{qid},{color}\n")

    logger.info("[NOVELTY] iTOL novelty strip → %s", output)
    return output


def compute_db_nearest_identities(mapped_derep: str, outdir: str, db, run_dataset: str, threads: Optional[int] = None):
    """Compute nearest-neighbour identities of run sequences against all non-run DB sequences.

    The comparison pool includes ALL datasets stored in the DB except *run_dataset*,
    which means every preload dataset (no matter how many) is automatically included.
    """
    db_preset_fasta = Path(outdir) / 'db_preload_seqs.fasta'
    try:
        with db.connect() as conn:
            cur = conn.cursor()
            cur.execute("SELECT id, sequence FROM sequences WHERE dataset IS NOT ? OR dataset IS NULL", (run_dataset,))
            rows = cur.fetchall()
        records = [(r[0], r[1]) for r in rows]
        if records:
            write_fasta(records, str(db_preset_fasta))
            # Log which datasets are in the pool
            try:
                with db.connect() as conn2:
                    cur2 = conn2.cursor()
                    cur2.execute(
                        "SELECT DISTINCT dataset FROM sequences WHERE dataset IS NOT ? OR dataset IS NULL",
                        (run_dataset,),
                    )
                    pool_ds = [r[0] for r in cur2.fetchall()]
                logger.info("[NOVELTY] Nearest-identity pool: %d seqs from datasets: %s", len(records), pool_ds)
            except Exception:
                pass
        else:
            return {}
    except Exception:
        return {}

    novel_results = {}
    if db_preset_fasta.exists():
        try:
            from branchmanager.utils.subprocess import run_cmd
            novelty_matches = Path(outdir) / 'novelty_matches.tsv'
            thread_flag = f" --threads {int(threads)}" if threads and int(threads) > 0 else ""
            cmd = f"vsearch --usearch_global {mapped_derep} --db {db_preset_fasta} --id 0.5 --strand both --blast6out {novelty_matches} --maxaccepts 1 --maxhits 1{thread_flag}"
            logger.info("[NOVELTY] Running vsearch against preload DB to compute novelty")
            run_cmd(cmd)
            if novelty_matches.exists():
                with open(novelty_matches) as nm:
                    for line in nm:
                        parts = line.strip().split('\t')
                        if len(parts) < 3:
                            continue
                        try:
                            novel_results[parts[0]] = float(parts[2])
                        except Exception:
                            novel_results[parts[0]] = None
        except Exception as e:
            logger.warning("[NOVELTY] vsearch novelty run failed: %s", e)
    return novel_results


def build_run_novelty_itol(outdir: str, run_ids, mapped_derep: str, db, run_dataset: str, orig_to_short: dict, threads: Optional[int] = None):
    """Build the run novelty iTOL strip from nearest-neighbour identities."""
    if not run_ids:
        return None

    from branchmanager.pipeline import itol

    nov_lines = [
        'DATASET_COLORSTRIP',
        'SEPARATOR COMMA',
        'DATASET_LABEL,Novelty (nearest identity)',
        'COLOR,#AAAAAA',
        'MARGIN,5',
        'SHOW_INTERNAL,0'
    ]
    try:
        percents = list(range(90, 101))
        colors = [itol._novelty_color_for_pct('<90')]
        colors += [itol._identity_to_color(p / 100.0, vmin=0.90, vmax=1.0) for p in percents]
        labels = ['<90'] + [f"{p}%" for p in percents]
        nov_lines.append('LEGEND_TITLE,Novelty (identity)')
        nov_lines.append('LEGEND_SHAPES,' + ','.join(['1'] * len(labels)))
        nov_lines.append('LEGEND_COLORS,' + ','.join(colors))
        nov_lines.append('LEGEND_LABELS,' + ','.join(labels))
    except Exception:
        pass

    novel_results = compute_db_nearest_identities(mapped_derep, outdir, db, run_dataset, threads=threads)

    def _candidate_forms_for_short(x):
        forms = []
        try:
            s = orig_to_short.get(x)
            if s:
                forms.append(s)
        except Exception:
            pass
        forms.append(x)
        if '|' in x:
            forms.append(x.split('|')[-1])
        forms.append(x.split()[-1] if x.split() else x)
        try:
            forms.append(forms[0].upper())
            forms.append(forms[0].lower())
        except Exception:
            pass
        import re as _re
        forms.append(_re.sub(r'_n?\d+$', '', forms[0] if forms else x))
        seen = set()
        out = []
        for f in forms:
            if not f:
                continue
            if f not in seen:
                seen.add(f)
                out.append(f)
        return out

    nov_lines.append('DATA')
    for iid in run_ids:
        short_id = None
        if novel_results:
            for cand in _candidate_forms_for_short(iid):
                if cand in novel_results:
                    short_id = cand
                    break
        if short_id is None:
            short_id = iid
        color = '#cccccc'
        if novel_results and short_id in novel_results:
            v = novel_results.get(short_id)
            if v is not None:
                try:
                    color = itol._identity_to_color(float(v), vmin=0.90, vmax=1.0)
                except Exception:
                    color = '#cccccc'
            else:
                try:
                    with db.connect() as conn:
                        cur = conn.cursor()
                        cur.execute('SELECT identity FROM distances WHERE id = ? LIMIT 1', (iid,))
                        r = cur.fetchone()
                        if r:
                            color = itol._identity_to_color(r[0], vmin=0.90, vmax=1.0)
                except Exception:
                    color = '#cccccc'
        nov_lines.append(f'{iid},{color}')

    nov_path = Path(outdir) / 'itol_novelty.itol'
    nov_path.write_text('\n'.join(nov_lines) + '\n')
    logger.info('[ITOL] Wrote novelty ITOL colorstrip to %s', nov_path)
    return str(nov_path)


def _load_novelty_summary(novelty_tsv: str):
    summary = {}
    with open(novelty_tsv) as fh:
        next(fh, None)
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4:
                continue
            try:
                nearest_identity = float(parts[1]) if parts[1] not in ('', 'NA') else 0.0
            except Exception:
                nearest_identity = 0.0
            summary[parts[0]] = {
                'nearest_identity': nearest_identity,
                'nearest_hit': parts[2] if len(parts) > 2 else 'None',
                'novel': parts[3] if len(parts) > 3 else 'True',
            }
    return summary


def _novelty_priority_from_metrics(nearest_identity: float, matches_ge_99: int, matches_ge_97: int):
    if nearest_identity < 97.0 and matches_ge_99 <= 1 and matches_ge_97 <= 3:
        return 'HIGH'
    if nearest_identity < 98.5 and matches_ge_99 <= 2 and matches_ge_97 <= 10:
        return 'MEDIUM'
    return 'LOW'


def _novelty_score_from_metrics(nearest_identity: float, matches_ge_99: int, matches_ge_97: int, matches_ge_95: int) -> float:
    distance_component = max(0.0, min(60.0, (100.0 - float(nearest_identity)) * 3.0))
    density_component = 0.0
    if matches_ge_99 <= 1:
        density_component += 20.0
    elif matches_ge_99 <= 5:
        density_component += 10.0
    if matches_ge_97 <= 3:
        density_component += 15.0
    elif matches_ge_97 <= 10:
        density_component += 7.5
    if matches_ge_95 <= 10:
        density_component += 5.0
    return round(min(100.0, distance_component + density_component), 2)


def _normalise_dataset_names(dataset_names) -> list[str]:
    if not dataset_names:
        return []
    if isinstance(dataset_names, str):
        raw = [dataset_names]
    else:
        raw = list(dataset_names)
    names: list[str] = []
    for item in raw:
        for part in str(item).split(','):
            name = part.strip()
            if name and name not in names:
                names.append(name)
    return names


def _query_db_pool(db, *, run_dataset: Optional[str] = None, include_datasets=None):
    include = _normalise_dataset_names(include_datasets)
    clauses = ["sequence IS NOT NULL", "sequence != ''"]
    params: list[object] = []
    if run_dataset:
        clauses.append("(dataset != ? OR dataset IS NULL)")
        params.append(run_dataset)
    if include:
        placeholders = ','.join('?' for _ in include)
        clauses.append(f"dataset IN ({placeholders})")
        params.extend(include)
    where = ' AND '.join(clauses)
    with db.connect() as conn:
        cur = conn.cursor()
        cur.execute(f"SELECT id, sequence, dataset FROM sequences WHERE {where}", tuple(params))
        rows = cur.fetchall()
    records = [(str(r[0]), str(r[1])) for r in rows if r[0] and r[1]]
    datasets = sorted({str(r[2]) for r in rows if r[2] not in (None, '')})
    return records, datasets


def _write_pool_fasta(records, path: Path) -> Optional[str]:
    if not records:
        return None
    path.parent.mkdir(parents=True, exist_ok=True)
    write_fasta(records, str(path))
    return str(path)


def _build_db_pool_fasta(
    *,
    db,
    outdir: str,
    name: str,
    run_dataset: Optional[str] = None,
    include_datasets=None,
) -> tuple[Optional[str], str, list[str], int]:
    if db is None:
        return None, 'none', [], 0
    try:
        records, datasets = _query_db_pool(
            db,
            run_dataset=run_dataset,
            include_datasets=include_datasets,
        )
    except Exception as e:
        logger.warning('[NOVELTY] Failed to query %s comparison pool from DB: %s', name, e)
        return None, 'none', [], 0
    path = _write_pool_fasta(records, Path(outdir) / f'novelty_{name}_pool.fasta')
    if not path:
        return None, 'none', datasets, 0
    return path, name, datasets, len(records)


def _run_nearest_search(input_fasta: str, db_fasta: str, matches_path: Path, threads: Optional[int] = None):
    thread_flag = f" --threads {int(threads)}" if threads and int(threads) > 0 else ""
    cmd = (
        f"vsearch --usearch_global {input_fasta} "
        f"--db {db_fasta} "
        f"--id 0.5 "
        f"--strand both "
        f"--blast6out {matches_path} "
        f"--maxaccepts 1 --maxhits 1 "
        f"--query_cov 0.7{thread_flag}"
    )
    run_cmd(cmd)
    return _parse_best_hits(str(matches_path))


def _run_density_search(input_fasta: str, db_fasta: str, matches_path: Path, threads: Optional[int] = None):
    thread_flag = f" --threads {int(threads)}" if threads and int(threads) > 0 else ""
    cmd = (
        f"vsearch --usearch_global {input_fasta} "
        f"--db {db_fasta} "
        f"--id 0.95 "
        f"--strand both "
        f"--blast6out {matches_path} "
        f"--maxaccepts 0 --maxhits 0 "
        f"--query_cov 0.7{thread_flag}"
    )
    run_cmd(cmd)

    density_counts = {}
    if matches_path.exists():
        with open(matches_path) as fh:
            for line in fh:
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 3:
                    continue
                qid = parts[0]
                try:
                    pct = float(parts[2])
                except Exception:
                    continue
                bucket = density_counts.setdefault(qid, {'ge_99': 0, 'ge_97': 0, 'ge_95': 0})
                if pct >= 95.0:
                    bucket['ge_95'] += 1
                if pct >= 97.0:
                    bucket['ge_97'] += 1
                if pct >= 99.0:
                    bucket['ge_99'] += 1
    return density_counts


def _crowding_from_counts(counts: dict) -> str:
    if counts['ge_99'] <= 1 and counts['ge_97'] <= 1:
        return 'isolated'
    if counts['ge_97'] <= 3:
        return 'sparse'
    if counts['ge_97'] <= 10:
        return 'moderate'
    return 'crowded'


def _metrics_from_nearest_and_density(qids, nearest_hits, density_counts, *, source: str, id_threshold: float):
    novel_threshold_pct = id_threshold * 100.0
    rows = {}
    for qid in qids:
        pct_id, hit_id = nearest_hits.get(qid, (0.0, 'None'))
        counts = density_counts.get(qid, {'ge_99': 0, 'ge_97': 0, 'ge_95': 0})
        crowding = _crowding_from_counts(counts)
        score = _novelty_score_from_metrics(pct_id, counts['ge_99'], counts['ge_97'], counts['ge_95'])
        priority = _novelty_priority_from_metrics(pct_id, counts['ge_99'], counts['ge_97'])
        rows[qid] = {
            'nearest_identity': pct_id,
            'nearest_hit': hit_id,
            'novel': str(pct_id < novel_threshold_pct),
            'matches_ge_99': counts['ge_99'],
            'matches_ge_97': counts['ge_97'],
            'matches_ge_95': counts['ge_95'],
            'novelty_score': score,
            'crowding': crowding,
            'sequencing_priority': priority,
            'density_source': source,
        }
    return rows


def _empty_pool_metrics(qids, *, source: str = 'none'):
    return {
        qid: {
            'nearest_identity': 0.0,
            'nearest_hit': 'None',
            'novel': 'True',
            'matches_ge_99': 0,
            'matches_ge_97': 0,
            'matches_ge_95': 0,
            'novelty_score': 0.0,
            'crowding': 'sparse',
            'sequencing_priority': 'HIGH',
            'density_source': source,
        }
        for qid in qids
    }


def _empty_sequencing_context(qids):
    return {
        qid: {
            'partner_id': 'NA',
            'selected_for_wgs': 'NA',
            'nearest_selected_hit': 'None',
            'nearest_selected_identity': 0.0,
            'selected_ge_99': 0,
            'selected_ge_97': 0,
            'selected_ge_95': 0,
            'clade_already_selected': 'False',
            'adjusted_novelty_score': 'NA',
            'adjusted_priority': 'NA',
            'source': 'none',
        }
        for qid in qids
    }


def _compute_sequencing_context(input_fasta: str, outdir: str, db, qids, primary_metrics, threads: Optional[int] = None):
    if db is None:
        return _empty_sequencing_context(qids)

    try:
        metadata = db.get_sequencing_metadata_for_ids(qids)
    except Exception:
        metadata = {}

    selected_records = []
    try:
        with db.connect() as conn:
            cur = conn.cursor()
            cur.execute(
                "SELECT s.id, s.sequence FROM sequences s "
                "JOIN sequencing_metadata m ON s.id = m.id "
                "WHERE m.selected_for_wgs = 1 AND s.sequence IS NOT NULL AND s.sequence != ''"
            )
            selected_records = [(str(r[0]), str(r[1])) for r in cur.fetchall() if r[0] and r[1]]
    except Exception as e:
        logger.warning('[NOVELTY] Failed to load selected-for-WGS pool from DB: %s', e)
        selected_records = []

    selected_hits = {qid: [] for qid in qids}
    if selected_records:
        selected_pool = Path(outdir) / 'novelty_selected_for_wgs_pool.fasta'
        write_fasta(selected_records, str(selected_pool))
        matches_path = Path(outdir) / 'novelty_selected_for_wgs_matches.tsv'
        try:
            _run_density_search(input_fasta, str(selected_pool), matches_path, threads=threads)
            if matches_path.exists():
                with open(matches_path) as fh:
                    for line in fh:
                        parts = line.rstrip('\n').split('\t')
                        if len(parts) < 3:
                            continue
                        qid, hit_id = parts[0], parts[1]
                        if qid not in selected_hits or hit_id == qid:
                            continue
                        try:
                            pct = float(parts[2])
                        except Exception:
                            continue
                        selected_hits[qid].append((pct, hit_id))
        except Exception as e:
            logger.warning('[NOVELTY] Failed to compute selected-for-WGS neighbourhood: %s', e)

    context = {}
    for qid in qids:
        hits = selected_hits.get(qid, [])
        ge_99 = sum(1 for pct, _ in hits if pct >= 99.0)
        ge_97 = sum(1 for pct, _ in hits if pct >= 97.0)
        ge_95 = sum(1 for pct, _ in hits if pct >= 95.0)
        nearest_pct, nearest_hit = (0.0, 'None')
        if hits:
            nearest_pct, nearest_hit = max(hits, key=lambda item: item[0])
        meta = metadata.get(qid, {})
        current_selected = meta.get('selected_for_wgs')
        base = primary_metrics.get(qid, {})
        try:
            base_score = float(base.get('novelty_score', 0.0))
        except Exception:
            base_score = 0.0
        base_priority = str(base.get('sequencing_priority', 'NA'))
        clade_already_selected = ge_97 > 0
        if current_selected:
            adjusted_priority = 'SELECTED'
            adjusted_score = base_score
        elif clade_already_selected:
            adjusted_priority = 'LOW_ALREADY_SELECTED_CLADE'
            adjusted_score = max(0.0, base_score - 25.0)
        elif ge_95 > 0 and base_priority == 'HIGH':
            adjusted_priority = 'MEDIUM_SELECTED_NEARBY'
            adjusted_score = max(0.0, base_score - 10.0)
        else:
            adjusted_priority = base_priority
            adjusted_score = base_score
        context[qid] = {
            'partner_id': meta.get('partner_id', 'NA'),
            'selected_for_wgs': 'True' if current_selected is True else ('False' if current_selected is False else 'NA'),
            'nearest_selected_hit': nearest_hit,
            'nearest_selected_identity': nearest_pct,
            'selected_ge_99': ge_99,
            'selected_ge_97': ge_97,
            'selected_ge_95': ge_95,
            'clade_already_selected': str(clade_already_selected),
            'adjusted_novelty_score': adjusted_score,
            'adjusted_priority': adjusted_priority,
            'source': 'sequencing_metadata' if metadata or selected_records else 'none',
        }
    return context


def _write_novelty_metrics_table(output: Path, qids, primary_metrics, baseline_metrics, all_known_metrics, reference_metrics, sequencing_context=None):
    headers = [
        'ID',
        'NearestIdentity', 'NearestHit', 'Novel', 'MatchesGE99', 'MatchesGE97', 'MatchesGE95',
        'NoveltyScore', 'Crowding', 'SequencingPriority', 'DensitySource',
        'BaselineNearestIdentity', 'BaselineNearestHit', 'BaselineNovel',
        'BaselineMatchesGE99', 'BaselineMatchesGE97', 'BaselineMatchesGE95',
        'BaselineNoveltyScore', 'BaselineCrowding', 'BaselineSequencingPriority', 'BaselineDensitySource',
        'AllKnownNearestIdentity', 'AllKnownNearestHit', 'AllKnownNovel',
        'AllKnownMatchesGE99', 'AllKnownMatchesGE97', 'AllKnownMatchesGE95',
        'AllKnownNoveltyScore', 'AllKnownCrowding', 'AllKnownSequencingPriority', 'AllKnownDensitySource',
        'ReferenceNearestIdentity', 'ReferenceNearestHit', 'ReferenceNovel',
        'ReferenceMatchesGE99', 'ReferenceMatchesGE97', 'ReferenceMatchesGE95',
        'ReferenceNoveltyScore', 'ReferenceCrowding', 'ReferenceSequencingPriority', 'ReferenceDensitySource',
        'PartnerID', 'SelectedForGenomeSequencing',
        'NearestSelectedGenomeHit', 'NearestSelectedGenomeIdentity',
        'SelectedGenomeMatchesGE99', 'SelectedGenomeMatchesGE97', 'SelectedGenomeMatchesGE95',
        'CladeAlreadySelectedForGenomeSequencing',
        'GenomeSequencingAdjustedNoveltyScore', 'GenomeSequencingAdjustedPriority',
        'GenomeSequencingMetadataSource',
    ]

    def vals(metrics, qid):
        row = metrics.get(qid) or _empty_pool_metrics([qid])[qid]
        return [
            f"{float(row['nearest_identity']):.2f}",
            row['nearest_hit'],
            row['novel'],
            str(row['matches_ge_99']),
            str(row['matches_ge_97']),
            str(row['matches_ge_95']),
            f"{float(row['novelty_score']):.2f}",
            row['crowding'],
            row['sequencing_priority'],
            row['density_source'],
        ]

    with open(output, 'w') as fh:
        fh.write('\t'.join(headers) + '\n')
        sequencing_context = sequencing_context or _empty_sequencing_context(qids)
        for qid in qids:
            primary = vals(primary_metrics, qid)
            baseline = vals(baseline_metrics, qid)
            all_known = vals(all_known_metrics, qid)
            reference = vals(reference_metrics, qid)
            seq_ctx = sequencing_context.get(qid) or _empty_sequencing_context([qid])[qid]
            sequencing = [
                seq_ctx['partner_id'],
                seq_ctx['selected_for_wgs'],
                seq_ctx['nearest_selected_hit'],
                f"{float(seq_ctx['nearest_selected_identity']):.2f}",
                str(seq_ctx['selected_ge_99']),
                str(seq_ctx['selected_ge_97']),
                str(seq_ctx['selected_ge_95']),
                seq_ctx['clade_already_selected'],
                f"{float(seq_ctx['adjusted_novelty_score']):.2f}" if seq_ctx['adjusted_novelty_score'] != 'NA' else 'NA',
                seq_ctx['adjusted_priority'],
                seq_ctx['source'],
            ]
            fh.write('\t'.join([qid] + primary + baseline + all_known + reference + sequencing) + '\n')
    return str(output)


def build_reference_novelty_metrics(
    input_fasta: str,
    ref_fasta: str,
    novelty_tsv: str,
    outdir: str,
    threads: Optional[int] = None,
    db=None,
    run_dataset: Optional[str] = None,
    target_fasta: Optional[str] = None,
    baseline_datasets=None,
):
    """Write per-sequence novelty metrics comparing against user-submitted sequences.

    Density (neighbourhood crowding) is computed against previously submitted
    sequences (preload + all other run datasets stored in the DB), NOT against
    the external reference database.  This ensures that novelty is always
    expressed relative to the sequences *the user* has already characterised,
    which is the core BranchManager design goal.

    If ``target_fasta`` is provided it overrides both the DB lookup and the
    ``ref_fasta`` fallback — only those sequences are used as the density
    comparison pool.

    Falls back to ``ref_fasta`` only when no DB sequences and no
    ``target_fasta`` are available.

    Parameters
    ----------
    input_fasta  : Query sequences (dereplicated, QC-passed).
    ref_fasta    : External reference FASTA — used as a last resort fallback only.
    novelty_tsv  : Path to novelty.tsv produced by run_novelty.
    outdir       : Output directory.
    threads      : VSEARCH thread count.
    db           : Database instance (used to pull previously submitted seqs).
    run_dataset  : Current run dataset name (excluded from density DB).
    target_fasta : Explicit FASTA of sequences to compare against; overrides DB/ref.
    baseline_datasets : DB dataset label(s) to treat as the cultured/reference
                        baseline pool, e.g. Hungate. Can be a string,
                        comma-separated string, or iterable.

    If the DB contains ``sequencing_metadata`` rows, the output also reports
    whether each query or any close 16S neighbour (>=97% identity) has already
    been selected for WGS/full-genome sequencing.
    """
    summary = _load_novelty_summary(novelty_tsv)
    output = Path(outdir) / 'novelty_metrics.tsv'
    qids = list(summary.keys())
    if not qids:
        qids = [
            h for h, _ in read_fasta(input_fasta)
            if not tree_pipeline.is_ref_anchor(h)
        ]

    all_pool_path: Optional[str] = None
    all_source = 'none'
    all_datasets: list[str] = []
    all_count = 0
    if target_fasta and Path(target_fasta).exists():
        all_pool_path = target_fasta
        all_source = 'target_fasta'
        logger.info('[NOVELTY] Using target FASTA for all-known novelty metrics: %s', target_fasta)
    elif db is not None:
        all_pool_path, all_source, all_datasets, all_count = _build_db_pool_fasta(
            db=db,
            outdir=outdir,
            name='all_known',
            run_dataset=run_dataset,
        )
        if all_pool_path:
            logger.info(
                '[NOVELTY] All-known pool: %d sequences from %d dataset(s): %s',
                all_count, len(all_datasets), all_datasets,
            )

    if all_pool_path is None and ref_fasta:
        all_pool_path = _ensure_uncompressed(ref_fasta, outdir)
        all_source = 'reference_fasta'
        logger.info(
            '[NOVELTY] WARNING: No user-submitted all-known pool available. '
            'Falling back to external reference (%s).',
            ref_fasta,
        )

    if all_pool_path:
        all_nearest = {
            qid: (meta['nearest_identity'], meta['nearest_hit'])
            for qid, meta in summary.items()
        }
        all_density = _run_density_search(
            input_fasta,
            all_pool_path,
            Path(outdir) / 'novelty_all_known_density_matches.tsv',
            threads=threads,
        )
        all_metrics = _metrics_from_nearest_and_density(
            qids,
            all_nearest,
            all_density,
            source=all_source,
            id_threshold=0.97,
        )
    else:
        logger.warning('[NOVELTY] No all-known comparison pool available; all-known metrics will be empty')
        all_metrics = _empty_pool_metrics(qids, source='none')

    reference_pool_path: Optional[str] = None
    reference_source = 'none'
    if ref_fasta:
        try:
            reference_pool_path = _ensure_uncompressed(ref_fasta, outdir)
            reference_source = 'reference_fasta'
        except Exception as e:
            logger.warning('[NOVELTY] Could not prepare reference novelty pool %s: %s', ref_fasta, e)
            reference_pool_path = None

    if reference_pool_path:
        reference_nearest = _run_nearest_search(
            input_fasta,
            reference_pool_path,
            Path(outdir) / 'novelty_reference_nearest_matches.tsv',
            threads=threads,
        )
        reference_density = _run_density_search(
            input_fasta,
            reference_pool_path,
            Path(outdir) / 'novelty_reference_density_matches.tsv',
            threads=threads,
        )
        reference_metrics = _metrics_from_nearest_and_density(
            qids,
            reference_nearest,
            reference_density,
            source=reference_source,
            id_threshold=0.97,
        )
    else:
        reference_metrics = _empty_pool_metrics(qids, source='none')

    baseline_names = _normalise_dataset_names(baseline_datasets)
    baseline_pool_path: Optional[str] = None
    baseline_source = 'none'
    baseline_datasets_found: list[str] = []
    baseline_count = 0
    if baseline_names and db is not None:
        baseline_pool_path, baseline_source, baseline_datasets_found, baseline_count = _build_db_pool_fasta(
            db=db,
            outdir=outdir,
            name='baseline',
            run_dataset=run_dataset,
            include_datasets=baseline_names,
        )
        if baseline_pool_path:
            baseline_source = 'baseline:' + ','.join(baseline_datasets_found or baseline_names)
            logger.info(
                '[NOVELTY] Baseline pool: %d sequences from dataset(s): %s',
                baseline_count, baseline_datasets_found or baseline_names,
            )
        else:
            logger.info('[NOVELTY] No baseline sequences found for dataset(s): %s', baseline_names)

    if baseline_pool_path:
        baseline_nearest = _run_nearest_search(
            input_fasta,
            baseline_pool_path,
            Path(outdir) / 'novelty_baseline_nearest_matches.tsv',
            threads=threads,
        )
        baseline_density = _run_density_search(
            input_fasta,
            baseline_pool_path,
            Path(outdir) / 'novelty_baseline_density_matches.tsv',
            threads=threads,
        )
        baseline_metrics = _metrics_from_nearest_and_density(
            qids,
            baseline_nearest,
            baseline_density,
            source=baseline_source,
            id_threshold=0.97,
        )
    else:
        baseline_metrics = _empty_pool_metrics(qids, source='none')

    primary_metrics = baseline_metrics if baseline_pool_path else all_metrics
    sequencing_context = _compute_sequencing_context(
        input_fasta,
        outdir,
        db,
        qids,
        primary_metrics,
        threads=threads,
    )
    _write_novelty_metrics_table(
        output,
        qids,
        primary_metrics,
        baseline_metrics,
        all_metrics,
        reference_metrics,
        sequencing_context,
    )
    logger.info(
        '[NOVELTY] Wrote novelty metrics to %s (primary=%s; baseline=%s; all-known=%s; reference=%s)',
        output,
        'baseline' if baseline_pool_path else all_source,
        baseline_source,
        all_source,
        reference_source,
    )
    return str(output)
