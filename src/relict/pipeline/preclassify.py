"""Pre-classification pipeline for Relict.

This module classifies one or more FASTA files (reference collections such as
Hungate 16S, SILVA, RDP etc.) against a reference database using vsearch and
writes taxonomy outputs that the main Relict pipeline can consume without
requiring on-the-fly classification at each run.

Output files written to *outdir*
---------------------------------
  {dataset}_classification.tsv   Full vsearch matches (ID, BestHit, Identity,
                                  Taxon, Confidence).  Human readable.
  {dataset}_taxonomy.tsv         Condensed (ID, Taxon, Confidence).  Passable
                                  directly to ``relict preload --taxa-assignments``.
  combined_taxonomy.tsv          All datasets merged with Dataset column.
  pipeline_taxonomy.tsv          All datasets merged without Dataset column;
                                  passable directly to ``relict preload`` or
                                  ``relict run`` via ``--taxa-assignments``.
  preclassify_summary.txt        Plain-text human-readable classification
                                  summary (counts, dataset labels, top hits).
"""

from __future__ import annotations

import gzip
import logging
import os
import re
import textwrap
from pathlib import Path
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Known 16S dataset registry
# ---------------------------------------------------------------------------
KNOWN_DATASETS: Dict[str, dict] = {
    "hungate16s": {
        "label": "Hungate16S",
        "description": "Hungate1000 project – curated full-length 16S rRNA "
                       "gene sequences from rumen bacteria and archaea",
        "url": "https://www.nature.com/articles/nbt.4387",
        "typical_min_len": 1200,
        "color": "#e6ab02",
    },
    "hungate": {
        "label": "Hungate16S",
        "description": "Alias for hungate16s",
        "url": "https://www.nature.com/articles/nbt.4387",
        "typical_min_len": 1200,
        "color": "#e6ab02",
    },
    "silva": {
        "label": "SILVA",
        "description": "SILVA ribosomal RNA gene database project",
        "url": "https://www.arb-silva.de",
        "typical_min_len": 1000,
        "color": "#1b9e77",
    },
    "rdp": {
        "label": "RDP",
        "description": "Ribosomal Database Project",
        "url": "https://rdp.cme.msu.edu",
        "typical_min_len": 1000,
        "color": "#7570b3",
    },
    "homd": {
        "label": "HOMD",
        "description": "Human Oral Microbiome Database",
        "url": "https://www.homd.org",
        "typical_min_len": 900,
        "color": "#d95f02",
    },
    "greengenes2": {
        "label": "GreenGenes2",
        "description": "GreenGenes2 16S rRNA reference database",
        "url": "https://greengenes2.ucsd.edu",
        "typical_min_len": 1000,
        "color": "#66a61e",
    },
    "gg2": {
        "label": "GreenGenes2",
        "description": "Alias for greengenes2",
        "url": "https://greengenes2.ucsd.edu",
        "typical_min_len": 1000,
        "color": "#66a61e",
    },
    "ncbi16s": {
        "label": "NCBI16S",
        "description": "NCBI 16S rRNA RefSeq collection",
        "url": "https://www.ncbi.nlm.nih.gov/refseq/targetedloci/",
        "typical_min_len": 900,
        "color": "#a6761d",
    },
    "gtdb": {
        "label": "GTDB",
        "description": "Genome Taxonomy Database 16S rRNA representative sequences",
        "url": "https://gtdb.ecogenomic.org",
        "typical_min_len": 1000,
        "color": "#386cb0",
    },
}


def resolve_dataset_meta(name: str) -> dict:
    """Return metadata for *name* (case-insensitive lookup against registry).

    Falls back to a generic entry so callers never have to handle ``None``.
    """
    key = name.lower().replace("-", "").replace("_", "").replace(" ", "")
    meta = KNOWN_DATASETS.get(key)
    if meta is None:
        for k, v in KNOWN_DATASETS.items():
            if k in key or key in k:
                meta = v
                break
    if meta is None:
        meta = {
            "label": name,
            "description": f"User-supplied dataset '{name}'",
            "url": "",
            "typical_min_len": 900,
            "color": _name_to_color(name),
        }
    return meta


def _name_to_color(name: str) -> str:
    """Deterministic hex color derived from *name*."""
    import hashlib
    h = hashlib.md5(name.encode()).hexdigest()
    r = int(h[0:2], 16) | 0x40
    g = int(h[2:4], 16) | 0x40
    b = int(h[4:6], 16) | 0x40
    return f"#{r:02x}{g:02x}{b:02x}"


# ---------------------------------------------------------------------------
# ID normalisation helpers  (mirrors logic in pipeline/classify.py)
# ---------------------------------------------------------------------------

def _norm_id(x: str) -> str:
    x = x.split()[0]
    if "|" in x:
        x = x.split("|")[-1]
    if "#" in x:
        x = x.split("#")[0]
    if ":" in x:
        x = x.split(":")[0]
    x = x.strip().strip("()[]{}")
    x = re.sub(r"[^0-9A-Za-z]+", "_", x).strip("_")
    return x


def _canon_id(x: Optional[str]) -> Optional[str]:
    if x is None:
        return x
    y = str(x).split()[0]
    if "#" in y:
        y = y.split("#")[0]
    y = re.sub(r":\d+-\d+\(.*\)$", "", y)
    y = re.sub(r":\d+-\d+$", "", y)
    if "|" in y:
        y = y.split("|")[-1]
    y = re.sub(r"[^0-9A-Za-z]+", "_", y).strip("_")
    return y


# ---------------------------------------------------------------------------
# Core classification logic
# ---------------------------------------------------------------------------

def _load_taxa_map(taxa_tsv: Optional[str]) -> dict:
    """Load a FeatureID→(Taxon, Confidence) mapping from *taxa_tsv*.

    Accepts gzipped files.  Stores three keys per record (original, normalised,
    canonical) to maximise hit rate when resolving vsearch best-hit IDs.
    """
    taxa_map: dict = {}
    if not taxa_tsv:
        return taxa_map
    open_fn = gzip.open if str(taxa_tsv).endswith(".gz") else open
    try:
        with open_fn(taxa_tsv, "rt") as fh:  # type: ignore[call-overload]
            first = fh.readline()
            if not any(kw in first.lower() for kw in ("feature", "taxon", "taxonomy", "id\t")):
                parts = first.strip().split("\t")
                if len(parts) >= 2:
                    fid, tax = parts[0], parts[1]
                    try:
                        conf: Optional[float] = float(parts[2]) if len(parts) > 2 else None
                    except Exception:
                        conf = None
                    for k in (fid, _norm_id(fid), _canon_id(fid)):
                        if k:
                            taxa_map[k] = (tax, conf)
            for line in fh:
                parts = line.strip().split("\t")
                if len(parts) < 2:
                    continue
                fid, tax = parts[0], parts[1]
                try:
                    conf = float(parts[2]) if len(parts) > 2 else None
                except Exception:
                    conf = None
                for k in (fid, _norm_id(fid), _canon_id(fid)):
                    if k:
                        taxa_map[k] = (tax, conf)
        logger.info("[PRECLASSIFY] Loaded %d taxa mappings from %s", len(taxa_map), taxa_tsv)
    except Exception as exc:
        logger.warning("[PRECLASSIFY] Failed to load taxa_tsv %s: %s", taxa_tsv, exc)
    return taxa_map


def _resolve_taxon(sid: str, taxa_map: dict) -> Tuple[Optional[str], Optional[float]]:
    """Look up the taxonomy for vsearch best-hit *sid* in *taxa_map*."""
    for key in (sid, _norm_id(sid), _canon_id(sid)):
        if key in taxa_map:
            return taxa_map[key]
    tokens = [t for t in re.split(r"[|_:\-.]+", sid) if t]
    for tok in tokens:
        for k in (tok, _norm_id(tok), _canon_id(tok)):
            if k in taxa_map:
                return taxa_map[k]
    return None, None


def _is_fasta_file(path: str) -> bool:
    """Return True if *path* looks like a FASTA file (not a TSV taxonomy file)."""
    p = str(path).lower()
    if p.endswith('.gz'):
        p = p[:-3]
    return p.endswith(('.fasta', '.fa', '.fna', '.ffn', '.faa', '.frn'))


# ---------------------------------------------------------------------------
# Retry cascade configuration
# ---------------------------------------------------------------------------
# Sequences unclassified at round 1 are re-searched with progressively relaxed
# parameters.  The main cause of NCBI BLAST succeeding where vsearch fails is
# the --query_cov filter (NCBI has none by default) and vsearch's --maxrejects
# heuristic giving up too early on N-rich / divergent sequences.
#
# Each round entry: (query_cov, maxrejects_multiplier, label)
#   query_cov           – 0.0 means no --query_cov flag is passed at all
#   maxrejects_mult     – multiplied against the initial maxrejects value
#   label               – used in log messages and TSV annotations
_RETRY_ROUNDS: List[Tuple[float, int, str]] = [
    (0.5, 2,  "retry1_cov50"),   # lower coverage, 2× rejects
    (0.3, 4,  "retry2_cov30"),   # even lower coverage, 4× rejects
    (0.0, 8,  "retry3_nocov"),   # no coverage filter at all, 8× rejects
]


def _run_vsearch_pass(
    query_fasta: str,
    ref_fasta: str,
    matches_path: str,
    min_identity: float,
    max_hits: int,
    max_rejects: int,
    query_cov: float,
    threads: int,
    dataset_name: str,
    round_label: str,
) -> Dict[str, Tuple[str, float]]:
    """Run one vsearch --usearch_global pass and return best hit per query.

    Returns ``{qid: (sid, pct_identity)}`` containing only the highest-identity
    hit for each query that appears in *matches_path*.
    """
    import shlex
    import subprocess

    cmd_list = [
        "vsearch",
        "--usearch_global", query_fasta,
        "--db", ref_fasta,
        "--id", str(min_identity),
        "--strand", "both",
        "--blast6out", matches_path,
        "--maxaccepts", str(int(max_hits)),
        "--maxhits",    str(int(max_hits)),
        "--maxrejects", str(int(max_rejects)),
    ]
    if query_cov > 0.0:
        cmd_list += ["--query_cov", str(query_cov)]
    if threads and int(threads) > 0:
        cmd_list += ["--threads", str(int(threads))]

    logger.info(
        "[PRECLASSIFY] vsearch %s for '%s' (qcov=%s, maxrejects=%d, id≥%.2f)",
        round_label, dataset_name,
        f"{query_cov:.0%}" if query_cov > 0 else "none",
        max_rejects, min_identity,
    )
    logger.debug("[PRECLASSIFY] cmd: %s", " ".join(shlex.quote(a) for a in cmd_list))

    result = subprocess.run(cmd_list, shell=False)
    if result.returncode != 0:
        raise RuntimeError(
            f"vsearch failed ({round_label}) for dataset '{dataset_name}' "
            f"(exit {result.returncode}). "
            f"Command: {' '.join(shlex.quote(a) for a in cmd_list)}"
        )

    best: Dict[str, Tuple[str, float]] = {}
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
                    pct = 0.0
                if qid not in best or pct > best[qid][1]:
                    best[qid] = (sid, pct)
    except FileNotFoundError:
        pass
    return best


def _write_fasta_subset(source_fasta: str, target_ids: set, out_path: str) -> int:
    """Write sequences from *source_fasta* whose IDs are in *target_ids* to *out_path*.

    Does a plain-text scan so it works even when relict.utils.fasta is unavailable.
    Returns the number of sequences written.
    """
    written = 0
    open_fn = gzip.open if str(source_fasta).endswith(".gz") else open
    inside = False
    buf: List[str] = []
    try:
        with open_fn(source_fasta, "rt") as fh:  # type: ignore[call-overload]
            with open(out_path, "w") as out:
                for raw in fh:
                    raw = raw.rstrip("\n\r")
                    if raw.startswith(">"):
                        rid = raw[1:].split()[0]
                        inside = rid in target_ids
                        if inside:
                            out.write(raw + "\n")
                            written += 1
                    elif inside:
                        out.write(raw + "\n")
    except Exception as exc:
        logger.warning("[PRECLASSIFY] _write_fasta_subset failed: %s", exc)
    return written


def classify_fasta(
    fasta_path: str,
    ref_fasta: str,
    outdir: str,
    dataset_name: str,
    taxa_tsv: Optional[str] = None,
    threads: int = 4,
    min_identity: float = 0.80,
    low_confidence_threshold: float = 0.97,
    max_hits: int = 10,
    max_rejects: int = 256,
) -> dict:
    """Classify sequences in *fasta_path* against *ref_fasta* using vsearch.

    vsearch collects up to *max_hits* candidate hits per query sequence
    (``--maxaccepts max_hits --maxhits max_hits --maxrejects max_rejects``).
    After all hits are collected, only the **best hit per query** (highest
    % identity) is kept for taxonomy assignment.  Using multiple hits + a
    post-filter is far more robust than ``--maxhits 1``, which stops at the
    first acceptable hit regardless of whether a better match exists nearby
    in the k-mer index.  *max_rejects* is also raised from vsearch's default
    of 32 so that ambiguous sequences (e.g. those containing many N's) are
    not abandoned prematurely.

    Writes four output files:
    - ``{outdir}/{dataset_name}_classification.tsv``  (full, human-readable)
    - ``{outdir}/{dataset_name}_taxonomy.tsv``        (condensed, pipeline-ready)
    - ``{outdir}/{dataset_name}_unclassified.tsv``    (sequences with no vsearch hit)
    - ``{outdir}/{dataset_name}_low_confidence.tsv``  (hits below *low_confidence_threshold*)

    Returns a dict with keys: classification_tsv, taxonomy_tsv, unclassified_tsv,
    low_confidence_tsv, n_input, n_classified, n_unclassified, n_low_confidence.
    """
    from relict.utils.fasta import read_fasta, write_fasta

    Path(outdir).mkdir(parents=True, exist_ok=True)

    # ── Normalise query FASTA ────────────────────────────────────────────────
    norm_query = os.path.join(outdir, f"{dataset_name}_input_norm.fasta")
    input_ids: Dict[str, int] = {}  # id → sequence length
    try:
        records = list(read_fasta(fasta_path))
        if not records:
            raise ValueError(f"No sequences found in {fasta_path}")
        write_fasta(records, norm_query)
        for rec in records:
            # read_fasta typically yields (header, seq) or objects; handle both
            if isinstance(rec, tuple):
                hdr, seq = rec[0], rec[1]
            else:
                hdr, seq = rec.id, str(rec.seq)
            rid = hdr.split()[0].lstrip(">")
            input_ids[rid] = len(seq)
        logger.info(
            "[PRECLASSIFY] Normalised %d sequences from '%s' → %s",
            len(records), dataset_name, norm_query,
        )
    except Exception as exc:
        logger.warning(
            "[PRECLASSIFY] Could not normalise query FASTA (%s); using original: %s",
            exc, fasta_path,
        )
        norm_query = fasta_path
        # Fall back: collect IDs by scanning original file
        try:
            open_fn = gzip.open if str(fasta_path).endswith(".gz") else open
            with open_fn(fasta_path, "rt") as fh:  # type: ignore[call-overload]
                cur_id: Optional[str] = None
                cur_len = 0
                for raw in fh:
                    raw = raw.rstrip()
                    if raw.startswith(">"):
                        if cur_id is not None:
                            input_ids[cur_id] = cur_len
                        cur_id = raw[1:].split()[0]
                        cur_len = 0
                    else:
                        cur_len += len(raw)
                if cur_id is not None:
                    input_ids[cur_id] = cur_len
        except Exception:
            pass

    n_input = len(input_ids)

    # ── Decompress reference if gzipped ─────────────────────────────────────
    ref_to_use = ref_fasta
    if str(ref_fasta).endswith(".gz"):
        ref_unc = os.path.join(outdir, f"{dataset_name}_ref_uncompressed.fasta")
        if not os.path.exists(ref_unc):
            logger.info("[PRECLASSIFY] Decompressing reference %s → %s", ref_fasta, ref_unc)
            recs = list(read_fasta(ref_fasta))
            write_fasta(recs, ref_unc)
        ref_to_use = ref_unc

    # ── Round 1: initial vsearch pass ───────────────────────────────────────
    matches_r1 = os.path.join(outdir, f"{dataset_name}_matches.tsv")
    best_hit: Dict[str, Tuple[str, float]] = _run_vsearch_pass(
        query_fasta=norm_query,
        ref_fasta=ref_to_use,
        matches_path=matches_r1,
        min_identity=min_identity,
        max_hits=max_hits,
        max_rejects=max_rejects,
        query_cov=0.7,
        threads=threads,
        dataset_name=dataset_name,
        round_label="round1",
    )
    # Track which round classified each sequence (for reporting)
    classify_round: Dict[str, str] = {qid: "round1" for qid in best_hit}

    logger.info(
        "[PRECLASSIFY] Round 1: %d/%d sequences classified",
        len(best_hit), n_input,
    )

    # ── Retry cascade: re-run on unclassified sequences with relaxed params ─
    # Each retry uses a lower --query_cov and higher --maxrejects.
    # This mirrors NCBI BLAST's behaviour (no query-coverage filter by default)
    # and gives N-rich / partial sequences more chances to find a valid hit.
    for retry_qcov, reject_mult, round_label in _RETRY_ROUNDS:
        still_unclassified = {rid for rid in input_ids if rid not in best_hit}
        if not still_unclassified:
            break

        subset_fasta = os.path.join(outdir, f"{dataset_name}_{round_label}_input.fasta")
        n_sub = _write_fasta_subset(norm_query, still_unclassified, subset_fasta)
        if n_sub == 0:
            break

        retry_matches = os.path.join(outdir, f"{dataset_name}_{round_label}_matches.tsv")
        retry_rejects  = max_rejects * reject_mult

        try:
            new_hits = _run_vsearch_pass(
                query_fasta=subset_fasta,
                ref_fasta=ref_to_use,
                matches_path=retry_matches,
                min_identity=min_identity,
                max_hits=max_hits,
                max_rejects=retry_rejects,
                query_cov=retry_qcov,
                threads=threads,
                dataset_name=dataset_name,
                round_label=round_label,
            )
        except Exception as exc:
            logger.warning("[PRECLASSIFY] %s failed: %s — skipping", round_label, exc)
            continue

        newly_found = 0
        for qid, hit in new_hits.items():
            if qid not in best_hit or hit[1] > best_hit[qid][1]:
                best_hit[qid] = hit
                classify_round[qid] = round_label
                newly_found += 1

        logger.info(
            "[PRECLASSIFY] %s: %d/%d newly classified (qcov=%s, maxrejects=%d)",
            round_label, newly_found, n_sub,
            f"{retry_qcov:.0%}" if retry_qcov > 0 else "none",
            retry_rejects,
        )

    logger.info(
        "[PRECLASSIFY] All rounds complete: %d/%d classified, %d remain unclassified",
        len(best_hit), n_input, n_input - len(best_hit),
    )

    # ── Parse vsearch matches — select best hit per query ───────────────────
    if taxa_tsv and _is_fasta_file(taxa_tsv):
        logger.info(
            "[PRECLASSIFY] --taxa points to a FASTA file (%s); "
            "taxonomy will be parsed from reference FASTA headers instead.",
            taxa_tsv,
        )
        effective_taxa_tsv = None
    else:
        effective_taxa_tsv = taxa_tsv

    taxa_map = _load_taxa_map(effective_taxa_tsv)

    if not taxa_map:
        try:
            from relict.pipeline.classify import _load_taxa_map_from_reference_fasta
            taxa_map, _warn = _load_taxa_map_from_reference_fasta(ref_to_use)
            if taxa_map:
                logger.info(
                    "[PRECLASSIFY] Parsed taxonomy from reference FASTA headers (%d entries)",
                    len(taxa_map) // 3,
                )
        except Exception as exc:
            logger.warning("[PRECLASSIFY] Could not parse taxa from FASTA headers: %s", exc)

    # best_hit is already fully populated by the retry cascade above.
    # Resolve taxonomy for each best hit.
    classification_tsv = os.path.join(outdir, f"{dataset_name}_classification.tsv")
    taxonomy_tsv       = os.path.join(outdir, f"{dataset_name}_taxonomy.tsv")
    unclassified_tsv   = os.path.join(outdir, f"{dataset_name}_unclassified.tsv")
    low_conf_tsv       = os.path.join(outdir, f"{dataset_name}_low_confidence.tsv")

    # Resolve taxonomy for each best hit
    full_rows: List[Tuple[str, str, float, Optional[str], Optional[float]]] = []
    tax_rows:  List[Tuple[str, Optional[str], Optional[float]]] = []
    classified_ids: set = set()

    for qid, (sid, pct) in best_hit.items():
        tax, conf = _resolve_taxon(sid, taxa_map)
        if conf is None and pct:
            conf = round(pct / 100.0, 4)
        classified_ids.add(qid)
        full_rows.append((qid, sid, pct, tax, conf))
        tax_rows.append((qid, tax, conf))

    written = len(full_rows)

    # ── Write full classification TSV ────────────────────────────────────────
    with open(classification_tsv, "w", newline="") as cf:
        cf.write("ID\tBestHit\tIdentity\tTaxon\tConfidence\tClassifiedRound\n")
        for qid, sid, pct, tax, conf in full_rows:
            cf.write(
                f"{qid}\t{sid}\t{pct}\t"
                f"{tax if tax is not None else 'NA'}\t"
                f"{conf if conf is not None else 'NA'}\t"
                f"{classify_round.get(qid, 'round1')}\n"
            )
    logger.info(
        "[PRECLASSIFY] Wrote %d classification rows for dataset '%s' → %s",
        written, dataset_name, classification_tsv,
    )

    # ── Write condensed taxonomy TSV ─────────────────────────────────────────
    with open(taxonomy_tsv, "w", newline="") as tf:
        tf.write("ID\tTaxon\tConfidence\n")
        for qid, tax, conf in tax_rows:
            tf.write(
                f"{qid}\t{tax if tax is not None else 'NA'}\t"
                f"{conf if conf is not None else 'NA'}\n"
            )
    logger.info("[PRECLASSIFY] Wrote condensed taxonomy → %s", taxonomy_tsv)

    # ── Write unclassified sequences TSV ────────────────────────────────────
    # Sequences with NO hit across ALL retry rounds.
    n_retry_rounds = len(_RETRY_ROUNDS)
    unclassified_ids = [
        (rid, input_ids.get(rid, 0))
        for rid in input_ids
        if rid not in classified_ids
    ]
    with open(unclassified_tsv, "w", newline="") as uf:
        uf.write(f"# Sequences from dataset '{dataset_name}' with no vsearch hit\n")
        uf.write(
            f"# Searched {1 + n_retry_rounds} rounds: "
            f"qcov 70%→50%→30%→none, "
            f"maxrejects {max_rejects}→{max_rejects*2}→{max_rejects*4}→{max_rejects*8}\n"
        )
        uf.write(f"# min_identity={min_identity}\n")
        uf.write(
            "# Possible causes: chimeric sequence, non-16S contamination, "
            "extreme divergence from reference, or excessive ambiguous bases.\n"
        )
        uf.write("ID\tLength\tNote\n")
        for rid, rlen in sorted(unclassified_ids, key=lambda x: x[0]):
            uf.write(
                f"{rid}\t{rlen}\t"
                f"No hit after {1 + n_retry_rounds} vsearch rounds "
                f"(id≥{min_identity:.0%}, qcov relaxed to none)\n"
            )
    n_unclassified = len(unclassified_ids)
    logger.info(
        "[PRECLASSIFY] %d/%d sequences unclassified after all rounds → %s",
        n_unclassified, n_input, unclassified_tsv,
    )

    # ── Write low-confidence sequences TSV ──────────────────────────────────
    low_conf_rows = [
        (qid, sid, pct, tax, conf)
        for qid, sid, pct, tax, conf in full_rows
        if pct < (low_confidence_threshold * 100)
    ]
    with open(low_conf_tsv, "w", newline="") as lf:
        lf.write(
            f"# Sequences from dataset '{dataset_name}' classified below "
            f"{low_confidence_threshold * 100:.1f}% identity\n"
        )
        lf.write(
            "# These hits are ambiguous at species level and may require "
            "manual review or lower-rank assignment.\n"
        )
        lf.write("ID\tBestHit\tIdentity\tTaxon\tConfidence\tClassifiedRound\n")
        for qid, sid, pct, tax, conf in sorted(low_conf_rows, key=lambda x: x[2]):
            lf.write(
                f"{qid}\t{sid}\t{pct}\t"
                f"{tax if tax is not None else 'NA'}\t"
                f"{conf if conf is not None else 'NA'}\t"
                f"{classify_round.get(qid, 'round1')}\n"
            )
    n_low_conf = len(low_conf_rows)
    logger.info(
        "[PRECLASSIFY] %d sequences have low-confidence hits (<%d%%) → %s",
        n_low_conf, int(low_confidence_threshold * 100), low_conf_tsv,
    )

    return {
        "classification_tsv": classification_tsv,
        "taxonomy_tsv":       taxonomy_tsv,
        "unclassified_tsv":   unclassified_tsv,
        "low_confidence_tsv": low_conf_tsv,
        "n_input":            n_input,
        "n_classified":       written,
        "n_unclassified":     n_unclassified,
        "n_low_confidence":   n_low_conf,
        "low_confidence_threshold": low_confidence_threshold,
        "min_identity":       min_identity,
        "max_hits":           max_hits,
        "max_rejects":        max_rejects,
    }


# ---------------------------------------------------------------------------
# Multi-dataset orchestration
# ---------------------------------------------------------------------------

DatasetSpec = Tuple[str, str]  # (dataset_name, fasta_path)


def run_preclassify(
    datasets: List[DatasetSpec],
    ref_fasta: str,
    outdir: str,
    taxa_tsv: Optional[str] = None,
    threads: int = 4,
    min_identity: float = 0.80,
    low_confidence_threshold: float = 0.97,
    max_hits: int = 10,
    max_rejects: int = 256,
) -> str:
    """Classify multiple FASTA datasets and produce merged outputs.

    Parameters
    ----------
    datasets:
        List of ``(dataset_name, fasta_path)`` pairs.
    ref_fasta:
        Path to the reference FASTA used for vsearch.
    outdir:
        Directory where all outputs are written.
    taxa_tsv:
        Optional taxonomy annotation file (FeatureID\\tTaxon\\tConfidence).
    threads:
        vsearch thread count.
    min_identity:
        Minimum vsearch identity threshold (0–1, default 0.80 = 80 %).
    low_confidence_threshold:
        Identity threshold below which a hit is flagged as low-confidence
        (0–1, default 0.97 = 97 %).  Sequences below this are written to
        ``{dataset}_low_confidence.tsv`` for manual review.
    max_hits:
        Maximum candidate hits per query that vsearch collects
        (``--maxaccepts`` / ``--maxhits``).  The best hit (highest %
        identity) is selected after collection.  Default: 10.
    max_rejects:
        ``--maxrejects`` passed to vsearch.  Raising this above vsearch's
        default of 32 prevents premature abandonment of ambiguous sequences
        (e.g. those with many N's).  Default: 256.

    Returns
    -------
    str
        Path to the pipeline-ready combined taxonomy TSV.
    """
    Path(outdir).mkdir(parents=True, exist_ok=True)

    results: Dict[str, dict] = {}  # name → classify_fasta result dict
    dataset_metas: Dict[str, dict] = {}

    for ds_name, fasta_path in datasets:
        meta = resolve_dataset_meta(ds_name)
        dataset_metas[ds_name] = meta
        logger.info(
            "[PRECLASSIFY] Processing dataset '%s' (%s) from %s",
            ds_name, meta["label"], fasta_path,
        )
        try:
            res = classify_fasta(
                fasta_path=fasta_path,
                ref_fasta=ref_fasta,
                outdir=outdir,
                dataset_name=ds_name,
                taxa_tsv=taxa_tsv,
                threads=threads,
                min_identity=min_identity,
                low_confidence_threshold=low_confidence_threshold,
                max_hits=max_hits,
                max_rejects=max_rejects,
            )
            results[ds_name] = res
        except Exception as exc:
            logger.error("[PRECLASSIFY] Failed to classify dataset '%s': %s", ds_name, exc)

    # Build combined taxonomy TSV (with Dataset column — human-readable)
    combined_path = os.path.join(outdir, "combined_taxonomy.tsv")
    combined_rows: List[Tuple[str, str, str, str]] = []  # id, dataset, taxon, conf
    seen_ids: set = set()
    for ds_name, res in results.items():
        tax_tsv = res["taxonomy_tsv"]
        try:
            with open(tax_tsv) as fh:
                next(fh, None)  # skip header
                for line in fh:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 1 or not parts[0]:
                        continue
                    qid = parts[0]
                    tax = parts[1] if len(parts) > 1 else "NA"
                    conf = parts[2] if len(parts) > 2 else "NA"
                    if qid not in seen_ids:
                        seen_ids.add(qid)
                        combined_rows.append((qid, ds_name, tax, conf))
        except Exception as exc:
            logger.warning("[PRECLASSIFY] Could not read taxonomy for dataset '%s': %s", ds_name, exc)

    with open(combined_path, "w", newline="") as fh:
        fh.write("ID\tDataset\tTaxon\tConfidence\n")
        for qid, ds_name, tax, conf in combined_rows:
            fh.write(f"{qid}\t{ds_name}\t{tax}\t{conf}\n")
    logger.info("[PRECLASSIFY] Wrote combined taxonomy (%d rows) → %s", len(combined_rows), combined_path)

    # Pipeline-compatible taxonomy (without Dataset column)
    # Pass directly as --taxa-assignments to ``relict preload`` or ``relict run``.
    pipeline_combined_path = os.path.join(outdir, "pipeline_taxonomy.tsv")
    with open(pipeline_combined_path, "w", newline="") as fh:
        fh.write("ID\tTaxon\tConfidence\n")
        for qid, _, tax, conf in combined_rows:
            fh.write(f"{qid}\t{tax}\t{conf}\n")
    logger.info("[PRECLASSIFY] Wrote pipeline-ready taxonomy → %s", pipeline_combined_path)

    # Dataset summary CSV
    summary_csv_path = os.path.join(outdir, "classification_summary.csv")
    _write_dataset_summary_csv(
        csv_path=summary_csv_path,
        datasets=datasets,
        results=results,
        dataset_metas=dataset_metas,
        low_confidence_threshold=low_confidence_threshold,
    )
    logger.info("[PRECLASSIFY] Wrote dataset summary CSV → %s", summary_csv_path)

    # Human-readable summary
    summary_path = os.path.join(outdir, "preclassify_summary.txt")
    _write_summary(
        summary_path=summary_path,
        datasets=datasets,
        results=results,
        combined_rows=combined_rows,
        dataset_metas=dataset_metas,
        ref_fasta=ref_fasta,
        taxa_tsv=taxa_tsv,
        combined_path=combined_path,
        pipeline_combined_path=pipeline_combined_path,
        summary_csv_path=summary_csv_path,
        low_confidence_threshold=low_confidence_threshold,
    )
    logger.info("[PRECLASSIFY] Wrote human-readable summary → %s", summary_path)

    return pipeline_combined_path


# ---------------------------------------------------------------------------
# Dataset summary CSV writer
# ---------------------------------------------------------------------------

# Confidence tier boundaries (identity %)
# High   : ≥ high_threshold   (default 97 %) – reliable species-level assignment
# Medium : 90 % ≤ x < 97 %   – genus-level confidence
# Low    : min_id ≤ x < 90 % – family/order level, ambiguous
# No hit : no vsearch match after all retry rounds
_MEDIUM_CONF_FLOOR = 90.0  # lower boundary of the medium tier (%)


def _write_dataset_summary_csv(
    csv_path: str,
    datasets: List[DatasetSpec],
    results: Dict[str, dict],
    dataset_metas: Dict[str, dict],
    low_confidence_threshold: float = 0.97,
) -> None:
    """Write ``classification_summary.csv`` — one row per dataset.

    Columns
    -------
    Dataset, Label, InputFASTA, TotalInput,
    Classified, Unclassified,
    HighConfidence, MediumConfidence, LowConfidence,
    PctHighConfidence, PctMediumConfidence, PctLowConfidence,
    PctClassified, PctUnclassified,
    HighConfThreshold_pct, MediumConfFloor_pct, MinIdentity_pct,
    ClassificationTSV, UnclassifiedTSV, LowConfidenceTSV
    """
    import csv as _csv

    high_floor  = low_confidence_threshold * 100   # e.g. 97.0
    med_floor   = _MEDIUM_CONF_FLOOR               # 90.0

    header = [
        "Dataset", "Label", "InputFASTA",
        "TotalInput", "Classified", "Unclassified",
        "HighConfidence", "MediumConfidence", "LowConfidence",
        "PctHighConfidence", "PctMediumConfidence", "PctLowConfidence",
        "PctClassified", "PctUnclassified",
        "HighConfThreshold_pct", "MediumConfFloor_pct", "MinIdentity_pct",
        "ClassificationTSV", "UnclassifiedTSV", "LowConfidenceTSV",
    ]

    def _p(n: int, total: int) -> str:
        return f"{n / total * 100:.1f}" if total > 0 else "0.0"

    rows = []
    for ds_name, fasta_path in datasets:
        meta = dataset_metas.get(ds_name, {})
        label = meta.get("label", ds_name)

        if ds_name not in results:
            rows.append(
                [ds_name, label, fasta_path] +
                ["FAILED"] * (len(header) - 3)
            )
            continue

        res = results[ds_name]
        n_input   = int(res.get("n_input", 0) or 0)
        n_unclass = int(res.get("n_unclassified", 0) or 0)
        min_id    = float(res.get("min_identity", 0.80) or 0.80)
        class_tsv = res.get("classification_tsv", "")
        unc_tsv   = res.get("unclassified_tsv", "")
        lc_tsv    = res.get("low_confidence_tsv", "")

        # Count per-tier by reading the classification TSV identity column
        n_high = n_med = n_low = 0
        try:
            with open(str(class_tsv)) as fh:
                reader = _csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    try:
                        pct = float(row.get("Identity", 0) or 0)
                    except (ValueError, TypeError):
                        pct = 0.0
                    if pct >= high_floor:
                        n_high += 1
                    elif pct >= med_floor:
                        n_med += 1
                    else:
                        n_low += 1
        except Exception as exc:
            logger.warning("[PRECLASSIFY] summary CSV: could not read %s: %s", class_tsv, exc)

        n_classified = n_high + n_med + n_low

        rows.append([
            ds_name, label, fasta_path,
            n_input, n_classified, n_unclass,
            n_high, n_med, n_low,
            _p(n_high, n_input), _p(n_med, n_input), _p(n_low, n_input),
            _p(n_classified, n_input), _p(n_unclass, n_input),
            f"{high_floor:.0f}", f"{med_floor:.0f}", f"{min_id * 100:.0f}",
            class_tsv, unc_tsv, lc_tsv,
        ])

    with open(csv_path, "w", newline="") as fh:
        writer = _csv.writer(fh)
        writer.writerow(header)
        writer.writerows(rows)


# ---------------------------------------------------------------------------
# Summary writer
# ---------------------------------------------------------------------------

def _write_summary(
    summary_path: str,
    datasets: List[DatasetSpec],
    results: Dict[str, dict],
    combined_rows: list,
    dataset_metas: Dict[str, dict],
    ref_fasta: str,
    taxa_tsv: Optional[str],
    combined_path: str,
    pipeline_combined_path: str,
    summary_csv_path: str = "",
    low_confidence_threshold: float = 0.97,
) -> None:
    lines: List[str] = []
    lines.append("=" * 70)
    lines.append("Relict  –  Pre-classification summary")
    lines.append("=" * 70)
    lines.append(f"\nReference FASTA : {ref_fasta}")
    if taxa_tsv:
        lines.append(f"Taxa TSV        : {taxa_tsv}")
    # Pull vsearch settings from first successful result for display
    _first = next((r for r in results.values() if isinstance(r, dict)), {})
    lines.append(f"Candidate hits  : up to {_first.get('max_hits', 10)} per query (best selected)")
    lines.append(f"Max rejects     : {_first.get('max_rejects', 256)}")
    lines.append(f"Low-conf cutoff : <{low_confidence_threshold * 100:.0f}% identity")
    lines.append(f"\nDatasets processed: {len(datasets)}")
    lines.append("")

    for ds_name, fasta_path in datasets:
        meta = dataset_metas.get(ds_name, {})
        label = meta.get("label", ds_name)
        desc = meta.get("description", "")
        url = meta.get("url", "")
        lines.append(f"  ┌─ Dataset : {ds_name}  ({label})")
        if desc:
            lines.append(f"  │  {textwrap.fill(desc, width=64, subsequent_indent='  │  ')}")
        if url:
            lines.append(f"  │  URL: {url}")
        lines.append(f"  │  Input FASTA : {fasta_path}")

        if ds_name in results:
            res = results[ds_name]
            n_input      = res.get("n_input", 0)
            n_classified = res.get("n_classified", 0)
            n_unclass    = res.get("n_unclassified", 0)
            n_low        = res.get("n_low_confidence", 0)
            class_tsv    = res.get("classification_tsv", "")
            tax_tsv_path = res.get("taxonomy_tsv", "")
            unclass_tsv  = res.get("unclassified_tsv", "")
            low_conf_tsv = res.get("low_confidence_tsv", "")
            min_id_pct   = int(float(res.get("min_identity") or 0.80) * 100)
            lc_pct       = int(low_confidence_threshold * 100)

            # Derived counts — all percentages use n_input as denominator
            n_high = (n_classified - n_low) if isinstance(n_classified, int) and isinstance(n_low, int) else "?"
            _ni = n_input if isinstance(n_input, int) and n_input > 0 else None

            def _pct(n: object) -> str:
                if _ni and isinstance(n, int):
                    return f"{n / _ni * 100:.1f}%"
                return "?"

            lines.append(f"  │")
            lines.append(f"  │  Input sequences                    : {n_input}")
            lines.append(f"  │")

            # High-confidence — the "good" classifications
            if isinstance(n_high, int) and isinstance(n_classified, int):
                lines.append(
                    f"  │  ✓  High-confidence   (≥{lc_pct}% identity)  : "
                    f"{n_high}/{n_input}  ({_pct(n_high)})"
                )
            else:
                lines.append(f"  │  High-confidence (≥{lc_pct}% identity) : ?")

            # Medium-confidence (90–97 %)
            n_med_raw = res.get("n_medium_confidence")  # may be absent in older runs
            # Fall back to reading classification TSV if needed
            if n_med_raw is None:
                n_med_raw = 0
                try:
                    import csv as _csv2
                    with open(str(res.get("classification_tsv", ""))) as _cf:
                        for _row in _csv2.DictReader(_cf, delimiter="\t"):
                            try:
                                _p2 = float(_row.get("Identity", 0) or 0)
                            except (ValueError, TypeError):
                                _p2 = 0.0
                            if _MEDIUM_CONF_FLOOR <= _p2 < (lc_pct):
                                n_med_raw += 1
                except Exception:
                    pass
            n_med = n_med_raw
            if isinstance(n_med, int) and n_med > 0:
                lines.append(
                    f"  │  ~  Medium-confidence ({_MEDIUM_CONF_FLOOR:.0f}–{lc_pct}% identity): "
                    f"{n_med}/{n_input}  ({_pct(n_med)})  → {low_conf_tsv}"
                )
            elif isinstance(n_med, int):
                lines.append(
                    f"  │     Medium-confidence ({_MEDIUM_CONF_FLOOR:.0f}–{lc_pct}% identity): 0"
                )

            # Low-confidence (min_id – 90 %)
            n_low_strict = (n_low - n_med) if isinstance(n_low, int) and isinstance(n_med, int) else n_low
            if isinstance(n_low_strict, int) and n_low_strict > 0:
                lines.append(
                    f"  │  ⚠  Low-confidence    ({min_id_pct}–{_MEDIUM_CONF_FLOOR:.0f}% identity): "
                    f"{n_low_strict}/{n_input}  ({_pct(n_low_strict)})  → {low_conf_tsv}"
                )
            elif isinstance(n_low_strict, int):
                lines.append(
                    f"  │     Low-confidence    ({min_id_pct}–{_MEDIUM_CONF_FLOOR:.0f}% identity): 0"
                )

            # Unclassified — no hit at all across all retry rounds
            if isinstance(n_unclass, int) and n_unclass > 0:
                lines.append(
                    f"  │  ✗  Unclassified (no hit)               : "
                    f"{n_unclass}/{n_input}  ({_pct(n_unclass)})  → {unclass_tsv}"
                )
            elif isinstance(n_unclass, int):
                lines.append(
                    f"  │     Unclassified (no hit)               : 0"
                )

            lines.append(f"  │")
            lines.append(f"  │  Classification TSV   : {class_tsv}")
            lines.append(f"  │  Taxonomy TSV         : {tax_tsv_path}")
            if isinstance(n_unclass, int) and n_unclass > 0:
                lines.append(f"  │  Unclassified TSV     : {unclass_tsv}")
            if isinstance(n_low, int) and n_low > 0:
                lines.append(f"  │  Low-confidence TSV   : {low_conf_tsv}")
        else:
            lines.append("  │  [classification FAILED – see log]")

        lines.append("  └─")
        lines.append("")

    # Overall totals — all percentages over total input
    total_input   = sum(int(r.get("n_input", 0))          for r in results.values() if isinstance(r, dict))
    total_class   = sum(int(r.get("n_classified", 0))     for r in results.values() if isinstance(r, dict))
    total_unclass = sum(int(r.get("n_unclassified", 0))   for r in results.values() if isinstance(r, dict))
    total_low_all = sum(int(r.get("n_low_confidence", 0)) for r in results.values() if isinstance(r, dict))

    # Compute medium/low split from classification TSVs
    total_high = 0
    total_med  = 0
    total_low_strict = 0
    high_floor_pct = low_confidence_threshold * 100
    for r in results.values():
        if not isinstance(r, dict):
            continue
        try:
            import csv as _csv3
            with open(str(r.get("classification_tsv", ""))) as _cf:
                for _row in _csv3.DictReader(_cf, delimiter="\t"):
                    try:
                        _p3 = float(_row.get("Identity", 0) or 0)
                    except (ValueError, TypeError):
                        _p3 = 0.0
                    if _p3 >= high_floor_pct:
                        total_high += 1
                    elif _p3 >= _MEDIUM_CONF_FLOOR:
                        total_med += 1
                    else:
                        total_low_strict += 1
        except Exception:
            pass

    lc_pct_global = int(low_confidence_threshold * 100)

    def _gpct(n: int) -> str:
        return f"{n / total_input * 100:.1f}%" if total_input > 0 else "?"

    lines.append("─" * 70)
    lines.append("Overall totals across all datasets")
    lines.append("─" * 70)
    lines.append(f"  Input sequences                              : {total_input}")
    lines.append(
        f"  ✓  High-confidence   (≥{lc_pct_global}% identity)     : "
        f"{total_high}  ({_gpct(total_high)})"
    )
    lines.append(
        f"  ~  Medium-confidence ({_MEDIUM_CONF_FLOOR:.0f}–{lc_pct_global}% identity) : "
        f"{total_med}  ({_gpct(total_med)})"
        + ("  ← review *_low_confidence.tsv" if total_med else "")
    )
    lines.append(
        f"  ⚠  Low-confidence    (<{_MEDIUM_CONF_FLOOR:.0f}% identity)     : "
        f"{total_low_strict}  ({_gpct(total_low_strict)})"
        + ("  ← review *_low_confidence.tsv" if total_low_strict else "")
    )
    if total_unclass:
        lines.append(
            f"  ✗  Unclassified     (no hit, all rounds)     : "
            f"{total_unclass}  ({_gpct(total_unclass)})  ← review *_unclassified.tsv"
        )
    else:
        lines.append(f"  ✓  Unclassified     (no hit, all rounds)     : 0")
    lines.append("")

    lines.append(f"Combined taxonomy (with Dataset column)  : {combined_path}")
    lines.append(f"Pipeline-ready taxonomy (ID,Taxon,Conf)  : {pipeline_combined_path}")
    lines.append(f"Dataset summary CSV                      : {summary_csv_path}")
    lines.append(f"Total rows in combined taxonomy          : {len(combined_rows)}")
    lines.append("")
    lines.append("─" * 70)
    lines.append("How to use these outputs with Relict")
    lines.append("─" * 70)
    lines.append(
        "\n1. Preload a classified dataset (taxonomy pre-computed, no re-classification):\n"
        "   relict preload \\\n"
        "     --fasta <original_fasta.fasta> \\\n"
        f"     --taxa-assignments {pipeline_combined_path} \\\n"
        "     --db my_project.db \\\n"
        "     --dataset <dataset_name> \\\n"
        "     -o preload_out/\n"
    )
    lines.append(
        "2. Run the main pipeline using the pre-classified taxonomy:\n"
        "   relict run \\\n"
        "     --input my_samples.fasta \\\n"
        "     --ref <ref.fasta> \\\n"
        f"     --taxa-assignments {pipeline_combined_path} \\\n"
        "     --db my_project.db \\\n"
        "     --dataset Batch1 \\\n"
        "     -o run_out/\n"
        "   (The --taxa-assignments flag prevents on-the-fly re-classification.)\n"
    )
    lines.append("=" * 70)

    with open(summary_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")

