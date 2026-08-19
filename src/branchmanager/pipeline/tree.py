"""
tree.py — Phylogenetic tree construction and incremental updating for BranchManager.

Architecture
------------
The tree is built from THREE layers of sequences:

  1. Reference anchors  (reference_anchors.fasta, bundled with BranchManager)
     Curated NCBI RefSeq 16S sequences spanning major rumen/gut phyla.
     These constrain the topology so the tree reflects real taxonomy.
     They are NEVER shown in outputs - they are scaffolding only.
     Identified by the BRANCHMANAGER_REF_ prefix in their headers.

  2. Filing Cabinet sequences  (e.g. Hungate1000)
     User-supplied baseline dataset registered via `branchmanager filing-cabinet`.
     Stored in the DB with dataset label. Shown in iTOL.

  3. Performance Review sequences
     User query sequences submitted via `branchmanager performance-review`. Shown in iTOL.

Tree building strategy
----------------------
  - Initial build: align all three layers together with mafft --auto
  - Incremental update: add new run sequences via mafft --addfragments
    against the existing combined alignment (which already contains anchors)
  - Every build labels internal nodes as NODE0000, NODE0001, ...
    (bootstrap values are stripped from node names to keep them clean)

Reference anchor FASTA
----------------------
The bundled anchor file lives at:
    <package_root>/data/reference_anchors.fasta

Each sequence header must be:
    >BRANCHMANAGER_REF_<PhylumName> accession=<ACC> source=<SOURCE>

Sequences in that file are excluded from all result outputs, novelty
scoring, and iTOL files. They exist purely to constrain tree topology.
"""

from __future__ import annotations

import hashlib
import logging
import shlex
import os as _os
import re
import shutil
import shutil as _shutil
import sys as _sys
from collections import defaultdict
from pathlib import Path
from typing import Optional

from branchmanager.utils.fasta import read_fasta, reverse_complement, write_fasta
from branchmanager.utils.subprocess import run_cmd

logger = logging.getLogger(__name__)


# Prefix that marks a sequence as a reference anchor.
REF_ANCHOR_PREFIX = "BRANCHMANAGER_REF_"

# Path to the bundled anchor file, relative to this module's location
_MODULE_DIR = Path(__file__).resolve().parent
BUILTIN_ANCHOR_FILE = _MODULE_DIR.parent / "data" / "reference_anchors.fasta"


# Public helpers

def is_ref_anchor(seq_id: str) -> bool:
    """Return True if a sequence ID belongs to a reference anchor."""
    return str(seq_id).startswith(REF_ANCHOR_PREFIX)


def filter_anchors_from_output(pairs: list[tuple]) -> list[tuple]:
    """Remove reference anchor entries from any (id, ...) pair list."""
    return [p for p in pairs if not is_ref_anchor(str(p[0]))]


# Reference anchor management

def get_anchor_file(custom_anchor_file: Optional[str] = None) -> Optional[Path]:
    """
    Return the path to the anchor FASTA to use.
    Priority:
      1. custom_anchor_file argument (user-supplied)
      2. BUILTIN_ANCHOR_FILE (bundled with BranchManager)
    Returns None if neither exists.
    """
    if custom_anchor_file:
        p = Path(custom_anchor_file)
        if p.exists():
            return p
        logger.warning("[ANCHORS] Custom anchor file not found: %s", p)

    if BUILTIN_ANCHOR_FILE.exists():
        return BUILTIN_ANCHOR_FILE

    logger.warning(
        "[ANCHORS] No reference anchor file found. Tree topology may not "
        "reflect true taxonomy. Consider running: branchmanager download-anchors"
    )
    return None


def load_anchor_sequences(anchor_file: Optional[Path]) -> list[tuple[str, str]]:
    """Load anchor sequences, adding the canonical anchor prefix when absent."""
    if anchor_file is None:
        return []
    records = []
    for header, seq in read_fasta(str(anchor_file)):
        if not is_ref_anchor(header):
            # Auto-prefix headers that are missing it so the file is forgiving
            header = f"{REF_ANCHOR_PREFIX}{header}"
        records.append((header, seq))
    logger.info("[ANCHORS] Loaded %d reference anchor sequences from %s", len(records), anchor_file)
    return records


def write_anchor_fasta(outdir: Path, anchor_file: Optional[Path]) -> Optional[Path]:
    """Write anchor sequences to outdir/ref_anchors.fasta. Returns path or None."""
    records = load_anchor_sequences(anchor_file)
    if not records:
        return None
    out_path = outdir / "ref_anchors.fasta"
    write_fasta(records, str(out_path))
    return out_path


# Combined FASTA construction

def build_combined_fasta(
    user_fasta: str,
    outdir: Path,
    anchor_file: Optional[Path] = None,
    db=None,
    db_dataset: Optional[str] = None,
    ref_fasta: Optional[str] = None,
    threads: int = 4,
    orientation_summary_path: Optional[Path] = None,
    build_mode: str = 'initial',
) -> Path:
    """
    Build a combined FASTA containing:
      - Reference anchor sequences (BRANCHMANAGER_REF_* prefixed)
      - All sequences from the DB that passed QC (if db provided)
      - User/query sequences from user_fasta

    Deduplicates by sequence content (MD5). The anchor sequences always
    take priority — if an identical sequence appears in user data it is
    dropped in favour of the anchor version.

    Returns path to combined FASTA.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    combined_path = outdir / "combined_input.fasta"

    seen_hashes: set[str] = set()
    all_records: list[tuple[str, str]] = []

    # 1. Anchors first — they set the reference topology
    anchors = load_anchor_sequences(anchor_file)
    for h, s in anchors:
        h_md5 = hashlib.md5(s.upper().encode()).hexdigest()
        seen_hashes.add(h_md5)
        all_records.append((h, s))

    # 2. DB sequences (Filing Cabinet and previous reviews)
    if db is not None:
        db_fasta = outdir / "db_sequences.fasta"
        try:
            got = db.get_sequences_fasta(str(db_fasta), dataset=db_dataset)
            if got and db_fasta.exists():
                logger.info("[TREE] Successfully exported %d DB sequences", sum(1 for _ in read_fasta(str(db_fasta))))
        except Exception as e:
            logger.warning("[TREE] Failed to get DB sequences: %s", e)
            got = False
        if got and db_fasta.exists():
            db_fasta_to_use, db_orientation_rows = _orient_tree_input_fasta(
                str(db_fasta),
                ref_fasta=ref_fasta,
                anchor_fasta=str(anchor_file) if anchor_file else None,
                outdir=outdir,
                label='db_sequences',
                build_mode=build_mode,
                threads=threads,
            )
            if orientation_summary_path and db_orientation_rows:
                _write_tree_orientation_summary(orientation_summary_path, db_orientation_rows, append=orientation_summary_path.exists())
            for h, s in read_fasta(str(db_fasta_to_use)):
                h_md5 = hashlib.md5(s.upper().encode()).hexdigest()
                if h_md5 not in seen_hashes:
                    seen_hashes.add(h_md5)
                    all_records.append((h, s))

    # 3. User sequences
    user_fasta_to_use, user_orientation_rows = _orient_tree_input_fasta(
        user_fasta,
        ref_fasta=ref_fasta,
        anchor_fasta=str(anchor_file) if anchor_file else None,
        outdir=outdir,
        label='user_sequences',
        build_mode=build_mode,
        threads=threads,
    )
    user_seq_count = sum(1 for _ in read_fasta(str(user_fasta_to_use)))
    logger.info("[TREE] User sequences after orientation: %d sequences", user_seq_count)
    if orientation_summary_path and user_orientation_rows:
        _write_tree_orientation_summary(orientation_summary_path, user_orientation_rows, append=orientation_summary_path.exists())
    for h, s in read_fasta(user_fasta_to_use):
        h_md5 = hashlib.md5(s.upper().encode()).hexdigest()
        if h_md5 not in seen_hashes:
            seen_hashes.add(h_md5)
            all_records.append((h, s))

    # Count the sources for logging
    all_records_map = [(h, is_ref_anchor(h)) for h, _ in all_records]
    n_anchors = sum(1 for _, is_anchor in all_records_map if is_anchor)
    n_data = len(all_records_map) - n_anchors

    write_fasta(all_records, str(combined_path))
    logger.info(
        "[TREE] Combined FASTA: %d sequences total (%d anchors, %d data sequences [DB+user])",
        len(all_records), n_anchors, n_data,
    )
    return combined_path


# Duplicate header handling

def _sanitise_name(x: str) -> str:
    token = str(x).split()[0]
    if "|" in token:
        token = token.split("|")[-1]
    token = re.sub(r"[^0-9A-Za-z_.\-]+", "_", token).strip("_")
    if not token:
        token = hashlib.md5(str(x).encode()).hexdigest()[:8]
    return token


def _make_unique_fasta(fasta_path: str, outdir: Path) -> tuple[str, dict]:
    """
    Ensure all headers in fasta_path are unique and Newick-safe.
    Returns (path_to_use, mapping) where mapping is {unique_id: original_id}.
    If no duplicates exist returns (fasta_path, {}).
    """
    records = list(read_fasta(fasta_path))
    seen: dict[str, int] = defaultdict(int)
    mapping: dict[str, str] = {}
    new_records: list[tuple[str, str]] = []
    has_duplicates = False

    for h, seq in records:
        base = _sanitise_name(h)
        if seen[base] == 0:
            new_h = base
        else:
            has_duplicates = True
            new_h = f"{base}__dup{seen[base]}"
        seen[base] += 1
        mapping[new_h] = h
        new_records.append((new_h, seq))

    if not has_duplicates:
        return fasta_path, {}

    outdir.mkdir(parents=True, exist_ok=True)
    unique_path = outdir / (Path(fasta_path).stem + "_unique.fasta")
    write_fasta(new_records, str(unique_path))

    with open(outdir / "id_map.tsv", "w") as fh:
        fh.write("unique_id\toriginal_id\n")
        for new_id, orig in mapping.items():
            fh.write(f"{new_id}\t{orig}\n")

    return str(unique_path), mapping


# Internal node labelling

def _label_internal_nodes(newick_text: str, prefix: str = "NODE") -> str:
    """
    Replace unlabelled internal nodes (and bootstrap values) with clean
    deterministic labels: NODE0000, NODE0001, ...

    FastTree writes bootstrap values after closing parentheses like:
        (...) 0.963 :0.05
    This function replaces those numeric-only labels with clean NODE ids
    so that iTOL can reference them unambiguously.
    """
    existing = set(re.findall(r"\b" + re.escape(prefix) + r"\d+\b", newick_text))
    counter = 0

    def _next_name() -> str:
        nonlocal counter
        while True:
            name = f"{prefix}{counter:04d}"
            counter += 1
            if name not in existing:
                existing.add(name)
                return name

    # Match ')' then optional numeric bootstrap value, then consume the
    # next delimiter once to avoid duplicating it in the output.
    def _replace(m: re.Match) -> str:
        delim = m.group(1)
        return ")" + _next_name() + delim

    pattern = re.compile(r"\)\s*(?:[0-9]+(?:\.[0-9]+)?)?\s*([:\),;])")

    try:
        return pattern.sub(_replace, newick_text)
    except Exception:
        return newick_text


def _repair_internal_node_label_delimiters(newick_text: Optional[str]) -> str:
    """Repair malformed node labels such as NODE0001::0.03 (duplicated delimiter)."""
    if not newick_text:
        return ''
    return re.sub(r"\b(NODE\d+)(::+)", r"\1:", newick_text)


# Normalisation helpers

def _norm_id(x: str) -> str:
    """Return the persistent FASTA ID without older shortening heuristics."""
    return str(x).strip() if x is not None else ''


def _new_sequences_only(
    user_fasta: str,
    base_aln: Path,
    db,
    outdir: Path,
) -> list[tuple[str, str]]:
    """
    Filter user_fasta to only sequences not already in base_aln (the tree).

    By default, all sequences not yet in the tree are added, even if they
    exist in the DB from previous runs. This ensures new submissions are
    always visible in the tree.

    Anchor sequences in base_aln are always excluded from this comparison
    (they are added separately during build_combined_fasta).
    """
    base_ids_norm: set[str] = set()
    try:
        for h, _ in read_fasta(str(base_aln)):
            if not is_ref_anchor(h):
                base_ids_norm.add(_norm_id(h))
    except Exception:
        pass

    new_records = []
    for h, s in read_fasta(user_fasta):
        if _norm_id(h) not in base_ids_norm:
            new_records.append((h, s))

    return new_records


def _ensure_uncompressed_reference(ref_fasta: Optional[str], outdir: Path) -> Optional[str]:
    if not ref_fasta:
        return None
    if not str(ref_fasta).endswith('.gz'):
        return ref_fasta
    out = outdir / 'tree_orientation_ref.fasta'
    if not out.exists():
        write_fasta(list(read_fasta(ref_fasta)), str(out))
    return str(out)


def _parse_orientation_by_query(matches_path: str):
    orientations = {}
    best_identity = {}
    try:
        with open(matches_path) as handle:
            for line in handle:
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 8:
                    continue
                qid = parts[0]
                try:
                    pct = float(parts[2])
                    qstart = int(parts[6])
                    qend = int(parts[7])
                except Exception:
                    continue
                if qid not in best_identity or pct > best_identity[qid]:
                    best_identity[qid] = pct
                    orientations[qid] = {
                        'orientation': 'reverse' if qstart > qend else 'forward',
                        'best_hit': parts[1] if len(parts) > 1 else 'NA',
                        'best_identity': pct,
                        'qstart': qstart,
                        'qend': qend,
                    }
    except FileNotFoundError:
        return {}
    return orientations


def _write_tree_orientation_summary(path: Path, rows, *, append: bool = False):
    mode = 'a' if append and path.exists() else 'w'
    with open(path, mode) as handle:
        if mode == 'w':
            handle.write(
                'BuildMode\tSourceGroup\tSequenceID\tOrientationCall\tReverseComplemented\t'
                'BestHit\tBestIdentity\tQueryStart\tQueryEnd\tReferenceSource\tStatus\n'
            )
        for row in rows:
            handle.write(
                f"{row.get('build_mode', '')}\t{row.get('source_group', '')}\t{row.get('sequence_id', '')}\t"
                f"{row.get('orientation_call', '')}\t{row.get('reverse_complemented', '')}\t"
                f"{row.get('best_hit', 'NA')}\t{row.get('best_identity', 'NA')}\t"
                f"{row.get('qstart', 'NA')}\t{row.get('qend', 'NA')}\t"
                f"{row.get('reference_source', '')}\t{row.get('status', '')}\n"
            )
    return str(path)


def _orient_tree_input_fasta(
    input_fasta: str,
    *,
    ref_fasta: Optional[str],
    anchor_fasta: Optional[str],
    outdir: Path,
    label: str,
    build_mode: str = 'manual',
    threads: int = 4,
) -> tuple[str, list[dict[str, object]]]:
    """Orient tree-input sequences against the reference/anchor database.

    This is tree-only normalisation: stored DB sequences are not mutated.
    Sequences without confident orientation evidence are left unchanged.
    """
    ref_to_use = _ensure_uncompressed_reference(ref_fasta, outdir)
    if not ref_to_use and anchor_fasta:
        ref_to_use = str(anchor_fasta)
    records = list(read_fasta(input_fasta))
    data_records = [(h, s) for h, s in records if not is_ref_anchor(h)]
    if not data_records:
        return input_fasta, []

    if not ref_to_use:
        return input_fasta, [
            {
                'build_mode': build_mode,
                'source_group': label,
                'sequence_id': h,
                'orientation_call': 'unknown',
                'reverse_complemented': 'False',
                'best_hit': 'NA',
                'best_identity': 'NA',
                'qstart': 'NA',
                'qend': 'NA',
                'reference_source': 'none',
                'status': 'no_reference',
            }
            for h, _ in data_records
        ]

    tmp_matches = outdir / f'{label}_orientation_matches.tsv'
    thread_flag = f' --threads {threads}' if threads and int(threads) > 0 else ''
    cmd = (
        f'vsearch --usearch_global {shlex.quote(str(input_fasta))} '
        f'--db {shlex.quote(str(ref_to_use))} '
        f'--id 0.5 '
        f'--strand both '
        f'--blast6out {shlex.quote(str(tmp_matches))} '
        f'--maxaccepts 1 --maxhits 1 '
        f'--query_cov 0.5{thread_flag}'
    )

    try:
        run_cmd(cmd)
        orientations = _parse_orientation_by_query(str(tmp_matches))
    except Exception as e:
        raise RuntimeError(
            f'Tree-input orientation search failed for {input_fasta}: {e}'
        ) from e
    finally:
        try:
            if tmp_matches.exists():
                tmp_matches.unlink()
        except Exception:
            pass

    flipped = 0
    oriented_records = []
    summary_rows = []
    for header, seq in records:
        if is_ref_anchor(header):
            oriented_records.append((header, seq))
            continue
        hit = orientations.get(header)
        if hit and hit.get('orientation') == 'reverse':
            oriented_records.append((header, reverse_complement(seq)))
            flipped += 1
            status = 'reverse_flipped'
            orientation_call = 'reverse'
            reverse_comped = 'True'
        else:
            oriented_records.append((header, seq))
            if hit:
                status = 'forward_kept'
                orientation_call = 'forward'
            else:
                status = 'no_hit'
                orientation_call = 'unknown'
            reverse_comped = 'False'
        summary_rows.append({
            'build_mode': build_mode,
            'source_group': label,
            'sequence_id': header,
            'orientation_call': orientation_call,
            'reverse_complemented': reverse_comped,
            'best_hit': hit.get('best_hit', 'NA') if hit else 'NA',
            'best_identity': f"{hit.get('best_identity'):.2f}" if hit and hit.get('best_identity') is not None else 'NA',
            'qstart': hit.get('qstart', 'NA') if hit else 'NA',
            'qend': hit.get('qend', 'NA') if hit else 'NA',
            'reference_source': Path(ref_to_use).name,
            'status': status,
        })

    if flipped == 0:
        logger.info('[TREE] Orientation check for %s found no reverse-strand sequences', input_fasta)
        return input_fasta, summary_rows

    oriented_path = outdir / f'{Path(input_fasta).stem}_oriented.fasta'
    write_fasta(oriented_records, str(oriented_path))
    logger.info(
        '[TREE] Oriented %d/%d sequences for tree input %s using %s',
        flipped,
        len(data_records),
        input_fasta,
        ref_to_use,
    )
    return str(oriented_path), summary_rows


# Tree building

def _run_fasttree(aln_path: Path, tree_path: Path) -> bool:
    """Run FastTree -nt -gtr on aln_path, write to tree_path. Returns success."""
    cmd = f"FastTree -nt -gtr {shlex.quote(str(aln_path))} > {shlex.quote(str(tree_path))}"
    logger.info("[TREE] Running FastTree: %s", cmd)
    try:
        run_cmd(cmd)
        logger.info("[TREE] FastTree complete → %s", tree_path)
        return True
    except RuntimeError as e:
        logger.warning(
            "[TREE] FastTree failed: %s\n"
            "Install via: conda install -c bioconda fasttree", e
        )
        return False


def _resolve_iqtree_binary() -> Optional[str]:
    """Find the IQ-TREE binary, trying multiple candidate names.

    Conda environments sometimes install IQ-TREE as ``iqtree`` rather than
    ``iqtree2``, and ``/bin/sh`` may not inherit the activated conda PATH.
    This function:
      1. Searches all candidates via ``shutil.which`` (honours PATH).
      2. Falls back to ``dirname(sys.executable)`` — the conda env's bin dir —
         so the correct binary is always found when BranchManager is run from within
         the conda environment.
    """
    candidates = ["iqtree2", "iqtree", "IQ-TREE", "iqtree-omp"]
    for name in candidates:
        found = _shutil.which(name)
        if found:
            logger.debug("[TREE] Resolved IQ-TREE binary via PATH: %s", found)
            return found
    # Fall back to the bin directory that contains the current Python interpreter
    env_bin = _os.path.dirname(_sys.executable)
    for name in candidates:
        candidate = _os.path.join(env_bin, name)
        if _os.path.isfile(candidate) and _os.access(candidate, _os.X_OK):
            logger.debug("[TREE] Resolved IQ-TREE binary via env bin: %s", candidate)
            return candidate
    return None


def _resolve_tree_builder_binary(tree_method: str) -> tuple[Optional[str], tuple[str, ...]]:
    if tree_method in ("iqtree", "iqtree-fast"):
        return _resolve_iqtree_binary(), ("iqtree2", "iqtree", "IQ-TREE", "iqtree-omp")
    for candidate in ("FastTree", "fasttree"):
        found = _shutil.which(candidate)
        if found:
            return found, ("FastTree", "fasttree")
    return None, ("FastTree", "fasttree")


def preflight_tree_tools(tree_method: str = "fasttree", *, require_mafft: bool = True) -> None:
    """Fail fast when tree-building dependencies are unavailable."""
    details: list[str] = []

    if require_mafft and not _shutil.which("mafft"):
        details.append("mafft (mafft)")

    tree_binary, tree_candidates = _resolve_tree_builder_binary(tree_method)
    if not tree_binary:
        tool = "iqtree" if tree_method in ("iqtree", "iqtree-fast") else "FastTree"
        details.append(f"{tool} ({', '.join(tree_candidates)})")

    if details:
        raise RuntimeError(
            "Missing required external tools for tree building: "
            + "; ".join(details)
            + ". Install them and rerun the command."
        )


def _run_iqtree(
    aln_path: Path,
    tree_path: Path,
    threads: int = 4,
    model: str = "GTR+G+I",
    fast: bool = False,
) -> bool:
    """Run IQ-TREE 2/3 on *aln_path*, write the best tree to *tree_path*.

    IQ-TREE produces several output files under a common prefix; this
    function copies the ``.treefile`` to *tree_path* so the rest of the
    pipeline is unaffected.

    Parameters
    ----------
    fast  : When True, pass ``-fast`` for a quicker (but still much more
            accurate than FastTree) search.  Recommended for incremental runs.
    """
    binary = _resolve_iqtree_binary()
    if not binary:
        logger.warning(
            "[TREE] IQ-TREE binary not found in PATH or conda env bin dir.\n"
            "Install via: conda install -c bioconda iqtree"
        )
        return False

    prefix = str(aln_path.parent / aln_path.stem) + "_iqtree"
    thread_flag = f"-T {threads}" if threads and int(threads) > 0 else "-T AUTO"
    fast_flag = "-fast" if fast else ""
    cmd = (
        f"{shlex.quote(str(binary))} -s {shlex.quote(str(aln_path))} -m {model} {thread_flag} "
        f"-pre {shlex.quote(str(prefix))} --redo -quiet {fast_flag}"
    ).strip()
    logger.info("[TREE] Running IQ-TREE%s: %s", " (fast)" if fast else "", cmd)
    try:
        run_cmd(cmd)
    except RuntimeError as e:
        logger.warning(
            "[TREE] IQ-TREE failed: %s\n"
            "Install via: conda install -c bioconda iqtree", e,
        )
        return False

    treefile = Path(prefix + ".treefile")
    if not treefile.exists():
        logger.warning("[TREE] IQ-TREE did not produce a treefile at %s", treefile)
        return False

    _shutil.copyfile(treefile, tree_path)
    logger.info("[TREE] IQ-TREE complete → %s", tree_path)
    return True


def _build_tree(
    fasta_for_tree: str,
    tree_path: Path,
    threads: int = 4,
    method: str = "fasttree",
) -> bool:
    """Dispatch to the chosen tree-building backend.

    Parameters
    ----------
    method : ``'fasttree'``   — approximate ML, GTR+CAT (default, fast)
             ``'iqtree'``     — full ML, GTR+G+I, SPR+NNI (more accurate)
             ``'iqtree-fast'``— full ML with IQ-TREE ``-fast`` flag (compromise)
    """
    if method in ("iqtree", "iqtree-fast"):
        fast = method == "iqtree-fast"
        return _run_iqtree(Path(fasta_for_tree), tree_path, threads=threads, fast=fast)
    return _run_fasttree(Path(fasta_for_tree), tree_path)


def _run_mafft_full(input_fasta: Path, output_fasta: Path, threads: int = 4) -> bool:
    """Full de-novo MAFFT alignment (--auto mode)."""
    thread_flag = f"--thread {threads}" if threads > 0 else ""
    cmd = (
        f"mafft {thread_flag} --auto {shlex.quote(str(input_fasta))} "
        f"> {shlex.quote(str(output_fasta))}"
    ).strip()
    logger.info("[TREE] Running MAFFT (full): %s", cmd)
    try:
        run_cmd(cmd)
        return True
    except RuntimeError as e:
        logger.warning("[TREE] MAFFT failed: %s", e)
        return False


def _run_mafft_add(
    new_fasta: Path,
    backbone_aln: Path,
    output_fasta: Path,
    threads: int = 4,
) -> bool:
    """Add near-full-length sequences to an existing alignment via --add."""
    thread_flag = f"--thread {threads}" if threads > 0 else ""
    cmd = (
        f"mafft {thread_flag} --add {shlex.quote(str(new_fasta))} "
        f"{shlex.quote(str(backbone_aln))} > {shlex.quote(str(output_fasta))}"
    ).strip()
    logger.info("[TREE] Running MAFFT (add): %s", cmd)
    try:
        run_cmd(cmd)
        return True
    except RuntimeError as e:
        logger.warning("[TREE] MAFFT add failed: %s", e)
        return False


def _run_mafft_addfragments(
    new_fasta: Path,
    backbone_aln: Path,
    output_fasta: Path,
    threads: int = 4,
) -> bool:
    """Add new sequences to an existing alignment via --addfragments."""
    thread_flag = f"--thread {threads}" if threads > 0 else ""
    cmd = (
        f"mafft {thread_flag} --addfragments {shlex.quote(str(new_fasta))} "
        f"{shlex.quote(str(backbone_aln))} > {shlex.quote(str(output_fasta))}"
    ).strip()
    logger.info("[TREE] Running MAFFT (addfragments): %s", cmd)
    try:
        run_cmd(cmd)
        return True
    except RuntimeError as e:
        logger.warning("[TREE] MAFFT addfragments failed: %s", e)
        return False


def _choose_mafft_incremental_mode(new_records: list[tuple[str, str]], backbone_aln: Path) -> str:
    """Choose MAFFT incremental mode based on sequence lengths.

    Use --add for mostly full-length sequences so existing topology/columns are
    preserved more faithfully. Use --addfragments for shorter/fragmentary input.
    """
    if not new_records:
        return 'addfragments'
    try:
        backbone_lengths = [
            len(seq.replace('-', ''))
            for header, seq in read_fasta(str(backbone_aln))
            if seq and not is_ref_anchor(header)
        ]
    except Exception:
        backbone_lengths = []
    new_lengths = [len(seq) for _, seq in new_records if seq]
    if not new_lengths:
        return 'addfragments'
    median_backbone = sorted(backbone_lengths)[len(backbone_lengths) // 2] if backbone_lengths else 0
    median_new = sorted(new_lengths)[len(new_lengths) // 2]
    if min(new_lengths) >= 1000 and (not median_backbone or median_new >= 0.85 * median_backbone):
        return 'add'
    return 'addfragments'


def _seed_backbone_from_previous_review(previous_review: Optional[str], outdir: Path) -> bool:
    """Copy a previous-review backbone alignment/tree into a fresh run outdir if present."""
    if not previous_review:
        return False
    p = Path(previous_review)
    previous_alignment = p / 'current_alignment.fasta'
    if not previous_alignment.exists():
        previous_alignment = p / 'tree' / 'current_alignment.fasta'
    if not previous_alignment.exists():
        return False
    outdir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(previous_alignment, outdir / 'current_alignment.fasta')
    previous_tree = p / 'current_tree.nwk'
    if not previous_tree.exists():
        previous_tree = p / 'tree' / 'current_tree.nwk'
    if previous_tree.exists():
        try:
            shutil.copyfile(previous_tree, outdir / 'current_tree.nwk')
        except Exception:
            pass
    previous_orientation = p / 'tree_orientation_summary.tsv'
    if not previous_orientation.exists():
        previous_orientation = p / 'tree' / 'tree_orientation_summary.tsv'
    if previous_orientation.exists():
        try:
            shutil.copyfile(previous_orientation, outdir / 'tree_orientation_summary.tsv')
        except Exception:
            pass
    logger.info('[TREE] Seeded backbone alignment from previous review %s', p)
    return True


def collect_tree_build_warnings(user_fasta: str, anchor_file: Optional[str] = None, db=None, db_dataset: Optional[str] = None):
    warnings = []
    try:
        user_records = list(read_fasta(user_fasta))
    except Exception as e:
        return [{'category': 'USER_FASTA_READ_FAILED', 'subject_id': str(user_fasta), 'detail': str(e)}]

    if len(user_records) < 3:
        warnings.append({'category': 'LOW_SEQUENCE_COUNT', 'subject_id': str(len(user_records)), 'detail': 'Very few user sequences are available; topology will be weak and novelty calls unstable.'})

    short_records = [h for h, s in user_records if len(s) < 1200]
    if short_records:
        warnings.append({'category': 'PARTIAL_16S_SEQUENCES', 'subject_id': str(len(short_records)), 'detail': 'Sequences shorter than 1200 bp reduce phylogenetic signal and can produce unstable placement.'})

    n_rich = [h for h, s in user_records if s and (s.upper().count('N') / float(len(s))) > 0.01]
    if n_rich:
        warnings.append({'category': 'HIGH_N_CONTENT', 'subject_id': str(len(n_rich)), 'detail': 'Sequences with high ambiguous-base content can produce unreliable alignment and branch lengths.'})

    anchors = load_anchor_sequences(get_anchor_file(anchor_file))
    if len(anchors) < 2:
        warnings.append({'category': 'MISSING_REFERENCE_ANCHORS', 'subject_id': '0', 'detail': 'Tree was built without a meaningful anchor backbone; clades may be poorly constrained.'})

    return warnings


def summarise_alignment_quality(alignment_fasta: str):
    warnings = []
    p = Path(alignment_fasta)
    if not p.exists():
        return warnings
    try:
        records = list(read_fasta(str(p)))
    except Exception as e:
        return [{'category': 'ALIGNMENT_READ_FAILED', 'subject_id': str(alignment_fasta), 'detail': str(e)}]
    if not records:
        return [{'category': 'EMPTY_ALIGNMENT', 'subject_id': str(alignment_fasta), 'detail': 'Alignment file is empty.'}]

    all_gap = [h for h, s in records if s and set(s) <= {'-'}]
    if all_gap:
        warnings.append({'category': 'ALL_GAP_ALIGNMENT_ROWS', 'subject_id': str(len(all_gap)), 'detail': 'Some aligned rows are entirely gaps, indicating failed fragment placement.'})

    ungapped_lengths = [len(s.replace('-', '')) for _, s in records if s]
    if ungapped_lengths and min(ungapped_lengths) < 300:
        warnings.append({'category': 'VERY_SHORT_ALIGNED_FRAGMENTS', 'subject_id': str(min(ungapped_lengths)), 'detail': 'Some aligned sequences retain very little ungapped sequence and may place unreliably.'})
    return warnings


def write_tree_warning_tsv(outdir: str, warning_rows):
    path = Path(outdir) / 'tree_build_warnings.tsv'
    with open(path, 'w') as fh:
        fh.write('Category\tSubjectID\tDetail\n')
        for row in warning_rows:
            fh.write(f"{row.get('category', 'WARNING')}\t{row.get('subject_id', '')}\t{row.get('detail', '')}\n")
    return str(path)


def _prune_anchor_leaves(newick: str) -> str:
    """Remove all BranchManager reference-anchor leaf nodes from a Newick string.

    Anchor sequences constrain tree topology at *build time* (during MAFFT
    alignment and FastTree inference).  Once the tree has been written the
    topology of the data sequences is fixed — the anchors have done their job
    and can be stripped from the stored newick without any loss of information.

    Removing them here means:
      • iTOL displays only the data sequences — no mystery unlabelled tips.
      • All downstream tools that read current_tree.nwk see a clean tree.
      • The ``current_alignment.fasta`` **retains** the anchors so that
        incremental ``mafft --addfragments`` runs still benefit from the
        anchored backbone on the next update.

    Algorithm
    ---------
    Iteratively applies three regex passes until the string stabilises:
      1. ``,ANCHOR`` — anchor that is a non-first sibling
      2. ``ANCHOR,`` — anchor that is the first sibling (its comma is trailing)
      3. ``ANCHOR``  — sole child of a clade (degenerate after passes 1-2)
    Then cleans up any ``(,``, ``,)``, ``()``, or ``,,`` artefacts produced
    by the removals.

    In practice, with 26 anchors spread across hundreds of data sequences,
    the degenerate single-child case essentially never occurs.
    """
    _ANCHOR_PAT = r'BRANCHMANAGER_REF_[^,:()\s;]+(?::[0-9Ee.+\-]+)?'

    for _ in range(100):          # safety iteration cap
        before = newick
        newick = re.sub(r',\s*' + _ANCHOR_PAT, '', newick)   # ,ANCHOR
        newick = re.sub(_ANCHOR_PAT + r'\s*,', '', newick)   # ANCHOR,
        newick = re.sub(_ANCHOR_PAT, '', newick)              # lone ANCHOR
        # Clean artefacts.
        # IMPORTANT: remove empty clades *including* any trailing internal-node
        # label and branch-length that FastTree attaches to the closing ')'.
        # Without this, a sole-anchor clade like (BRANCHMANAGER_REF_X:0.5)1.000:0.09
        # becomes ()1.000:0.09 then just 1.000:0.09 which looks like a leaf.
        newick = re.sub(                                      # ()label:len or ()label
            r'\(\s*\)(?:[^,():;\s]+(?::[^,():;\s]+)?)?', '', newick
        )
        newick = re.sub(r'\(,', '(', newick)                  # (,X
        newick = re.sub(r',\)', ')', newick)                  # X,)
        newick = re.sub(r',,+', ',', newick)                  # ,,
        if newick == before:
            break

    return newick.strip()


def _finalise_tree(out: Path, id_map: dict) -> None:
    """Remap IDs, prune anchor leaves, and label internal nodes."""
    tree_path = out / "current_tree.nwk"
    if not tree_path.exists():
        return

    newick = tree_path.read_text()

    # Remap unique IDs back to original IDs if we created them
    if id_map:
        for new_id, orig_id in id_map.items():
            newick = newick.replace(new_id, orig_id)

    # Strip reference anchor leaves — they have done their job constraining
    # the topology and must not appear in iTOL or other visualisations.
    before_prune = newick.count(REF_ANCHOR_PREFIX)
    newick = _prune_anchor_leaves(newick)
    after_prune = newick.count(REF_ANCHOR_PREFIX)
    if before_prune:
        logger.info("[TREE] Pruned %d anchor leaf/leaves from stored newick "
                    "(%d unexpected occurrences remaining)",
                    before_prune - after_prune, after_prune)

    # Repair malformed labels first, then relabel internal nodes cleanly.
    newick = _repair_internal_node_label_delimiters(newick)
    newick = _label_internal_nodes(newick)
    tree_path.write_text(newick)
    logger.info("[TREE] Internal nodes labelled; tree finalised at %s", tree_path)


# Main entry point

def initialise_or_update_tree(
    ref_fasta: str,
    user_fasta: str,
    outdir: str,
    db=None,
    db_dataset: Optional[str] = None,
    threads: int = 4,
    anchor_file: Optional[str] = None,
    force_rebuild: bool = False,
    previous_review: Optional[str] = None,
    tree_method: str = "fasttree",
):
    """
    Initialise or incrementally update the BranchManager phylogenetic tree.

    Parameters
    ----------
    ref_fasta      : Path to reference FASTA used for classification (not used
                     for tree building; accepted for older callers).
    user_fasta     : FASTA of new query sequences to add.
    outdir         : Output directory. Persistent files:
                       current_alignment.fasta — the growing backbone
                       current_tree.nwk        — the current tree
    db             : Optional Database instance to pull existing sequences from.
    db_dataset     : If set, only pull this dataset from the DB.
    threads        : Number of threads for MAFFT and the tree builder.
    anchor_file    : Path to custom reference anchor FASTA.
                     Falls back to the bundled data/reference_anchors.fasta.
    force_rebuild  : If True, rebuild the entire tree from scratch even if
                     current_alignment.fasta already exists.
    tree_method    : Tree-building backend to use.
                     ``'fasttree'``    — approximate ML, GTR+CAT (default, fast).
                     ``'iqtree'``      — full ML, GTR+G+I + SPR/NNI (more accurate,
                                         slower; recommended for publication figures).
                     ``'iqtree-fast'`` — IQ-TREE 2 with ``-fast`` flag; a good
                                         compromise for incremental exploratory runs.

    Strategy
    --------
    FIRST RUN (no current_alignment.fasta):
      1. Load reference anchors
      2. Pull DB sequences
      3. Load user sequences
      4. Combine all three → combined_input.fasta
      5. MAFFT --auto on combined_input.fasta → current_alignment.fasta
      6. Tree builder → current_tree.nwk

    SUBSEQUENT RUNS (current_alignment.fasta exists):
      1. Find sequences in user_fasta not yet in alignment or DB
      2. If new sequences exist:
         a. MAFFT --addfragments new_seqs onto current_alignment.fasta
         b. Tree builder → current_tree.nwk
         c. Replace current_alignment.fasta with the new combined alignment
      3. If no new sequences: nothing to do (tree is already up to date)

    NOTE: Reference anchors are always present in the backbone alignment.
    On subsequent runs they are already in current_alignment.fasta so
    they don't need to be added again.
    """
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    current_aln = out / "current_alignment.fasta"
    current_tree = out / "current_tree.nwk"
    threads = int(threads) if threads and int(threads) > 0 else 4

    preflight_tree_tools(tree_method, require_mafft=True)

    resolved_anchor_file = get_anchor_file(anchor_file)
    orientation_summary = out / 'tree_orientation_summary.tsv'

    with open(user_fasta) as handle:
        user_count = sum(1 for line in handle if line.startswith(">"))
    logger.info("[TREE] Input: %d user sequences, previous_review=%s, force_rebuild=%s, outdir=%s",
                user_count, previous_review, force_rebuild, outdir)

    if not current_aln.exists() and not force_rebuild:
        seeded = _seed_backbone_from_previous_review(previous_review, out)
        if not seeded:
            # Performance Review organises persistent tree files under out/tree at the
            # end of a run. Reuse that backbone automatically when the same
            # output directory is supplied again.
            _seed_backbone_from_previous_review(str(out), out)

    if not current_aln.exists() or force_rebuild:
        logger.info("[TREE] %s — building from scratch",
                    "Force rebuild requested" if force_rebuild else "No backbone alignment found")
        try:
            if orientation_summary.exists():
                orientation_summary.unlink()
        except Exception:
            pass

        combined_fasta = build_combined_fasta(
            user_fasta=user_fasta,
            outdir=out,
            anchor_file=resolved_anchor_file,
            db=db,
            db_dataset=db_dataset,
            ref_fasta=ref_fasta,
            threads=threads,
            orientation_summary_path=orientation_summary,
            build_mode='rebuild' if force_rebuild else 'initial',
        )

        aligned_fasta = out / "current_alignment.fasta"
        if not _run_mafft_full(combined_fasta, aligned_fasta, threads):
            logger.error("[TREE] Alignment failed — cannot build tree")
            return

        fasta_for_tree, id_map = _make_unique_fasta(str(aligned_fasta), out)
        if not _build_tree(fasta_for_tree, current_tree, threads=threads, method=tree_method):
            return

        _finalise_tree(out, id_map)
        logger.info("[TREE] Initial tree built with %d sequences (%d anchors)",
                    sum(1 for _ in read_fasta(str(combined_fasta))),
                    len(load_anchor_sequences(resolved_anchor_file)))
        return

    if user_count == 0:
        if not current_tree.exists():
            logger.info("[TREE] No user sequences; building tree from existing alignment")
            fasta_for_tree, id_map = _make_unique_fasta(str(current_aln), out)
            if _build_tree(fasta_for_tree, current_tree, threads=threads, method=tree_method):
                _finalise_tree(out, id_map)
        else:
            logger.info("[TREE] No new sequences and tree exists — nothing to do")
        return

    new_records = _new_sequences_only(user_fasta, current_aln, db, out)
    logger.info("[TREE] _new_sequences_only: %d/%d user sequences are truly new (not in current alignment)",
                len(new_records), user_count)

    if not new_records:
        logger.info("[TREE] All user sequences already in alignment — nothing to add")
        if not current_tree.exists():
            fasta_for_tree, id_map = _make_unique_fasta(str(current_aln), out)
            if _build_tree(fasta_for_tree, current_tree, threads=threads, method=tree_method):
                _finalise_tree(out, id_map)
        return

    logger.info("[TREE] %d new sequences to add to existing alignment", len(new_records))

    new_fasta_path = out / "new_sequences.fasta"
    write_fasta(new_records, str(new_fasta_path))
    new_fasta_to_add, orientation_rows = _orient_tree_input_fasta(
        str(new_fasta_path),
        ref_fasta=ref_fasta,
        anchor_fasta=str(resolved_anchor_file) if resolved_anchor_file else None,
        outdir=out,
        label='new_sequences',
        build_mode='incremental',
        threads=threads,
    )
    if orientation_rows:
        _write_tree_orientation_summary(orientation_summary, orientation_rows, append=True)

    combined_aln = out / "combined_aln.fasta"
    mode = _choose_mafft_incremental_mode(new_records, current_aln)
    if mode == 'add':
        ok = _run_mafft_add(Path(new_fasta_to_add), current_aln, combined_aln, threads)
    else:
        ok = _run_mafft_addfragments(Path(new_fasta_to_add), current_aln, combined_aln, threads)
    if not ok:
        logger.error("[TREE] addfragments failed — tree not updated")
        return

    fasta_for_tree, id_map = _make_unique_fasta(str(combined_aln), out)
    if not _build_tree(fasta_for_tree, current_tree, threads=threads, method=tree_method):
        return

    _finalise_tree(out, id_map)

    # Promote combined alignment to current
    current_aln.write_text(combined_aln.read_text())
    logger.info("[TREE] Tree updated → %s", current_tree)
