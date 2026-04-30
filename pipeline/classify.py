
"""Classification helpers for phylo16s.

This module wraps a simple vsearch-based best-hit classifier and converts
the blast6-style output into a small taxonomy table. It logs useful
diagnostic information to help debugging classification and mapping issues.
"""

from utils.subprocess import run_cmd
import logging
import os
import gzip
import re
from utils.fasta import read_fasta, write_fasta

logger = logging.getLogger(__name__)


def run_classification(input_fasta, outdir, ref_fasta=None, taxa_tsv=None, threads=None):
    """Classify sequences by searching the provided reference fasta.

    This uses vsearch --usearch_global to find the best hit in the reference
    fasta and writes a simple taxonomy table with columns:
    ID\tBestHit\tIdentity\tTaxon\tConfidence
    """
    output = f"{outdir}/taxonomy.tsv"

    if ref_fasta is None:
        raise ValueError("ref_fasta must be provided to run_classification")

    # if ref_fasta is gzipped, create an uncompressed copy in outdir
    ref_to_use = ref_fasta
    if str(ref_fasta).endswith('.gz'):
        ref_unc = os.path.join(outdir, 'ref_uncompressed.fasta')
        if not os.path.exists(ref_unc):
            logger.info("Uncompressing reference fasta %s -> %s", ref_fasta, ref_unc)
            records = list(read_fasta(ref_fasta))
            write_fasta(records, ref_unc)
        ref_to_use = ref_unc

    # run similarity search; keep a single best hit per query
    thread_flag = f" --threads {int(threads)}" if threads and int(threads) > 0 else ""
    cmd = (
        f"vsearch --usearch_global {input_fasta} --db {ref_to_use} --id 0.8 "
        f"--blast6out {outdir}/matches.tsv --maxaccepts 1 --maxhits 1{thread_flag}"
    )

    logger.debug("Running vsearch command: %s", cmd)
    logger.info("[CLASSIFY] Running vsearch for classification (input=%s, db=%s)", input_fasta, ref_to_use)
    run_cmd(cmd)
    logger.info("[CLASSIFY] vsearch finished; matches written to %s/matches.tsv", outdir)

    # helper to normalise ids (available even if no taxa_tsv provided)
    def _norm_id(x: str) -> str:
        # take first token (drop whitespace-separated metadata)
        x = x.split()[0]
        # if pipes present (e.g. ...|...|...|ID:coords), prefer the last field
        if '|' in x:
            x = x.split('|')[-1]
        # drop everything after a '#' (common in names like name#_NODE_...)
        if '#' in x:
            x = x.split('#')[0]
        # drop coordinate suffixes like :123-456(+) or :0-123(-)
        if ':' in x:
            x = x.split(':')[0]
        # remove trailing parentheses or stray punctuation
        x = x.strip().strip('()[]{}')
        # collapse repeated non-alphanumeric runs to single underscore for stability
        x = re.sub(r'[^0-9A-Za-z]+', '_', x).strip('_')
        return x

    def _canon_id(x: str) -> str:
        """More aggressive canonicalization: remove node fragments (#...), coords (:..), parentheses, and collapse non-alnum to underscore."""
        if x is None:
            return x
        y = str(x).split()[0]
        # remove everything after '#' (node annotations)
        if '#' in y:
            y = y.split('#')[0]
        # remove coordinate suffixes like :123-456(+) or :0-123(-)
        y = re.sub(r':\d+-\d+\(.*\)$', '', y)
        y = re.sub(r':\d+-\d+$', '', y)
        # prefer last pipe field
        if '|' in y:
            y = y.split('|')[-1]
        # collapse non-alnum to underscore
        y = re.sub(r'[^0-9A-Za-z]+', '_', y).strip('_')
        return y

    # optionally load taxa mapping (FeatureID -> (Taxon, Confidence))
    taxa_map = {}
    if taxa_tsv:
        open_fn = gzip.open if str(taxa_tsv).endswith('.gz') else open
        try:
            with open_fn(taxa_tsv, 'rt') as t:
                # attempt to detect and skip header
                first = t.readline()
                if not ("Feature" in first or "Taxon" in first):
                    parts = first.strip().split("\t")
                    if len(parts) >= 2:
                        try:
                            conf = float(parts[2]) if len(parts) > 2 else None
                        except Exception:
                            conf = None
                        fid = parts[0]
                        tax = parts[1]
                        taxa_map[fid] = (tax, conf)
                        taxa_map[_norm_id(fid)] = (tax, conf)
                        taxa_map[_canon_id(fid)] = (tax, conf)
                for line in t:
                    parts = line.strip().split("\t")
                    if len(parts) < 2:
                        continue
                    fid = parts[0]
                    tax = parts[1]
                    try:
                        conf = float(parts[2]) if len(parts) > 2 else None
                    except Exception:
                        conf = None
                    taxa_map[fid] = (tax, conf)
                    taxa_map[_norm_id(fid)] = (tax, conf)
                    taxa_map[_canon_id(fid)] = (tax, conf)
            logger.info("Loaded %d taxa mappings from %s", len(taxa_map)//2, taxa_tsv)
        except Exception as e:
            logger.warning("Failed to load taxa mapping %s: %s", taxa_tsv, e)

    # parse matches.tsv and write a small taxonomy table
    matches_path = os.path.join(outdir, 'matches.tsv')
    written = 0
    with open(output, 'w') as out:
        out.write("ID\tBestHit\tIdentity\tTaxon\tConfidence\n")
        try:
            with open(matches_path) as m:
                for line in m:
                    if not line.strip():
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) < 3:
                        continue
                    qid = parts[0]
                    sid = parts[1]
                    identity = parts[2]
                    tax = None
                    conf = None
                    if sid in taxa_map:
                        tax, conf = taxa_map[sid]
                    else:
                        # try normalized and canonical forms
                        tax, conf = taxa_map.get(_norm_id(sid), (None, None))
                        if tax is None:
                            tax, conf = taxa_map.get(_canon_id(sid), (None, None))
                        # fallback: try matching on tokens from the sid (split by common separators)
                        # Many reference IDs can include prefixes/suffixes (e.g. RS-GCF-...-NZ-...);
                        # attempt to find any token that maps to a taxa entry.
                        if tax is None:
                            import re as _re
                            tokens = [t for t in _re.split(r"[\|_:\-.]+", sid) if t]
                            for tok in tokens:
                                if tok in taxa_map:
                                    tax, conf = taxa_map[tok]
                                    break
                                nt = _norm_id(tok)
                                if nt in taxa_map:
                                    tax, conf = taxa_map[nt]
                                    break
                                ct = _canon_id(tok)
                                if ct in taxa_map:
                                    tax, conf = taxa_map[ct]
                                    break
                    out.write(f"{qid}\t{sid}\t{identity}\t{tax if tax is not None else 'NA'}\t{conf if conf is not None else 'NA'}\n")
                    written += 1
        except FileNotFoundError:
            logger.warning("matches file %s not found; taxonomy output will contain only header", matches_path)
        except Exception as e:
            logger.warning("Error while parsing matches file %s: %s", matches_path, e)

    logger.info("Wrote %d taxonomy entries to %s", written, output)
    return output

