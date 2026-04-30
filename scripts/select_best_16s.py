#!/usr/bin/env python3
"""Select best 16S per genome and assign taxonomy from a greengenes-style DB.

This script expects a directory where each file is a FASTA containing
predicted 16S sequences for a single genome. For each genome file it:

- searches all sequences against the provided reference FASTA (gg2) with vsearch
- selects the sequence with the highest percent identity to the reference
- assigns taxonomy from a taxa TSV (FeatureID\tTaxon\tConfidence) if provided
- writes a combined FASTA of selected sequences and a summary TSV

The script supports gzipped reference FASTA and gzipped taxa TSVs.
"""
import argparse
from pathlib import Path
# os not required
import gzip
import re
from utils.subprocess import run_cmd
from utils.fasta import read_fasta, write_fasta


def _norm_id(x: str) -> str:
    x = x.split()[0]
    if '|' in x:
        x = x.split('|')[-1]
    x = re.sub(r'_[0-9]+$', '', x)
    return x


def load_taxa(taxa_tsv):
    taxa_map = {}
    if not taxa_tsv:
        return taxa_map
    open_fn = gzip.open if str(taxa_tsv).endswith('.gz') else open
    with open_fn(taxa_tsv, 'rt') as f:
        first = f.readline()
        if 'Feature' in first or 'Taxon' in first:
            # header present, continue
            pass
        else:
            parts = first.strip().split('\t')
            if len(parts) >= 2:
                conf = float(parts[2]) if len(parts) > 2 else None
                taxa_map[parts[0]] = (parts[1], conf)
                taxa_map[_norm_id(parts[0])] = (parts[1], conf)
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            fid = parts[0]
            tax = parts[1]
            conf = float(parts[2]) if len(parts) > 2 else None
            taxa_map[fid] = (tax, conf)
            taxa_map[_norm_id(fid)] = (tax, conf)
    return taxa_map


def ensure_ref_uncompressed(ref, outdir):
    # if gz, write uncompressed copy into outdir
    if str(ref).endswith('.gz'):
        out_unc = Path(outdir) / 'ref_uncompressed.fasta'
        if not out_unc.exists():
            print(f"[REF] Decompressing {ref} -> {out_unc}")
            records = list(read_fasta(ref))
            write_fasta(records, str(out_unc))
        return str(out_unc)
    return str(ref)


def select_best_for_genome(genome_fasta, ref_fasta, outdir, verbose=False):
    # run vsearch best-hit per query
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    # use filename without .gz and without fasta extension to name matches file
    gf_name = Path(genome_fasta).name
    if gf_name.endswith('.gz'):
        gf_base = gf_name[:-3]
    else:
        gf_base = gf_name
    matches = outdir / f"matches_{Path(gf_base).stem}.tsv"
    cmd = f"vsearch --usearch_global {genome_fasta} --db {ref_fasta} --id 0.0 --blast6out {matches} --maxaccepts 1 --maxhits 1"
    run_cmd(cmd)

    best_hits = {}
    if matches.exists():
        linecount = 0
        with open(matches) as f:
            for line in f:
                linecount += 1
                parts = line.strip().split('\t')
                if len(parts) < 3:
                    continue
                qid = parts[0]
                sid = parts[1]
                try:
                    ident = float(parts[2])
                except Exception:
                    ident = 0.0
                best_hits[qid] = (sid, ident)
        if verbose:
            print(f"[VSEARCH] matches found: {linecount} lines -> {len(best_hits)} unique queries ({matches})")
            # print first few lines for debugging
            try:
                with open(matches) as f:
                    for i, L in enumerate(f):
                        if i >= 5:
                            break
                        print("[VSEARCH] sample:", L.strip())
            except Exception:
                pass
    else:
        if verbose:
            print(f"[VSEARCH] matches file not found: {matches}")

    # choose best sequence: highest identity; if no matches, fallback to longest
    best_qid = None
    best_ident = -1.0
    best_sid = None
    # build seqs map with multiple keys to improve matching with vsearch qids
    seqs = {}
    for h, s in read_fasta(genome_fasta):
        seqs[h] = s
        # also map first token (vsearch reports qid up to first whitespace)
        first = h.split()[0]
        if first not in seqs:
            seqs[first] = s
        # also map normalized id (strip trailing _N and pipe prefixes)
        try:
            norm = _norm_id(h)
            if norm not in seqs:
                seqs[norm] = s
        except Exception:
            pass
    if verbose:
        # count unique original headers (approx): count headers returned by read_fasta
        unique_headers = set()
        for h, _ in read_fasta(genome_fasta):
            unique_headers.add(h)
        print(f"[READ] {genome_fasta} -> {len(unique_headers)} sequences read (mapped keys: {len(seqs)})")

    for qid, seq in seqs.items():
        if qid in best_hits:
            sid, ident = best_hits[qid]
            if ident > best_ident:
                best_ident = ident
                best_qid = qid
                best_sid = sid

    if best_qid is None:
        # no hits at all, pick longest sequence
        best_qid = max(seqs.keys(), key=lambda k: len(seqs[k])) if seqs else None
        best_sid = None
        best_ident = 0.0

    best_seq = seqs.get(best_qid) if best_qid else None
    if verbose:
        if best_seq:
            seq_len = len(best_seq) if isinstance(best_seq, (str, bytes, list, tuple)) else 0
            print(f"[SELECT] chosen {best_qid} (len={seq_len}) ident={best_ident} sid={best_sid}")
        else:
            print(f"[SELECT] no sequence chosen for {genome_fasta}; seqs_present={bool(seqs)}")
    return best_qid, best_seq, best_sid, best_ident


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--input-dir', required=True)
    p.add_argument('--ref', required=True, help='Reference fasta (gg2) - may be gzipped')
    p.add_argument('--taxa', required=False, help='Taxa TSV (FeatureID\tTaxon\tConfidence) - may be gzipped')
    p.add_argument('--out', required=True)
    p.add_argument('--min-identity', dest='min_id', type=float, default=0.0, help='Minimum percent identity required to accept a selected 16S (default 0.0)')
    p.add_argument('--min-len', dest='min_len', type=int, default=100, help='Minimum sequence length to consider (default: 100)')
    p.add_argument('--max-n', dest='max_n', type=int, default=3, help='Maximum number of Ns allowed in sequence (default: 3)')
    p.add_argument('--max-per-genome', dest='max_per_genome', type=int, default=1, help='Maximum sequences to select per genome (default: 1)')
    p.add_argument('--verbose', action='store_true', help='Enable verbose logging')
    args = p.parse_args()

    inp = Path(args.input_dir)
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)

    ref_to_use = ensure_ref_uncompressed(args.ref, outdir)
    taxa_map = load_taxa(args.taxa) if args.taxa else {}

    # We'll do a single-pass search: concatenate filtered sequences from all genomes
    # into one query FASTA, run vsearch once, then choose the best sequence per genome.
    selected = []
    summary = []

    valid_exts = ('.fa', '.fasta', '.fna', '.ffn')
    combined_queries = outdir / 'combined_queries.fasta'
    matches = outdir / 'matches_combined.tsv'

    # collect filtered records and mapping genome->qids
    records = []
    genome_qids = {}
    seq_map = {}

    for fasta in sorted(inp.iterdir()):
        if not fasta.is_file():
            continue
        name = fasta.name
        lower = name.lower()
        ok = False
        for ext in valid_exts:
            if lower.endswith(ext) or lower.endswith(ext + '.gz'):
                ok = True
                break
        if not ok:
            if getattr(args, 'verbose', False):
                print(f"[SKIP] Ignoring non-FASTA file: {fasta}")
            continue
        genome_id = fasta.stem
        if getattr(args, 'verbose', False):
            print(f"[COLLECT] Scanning {genome_id} ({fasta})")
        for h, s in read_fasta(str(fasta)):
            # apply simple length/N filters
            if args.min_len and len(s) < args.min_len:
                continue
            if args.max_n is not None and s.count('N') > args.max_n:
                continue
            # make a query id that encodes genome and original id; vsearch reports qid up to first whitespace
            orig_token = h.split()[0]
            qid = f"{genome_id}|{orig_token}"
            records.append((qid, s))
            genome_qids.setdefault(genome_id, []).append(qid)
            seq_map[qid] = s

    # write combined queries
    best_hits = {}
    if records:
        write_fasta(records, str(combined_queries))
        # run vsearch once on combined queries
        cmd = f"vsearch --usearch_global {combined_queries} --db {ref_to_use} --id 0.0 --blast6out {matches} --maxaccepts 1 --maxhits 1"
        run_cmd(cmd)

        # parse matches into best_hits per qid
        best_hits = {}
        if Path(matches).exists():
            with open(matches) as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) < 3:
                        continue
                    qid = parts[0]
                    sid = parts[1]
                    try:
                        ident = float(parts[2])
                    except Exception:
                        ident = 0.0
                    best_hits[qid] = (sid, ident)
        else:
            if getattr(args, 'verbose', False):
                print(f"[VSEARCH] combined matches file not found: {matches}")

    # select top-N per genome from parsed results (default N=1)
    max_per = max(1, int(getattr(args, 'max_per_genome', 1)))
    for fasta in sorted(inp.iterdir()):
        if not fasta.is_file():
            continue
        name = fasta.name
        lower = name.lower()
        ok = False
        for ext in valid_exts:
            if lower.endswith(ext) or lower.endswith(ext + '.gz'):
                ok = True
                break
        if not ok:
            continue
        genome_id = fasta.stem
        qids = genome_qids.get(genome_id, [])

        # build candidate list: (qid, identity, sid)
        candidates = []
        for qid in qids:
            if qid in best_hits:
                sid, ident = best_hits[qid]
                if ident is None:
                    ident = 0.0
                candidates.append((qid, float(ident), sid))
            else:
                # no hit -> treat as identity 0.0 (still eligible if min_id==0)
                candidates.append((qid, 0.0, None))

        if not candidates:
            if getattr(args, 'verbose', False):
                print(f"[SELECT] no candidates for {genome_id}")
            continue

        # sort candidates by identity desc, then by sequence length desc as tiebreaker
        candidates.sort(key=lambda t: (t[1], len(seq_map.get(t[0], ''))), reverse=True)

        picks = candidates[:max_per]
        for qid, ident, sid in picks:
            # apply min identity threshold
            if args.min_id > 0.0 and ident < args.min_id:
                if getattr(args, 'verbose', False):
                    print(f"[INFO] Skipping {qid} for {genome_id}: ident {ident} < min_id {args.min_id}")
                continue
            seq = seq_map.get(qid)
            if not seq:
                continue
            # selected id should preserve original header token
            orig = qid.split('|', 1)[1] if '|' in qid else qid
            selected_id = f"{genome_id}|{orig}"
            selected.append((selected_id, seq))

            tax = 'NA'
            conf = None
            if sid and taxa_map:
                entry = taxa_map.get(sid)
                if entry is None:
                    entry = taxa_map.get(_norm_id(str(sid)))
                if entry:
                    tax, conf = entry

            summary.append((genome_id, orig, sid or 'NA', f"{ident}", tax or 'NA', str(conf) if conf is not None else 'NA'))

    # write combined selected fasta and summary
    combined_fasta = outdir / 'selected_16s.fasta'
    if selected:
        write_fasta(selected, str(combined_fasta))
    else:
        # write an empty file and warn
        open(combined_fasta, 'w').close()
        print('[WARN] No sequences were selected; combined FASTA is empty')

    summary_path = outdir / 'selected_16s_summary.tsv'
    with open(summary_path, 'w') as s:
        s.write('genome_id\tbest_seq_id\tbest_hit\tidentity\ttaxon\tconfidence\n')
        for row in summary:
            s.write('\t'.join(row) + '\n')

    print(f"[DONE] Wrote {len(selected)} selected sequences to {combined_fasta}")
    print(f"[DONE] Summary written to {summary_path}")


if __name__ == '__main__':
    main()

