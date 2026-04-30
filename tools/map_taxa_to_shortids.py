#!/usr/bin/env python3
"""Map combined taxonomy IDs to tree short IDs and write corrected iTOL files.

Usage:
  python tools/map_taxa_to_shortids.py \
    --combined /path/to/combined_taxonomy.tsv \
    --user-map /path/to/user_id_map.tsv \
    --db /path/to/Current_DB \
    --preload-itol /path/to/Current_Preload/itol_dataset_preload.itol \
    --outdir /path/to/Current_New

This script attempts several heuristics to map taxonomy rows (which often refer
to original/reference ids) to the short IDs used as tree tip labels. It then
writes a corrected `itol_phylum_colors.itol` (and `itol_membership.itol`) into
the provided outdir so colours will be applied to the tree tips.
"""
import argparse
from pathlib import Path
import re
import sqlite3
import csv
import sys


def canonicalize(x: str) -> str:
    if x is None:
        return ''
    x = str(x).split()[0]
    if '|' in x:
        x = x.split('|')[-1]
    x = re.sub(r'_[0-9]+$', '', x)
    return x


def read_user_map(path: Path):
    # user_id_map: short_id\toriginal_header
    orig_to_short = {}
    short_set = set()
    if not path.exists():
        return orig_to_short, short_set
    # use csv to robustly parse TSV and allow additional columns
    with open(path, newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        # optionally skip header if it looks like one
        try:
            first = next(reader)
        except StopIteration:
            return orig_to_short, short_set
        if len(first) >= 2 and (first[0].lower() in ('short', 'short_id') or first[1].lower().startswith('orig')):
            pass  # header consumed
        else:
            # treat first row as data
            parts = first
            if len(parts) >= 2:
                short, orig = parts[0].strip(), parts[1].strip()
                if orig:
                    orig_to_short[orig] = short
                    orig_to_short[canonicalize(orig)] = short
                    short_set.add(short)
        for parts in reader:
            if not parts:
                continue
            if len(parts) >= 2:
                short, orig = parts[0].strip(), parts[1].strip()
                if orig:
                    orig_to_short[orig] = short
                    orig_to_short[canonicalize(orig)] = short
                    short_set.add(short)
    return orig_to_short, short_set


def read_preload_itol(itol_path: Path):
    pre_ids = set()
    pre_color = None
    if not itol_path.exists():
        return pre_ids, pre_color
    with open(itol_path) as pf:
        for ln in pf:
            ln = ln.strip()
            if not ln:
                continue
            if ln.startswith('COLOR,') and pre_color is None:
                try:
                    pre_color = ln.split(',', 1)[1]
                except Exception:
                    pre_color = None
            # data lines are id,color (avoid header lines)
            if ',' in ln and not any(ln.startswith(x) for x in ('DATASET_COLORSTRIP','SEPARATOR','DATASET_LABEL','COLOR','MARGIN','SHOW_INTERNAL','DATA')):
                a, b = ln.split(',', 1)
                pre_ids.add(a)
    return pre_ids, pre_color


def query_db_aliases(db_path: Path):
    # build mapping original_header -> canonical/short id stored in seq_aliases
    mapping = {}
    try:
        conn = sqlite3.connect(str(db_path))
        cur = conn.cursor()
        cur.execute('SELECT canonical_id, original_header FROM seq_aliases')
        for cid, orig in cur.fetchall():
            if orig:
                mapping[orig] = cid
                mapping[canonicalize(orig)] = cid
        conn.close()
    except Exception:
        pass
    return mapping


def parse_combined(path: Path):
    taxa = {}
    # use csv to robustly parse TSV; allow taxonomy to contain tabs by joining remaining columns
    with open(path, newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        try:
            first = next(reader)
        except StopIteration:
            return taxa
        # if header-like (first cell equals 'qid' or 'id'), skip it
        if len(first) >= 1 and first[0].lower() in ('qid', 'id', 'query'):
            pass
        else:
            # treat first as data
            if len(first) >= 1:
                qid = first[0].strip()
                tax = '\t'.join(first[1:]).strip() if len(first) > 1 else 'NA'
                taxa[qid] = tax
        for parts in reader:
            if not parts:
                continue
            qid = parts[0].strip()
            tax = '\t'.join(parts[1:]).strip() if len(parts) > 1 else 'NA'
            taxa[qid] = tax
    return taxa


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--combined', required=True)
    p.add_argument('--user-map', required=False)
    p.add_argument('--preload-itol', required=False)
    p.add_argument('--db', required=False)
    p.add_argument('--outdir', required=True)
    args = p.parse_args()

    combined = Path(args.combined)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    taxa = parse_combined(combined)

    orig_to_short = {}
    short_set = set()
    if args.user_map:
        om, ss = read_user_map(Path(args.user_map))
        orig_to_short.update(om)
        short_set.update(ss)

    db_aliases = {}
    if args.db:
        db_aliases = query_db_aliases(Path(args.db))
        # also try mapping canonical ids present in seq_aliases
        for k, v in list(db_aliases.items()):
            if v not in short_set:
                short_set.add(v)

    pre_ids = set()
    pre_color = None
    if args.preload_itol:
        pre_ids, pre_color = read_preload_itol(Path(args.preload_itol))

    # avoid importing the package to keep this script runnable from repo root
    # reimplement minimal functions from phylo16s.pipeline.itol used below
    def _hash_to_hue(s: str) -> float:
        import hashlib
        h = hashlib.md5(s.encode('utf-8')).hexdigest()
        val = int(h[:8], 16)
        return (val % 360) / 360.0

    def _hsv_to_hex(h, s=0.65, v=0.95):
        import colorsys
        r, g, b = colorsys.hsv_to_rgb(h, s, v)
        return '#{:02x}{:02x}{:02x}'.format(int(r*255), int(g*255), int(b*255))

    def _name_to_color_by_rank(name: str, rank: str = None) -> str:
        h = _hash_to_hue(name)
        rank_offsets = {'phylum': 0.0, 'family': 0.33, 'genus': 0.66}
        offset = rank_offsets.get(rank, 0.0)
        h = (h + offset) % 1.0
        sv_map = {'phylum': (0.85, 0.97), 'family': (0.65, 0.90), 'genus': (0.45, 0.80)}
        s, v = sv_map.get(rank, (0.65, 0.90))
        return _hsv_to_hex(h, s=s, v=v)

    def parse_taxon_string(taxon: str):
        parts = [p.strip() for p in str(taxon).split(';')]
        res = {}
        for p in parts:
            if '__' in p:
                k, v = p.split('__', 1)
                k = k.strip()
                v = v.strip()
                if v:
                    res[k] = v
        return res
    # use these local functions as itol_mod
    class itol_mod:
        @staticmethod
        def _name_to_color_by_rank(name, rank=None):
            return _name_to_color_by_rank(name, rank)
        @staticmethod
        def parse_taxon_string(t):
            return parse_taxon_string(t)
    # heuristics mapping
    mapped = {}
    unmapped = []
    for qid, tax in taxa.items():
        tgt = None
        # 1) already short id (three uppercase letters + two or more digits)
        if re.match(r'^[A-Z]{3}\d{2,}$', qid):
            tgt = qid
        # 2) direct match in user-provided maps
        if tgt is None and qid in orig_to_short:
            tgt = orig_to_short[qid]
        # 3) exact match in db alias original headers
        if tgt is None and qid in db_aliases:
            tgt = db_aliases[qid]
        # 4) canonicalize and test
        if tgt is None:
            c = canonicalize(qid)
            if c in orig_to_short:
                tgt = orig_to_short[c]
            elif c in db_aliases:
                tgt = db_aliases[c]
        if tgt:
            mapped[qid] = (tgt, tax)
        else:
            unmapped.append(qid)

    print(f"Total combined ids: {len(taxa)}; mapped to short ids: {len(mapped)}; unmapped: {len(unmapped)}")
    if unmapped:
        print('Examples of unmapped ids:', unmapped[:20])

    # Build phylum color dataset for mapped short ids
    ph_pairs = []
    for orig_id, (short, tax) in mapped.items():
        parsed = itol_mod.parse_taxon_string(tax)
        ph = parsed.get('p', 'unknown')
        color = itol_mod._name_to_color_by_rank(ph, 'phylum')
        ph_pairs.append((short, color))

    # write phylum colorstrip
    itol_path = outdir / 'itol_phylum_colors_fixed.itol'
    with open(itol_path, 'w') as f:
        f.write('DATASET_COLORSTRIP\n')
        f.write('SEPARATOR COMMA\n')
        f.write('DATASET_LABEL,Phylum colors\n')
        f.write('COLOR,#AAAAAA\n')
        f.write('MARGIN,5\n')
        f.write('SHOW_INTERNAL,0\n')
        f.write('DATA\n')
        for sid, col in ph_pairs:
            f.write(f"{sid},{col}\n")

    # membership strip: preload vs run vs other
    run_color = '#ff3333'
    other_color = '#cccccc'
    if pre_color is None:
        pre_color = '#1f78b4'

    membership_path = outdir / 'itol_membership_fixed.itol'
    with open(membership_path, 'w') as f:
        f.write('DATASET_COLORSTRIP\n')
        f.write('SEPARATOR COMMA\n')
        f.write('DATASET_LABEL,Dataset membership\n')
        f.write('COLOR,#AAAAAA\n')
        f.write('MARGIN,5\n')
        f.write('SHOW_INTERNAL,0\n')
        f.write('DATA\n')
        for orig_id, (short, tax) in mapped.items():
            if short in pre_ids:
                col = pre_color
            elif short in short_set:
                col = run_color
            else:
                col = other_color
            f.write(f"{short},{col}\n")

    print('Wrote:', itol_path, membership_path)
    print('If these files are correct, upload them to iTOL with your tree or replace the existing itol files in your run outdir and re-upload to iTOL.')


if __name__ == '__main__':
    main()





