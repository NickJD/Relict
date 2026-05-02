#!/usr/bin/env python3
"""Safe converter: produce a branch-ready copy of an existing iTOL colorstrip
file without changing the dataset template. iTOL does not accept
`DATASET_BRANCHCOLORS` — instead upload the colorstrip and enable branch
display in the web UI. This script simply writes a file with `_branch` in
the name but preserves the original, valid header.

Usage:
  convert_to_branch_itol.py itol_phylum_colors.itol [out.itol]
"""
from pathlib import Path
import sys


def convert_to_branch_filename(infile, outfile=None):
    p = Path(infile)
    if not p.exists():
        raise FileNotFoundError(f"Input file not found: {p}")
    if outfile is None:
        outfile = p.with_name(p.stem + '_branch' + p.suffix)
    # preserve the file content exactly (do NOT change the DATASET header)
    txt = p.read_text()
    Path(outfile).write_text(txt)
    print(f"Wrote branch-ready dataset (same template): {outfile}")


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: convert_to_branch_itol.py itol_phylum_colors.itol [out.itol]")
        sys.exit(1)
    try:
        convert_to_branch_filename(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(2)
