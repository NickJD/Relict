#!/usr/bin/env python3
"""
build_anchor_fasta.py — Download and format reference anchor sequences for Relict.

Fetches one high-quality full-length 16S rRNA type strain sequence per major
bacterial and archaeal phylum from NCBI RefSeq, then writes them in the format
expected by Relict's tree builder:

    >RELICT_REF_<PhylumName> accession=<ACC> source=NCBI_RefSeq

Output: src/relict/data/reference_anchors.fasta

Usage
-----
    python scripts/build_anchor_fasta.py
    python scripts/build_anchor_fasta.py --out /path/to/custom_anchors.fasta
    python scripts/build_anchor_fasta.py --email your@email.com   # recommended by NCBI

Why this set?
-------------
The anchor set is designed to cover the major phyla found in gut and rumen
microbiome studies (emphasis on ruminants).  All sequences are NCBI RefSeq
type strain 16S entries (NR_ prefix) — the most stable, peer-reviewed, and
consistently updated 16S sequences available.

Anchors are topology scaffolds only — they never appear in your sequence
assessment, novelty metrics, or iTOL outputs.  The goal is to give FastTree
enough long-range signal to place the phylum-level branches correctly so that
your Hungate preload and run sequences always end up in the right neighbourhood.

Phyla covered
-------------
Bacteria
  Bacillota (Firmicutes)           — Ruminococcus, Clostridium, Bacillus,
                                     Lactobacillus, Streptococcus
  Bacteroidota                      — Bacteroides, Prevotella, Porphyromonas
  Pseudomonadota (Proteobacteria)   — Escherichia, Helicobacter, Campylobacter
  Actinomycetota (Actinobacteria)   — Bifidobacterium, Streptomyces
  Spirochaetota                     — Treponema, Borrelia
  Fibrobacterota                    — Fibrobacter succinogenes (key rumen fibre
                                       degrader — easily overlooked by primers)
  Fusobacteriota                    — Fusobacterium
  Planctomycetota                   — Planctomyces
  Verrucomicrobiota                 — Akkermansia muciniphila (gut mucosa)
  Chlorobiota (Chlorobi)            — Chlorobium
  Cyanobacteriota                   — Synechococcus
  Deinococcota                      — Deinococcus radiodurans
  Thermotogota                      — Thermotoga maritima
Archaea (outgroup / root)
  Methanobacteriota (Euryarchaeota) — Methanobrevibacter ruminantium (dominant
                                       rumen methanogen), Methanobacterium
  Thermoproteota (Crenarchaeota)    — Sulfolobus acidocaldarius (deep root)

NCBI accessions
---------------
These are RefSeq 16S rRNA gene entries (NR_*) for named type strains.
All are full-length or near-full-length (>1 400 bp).  Verify at:
    https://www.ncbi.nlm.nih.gov/nuccore/<ACCESSION>
"""

import argparse
import sys
import textwrap
import time
import urllib.error
import urllib.request
from pathlib import Path

# ---------------------------------------------------------------------------
# Anchor definitions
# Each entry: (relict_label, ncbi_accession, description)
# ---------------------------------------------------------------------------
ANCHORS = [
    # ── Bacillota (Firmicutes) ────────────────────────────────────────────
    ("Bacillota_Ruminococcus",   "NR_025930", "Ruminococcus albus 7 — rumen fibre degrader (Lachnospiraceae s.l.)"),
    ("Bacillota_Clostridium",    "NR_074545", "Clostridium butyricum ATCC 19398 — type genus Clostridia"),
    ("Bacillota_Bacillus",       "NR_027552", "Bacillus subtilis subsp. subtilis str. 168 — model Bacilli"),
    ("Bacillota_Lactobacillus",  "NR_075051", "Lactobacillus acidophilus ATCC 4356"),
    ("Bacillota_Streptococcus",  "NR_113594", "Streptococcus equinus NBRC 12553"),
    ("Bacillota_Butyrivibrio",   "NR_044858", "Butyrivibrio fibrisolvens D1 — rumen cellulolytic"),

    # ── Bacteroidota ──────────────────────────────────────────────────────
    ("Bacteroidota_Bacteroides",  "NR_041386", "Bacteroides fragilis ATCC 25285 — human/rumen gut"),
    ("Bacteroidota_Prevotella",   "NR_044825", "Prevotella bryantii B14 — dominant rumen Bacteroidota"),
    ("Bacteroidota_Porphyromonas","NR_040847", "Porphyromonas gingivalis ATCC 33277"),

    # ── Pseudomonadota (Proteobacteria) ───────────────────────────────────
    ("Pseudomonadota_Ecoli",       "NR_102804", "Escherichia coli K-12 — classic gamma-Proteobacteria reference"),
    ("Pseudomonadota_Helicobacter","NR_073694", "Helicobacter pylori 26695 — epsilon-Proteobacteria"),
    ("Pseudomonadota_Wolinella",   "NR_043184", "Wolinella succinogenes DSM 1740 — rumen epsilon-Proteobacteria"),

    # ── Actinomycetota (Actinobacteria) ───────────────────────────────────
    ("Actinomycetota_Bifidobacterium", "NR_040783", "Bifidobacterium longum NCC2705 — gut Actinobacteria"),
    ("Actinomycetota_Streptomyces",    "NR_043823", "Streptomyces griseus subsp. griseus ATCC 23345"),

    # ── Spirochaetota ─────────────────────────────────────────────────────
    ("Spirochaetota_Treponema", "NR_027243", "Treponema pallidum subsp. pallidum str. Nichols"),
    ("Spirochaetota_Borrelia",  "NR_025890", "Borrelia burgdorferi B31"),

    # ── Fibrobacterota ────────────────────────────────────────────────────
    # Critical: Fibrobacter succinogenes is the dominant rumen cellulose
    # degrader yet is missed by many universal primers.  Including it here
    # gives the tree a stable position to anchor these sequences when they
    # do appear.
    ("Fibrobacterota_Fibrobacter", "NR_041558", "Fibrobacter succinogenes subsp. succinogenes S85 — primary rumen fibre degrader"),

    # ── Fusobacteriota ────────────────────────────────────────────────────
    ("Fusobacteriota_Fusobacterium", "NR_026043", "Fusobacterium nucleatum subsp. nucleatum ATCC 25586"),

    # ── Planctomycetota ───────────────────────────────────────────────────
    ("Planctomycetota_Planctomyces", "NR_043399", "Planctomyces maris DSM 8797"),

    # ── Verrucomicrobiota ─────────────────────────────────────────────────
    ("Verrucomicrobiota_Akkermansia", "NR_042817", "Akkermansia muciniphila ATCC BAA-835 — gut mucosa"),

    # ── Thermotogota ──────────────────────────────────────────────────────
    ("Thermotogota_Thermotoga", "NR_043084", "Thermotoga maritima MSB8"),

    # ── Deinococcota ──────────────────────────────────────────────────────
    ("Deinococcota_Deinococcus", "NR_036779", "Deinococcus radiodurans R1"),

    # ── Cyanobacteriota ───────────────────────────────────────────────────
    ("Cyanobacteriota_Synechococcus", "NR_043317", "Synechococcus elongatus PCC 6301"),

    # ── Archaea — outgroup / root ─────────────────────────────────────────
    # Methanobrevibacter ruminantium is the predominant methanogen in cattle
    # rumen.  Its 16S is sufficiently divergent from Bacteria to function as
    # a meaningful outgroup for rumen bacterial trees.
    ("Archaea_Methanobrevibacter", "NR_044812", "Methanobrevibacter ruminantium M1 — dominant rumen methanogen"),
    ("Archaea_Methanobacterium",   "NR_028243", "Methanobacterium thermoautotrophicum str. Marburg"),
    ("Archaea_Sulfolobus",         "NR_043325", "Sulfolobus acidocaldarius DSM 639 — deep archaeal root (Crenarchaeota)"),
]

NCBI_EFETCH = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    "?db=nuccore&id={acc}&rettype=fasta&retmode=text"
)

# ---------------------------------------------------------------------------

def fetch_sequence(acc: str, email: str, retries: int = 3) -> str:
    """Fetch a single FASTA sequence from NCBI eutils.  Returns raw FASTA text."""
    url = NCBI_EFETCH.format(acc=acc)
    if email:
        url += f"&email={email}&tool=relict_build_anchors"
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=30) as resp:
                return resp.read().decode("utf-8")
        except urllib.error.HTTPError as e:
            if e.code == 429 or e.code == 503:
                wait = 5 * (attempt + 1)
                print(f"  Rate-limited ({e.code}); waiting {wait}s …", file=sys.stderr)
                time.sleep(wait)
            else:
                raise
        except Exception as e:
            if attempt < retries - 1:
                time.sleep(3)
            else:
                raise RuntimeError(f"Failed to fetch {acc}: {e}") from e
    raise RuntimeError(f"Exhausted retries for {acc}")


def parse_fasta_sequence(fasta_text: str) -> tuple[str, str]:
    """Extract (header, sequence) from a single-entry FASTA string."""
    lines = fasta_text.strip().splitlines()
    if not lines or not lines[0].startswith(">"):
        raise ValueError(f"Unexpected FASTA format: {fasta_text[:120]!r}")
    header = lines[0][1:].strip()
    seq = "".join(l.strip() for l in lines[1:] if not l.startswith(">"))
    return header, seq


def build_anchor_fasta(out_path: Path, email: str, skip_existing: bool = False):
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if skip_existing and out_path.exists():
        print(f"Anchor file already exists: {out_path}  (use --force to overwrite)")
        return

    records = []
    failed = []

    print(f"Fetching {len(ANCHORS)} anchor sequences from NCBI RefSeq …")
    for relict_label, acc, description in ANCHORS:
        print(f"  {acc:14s}  {relict_label} …", end=" ", flush=True)
        try:
            raw = fetch_sequence(acc, email)
            _orig_header, seq = parse_fasta_sequence(raw)
            if len(seq) < 800:
                print(f"WARN — sequence only {len(seq)} bp; keeping anyway")
            else:
                print(f"OK ({len(seq)} bp)")
            relict_header = f"RELICT_REF_{relict_label} accession={acc} source=NCBI_RefSeq desc={description!r}"
            records.append((relict_header, seq))
            time.sleep(0.4)  # stay well within NCBI rate limits
        except Exception as e:
            print(f"FAILED — {e}")
            failed.append((acc, relict_label, str(e)))

    # Write output
    with open(out_path, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            for chunk in textwrap.wrap(seq, 80):
                fh.write(chunk + "\n")

    print(f"\nWrote {len(records)} anchor sequences to {out_path}")
    if failed:
        print(f"\nFailed ({len(failed)}):")
        for acc, label, err in failed:
            print(f"  {acc}  {label}: {err}")
        print("\nRe-run the script to retry, or fetch these manually from")
        print("https://www.ncbi.nlm.nih.gov/nuccore/<ACCESSION>?report=fasta")
        print("and append them with headers like:")
        print("  >RELICT_REF_<PhylumName> accession=<ACC> source=NCBI_RefSeq")


# ---------------------------------------------------------------------------

def main():
    default_out = Path(__file__).resolve().parent.parent / "src" / "relict" / "data" / "reference_anchors.fasta"

    parser = argparse.ArgumentParser(
        description="Download NCBI RefSeq type-strain 16S sequences to use as Relict tree anchors."
    )
    parser.add_argument(
        "--out", default=str(default_out), metavar="PATH",
        help=f"Output FASTA path (default: {default_out})",
    )
    parser.add_argument(
        "--email", default="", metavar="EMAIL",
        help="Your email address — recommended by NCBI for E-utilities usage",
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Overwrite existing output file",
    )
    args = parser.parse_args()

    build_anchor_fasta(
        Path(args.out),
        email=args.email,
        skip_existing=not args.force,
    )


if __name__ == "__main__":
    main()


