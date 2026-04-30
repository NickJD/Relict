#!/usr/bin/env python3
"""Ad-hoc test runner for the pipeline changes (no pytest required)."""
import tempfile
from pathlib import Path

import sys
from pathlib import Path as _P
# ensure repository root is on sys.path so we can import pipeline modules
repo_root = _P(__file__).resolve().parents[1]
sys.path.insert(0, str(repo_root))

import pipeline.tree as tree_mod
import pipeline.classify as classify_mod
import pipeline.novelty as novelty_mod


def fake_run_tree(cmd):
    cmd = str(cmd)
    # extract outdir from redirection in command strings
    if "--auto" in cmd and ">" in cmd:
        parts = cmd.split(">")
        out = parts[1].strip()
        outp = Path(out)
        # create a simple alignment file
        outp.write_text(">ref1\nACTG\n")
    elif "--addfragments" in cmd and ">" in cmd:
        parts = cmd.split(">")
        out = parts[1].strip()
        Path(out).write_text(">ref1\nACTG\n>u1\nACTC\n")
    elif "FastTree" in cmd and ">" in cmd:
        parts = cmd.split(">")
        out = parts[1].strip()
        Path(out).write_text("(ref1,u1);")
    else:
        print("fake_run_tree got unexpected cmd:", cmd)


def fake_run_search(cmd):
    # ignore cmd contents; always write to out/matches.tsv
    # this function will be bound below where `out` is defined in outer scope
    out_dir = fake_run_search._outdir
    (out_dir / "matches.tsv").write_text("q1\tref1\t0.95\n")


def main():
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        out = td / "out"
        out.mkdir()

        ref = td / "ref.fasta"
        ref.write_text(">ref1\nACTGACTGACTG\n")

        user = td / "user.fasta"
        user.write_text(">u1\nACTGACTGACTC\n")

        # test tree
        tree_mod.run_cmd = fake_run_tree
        tree_mod.initialise_or_update_tree(str(ref), str(user), str(out))
        print("tree files:", list(out.iterdir()))

        # test classify/novelty
        fake_run_search._outdir = out
        classify_mod.run_cmd = fake_run_search
        novelty_mod.run_cmd = fake_run_search
        classify_mod.run_classification(str(user), str(out), ref_fasta=str(ref))
        novelty_mod.run_novelty(str(user), str(ref), str(out))
        print("taxonomy exists:", (out / "taxonomy.tsv").exists())
        print("novelty exists:", (out / "novelty.tsv").exists())


if __name__ == '__main__':
    main()




