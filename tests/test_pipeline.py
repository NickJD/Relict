import os
from pathlib import Path

import phylo16s.pipeline.tree as tree_mod
import phylo16s.pipeline.classify as classify_mod
import phylo16s.pipeline.novelty as novelty_mod


def test_tree_addition(tmp_path, monkeypatch):
    outdir = tmp_path / "out"
    outdir.mkdir()

    ref = tmp_path / "ref.fasta"
    ref.write_text(">ref1\nACTGACTGACTG\n>ref2\nACTGACTGACTA\n")

    user = tmp_path / "user.fasta"
    user.write_text(">u1\nACTGACTGACTC\n")

    # monkeypatch run_cmd to simulate external tools
    def fake_run(cmd):
        cmd = str(cmd)
        if "--auto" in cmd:
            # simulate mafft --auto: write ref_aln.fasta
            (outdir / "ref_aln.fasta").write_text(ref.read_text())
        elif "--addfragments" in cmd:
            # simulate mafft --addfragments: combine ref_aln and user
            combined = ref.read_text() + user.read_text()
            (outdir / "combined_aln.fasta").write_text(combined)
        elif "FastTree" in cmd:
            # simulate tree creation
            (outdir / "current_tree.nwk").write_text("(ref1,u1);")
        else:
            raise RuntimeError(f"Unexpected cmd in fake_run: {cmd}")

    monkeypatch.setattr(tree_mod, "run_cmd", fake_run)

    # Run the function
    tree_mod.initialise_or_update_tree(str(ref), str(user), str(outdir))

    # check outputs
    assert (outdir / "current_alignment.fasta").exists()
    assert (outdir / "current_tree.nwk").exists()


def test_classify_and_novelty(tmp_path, monkeypatch):
    outdir = tmp_path / "out"
    outdir.mkdir()

    ref = tmp_path / "ref.fasta"
    ref.write_text(">ref1\nACTGACTGACTG\n")

    user = tmp_path / "user.fasta"
    user.write_text(">q1\nACTGACTGACTC\n")

    # fake run_cmd to create matches.tsv
    def fake_run(cmd):
        cmd = str(cmd)
        if "--usearch_global" in cmd:
            # write a match line: qid\trefid\t0.95
            (outdir / "matches.tsv").write_text("q1\tref1\t0.95\n")
        else:
            raise RuntimeError("Unexpected cmd")

    monkeypatch.setattr(classify_mod, "run_cmd", fake_run)
    monkeypatch.setattr(novelty_mod, "run_cmd", fake_run)

    tax = classify_mod.run_classification(str(user), str(outdir), ref_fasta=str(ref))
    assert (outdir / "taxonomy.tsv").exists()
    with open(tax) as f:
        assert "q1" in f.read()

    nov = novelty_mod.run_novelty(str(user), str(ref), str(outdir), id_threshold=0.9)
    assert (outdir / "novelty.tsv").exists()
    with open(nov) as f:
        data = f.read()
        assert "q1" in data

