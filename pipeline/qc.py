from utils.fasta import read_fasta, write_fasta


def run_qc(input_fasta, outdir, min_len=1200, max_n=5):
    """Filter sequences by length and allowed Ns.

    Writes filtered sequences to outdir/qc.fasta and a small qc.stats file
    with counts and a short sample of kept IDs to help debugging when
    outputs are unexpectedly empty.
    """
    output = f"{outdir}/qc.fasta"

    filtered = []
    total = 0

    for header, seq in read_fasta(input_fasta):
        total += 1
        if len(seq) < min_len:
            continue
        if seq.upper().count("N") > max_n:
            continue

        filtered.append((header, seq))

    write_fasta(filtered, output)

    # write qc stats for easier debugging
    stats_path = f"{outdir}/qc.stats"
    with open(stats_path, "w") as sf:
        sf.write(f"total_input\t{total}\n")
        sf.write(f"kept\t{len(filtered)}\n")
        sf.write("sample_ids:\n")
        for h, _ in filtered[:10]:
            sf.write(h + "\n")

    print(f"[QC] Read {total} sequences, kept {len(filtered)} (min_len={min_len}, max_n={max_n})")
    print(f"[QC] Wrote filtered sequences to {output} and stats to {stats_path}")

    return output
