import logging

from relict.utils.fasta import read_fasta, write_fasta

logger = logging.getLogger(__name__)


def run_qc(input_fasta, outdir, min_len=1200, max_n=5):
    """Filter sequences by length and allowed Ns.

    Writes filtered sequences to outdir/qc.fasta and a qc.stats file
    with counts and per-rejection-reason breakdown to help debug when
    outputs are unexpectedly empty.
    """
    output = f"{outdir}/qc.fasta"

    filtered = []
    total = 0
    rejected_too_short = []
    rejected_too_many_n = []

    for header, seq in read_fasta(input_fasta):
        total += 1
        if len(seq) < min_len:
            rejected_too_short.append(header)
            continue
        if seq.upper().count("N") > max_n:
            rejected_too_many_n.append(header)
            continue
        filtered.append((header, seq))

    write_fasta(filtered, output)

    # write qc stats for easier debugging
    stats_path = f"{outdir}/qc.stats"
    with open(stats_path, "w") as sf:
        sf.write(f"total_input\t{total}\n")
        sf.write(f"kept\t{len(filtered)}\n")
        sf.write(f"rejected_too_short\t{len(rejected_too_short)}\n")
        sf.write(f"rejected_too_many_n\t{len(rejected_too_many_n)}\n")
        sf.write(f"min_len\t{min_len}\n")
        sf.write(f"max_n\t{max_n}\n")
        sf.write("sample_ids_kept:\n")
        for h, _ in filtered[:10]:
            sf.write(f"\t{h}\n")
        if rejected_too_short:
            sf.write("sample_ids_too_short:\n")
            for h in rejected_too_short[:5]:
                sf.write(f"\t{h}\n")
        if rejected_too_many_n:
            sf.write("sample_ids_too_many_n:\n")
            for h in rejected_too_many_n[:5]:
                sf.write(f"\t{h}\n")

    logger.info(
        "[QC] %d input → %d kept, %d too short (<=%dbp), %d high-N (>%d). Stats: %s",
        total, len(filtered), len(rejected_too_short), min_len,
        len(rejected_too_many_n), max_n, stats_path,
    )
    if rejected_too_short:
        logger.debug("[QC] Rejected (too short): %s", ", ".join(rejected_too_short[:10]))
    if rejected_too_many_n:
        logger.debug("[QC] Rejected (high N): %s", ", ".join(rejected_too_many_n[:10]))

    return output

