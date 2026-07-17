import logging

from branchmanager.utils.fasta import read_fasta, write_fasta

logger = logging.getLogger(__name__)


def run_qc(input_fasta, outdir, min_len=800, max_n_percent=5.0):
    """Filter sequences by length and percentage of ambiguous bases.

    Writes filtered sequences to outdir/qc.fasta, a qc.stats summary, and
    qc_rejections.tsv with length/N counts and per-sequence rejection reasons.
    """
    max_n_percent = float(max_n_percent)
    if not 0.0 <= max_n_percent <= 100.0:
        raise ValueError("max_n_percent must be between 0 and 100")
    output = f"{outdir}/qc.fasta"

    filtered = []
    total = 0
    rejected_too_short = []
    rejected_too_many_n = []
    rejected_details = []

    for header, seq in read_fasta(input_fasta):
        total += 1
        seq_len = len(seq)
        n_count = seq.upper().count("N")
        n_percent = (100.0 * n_count / seq_len) if seq_len else 0.0
        reasons = []
        if seq_len < min_len:
            reasons.append("too_short")
        if n_percent > max_n_percent:
            reasons.append("too_many_n")

        if reasons:
            if "too_short" in reasons:
                rejected_too_short.append(header)
            if "too_many_n" in reasons:
                rejected_too_many_n.append(header)
            rejected_details.append((header, seq_len, n_count, n_percent, ";".join(reasons)))
            continue
        filtered.append((header, seq))

    write_fasta(filtered, output)

    rejection_path = f"{outdir}/qc_rejections.tsv"
    with open(rejection_path, "w") as rf:
        rf.write("ID\tLength\tNCount\tNPercent\tReasons\tMinLength\tMaxNPercent\n")
        for header, seq_len, n_count, n_percent, reasons in rejected_details:
            rf.write(
                f"{header}\t{seq_len}\t{n_count}\t{n_percent:.6f}\t{reasons}\t"
                f"{min_len}\t{max_n_percent:g}\n"
            )

    # write qc stats for easier debugging
    stats_path = f"{outdir}/qc.stats"
    with open(stats_path, "w") as sf:
        sf.write(f"total_input\t{total}\n")
        sf.write(f"kept\t{len(filtered)}\n")
        sf.write(f"rejected_total\t{len(rejected_details)}\n")
        sf.write(f"rejected_too_short\t{len(rejected_too_short)}\n")
        sf.write(f"rejected_too_many_n\t{len(rejected_too_many_n)}\n")
        sf.write(f"min_len\t{min_len}\n")
        sf.write(f"max_n_percent\t{max_n_percent:g}\n")
        sf.write(f"rejection_details\t{rejection_path}\n")
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
        "[QC] %d input → %d kept, %d rejected; %d too short (<%dbp), %d high-N (>%.2f%%). Stats: %s; details: %s",
        total, len(filtered), len(rejected_details), len(rejected_too_short), min_len,
        len(rejected_too_many_n), max_n_percent, stats_path, rejection_path,
    )
    if rejected_too_short:
        logger.debug("[QC] Rejected (too short): %s", ", ".join(rejected_too_short[:10]))
    if rejected_too_many_n:
        logger.debug("[QC] Rejected (high N): %s", ", ".join(rejected_too_many_n[:10]))

    return output
