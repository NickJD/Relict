import logging

from branchmanager.utils.fasta import read_fasta, write_fasta

logger = logging.getLogger(__name__)


def run_qc(input_fasta, outdir, min_len=800, max_n=5):
    """Filter sequences by length and allowed Ns.

    Writes filtered sequences to outdir/qc.fasta, a qc.stats summary, and
    qc_rejections.tsv with length/N counts and per-sequence rejection reasons.
    """
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
        reasons = []
        if seq_len < min_len:
            reasons.append("too_short")
        if n_count > max_n:
            reasons.append("too_many_n")

        if reasons:
            if "too_short" in reasons:
                rejected_too_short.append(header)
            if "too_many_n" in reasons:
                rejected_too_many_n.append(header)
            rejected_details.append((header, seq_len, n_count, ";".join(reasons)))
            continue
        filtered.append((header, seq))

    write_fasta(filtered, output)

    rejection_path = f"{outdir}/qc_rejections.tsv"
    with open(rejection_path, "w") as rf:
        rf.write("ID\tLength\tNCount\tReasons\tMinLength\tMaxN\n")
        for header, seq_len, n_count, reasons in rejected_details:
            rf.write(f"{header}\t{seq_len}\t{n_count}\t{reasons}\t{min_len}\t{max_n}\n")

    # write qc stats for easier debugging
    stats_path = f"{outdir}/qc.stats"
    with open(stats_path, "w") as sf:
        sf.write(f"total_input\t{total}\n")
        sf.write(f"kept\t{len(filtered)}\n")
        sf.write(f"rejected_total\t{len(rejected_details)}\n")
        sf.write(f"rejected_too_short\t{len(rejected_too_short)}\n")
        sf.write(f"rejected_too_many_n\t{len(rejected_too_many_n)}\n")
        sf.write(f"min_len\t{min_len}\n")
        sf.write(f"max_n\t{max_n}\n")
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
        "[QC] %d input → %d kept, %d rejected; %d too short (<%dbp), %d high-N (>%d). Stats: %s; details: %s",
        total, len(filtered), len(rejected_details), len(rejected_too_short), min_len,
        len(rejected_too_many_n), max_n, stats_path, rejection_path,
    )
    if rejected_too_short:
        logger.debug("[QC] Rejected (too short): %s", ", ".join(rejected_too_short[:10]))
    if rejected_too_many_n:
        logger.debug("[QC] Rejected (high N): %s", ", ".join(rejected_too_many_n[:10]))

    return output
