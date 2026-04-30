from utils.subprocess import run_cmd
import logging
logger = logging.getLogger(__name__)


def run_novelty(input_fasta, ref_fasta, outdir, id_threshold=0.97, threads=None):
    """Mark sequences as novel if they have no hit above id_threshold
    against the reference fasta.
    """
    output = f"{outdir}/novelty.tsv"

    if ref_fasta is None:
        raise ValueError("ref_fasta must be provided to run_novelty")

    # if ref_fasta is gzipped, create an uncompressed copy in outdir
    import os
    from utils.fasta import read_fasta, write_fasta
    ref_to_use = ref_fasta
    if str(ref_fasta).endswith('.gz'):
        ref_unc = os.path.join(outdir, 'ref_uncompressed.fasta')
        if not os.path.exists(ref_unc):
            records = list(read_fasta(ref_fasta))
            write_fasta(records, ref_unc)
        ref_to_use = ref_unc

    # If a matches.tsv already exists (produced by classification), reuse it
    # to avoid running vsearch twice. If not present, run vsearch with the
    # requested id_threshold to produce matches.tsv.
    matches_path = f"{outdir}/matches.tsv"
    if not os.path.exists(matches_path):
        thread_flag = f" --threads {int(threads)}" if threads and int(threads) > 0 else ""
        cmd = f"""
        vsearch --usearch_global {input_fasta} \
                --db {ref_to_use} \
                --id {id_threshold} \
                --blast6out {matches_path}{thread_flag}
        """
        logger.info("[NOVELTY] Running vsearch for novelty detection (input=%s, db=%s, id_threshold=%s)", input_fasta, ref_to_use, id_threshold)
        run_cmd(cmd)
        logger.info("[NOVELTY] vsearch finished; matches written to %s", matches_path)

    matched = set()
    try:
        with open(matches_path) as f:
            any_line = False
            for line in f:
                parts = line.strip().split('\t')
                if not parts:
                    continue
                any_line = True
                q = parts[0]
                # BLAST6-like output: qid, sid, identity, ... -> identity at index 2
                try:
                    ident = float(parts[2]) if len(parts) > 2 else 0.0
                except Exception:
                    ident = 0.0
                if ident >= id_threshold:
                    matched.add(q)
            if not any_line:
                print(f"[NOVELTY] Warning: matches file {matches_path} is empty")
    except FileNotFoundError:
        print(f"[NOVELTY] Warning: matches file {matches_path} not found; assuming no matches")

    # read input fasta headers
    from utils.fasta import read_fasta

    with open(output, "w") as f:
        f.write("ID\tNovel\n")
        for header, _ in read_fasta(input_fasta):
            f.write(f"{header}\t{header not in matched}\n")

    return output

