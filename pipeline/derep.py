from collections import defaultdict
from utils.fasta import read_fasta, write_fasta


def run_derep(input_fasta, outdir):
    output = f"{outdir}/derep.fasta"

    seq_map = {}
    clusters = defaultdict(list)
    total = 0

    for header, seq in read_fasta(input_fasta):
        total += 1
        if seq not in seq_map:
            seq_map[seq] = header
        clusters[seq].append(header)

    derep_records = [(rep, seq) for seq, rep in seq_map.items()]
    write_fasta(derep_records, output)

    # Save mapping
    map_path = f"{outdir}/derep_map.tsv"
    with open(map_path, "w") as f:
        for seq, ids in clusters.items():
            rep = seq_map[seq]
            for i in ids:
                f.write(f"{i}\t{rep}\n")

    print(f"[DEREP] Read {total} sequences, wrote {len(derep_records)} dereplicated sequences to {output}")
    print(f"[DEREP] Derep mapping written to {map_path}")

    return output

