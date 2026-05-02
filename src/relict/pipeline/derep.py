import hashlib
import logging
from collections import defaultdict

from relict.utils.fasta import read_fasta, write_fasta

logger = logging.getLogger(__name__)


def run_derep(input_fasta, outdir):
    """Deduplicate sequences by exact sequence content.

    Uses MD5 hashing to avoid holding full sequences as dictionary keys,
    making this memory-efficient for large input datasets.
    """
    output = f"{outdir}/derep.fasta"

    # md5_hash -> (representative_header, sequence)
    seq_hash_to_rep: dict[str, tuple[str, str]] = {}
    # md5_hash -> [all headers with this sequence]
    clusters: dict[str, list[str]] = defaultdict(list)
    total = 0

    for header, seq in read_fasta(input_fasta):
        total += 1
        md5 = hashlib.md5(seq.upper().encode()).hexdigest()
        if md5 not in seq_hash_to_rep:
            seq_hash_to_rep[md5] = (header, seq)
        clusters[md5].append(header)

    derep_records = [rep_seq for rep_seq in seq_hash_to_rep.values()]
    write_fasta(derep_records, output)

    # Save mapping: member -> representative
    map_path = f"{outdir}/derep_map.tsv"
    with open(map_path, "w") as f:
        f.write("member\trepresentative\n")
        for md5, members in clusters.items():
            rep_header, _ = seq_hash_to_rep[md5]
            for member in members:
                f.write(f"{member}\t{rep_header}\n")

    n_clusters = len(derep_records)
    n_collapsed = total - n_clusters
    logger.info(
        "[DEREP] %d input → %d unique sequences (%d duplicates removed). Map: %s",
        total, n_clusters, n_collapsed, map_path,
    )
    return output
