import gzip


_RC_TABLE = str.maketrans({
    'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A',
    'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K',
    'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B', 'N': 'N',
    'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'u': 'a',
    'r': 'y', 'y': 'r', 's': 's', 'w': 'w', 'k': 'm', 'm': 'k',
    'b': 'v', 'd': 'h', 'h': 'd', 'v': 'b', 'n': 'n',
})


def _open_maybe_gzip(path, mode="rt"):
    """Open a file path, transparently handling .gz compressed files.

    mode should be text mode ('rt' or 'wt') when reading/writing strings.
    """
    if str(path).endswith('.gz'):
        return gzip.open(path, mode)
    return open(path, mode)


def read_fasta(path):
    with _open_maybe_gzip(path, 'rt') as f:
        header = None
        seq = []

        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    yield header, "".join(seq)
                header = line[1:]
                seq = []
            else:
                seq.append(line)

        if header:
            yield header, "".join(seq)


def write_fasta(records, path):
    # ensure parent dir exists
    from pathlib import Path
    p = Path(path)
    if p.parent and not p.parent.exists():
        p.parent.mkdir(parents=True, exist_ok=True)

    # choose gzip writer if path endswith .gz
    mode = 'wt'
    with _open_maybe_gzip(path, mode) as f:
        for h, s in records:
            f.write(f">{h}\n{s}\n")


def reverse_complement(seq: str) -> str:
    if seq is None:
        return ''
    return str(seq).translate(_RC_TABLE)[::-1]

