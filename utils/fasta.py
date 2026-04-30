import gzip


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
