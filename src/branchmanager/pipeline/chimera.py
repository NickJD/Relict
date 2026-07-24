"""Reference-based chimera screening for marker sequences."""

from __future__ import annotations

import shlex
from pathlib import Path

from branchmanager.pipeline.classify import _ensure_uncompressed
from branchmanager.utils.fasta import read_fasta
from branchmanager.utils.subprocess import run_cmd


def _parse_uchime_row(line: str) -> tuple[str, dict] | None:
    """Parse VSEARCH's score-first UCHIME output format."""
    parts = line.rstrip('\n').split('\t')
    if len(parts) < 2 or not parts[1]:
        return None
    call = str(parts[-1]).strip().upper()
    try:
        score = float(parts[0])
    except (TypeError, ValueError):
        score = None
    return parts[1], {
        'call': 'CHIMERA' if call == 'Y' else ('PASS' if call == 'N' else 'INDETERMINATE'),
        'score': score,
        'left_parent': parts[2] if len(parts) > 2 else '',
        'right_parent': parts[3] if len(parts) > 3 else '',
    }


def run_reference_screen(query_fasta: str, reference_fasta: str, outdir: str, *, threads: int = 4) -> tuple[str, dict[str, dict]]:
    output = Path(outdir)
    output.mkdir(parents=True, exist_ok=True)
    reference = _ensure_uncompressed(reference_fasta, outdir, out_name='chimera_reference.fasta')
    raw = output / 'chimera_uchime.tsv'
    chimeras = output / 'chimera_flagged.fasta'
    nonchimeras = output / 'chimera_passed.fasta'
    run_cmd(
        f'vsearch --uchime_ref {shlex.quote(str(query_fasta))} '
        f'--db {shlex.quote(str(reference))} --uchimeout {shlex.quote(str(raw))} '
        f'--chimeras {shlex.quote(str(chimeras))} '
        f'--nonchimeras {shlex.quote(str(nonchimeras))} --threads {max(1, int(threads))}'
    )
    parsed = {}
    if raw.is_file():
        with open(raw) as handle:
            for line in handle:
                parsed_row = _parse_uchime_row(line)
                if parsed_row is None:
                    continue
                sequence_id, result = parsed_row
                parsed[sequence_id] = result
    report = output / 'chimera_screen.tsv'
    with open(report, 'w') as handle:
        handle.write('SequenceID\tChimeraCall\tUCHIMEScore\tLeftParent\tRightParent\n')
        for sequence_id, _sequence in read_fasta(query_fasta):
            row = parsed.get(sequence_id, {'call': 'INDETERMINATE', 'score': None, 'left_parent': '', 'right_parent': ''})
            score = 'NA' if row['score'] is None else f"{row['score']:.6g}"
            handle.write(f"{sequence_id}\t{row['call']}\t{score}\t{row['left_parent']}\t{row['right_parent']}\n")
            parsed[sequence_id] = row
    return str(report), parsed
