"""Audited removal of sequences from active BranchManager project state."""

from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path

from branchmanager.utils.fasta import write_fasta


def _open_text(path: str | Path):
    return gzip.open(path, 'rt', newline='') if str(path).lower().endswith('.gz') else open(path, newline='')


def _normalise_header(value: object) -> str:
    return str(value or '').strip().lstrip('\ufeff').lower().replace('-', '_').replace(' ', '_')


def _row_value(row: dict[str, str], *names: str) -> str:
    for name in names:
        value = row.get(_normalise_header(name), '')
        if value:
            return str(value).strip()
    return ''


def load_exit_requests(
    sequence_ids=None, *, input_path: str | Path | None = None,
    default_reason: str = '',
) -> list[dict[str, str]]:
    """Load direct IDs and/or a CSV/TSV request table."""
    requests = [
        {
            'sequence_id': str(sequence_id or '').strip(),
            'reason': str(default_reason or '').strip(),
            'source_request': 'command_line',
        }
        for sequence_id in sequence_ids or []
    ]
    if input_path:
        source = Path(input_path).expanduser().resolve()
        with _open_text(source) as handle:
            sample = handle.read(4096)
            handle.seek(0)
            delimiter = '\t' if str(source).lower().endswith(('.tsv', '.tsv.gz', '.txt', '.txt.gz')) else ','
            try:
                delimiter = csv.Sniffer().sniff(sample, delimiters=',\t').delimiter
            except csv.Error:
                pass
            for raw_row in csv.DictReader(handle, delimiter=delimiter):
                row = {
                    _normalise_header(key): str(value or '').strip()
                    for key, value in raw_row.items() if key is not None
                }
                requests.append({
                    'sequence_id': _row_value(row, 'sequence_id', 'isolate_id', 'sample_id', 'id'),
                    'reason': _row_value(row, 'reason', 'removal_reason', 'detail', 'notes') or str(default_reason or '').strip(),
                    'source_request': str(source),
                })
    if not requests:
        raise ValueError('Provide at least one --sequence-id or an --input request table')
    missing_reason = [row.get('sequence_id') or '<blank>' for row in requests if not row.get('reason')]
    if missing_reason:
        raise ValueError(
            'Every Exit Interview requires a reason; missing for: ' + ', '.join(missing_reason)
        )
    return requests


def write_exit_interview_report(path: str | Path, rows) -> str:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = (
        'requested_id', 'sequence_id', 'dataset', 'dataset_role', 'partner_id',
        'sequence_length', 'sequence_sha256', 'taxonomy',
        'selected_for_genome_sequencing', 'genome_already_sequenced',
        'genome_records', 'status', 'reason', 'removed_records', 'downstream_action',
    )
    with open(output, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        for source in rows:
            row = dict(source)
            removed = row.get('removed_records', '')
            if isinstance(removed, dict):
                row['removed_records'] = json.dumps(removed, sort_keys=True, separators=(',', ':'))
            row['downstream_action'] = (
                'Rerun affected Performance Reviews before using tree placement, crowding, or selection rankings'
                if row.get('status') in {'READY', 'REMOVED'} else ''
            )
            writer.writerow(row)
    return str(output)


def write_departing_fasta(path: str | Path, rows) -> str:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    records = [
        (str(row['sequence_id']), str(row.get('_sequence') or ''))
        for row in rows if row.get('status') in {'READY', 'REMOVED'} and row.get('_sequence')
    ]
    write_fasta(records, str(output))
    return str(output)
