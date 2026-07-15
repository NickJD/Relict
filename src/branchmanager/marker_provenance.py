"""Bridge marker preparation QC into downstream candidate assessment."""

from __future__ import annotations

import csv
import gzip
from pathlib import Path
from typing import Optional

from branchmanager.run_manifest import sha256_file


APPROVED_REVIEW_VALUES = {'approved', 'approve', 'accepted', 'accept', 'pass', 'yes'}
REJECTED_REVIEW_VALUES = {'rejected', 'reject', 'failed', 'fail', 'no'}


def _open(path: str | Path):
    return gzip.open(path, 'rt', newline='') if str(path).lower().endswith('.gz') else open(path, newline='')


def _normalise(value: object) -> str:
    return str(value or '').strip().lstrip('\ufeff').lower().replace('-', '_').replace(' ', '_')


def _read_table(path: str | Path) -> list[dict[str, str]]:
    with _open(path) as handle:
        sample = handle.read(4096)
        handle.seek(0)
        delimiter = '\t' if str(path).lower().endswith(('.tsv', '.tsv.gz', '.txt', '.txt.gz')) else ','
        try:
            delimiter = csv.Sniffer().sniff(sample, delimiters=',\t').delimiter
        except Exception:
            pass
        return [
            {_normalise(key): str(value or '').strip() for key, value in row.items() if key is not None}
            for row in csv.DictReader(handle, delimiter=delimiter)
        ]


def discover_marker_qc(input_fasta: str | Path, explicit: Optional[str | Path] = None) -> Optional[Path]:
    if explicit:
        candidate = Path(explicit).expanduser()
        if not candidate.is_file():
            raise ValueError(f'Marker-QC sidecar not found: {candidate}')
        return candidate.resolve()
    parent = Path(input_fasta).expanduser().resolve().parent
    for name in ('assembly_report.tsv', 'marker_qc.tsv'):
        candidate = parent / name
        if candidate.is_file():
            return candidate
    return None


def _review_decisions(path: Optional[str | Path]) -> dict[str, dict[str, str]]:
    if not path:
        return {}
    decisions = {}
    for row in _read_table(path):
        sequence_id = row.get('sequence_id') or row.get('isolate_id') or row.get('id') or ''
        raw = row.get('decision') or row.get('review_status') or row.get('manual_review_status') or ''
        key = raw.strip().lower()
        if key in APPROVED_REVIEW_VALUES:
            status = 'APPROVED'
        elif key in REJECTED_REVIEW_VALUES:
            status = 'REJECTED'
        else:
            status = 'NOT_REVIEWED'
        if sequence_id:
            decisions[sequence_id] = {
                'status': status,
                'reviewer': row.get('reviewer', ''),
                'notes': row.get('notes') or row.get('detail') or '',
            }
    return decisions


def load_marker_provenance(
    input_fasta: str | Path,
    *,
    qc_path: Optional[str | Path] = None,
    review_path: Optional[str | Path] = None,
    id_map: Optional[dict[str, str]] = None,
    accept_unverified: bool = False,
) -> tuple[list[dict], Optional[str]]:
    """Return per-sequence provenance and the resolved QC sidecar path."""
    resolved = discover_marker_qc(input_fasta, qc_path)
    decisions = _review_decisions(review_path)
    mapping = id_map or {}
    input_hash = sha256_file(input_fasta)
    manifest = Path(input_fasta).resolve().parent / 'run_manifest.json'
    read_qc = Path(input_fasta).resolve().parent / 'read_qc.tsv'
    read_files_by_sequence: dict[str, list[str]] = {}
    if read_qc.is_file():
        for row in _read_table(read_qc):
            sequence_id = row.get('sequenceid') or row.get('sequence_id') or ''
            source = row.get('sourcefile') or row.get('source_file') or ''
            if sequence_id and source:
                read_files_by_sequence.setdefault(sequence_id, []).append(source)

    records = []
    if resolved:
        for row in _read_table(resolved):
            source_id = row.get('sequenceid') or row.get('sequence_id') or row.get('id') or ''
            sequence_id = mapping.get(source_id, source_id)
            if not sequence_id:
                continue
            qc_class = row.get('qcclass') or row.get('qc_class') or 'QUALITY_UNVERIFIED'
            recommendation = row.get('recommendation') or 'MANUAL_REVIEW'
            reasons = row.get('reasons') or ''
            review = decisions.get(source_id) or decisions.get(sequence_id) or {}
            records.append({
                'sequence_id': sequence_id,
                'source_manifest': str(manifest) if manifest.is_file() else '',
                'source_sequence_file': str(Path(input_fasta).resolve()),
                'source_sequence_sha256': input_hash,
                'source_read_ids': row.get('readids') or row.get('read_ids') or '',
                'source_read_files': ';'.join(sorted(set(read_files_by_sequence.get(source_id, [])))),
                'marker_qc_class': qc_class,
                'marker_qc_recommendation': recommendation,
                'marker_qc_reasons': reasons,
                'manual_review_status': review.get('status', 'NOT_REVIEWED'),
            })
    else:
        from branchmanager.utils.fasta import read_fasta
        status = 'APPROVED' if accept_unverified else 'NOT_REVIEWED'
        for source_id, _sequence in read_fasta(input_fasta):
            sequence_id = mapping.get(source_id, source_id)
            records.append({
                'sequence_id': sequence_id,
                'source_manifest': '',
                'source_sequence_file': str(Path(input_fasta).resolve()),
                'source_sequence_sha256': input_hash,
                'source_read_ids': '',
                'source_read_files': '',
                'marker_qc_class': 'QUALITY_UNVERIFIED',
                'marker_qc_recommendation': 'ACCEPTED_BY_USER' if accept_unverified else 'MANUAL_REVIEW',
                'marker_qc_reasons': 'no_marker_qc_sidecar',
                'manual_review_status': status,
            })
    return records, str(resolved) if resolved else None


def marker_qc_flag(record: Optional[dict]) -> str:
    if not record:
        return 'MARKER_QC_UNVERIFIED'
    qc_class = str(record.get('marker_qc_class') or '').upper()
    review = str(record.get('manual_review_status') or '').upper()
    if qc_class == 'FAIL_QC' or review == 'REJECTED':
        return 'MARKER_QC_FAILED'
    if qc_class == 'PASS_HIGH_CONFIDENCE':
        return ''
    if review == 'APPROVED':
        return 'MARKER_QC_REVIEW_APPROVED'
    return 'MARKER_QC_REVIEW_REQUIRED'


def write_marker_qc_bridge(path: str | Path, records: list[dict]) -> str:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        'sequence_id', 'marker_qc_class', 'marker_qc_recommendation', 'manual_review_status',
        'marker_qc_reasons', 'source_sequence_file', 'source_sequence_sha256',
        'source_read_ids', 'source_read_files', 'source_manifest', 'chimera_call', 'chimera_score',
    ]
    with open(output, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(records)
    return str(output)
