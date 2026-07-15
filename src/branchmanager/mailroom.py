"""AB1 batch inventory and supplier-metadata reconciliation before Onboarding."""

from __future__ import annotations

import csv
import json
import os
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Optional

from branchmanager.onboarding import _read_table
from branchmanager.pipeline.paper_trail import (
    DEFAULT_PRIMER_SEQUENCES,
    _decode_abif_text,
    _read_ab1_details,
    _read_abif_entries,
    _row_processing_mode,
    _row_value,
)


AB1_SUFFIXES = ('.ab1', '.abi', '.ab1.gz', '.abi.gz')
SUPPLIER_ID_COLUMNS = (
    'sequencing_id', 'supplier_read_id', 'supplier_id', 'trace_id',
    'read_id', 'ab1_id', 'file', 'filename', 'read_file', 'ab1_file',
)
DIRECTION_COLUMNS = ('direction', 'orientation', 'read', 'read_direction')
DATASET_COLUMNS = ('dataset', 'batch_id', 'batch', 'submission_id')


def _strip_ab1_suffix(value: str | Path) -> str:
    name = Path(str(value)).name
    lower = name.lower()
    for suffix in sorted(AB1_SUFFIXES, key=len, reverse=True):
        if lower.endswith(suffix):
            return name[:-len(suffix)]
    return Path(name).stem


def _normalise_supplier_id(value: str | Path) -> str:
    stem = _strip_ab1_suffix(value).strip()
    return re.sub(
        r'[-_](?:f|r|forward|reverse)$', '', stem, flags=re.IGNORECASE,
    ).casefold()


def _normalise_direction(value: object) -> str:
    key = re.sub(r'[^a-z]+', '', str(value or '').strip().lower())
    if key in {'f', 'fwd', 'forward'}:
        return 'forward'
    if key in {'r', 'rev', 'reverse'}:
        return 'reverse'
    return ''


def _relative_path(path: Path, output: Path) -> str:
    try:
        return os.path.relpath(path.resolve(), output.resolve())
    except (OSError, ValueError):
        return str(path.resolve())


def _embedded_primer(entries: dict[tuple[str, int], bytes]) -> str:
    names = sorted(DEFAULT_PRIMER_SEQUENCES, key=len, reverse=True)
    found = set()
    for (tag, _number), value in entries.items():
        if tag in {'PBAS', 'PCON', 'PLOC', 'DATA'}:
            continue
        text = _decode_abif_text(value)
        for name in names:
            if re.search(rf'(?<![A-Za-z0-9]){re.escape(name)}(?![A-Za-z0-9])', text, re.IGNORECASE):
                found.add(name)
    return next(iter(found)) if len(found) == 1 else ''


def _inspect_trace(path: Path, output: Path) -> dict[str, object]:
    row: dict[str, object] = {
        'read_file': _relative_path(path, output),
        'filename': path.name,
        'supplier_prefix': _normalise_supplier_id(path.name).split('_', 1)[0],
        'abif_sample': '',
        'abif_comment': '',
        'abif_run_name': '',
        'instrument_model': '',
        'raw_length': '',
        'quality_available': '',
        'trace_available': '',
        'embedded_primer': '',
        'parse_status': 'PASS',
        'parse_error': '',
        '_path': path.resolve(),
    }
    try:
        entries = _read_abif_entries(path)
        details = _read_ab1_details(path)
        row.update({
            'abif_sample': _decode_abif_text(entries.get(('SMPL', 1), b'')),
            'abif_comment': _decode_abif_text(entries.get(('CMNT', 1), b'')),
            'abif_run_name': _decode_abif_text(entries.get(('RunN', 1), b'')),
            'instrument_model': _decode_abif_text(entries.get(('MODL', 1), b'')),
            'raw_length': len(str(details.get('sequence') or '')),
            'quality_available': 'yes' if details.get('quality_available') else 'no',
            'trace_available': 'yes' if details.get('trace_available') else 'no',
            'embedded_primer': _embedded_primer(entries),
        })
    except Exception as exc:
        row['parse_status'] = 'FAIL'
        row['parse_error'] = str(exc)
    return row


def _collect_ab1_files(read_dir: str | Path, recursive: bool) -> list[Path]:
    root = Path(read_dir).expanduser().resolve()
    if not root.is_dir():
        raise ValueError(f'AB1 directory does not exist: {root}')
    iterator: Iterable[Path] = root.rglob('*') if recursive else root.iterdir()
    return sorted(
        child.resolve() for child in iterator
        if child.is_file() and str(child).lower().endswith(AB1_SUFFIXES)
    )


def _build_file_index(inventory: list[dict[str, object]]) -> dict[str, set[Path]]:
    index: dict[str, set[Path]] = defaultdict(set)
    for row in inventory:
        path = Path(row['_path'])
        stem = _strip_ab1_suffix(path)
        keys = {
            path.name.casefold(), stem.casefold(), _normalise_supplier_id(stem),
            stem.split('_', 1)[0].casefold(),
            str(row.get('abif_sample') or '').strip().casefold(),
        }
        for key in keys:
            if key:
                index[key].add(path)
    return index


def _resolve_supplier_read(
    value: str,
    read_dir: Path,
    index: dict[str, set[Path]],
) -> list[Path]:
    supplied = Path(str(value)).expanduser()
    direct = supplied if supplied.is_absolute() else read_dir / supplied
    if direct.is_file():
        return [direct.resolve()]
    stem = _strip_ab1_suffix(value)
    keys = {
        str(value).strip().casefold(), Path(str(value)).name.casefold(),
        stem.casefold(), _normalise_supplier_id(stem),
    }
    matches: set[Path] = set()
    for key in keys:
        matches.update(index.get(key, set()))
    return sorted(matches)


def prepare_ab1_map(
    read_dir: str | Path,
    submission_metadata: str | Path,
    outdir: str | Path,
    *,
    dataset: str,
    forward_primer: Optional[str] = None,
    reverse_primer: Optional[str] = None,
    processing_mode: str = 'auto',
    recursive: bool = True,
) -> dict[str, object]:
    """Inspect an AB1 batch and write a validated one-row-per-read map."""
    output = Path(outdir).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    read_root = Path(read_dir).expanduser().resolve()
    rows = _read_table(submission_metadata)
    files = _collect_ab1_files(read_root, recursive)
    inventory = [_inspect_trace(path, output) for path in files]
    inventory_by_path = {Path(row['_path']): row for row in inventory}
    file_index = _build_file_index(inventory)
    report: list[dict[str, object]] = []
    provisional: list[dict[str, object]] = []
    expected_reads = Counter()
    file_owner: dict[Path, str] = {}

    for line_number, row in enumerate(rows, start=2):
        sequence_id = _row_value(
            row, 'sequence_id', 'isolate_id', 'isolate_number', 'sample_id',
            'unique_id', 'sample', 'isolate',
        )
        supplier_id = _row_value(row, *SUPPLIER_ID_COLUMNS)
        direction = _normalise_direction(_row_value(row, *DIRECTION_COLUMNS))
        supplied_dataset = _row_value(row, *DATASET_COLUMNS)
        row_dataset = str(dataset).strip()
        supplied_primer = _row_value(row, 'primer', 'primer_name')
        explicit_mode = _row_processing_mode(row)

        if not sequence_id:
            report.append(_issue('ERROR', line_number, '', supplier_id, 'MISSING_SEQUENCE_ID', 'supplier metadata requires an isolate/sequence ID'))
            continue
        expected_reads[sequence_id] += 1
        if supplied_dataset and supplied_dataset.casefold() != row_dataset.casefold():
            report.append(_issue('ERROR', line_number, sequence_id, supplier_id, 'CONFLICTING_DATASET', f'supplier row uses {supplied_dataset!r}; Mailroom batch is {row_dataset!r}'))
        if not supplier_id:
            report.append(_issue('ERROR', line_number, sequence_id, '', 'MISSING_SUPPLIER_READ_ID', 'supplier metadata requires a sequencing/read/file identifier'))
            continue
        if not direction and supplied_primer:
            direction = 'reverse' if supplied_primer.upper().endswith('R') else ('forward' if supplied_primer.upper().endswith('F') else '')
        if not direction:
            report.append(_issue('ERROR', line_number, sequence_id, supplier_id, 'MISSING_DIRECTION', 'read direction must be forward or reverse'))

        matches = _resolve_supplier_read(supplier_id, read_root, file_index)
        if not matches:
            report.append(_issue('ERROR', line_number, sequence_id, supplier_id, 'READ_FILE_NOT_FOUND', f'no AB1 file matched {supplier_id!r}'))
            continue
        if len(matches) > 1:
            report.append(_issue('ERROR', line_number, sequence_id, supplier_id, 'AMBIGUOUS_READ_FILE', ';'.join(str(path) for path in matches)))
            continue
        path = matches[0]
        previous_owner = file_owner.get(path)
        if previous_owner:
            report.append(_issue('ERROR', line_number, sequence_id, supplier_id, 'READ_FILE_MULTIPLE_ASSIGNMENTS', f'{path} is already assigned to {previous_owner}'))
            continue
        file_owner[path] = sequence_id

        embedded_primer = str(inventory_by_path[path].get('embedded_primer') or '')
        configured_primer = forward_primer if direction == 'forward' else reverse_primer if direction == 'reverse' else None
        if supplied_primer:
            primer = supplied_primer
            primer_assignment = 'supplied_by_metadata'
        elif embedded_primer:
            primer = embedded_primer
            primer_assignment = 'embedded_in_abif'
        elif configured_primer:
            primer = str(configured_primer)
            primer_assignment = 'configured_for_batch'
        else:
            primer = 'unknown'
            primer_assignment = 'unresolved'
            report.append(_issue('WARNING', line_number, sequence_id, supplier_id, 'UNRESOLVED_PRIMER', 'AB1 contains no primer name; supply it in metadata or with --forward-primer/--reverse-primer'))

        primer_upper = primer.upper()
        if (
            primer != 'unknown'
            and ((primer_upper.endswith('R') and direction == 'forward')
                 or (primer_upper.endswith('F') and direction == 'reverse'))
        ):
            report.append(_issue('ERROR', line_number, sequence_id, supplier_id, 'PRIMER_DIRECTION_CONFLICT', f'primer {primer!r} conflicts with {direction!r} direction'))

        if inventory_by_path[path]['parse_status'] != 'PASS':
            report.append(_issue('ERROR', line_number, sequence_id, supplier_id, 'AB1_PARSE_FAILED', str(inventory_by_path[path]['parse_error'])))

        provisional.append({
            'sequence_id': sequence_id,
            'dataset': row_dataset,
            'read_file': _relative_path(path, output),
            'primer': primer,
            'direction': direction or 'unknown',
            'processing_mode': explicit_mode or '',
            'primer_assignment': primer_assignment,
        })

    grouped = defaultdict(list)
    for row in provisional:
        grouped[row['sequence_id']].append(row)
    for sequence_id, group in grouped.items():
        explicit_modes = {row['processing_mode'] for row in group if row['processing_mode']}
        if len(explicit_modes) > 1:
            report.append(_issue('ERROR', '', sequence_id, '', 'CONFLICTING_PROCESSING_MODE', ','.join(sorted(explicit_modes))))
        selected_mode = (
            processing_mode if processing_mode in {'assemble', 'best_read'}
            else next(iter(explicit_modes)) if len(explicit_modes) == 1
            else 'assemble' if expected_reads[sequence_id] > 1
            else 'best_read'
        )
        for row in group:
            row['processing_mode'] = selected_mode
        directions = {row['direction'] for row in group}
        if expected_reads[sequence_id] > 1 and directions != {'forward', 'reverse'}:
            report.append(_issue('WARNING', '', sequence_id, '', 'INCOMPLETE_DIRECTION_PAIR', f'expected forward+reverse reads; observed {",".join(sorted(directions))}'))

    mapped_paths = set(file_owner)
    for path in files:
        if path not in mapped_paths:
            report.append(_issue('WARNING', '', '', '', 'UNMAPPED_AB1_FILE', str(path)))

    map_path = output / 'ab1_map.tsv'
    inventory_path = output / 'ab1_inventory.tsv'
    report_path = output / 'mailroom_report.tsv'
    summary_path = output / 'mailroom_summary.json'
    _write_tsv(map_path, provisional, (
        'sequence_id', 'dataset', 'read_file', 'primer', 'direction',
        'processing_mode', 'primer_assignment',
    ))
    _write_tsv(inventory_path, inventory, (
        'read_file', 'filename', 'supplier_prefix', 'abif_sample', 'abif_comment',
        'abif_run_name', 'instrument_model', 'raw_length', 'quality_available',
        'trace_available', 'embedded_primer', 'parse_status', 'parse_error',
    ))
    _write_tsv(report_path, report, (
        'severity', 'line', 'sequence_id', 'supplier_read_id', 'code', 'detail',
    ))

    errors = sum(row['severity'] == 'ERROR' for row in report)
    warnings = sum(row['severity'] == 'WARNING' for row in report)
    unresolved = sum(row['primer_assignment'] == 'unresolved' for row in provisional)
    status = 'FAIL' if errors else 'REVIEW_REQUIRED' if unresolved else 'PASS'
    summary = {
        'status': status,
        'dataset': str(dataset),
        'metadata_rows': len(rows),
        'physical_ab1_files': len(files),
        'mapped_reads': len(provisional),
        'isolates': len(expected_reads),
        'mapped_isolates': len(grouped),
        'unresolved_primers': unresolved,
        'errors': errors,
        'warnings': warnings,
    }
    summary_path.write_text(json.dumps(summary, indent=2) + '\n')
    return {
        **summary,
        'ab1_map': str(map_path),
        'inventory': str(inventory_path),
        'report': str(report_path),
        'summary': str(summary_path),
    }


def _issue(
    severity: str,
    line: object,
    sequence_id: str,
    supplier_read_id: str,
    code: str,
    detail: str,
) -> dict[str, object]:
    return {
        'severity': severity,
        'line': line,
        'sequence_id': sequence_id,
        'supplier_read_id': supplier_read_id,
        'code': code,
        'detail': detail,
    }


def _write_tsv(path: Path, rows: list[dict[str, object]], fields: tuple[str, ...]) -> None:
    with open(path, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
