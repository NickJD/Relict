"""Preflight checks for a production BranchManager installation and project."""

from __future__ import annotations

import importlib
import json
import shutil
import sqlite3
import sys
import tempfile
from pathlib import Path

from branchmanager.utils.fasta import read_fasta


def _row(check: str, status: str, detail: str, required: bool = True) -> dict:
    return {'check': check, 'status': status, 'required': required, 'detail': detail}


def run_it_desk_checks(*, db_path=None, references=(), output_dir=None, tree_method='fasttree') -> list[dict]:
    rows = []
    rows.append(_row(
        'python_version', 'PASS' if sys.version_info >= (3, 10) else 'FAIL',
        sys.version.replace('\n', ' '),
    ))
    for module, required in (('Bio', True), ('matplotlib', True), ('openpyxl', False)):
        try:
            loaded = importlib.import_module(module)
            version = getattr(loaded, '__version__', 'available')
            rows.append(_row(f'python_module:{module}', 'PASS', str(version), required))
        except Exception as exc:
            rows.append(_row(f'python_module:{module}', 'FAIL' if required else 'WARN', str(exc), required))

    binaries = [('vsearch', ('vsearch',), True), ('mafft', ('mafft',), True)]
    if tree_method.startswith('iqtree'):
        binaries.append(('tree_builder', ('iqtree3', 'iqtree2', 'iqtree'), True))
    else:
        binaries.append(('tree_builder', ('FastTree', 'fasttree'), True))
    for label, candidates, required in binaries:
        found = next((shutil.which(item) for item in candidates if shutil.which(item)), None)
        rows.append(_row(label, 'PASS' if found else 'FAIL', found or f'not found: {", ".join(candidates)}', required))

    for index, reference in enumerate(references or (), start=1):
        item = Path(reference).expanduser()
        status = 'FAIL'
        detail = str(item)
        if item.is_file() and item.stat().st_size > 0:
            try:
                header, sequence = next(iter(read_fasta(str(item))))
                if not header or not sequence:
                    raise ValueError('first FASTA record is empty')
                invalid = set(sequence.upper()) - set('ACGTUNRYKMSWBDHV.-')
                if invalid:
                    raise ValueError(f'non-IUPAC symbols in first sequence: {"".join(sorted(invalid))}')
                status = 'PASS'
                detail = f'{item.resolve()} ({item.stat().st_size} bytes; first record={header})'
            except Exception as exc:
                detail = f'{item.resolve()}: could not read as nucleotide FASTA: {exc}'
        rows.append(_row(f'reference:{index}', status, detail))

    if db_path:
        item = Path(db_path).expanduser()
        if item.is_file():
            try:
                with sqlite3.connect(str(item)) as conn:
                    result = conn.execute('PRAGMA integrity_check').fetchone()[0]
                    tables = {
                        row[0] for row in conn.execute(
                            "SELECT name FROM sqlite_master WHERE type='table'"
                        )
                    }
                required_tables = {'sequences', 'taxonomy', 'sequencing_metadata'}
                missing = sorted(required_tables - tables)
                status = 'PASS' if result == 'ok' and not missing else 'FAIL'
                detail = str(result) if not missing else f'{result}; missing tables={",".join(missing)}'
                rows.append(_row('project_database', status, detail))
            except Exception as exc:
                rows.append(_row('project_database', 'FAIL', str(exc)))
        else:
            parent = item.resolve().parent
            rows.append(_row('project_database', 'PASS' if parent.is_dir() else 'FAIL', f'new database; parent={parent}'))

    if output_dir:
        parent = Path(output_dir).expanduser().resolve()
        probe_dir = parent if parent.is_dir() else parent.parent
        try:
            if not probe_dir.is_dir():
                raise FileNotFoundError(f'parent directory does not exist: {probe_dir}')
            with tempfile.NamedTemporaryFile(prefix='.branchmanager_write_test_', dir=probe_dir):
                pass
            rows.append(_row('output_directory', 'PASS', f'writable: {probe_dir}'))
        except Exception as exc:
            rows.append(_row('output_directory', 'FAIL', str(exc)))
    return rows


def write_it_desk_report(outdir: str | Path, rows: list[dict]) -> dict:
    output = Path(outdir)
    output.mkdir(parents=True, exist_ok=True)
    tsv = output / 'it_desk_report.tsv'
    with open(tsv, 'w') as handle:
        handle.write('Check\tStatus\tRequired\tDetail\n')
        for row in rows:
            detail = str(row['detail']).replace('\t', ' ').replace('\n', ' ')
            handle.write(f"{row['check']}\t{row['status']}\t{str(row['required']).lower()}\t{detail}\n")
    summary = {
        'status': 'FAIL' if any(row['required'] and row['status'] == 'FAIL' for row in rows) else 'PASS',
        'checks': rows,
    }
    json_path = output / 'it_desk_summary.json'
    json_path.write_text(json.dumps(summary, indent=2) + '\n')
    return {'tsv': str(tsv), 'json': str(json_path), 'status': summary['status']}
