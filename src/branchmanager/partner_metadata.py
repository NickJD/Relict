from __future__ import annotations

import csv
import gzip
from pathlib import Path
from typing import Iterable


ID_COLUMNS = {
    'id', 'sequenceid', 'sequence_id', 'sequence id', 'seqid', 'seq_id',
    'queryid', 'query_id', 'fastaid', 'fasta_id', 'isolateid',
    'isolate_id', 'isolate id', 'sampleid', 'sample_id',
}
PARTNER_ID_COLUMNS = {
    'partnerid', 'partner_id', 'partner id', 'partner', 'partneracronym',
    'partner_acronym', 'partner acronym',
}
SELECTED_COLUMNS = {
    'selected', 'selectedforwgs', 'selected_for_wgs', 'selected for wgs',
    'selectedforgenomesequencing', 'selected_for_genome_sequencing',
    'selected for genome sequencing', 'selectedforfullgenomesequencing',
    'selected_for_full_genome_sequencing',
    'selected for full genome sequencing', 'wgs', 'wgsselected',
    'wgs_selected', 'wgs selected', 'genomesequencing',
    'genome_sequencing', 'genome sequencing', 'fullgenomesequencing',
    'full_genome_sequencing', 'full genome sequencing',
}


def _normalise_column(name: object) -> str:
    return str(name or '').strip().lower().replace('-', '_')


def _normalise_for_matching(name: object) -> str:
    return _normalise_column(name).replace('_', '').replace(' ', '')


def _open_text(path: str | Path):
    if str(path).endswith('.gz'):
        return gzip.open(path, 'rt', newline='')
    return open(path, 'rt', newline='')


def _truthy_selected(value: object) -> bool:
    if value is None:
        return False
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return float(value) != 0.0
    text = str(value).strip().lower()
    if text in ('', '0', 'false', 'f', 'no', 'n', 'none', 'na', 'n/a', 'not selected', 'not sequenced'):
        return False
    if text.startswith('not ') or text.startswith('no '):
        return False
    if text in (
        '1', 'true', 't', 'yes', 'y', 'selected', 'select', 'wgs',
        'genome', 'genome sequencing', 'full genome sequencing',
        'sequenced', 'for sequencing',
    ):
        return True
    return any(marker in text for marker in ('selected', 'wgs', 'sequenc'))


def _find_column(fieldnames: Iterable[str], candidates: set[str]) -> str | None:
    fields = list(fieldnames or [])
    exact = {_normalise_column(c): c for c in candidates}
    compact = {_normalise_for_matching(c): c for c in candidates}
    for field in fields:
        norm = _normalise_column(field)
        if norm in exact:
            return field
    for field in fields:
        norm = _normalise_for_matching(field)
        if norm in compact:
            return field
    return None


def _rows_from_delimited(path: str | Path) -> list[dict]:
    suffixes = ''.join(Path(str(path).removesuffix('.gz')).suffixes).lower()
    delimiter = '\t' if suffixes.endswith('.tsv') or suffixes.endswith('.txt') else ','
    with _open_text(path) as handle:
        sample = handle.read(4096)
        handle.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=',\t;')
            delimiter = dialect.delimiter
        except Exception:
            pass
        reader = csv.DictReader(handle, delimiter=delimiter)
        return [row for row in reader if any(v not in (None, '') for v in row.values())]


def load_partner_sequencing_metadata(path: str | Path, sheet_name: str | None = None) -> list[dict]:
    """Read partner sequencing-selection metadata from CSV/TSV sidecar files.

    Returns rows with normalized keys:
      source_id, partner_id, selected_for_wgs, raw_selected_value
    """
    p = Path(path)
    lower = str(p).lower()
    if lower.endswith(('.csv', '.tsv', '.txt', '.csv.gz', '.tsv.gz', '.txt.gz')):
        raw_rows = _rows_from_delimited(p)
    else:
        raise ValueError(f"Unsupported partner metadata format: {path}. Use CSV/TSV sidecar metadata.")

    if not raw_rows:
        return []

    fieldnames = list(raw_rows[0].keys())
    id_col = _find_column(fieldnames, ID_COLUMNS)
    selected_col = _find_column(fieldnames, SELECTED_COLUMNS)
    partner_col = _find_column(fieldnames, PARTNER_ID_COLUMNS)

    if not id_col:
        raise ValueError(
            'Partner metadata must contain a sequence ID column matching FASTA IDs, such as sequence_id, isolate_id, sample_id, or ID.'
        )
    if not partner_col:
        raise ValueError(
            'Partner metadata must contain a partner acronym column such as partner_id, partner, or partner_acronym.'
        )
    if not selected_col:
        raise ValueError(
            'Partner metadata must contain a WGS-selection column such as selected_for_wgs or selected_for_genome_sequencing.'
        )

    rows = []
    for raw in raw_rows:
        source_id = str(raw.get(id_col) or '').strip()
        if not source_id:
            continue
        partner_id = str(raw.get(partner_col) or '').strip()
        if not partner_id:
            continue
        raw_selected = raw.get(selected_col)
        rows.append({
            'source_id': source_id,
            'partner_id': partner_id or source_id,
            'selected_for_wgs': _truthy_selected(raw_selected),
            'raw_selected_value': '' if raw_selected is None else str(raw_selected),
        })
    return rows
