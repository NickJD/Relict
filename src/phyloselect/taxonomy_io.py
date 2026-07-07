from __future__ import annotations

import csv
import gzip
import re
from pathlib import Path
from typing import Iterator, Optional


ID_COLUMNS = {
    'id',
    'featureid',
    'feature',
    'otuid',
    'otu',
    'asvid',
    'asv',
    'sequenceid',
    'seqid',
    'queryid',
    'accession',
    'accessionid',
}

TAXONOMY_COLUMNS = {
    'taxon',
    'taxonomy',
    'lineage',
    'classification',
    'taxonomiclineage',
}

CONFIDENCE_COLUMNS = {
    'confidence',
    'conf',
    'score',
    'probability',
}


def open_text_maybe_gzip(path: str | Path):
    if str(path).endswith('.gz'):
        return gzip.open(path, 'rt', newline='')
    return open(path, 'rt', newline='')


def _without_gz_suffix(path: str | Path) -> str:
    text = str(path).lower()
    if text.endswith('.gz'):
        text = text[:-3]
    return text


def _normalise_column_name(value: str) -> str:
    return re.sub(r'[^a-z0-9]+', '', str(value or '').strip().lower())


def _guess_delimiter(path: str | Path, sample_line: str) -> str:
    lower = _without_gz_suffix(path)
    if lower.endswith('.csv'):
        return ','
    if lower.endswith(('.tsv', '.tab', '.txt')):
        return '\t'
    return '\t' if '\t' in sample_line else ','


def _parse_fields(line: str, delimiter: str) -> list[str]:
    return next(csv.reader([line], delimiter=delimiter))


def _looks_like_header(fields: list[str]) -> bool:
    cols = {_normalise_column_name(f) for f in fields}
    return bool(cols & TAXONOMY_COLUMNS) or bool(cols & ID_COLUMNS and len(fields) > 1)


def _find_column(fields: list[str], candidates: set[str], fallback: Optional[int]) -> Optional[int]:
    for i, name in enumerate(fields):
        if _normalise_column_name(name) in candidates:
            return i
    return fallback if fallback is not None and fallback < len(fields) else None


def _parse_confidence(value: object) -> Optional[float]:
    if value is None:
        return None
    text = str(value).strip()
    if not text or text.upper() in {'NA', 'N/A', 'NONE', 'NULL'}:
        return None
    try:
        return float(text)
    except Exception:
        return None


def _clean_comment_header(line: str) -> str:
    stripped = line.strip()
    if stripped.startswith('#'):
        return stripped[1:].lstrip()
    return stripped


def iter_taxonomy_assignment_rows(path: str | Path) -> Iterator[dict[str, object]]:
    """Yield ID/taxonomy/confidence rows from TSV, CSV, TSV.GZ, or CSV.GZ.

    Supported table shapes include common QIIME-style
    ``FeatureID<TAB>Taxon<TAB>Confidence`` files, GTDB-style two-column
    ``accession<TAB>classification`` files, and equivalent comma-delimited CSVs.
    If a header row is present, column names decide which fields to use. If no
    header is present, the first three fields are interpreted as ID, taxonomy,
    and optional confidence.
    """
    with open_text_maybe_gzip(path) as handle:
        first_line = ''
        for raw in handle:
            stripped = raw.strip()
            if not stripped:
                continue
            candidate = _clean_comment_header(raw)
            lower = candidate.lower()
            if raw.lstrip().startswith('#') and not any(k in lower for k in ('taxon', 'taxonomy', 'lineage')):
                continue
            first_line = candidate
            break

        if not first_line:
            return

        delimiter = _guess_delimiter(path, first_line)
        first_fields = _parse_fields(first_line, delimiter)
        has_header = _looks_like_header(first_fields)
        if has_header:
            id_idx = _find_column(first_fields, ID_COLUMNS, 0)
            tax_idx = _find_column(first_fields, TAXONOMY_COLUMNS, 1)
            conf_idx = _find_column(first_fields, CONFIDENCE_COLUMNS, None)
        else:
            id_idx = 0
            tax_idx = 1 if len(first_fields) > 1 else None
            conf_idx = 2 if len(first_fields) > 2 else None

        def emit(fields: list[str]) -> Optional[dict[str, object]]:
            if id_idx is None or tax_idx is None:
                return None
            if len(fields) <= max(id_idx, tax_idx):
                return None
            fid = str(fields[id_idx]).strip()
            tax = str(fields[tax_idx]).strip()
            if not fid or not tax:
                return None
            conf = None
            if conf_idx is not None and conf_idx < len(fields):
                conf = _parse_confidence(fields[conf_idx])
            return {'id': fid, 'taxonomy': tax, 'confidence': conf}

        if not has_header:
            row = emit(first_fields)
            if row:
                yield row

        reader = csv.reader(handle, delimiter=delimiter)
        for fields in reader:
            if not fields:
                continue
            if len(fields) == 1 and not str(fields[0]).strip():
                continue
            if str(fields[0]).lstrip().startswith('#'):
                continue
            row = emit(fields)
            if row:
                yield row
