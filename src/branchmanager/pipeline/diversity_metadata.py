"""Metadata-driven diversity scoring helpers for selection overlays."""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List

from branchmanager.taxonomy import canonicalise_sequence_id


ID_COLUMNS = {
    'id', 'sequenceid', 'sequence_id', 'sequence id', 'seqid', 'seq_id',
    'queryid', 'query_id', 'fastaid', 'fasta_id', 'isolateid',
    'isolate_id', 'isolate id', 'sampleid', 'sample_id',
}
PARTNER_ID_COLUMNS = {
    'partnerid', 'partner_id', 'partner id', 'partner', 'partneracronym',
    'partner_acronym', 'partner acronym',
}
DIVERSITY_STRONG_THRESHOLD = 60.0
DIVERSITY_STRENGTH_RANKS = {
    'na': 0,
    'low': 0,
    'moderate': 1,
    'high': 2,
}


def _safe_text(value: object) -> str:
    if value is None:
        return ''
    return str(value).replace('\t', ' ').replace('\n', ' ').replace('\r', ' ').strip()


def _normalise_column(name: object) -> str:
    return str(name or '').strip().lower().replace('-', '_')


def _normalise_for_matching(name: object) -> str:
    return _normalise_column(name).replace('_', '').replace(' ', '')


def _normalise_value(value: object) -> str:
    return ' '.join(_safe_text(value).split()).casefold()


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


def _open_text(path: str | Path):
    if str(path).endswith('.gz'):
        import gzip
        return gzip.open(path, 'rt', newline='')
    return open(path, 'rt', newline='')


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


def load_diversity_metadata(path: str | Path) -> list[dict]:
    """Read isolate-level diversity metadata from CSV/TSV sidecar files."""
    p = Path(path)
    lower = str(p).lower()
    if lower.endswith(('.csv', '.tsv', '.txt', '.csv.gz', '.tsv.gz', '.txt.gz')):
        raw_rows = _rows_from_delimited(p)
    else:
        raise ValueError(f'Unsupported diversity metadata format: {path}. Use CSV/TSV sidecar metadata.')

    if not raw_rows:
        return []

    fieldnames = list(raw_rows[0].keys())
    id_col = _find_column(fieldnames, ID_COLUMNS)
    partner_col = _find_column(fieldnames, PARTNER_ID_COLUMNS)
    if not id_col:
        raise ValueError(
            'Diversity metadata must contain an isolate ID column matching FASTA/project IDs, such as sequence_id, isolate_id, sample_id, or ID.'
        )

    rows = []
    seen: dict[str, tuple[tuple[str, str], ...]] = {}
    for row_number, raw in enumerate(raw_rows, start=2):
        source_id = _safe_text(raw.get(id_col))
        if not source_id:
            continue
        partner_id = _safe_text(raw.get(partner_col)) if partner_col else ''
        metadata = {
            _safe_text(field): _safe_text(value)
            for field, value in raw.items()
            if field not in {id_col, partner_col} and _safe_text(value) not in {'', 'NA', 'None'}
        }
        signature = tuple(sorted((field, _normalise_value(value)) for field, value in metadata.items()))
        if source_id in seen:
            previous = seen[source_id]
            duplicate_kind = 'duplicate' if previous == signature else 'conflicting'
            raise ValueError(
                f'Diversity metadata contains a {duplicate_kind} isolate ID {source_id!r} on row {row_number}; the cumulative ledger requires exactly one row per isolate.'
            )
        seen[source_id] = signature
        rows.append({
            'source_id': source_id,
            'partner_id': partner_id or source_id,
            'metadata': metadata,
            'row_number': row_number,
        })
    return rows


def _id_key(value: object) -> str:
    text = _safe_text(value)
    if not text:
        return ''
    try:
        return canonicalise_sequence_id(text)
    except Exception:
        return text


def _score_field(total: int, count: int) -> float:
    if total <= 1 or count <= 0:
        return 0.0
    if count >= total:
        return 0.0
    return round(100.0 * (total - count) / (total - 1), 2)


def annotate_assessment_rows(rows: List[dict], metadata_rows: List[dict]) -> tuple[int, list[tuple[str, str]]]:
    """Annotate assessment rows with diversity cues and rarity scores."""
    if not rows:
        warnings = [(str(row.get('source_id') or ''), 'diversity_metadata_id_not_found_in_assessment_rows') for row in metadata_rows]
        return 0, warnings

    if not metadata_rows:
        for row in rows:
            row.update({
                'diversity_metadata_present': 'No',
                'diversity_metadata_strength': 'NA',
                'diversity_metadata_score': '0.00',
                'diversity_metadata_informative_fields': '0',
                'diversity_metadata_reason': 'no diversity metadata supplied',
            })
        return 0, []

    lookup: Dict[str, dict] = {}
    for row in rows:
        keys = {str(row.get('id') or '').strip(), _id_key(row.get('id'))}
        for key in {item for item in keys if item}:
            lookup[key] = row

    metadata_by_key: Dict[str, dict] = {}
    warnings: list[tuple[str, str]] = []
    for row in metadata_rows:
        source_id = _safe_text(row.get('source_id'))
        if not source_id:
            continue
        keys = {source_id, _id_key(source_id)}
        matched = False
        for key in {item for item in keys if item}:
            existing = metadata_by_key.get(key)
            if existing is not None and existing.get('metadata') != row.get('metadata'):
                raise ValueError(
                    f'Diversity metadata contains conflicting entries for {source_id!r}; the cumulative ledger requires one row per isolate.'
                )
            metadata_by_key[key] = row
            if key in lookup:
                matched = True
        if not matched:
            warnings.append((source_id, 'diversity_metadata_id_not_found_in_assessment_rows'))

    field_stats: dict[str, dict[str, object]] = {}
    for row in metadata_rows:
        record = row.get('metadata') or {}
        for field, value in record.items():
            if not value:
                continue
            norm = _normalise_value(value)
            if not norm:
                continue
            stats = field_stats.setdefault(field, {'total': 0, 'counts': Counter(), 'display': {}})
            stats['total'] = int(stats['total']) + 1
            counts: Counter = stats['counts']  # type: ignore[assignment]
            counts[norm] += 1
            display = stats['display']  # type: ignore[assignment]
            display.setdefault(norm, value)

    matched_count = 0
    for row in rows:
        candidate_keys = {str(row.get('id') or '').strip(), _id_key(row.get('id'))}
        metadata_row = None
        for key in {item for item in candidate_keys if item}:
            metadata_row = metadata_by_key.get(key)
            if metadata_row is not None:
                break
        if metadata_row is None:
            row.update({
                'diversity_metadata_present': 'No',
                'diversity_metadata_strength': 'NA',
                'diversity_metadata_score': '0.00',
                'diversity_metadata_informative_fields': '0',
                'diversity_metadata_reason': 'no diversity metadata supplied',
            })
            continue

        matched_count += 1
        record = metadata_row.get('metadata') or {}
        signals: list[tuple[float, str, str, int, int]] = []
        for field, value in record.items():
            if not value:
                continue
            norm = _normalise_value(value)
            if not norm:
                continue
            stats = field_stats.get(field)
            if not stats:
                continue
            total = int(stats.get('total') or 0)
            counts: Counter = stats.get('counts') or Counter()  # type: ignore[assignment]
            display = stats.get('display') or {}
            count = int(counts.get(norm, 0))
            score = _score_field(total, count)
            if score <= 0:
                continue
            signals.append((score, field, str(display.get(norm, value)), count, total))

        signals.sort(key=lambda item: (item[0], item[1]), reverse=True)
        informative = [signal for signal in signals if signal[0] >= DIVERSITY_STRONG_THRESHOLD]
        best_score = signals[0][0] if signals else 0.0
        if signals:
            if best_score >= 90.0 or len(informative) >= 2:
                strength = 'HIGH'
            elif best_score >= 60.0:
                strength = 'MODERATE'
            else:
                strength = 'LOW'
            reason = '; '.join(
                f'{field}={value} ({count}/{total})'
                for score, field, value, count, total in signals[:2]
                if score > 0
            ) or 'metadata supplied but not cohort-distinct'
        else:
            strength = 'LOW'
            reason = 'metadata supplied but not cohort-distinct'
        row.update({
            'diversity_metadata_present': 'Yes',
            'diversity_metadata_strength': strength,
            'diversity_metadata_score': f'{best_score:.2f}',
            'diversity_metadata_informative_fields': str(len(informative)),
            'diversity_metadata_reason': reason,
        })
    return matched_count, warnings
