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
GENOME_AVAILABLE_COLUMNS = {
    'alreadysequenced', 'already_sequenced', 'already sequenced',
    'genomeavailable', 'genome_available', 'genome available',
    'genomesequenced', 'genome_sequenced', 'genome sequenced',
    'wgsavailable', 'wgs_available', 'wgs available',
}
SELECTION_COLUMNS = {
    'selected', 'selectedforwgs', 'selected_for_wgs', 'selected for wgs',
    'selectedforgenomesequencing', 'selected_for_genome_sequencing',
    'selected for genome sequencing', 'selectedforfullgenomesequencing',
    'selected_for_full_genome_sequencing',
    'selected for full genome sequencing', 'wgs', 'wgsselected',
    'wgs_selected', 'wgs selected', 'genomesequencing',
    'genomesequencingselected', 'genome_sequencing_selected',
    'genome sequencing selected',
}


def _normalise_column(name: object) -> str:
    return str(name or '').strip().lower().replace('-', '_')


def _normalise_for_matching(name: object) -> str:
    return _normalise_column(name).replace('_', '').replace(' ', '')


def _open_text(path: str | Path):
    if str(path).endswith('.gz'):
        return gzip.open(path, 'rt', newline='')
    return open(path, 'rt', newline='')


def _parse_boolean(value: object, *, field: str) -> bool:
    if value is None:
        return False
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return float(value) != 0.0
    text = str(value).strip().lower()
    if text in ('0', 'false', 'f', 'no', 'n'):
        return False
    if text in ('1', 'true', 't', 'yes', 'y'):
        return True
    raise ValueError(f'{field} must be an explicit yes/no Boolean; received {value!r}')


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
    """Read partner genome-availability metadata from CSV/TSV sidecar files.

    Returns rows with normalised keys:
      source_id, partner_id, selected_for_sequencing, genome_available,
      selected_for_wgs (internal compatibility alias for genome_available)
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
    genome_col = _find_column(fieldnames, GENOME_AVAILABLE_COLUMNS)
    selection_col = _find_column(fieldnames, SELECTION_COLUMNS)
    partner_col = _find_column(fieldnames, PARTNER_ID_COLUMNS)

    if not id_col:
        raise ValueError(
            'Partner metadata must contain a sequence ID column matching FASTA IDs, such as sequence_id, isolate_id, sample_id, or ID.'
        )
    if not partner_col:
        raise ValueError(
            'Partner metadata must contain a partner acronym column such as partner_id, partner, or partner_acronym.'
        )
    if not genome_col:
        raise ValueError(
            'Partner metadata must contain an already-sequenced/genome-available column such as already_sequenced.'
        )

    rows = []
    seen: dict[str, tuple[str, bool, bool, int]] = {}
    for row_number, raw in enumerate(raw_rows, start=2):
        source_id = str(raw.get(id_col) or '').strip()
        if not source_id:
            continue
        partner_id = str(raw.get(partner_col) or '').strip()
        if not partner_id:
            continue
        raw_genome = raw.get(genome_col)
        raw_selection = raw.get(selection_col) if selection_col else 'no'
        genome_available = _parse_boolean(raw_genome, field='already_sequenced')
        selected_for_sequencing = _parse_boolean(
            raw_selection if raw_selection not in (None, '') else 'no',
            field='selected_for_genome_sequencing',
        )
        signature = (
            partner_id.casefold(), selected_for_sequencing, genome_available, row_number,
        )
        if source_id in seen:
            previous = seen[source_id]
            duplicate_kind = (
                'duplicate' if previous[:3] == signature[:3] else 'conflicting'
            )
            raise ValueError(
                f'Partner metadata contains a {duplicate_kind} sequence_id {source_id!r} '
                f'on rows {previous[3]} and {row_number}; the cumulative ledger requires '
                'exactly one row per isolate.'
            )
        seen[source_id] = signature
        rows.append({
            'source_id': source_id,
            'partner_id': partner_id or source_id,
            'selected_for_sequencing': selected_for_sequencing,
            'genome_available': genome_available,
            # selected_for_wgs is retained internally until the storage schema
            # can be renamed without ambiguity; it means genome available.
            'selected_for_wgs': genome_available,
            'raw_selected_value': '' if raw_genome is None else str(raw_genome),
            'raw_commitment_value': '' if raw_selection is None else str(raw_selection),
        })
    return rows
