from __future__ import annotations

import csv
import difflib
import gzip
import logging
import math
import re
import struct
from dataclasses import dataclass
from html import escape as html_escape
from pathlib import Path
from typing import Iterable, Optional

from branchmanager.utils.fasta import read_fasta, reverse_complement, write_fasta


logger = logging.getLogger(__name__)


DEFAULT_PRIMERS = (
    '27F', '341F', '515F', '519F', '785F',
    '534R', '806R', '907R', '926R', '1100R', '1392R', '1406R', '1492R',
)

SUPPORTED_READ_EXTENSIONS = (
    '.ab1', '.abi', '.ab1.gz', '.abi.gz',
    '.fasta', '.fa', '.fna', '.fasta.gz', '.fa.gz', '.fna.gz',
    '.fastq', '.fq', '.fastq.gz', '.fq.gz',
)

READ_FILE_COLUMNS = {
    'file', 'filename', 'read_file', 'read_filename', 'ab1_file', 'abi_file',
    'trace_file', 'chromatogram_file', 'path', 'filepath', 'file_path',
    'read_id', 'trace_id', 'sequencing_id',
}
READ_FILE_LIST_COLUMNS = {
    'files', 'filenames', 'read_files', 'read_filenames', 'ab1_files', 'abi_files',
    'trace_files', 'chromatogram_files', 'paths', 'filepaths', 'file_paths',
}
SEQUENCE_ID_COLUMNS = (
    'sequence_id', 'seq_id', 'id', 'isolate_id', 'sample_id',
    'unique_id', 'uniqueid', 'unique_sequence_id',
    'sample', 'isolate', 'isolate_number', 'sample_number',
)
PROCESSING_MODE_COLUMNS = (
    'mode', 'processing_mode', 'process_mode', 'assembly_mode', 'ab1_mode',
    'sanger_mode', 'read_mode', 'read_handling', 'read_strategy',
    'sequence_mode', 'action', 'operation', 'workflow',
)
PROCESSING_TAG_COLUMNS = ('tag', 'tags', 'flag', 'flags')
ROW_METADATA_COLUMNS = set(SEQUENCE_ID_COLUMNS) | {
    'primer', 'primer_name', 'direction', 'orientation', 'read', 'read_direction',
} | set(PROCESSING_MODE_COLUMNS) | set(PROCESSING_TAG_COLUMNS)


@dataclass
class SangerRead:
    read_id: str
    sequence_id: str
    primer: str
    direction: str
    source_file: str
    raw_sequence: str
    raw_qualities: list[int]
    quality_available: bool = True
    quality_source: str = 'provided'
    trimmed_sequence: str = ''
    trimmed_qualities: list[int] | None = None
    trim_start: int = 0
    trim_end: int = 0
    masked_sequence: str = ''
    masked_qualities: list[int] | None = None
    masked_bases: int = 0
    longest_low_quality_run: int = 0
    oriented_sequence: str = ''
    oriented_qualities: list[int] | None = None
    status: str = 'pending'
    warning: str = ''
    qc_class: str = 'PENDING'
    qc_reasons: list[str] | None = None
    processing_mode: str = 'assemble'
    processing_mode_explicit: bool = False


def _open_text(path: str | Path):
    if str(path).endswith('.gz'):
        return gzip.open(path, 'rt', newline='')
    return open(path, 'rt', newline='')


def _read_abif_entries(path: str | Path) -> dict[tuple[str, int], bytes]:
    if str(path).lower().endswith('.gz'):
        with gzip.open(path, 'rb') as handle:
            data = handle.read()
    else:
        data = Path(path).read_bytes()
    if data[:4] != b'ABIF':
        raise ValueError(f'Not an ABIF/AB1 file: {path}')
    if len(data) < 34:
        raise ValueError(f'Truncated ABIF/AB1 file: {path}')

    def parse_entry(offset: int):
        if offset + 28 > len(data):
            raise ValueError(f'Truncated ABIF directory entry in {path}')
        tag, tag_num, elem_type, elem_size, elem_count, data_size, data_offset, data_handle = struct.unpack(
            '>4sIHHIIII',
            data[offset:offset + 28],
        )
        inline = data[offset + 20:offset + 24]
        return {
            'tag': tag.decode('ascii', errors='replace'),
            'tag_num': tag_num,
            'elem_type': elem_type,
            'elem_size': elem_size,
            'elem_count': elem_count,
            'data_size': data_size,
            'data_offset': data_offset,
            'data_handle': data_handle,
            'inline': inline,
        }

    root = parse_entry(6)
    directory_offset = root['data_offset']
    directory_count = root['elem_count']
    entries: dict[tuple[str, int], bytes] = {}
    for idx in range(directory_count):
        entry_offset = directory_offset + idx * 28
        entry = parse_entry(entry_offset)
        if entry['data_size'] <= 4:
            value = entry['inline'][:entry['data_size']]
        else:
            start = entry['data_offset']
            end = start + entry['data_size']
            value = data[start:end]
        entries[(entry['tag'], entry['tag_num'])] = value
    return entries


def _read_ab1_with_quality_flag(path: str | Path) -> tuple[str, list[int], bool]:
    entries = _read_abif_entries(path)
    bases = entries.get(('PBAS', 2)) or entries.get(('PBAS', 1))
    qualities = entries.get(('PCON', 2)) or entries.get(('PCON', 1))
    if bases is None:
        raise ValueError(f'AB1 file has no PBAS base-call tag: {path}')
    seq = bases.decode('ascii', errors='ignore').replace('\x00', '').strip().upper()
    seq = ''.join(base if base in 'ACGTRYSWKMBDHVN' else 'N' for base in seq)
    if qualities is None:
        qual = [0] * len(seq)
        quality_available = False
    else:
        qual = [int(q) for q in qualities[:len(seq)]]
        if len(qual) < len(seq):
            qual.extend([0] * (len(seq) - len(qual)))
        quality_available = True
    return seq, qual, quality_available


def read_ab1(path: str | Path) -> tuple[str, list[int]]:
    """Read base calls and phred-like qualities from an AB1/ABIF file."""
    seq, qual, _quality_available = _read_ab1_with_quality_flag(path)
    return seq, qual


def read_fastq(path: str | Path) -> list[tuple[str, str, list[int]]]:
    records = []
    with _open_text(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline()
            plus = handle.readline()
            qual = handle.readline()
            if not qual:
                break
            if not header.startswith('@') or not plus.startswith('+'):
                raise ValueError(f'Malformed FASTQ near header {header.strip()} in {path}')
            rid = header[1:].strip()
            sequence = seq.strip().upper()
            qualities = [max(0, ord(ch) - 33) for ch in qual.rstrip('\n\r')]
            if len(qualities) < len(sequence):
                qualities.extend([0] * (len(sequence) - len(qualities)))
            records.append((rid, sequence, qualities[:len(sequence)]))
    return records


def trim_by_quality(
    sequence: str,
    qualities: list[int],
    *,
    min_quality: int = 20,
    window: int = 20,
) -> tuple[str, list[int], int, int]:
    """Return the highest-scoring Phred/Mott trim interval.

    A base contributes positively when its error probability is below the
    selected cutoff and negatively when it is above it. This preserves high-
    quality internal sequence while trimming unreliable ends, matching the
    usual Sanger trace-cleanup model more closely than a pure moving average.
    """
    seq = str(sequence or '').upper()
    qual = list(qualities or [])
    if len(qual) < len(seq):
        qual.extend([0] * (len(seq) - len(qual)))
    qual = qual[:len(seq)]
    n = len(seq)
    if n == 0:
        return '', [], 0, 0
    cutoff_error = _error_probability(min_quality)
    best_score = 0.0
    best_start = 0
    best_end = 0
    score = 0.0
    start = 0
    for idx, q in enumerate(qual):
        base_score = cutoff_error - _error_probability(q)
        if score <= 0:
            score = base_score
            start = idx
        else:
            score += base_score
        if score > best_score:
            best_score = score
            best_start = start
            best_end = idx + 1
    if best_end <= best_start:
        return '', [], 0, 0

    left = best_start
    right = best_end

    while left < right and seq[left] == 'N':
        left += 1
    while right > left and seq[right - 1] == 'N':
        right -= 1

    return seq[left:right], qual[left:right], left, right


def infer_read_info(path_or_id: str, primers: Iterable[str] = DEFAULT_PRIMERS) -> tuple[str, str, str]:
    stem = Path(str(path_or_id)).stem
    tokens = [tok for tok in re.split(r'[_\-. ]+', stem) if tok]
    primer = ''
    primer_index = None
    primer_lookup = {p.upper(): p for p in primers}
    for idx, token in enumerate(tokens):
        got = primer_lookup.get(token.upper())
        if got:
            primer = got
            primer_index = idx
            break
    if not primer:
        for item in primers:
            if re.search(rf'(^|[_\-. ]){re.escape(item)}($|[_\-. ])', stem, flags=re.IGNORECASE):
                primer = item
                break
    direction = 'reverse' if primer.upper().endswith('R') else 'forward'
    if primer_index is not None:
        seq_tokens = tokens[:primer_index] + tokens[primer_index + 1:]
        sequence_id = '_'.join(seq_tokens) or stem
    elif primer:
        sequence_id = re.sub(rf'(^|[_\-. ]){re.escape(primer)}($|[_\-. ])', '_', stem, flags=re.IGNORECASE).strip('_-. ')
        sequence_id = sequence_id or stem
    else:
        sequence_id = stem
    return sequence_id, primer or 'unknown', direction


def _normalise_column(name: object) -> str:
    return str(name or '').strip().lstrip('\ufeff').lower().replace('-', '_').replace(' ', '_')


def _normalise_row(row: dict) -> dict[str, str]:
    norm = {}
    for key, value in row.items():
        if key is None:
            continue
        norm[_normalise_column(key)] = '' if value is None else str(value).strip()
    return norm


def _row_value(row: dict[str, str], *names: str) -> str:
    for name in names:
        got = row.get(name)
        if got not in (None, ''):
            return str(got).strip()
    return ''


def normalise_processing_mode(value: object, default: Optional[str] = None) -> Optional[str]:
    """Map user-facing Sanger processing tags onto stable workflow modes."""
    text = str(value or '').strip().lower()
    if not text:
        return default
    key = re.sub(r'[^a-z0-9]+', '_', text).strip('_')
    best_markers = {
        'best', 'best_read', 'bestread', 'quality_best', 'best_quality',
        'highest_quality', 'select_best', 'choose_best', 'single_best',
        'pick_best', 'independent', 'independently', 'separate',
        'separate_reads', 'individual', 'individual_reads', 'no_assemble',
        'noassembly', 'no_assembly', 'dont_assemble', 'do_not_assemble',
        'unassembled', 'convert_only', 'fasta_only',
    }
    assemble_markers = {
        'assemble', 'assembled', 'assembly', 'merge', 'merged',
        'consensus', 'overlap', 'assemble_reads', 'merge_reads',
        'long_read', 'longer_sequence',
    }

    if key in best_markers:
        return 'best_read'
    if key in assemble_markers:
        return 'assemble'

    chunks = [
        re.sub(r'[^a-z0-9]+', '_', chunk).strip('_')
        for chunk in re.split(r'[;,|/]+', text)
        if chunk.strip()
    ]
    for chunk in chunks:
        if chunk in best_markers or any(marker in chunk for marker in best_markers):
            return 'best_read'
    for chunk in chunks:
        if chunk in assemble_markers or any(marker in chunk for marker in assemble_markers):
            return 'assemble'

    if any(marker in key for marker in best_markers):
        return 'best_read'
    if any(marker in key for marker in assemble_markers):
        return 'assemble'
    return default


def _row_processing_mode(row: dict[str, str]) -> Optional[str]:
    for column in PROCESSING_MODE_COLUMNS:
        mode = normalise_processing_mode(row.get(column), default=None)
        if mode:
            return mode
    for column in PROCESSING_TAG_COLUMNS:
        mode = normalise_processing_mode(row.get(column), default=None)
        if mode:
            return mode
    return None


def _split_file_cell(value: str) -> list[str]:
    if not value:
        return []
    cleaned = str(value).strip().strip('"').strip("'")
    if not cleaned:
        return []
    return [part.strip().strip('"').strip("'") for part in re.split(r'\s*[;,|]\s*', cleaned) if part.strip()]


def _is_supported_read_file(path: str | Path) -> bool:
    lower = str(path).lower()
    return lower.endswith(SUPPORTED_READ_EXTENSIONS)


def _primer_from_column(column: str, primers: Iterable[str]) -> str:
    normalised = _normalise_column(column)
    primer_lookup = {_normalise_column(p): p for p in primers}
    for key, primer in sorted(primer_lookup.items(), key=lambda item: -len(item[0])):
        if re.search(rf'(^|_){re.escape(key)}($|_)', normalised):
            return primer
    return ''


def _metadata_file_columns(row: dict[str, str], primers: Iterable[str]) -> list[tuple[str, str]]:
    columns: list[tuple[str, str]] = []
    for column, value in row.items():
        if not value or column in ROW_METADATA_COLUMNS:
            continue
        primer_from_column = _primer_from_column(column, primers)
        is_file_column = (
            column in READ_FILE_COLUMNS
            or column in READ_FILE_LIST_COLUMNS
            or column.startswith(('file_', 'read_file_', 'ab1_file_', 'abi_file_', 'trace_file_', 'chromatogram_file_'))
            or column.endswith(('_file', '_files', '_path', '_paths'))
            or bool(primer_from_column)
        )
        if not is_file_column and not _is_supported_read_file(value):
            continue
        for file_id in _split_file_cell(value):
            columns.append((file_id, primer_from_column))
    return columns


def _resolve_metadata_path(value: str, base_dir: Path) -> Path:
    p = Path(value).expanduser()
    if not p.is_absolute():
        p = base_dir / p
    return p


def _metadata_key_variants(value: str | Path) -> list[str]:
    raw = str(value or '').strip()
    if not raw:
        return []
    path = Path(raw)
    candidates = [raw, path.name, path.stem]
    stripped = []
    for item in candidates:
        got = re.sub(r'([-_](f|r|forward|reverse))$', '', str(item).strip(), flags=re.IGNORECASE)
        stripped.append(got)
    seen = set()
    variants = []
    for item in candidates + stripped:
        if item and item not in seen:
            variants.append(item)
            seen.add(item)
    return variants


def _load_read_metadata(
    path: Optional[str | Path],
    primers: Iterable[str] = DEFAULT_PRIMERS,
) -> tuple[dict[str, dict[str, str]], list[Path]]:
    if not path:
        return {}, []
    metadata_path = Path(path).expanduser()
    base_dir = metadata_path.parent
    with _open_text(path) as handle:
        sample = handle.read(4096)
        handle.seek(0)
        delimiter = '\t' if str(path).endswith(('.tsv', '.tsv.gz', '.txt', '.txt.gz')) else ','
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=',\t')
            delimiter = dialect.delimiter
        except Exception:
            pass
        rows = [_normalise_row(row) for row in csv.DictReader(handle, delimiter=delimiter)]

    meta = {}
    listed_files: list[Path] = []
    for row in rows:
        sequence_id = _row_value(row, *SEQUENCE_ID_COLUMNS)
        row_primer = _row_value(row, 'primer', 'primer_name')
        row_direction = _row_value(row, 'direction', 'orientation', 'read', 'read_direction').lower()
        row_mode = _row_processing_mode(row)
        for file_id, column_primer in _metadata_file_columns(row, primers):
            resolved = _resolve_metadata_path(file_id, base_dir)
            inferred_sequence_id, inferred_primer, inferred_direction = infer_read_info(Path(file_id).name, primers)
            primer = row_primer or column_primer or inferred_primer
            if row_direction:
                direction = row_direction
            elif row_primer or column_primer:
                direction = 'reverse' if str(primer).upper().endswith('R') else 'forward'
            else:
                direction = inferred_direction
            if direction not in ('forward', 'reverse'):
                direction = 'reverse' if str(primer).upper().endswith('R') else 'forward'
            entry = {
                'sequence_id': sequence_id or inferred_sequence_id,
                'primer': primer,
                'direction': direction,
            }
            if row_mode:
                entry['processing_mode'] = row_mode
            listed_files.append(resolved)
            for key in _metadata_key_variants(file_id) + _metadata_key_variants(resolved):
                if key:
                    meta[key] = entry
    return meta, listed_files


def _collect_input_files(paths: Iterable[str | Path], *, recursive: bool = True) -> list[Path]:
    files: list[Path] = []
    for raw in paths:
        p = Path(raw)
        if p.is_dir():
            iterator = p.rglob('*') if recursive else p.glob('*')
            for child in iterator:
                if child.is_file() and _is_supported_read_file(child):
                    files.append(child)
        elif p.is_file() and _is_supported_read_file(p):
            files.append(p)
        else:
            logger.warning('[SANGER] Ignoring unsupported or missing input path: %s', raw)
    seen = set()
    unique = []
    for path in sorted(files):
        key = str(path)
        if key in seen:
            continue
        seen.add(key)
        unique.append(path)
    return unique


def _read_records_from_file(path: Path) -> list[tuple[str, str, list[int], bool, str]]:
    lower = path.name.lower()
    if lower.endswith(('.ab1', '.abi', '.ab1.gz', '.abi.gz')):
        seq, qual, quality_available = _read_ab1_with_quality_flag(path)
        read_id = Path(path.stem).stem if lower.endswith('.gz') else path.stem
        source = 'ab1_pcon' if quality_available else 'ab1_missing_pcon'
        return [(read_id, seq, qual, quality_available, source)]
    if lower.endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz')):
        return [(rid, seq, qual, True, 'fastq') for rid, seq, qual in read_fastq(path)]
    if lower.endswith(('.fasta', '.fa', '.fna', '.fasta.gz', '.fa.gz', '.fna.gz')):
        return [
            (h, s.upper(), [40] * len(s), False, 'fasta_assumed_q40')
            for h, s in read_fasta(str(path))
        ]
    raise ValueError(f'Unsupported input file: {path}')


def _identity_for_offset(seq_a: str, seq_b: str, offset: int) -> tuple[int, int, float]:
    start = max(0, offset)
    end = min(len(seq_a), offset + len(seq_b))
    overlap = max(0, end - start)
    if overlap <= 0:
        return 0, 0, 0.0
    matches = 0
    compared = 0
    for pos_a in range(start, end):
        pos_b = pos_a - offset
        a = seq_a[pos_a]
        b = seq_b[pos_b]
        if a == 'N' or b == 'N':
            continue
        compared += 1
        if a == b:
            matches += 1
    identity = matches / compared if compared else 0.0
    return overlap, compared, identity


def find_best_overlap(seq_a: str, seq_b: str, *, min_overlap: int = 40, min_identity: float = 0.85):
    best = None
    for offset in range(-len(seq_b) + min_overlap, len(seq_a) - min_overlap + 1):
        overlap, compared, identity = _identity_for_offset(seq_a, seq_b, offset)
        if overlap < min_overlap or compared < min_overlap:
            continue
        if identity < min_identity:
            continue
        candidate = (identity, overlap, compared, offset)
        if best is None or candidate > best:
            best = candidate
    if best is None:
        return None
    identity, overlap, compared, offset = best
    return {
        'offset': offset,
        'overlap': overlap,
        'compared': compared,
        'identity': identity,
    }


def merge_pair(
    seq_a: str,
    qual_a: list[int],
    seq_b: str,
    qual_b: list[int],
    *,
    min_overlap: int = 40,
    min_identity: float = 0.85,
    quality_difference: int = 10,
):
    overlap = find_best_overlap(seq_a, seq_b, min_overlap=min_overlap, min_identity=min_identity)
    if overlap is None:
        return None
    offset = overlap['offset']
    start = min(0, offset)
    end = max(len(seq_a), offset + len(seq_b))
    consensus = []
    consensus_qual = []
    conflicts = 0
    ambiguous_conflicts = 0
    for pos in range(start, end):
        idx_a = pos
        idx_b = pos - offset
        has_a = 0 <= idx_a < len(seq_a)
        has_b = 0 <= idx_b < len(seq_b)
        if has_a and has_b:
            base_a = seq_a[idx_a]
            base_b = seq_b[idx_b]
            q_a = qual_a[idx_a] if idx_a < len(qual_a) else 0
            q_b = qual_b[idx_b] if idx_b < len(qual_b) else 0
            c, ac = _append_consensus_base(
                consensus,
                consensus_qual,
                base_a,
                q_a,
                base_b,
                q_b,
                quality_difference=quality_difference,
            )
            conflicts += c
            ambiguous_conflicts += ac
        elif has_a:
            consensus.append(seq_a[idx_a])
            consensus_qual.append(qual_a[idx_a] if idx_a < len(qual_a) else 0)
        elif has_b:
            consensus.append(seq_b[idx_b])
            consensus_qual.append(qual_b[idx_b] if idx_b < len(qual_b) else 0)
    return ''.join(consensus), consensus_qual, {
        **overlap,
        'conflicts': conflicts,
        'ambiguous_conflicts': ambiguous_conflicts,
        'method': 'ungapped_overlap',
    }


def _append_consensus_base(
    consensus: list[str],
    consensus_qual: list[int],
    base_a: str,
    q_a: int,
    base_b: str,
    q_b: int,
    *,
    quality_difference: int = 10,
) -> tuple[int, int]:
    conflicts = 0
    ambiguous_conflicts = 0
    if base_a == base_b:
        consensus.append(base_a)
        consensus_qual.append(min(60, max(q_a, 0) + max(q_b, 0)))
    elif base_a == 'N':
        consensus.append(base_b)
        consensus_qual.append(q_b)
    elif base_b == 'N':
        consensus.append(base_a)
        consensus_qual.append(q_a)
    else:
        conflicts = 1
        if abs(q_a - q_b) >= quality_difference:
            if q_a > q_b:
                base, q = base_a, max(2, q_a - q_b)
            else:
                base, q = base_b, max(2, q_b - q_a)
        else:
            base, q = 'N', min(q_a, q_b)
            ambiguous_conflicts = 1
        consensus.append(base)
        consensus_qual.append(q)
    return conflicts, ambiguous_conflicts


def _merge_pair_gapped_one_way(
    seq_a: str,
    qual_a: list[int],
    seq_b: str,
    qual_b: list[int],
    *,
    min_overlap: int = 40,
    min_identity: float = 0.85,
    quality_difference: int = 10,
):
    matcher = difflib.SequenceMatcher(None, seq_a, seq_b, autojunk=False)
    blocks = [block for block in matcher.get_matching_blocks() if block.size]
    if not blocks:
        return None
    first = blocks[0]
    last = blocks[-1]
    matches = sum(block.size for block in blocks)
    span_a = (last.a + last.size) - first.a
    span_b = (last.b + last.size) - first.b
    overlap_length = max(span_a, span_b)
    identity = matches / overlap_length if overlap_length else 0.0
    if matches < min_overlap or identity < min_identity:
        return None

    consensus: list[str] = []
    consensus_qual: list[int] = []
    conflicts = 0
    ambiguous_conflicts = 0
    for tag, a0, a1, b0, b1 in matcher.get_opcodes():
        if tag == 'equal':
            for idx_a, idx_b in zip(range(a0, a1), range(b0, b1)):
                q_a = qual_a[idx_a] if idx_a < len(qual_a) else 0
                q_b = qual_b[idx_b] if idx_b < len(qual_b) else 0
                c, ac = _append_consensus_base(
                    consensus,
                    consensus_qual,
                    seq_a[idx_a],
                    q_a,
                    seq_b[idx_b],
                    q_b,
                    quality_difference=quality_difference,
                )
                conflicts += c
                ambiguous_conflicts += ac
        elif tag == 'replace':
            len_a = a1 - a0
            len_b = b1 - b0
            paired = min(len_a, len_b)
            for step in range(paired):
                idx_a = a0 + step
                idx_b = b0 + step
                q_a = qual_a[idx_a] if idx_a < len(qual_a) else 0
                q_b = qual_b[idx_b] if idx_b < len(qual_b) else 0
                c, ac = _append_consensus_base(
                    consensus,
                    consensus_qual,
                    seq_a[idx_a],
                    q_a,
                    seq_b[idx_b],
                    q_b,
                    quality_difference=quality_difference,
                )
                conflicts += c
                ambiguous_conflicts += ac
            for idx_a in range(a0 + paired, a1):
                consensus.append(seq_a[idx_a])
                consensus_qual.append(qual_a[idx_a] if idx_a < len(qual_a) else 0)
                conflicts += 1
            for idx_b in range(b0 + paired, b1):
                consensus.append(seq_b[idx_b])
                consensus_qual.append(qual_b[idx_b] if idx_b < len(qual_b) else 0)
                conflicts += 1
        elif tag == 'delete':
            for idx_a in range(a0, a1):
                consensus.append(seq_a[idx_a])
                consensus_qual.append(qual_a[idx_a] if idx_a < len(qual_a) else 0)
        elif tag == 'insert':
            for idx_b in range(b0, b1):
                consensus.append(seq_b[idx_b])
                consensus_qual.append(qual_b[idx_b] if idx_b < len(qual_b) else 0)

    return ''.join(consensus), consensus_qual, {
        'offset': 'gapped',
        'overlap': overlap_length,
        'compared': overlap_length,
        'identity': identity,
        'conflicts': conflicts,
        'ambiguous_conflicts': ambiguous_conflicts,
        'method': 'gapped_anchor',
    }


def merge_pair_gapped(
    seq_a: str,
    qual_a: list[int],
    seq_b: str,
    qual_b: list[int],
    *,
    min_overlap: int = 40,
    min_identity: float = 0.85,
    quality_difference: int = 10,
):
    first = _merge_pair_gapped_one_way(
        seq_a,
        qual_a,
        seq_b,
        qual_b,
        min_overlap=min_overlap,
        min_identity=min_identity,
        quality_difference=quality_difference,
    )
    second = _merge_pair_gapped_one_way(
        seq_b,
        qual_b,
        seq_a,
        qual_a,
        min_overlap=min_overlap,
        min_identity=min_identity,
        quality_difference=quality_difference,
    )
    if first is None:
        return second
    if second is None:
        return first
    first_stats = first[2]
    second_stats = second[2]
    first_key = (
        float(first_stats.get('identity', 0.0)),
        int(first_stats.get('overlap', 0) or 0),
        int(first_stats.get('compared', 0) or 0),
    )
    second_key = (
        float(second_stats.get('identity', 0.0)),
        int(second_stats.get('overlap', 0) or 0),
        int(second_stats.get('compared', 0) or 0),
    )
    return second if second_key > first_key else first


def assemble_read_group(
    reads: list[SangerRead],
    *,
    min_overlap: int = 40,
    min_identity: float = 0.85,
    quality_difference: int = 10,
):
    usable = [
        read for read in reads
        if read.status == 'kept' and read.oriented_sequence
    ]
    if not usable:
        return '', {
            'status': 'failed_no_reads',
            'read_count': len(reads),
            'used_reads': 0,
            'overlap_identity': 'NA',
            'overlap_length': 'NA',
            'conflicts': 'NA',
            'unmerged_reads': len(reads),
            'method': 'failed_no_reads',
            'qualities': [],
        }
    usable.sort(key=lambda r: (r.direction != 'forward', -len(r.oriented_sequence), r.read_id))
    consensus = usable[0].oriented_sequence
    qualities = usable[0].oriented_qualities or [40] * len(consensus)
    used = 1
    merge_stats = []
    unmerged = []
    for read in usable[1:]:
        got = merge_pair(
            consensus,
            qualities,
            read.oriented_sequence,
            read.oriented_qualities or [40] * len(read.oriented_sequence),
            min_overlap=min_overlap,
            min_identity=min_identity,
            quality_difference=quality_difference,
        )
        if got is None:
            got = merge_pair_gapped(
                consensus,
                qualities,
                read.oriented_sequence,
                read.oriented_qualities or [40] * len(read.oriented_sequence),
                min_overlap=min_overlap,
                min_identity=min_identity,
                quality_difference=quality_difference,
            )
        if got is None:
            unmerged.append(read.read_id)
            continue
        consensus, qualities, stats = got
        merge_stats.append(stats)
        used += 1
    if merge_stats:
        status = 'assembled'
        overlap_identity = min(stat['identity'] for stat in merge_stats)
        overlap_length = min(stat['overlap'] for stat in merge_stats)
        conflicts = sum(stat['conflicts'] for stat in merge_stats)
        ambiguous = sum(stat['ambiguous_conflicts'] for stat in merge_stats)
        method = ';'.join(sorted({str(stat.get('method', 'unknown')) for stat in merge_stats}))
    elif len(usable) == 1:
        status = 'single_read'
        overlap_identity = 'NA'
        overlap_length = 'NA'
        conflicts = 0
        ambiguous = 0
        method = 'single_read'
    else:
        status = 'partial_no_overlap'
        overlap_identity = 'NA'
        overlap_length = 'NA'
        conflicts = 0
        ambiguous = 0
        method = 'unmerged'
    return consensus, {
        'status': status,
        'read_count': len(reads),
        'used_reads': used,
        'overlap_identity': overlap_identity,
        'overlap_length': overlap_length,
        'conflicts': conflicts,
        'ambiguous_conflicts': ambiguous,
        'unmerged_reads': ';'.join(unmerged),
        'mean_quality': round(sum(qualities) / len(qualities), 2) if qualities else 0.0,
        'method': method,
        'qualities': qualities,
    }


def _resolve_group_processing_mode(
    reads: list[SangerRead],
    default_mode: str,
) -> tuple[str, str]:
    explicit_modes = sorted({read.processing_mode for read in reads if read.processing_mode_explicit})
    if not explicit_modes:
        return default_mode, ''
    if 'best_read' in explicit_modes:
        mode = 'best_read'
    elif 'assemble' in explicit_modes:
        mode = 'assemble'
    else:
        mode = default_mode
    warning = ''
    if len(explicit_modes) > 1:
        warning = 'conflicting_processing_modes_' + '_'.join(explicit_modes)
    return mode, warning


def _read_quality_key(read: SangerRead) -> tuple[float, float, int, int, str]:
    qualities = read.oriented_qualities or read.trimmed_qualities or []
    sequence = read.oriented_sequence or read.trimmed_sequence
    mean_q = _mean_quality(qualities)
    mean_error = _mean_error_probability(qualities)
    n_count = sequence.upper().count('N')
    return (mean_q, -mean_error, len(sequence), -n_count, read.read_id)


def select_best_read_group(reads: list[SangerRead]) -> tuple[str, dict[str, object]]:
    kept = [read for read in reads if read.status == 'kept' and read.oriented_sequence]
    if not kept:
        return '', {
            'status': 'failed_no_reads',
            'read_count': len(reads),
            'used_reads': 0,
            'overlap_identity': 'NA',
            'overlap_length': 'NA',
            'conflicts': 'NA',
            'ambiguous_conflicts': 'NA',
            'unmerged_reads': len(reads),
            'mean_quality': 0.0,
            'selected_read_id': '',
            'method': 'best_read',
            'qualities': [],
        }
    best = max(kept, key=_read_quality_key)
    qualities = best.oriented_qualities or []
    unselected = [read.read_id for read in kept if read.read_id != best.read_id]
    return best.oriented_sequence, {
        'status': 'best_read',
        'read_count': len(reads),
        'used_reads': 1,
        'overlap_identity': 'NA',
        'overlap_length': 'NA',
        'conflicts': 'NA',
        'ambiguous_conflicts': 'NA',
        'unmerged_reads': ';'.join(unselected),
        'mean_quality': round(_mean_quality(qualities), 2),
        'selected_read_id': best.read_id,
        'method': 'best_read',
        'qualities': qualities,
    }


def _classify_read_qc(
    read: SangerRead,
    *,
    read_min_length: int,
    min_mean_quality: float,
    max_read_expected_errors: float,
    max_n_percent: float,
    max_internal_low_quality_run: int,
    allow_missing_quality: bool,
) -> None:
    fail_reasons: list[str] = []
    warn_reasons: list[str] = []
    sequence = read.masked_sequence or read.trimmed_sequence
    qualities = read.masked_qualities or read.trimmed_qualities or []
    expected_errors = _expected_errors(qualities, sequence)
    mean_q = _mean_quality(qualities)
    n_pct = _n_percent(sequence)

    if read.quality_source == 'ab1_missing_pcon' and not allow_missing_quality:
        fail_reasons.append('quality_scores_missing_from_ab1')
    elif read.quality_source == 'ab1_missing_pcon':
        warn_reasons.append('quality_scores_missing_from_ab1_allowed')
    elif read.quality_source == 'fasta_assumed_q40':
        warn_reasons.append('quality_scores_unavailable_assumed_q40')

    if len(sequence) < read_min_length:
        fail_reasons.append(f'trimmed_length_lt_{read_min_length}')
    if sequence and mean_q < min_mean_quality:
        fail_reasons.append(f'mean_q_lt_{min_mean_quality:g}')
    if expected_errors > max_read_expected_errors:
        fail_reasons.append(f'expected_errors_gt_{max_read_expected_errors:g}')
    if n_pct > max_n_percent:
        fail_reasons.append(f'n_percent_gt_{max_n_percent:g}')
    if read.longest_low_quality_run > max_internal_low_quality_run:
        fail_reasons.append(f'internal_low_quality_run_gt_{max_internal_low_quality_run}')
    elif read.longest_low_quality_run > 0:
        warn_reasons.append(f'internal_low_quality_masked_run_{read.longest_low_quality_run}')
    if read.masked_bases:
        warn_reasons.append(f'masked_low_quality_bases_{read.masked_bases}')

    read.qc_class = _qc_class_from_reasons(fail_reasons, warn_reasons)
    read.qc_reasons = fail_reasons + warn_reasons
    if read.qc_class == 'FAIL_QC':
        read.status = 'filtered'
    else:
        read.status = 'kept'
    read.warning = _join_reasons(read.qc_reasons)


def _classify_output_qc(
    consensus: str,
    qualities: list[int],
    stats: dict[str, object],
    group: list[SangerRead],
    *,
    final_min_length: int,
    min_mean_quality: float,
    max_output_expected_errors: float,
    max_n_percent: float,
    max_conflict_density: float,
    processing_mode: str,
    mode_warning: str,
) -> dict[str, object]:
    fail_reasons: list[str] = []
    warn_reasons: list[str] = []
    output_length = len(consensus or '')
    mean_q = _mean_quality(qualities)
    expected_errors = _expected_errors(qualities, consensus)
    n_pct = _n_percent(consensus)
    conflicts_raw = stats.get('conflicts', 0)
    ambiguous_raw = stats.get('ambiguous_conflicts', 0)
    try:
        conflicts = int(conflicts_raw)
    except Exception:
        conflicts = 0
    try:
        ambiguous = int(ambiguous_raw)
    except Exception:
        ambiguous = 0
    conflict_density = (conflicts / output_length) * 100.0 if output_length else 0.0

    if not consensus:
        fail_reasons.append('no_output_sequence')
    if output_length and output_length < final_min_length:
        fail_reasons.append(f'output_length_lt_{final_min_length}')
    if consensus and mean_q < min_mean_quality:
        fail_reasons.append(f'output_mean_q_lt_{min_mean_quality:g}')
    if expected_errors > max_output_expected_errors:
        fail_reasons.append(f'output_expected_errors_gt_{max_output_expected_errors:g}')
    if n_pct > max_n_percent:
        fail_reasons.append(f'output_n_percent_gt_{max_n_percent:g}')
    if conflict_density > max_conflict_density:
        fail_reasons.append(f'overlap_conflict_density_gt_{max_conflict_density:g}')

    status = str(stats.get('status', 'unknown'))
    if status in {'failed_no_reads', 'partial_no_overlap'}:
        fail_reasons.append(status)
    elif status in {'single_read'} and len(group) > 1 and processing_mode == 'assemble':
        warn_reasons.append('only_one_primer_read_passed_qc')
    if ambiguous:
        warn_reasons.append(f'ambiguous_overlap_conflicts_{ambiguous}')
    if mode_warning:
        warn_reasons.append(mode_warning)
    failed_reads = [read.read_id for read in group if read.qc_class == 'FAIL_QC']
    if failed_reads:
        warn_reasons.append('failed_reads_' + str(len(failed_reads)))
    if any(read.qc_class == 'PASS_WITH_WARNINGS' for read in group):
        warn_reasons.append('read_level_warnings')
    if 'gapped_anchor' in str(stats.get('method', '')):
        warn_reasons.append('gapped_overlap_fallback_used')

    qc_class = _qc_class_from_reasons(fail_reasons, warn_reasons)
    recommendation = _recommendation_from_qc_class(qc_class)
    if recommendation == 'RESEQUENCE':
        suggested_action = 'Repeat Sanger sequencing for this isolate; inspect chromatograms and primer setup before reuse.'
    elif recommendation == 'MANUAL_REVIEW':
        suggested_action = 'Inspect chromatograms/overlap manually; accept only if trace evidence supports the consensus.'
    else:
        suggested_action = 'Accept for downstream BranchManager evaluate.'
    reasons = fail_reasons + warn_reasons
    return {
        'qc_class': qc_class,
        'recommendation': recommendation,
        'reasons': reasons,
        'suggested_action': suggested_action,
        'mean_quality': round(mean_q, 2),
        'expected_errors': round(expected_errors, 4),
        'n_percent': round(n_pct, 3),
        'conflict_density': round(conflict_density, 3),
        'failed_read_ids': ';'.join(failed_reads),
    }


def _clean_field(value: object) -> str:
    return str(value if value is not None else '').replace('\t', ' ').replace('\r', ' ').replace('\n', ' ')


def _mean_quality(qualities: list[int] | None) -> float:
    return sum(qualities) / len(qualities) if qualities else 0.0


def _mean_error_probability(qualities: list[int] | None) -> float:
    return sum(10 ** (-q / 10) for q in qualities) / len(qualities) if qualities else 0.0


def _error_probability(q: int) -> float:
    return 10 ** (-max(0, int(q)) / 10)


def _expected_errors(qualities: list[int] | None, sequence: str | None = None) -> float:
    if not qualities:
        return 0.0
    total = 0.0
    for idx, q in enumerate(qualities):
        if sequence is not None and idx < len(sequence) and sequence[idx].upper() not in {'A', 'C', 'G', 'T'}:
            continue
        total += _error_probability(q)
    return total


def _n_percent(sequence: str) -> float:
    if not sequence:
        return 0.0
    return (sequence.upper().count('N') / len(sequence)) * 100.0


def _longest_true_run(values: Iterable[bool]) -> int:
    longest = 0
    current = 0
    for value in values:
        if value:
            current += 1
            longest = max(longest, current)
        else:
            current = 0
    return longest


def mask_low_quality_bases(
    sequence: str,
    qualities: list[int],
    *,
    mask_quality: int = 20,
) -> tuple[str, list[int], int, int]:
    masked: list[str] = []
    out_qualities: list[int] = []
    masked_count = 0
    low_quality_flags: list[bool] = []
    for idx, base in enumerate(str(sequence or '').upper()):
        q = qualities[idx] if idx < len(qualities) else 0
        low_quality = q < mask_quality
        ambiguous = base not in {'A', 'C', 'G', 'T'}
        low_quality_flags.append(low_quality or ambiguous)
        if low_quality or ambiguous:
            masked.append('N')
            out_qualities.append(min(q, 2))
            if base != 'N' or low_quality:
                masked_count += 1
        else:
            masked.append(base)
            out_qualities.append(q)
    return ''.join(masked), out_qualities, masked_count, _longest_true_run(low_quality_flags)


def _join_reasons(reasons: Iterable[str]) -> str:
    return ';'.join(str(reason) for reason in reasons if str(reason))


def _qc_class_from_reasons(fail_reasons: list[str], warn_reasons: list[str]) -> str:
    if fail_reasons:
        return 'FAIL_QC'
    if warn_reasons:
        return 'PASS_WITH_WARNINGS'
    return 'PASS_HIGH_CONFIDENCE'


def _recommendation_from_qc_class(qc_class: str) -> str:
    if qc_class == 'FAIL_QC':
        return 'RESEQUENCE'
    if qc_class == 'PASS_WITH_WARNINGS':
        return 'MANUAL_REVIEW'
    return 'ACCEPT'


def _write_per_base_error_tsv(path: Path, reads: list[SangerRead]) -> None:
    with open(path, 'w') as handle:
        handle.write(
            'ReadID\tSequenceID\tPrimer\tDirection\tSourceFile\tPosition\tBase\tQuality\t'
            'ErrorProbability\tTrimRegion\tMaskedLowQuality\tBaseAfterMask\t'
            'KeptAfterTrim\tOrientedPosition\n'
        )
        for read in reads:
            for idx, base in enumerate(read.raw_sequence):
                q = read.raw_qualities[idx] if idx < len(read.raw_qualities) else 0
                masked_low_quality = 'no'
                base_after_mask = ''
                if idx < read.trim_start:
                    trim_region = 'left_trimmed'
                elif idx >= read.trim_end:
                    trim_region = 'right_trimmed'
                else:
                    trim_region = 'kept'
                    trim_idx = idx - read.trim_start
                    if trim_idx < len(read.masked_sequence):
                        base_after_mask = read.masked_sequence[trim_idx]
                        if base_after_mask == 'N' and base.upper() != 'N':
                            masked_low_quality = 'yes'
                kept = read.status == 'kept' and trim_region == 'kept'
                if kept and read.direction == 'reverse':
                    oriented_position = read.trim_end - idx
                elif kept:
                    oriented_position = idx - read.trim_start + 1
                else:
                    oriented_position = ''
                handle.write(
                    f'{_clean_field(read.read_id)}\t{_clean_field(read.sequence_id)}\t'
                    f'{_clean_field(read.primer)}\t{_clean_field(read.direction)}\t'
                    f'{_clean_field(read.source_file)}\t{idx + 1}\t{base}\t{q}\t'
                    f'{_error_probability(q):.6g}\t{trim_region}\t'
                    f'{masked_low_quality}\t{base_after_mask}\t'
                    f'{"yes" if kept else "no"}\t{oriented_position}\n'
                )


def _quality_points(qualities: list[int], x0: float, y0: float, width: float, height: float) -> str:
    if not qualities:
        return ''
    max_points = 700
    step = max(1, len(qualities) // max_points)
    binned: list[tuple[int, float]] = []
    for start in range(0, len(qualities), step):
        chunk = qualities[start:start + step]
        binned.append((start + len(chunk) // 2, _mean_quality(chunk)))
    denom = max(1, len(qualities) - 1)
    points = []
    for pos, q in binned:
        x = x0 + (pos / denom) * width
        y = y0 + height - (min(45.0, max(0.0, q)) / 45.0) * height
        points.append(f'{x:.1f},{y:.1f}')
    return ' '.join(points)


def _write_read_error_profile_svg(path: Path, reads: list[SangerRead]) -> None:
    width = 1180
    panel_height = 112
    left = 190
    right = 38
    graph_width = width - left - right
    height = 72 + max(1, len(reads)) * panel_height
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<text x="24" y="32" font-family="Arial, sans-serif" font-size="20" font-weight="700" fill="#1f2937">Sanger read error and trimming profiles</text>',
        '<text x="24" y="54" font-family="Arial, sans-serif" font-size="12" fill="#4b5563">Blue line shows Phred quality; lower quality means higher expected base-call error. Green band is the retained trim window.</text>',
    ]
    if not reads:
        lines.append('<text x="24" y="92" font-family="Arial, sans-serif" font-size="14" fill="#6b7280">No reads parsed.</text>')
    for idx, read in enumerate(reads):
        top = 76 + idx * panel_height
        graph_top = top + 32
        graph_height = 48
        raw_len = max(1, len(read.raw_sequence))
        trim_start_x = left + (read.trim_start / raw_len) * graph_width
        trim_end_x = left + (read.trim_end / raw_len) * graph_width
        points = _quality_points(read.raw_qualities, left, graph_top, graph_width, graph_height)
        raw_mean = _mean_quality(read.raw_qualities)
        trim_mean = _mean_quality(read.trimmed_qualities or [])
        masked_mean = _mean_quality(read.masked_qualities or [])
        status_colour = '#047857' if read.qc_class == 'PASS_HIGH_CONFIDENCE' else ('#b45309' if read.status == 'kept' else '#b91c1c')
        q20_y = graph_top + graph_height - (20 / 45.0) * graph_height
        q30_y = graph_top + graph_height - (30 / 45.0) * graph_height
        lines.extend([
            f'<rect x="18" y="{top - 12}" width="{width - 36}" height="{panel_height - 12}" rx="6" fill="#f9fafb" stroke="#e5e7eb"/>',
            f'<text x="30" y="{top + 10}" font-family="Arial, sans-serif" font-size="13" font-weight="700" fill="#111827">{html_escape(read.read_id)}</text>',
            f'<text x="30" y="{top + 30}" font-family="Arial, sans-serif" font-size="11" fill="#4b5563">isolate {html_escape(read.sequence_id)} | {html_escape(read.primer)} | {html_escape(read.direction)} | mode {html_escape(read.processing_mode)}</text>',
            f'<text x="30" y="{top + 50}" font-family="Arial, sans-serif" font-size="11" fill="{status_colour}">{html_escape(read.status)} | {html_escape(read.qc_class)} | {html_escape(read.warning)}</text>',
            f'<text x="30" y="{top + 70}" font-family="Arial, sans-serif" font-size="11" fill="#4b5563">raw {len(read.raw_sequence)} bp, trimmed {len(read.trimmed_sequence)} bp, masked {read.masked_bases} bp, mean Q {raw_mean:.1f} -> {trim_mean:.1f} -> {masked_mean:.1f}</text>',
            f'<line x1="{left}" y1="{graph_top + graph_height}" x2="{left + graph_width}" y2="{graph_top + graph_height}" stroke="#d1d5db"/>',
            f'<line x1="{left}" y1="{graph_top}" x2="{left + graph_width}" y2="{graph_top}" stroke="#e5e7eb"/>',
            f'<line x1="{left}" y1="{q20_y:.1f}" x2="{left + graph_width}" y2="{q20_y:.1f}" stroke="#f59e0b" stroke-width="0.8" stroke-dasharray="4 3"/>',
            f'<line x1="{left}" y1="{q30_y:.1f}" x2="{left + graph_width}" y2="{q30_y:.1f}" stroke="#10b981" stroke-width="0.8" stroke-dasharray="4 3"/>',
            f'<rect x="{trim_start_x:.1f}" y="{graph_top}" width="{max(0.0, trim_end_x - trim_start_x):.1f}" height="{graph_height}" fill="#bbf7d0" opacity="0.55"/>',
        ])
        mask_start = None
        for offset, masked_base in enumerate(read.masked_sequence + 'A'):
            original_base = read.trimmed_sequence[offset] if offset < len(read.trimmed_sequence) else 'A'
            is_masked = offset < len(read.masked_sequence) and masked_base == 'N' and original_base != 'N'
            if is_masked and mask_start is None:
                mask_start = offset
            elif not is_masked and mask_start is not None:
                seg_start = read.trim_start + mask_start
                seg_end = read.trim_start + offset
                x0 = left + (seg_start / raw_len) * graph_width
                x1 = left + (seg_end / raw_len) * graph_width
                lines.append(
                    f'<rect x="{x0:.1f}" y="{graph_top}" width="{max(1.0, x1 - x0):.1f}" height="{graph_height}" fill="#f97316" opacity="0.35"/>'
                )
                mask_start = None
        if points:
            lines.append(f'<polyline points="{points}" fill="none" stroke="#2563eb" stroke-width="1.5"/>')
        lines.extend([
            f'<text x="{left}" y="{graph_top + graph_height + 18}" font-family="Arial, sans-serif" font-size="10" fill="#6b7280">1</text>',
            f'<text x="{left + graph_width - 44}" y="{graph_top + graph_height + 18}" font-family="Arial, sans-serif" font-size="10" fill="#6b7280">{len(read.raw_sequence)} bases</text>',
        ])
    lines.append('</svg>')
    path.write_text('\n'.join(lines))


def _write_assembly_overview_svg(
    path: Path,
    groups: dict[str, list[SangerRead]],
    assembly_rows: list[dict[str, object]],
) -> None:
    width = 1180
    left = 260
    right = 40
    graph_width = width - left - right
    row_lookup = {str(row.get('SequenceID', '')): row for row in assembly_rows}
    group_heights = {sid: 82 + max(1, len(reads)) * 22 for sid, reads in groups.items()}
    height = 74 + sum(group_heights.get(sid, 100) for sid in sorted(groups))
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<text x="24" y="32" font-family="Arial, sans-serif" font-size="20" font-weight="700" fill="#1f2937">Sanger isolate assembly overview</text>',
        '<text x="24" y="54" font-family="Arial, sans-serif" font-size="12" fill="#4b5563">Consensus length, read contribution, overlap identity, conflicts, and unmerged reads per isolate.</text>',
    ]
    y = 82
    if not groups:
        lines.append('<text x="24" y="92" font-family="Arial, sans-serif" font-size="14" fill="#6b7280">No read groups parsed.</text>')
    for sequence_id in sorted(groups):
        reads = sorted(groups[sequence_id], key=lambda r: (r.status != 'kept', r.direction != 'forward', r.read_id))
        row = row_lookup.get(sequence_id, {})
        group_height = group_heights.get(sequence_id, 100)
        output_length = int(row.get('OutputLength') or 0)
        max_len = max([output_length] + [len(r.oriented_sequence or r.trimmed_sequence or r.raw_sequence) for r in reads] + [1])
        status = str(row.get('Status', 'unknown'))
        qc_class = str(row.get('QCClass', ''))
        if qc_class == 'PASS_HIGH_CONFIDENCE':
            status_colour = '#047857'
        elif qc_class == 'FAIL_QC' or status.startswith('filtered') or status.startswith('failed'):
            status_colour = '#b91c1c'
        else:
            status_colour = '#b45309'
        mode = str(row.get('ProcessingMode', 'assemble'))
        reasons = str(row.get('Reasons', ''))
        selected = str(row.get('SelectedReadID', ''))
        mode_warning = str(row.get('ModeWarning', ''))
        extra_bits = []
        if selected:
            extra_bits.append(f'selected {selected}')
        if reasons:
            extra_bits.append(reasons)
        if mode_warning:
            extra_bits.append(mode_warning)
        extra = '; ' + '; '.join(extra_bits) if extra_bits else ''
        details = (
            f'{status}; {qc_class}; {row.get("Recommendation", "")}; mode {mode}; output {output_length} bp; reads {row.get("ReadCount", len(reads))}; '
            f'used {row.get("UsedReads", 0)}; overlap {row.get("OverlapIdentity", "NA")} / {row.get("OverlapLength", "NA")}; '
            f'conflicts {row.get("Conflicts", "NA")}; unmerged {row.get("UnmergedReads", "")}; '
            f'method {row.get("MergeMethod", "NA")}{extra}'
        )
        lines.extend([
            f'<rect x="18" y="{y - 20}" width="{width - 36}" height="{group_height - 10}" rx="6" fill="#f9fafb" stroke="#e5e7eb"/>',
            f'<text x="30" y="{y + 2}" font-family="Arial, sans-serif" font-size="14" font-weight="700" fill="#111827">{html_escape(sequence_id)}</text>',
            f'<text x="30" y="{y + 22}" font-family="Arial, sans-serif" font-size="11" fill="{status_colour}">{html_escape(details)}</text>',
        ])
        consensus_width = (output_length / max_len) * graph_width if output_length else 0
        lines.extend([
            f'<text x="30" y="{y + 48}" font-family="Arial, sans-serif" font-size="11" fill="#4b5563">consensus</text>',
            f'<rect x="{left}" y="{y + 38}" width="{max(1.0, consensus_width):.1f}" height="10" rx="2" fill="#111827" opacity="0.82"/>',
        ])
        read_y = y + 66
        for read in reads:
            read_len = len(read.oriented_sequence or read.trimmed_sequence or read.raw_sequence)
            bar_width = (read_len / max_len) * graph_width
            colour = '#2563eb' if read.status == 'kept' else '#b91c1c'
            label = f'{read.read_id} | {read.primer} | {read.direction} | {read.status} | {read_len} bp'
            lines.extend([
                f'<text x="30" y="{read_y + 9}" font-family="Arial, sans-serif" font-size="10" fill="#4b5563">{html_escape(label)}</text>',
                f'<rect x="{left}" y="{read_y}" width="{max(1.0, bar_width):.1f}" height="10" rx="2" fill="{colour}" opacity="0.72"/>',
            ])
            read_y += 22
        y += group_height
    lines.append('</svg>')
    path.write_text('\n'.join(lines))


def run_sanger(
    inputs: Iterable[str | Path],
    outdir: str | Path,
    *,
    read_metadata: Optional[str | Path] = None,
    sample_map: Optional[str | Path] = None,
    primers: Iterable[str] = DEFAULT_PRIMERS,
    min_quality: int = 20,
    window: int = 20,
    min_length: int = 800,
    min_read_length: Optional[int] = None,
    min_mean_quality: float = 25.0,
    mask_quality: int = 20,
    max_read_expected_errors: float = 8.0,
    max_output_expected_errors: float = 5.0,
    max_n_percent: float = 1.0,
    max_internal_low_quality_run: int = 20,
    max_conflict_density: float = 1.0,
    quality_difference: int = 10,
    allow_missing_quality: bool = False,
    min_overlap: int = 40,
    min_overlap_identity: float = 0.85,
    assemble: bool = True,
    recursive: bool = True,
) -> dict[str, str]:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    final_min_length = int(min_length or 800)
    read_min_length = int(min_read_length if min_read_length is not None else final_min_length)
    default_mode = 'assemble' if assemble else 'best_read'
    metadata, listed_files = _load_read_metadata(read_metadata, primers)
    sample_metadata, sample_files = _load_read_metadata(sample_map, primers)
    metadata.update(sample_metadata)
    listed_files.extend(sample_files)
    files = _collect_input_files(inputs or [], recursive=recursive)
    seen_files = {str(path) for path in files}
    for path in listed_files:
        if path.is_file() and _is_supported_read_file(path) and str(path) not in seen_files:
            files.append(path)
            seen_files.add(str(path))
    files = sorted(files)
    if not files:
        raise ValueError('No AB1/FASTA/FASTQ input files found; provide --input and/or --sample-map')

    reads: list[SangerRead] = []
    raw_records = []
    trimmed_records = []

    for path in files:
        for record_id, sequence, qualities, quality_available, quality_source in _read_records_from_file(path):
            stem_prefix = path.stem.split('_', 1)[0]
            record_prefix = str(record_id).split('_', 1)[0]
            key_candidates = [str(path), path.name, path.stem, record_id, stem_prefix, record_prefix]
            meta = {}
            for key in key_candidates:
                if key in metadata:
                    meta = metadata[key]
                    break
            inferred_seq_id, inferred_primer, inferred_direction = infer_read_info(record_id if record_id != path.stem else path.name, primers)
            sequence_id = meta.get('sequence_id') or inferred_seq_id
            primer = meta.get('primer') or inferred_primer
            direction = meta.get('direction') or inferred_direction
            if direction not in ('forward', 'reverse'):
                direction = 'reverse' if str(primer).upper().endswith('R') else 'forward'
            meta_mode = normalise_processing_mode(meta.get('processing_mode'), default=None)
            read = SangerRead(
                read_id=record_id,
                sequence_id=sequence_id,
                primer=primer,
                direction=direction,
                source_file=str(path),
                raw_sequence=sequence,
                raw_qualities=qualities,
                quality_available=quality_available,
                quality_source=quality_source,
                processing_mode=meta_mode or default_mode,
                processing_mode_explicit=bool(meta_mode),
            )
            raw_records.append((read.read_id, read.raw_sequence))
            trimmed, trimmed_q, start, end = trim_by_quality(
                sequence,
                qualities,
                min_quality=min_quality,
                window=window,
            )
            read.trimmed_sequence = trimmed
            read.trimmed_qualities = trimmed_q
            read.trim_start = start
            read.trim_end = end
            masked, masked_q, masked_count, lowq_run = mask_low_quality_bases(
                trimmed,
                trimmed_q,
                mask_quality=mask_quality,
            )
            read.masked_sequence = masked
            read.masked_qualities = masked_q
            read.masked_bases = masked_count
            read.longest_low_quality_run = lowq_run
            _classify_read_qc(
                read,
                read_min_length=read_min_length,
                min_mean_quality=float(min_mean_quality),
                max_read_expected_errors=float(max_read_expected_errors),
                max_n_percent=float(max_n_percent),
                max_internal_low_quality_run=int(max_internal_low_quality_run),
                allow_missing_quality=bool(allow_missing_quality),
            )
            if read.status == 'kept':
                if direction == 'reverse':
                    read.oriented_sequence = reverse_complement(masked)
                    read.oriented_qualities = list(reversed(masked_q))
                else:
                    read.oriented_sequence = masked
                    read.oriented_qualities = masked_q
                trimmed_records.append((read.read_id, read.oriented_sequence))
            reads.append(read)

    raw_fasta = out / 'raw_reads.fasta'
    trimmed_fasta = out / 'trimmed_oriented_reads.fasta'
    assembled_fasta = out / 'assembled.fasta'
    read_qc_tsv = out / 'read_qc.tsv'
    per_base_error_tsv = out / 'per_base_error.tsv'
    assembly_tsv = out / 'assembly_report.tsv'
    recommendations_tsv = out / 'resequence_recommendations.tsv'
    qc_policy_tsv = out / 'sanger_qc_policy.tsv'
    read_error_svg = out / 'read_error_profiles.svg'
    assembly_svg = out / 'assembly_overview.svg'
    summary_txt = out / 'sanger_summary.txt'

    write_fasta(raw_records, str(raw_fasta))
    write_fasta(trimmed_records, str(trimmed_fasta))
    _write_per_base_error_tsv(per_base_error_tsv, reads)
    _write_read_error_profile_svg(read_error_svg, reads)

    with open(read_qc_tsv, 'w') as handle:
        handle.write(
            'ReadID\tSequenceID\tPrimer\tDirection\tSourceFile\tRawLength\tTrimmedLength\t'
            'MaskedLength\tTrimStart\tTrimEnd\tLeftTrimmedBases\tRightTrimmedBases\t'
            'MeanRawQuality\tMeanTrimmedQuality\tMeanMaskedQuality\t'
            'RawExpectedErrors\tTrimmedExpectedErrors\tMaskedExpectedErrors\t'
            'MeanRawErrorProbability\tMeanTrimmedErrorProbability\tMaskedNPercent\t'
            'MaskedBases\tLongestLowQualityRun\tQualityAvailable\tQualitySource\t'
            'Status\tQCClass\tReasons\tProcessingMode\tProcessingModeExplicit\n'
        )
        for read in reads:
            raw_mean = _mean_quality(read.raw_qualities)
            trimmed_q = read.trimmed_qualities or []
            trim_mean = _mean_quality(trimmed_q)
            masked_q = read.masked_qualities or []
            masked_mean = _mean_quality(masked_q)
            left_trimmed = read.trim_start
            right_trimmed = max(0, len(read.raw_sequence) - read.trim_end)
            handle.write(
                f'{_clean_field(read.read_id)}\t{_clean_field(read.sequence_id)}\t'
                f'{_clean_field(read.primer)}\t{_clean_field(read.direction)}\t{_clean_field(read.source_file)}\t'
                f'{len(read.raw_sequence)}\t{len(read.trimmed_sequence)}\t{len(read.masked_sequence)}\t'
                f'{read.trim_start}\t{read.trim_end}\t{left_trimmed}\t{right_trimmed}\t'
                f'{raw_mean:.2f}\t{trim_mean:.2f}\t{masked_mean:.2f}\t'
                f'{_expected_errors(read.raw_qualities, read.raw_sequence):.4f}\t'
                f'{_expected_errors(trimmed_q, read.trimmed_sequence):.4f}\t'
                f'{_expected_errors(masked_q, read.masked_sequence):.4f}\t'
                f'{_mean_error_probability(read.raw_qualities):.6g}\t{_mean_error_probability(trimmed_q):.6g}\t'
                f'{_n_percent(read.masked_sequence):.3f}\t{read.masked_bases}\t'
                f'{read.longest_low_quality_run}\t{"yes" if read.quality_available else "no"}\t'
                f'{_clean_field(read.quality_source)}\t{read.status}\t{read.qc_class}\t'
                f'{_clean_field(read.warning)}\t{read.processing_mode}\t'
                f'{"yes" if read.processing_mode_explicit else "no"}\n'
            )

    groups: dict[str, list[SangerRead]] = {}
    for read in reads:
        groups.setdefault(read.sequence_id, []).append(read)

    assembled_records = []
    assembly_rows: list[dict[str, object]] = []
    recommendation_rows: list[dict[str, object]] = []
    with open(assembly_tsv, 'w') as handle:
        handle.write(
            'SequenceID\tStatus\tReadCount\tUsedReads\tOutputLength\tMeanQuality\t'
            'OutputExpectedErrors\tOutputNPercent\tConflictDensity\tQCClass\tRecommendation\t'
            'OverlapIdentity\tOverlapLength\tConflicts\tAmbiguousConflicts\tUnmergedReads\t'
            'ReadIDs\tKeptReadIDs\tMergeMethod\tProcessingMode\tSelectedReadID\t'
            'PassLengthQC\tReasons\tSuggestedAction\tFailedReadIDs\tModeWarning\n'
        )
        for sequence_id, group in sorted(groups.items()):
            processing_mode, mode_warning = _resolve_group_processing_mode(group, default_mode)
            if processing_mode == 'assemble':
                consensus, stats = assemble_read_group(
                    group,
                    min_overlap=min_overlap,
                    min_identity=min_overlap_identity,
                    quality_difference=quality_difference,
                )
            else:
                consensus, stats = select_best_read_group(group)
            consensus_qualities = list(stats.get('qualities') or [])
            final_qc = _classify_output_qc(
                consensus,
                consensus_qualities,
                stats,
                group,
                final_min_length=final_min_length,
                min_mean_quality=float(min_mean_quality),
                max_output_expected_errors=float(max_output_expected_errors),
                max_n_percent=float(max_n_percent),
                max_conflict_density=float(max_conflict_density),
                processing_mode=processing_mode,
                mode_warning=mode_warning,
            )
            pass_length_qc = 'no' if any(
                str(reason).startswith('output_length_lt_')
                for reason in final_qc.get('reasons', [])
            ) or not consensus else 'yes'
            if consensus and pass_length_qc == 'no':
                stats['status'] = f'filtered_output_length_lt_{final_min_length}'
            if consensus and final_qc['qc_class'] != 'FAIL_QC':
                assembled_records.append((sequence_id, consensus))
            overlap_identity = stats['overlap_identity']
            if isinstance(overlap_identity, float):
                overlap_identity = f'{overlap_identity:.4f}'
            read_ids = ';'.join(_clean_field(read.read_id) for read in group)
            kept_read_ids = ';'.join(_clean_field(read.read_id) for read in group if read.status == 'kept')
            reasons = _join_reasons(final_qc.get('reasons', []))
            row = {
                'SequenceID': sequence_id,
                'Status': stats['status'],
                'ReadCount': stats['read_count'],
                'UsedReads': stats['used_reads'],
                'OutputLength': len(consensus),
                'MeanQuality': final_qc.get('mean_quality', stats.get('mean_quality', 0.0)),
                'OutputExpectedErrors': final_qc.get('expected_errors', 0.0),
                'OutputNPercent': final_qc.get('n_percent', 0.0),
                'ConflictDensity': final_qc.get('conflict_density', 0.0),
                'QCClass': final_qc.get('qc_class', 'FAIL_QC'),
                'Recommendation': final_qc.get('recommendation', 'RESEQUENCE'),
                'OverlapIdentity': overlap_identity,
                'OverlapLength': stats['overlap_length'],
                'Conflicts': stats['conflicts'],
                'AmbiguousConflicts': stats.get('ambiguous_conflicts', 'NA'),
                'UnmergedReads': stats['unmerged_reads'],
                'ReadIDs': read_ids,
                'KeptReadIDs': kept_read_ids,
                'MergeMethod': stats.get('method', 'NA'),
                'ProcessingMode': processing_mode,
                'SelectedReadID': stats.get('selected_read_id', ''),
                'PassLengthQC': pass_length_qc if consensus else 'no',
                'Reasons': reasons,
                'SuggestedAction': final_qc.get('suggested_action', ''),
                'FailedReadIDs': final_qc.get('failed_read_ids', ''),
                'ModeWarning': mode_warning,
            }
            assembly_rows.append(row)
            recommendation_rows.append(row)
            handle.write(
                f'{_clean_field(sequence_id)}\t{stats["status"]}\t{stats["read_count"]}\t{stats["used_reads"]}\t'
                f'{len(consensus)}\t{row["MeanQuality"]}\t{row["OutputExpectedErrors"]}\t'
                f'{row["OutputNPercent"]}\t{row["ConflictDensity"]}\t{row["QCClass"]}\t'
                f'{row["Recommendation"]}\t{overlap_identity}\t'
                f'{stats["overlap_length"]}\t{stats["conflicts"]}\t{stats.get("ambiguous_conflicts", "NA")}\t'
                f'{_clean_field(stats["unmerged_reads"])}\t{read_ids}\t{kept_read_ids}\t'
                f'{_clean_field(stats.get("method", "NA"))}\t{processing_mode}\t'
                f'{_clean_field(stats.get("selected_read_id", ""))}\t'
                f'{pass_length_qc if consensus else "no"}\t{_clean_field(reasons)}\t'
                f'{_clean_field(row["SuggestedAction"])}\t{_clean_field(row["FailedReadIDs"])}\t'
                f'{_clean_field(mode_warning)}\n'
            )

    write_fasta(assembled_records, str(assembled_fasta))
    _write_assembly_overview_svg(assembly_svg, groups, assembly_rows)
    with open(recommendations_tsv, 'w') as handle:
        handle.write(
            'SequenceID\tRecommendation\tQCClass\tReasons\tSuggestedAction\t'
            'Status\tProcessingMode\tOutputLength\tMeanQuality\tOutputExpectedErrors\t'
            'OutputNPercent\tConflictDensity\tReadCount\tUsedReads\tReadIDs\t'
            'KeptReadIDs\tFailedReadIDs\tSelectedReadID\tMergeMethod\n'
        )
        for row in recommendation_rows:
            handle.write(
                f'{_clean_field(row.get("SequenceID", ""))}\t{_clean_field(row.get("Recommendation", ""))}\t'
                f'{_clean_field(row.get("QCClass", ""))}\t{_clean_field(row.get("Reasons", ""))}\t'
                f'{_clean_field(row.get("SuggestedAction", ""))}\t{_clean_field(row.get("Status", ""))}\t'
                f'{_clean_field(row.get("ProcessingMode", ""))}\t{row.get("OutputLength", 0)}\t'
                f'{row.get("MeanQuality", 0.0)}\t{row.get("OutputExpectedErrors", 0.0)}\t'
                f'{row.get("OutputNPercent", 0.0)}\t{row.get("ConflictDensity", 0.0)}\t'
                f'{row.get("ReadCount", 0)}\t{row.get("UsedReads", 0)}\t'
                f'{_clean_field(row.get("ReadIDs", ""))}\t{_clean_field(row.get("KeptReadIDs", ""))}\t'
                f'{_clean_field(row.get("FailedReadIDs", ""))}\t{_clean_field(row.get("SelectedReadID", ""))}\t'
                f'{_clean_field(row.get("MergeMethod", ""))}\n'
            )
    qc_policy_tsv.write_text(
        '\n'.join([
            'Metric\tDefault\tMeaning',
            f'min_quality\t{min_quality}\tPhred cutoff used for Mott-style end trimming',
            f'mask_quality\t{mask_quality}\tBases below this Phred score inside the retained read are masked to N',
            f'min_read_length\t{read_min_length}\tMinimum retained read length before assembly or best-read selection',
            f'min_final_length\t{final_min_length}\tMinimum final isolate sequence length written to assembled.fasta',
            f'min_mean_quality\t{float(min_mean_quality):g}\tMinimum mean Phred score after masking for read/final QC',
            f'max_read_expected_errors\t{float(max_read_expected_errors):g}\tMaximum expected base-call errors allowed per retained read',
            f'max_output_expected_errors\t{float(max_output_expected_errors):g}\tMaximum expected base-call errors allowed in the final isolate sequence',
            f'max_n_percent\t{float(max_n_percent):g}\tMaximum percent N in read/final sequence after internal masking',
            f'max_internal_low_quality_run\t{int(max_internal_low_quality_run)}\tLongest internal low-quality/ambiguous run allowed before read failure',
            f'max_conflict_density\t{float(max_conflict_density):g}\tMaximum overlap conflicts per 100 final bases before final failure',
            f'quality_difference\t{int(quality_difference)}\tMinimum Phred difference required to choose one conflicting base over another',
            f'allow_missing_quality\t{bool(allow_missing_quality)}\tWhether AB1 reads missing PCON quality scores may pass QC',
            '',
        ])
    )
    kept_reads = sum(1 for r in reads if r.status == 'kept')
    assembled_count = len(assembled_records)
    filtered_outputs = sum(
        1 for row in assembly_rows
        if str(row.get('QCClass', 'FAIL_QC')) == 'FAIL_QC'
    )
    resequence_count = sum(1 for row in assembly_rows if row.get('Recommendation') == 'RESEQUENCE')
    review_count = sum(1 for row in assembly_rows if row.get('Recommendation') == 'MANUAL_REVIEW')
    accept_count = sum(1 for row in assembly_rows if row.get('Recommendation') == 'ACCEPT')
    mode_counts: dict[str, int] = {}
    for row in assembly_rows:
        mode = str(row.get('ProcessingMode') or 'unknown')
        mode_counts[mode] = mode_counts.get(mode, 0) + 1
    summary_txt.write_text(
        '\n'.join([
            'BranchManager Sanger/AB1 Processing Summary',
            f'Input files: {len(files)}',
            f'Reads parsed: {len(reads)}',
            f'Minimum trimmed read length: {read_min_length}',
            f'Minimum final sequence length: {final_min_length}',
            f'Reads kept after trimming: {kept_reads}',
            f'Output assembled sequences: {assembled_count}',
            f'Final QC accept: {accept_count}',
            f'Final QC manual review: {review_count}',
            f'Final QC resequence: {resequence_count}',
            f'Output groups failing final QC: {filtered_outputs}',
            'Processing modes: ' + ', '.join(f'{mode}={count}' for mode, count in sorted(mode_counts.items())),
            f'Final FASTA: {assembled_fasta}',
            f'Read QC: {read_qc_tsv}',
            f'Per-base error profile: {per_base_error_tsv}',
            f'Assembly report: {assembly_tsv}',
            f'Resequencing recommendations: {recommendations_tsv}',
            f'QC policy: {qc_policy_tsv}',
            f'Read profile visual: {read_error_svg}',
            f'Assembly visual: {assembly_svg}',
            '',
        ])
    )
    return {
        'raw_fasta': str(raw_fasta),
        'trimmed_fasta': str(trimmed_fasta),
        'assembled_fasta': str(assembled_fasta),
        'read_qc_tsv': str(read_qc_tsv),
        'per_base_error_tsv': str(per_base_error_tsv),
        'assembly_tsv': str(assembly_tsv),
        'recommendations_tsv': str(recommendations_tsv),
        'qc_policy_tsv': str(qc_policy_tsv),
        'read_error_svg': str(read_error_svg),
        'assembly_svg': str(assembly_svg),
        'summary': str(summary_txt),
    }
