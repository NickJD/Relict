from __future__ import annotations

import csv
import gzip
import io
import logging
import math
import re
import struct
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

from branchmanager.utils.fasta import read_fasta, reverse_complement, write_fasta


logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class PaperTrailQCPolicy:
    """Single source of truth for scientifically conservative Sanger QC defaults."""

    min_quality: int = 20
    min_read_length: int = 300
    min_final_length: int = 800
    min_mean_quality: float = 20.0
    mask_quality: int = 20
    max_read_expected_errors: float = 8.0
    max_output_expected_errors: float = 5.0
    warn_n_percent: float = 3.0
    max_n_percent: float = 5.0
    warn_internal_low_quality_run: int = 5
    max_internal_low_quality_run: int = 20
    max_conflict_density: float = 1.0
    secondary_peak_ratio: float = 0.33
    max_mixed_peak_percent: float = 15.0
    mixed_peak_min_quality: int = 20
    quality_difference: int = 10
    min_overlap: int = 40
    min_overlap_identity: float = 0.95


PAPER_TRAIL_QC_POLICY_VERSION = '2.0'
DEFAULT_QC_POLICY = PaperTrailQCPolicy()


DEFAULT_PRIMERS = (
    # 16S rRNA
    '8F', '27F', '63F', '338F', '341F', '357F', '515F', '519F', '785F', '1055F',
    '534R', '806R', '907R', '926R', '1100R', '1392R', '1406R', '1492R', '1525R',
    # 18S rRNA
    'NS1', 'NS2', 'NS3', 'NS4', 'NS5', 'NS6', 'NS7', 'NS8', '528F', '1391F', 'EukBr',
    # ITS
    'ITS1', 'ITS1F', 'ITS2', 'ITS3', 'ITS4', 'ITS4B', 'ITS5',
    # 28S / LSU rRNA
    'LR0R', 'LR3', 'LR5', 'LR6', 'NL1', 'NL4',
)

DEFAULT_MAX_REPORT_IMAGE_HEIGHT = 2400
MIN_REPORT_IMAGE_HEIGHT = 600
READ_ERROR_HEADER_HEIGHT = 76
READ_ERROR_PANEL_HEIGHT = 132
CHROMATOGRAM_HEADER_HEIGHT = 82
CHROMATOGRAM_PANEL_HEIGHT = 178
ASSEMBLY_HEADER_HEIGHT = 104


DEFAULT_PRIMER_SEQUENCES = {
    # 16S rRNA
    '8F': 'AGAGTTTGATCCTGGCTCAG',          # Near full-length forward, non-degenerate variant of 27F
    '27F': 'AGAGTTTGATCMTGGCTCAG',          # Standard near full-length forward
    '63F': 'CAGGCCTAACACATGCAAGTC',          # Universal forward
    '338F': 'ACTCCTACGGGAGGCAGCAG',          # V3 forward (Muyzer et al. 1993)
    '341F': 'CCTACGGGNGGCWGCAG',             # V3-V4 forward
    '357F': 'CTCCTACGGGAGGCAGCAG',           # V3 forward variant
    '515F': 'GTGYCAGCMGCCGCGGTAA',           # V4 forward (Earth Microbiome Project)
    '519F': 'CAGCMGCCGCGGTAANWC',            # V4 forward alternative
    '785F': 'GGATTAGATACCCBDGTAGTC',         # V5 forward
    '1055F': 'ATGGCTGTCGTCAGCT',             # V6 forward
    '534R': 'ATTACCGCGGCTGCTGG',             # V1-V3 reverse
    '806R': 'GGACTACNVGGGTWTCTAAT',          # V4 reverse (Earth Microbiome Project)
    '907R': 'CCGTCAATTCMTTTRAGTTT',          # V5 reverse
    '926R': 'CCGYCAATTYMTTTRAGTTT',          # V4-V5 reverse
    '1100R': 'GGGTTGCGCTCGTTG',              # V6 reverse
    '1392R': 'ACGGGCGGTGTGTRC',              # Near full-length reverse
    '1406R': 'ACGGGCGGTGTGTRCAA',            # Extended near full-length reverse
    '1492R': 'TACGGYTACCTTGTTACGACTT',       # Standard full-length reverse
    '1525R': 'AAGGAGGTGWTCCARCC',            # Full-length reverse
    # ITS (Internal Transcribed Spacer)
    'ITS1': 'TCCGTAGGTGAACCTGCGG',           # Universal ITS1 forward (White et al. 1990)
    'ITS1F': 'CTTGGTCATTTAGAGGAAGTAA',       # Fungal-specific ITS1 forward (Gardes & Bruns 1993)
    'ITS2': 'GCTGCGTTCTTCATCGATGC',          # ITS1/5.8S boundary reverse (White et al. 1990)
    'ITS3': 'GCATCGATGAAGAACGCAGC',          # 5.8S/ITS2 boundary forward (White et al. 1990)
    'ITS4': 'TCCTCCGCTTATTGATATGC',          # Universal ITS2 reverse (White et al. 1990)
    'ITS4B': 'CAGGAGACTTGTACACGGTCCAG',      # Basidiomycete-specific ITS2 reverse (Gardes & Bruns 1993)
    'ITS5': 'GGAAGTAAAAGTCGTAACAAGG',        # Upstream of ITS1 forward (White et al. 1990)
    # 18S rRNA
    'NS1': 'GTAGTCATATGCTTGTCTC',            # Full 18S forward (White et al. 1990)
    'NS2': 'GGCTGCTGGCACCAGACTTGC',          # Full 18S reverse (White et al. 1990)
    'NS3': 'GCAAGTCTGGTGCCAGCAGCC',          # Full 18S forward (White et al. 1990)
    'NS4': 'CTTCCGTCAATTCCTTTAAG',           # Full 18S reverse (White et al. 1990)
    'NS5': 'AACTTAAAGGAATTGACGGAAG',         # Full 18S forward (White et al. 1990)
    'NS6': 'GCATCACAGACCTGTTATTGCCTC',       # Full 18S reverse (White et al. 1990)
    'NS7': 'GAGGCAATAACAGGTCTGTGATGC',       # Full 18S forward (White et al. 1990)
    'NS8': 'TCCGCAGGTTCACCTACGGA',           # Full 18S reverse (White et al. 1990)
    '528F': 'GCGGTAATTCCAGCTCCAA',           # V4 region forward (Stoeck et al. 2010)
    '1391F': 'GTACACACCGCCCGTC',             # V9 region forward (Earth Microbiome Project)
    'EukBr': 'TGATCCTTCTGCAGGTTCACCTAC',     # V9 region reverse (Earth Microbiome Project)
    # 28S / LSU rRNA
    'LR0R': 'ACCCGCTGAACTTAAGC',             # Universal; primes into LSU from 18S (Vilgalys & Hester 1990)
    'LR3': 'CCGTGTTTCAAGACGGG',              # D1 domain reverse (Vilgalys & Hester 1990)
    'LR5': 'TCCTGAGGGAAACTTCG',              # D2 domain reverse (Vilgalys & Hester 1990)
    'LR6': 'CGCCAGTTCTGCTTACC',              # D3 domain reverse (Vilgalys & Hester 1990)
    'NL1': 'GCATATCAATAAGCGGAGGAAAAG',       # D1/D2 domain forward (O'Donnell 1993)
    'NL4': 'GGTCCGTGTTTCAAGACGG',            # D1/D2 domain reverse (O'Donnell 1993)
}

IUPAC_BASES = {
    'A': {'A'}, 'C': {'C'}, 'G': {'G'}, 'T': {'T'},
    'R': {'A', 'G'}, 'Y': {'C', 'T'}, 'S': {'G', 'C'}, 'W': {'A', 'T'},
    'K': {'G', 'T'}, 'M': {'A', 'C'}, 'B': {'C', 'G', 'T'},
    'D': {'A', 'G', 'T'}, 'H': {'A', 'C', 'T'}, 'V': {'A', 'C', 'G'},
    'N': {'A', 'C', 'G', 'T'},
}

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
    'read_mode', 'read_handling', 'read_strategy',
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
    trace_available: bool = False
    trace_order: str = ''
    peak_locations: list[int] | None = None
    trace_channels: dict[str, list[int]] | None = None
    secondary_peak_ratios: list[float] | None = None
    secondary_peak_bases: list[str] | None = None
    mixed_peak_count: int = 0
    mixed_peak_percent: float = 0.0
    primer_trimmed_bases: int = 0
    taxonomy_screen: str = ''


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


def _decode_abif_text(value: object) -> str:
    if isinstance(value, bytes):
        return value.decode('ascii', errors='ignore').replace('\x00', '').strip()
    return str(value or '').replace('\x00', '').strip()


def _integer_list(value: object) -> list[int]:
    if value is None:
        return []
    if isinstance(value, bytes):
        return [int(item) for item in value]
    if isinstance(value, (list, tuple)):
        return [int(item) for item in value]
    try:
        return [int(value)]
    except (TypeError, ValueError):
        return []


def _trace_peak_evidence(
    sequence: str,
    peak_locations: list[int],
    channels: dict[str, list[int]],
) -> tuple[list[float], list[str]]:
    ratios = []
    secondary_bases = []
    for idx, called in enumerate(sequence):
        if idx >= len(peak_locations):
            ratios.append(0.0)
            secondary_bases.append('')
            continue
        peak = peak_locations[idx]
        signals = {
            base: values[peak] if 0 <= peak < len(values) else 0
            for base, values in channels.items()
            if base in {'A', 'C', 'G', 'T'}
        }
        if not signals:
            ratios.append(0.0)
            secondary_bases.append('')
            continue
        ordered = sorted(signals.items(), key=lambda item: (item[1], item[0]), reverse=True)
        primary_signal = max(1, signals.get(called, ordered[0][1]))
        secondary = max(
            ((base, signal) for base, signal in signals.items() if base != called),
            key=lambda item: (item[1], item[0]),
            default=('', 0),
        )
        ratios.append(float(secondary[1]) / float(primary_signal))
        secondary_bases.append(secondary[0])
    return ratios, secondary_bases


def _read_ab1_details(path: str | Path) -> dict[str, object]:
    """Read calls, Phred values, peak positions, and four dye channels from ABIF."""
    try:
        from Bio import SeqIO
    except ImportError as exc:
        raise RuntimeError(
            'Biopython is required for scientific AB1 parsing. Install it in the active environment.'
        ) from exc

    try:
        if str(path).lower().endswith('.gz'):
            with gzip.open(path, 'rb') as handle:
                binary = io.BytesIO(handle.read())
                record = SeqIO.read(binary, 'abi')
        else:
            with open(path, 'rb') as handle:
                record = SeqIO.read(handle, 'abi')
    except Exception as biopython_error:
        # Some legacy/minimal ABIF writers use non-standard PCON encodings that
        # Biopython rejects. Preserve their calls and qualities for review, but
        # do not claim chromatogram-level evidence without decoded trace arrays.
        entries = _read_abif_entries(path)
        called = entries.get(('PBAS', 2)) or entries.get(('PBAS', 1))
        if not called:
            raise ValueError(f'AB1 contains no readable PBAS base calls: {path}') from biopython_error
        sequence = _decode_abif_text(called).upper()
        sequence = ''.join(base if base in IUPAC_BASES else 'N' for base in sequence)
        quality_bytes = entries.get(('PCON', 2)) or entries.get(('PCON', 1)) or b''
        qualities = _integer_list(quality_bytes)[:len(sequence)]
        quality_available = bool(qualities)
        if len(qualities) < len(sequence):
            qualities.extend([0] * (len(sequence) - len(qualities)))
        return {
            'sequence': sequence,
            'qualities': qualities,
            'quality_available': quality_available,
            'trace_available': False,
            'trace_order': '',
            'peak_locations': [],
            'trace_channels': {},
            'secondary_peak_ratios': [0.0] * len(sequence),
            'secondary_peak_bases': [''] * len(sequence),
            'parser_warning': 'legacy_abif_calls_only',
        }
    sequence = str(record.seq).strip().upper()
    sequence = ''.join(base if base in IUPAC_BASES else 'N' for base in sequence)
    qualities = [int(value) for value in record.letter_annotations.get('phred_quality', [])]
    quality_available = bool(qualities)
    if len(qualities) < len(sequence):
        qualities.extend([0] * (len(sequence) - len(qualities)))
    qualities = qualities[:len(sequence)]

    raw = record.annotations.get('abif_raw', {}) or {}
    peak_locations = _integer_list(raw.get('PLOC2') or raw.get('PLOC1'))
    order = _decode_abif_text(raw.get('FWO_1') or raw.get('FWO_2') or 'GATC').upper()
    if len(order) != 4 or set(order) != {'A', 'C', 'G', 'T'}:
        order = 'GATC'
    channels = {}
    for base, tag in zip(order, ('DATA9', 'DATA10', 'DATA11', 'DATA12')):
        values = _integer_list(raw.get(tag))
        if values:
            channels[base] = values
    trace_available = bool(peak_locations and len(channels) == 4)
    if trace_available:
        ratios, secondary = _trace_peak_evidence(sequence, peak_locations, channels)
    else:
        ratios, secondary = [0.0] * len(sequence), [''] * len(sequence)
    return {
        'sequence': sequence,
        'qualities': qualities,
        'quality_available': quality_available,
        'trace_available': trace_available,
        'trace_order': order,
        'peak_locations': peak_locations,
        'trace_channels': channels,
        'secondary_peak_ratios': ratios,
        'secondary_peak_bases': secondary,
    }


def _read_ab1_with_quality_flag(path: str | Path) -> tuple[str, list[int], bool]:
    details = _read_ab1_details(path)
    return (
        str(details['sequence']),
        list(details['qualities']),
        bool(details['quality_available']),
    )


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


def _iupac_matches(left: str, right: str) -> bool:
    return bool(IUPAC_BASES.get(left.upper(), {left.upper()}) & IUPAC_BASES.get(right.upper(), {right.upper()}))


def trim_leading_primer(
    sequence: str,
    qualities: list[int],
    primer_name: str,
    *,
    primer_sequences: Optional[dict[str, str]] = None,
    min_identity: float = 0.80,
    max_offset: int = 5,
) -> tuple[str, list[int], int]:
    """Remove a confidently observed sequencing-primer sequence at the read start."""
    lookup = primer_sequences or DEFAULT_PRIMER_SEQUENCES
    primer = str(lookup.get(str(primer_name).upper()) or lookup.get(str(primer_name)) or '').upper()
    seq = str(sequence or '').upper()
    if not primer or len(primer) < 12 or len(seq) < 12:
        return seq, list(qualities), 0
    best = None
    for offset in range(0, min(max_offset, max(0, len(seq) - 12)) + 1):
        compared = min(len(primer), len(seq) - offset)
        if compared < 12:
            continue
        matches = sum(_iupac_matches(seq[offset + idx], primer[idx]) for idx in range(compared))
        identity = matches / compared
        candidate = (identity, compared, -offset)
        if identity >= min_identity and (best is None or candidate > best[0]):
            best = (candidate, offset + compared)
    if best is None:
        return seq, list(qualities), 0
    trim_end = best[1]
    return seq[trim_end:], list(qualities)[trim_end:], trim_end


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


def _load_sample_map(
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
            logger.warning('[PAPER TRAIL] Ignoring unsupported or missing input path: %s', raw)
    seen = set()
    unique = []
    for path in sorted(files):
        key = str(path)
        if key in seen:
            continue
        seen.add(key)
        unique.append(path)
    return unique


def _read_records_from_file(path: Path) -> list[tuple[str, str, list[int], bool, str, dict]]:
    lower = path.name.lower()
    if lower.endswith(('.ab1', '.abi', '.ab1.gz', '.abi.gz')):
        details = _read_ab1_details(path)
        seq = str(details['sequence'])
        qual = list(details['qualities'])
        quality_available = bool(details['quality_available'])
        read_id = Path(path.stem).stem if lower.endswith('.gz') else path.stem
        source = 'ab1_pcon' if quality_available else 'ab1_missing_pcon'
        return [(read_id, seq, qual, quality_available, source, details)]
    if lower.endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz')):
        return [(rid, seq, qual, True, 'fastq', {}) for rid, seq, qual in read_fastq(path)]
    if lower.endswith(('.fasta', '.fa', '.fna', '.fasta.gz', '.fa.gz', '.fna.gz')):
        return [
            (h, s.upper(), [0] * len(s), False, 'fasta_quality_unverified', {})
            for h, s in read_fasta(str(path))
        ]
    raise ValueError(f'Unsupported input file: {path}')


def _phred_from_error(error_probability: float, *, maximum: int = 60) -> int:
    error = min(1.0, max(10 ** (-maximum / 10), float(error_probability)))
    return min(maximum, max(0, int(round(-10.0 * math.log10(error)))))


def _posterior_consensus(
    base_a: str,
    q_a: int,
    base_b: str,
    q_b: int,
    *,
    quality_difference: int = DEFAULT_QC_POLICY.quality_difference,
) -> tuple[str, int, bool]:
    """Infer one consensus base under an independent Phred observation model."""
    a = str(base_a).upper()
    b = str(base_b).upper()
    if a == 'N':
        return b, int(q_b), False
    if b == 'N':
        return a, int(q_a), False
    p_a = min(0.75, _error_probability(q_a))
    p_b = min(0.75, _error_probability(q_b))
    if a == b:
        correct = (1.0 - p_a) * (1.0 - p_b)
        wrong = (p_a * p_b) / 3.0
        posterior_error = wrong / max(correct + wrong, 1e-12)
        return a, _phred_from_error(posterior_error), False

    likelihoods = {
        a: (1.0 - p_a) * (p_b / 3.0),
        b: (p_a / 3.0) * (1.0 - p_b),
    }
    for alternative in {'A', 'C', 'G', 'T'} - {a, b}:
        likelihoods[alternative] = (p_a / 3.0) * (p_b / 3.0)
    total = max(sum(likelihoods.values()), 1e-12)
    winner, winner_likelihood = max(likelihoods.items(), key=lambda item: (item[1], item[0]))
    posterior = winner_likelihood / total
    quality = _phred_from_error(1.0 - posterior)
    odds = 10 ** (max(0, int(quality_difference)) / 10.0)
    minimum_posterior = odds / (1.0 + odds)
    if posterior < minimum_posterior:
        return 'N', min(int(q_a), int(q_b), quality), True
    return winner, quality, False


def _pairwise_aligner():
    try:
        from Bio import Align
    except ImportError as exc:
        raise RuntimeError('Biopython is required for quality-aware Sanger assembly.') from exc
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2.0
    aligner.mismatch_score = -3.0
    aligner.open_gap_score = -6.0
    aligner.extend_gap_score = -1.0
    aligner.end_gap_score = 0.0
    aligner.wildcard = 'N'
    return aligner


def _pairwise_alignment_columns(seq_a: str, qual_a: list[int], seq_b: str, qual_b: list[int]):
    aligner = _pairwise_aligner()
    alignments = aligner.align(seq_a, seq_b)
    if not alignments:
        return []
    coordinates = alignments[0].coordinates
    columns = []
    for segment in range(coordinates.shape[1] - 1):
        a0, a1 = int(coordinates[0, segment]), int(coordinates[0, segment + 1])
        b0, b1 = int(coordinates[1, segment]), int(coordinates[1, segment + 1])
        advance_a = a1 - a0
        advance_b = b1 - b0
        if advance_a and advance_b:
            if advance_a != advance_b:
                raise RuntimeError('Unexpected unequal aligned block from Biopython PairwiseAligner')
            for offset in range(advance_a):
                ia, ib = a0 + offset, b0 + offset
                columns.append((seq_a[ia], qual_a[ia] if ia < len(qual_a) else 0,
                                seq_b[ib], qual_b[ib] if ib < len(qual_b) else 0))
        elif advance_a:
            for ia in range(a0, a1):
                columns.append((seq_a[ia], qual_a[ia] if ia < len(qual_a) else 0, '-', None))
        elif advance_b:
            for ib in range(b0, b1):
                columns.append(('-', None, seq_b[ib], qual_b[ib] if ib < len(qual_b) else 0))
    return columns


def _align_read_to_consensus(consensus: str, sequence: str) -> list[tuple[int, int, int, int]]:
    """Return aligned consensus/read blocks as zero-based half-open coordinates."""
    if not consensus or not sequence:
        return []
    alignments = _pairwise_aligner().align(consensus, sequence)
    if not alignments:
        return []
    consensus_blocks, read_blocks = alignments[0].aligned
    return [
        (int(consensus_start), int(consensus_end), int(read_start), int(read_end))
        for (consensus_start, consensus_end), (read_start, read_end)
        in zip(consensus_blocks, read_blocks)
        if int(consensus_end) > int(consensus_start) and int(read_end) > int(read_start)
    ]


def _read_sequence_for_placement(read: SangerRead) -> str:
    if read.oriented_sequence:
        return read.oriented_sequence
    sequence = read.masked_sequence or read.trimmed_sequence or read.raw_sequence
    return reverse_complement(sequence) if read.direction == 'reverse' else sequence


def _build_read_placements(
    consensus: str,
    reads: list[SangerRead],
    used_read_ids: Iterable[str],
) -> dict[str, dict[str, object]]:
    used = {str(read_id) for read_id in used_read_ids}
    placements = {}
    for read in reads:
        sequence = _read_sequence_for_placement(read)
        segments = _align_read_to_consensus(consensus, sequence)
        placements[read.read_id] = {
            'segments': segments,
            'consensus_start': min((segment[0] for segment in segments), default=None),
            'consensus_end': max((segment[1] for segment in segments), default=None),
            'aligned_bases': sum(min(segment[1] - segment[0], segment[3] - segment[2]) for segment in segments),
            'read_length': len(sequence),
            'contributes': read.read_id in used,
        }
    return placements


def find_best_overlap(
    seq_a: str,
    seq_b: str,
    *,
    min_overlap: int = DEFAULT_QC_POLICY.min_overlap,
    min_identity: float = DEFAULT_QC_POLICY.min_overlap_identity,
):
    columns = _pairwise_alignment_columns(seq_a, [30] * len(seq_a), seq_b, [30] * len(seq_b))
    paired = [idx for idx, (a, _qa, b, _qb) in enumerate(columns) if a != '-' and b != '-']
    if not paired:
        return None
    first, last = paired[0], paired[-1]
    matches = 0
    compared = 0
    indel_columns = 0
    for idx in range(first, last + 1):
        a, _qa, b, _qb = columns[idx]
        if a == '-' or b == '-':
            indel_columns += 1
        elif a != 'N' and b != 'N':
            compared += 1
            matches += int(a == b)
    denominator = compared + indel_columns
    identity = matches / denominator if denominator else 0.0
    if compared < min_overlap or identity < min_identity:
        return None
    return {
        'offset': 'pairwise',
        'overlap': last - first + 1,
        'compared': compared,
        'identity': identity,
        'indel_columns': indel_columns,
        'columns': columns,
        'overlap_start': first,
        'overlap_end': last,
    }


def merge_pair(
    seq_a: str,
    qual_a: list[int],
    seq_b: str,
    qual_b: list[int],
    *,
    min_overlap: int = DEFAULT_QC_POLICY.min_overlap,
    min_identity: float = DEFAULT_QC_POLICY.min_overlap_identity,
    quality_difference: int = DEFAULT_QC_POLICY.quality_difference,
):
    overlap = find_best_overlap(seq_a, seq_b, min_overlap=min_overlap, min_identity=min_identity)
    if overlap is None:
        return None
    consensus = []
    consensus_qual = []
    conflicts = 0
    ambiguous_conflicts = 0
    first = int(overlap['overlap_start'])
    last = int(overlap['overlap_end'])
    for idx, (base_a, q_a, base_b, q_b) in enumerate(overlap['columns']):
        if base_a != '-' and base_b != '-':
            base, quality, ambiguous = _posterior_consensus(
                base_a,
                int(q_a or 0),
                base_b,
                int(q_b or 0),
                quality_difference=quality_difference,
            )
            consensus.append(base)
            consensus_qual.append(quality)
            if base_a != base_b and base_a != 'N' and base_b != 'N':
                conflicts += 1
                ambiguous_conflicts += int(ambiguous)
        elif base_a != '-':
            consensus.append(base_a)
            if first <= idx <= last:
                consensus_qual.append(min(10, int(q_a or 0)))
                conflicts += 1
                ambiguous_conflicts += 1
            else:
                consensus_qual.append(int(q_a or 0))
        elif base_b != '-':
            consensus.append(base_b)
            if first <= idx <= last:
                consensus_qual.append(min(10, int(q_b or 0)))
                conflicts += 1
                ambiguous_conflicts += 1
            else:
                consensus_qual.append(int(q_b or 0))
    return ''.join(consensus), consensus_qual, {
        'offset': overlap['offset'],
        'overlap': overlap['overlap'],
        'compared': overlap['compared'],
        'identity': overlap['identity'],
        'indel_columns': overlap['indel_columns'],
        'conflicts': conflicts,
        'ambiguous_conflicts': ambiguous_conflicts,
        'method': 'biopython_pairwise_overlap',
    }


def _append_consensus_base(
    consensus: list[str],
    consensus_qual: list[int],
    base_a: str,
    q_a: int,
    base_b: str,
    q_b: int,
    *,
    quality_difference: int = DEFAULT_QC_POLICY.quality_difference,
) -> tuple[int, int]:
    conflicts = 0
    ambiguous_conflicts = 0
    base, quality, ambiguous = _posterior_consensus(
        base_a,
        q_a,
        base_b,
        q_b,
        quality_difference=quality_difference,
    )
    if base_a != base_b and base_a != 'N' and base_b != 'N':
        conflicts = 1
        ambiguous_conflicts = int(ambiguous)
    consensus.append(base)
    consensus_qual.append(quality)
    return conflicts, ambiguous_conflicts


def merge_pair_gapped(*args, **kwargs):
    """Compatibility wrapper; the primary merger now models gaps directly."""
    return merge_pair(*args, **kwargs)


def assemble_read_group(
    reads: list[SangerRead],
    *,
    min_overlap: int = DEFAULT_QC_POLICY.min_overlap,
    min_identity: float = DEFAULT_QC_POLICY.min_overlap_identity,
    quality_difference: int = DEFAULT_QC_POLICY.quality_difference,
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
            'unmerged_reads': ';'.join(read.read_id for read in reads),
            'method': 'failed_no_reads',
            'qualities': [],
            'used_read_ids': [],
        }
    usable.sort(key=lambda r: (r.direction != 'forward', -len(r.oriented_sequence), r.read_id))
    if len(usable) == 1:
        only = usable[0]
        consensus = only.oriented_sequence
        qualities = only.oriented_qualities or [0] * len(consensus)
        return consensus, {
            'status': 'single_read', 'read_count': len(reads), 'used_reads': 1,
            'overlap_identity': 'NA', 'overlap_length': 'NA', 'conflicts': 0,
            'ambiguous_conflicts': 0, 'unmerged_reads': '',
            'mean_quality': round(_mean_quality(qualities), 2),
            'method': 'single_read', 'qualities': qualities,
            'used_read_ids': [only.read_id],
        }

    # Choose the strongest pair first. This avoids making the consensus depend
    # on filename or primer order when three or more reads are supplied.
    pair_candidates = []
    for left_idx in range(len(usable)):
        for right_idx in range(left_idx + 1, len(usable)):
            left = usable[left_idx]
            right = usable[right_idx]
            got = merge_pair(
                left.oriented_sequence,
                left.oriented_qualities or [0] * len(left.oriented_sequence),
                right.oriented_sequence,
                right.oriented_qualities or [0] * len(right.oriented_sequence),
                min_overlap=min_overlap,
                min_identity=min_identity,
                quality_difference=quality_difference,
            )
            if got is None:
                continue
            stats = got[2]
            key = (
                float(stats.get('identity', 0.0)), int(stats.get('compared', 0)),
                int(stats.get('overlap', 0)), len(got[0]),
                left.read_id, right.read_id,
            )
            pair_candidates.append((key, left_idx, right_idx, got))
    if not pair_candidates:
        seed = usable[0]
        return seed.oriented_sequence, {
            'status': 'partial_no_overlap', 'read_count': len(reads), 'used_reads': 1,
            'overlap_identity': 'NA', 'overlap_length': 'NA', 'conflicts': 0,
            'ambiguous_conflicts': 0,
            'unmerged_reads': ';'.join(read.read_id for read in usable[1:]),
            'mean_quality': round(_mean_quality(seed.oriented_qualities or []), 2),
            'method': 'unmerged', 'qualities': seed.oriented_qualities or [],
            'used_read_ids': [seed.read_id],
        }

    _key, left_idx, right_idx, initial = max(pair_candidates, key=lambda item: item[0])
    consensus, qualities, first_stats = initial
    merge_stats = [first_stats]
    used_indexes = {left_idx, right_idx}
    remaining = {idx for idx in range(len(usable)) if idx not in used_indexes}
    while remaining:
        next_candidates = []
        for idx in sorted(remaining):
            read = usable[idx]
            got = merge_pair(
                consensus,
                qualities,
                read.oriented_sequence,
                read.oriented_qualities or [0] * len(read.oriented_sequence),
                min_overlap=min_overlap,
                min_identity=min_identity,
                quality_difference=quality_difference,
            )
            if got is None:
                continue
            stats = got[2]
            key = (
                float(stats.get('identity', 0.0)), int(stats.get('compared', 0)),
                int(stats.get('overlap', 0)), len(got[0]), read.read_id,
            )
            next_candidates.append((key, idx, got))
        if not next_candidates:
            break
        _key, idx, got = max(next_candidates, key=lambda item: item[0])
        consensus, qualities, stats = got
        merge_stats.append(stats)
        used_indexes.add(idx)
        remaining.remove(idx)
    unmerged = [usable[idx].read_id for idx in sorted(remaining)]
    used = len(used_indexes)
    if merge_stats:
        status = 'assembled'
        overlap_identity = min(stat['identity'] for stat in merge_stats)
        overlap_length = min(stat['overlap'] for stat in merge_stats)
        conflicts = sum(stat['conflicts'] for stat in merge_stats)
        ambiguous = sum(stat['ambiguous_conflicts'] for stat in merge_stats)
        method = ';'.join(sorted({str(stat.get('method', 'unknown')) for stat in merge_stats}))
    else:
        raise RuntimeError('Internal assembly error: no merge statistics after selecting a seed pair')
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
        'used_read_ids': [usable[idx].read_id for idx in sorted(used_indexes)],
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
            'unmerged_reads': ';'.join(read.read_id for read in reads),
            'mean_quality': 0.0,
            'selected_read_id': '',
            'method': 'best_read',
            'qualities': [],
            'used_read_ids': [],
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
        'used_read_ids': [best.read_id],
    }


def _classify_read_qc(
    read: SangerRead,
    *,
    read_min_length: int,
    min_mean_quality: float,
    max_read_expected_errors: float,
    warn_n_percent: float,
    max_n_percent: float,
    warn_internal_low_quality_run: int,
    max_internal_low_quality_run: int,
    allow_missing_quality: bool,
    max_mixed_peak_percent: float,
) -> None:
    fail_reasons: list[str] = []
    warn_reasons: list[str] = []
    sequence = read.masked_sequence or read.trimmed_sequence
    qualities = read.masked_qualities or read.trimmed_qualities or []
    expected_errors = _expected_errors(qualities, sequence)
    mean_q = _mean_quality(qualities)
    n_pct = _n_percent(sequence)

    quality_verified = bool(read.quality_available)
    if read.quality_source == 'fasta_quality_unverified':
        warn_reasons.append('quality_scores_unavailable_fasta_not_assumed')
    elif not quality_verified and not allow_missing_quality:
        fail_reasons.append(f'quality_scores_missing_{read.quality_source}')
    elif not quality_verified:
        warn_reasons.append(f'quality_scores_missing_allowed_{read.quality_source}')
    if read.quality_source.startswith('ab1_') and not read.trace_available:
        warn_reasons.append('chromatogram_trace_channels_unavailable')
    if read.mixed_peak_percent > float(max_mixed_peak_percent):
        fail_reasons.append(f'mixed_peak_percent_gt_{float(max_mixed_peak_percent):g}')
    elif read.mixed_peak_percent > max(1.0, float(max_mixed_peak_percent) / 2.0):
        warn_reasons.append(f'mixed_peak_percent_{read.mixed_peak_percent:.2f}')
    if len(sequence) < read_min_length:
        fail_reasons.append(f'trimmed_length_lt_{read_min_length}')
    if quality_verified and sequence and mean_q < min_mean_quality:
        fail_reasons.append(f'mean_q_lt_{min_mean_quality:g}')
    if quality_verified and expected_errors > max_read_expected_errors:
        fail_reasons.append(f'expected_errors_gt_{max_read_expected_errors:g}')
    if n_pct > max_n_percent:
        fail_reasons.append(f'n_percent_gt_{max_n_percent:g}')
    elif n_pct > warn_n_percent:
        warn_reasons.append(f'n_percent_gt_{warn_n_percent:g}')
    if quality_verified and read.longest_low_quality_run > max_internal_low_quality_run:
        fail_reasons.append(f'internal_low_quality_run_gt_{max_internal_low_quality_run}')
    elif read.longest_low_quality_run > warn_internal_low_quality_run:
        warn_reasons.append(
            f'internal_low_quality_run_gt_{warn_internal_low_quality_run}'
        )

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
    warn_n_percent: float,
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
    quality_verified = any(read.quality_available for read in group if read.status == 'kept')
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
    if quality_verified and consensus and mean_q < min_mean_quality:
        fail_reasons.append(f'output_mean_q_lt_{min_mean_quality:g}')
    if quality_verified and expected_errors > max_output_expected_errors:
        fail_reasons.append(f'output_expected_errors_gt_{max_output_expected_errors:g}')
    if n_pct > max_n_percent:
        fail_reasons.append(f'output_n_percent_gt_{max_n_percent:g}')
    elif n_pct > warn_n_percent:
        warn_reasons.append(f'output_n_percent_gt_{warn_n_percent:g}')
    if conflict_density > max_conflict_density:
        fail_reasons.append(f'overlap_conflict_density_gt_{max_conflict_density:g}')

    status = str(stats.get('status', 'unknown'))
    if status in {'failed_no_reads', 'partial_no_overlap'}:
        fail_reasons.append(status)
    elif status in {'single_read'} and len(group) > 1 and processing_mode == 'assemble':
        warn_reasons.append('only_one_primer_read_passed_qc')
    if ambiguous:
        warn_reasons.append(f'ambiguous_overlap_conflicts_{ambiguous}')
    unmerged = str(stats.get('unmerged_reads') or '').strip()
    if unmerged and processing_mode == 'assemble':
        warn_reasons.append('unmerged_primer_reads_' + str(len([value for value in unmerged.split(';') if value])))
    if stats.get('taxonomy_conflict'):
        fail_reasons.append('primer_read_taxonomic_conflict')
    if mode_warning:
        warn_reasons.append(mode_warning)
    failed_reads = [read.read_id for read in group if read.qc_class == 'FAIL_QC']
    if failed_reads and (processing_mode == 'assemble' or not consensus):
        warn_reasons.append('failed_reads_' + str(len(failed_reads)))
    used_read_ids = {str(value) for value in stats.get('used_read_ids', [])}
    used_reads = [read for read in group if read.read_id in used_read_ids]
    for reason in sorted({
        reason
        for read in used_reads
        if read.qc_class == 'PASS_WITH_WARNINGS'
        for reason in (read.qc_reasons or [])
        if not str(reason).startswith('n_percent_gt_')
    }):
        warn_reasons.append(f'read_{reason}')
    if 'gapped_anchor' in str(stats.get('method', '')):
        warn_reasons.append('gapped_overlap_fallback_used')

    qc_class = _qc_class_from_reasons(fail_reasons, warn_reasons)
    recommendation = _recommendation_from_qc_class(qc_class)
    if recommendation == 'RESEQUENCE':
        suggested_action = 'Repeat Sanger sequencing for this isolate; inspect chromatograms and primer setup before reuse.'
    elif recommendation == 'MANUAL_REVIEW':
        suggested_action = 'Inspect chromatograms/overlap manually; accept only if trace evidence supports the consensus.'
    else:
        suggested_action = 'Accept for downstream BranchManager Performance Review.'
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


def _clean_fasta_token(value: object) -> str:
    text = _clean_field(value).strip()
    return re.sub(r'\s+', '_', text) if text else 'unknown'


def _review_sequence_for_failed_read(read: SangerRead) -> str:
    seq = read.oriented_sequence or read.masked_sequence or read.trimmed_sequence or read.raw_sequence
    if not seq:
        return ''
    if not read.oriented_sequence and read.direction == 'reverse':
        return reverse_complement(seq).upper()
    return str(seq).upper()


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
            'KeptAfterTrim\tOrientedPosition\tPeakLocation\tSignalA\tSignalC\tSignalG\t'
            'SignalT\tSecondaryPeakBase\tSecondaryPeakRatio\tMixedPeakFlag\n'
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
                peak_locations = read.peak_locations or []
                peak = peak_locations[idx] if idx < len(peak_locations) else None
                channels = read.trace_channels or {}
                signals = {
                    called: (
                        channels.get(called, [])[peak]
                        if peak is not None and 0 <= peak < len(channels.get(called, [])) else 0
                    )
                    for called in ('A', 'C', 'G', 'T')
                }
                ratios = read.secondary_peak_ratios or []
                secondary = read.secondary_peak_bases or []
                ratio = ratios[idx] if idx < len(ratios) else 0.0
                mixed = ratio >= 0.33 and q >= 20
                handle.write(
                    f'{_clean_field(read.read_id)}\t{_clean_field(read.sequence_id)}\t'
                    f'{_clean_field(read.primer)}\t{_clean_field(read.direction)}\t'
                    f'{_clean_field(read.source_file)}\t{idx + 1}\t{base}\t{q}\t'
                    f'{_error_probability(q):.6g}\t{trim_region}\t'
                    f'{masked_low_quality}\t{base_after_mask}\t'
                    f'{"yes" if kept else "no"}\t{oriented_position}\t'
                    f'{peak + 1 if peak is not None else ""}\t{signals["A"]}\t{signals["C"]}\t'
                    f'{signals["G"]}\t{signals["T"]}\t'
                    f'{secondary[idx] if idx < len(secondary) else ""}\t{ratio:.4f}\t'
                    f'{"yes" if mixed else "no"}\n'
                )


def _load_report_font(size: int, *, bold: bool = False):
    try:
        from PIL import ImageFont
    except ImportError as exc:
        raise RuntimeError(
            'PNG report figures require Pillow. Install it with `conda install pillow`.'
        ) from exc
    candidates = [
        '/System/Library/Fonts/Supplemental/Arial Bold.ttf' if bold else '/System/Library/Fonts/Supplemental/Arial.ttf',
        '/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf' if bold else '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',
        'DejaVuSans-Bold.ttf' if bold else 'DejaVuSans.ttf',
    ]
    for candidate in candidates:
        try:
            return ImageFont.truetype(candidate, size=size)
        except (OSError, ValueError):
            continue
    return ImageFont.load_default()


def _fit_report_text(draw, value: object, font, max_width: int) -> str:
    text = _clean_field(value)
    if draw.textlength(text, font=font) <= max_width:
        return text
    suffix = '...'
    while text and draw.textlength(text + suffix, font=font) > max_width:
        text = text[:-1]
    return text + suffix


def _draw_dashed_horizontal(draw, x0: int, x1: int, y: int, *, fill: str, width: int = 1) -> None:
    dash = 6
    gap = 4
    x = x0
    while x < x1:
        draw.line((x, y, min(x + dash, x1), y), fill=fill, width=width)
        x += dash + gap


def _quality_points(
    qualities: list[int], x0: float, y0: float, width: float, height: float,
) -> list[tuple[float, float]]:
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
        points.append((x, y))
    return points


def _write_read_error_profile_png(
    path: Path,
    reads: list[SangerRead],
    *,
    page_note: str = '',
) -> dict[str, int]:
    try:
        from PIL import Image, ImageDraw
    except ImportError as exc:
        raise RuntimeError(
            'PNG report figures require Pillow. Install it with `conda install pillow`.'
        ) from exc
    width = 1280
    panel_height = READ_ERROR_PANEL_HEIGHT
    left = 340
    right = 34
    graph_width = width - left - right
    height = READ_ERROR_HEADER_HEIGHT + max(1, len(reads)) * panel_height
    image = Image.new('RGB', (width, height), '#ffffff')
    draw = ImageDraw.Draw(image)
    title_font = _load_report_font(21, bold=True)
    body_font = _load_report_font(12)
    label_font = _load_report_font(12)
    label_bold_font = _load_report_font(14, bold=True)
    small_font = _load_report_font(10)
    draw.text((24, 18), 'Sanger read error and trimming profiles', fill='#1f2937', font=title_font)
    draw.text(
        (24, 48),
        'Blue: Phred quality. Green: retained trim window. Orange: internally masked low-quality bases.',
        fill='#4b5563', font=body_font,
    )
    if page_note:
        draw.text(
            (width - 24 - draw.textlength(page_note, font=body_font), 22),
            page_note, fill='#4b5563', font=body_font,
        )
    if not reads:
        draw.text((24, 88), 'No reads parsed.', fill='#6b7280', font=label_font)
    for idx, read in enumerate(reads):
        top = 78 + idx * panel_height
        graph_top = top + 36
        graph_height = 64
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
        draw.rounded_rectangle(
            (18, top - 10, width - 18, top + panel_height - 16),
            radius=6, fill='#f9fafb', outline='#e5e7eb', width=1,
        )
        text_width = left - 60
        draw.text((30, top + 4), _fit_report_text(draw, read.read_id, label_bold_font, width - 64), fill='#111827', font=label_bold_font)
        draw.text(
            (30, top + 27),
            _fit_report_text(draw, f'isolate {read.sequence_id} | {read.primer} | {read.direction} | mode {read.processing_mode}', label_font, text_width),
            fill='#4b5563', font=label_font,
        )
        draw.text(
            (30, top + 50),
            _fit_report_text(draw, f'{read.status} | {read.qc_class} | {read.warning}', label_font, text_width),
            fill=status_colour, font=label_font,
        )
        draw.text(
            (30, top + 73),
            _fit_report_text(
                draw,
                f'raw {len(read.raw_sequence)} bp | trimmed {len(read.trimmed_sequence)} bp | masked {read.masked_bases} bp',
                label_font, text_width,
            ),
            fill='#4b5563', font=label_font,
        )
        draw.text(
            (30, top + 96),
            _fit_report_text(draw, f'mean Q {raw_mean:.1f} -> {trim_mean:.1f} -> {masked_mean:.1f}', label_font, text_width),
            fill='#4b5563', font=label_font,
        )
        draw.rectangle(
            (round(trim_start_x), graph_top, round(trim_end_x), graph_top + graph_height),
            fill='#dcfce7',
        )
        draw.line((left, graph_top, left + graph_width, graph_top), fill='#e5e7eb', width=1)
        draw.line((left, graph_top + graph_height, left + graph_width, graph_top + graph_height), fill='#d1d5db', width=1)
        _draw_dashed_horizontal(draw, left, left + graph_width, round(q20_y), fill='#f59e0b')
        _draw_dashed_horizontal(draw, left, left + graph_width, round(q30_y), fill='#10b981')
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
                draw.rectangle((round(x0), graph_top, max(round(x0) + 1, round(x1)), graph_top + graph_height), fill='#fed7aa')
                mask_start = None
        if points:
            draw.line(points, fill='#2563eb', width=2)
        draw.text((left, graph_top + graph_height + 6), '1', fill='#6b7280', font=small_font)
        end_label = f'{len(read.raw_sequence)} bases'
        draw.text((left + graph_width - draw.textlength(end_label, font=small_font), graph_top + graph_height + 6), end_label, fill='#6b7280', font=small_font)
    image.save(path, format='PNG', optimize=True, dpi=(144, 144))
    return {'width': width, 'height': height}


def _write_chromatogram_png(
    path: Path,
    reads: list[SangerRead],
    *,
    secondary_peak_ratio: float = 0.33,
    mixed_peak_min_quality: int = 20,
    page_note: str = '',
) -> dict[str, int]:
    """Write downsampled four-channel chromatograms with trim and mixture evidence."""
    try:
        from PIL import Image, ImageDraw
    except ImportError as exc:
        raise RuntimeError(
            'PNG report figures require Pillow. Install it with `conda install pillow`.'
        ) from exc
    width = 1280
    panel_height = CHROMATOGRAM_PANEL_HEIGHT
    left = 300
    right = 36
    graph_width = width - left - right
    height = CHROMATOGRAM_HEADER_HEIGHT + max(1, len(reads)) * panel_height
    colours = {'A': '#16a34a', 'C': '#2563eb', 'G': '#111827', 'T': '#dc2626'}
    image = Image.new('RGB', (width, height), '#ffffff')
    draw = ImageDraw.Draw(image)
    title_font = _load_report_font(21, bold=True)
    body_font = _load_report_font(12)
    label_font = _load_report_font(12)
    label_bold_font = _load_report_font(14, bold=True)
    small_font = _load_report_font(10)
    draw.text((24, 18), 'Paper Trail: chromatogram evidence', fill='#1f2937', font=title_font)
    draw.text((24, 49), 'A/C/G/T dye channels; green is retained sequence and orange marks high secondary peaks.', fill='#4b5563', font=body_font)
    if page_note:
        draw.text(
            (width - 24 - draw.textlength(page_note, font=body_font), 52),
            page_note, fill='#4b5563', font=body_font,
        )
    for base, colour in colours.items():
        x = 870 + ('ACGT'.index(base) * 76)
        draw.line((x, 29, x + 22, 29), fill=colour, width=2)
        draw.text((x + 28, 22), base, fill='#374151', font=label_font)
    for index, read in enumerate(reads):
        top = 78 + index * panel_height
        graph_top = top + 46
        graph_height = 92
        channels = read.trace_channels or {}
        trace_length = max([len(values) for values in channels.values()] + [0])
        draw.rounded_rectangle((18, top - 8, width - 18, top + panel_height - 16), radius=6, fill='#f9fafb', outline='#e5e7eb')
        draw.text((30, top + 5), _fit_report_text(draw, read.read_id, label_bold_font, width - 64), fill='#111827', font=label_bold_font)
        draw.text(
            (30, top + 29),
            _fit_report_text(draw, f'{read.sequence_id} | {read.primer} | mixed peaks {read.mixed_peak_count} ({read.mixed_peak_percent:.2f}%)', label_font, width - 64),
            fill='#4b5563', font=label_font,
        )
        if not read.trace_available or not trace_length:
            draw.text((left, graph_top + 40), 'Trace channels unavailable; base calls and Phred values only.', fill='#b45309', font=label_font)
            continue
        peaks = read.peak_locations or []
        retained_start = peaks[read.trim_start] if read.trim_start < len(peaks) else 0
        retained_end_index = min(max(read.trim_start, read.trim_end - 1), len(peaks) - 1)
        retained_end = peaks[retained_end_index] if peaks else trace_length
        x_start = left + (retained_start / max(1, trace_length - 1)) * graph_width
        x_end = left + (retained_end / max(1, trace_length - 1)) * graph_width
        draw.rectangle((round(x_start), graph_top, max(round(x_start) + 1, round(x_end)), graph_top + graph_height), fill='#dcfce7')
        pooled = sorted(value for values in channels.values() for value in values if value >= 0)
        scale = pooled[min(len(pooled) - 1, int(len(pooled) * 0.99))] if pooled else 1
        scale = max(1, scale)
        step = max(1, trace_length // 1000)
        for base in 'ACGT':
            values = channels.get(base, [])
            points: list[tuple[float, float]] = []
            for pos in range(0, len(values), step):
                x = left + (pos / max(1, trace_length - 1)) * graph_width
                y = graph_top + graph_height - (min(scale, values[pos]) / scale) * graph_height
                points.append((x, y))
            if points:
                draw.line(points, fill=colours[base], width=1)
        ratios = read.secondary_peak_ratios or []
        for base_index in range(read.trim_start, min(read.trim_end, len(peaks), len(ratios), len(read.raw_qualities))):
            if ratios[base_index] < secondary_peak_ratio or read.raw_qualities[base_index] < mixed_peak_min_quality:
                continue
            x = left + (peaks[base_index] / max(1, trace_length - 1)) * graph_width
            draw.ellipse((round(x) - 3, graph_top + 5, round(x) + 3, graph_top + 11), fill='#f97316')
        draw.text((left, graph_top + graph_height + 6), 'trace sample 1', fill='#6b7280', font=small_font)
        end_label = f'{trace_length} samples'
        draw.text((left + graph_width - draw.textlength(end_label, font=small_font), graph_top + graph_height + 6), end_label, fill='#6b7280', font=small_font)
    image.save(path, format='PNG', optimize=True, dpi=(144, 144))
    return {'width': width, 'height': height}


def _assembly_group_height(reads: list[SangerRead]) -> int:
    return 118 + max(1, len(reads)) * 30


def _write_assembly_overview_png(
    path: Path,
    groups: dict[str, list[SangerRead]],
    assembly_rows: list[dict[str, object]],
    *,
    page_note: str = '',
) -> dict[str, int]:
    try:
        from PIL import Image, ImageDraw
    except ImportError as exc:
        raise RuntimeError(
            'PNG report figures require Pillow. Install it with `conda install pillow`.'
        ) from exc
    width = 1280
    left = 390
    right = 40
    graph_width = width - left - right
    row_lookup = {str(row.get('SequenceID', '')): row for row in assembly_rows}
    group_heights = {sid: _assembly_group_height(reads) for sid, reads in groups.items()}
    height = max(170, ASSEMBLY_HEADER_HEIGHT + sum(group_heights.get(sid, 148) for sid in sorted(groups)))
    image = Image.new('RGB', (width, height), '#ffffff')
    draw = ImageDraw.Draw(image)
    title_font = _load_report_font(21, bold=True)
    body_font = _load_report_font(12)
    label_font = _load_report_font(11)
    label_bold_font = _load_report_font(14, bold=True)
    small_font = _load_report_font(10)
    draw.text((24, 18), 'Merge Meeting: Sanger isolate assembly', fill='#1f2937', font=title_font)
    draw.text((24, 49), 'Reads are aligned to consensus coordinates; colour shows whether each trace contributed.', fill='#4b5563', font=body_font)
    legend = (
        ('#374151', 'consensus'),
        ('#2563eb', 'contributing read'),
        ('#d97706', 'kept but not used'),
        ('#b91c1c', 'failed QC'),
    )
    legend_x = 24
    for colour, label in legend:
        draw.rounded_rectangle((legend_x, 73, legend_x + 18, 82), radius=2, fill=colour)
        draw.text((legend_x + 24, 69), label, fill='#4b5563', font=small_font)
        legend_x += 30 + round(draw.textlength(label, font=small_font))
    if page_note:
        draw.text(
            (width - 24 - draw.textlength(page_note, font=body_font), 22),
            page_note, fill='#4b5563', font=body_font,
        )
    y = ASSEMBLY_HEADER_HEIGHT
    if not groups:
        draw.text((24, 92), 'No read groups parsed.', fill='#6b7280', font=body_font)
    for sequence_id in sorted(groups):
        reads = sorted(groups[sequence_id], key=lambda r: (r.status != 'kept', r.direction != 'forward', r.read_id))
        row = row_lookup.get(sequence_id, {})
        group_height = group_heights.get(sequence_id, 100)
        output_length = int(row.get('OutputLength') or 0)
        max_len = max([output_length] + [len(_read_sequence_for_placement(r)) for r in reads] + [1])
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
        draw.rounded_rectangle((18, y, width - 18, y + group_height - 10), radius=6, fill='#f9fafb', outline='#d1d5db')
        draw.text((30, y + 12), _fit_report_text(draw, sequence_id, label_bold_font, width - 64), fill='#111827', font=label_bold_font)
        draw.text((30, y + 36), _fit_report_text(draw, details, label_font, width - 64), fill=status_colour, font=label_font)
        draw.text((30, y + 65), f'consensus | {output_length} bp', fill='#4b5563', font=label_font)
        draw.rounded_rectangle((left, y + 65, left + graph_width, y + 77), radius=2, fill='#e5e7eb')
        consensus_width = graph_width if output_length else 1
        draw.rounded_rectangle((left, y + 65, left + consensus_width, y + 77), radius=2, fill='#374151')
        if output_length:
            draw.text((left, y + 80), '1', fill='#6b7280', font=small_font)
            end_label = f'{output_length} bp'
            draw.text((left + graph_width - draw.textlength(end_label, font=small_font), y + 80), end_label, fill='#6b7280', font=small_font)
        placements = row.get('_ReadPlacements') or {}
        read_y = y + 104
        for read in reads:
            placement = placements.get(read.read_id, {})
            read_len = int(placement.get('read_length') or len(_read_sequence_for_placement(read)))
            contributes = bool(placement.get('contributes'))
            colour = '#2563eb' if contributes else ('#d97706' if read.status == 'kept' else '#b91c1c')
            start = placement.get('consensus_start')
            end = placement.get('consensus_end')
            coordinate = f'{int(start) + 1}-{int(end)}' if start is not None and end is not None else 'not aligned'
            role = 'used' if contributes else 'not used'
            direction = 'R' if read.direction == 'reverse' else ('F' if read.direction == 'forward' else '?')
            primer = '' if str(read.primer).lower() in {'', 'unknown'} else f'{read.primer}/'
            label = (
                f'{read.read_id} | {role} | cons {coordinate} | '
                f'{primer}{direction} | {read.status} | {read_len} bp'
            )
            draw.text((30, read_y - 1), _fit_report_text(draw, label, small_font, left - 44), fill='#4b5563', font=small_font)
            draw.rounded_rectangle((left, read_y, left + graph_width, read_y + 11), radius=2, fill='#e5e7eb')
            segments = list(placement.get('segments') or [])
            if output_length and segments:
                for consensus_start, consensus_end, _read_start, _read_end in segments:
                    x_start = left + (consensus_start / output_length) * graph_width
                    x_end = left + (consensus_end / output_length) * graph_width
                    draw.rounded_rectangle((round(x_start), read_y, max(round(x_start) + 1, round(x_end)), read_y + 11), radius=2, fill=colour)
            else:
                bar_width = (read_len / max_len) * graph_width
                draw.rounded_rectangle((left, read_y, left + max(1, round(bar_width)), read_y + 11), radius=2, fill=colour)
            read_y += 30
        y += group_height
    image.save(path, format='PNG', optimize=True, dpi=(144, 144))
    return {'width': width, 'height': height}


def _paginate_fixed_height(
    records: list,
    *,
    header_height: int,
    panel_height: int,
    max_image_height: int,
) -> list[tuple[int, int, list]]:
    capacity = max(1, (max_image_height - header_height) // panel_height)
    if not records:
        return [(0, 0, [])]
    return [
        (start, min(len(records), start + capacity), records[start:start + capacity])
        for start in range(0, len(records), capacity)
    ]


def _write_paginated_read_visuals(
    visual_root: Path,
    reads: list[SangerRead],
    *,
    max_image_height: int,
    secondary_peak_ratio: float,
    mixed_peak_min_quality: int,
) -> tuple[list[str], list[str], list[dict[str, object]]]:
    report_specs = (
        (
            'read_error_profiles', READ_ERROR_HEADER_HEIGHT, READ_ERROR_PANEL_HEIGHT,
            lambda path, page_reads, note: _write_read_error_profile_png(
                path, page_reads, page_note=note,
            ),
        ),
        (
            'trace_chromatograms', CHROMATOGRAM_HEADER_HEIGHT, CHROMATOGRAM_PANEL_HEIGHT,
            lambda path, page_reads, note: _write_chromatogram_png(
                path,
                page_reads,
                secondary_peak_ratio=secondary_peak_ratio,
                mixed_peak_min_quality=mixed_peak_min_quality,
                page_note=note,
            ),
        ),
    )
    output_paths: dict[str, list[str]] = {}
    manifest_rows: list[dict[str, object]] = []
    for report_name, header_height, panel_height, writer in report_specs:
        report_dir = visual_root / report_name
        report_dir.mkdir(parents=True, exist_ok=True)
        for stale_page in report_dir.glob(f'{report_name}_page_*.png'):
            stale_page.unlink()
        pages = _paginate_fixed_height(
            reads,
            header_height=header_height,
            panel_height=panel_height,
            max_image_height=max_image_height,
        )
        output_paths[report_name] = []
        for page_number, (start, end, page_reads) in enumerate(pages, start=1):
            filename = f'{report_name}_page_{page_number:03d}.png'
            path = report_dir / filename
            note = f'Page {page_number} of {len(pages)}'
            dimensions = writer(path, page_reads, note)
            output_paths[report_name].append(str(path))
            manifest_rows.append({
                'Report': report_name,
                'Page': page_number,
                'PageCount': len(pages),
                'File': str(path.relative_to(visual_root.parent)),
                'FirstRecord': page_reads[0].read_id if page_reads else 'NA',
                'LastRecord': page_reads[-1].read_id if page_reads else 'NA',
                'RecordCount': len(page_reads),
                'WidthPixels': dimensions['width'],
                'HeightPixels': dimensions['height'],
                'MaxHeightPixels': max_image_height,
                'PageNotes': '',
            })
    return (
        output_paths['read_error_profiles'],
        output_paths['trace_chromatograms'],
        manifest_rows,
    )


def _assembly_pages(
    groups: dict[str, list[SangerRead]],
    *,
    max_image_height: int,
) -> list[dict[str, object]]:
    body_height = max_image_height - ASSEMBLY_HEADER_HEIGHT
    max_reads_per_chunk = max(1, (body_height - 118) // 30)
    entries: list[dict[str, object]] = []
    for sequence_id in sorted(groups):
        reads = sorted(
            groups[sequence_id],
            key=lambda read: (read.status != 'kept', read.direction != 'forward', read.read_id),
        )
        chunks = [
            reads[start:start + max_reads_per_chunk]
            for start in range(0, len(reads), max_reads_per_chunk)
        ] or [[]]
        for chunk_number, chunk in enumerate(chunks, start=1):
            entries.append({
                'sequence_id': sequence_id,
                'reads': chunk,
                'height': _assembly_group_height(chunk),
                'chunk_number': chunk_number,
                'chunk_count': len(chunks),
            })

    if not entries:
        return [{'groups': {}, 'sequence_ids': [], 'continuations': []}]
    pages: list[dict[str, object]] = []
    current_groups: dict[str, list[SangerRead]] = {}
    current_ids: list[str] = []
    current_continuations: list[str] = []
    current_height = ASSEMBLY_HEADER_HEIGHT
    for entry in entries:
        entry_height = int(entry['height'])
        if current_groups and (
            current_height + entry_height > max_image_height
            or str(entry['sequence_id']) in current_groups
        ):
            pages.append({
                'groups': current_groups,
                'sequence_ids': current_ids,
                'continuations': current_continuations,
            })
            current_groups = {}
            current_ids = []
            current_continuations = []
            current_height = ASSEMBLY_HEADER_HEIGHT
        sequence_id = str(entry['sequence_id'])
        current_groups[sequence_id] = list(entry['reads'])
        current_ids.append(sequence_id)
        if int(entry['chunk_count']) > 1:
            current_continuations.append(
                f'{sequence_id} read page {entry["chunk_number"]}/{entry["chunk_count"]}'
            )
        current_height += entry_height
    if current_groups:
        pages.append({
            'groups': current_groups,
            'sequence_ids': current_ids,
            'continuations': current_continuations,
        })
    return pages


def _write_paginated_assembly_visuals(
    visual_root: Path,
    groups: dict[str, list[SangerRead]],
    assembly_rows: list[dict[str, object]],
    *,
    max_image_height: int,
) -> tuple[list[str], list[dict[str, object]]]:
    report_name = 'assembly_overviews'
    report_dir = visual_root / report_name
    report_dir.mkdir(parents=True, exist_ok=True)
    for stale_page in report_dir.glob('assembly_overview_page_*.png'):
        stale_page.unlink()
    pages = _assembly_pages(groups, max_image_height=max_image_height)
    output_paths: list[str] = []
    manifest_rows: list[dict[str, object]] = []
    for page_number, page in enumerate(pages, start=1):
        path = report_dir / f'assembly_overview_page_{page_number:03d}.png'
        page_note = f'Page {page_number} of {len(pages)}'
        continuations = list(page['continuations'])
        if continuations:
            page_note += ' | continued isolate reads'
        dimensions = _write_assembly_overview_png(
            path,
            dict(page['groups']),
            assembly_rows,
            page_note=page_note,
        )
        sequence_ids = list(page['sequence_ids'])
        output_paths.append(str(path))
        manifest_rows.append({
            'Report': report_name,
            'Page': page_number,
            'PageCount': len(pages),
            'File': str(path.relative_to(visual_root.parent)),
            'FirstRecord': sequence_ids[0] if sequence_ids else 'NA',
            'LastRecord': sequence_ids[-1] if sequence_ids else 'NA',
            'RecordCount': len(set(sequence_ids)),
            'WidthPixels': dimensions['width'],
            'HeightPixels': dimensions['height'],
            'MaxHeightPixels': max_image_height,
            'PageNotes': '; '.join(continuations),
        })
    return output_paths, manifest_rows


def _write_visual_report_manifest(path: Path, rows: list[dict[str, object]]) -> None:
    fields = [
        'Report', 'Page', 'PageCount', 'File', 'FirstRecord', 'LastRecord',
        'RecordCount', 'WidthPixels', 'HeightPixels', 'MaxHeightPixels',
        'PageNotes',
    ]
    with open(path, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def _screen_primer_read_taxonomy(
    reads: list[SangerRead],
    trimmed_fasta: Path,
    outdir: Path,
    *,
    ref_fasta: Optional[str | Path],
    taxa_tsv: Optional[str | Path],
    threads: int,
) -> tuple[dict[str, bool], Optional[str]]:
    """Flag isolate groups whose independently classified primer reads disagree."""
    if not ref_fasta:
        return {}, None
    from branchmanager.pipeline import classify
    from branchmanager.taxonomy import parse_taxon_string

    screen_dir = outdir / 'read_taxonomy_screen'
    screen_dir.mkdir(parents=True, exist_ok=True)
    taxonomy_path = classify.run_classification(
        str(trimmed_fasta), str(screen_dir), ref_fasta=str(ref_fasta),
        taxa_tsv=str(taxa_tsv) if taxa_tsv else None, threads=threads,
    )
    by_read = {}
    with open(taxonomy_path) as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            by_read[str(row.get('ID') or '')] = row
    for read in reads:
        row = by_read.get(read.read_id, {})
        read.taxonomy_screen = str(row.get('Taxon') or 'NA')

    conflicts = {}
    groups: dict[str, list[SangerRead]] = {}
    for read in reads:
        groups.setdefault(read.sequence_id, []).append(read)
    report_path = outdir / 'read_taxonomy_concordance.tsv'
    with open(report_path, 'w') as handle:
        handle.write('SequenceID\tReadIDs\tTaxonomies\tComparedRank\tConcordant\tReason\n')
        for sequence_id, group in sorted(groups.items()):
            resolved = []
            for read in group:
                parsed = parse_taxon_string(read.taxonomy_screen)
                rank = str(parsed.get('g') or parsed.get('f') or '').strip()
                if rank:
                    resolved.append((read.read_id, rank, read.taxonomy_screen))
            distinct = {rank for _read_id, rank, _taxonomy in resolved}
            conflict = len(distinct) > 1
            conflicts[sequence_id] = conflict
            reason = 'different genus/family assignments across primer reads' if conflict else (
                'concordant' if resolved else 'insufficient classified primer reads'
            )
            handle.write(
                f'{_clean_field(sequence_id)}\t'
                f'{_clean_field(";".join(read.read_id for read in group))}\t'
                f'{_clean_field(" | ".join(read.taxonomy_screen for read in group))}\t'
                f'genus_then_family\t{"no" if conflict else "yes"}\t{reason}\n'
            )
    return conflicts, str(report_path)


def _validate_qc_policy(
    *,
    min_quality: int,
    min_read_length: int,
    min_final_length: int,
    min_mean_quality: float,
    mask_quality: int,
    max_read_expected_errors: float,
    max_output_expected_errors: float,
    warn_n_percent: float,
    max_n_percent: float,
    warn_internal_low_quality_run: int,
    max_internal_low_quality_run: int,
    max_conflict_density: float,
    secondary_peak_ratio: float,
    max_mixed_peak_percent: float,
    mixed_peak_min_quality: int,
    quality_difference: int,
    min_overlap: int,
    min_overlap_identity: float,
) -> None:
    positive_lengths = {
        'min_read_length': min_read_length,
        'min_final_length': min_final_length,
        'min_overlap': min_overlap,
    }
    for name, value in positive_lengths.items():
        if int(value) < 1:
            raise ValueError(f'{name} must be at least 1')
    non_negative = {
        'min_quality': min_quality,
        'min_mean_quality': min_mean_quality,
        'mask_quality': mask_quality,
        'max_read_expected_errors': max_read_expected_errors,
        'max_output_expected_errors': max_output_expected_errors,
        'warn_internal_low_quality_run': warn_internal_low_quality_run,
        'max_internal_low_quality_run': max_internal_low_quality_run,
        'max_conflict_density': max_conflict_density,
        'mixed_peak_min_quality': mixed_peak_min_quality,
        'quality_difference': quality_difference,
    }
    for name, value in non_negative.items():
        if float(value) < 0:
            raise ValueError(f'{name} cannot be negative')
    for name, value in {
        'warn_n_percent': warn_n_percent,
        'max_n_percent': max_n_percent,
        'max_mixed_peak_percent': max_mixed_peak_percent,
    }.items():
        if not 0 <= float(value) <= 100:
            raise ValueError(f'{name} must be between 0 and 100')
    if float(warn_n_percent) > float(max_n_percent):
        raise ValueError('warn_n_percent cannot exceed max_n_percent')
    if int(warn_internal_low_quality_run) > int(max_internal_low_quality_run):
        raise ValueError(
            'warn_internal_low_quality_run cannot exceed max_internal_low_quality_run'
        )
    if float(secondary_peak_ratio) <= 0:
        raise ValueError('secondary_peak_ratio must be greater than 0')
    if not 0 <= float(min_overlap_identity) <= 1:
        raise ValueError('min_overlap_identity must be between 0 and 1')


def run_paper_trail(
    inputs: Iterable[str | Path],
    outdir: str | Path,
    *,
    sample_map: Optional[str | Path] = None,
    primers: Iterable[str] = DEFAULT_PRIMERS,
    primer_sequences: Optional[dict[str, str]] = None,
    trim_primers: bool = True,
    min_quality: int = DEFAULT_QC_POLICY.min_quality,
    min_length: int = DEFAULT_QC_POLICY.min_final_length,
    min_read_length: Optional[int] = None,
    min_mean_quality: float = DEFAULT_QC_POLICY.min_mean_quality,
    mask_quality: int = DEFAULT_QC_POLICY.mask_quality,
    max_read_expected_errors: float = DEFAULT_QC_POLICY.max_read_expected_errors,
    max_output_expected_errors: float = DEFAULT_QC_POLICY.max_output_expected_errors,
    warn_n_percent: float = DEFAULT_QC_POLICY.warn_n_percent,
    max_n_percent: float = DEFAULT_QC_POLICY.max_n_percent,
    warn_internal_low_quality_run: int = DEFAULT_QC_POLICY.warn_internal_low_quality_run,
    max_internal_low_quality_run: int = DEFAULT_QC_POLICY.max_internal_low_quality_run,
    max_conflict_density: float = DEFAULT_QC_POLICY.max_conflict_density,
    secondary_peak_ratio: float = DEFAULT_QC_POLICY.secondary_peak_ratio,
    max_mixed_peak_percent: float = DEFAULT_QC_POLICY.max_mixed_peak_percent,
    mixed_peak_min_quality: int = DEFAULT_QC_POLICY.mixed_peak_min_quality,
    quality_difference: int = DEFAULT_QC_POLICY.quality_difference,
    allow_missing_quality: bool = False,
    min_overlap: int = DEFAULT_QC_POLICY.min_overlap,
    min_overlap_identity: float = DEFAULT_QC_POLICY.min_overlap_identity,
    screen_ref: Optional[str | Path] = None,
    screen_taxa: Optional[str | Path] = None,
    threads: int = 4,
    assemble: bool = True,
    recursive: bool = True,
    max_report_image_height: int = DEFAULT_MAX_REPORT_IMAGE_HEIGHT,
) -> dict[str, object]:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    final_min_length = int(min_length)
    read_min_length = int(
        min_read_length
        if min_read_length is not None
        else min(DEFAULT_QC_POLICY.min_read_length, final_min_length)
    )
    _validate_qc_policy(
        min_quality=min_quality,
        min_read_length=read_min_length,
        min_final_length=final_min_length,
        min_mean_quality=min_mean_quality,
        mask_quality=mask_quality,
        max_read_expected_errors=max_read_expected_errors,
        max_output_expected_errors=max_output_expected_errors,
        warn_n_percent=warn_n_percent,
        max_n_percent=max_n_percent,
        warn_internal_low_quality_run=warn_internal_low_quality_run,
        max_internal_low_quality_run=max_internal_low_quality_run,
        max_conflict_density=max_conflict_density,
        secondary_peak_ratio=secondary_peak_ratio,
        max_mixed_peak_percent=max_mixed_peak_percent,
        mixed_peak_min_quality=mixed_peak_min_quality,
        quality_difference=quality_difference,
        min_overlap=min_overlap,
        min_overlap_identity=min_overlap_identity,
    )
    max_report_image_height = int(max_report_image_height)
    if max_report_image_height < MIN_REPORT_IMAGE_HEIGHT:
        raise ValueError(
            f'max_report_image_height must be at least {MIN_REPORT_IMAGE_HEIGHT} pixels'
        )
    default_mode = 'assemble' if assemble else 'best_read'
    effective_primer_sequences = dict(DEFAULT_PRIMER_SEQUENCES)
    for name, sequence in (primer_sequences or {}).items():
        effective_primer_sequences[str(name).upper()] = str(sequence).upper()
    metadata, listed_files = _load_sample_map(sample_map, primers)
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
        for record_id, sequence, qualities, quality_available, quality_source, trace_details in _read_records_from_file(path):
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
                trace_available=bool(trace_details.get('trace_available')),
                trace_order=str(trace_details.get('trace_order') or ''),
                peak_locations=list(trace_details.get('peak_locations') or []),
                trace_channels=dict(trace_details.get('trace_channels') or {}),
                secondary_peak_ratios=list(trace_details.get('secondary_peak_ratios') or []),
                secondary_peak_bases=list(trace_details.get('secondary_peak_bases') or []),
            )
            raw_records.append((read.read_id, read.raw_sequence))
            if quality_available:
                trimmed, trimmed_q, start, end = trim_by_quality(
                    sequence,
                    qualities,
                    min_quality=min_quality,
                )
            else:
                start = 0
                end = len(sequence)
                while start < end and sequence[start] not in {'A', 'C', 'G', 'T'}:
                    start += 1
                while end > start and sequence[end - 1] not in {'A', 'C', 'G', 'T'}:
                    end -= 1
                trimmed = sequence[start:end]
                trimmed_q = list(qualities[start:end])
            if trim_primers:
                trimmed, trimmed_q, primer_bases = trim_leading_primer(
                    trimmed,
                    trimmed_q,
                    primer,
                    primer_sequences=effective_primer_sequences,
                )
                read.primer_trimmed_bases = primer_bases
                start += primer_bases
            read.trimmed_sequence = trimmed
            read.trimmed_qualities = trimmed_q
            read.trim_start = start
            read.trim_end = end
            ratios = read.secondary_peak_ratios or []
            retained_positions = [
                idx for idx in range(start, min(end, len(ratios), len(qualities)))
                if qualities[idx] >= int(mixed_peak_min_quality)
            ]
            read.mixed_peak_count = sum(
                ratios[idx] >= float(secondary_peak_ratio) for idx in retained_positions
            )
            read.mixed_peak_percent = (
                (read.mixed_peak_count / len(retained_positions)) * 100.0
                if retained_positions else 0.0
            )
            if quality_available:
                masked, masked_q, masked_count, lowq_run = mask_low_quality_bases(
                    trimmed,
                    trimmed_q,
                    mask_quality=mask_quality,
                )
            else:
                masked = ''.join(base if base in {'A', 'C', 'G', 'T'} else 'N' for base in trimmed)
                masked_q = list(trimmed_q)
                masked_count = sum(base == 'N' for base in masked)
                lowq_run = _longest_true_run(base == 'N' for base in masked)
            read.masked_sequence = masked
            read.masked_qualities = masked_q
            read.masked_bases = masked_count
            read.longest_low_quality_run = lowq_run
            _classify_read_qc(
                read,
                read_min_length=read_min_length,
                min_mean_quality=float(min_mean_quality),
                max_read_expected_errors=float(max_read_expected_errors),
                warn_n_percent=float(warn_n_percent),
                max_n_percent=float(max_n_percent),
                warn_internal_low_quality_run=int(warn_internal_low_quality_run),
                max_internal_low_quality_run=int(max_internal_low_quality_run),
                allow_missing_quality=bool(allow_missing_quality),
                max_mixed_peak_percent=float(max_mixed_peak_percent),
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
    assembly_placements_tsv = out / 'assembly_read_placements.tsv'
    recommendations_tsv = out / 'resequence_recommendations.tsv'
    marker_review_template_tsv = out / 'marker_review_template.tsv'
    qc_policy_tsv = out / 'paper_trail_qc_policy.tsv'
    visual_root = out / 'visual_reports'
    visual_manifest_tsv = visual_root / 'visual_report_manifest.tsv'
    for obsolete_name in (
        'read_error_profiles.png', 'trace_chromatograms.png', 'assembly_overview.png',
    ):
        obsolete_path = out / obsolete_name
        if obsolete_path.exists():
            obsolete_path.unlink()
    summary_txt = out / 'paper_trail_summary.txt'
    failed_qc_dir = out / 'failed_qc_sequences'
    failed_final_fasta = failed_qc_dir / 'failed_final_sequences.fasta'
    failed_read_fasta = failed_qc_dir / 'failed_read_sequences.fasta'
    failed_manifest_tsv = failed_qc_dir / 'failed_qc_manifest.tsv'
    failed_read_manifest_tsv = failed_qc_dir / 'failed_read_manifest.tsv'
    failed_qc_guide = failed_qc_dir / 'failed_qc_guide.txt'

    write_fasta(raw_records, str(raw_fasta))
    write_fasta(trimmed_records, str(trimmed_fasta))
    taxonomy_conflicts, taxonomy_screen_path = _screen_primer_read_taxonomy(
        reads,
        trimmed_fasta,
        out,
        ref_fasta=screen_ref,
        taxa_tsv=screen_taxa,
        threads=int(threads or 1),
    )
    _write_per_base_error_tsv(per_base_error_tsv, reads)
    read_error_pngs, chromatogram_pngs, visual_manifest_rows = _write_paginated_read_visuals(
        visual_root,
        reads,
        max_image_height=max_report_image_height,
        secondary_peak_ratio=float(secondary_peak_ratio),
        mixed_peak_min_quality=int(mixed_peak_min_quality),
    )

    with open(read_qc_tsv, 'w') as handle:
        handle.write(
            'ReadID\tSequenceID\tPrimer\tDirection\tSourceFile\tRawLength\tTrimmedLength\t'
            'MaskedLength\tTrimStart\tTrimEnd\tLeftTrimmedBases\tRightTrimmedBases\t'
            'MeanRawQuality\tMeanTrimmedQuality\tMeanMaskedQuality\t'
            'RawExpectedErrors\tTrimmedExpectedErrors\tMaskedExpectedErrors\t'
            'MeanRawErrorProbability\tMeanTrimmedErrorProbability\tMaskedNPercent\t'
            'MaskedBases\tLongestLowQualityRun\tQualityAvailable\tQualitySource\t'
            'TraceAvailable\tTraceOrder\tMixedPeakCount\tMixedPeakPercent\tPrimerTrimmedBases\t'
            'Status\tQCClass\tReasons\tProcessingMode\tProcessingModeExplicit\tTaxonomyScreen\n'
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
                f'{_clean_field(read.quality_source)}\t{"yes" if read.trace_available else "no"}\t'
                f'{_clean_field(read.trace_order)}\t{read.mixed_peak_count}\t{read.mixed_peak_percent:.3f}\t'
                f'{read.primer_trimmed_bases}\t{read.status}\t{read.qc_class}\t'
                f'{_clean_field(read.warning)}\t{read.processing_mode}\t'
                f'{"yes" if read.processing_mode_explicit else "no"}\t{_clean_field(read.taxonomy_screen)}\n'
            )

    groups: dict[str, list[SangerRead]] = {}
    for read in reads:
        groups.setdefault(read.sequence_id, []).append(read)

    assembled_records = []
    failed_final_records = []
    failed_read_records = []
    failed_read_manifest_rows: list[dict[str, object]] = []
    failed_isolate_manifest_rows: list[dict[str, object]] = []
    for read in reads:
        if read.qc_class != 'FAIL_QC':
            continue
        failed_seq = _review_sequence_for_failed_read(read)
        if failed_seq:
            failed_header = (
                f'{_clean_fasta_token(read.sequence_id)}'
                f'|read={_clean_fasta_token(read.read_id)}'
                f'|qc={_clean_fasta_token(read.qc_class)}'
                f'|status={_clean_fasta_token(read.status)}'
            )
            failed_read_records.append((failed_header, failed_seq))
        failed_read_manifest_rows.append({
            'SequenceID': read.sequence_id,
            'ReadID': read.read_id,
            'SourceFile': read.source_file,
            'QCClass': read.qc_class,
            'Recommendation': 'RESEQUENCE',
            'Status': read.status,
            'Reasons': read.warning,
            'OutputLength': len(failed_seq),
            'MeanQuality': _mean_quality(read.masked_qualities or read.trimmed_qualities or read.raw_qualities),
            'OutputExpectedErrors': _expected_errors(read.masked_qualities or read.trimmed_qualities or read.raw_qualities, failed_seq),
            'OutputNPercent': _n_percent(failed_seq),
            'ProcessingMode': read.processing_mode,
        })
    assembly_rows: list[dict[str, object]] = []
    recommendation_rows: list[dict[str, object]] = []
    with open(assembly_tsv, 'w') as handle:
        handle.write(
            'SequenceID\tStatus\tReadCount\tUsedReads\tOutputLength\tMeanQuality\t'
            'OutputExpectedErrors\tOutputNPercent\tConflictDensity\tQCClass\tRecommendation\t'
            'OverlapIdentity\tOverlapLength\tConflicts\tAmbiguousConflicts\tUnmergedReads\t'
            'ReadIDs\tKeptReadIDs\tUsedReadIDs\tMergeMethod\tProcessingMode\tSelectedReadID\t'
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
            stats['taxonomy_conflict'] = bool(taxonomy_conflicts.get(sequence_id, False))
            consensus_qualities = list(stats.get('qualities') or [])
            final_qc = _classify_output_qc(
                consensus,
                consensus_qualities,
                stats,
                group,
                final_min_length=final_min_length,
                min_mean_quality=float(min_mean_quality),
                max_output_expected_errors=float(max_output_expected_errors),
                warn_n_percent=float(warn_n_percent),
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
            elif consensus:
                failed_final_records.append((sequence_id, consensus))
            overlap_identity = stats['overlap_identity']
            if isinstance(overlap_identity, float):
                overlap_identity = f'{overlap_identity:.4f}'
            read_ids = ';'.join(_clean_field(read.read_id) for read in group)
            kept_read_ids = ';'.join(_clean_field(read.read_id) for read in group if read.status == 'kept')
            used_read_ids = [str(read_id) for read_id in stats.get('used_read_ids', [])]
            read_placements = _build_read_placements(consensus, group, used_read_ids)
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
                'UsedReadIDs': ';'.join(used_read_ids),
                'MergeMethod': stats.get('method', 'NA'),
                'ProcessingMode': processing_mode,
                'SelectedReadID': stats.get('selected_read_id', ''),
                'PassLengthQC': pass_length_qc if consensus else 'no',
                'Reasons': reasons,
                'SuggestedAction': final_qc.get('suggested_action', ''),
                'FailedReadIDs': final_qc.get('failed_read_ids', ''),
                'ModeWarning': mode_warning,
                '_ConsensusSequence': consensus,
                '_ReadPlacements': read_placements,
            }
            assembly_rows.append(row)
            recommendation_rows.append(row)
            if row['QCClass'] == 'FAIL_QC':
                failed_group_reads = [read for read in group if read.qc_class == 'FAIL_QC']
                recoverable_reads = [
                    read for read in group if _review_sequence_for_failed_read(read)
                ]
                best_recoverable = (
                    max(recoverable_reads, key=_read_quality_key)
                    if recoverable_reads else None
                )
                best_recoverable_sequence = (
                    _review_sequence_for_failed_read(best_recoverable)
                    if best_recoverable else ''
                )
                failed_isolate_manifest_rows.append({
                    'SequenceID': sequence_id,
                    'QCClass': row['QCClass'],
                    'Recommendation': row['Recommendation'],
                    'Status': row['Status'],
                    'Reasons': row['Reasons'],
                    'ReadCount': row['ReadCount'],
                    'KeptReadCount': len([read for read in group if read.status == 'kept']),
                    'FailedReadCount': len(failed_group_reads),
                    'ReadIDs': row['ReadIDs'],
                    'KeptReadIDs': row['KeptReadIDs'],
                    'FailedReadIDs': row['FailedReadIDs'],
                    'SourceFiles': ';'.join(sorted({read.source_file for read in group})),
                    'ReadFailureReasons': ' | '.join(
                        f'{read.read_id}[{read.warning}]' for read in failed_group_reads
                    ),
                    'OutputLength': row['OutputLength'],
                    'MeanQuality': row['MeanQuality'],
                    'OutputExpectedErrors': row['OutputExpectedErrors'],
                    'OutputNPercent': row['OutputNPercent'],
                    'ProcessingMode': row['ProcessingMode'],
                    'SelectedReadID': row['SelectedReadID'],
                    'MergeMethod': row['MergeMethod'],
                    'BestRecoverableReadID': best_recoverable.read_id if best_recoverable else '',
                    'BestRecoverableLength': len(best_recoverable_sequence),
                    'BestRecoverableNPercent': _n_percent(best_recoverable_sequence),
                })
            handle.write(
                f'{_clean_field(sequence_id)}\t{stats["status"]}\t{stats["read_count"]}\t{stats["used_reads"]}\t'
                f'{len(consensus)}\t{row["MeanQuality"]}\t{row["OutputExpectedErrors"]}\t'
                f'{row["OutputNPercent"]}\t{row["ConflictDensity"]}\t{row["QCClass"]}\t'
                f'{row["Recommendation"]}\t{overlap_identity}\t'
                f'{stats["overlap_length"]}\t{stats["conflicts"]}\t{stats.get("ambiguous_conflicts", "NA")}\t'
                f'{_clean_field(stats["unmerged_reads"])}\t{read_ids}\t{kept_read_ids}\t{_clean_field(row["UsedReadIDs"])}\t'
                f'{_clean_field(stats.get("method", "NA"))}\t{processing_mode}\t'
                f'{_clean_field(stats.get("selected_read_id", ""))}\t'
                f'{pass_length_qc if consensus else "no"}\t{_clean_field(reasons)}\t'
                f'{_clean_field(row["SuggestedAction"])}\t{_clean_field(row["FailedReadIDs"])}\t'
                f'{_clean_field(mode_warning)}\n'
            )

    with open(assembly_placements_tsv, 'w') as handle:
        handle.write(
            'SequenceID\tReadID\tPrimer\tDirection\tReadStatus\tContributesToConsensus\t'
            'ConsensusStart\tConsensusEnd\tAlignedBases\tReadLength\tAlignedSegments\n'
        )
        for row in assembly_rows:
            sequence_id = str(row.get('SequenceID', ''))
            reads_by_id = {read.read_id: read for read in groups.get(sequence_id, [])}
            for read_id, placement in (row.get('_ReadPlacements') or {}).items():
                read = reads_by_id.get(read_id)
                start = placement.get('consensus_start')
                end = placement.get('consensus_end')
                segments = ';'.join(
                    f'{consensus_start + 1}-{consensus_end}'
                    for consensus_start, consensus_end, _read_start, _read_end
                    in placement.get('segments', [])
                )
                handle.write(
                    f'{_clean_field(sequence_id)}\t{_clean_field(read_id)}\t'
                    f'{_clean_field(read.primer if read else "")}\t'
                    f'{_clean_field(read.direction if read else "")}\t'
                    f'{_clean_field(read.status if read else "")}\t'
                    f'{"yes" if placement.get("contributes") else "no"}\t'
                    f'{int(start) + 1 if start is not None else ""}\t'
                    f'{int(end) if end is not None else ""}\t'
                    f'{placement.get("aligned_bases", 0)}\t{placement.get("read_length", 0)}\t'
                    f'{segments}\n'
                )

    write_fasta(assembled_records, str(assembled_fasta))
    failed_qc_dir.mkdir(parents=True, exist_ok=True)
    write_fasta(failed_final_records, str(failed_final_fasta))
    write_fasta(failed_read_records, str(failed_read_fasta))
    with open(failed_manifest_tsv, 'w') as handle:
        handle.write(
            'SequenceID\tQCClass\tRecommendation\tStatus\tReasons\tReadCount\t'
            'KeptReadCount\tFailedReadCount\tReadIDs\tKeptReadIDs\tFailedReadIDs\t'
            'SourceFiles\tReadFailureReasons\tOutputLength\tMeanQuality\t'
            'OutputExpectedErrors\tOutputNPercent\tProcessingMode\tSelectedReadID\t'
            'MergeMethod\tBestRecoverableReadID\tBestRecoverableLength\t'
            'BestRecoverableNPercent\n'
        )
        for row in failed_isolate_manifest_rows:
            handle.write(
                f'{_clean_field(row.get("SequenceID", ""))}\t{_clean_field(row.get("QCClass", ""))}\t'
                f'{_clean_field(row.get("Recommendation", ""))}\t'
                f'{_clean_field(row.get("Status", ""))}\t{_clean_field(row.get("Reasons", ""))}\t'
                f'{row.get("ReadCount", 0)}\t{row.get("KeptReadCount", 0)}\t'
                f'{row.get("FailedReadCount", 0)}\t{_clean_field(row.get("ReadIDs", ""))}\t'
                f'{_clean_field(row.get("KeptReadIDs", ""))}\t'
                f'{_clean_field(row.get("FailedReadIDs", ""))}\t'
                f'{_clean_field(row.get("SourceFiles", ""))}\t'
                f'{_clean_field(row.get("ReadFailureReasons", ""))}\t'
                f'{row.get("OutputLength", 0)}\t{row.get("MeanQuality", 0.0)}\t'
                f'{row.get("OutputExpectedErrors", 0.0)}\t{row.get("OutputNPercent", 0.0)}\t'
                f'{_clean_field(row.get("ProcessingMode", ""))}\t'
                f'{_clean_field(row.get("SelectedReadID", ""))}\t{_clean_field(row.get("MergeMethod", ""))}\t'
                f'{_clean_field(row.get("BestRecoverableReadID", ""))}\t'
                f'{row.get("BestRecoverableLength", 0)}\t'
                f'{row.get("BestRecoverableNPercent", 0.0):.3f}\n'
            )
    with open(failed_read_manifest_tsv, 'w') as handle:
        handle.write(
            'SequenceID\tReadID\tSourceFile\tQCClass\tRecommendation\tStatus\tReasons\t'
            'RecoverySequenceLength\tMeanQuality\tExpectedErrors\tNPercent\tProcessingMode\n'
        )
        for row in failed_read_manifest_rows:
            handle.write(
                f'{_clean_field(row.get("SequenceID", ""))}\t{_clean_field(row.get("ReadID", ""))}\t'
                f'{_clean_field(row.get("SourceFile", ""))}\t{_clean_field(row.get("QCClass", ""))}\t'
                f'{_clean_field(row.get("Recommendation", ""))}\t{_clean_field(row.get("Status", ""))}\t'
                f'{_clean_field(row.get("Reasons", ""))}\t{row.get("OutputLength", 0)}\t'
                f'{float(row.get("MeanQuality", 0.0)):.2f}\t'
                f'{float(row.get("OutputExpectedErrors", 0.0)):.4f}\t'
                f'{float(row.get("OutputNPercent", 0.0)):.3f}\t'
                f'{_clean_field(row.get("ProcessingMode", ""))}\n'
            )
    failed_qc_guide.write_text(
        '\n'.join([
            'BranchManager failed-QC guide',
            '',
            'failed_qc_manifest.tsv has exactly one row per failed isolate/final marker.',
            'failed_read_manifest.tsv has one row per failed physical read and supplies the read-level evidence.',
            'failed_final_sequences.fasta and failed_read_sequences.fasta are recovery/review files; their records did not pass downstream QC.',
            '',
            'Common reason codes:',
            f'- trimmed_length_lt_{read_min_length}: fewer than {read_min_length} usable bases remained after Q{min_quality} Mott trimming.',
            f'- output_length_lt_{final_min_length}: a readable primer trace remained, but the final marker was shorter than {final_min_length} bp.',
            f'- n_percent_gt_{float(max_n_percent):g}: more than {float(max_n_percent):g}% of retained positions were ambiguous after Q{mask_quality} masking.',
            f'- mixed_peak_percent_gt_{float(max_mixed_peak_percent):g}: the chromatogram showed excessive secondary-peak evidence consistent with a mixed template or unresolved trace.',
            f'- mean_q_lt_{float(min_mean_quality):g}: retained base calls had mean Phred quality below {float(min_mean_quality):g}.',
            f'- expected_errors_gt_{float(max_read_expected_errors):g}: summed Phred error probabilities exceeded the read-level limit.',
            f'- internal_low_quality_run_gt_{int(max_internal_low_quality_run)}: a retained ambiguous/low-quality run exceeded {int(max_internal_low_quality_run)} bases.',
            '- failed_no_reads: no physical read passed read-level QC for final-marker construction.',
            '- primer_read_taxonomic_conflict: separate primer reads classified to different genus/family contexts.',
            '',
        ])
    )
    assembly_pngs, assembly_visual_rows = _write_paginated_assembly_visuals(
        visual_root,
        groups,
        assembly_rows,
        max_image_height=max_report_image_height,
    )
    visual_manifest_rows.extend(assembly_visual_rows)
    _write_visual_report_manifest(visual_manifest_tsv, visual_manifest_rows)
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
    with open(marker_review_template_tsv, 'w') as handle:
        handle.write(
            'sequence_id\tdecision\treviewer\tnotes\tqc_reasons\toutput_length\tn_percent\n'
        )
        for row in recommendation_rows:
            if row.get('Recommendation') != 'MANUAL_REVIEW':
                continue
            handle.write(
                f'{_clean_field(row.get("SequenceID", ""))}\t\t\t\t'
                f'{_clean_field(row.get("Reasons", ""))}\t{row.get("OutputLength", 0)}\t'
                f'{row.get("OutputNPercent", 0.0)}\n'
            )
    qc_policy_tsv.write_text(
        '\n'.join([
            'Metric\tValue\tMeaning',
            f'policy_version\t{PAPER_TRAIL_QC_POLICY_VERSION}\tVersioned BranchManager Paper Trail QC policy',
            'trim_algorithm\tmodified_mott\tMaximum-scoring trim interval calculated from Phred error probabilities',
            f'min_quality\t{min_quality}\tPhred cutoff used for Mott-style end trimming',
            f'mask_quality\t{mask_quality}\tBases below this Phred score inside the retained read are masked to N',
            f'min_read_length\t{read_min_length}\tMinimum retained read length before assembly or best-read selection',
            f'min_final_length\t{final_min_length}\tMinimum final isolate sequence length written to assembled.fasta',
            f'min_mean_quality\t{float(min_mean_quality):g}\tMinimum mean Phred score after masking for read/final QC',
            f'max_read_expected_errors\t{float(max_read_expected_errors):g}\tMaximum expected errors among called A/C/G/T bases per retained read',
            f'max_output_expected_errors\t{float(max_output_expected_errors):g}\tMaximum expected errors among called A/C/G/T bases in the final isolate sequence',
            f'warn_n_percent\t{float(warn_n_percent):g}\tPercent N above which a passing read/final sequence requires manual review',
            f'max_n_percent\t{float(max_n_percent):g}\tMaximum percent N allowed before read/final QC failure',
            f'warn_internal_low_quality_run\t{int(warn_internal_low_quality_run)}\tInternal low-quality/ambiguous run above which manual review is required',
            f'max_internal_low_quality_run\t{int(max_internal_low_quality_run)}\tLongest internal low-quality/ambiguous run allowed before read failure',
            f'max_conflict_density\t{float(max_conflict_density):g}\tMaximum overlap conflicts per 100 final bases before final failure',
            f'secondary_peak_ratio\t{float(secondary_peak_ratio):g}\tSecondary dye peak divided by called-base peak used to flag mixed-template evidence',
            f'mixed_peak_review_percent\t{max(1.0, float(max_mixed_peak_percent) / 2.0):g}\tMixed-peak percent above which a passing read requires manual review',
            f'max_mixed_peak_percent\t{float(max_mixed_peak_percent):g}\tMaximum retained high-quality bases with secondary-peak evidence',
            f'mixed_peak_min_quality\t{int(mixed_peak_min_quality)}\tMinimum called-base Phred score considered for mixed-peak screening',
            f'min_overlap\t{int(min_overlap)}\tMinimum compared bases required to assemble primer reads',
            f'min_overlap_identity\t{float(min_overlap_identity):g}\tMinimum overlap identity required to assemble primer reads',
            f'trim_primers\t{bool(trim_primers)}\tRemove a confidently observed known primer sequence at the read start',
            f'quality_difference\t{int(quality_difference)}\tMinimum Phred-scaled posterior odds required to resolve a conflicting overlap base',
            f'allow_missing_quality\t{bool(allow_missing_quality)}\tWhether AB1 reads missing PCON quality scores may pass QC',
            f'default_processing_mode\t{default_mode}\tFallback handling when a sample-map row does not specify assemble or best_read',
            f'max_report_image_height\t{max_report_image_height}\tMaximum PNG page height before visual reports are paginated',
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
            'BranchManager Paper Trail / Merge Meeting Summary',
            f'Input files: {len(files)}',
            f'Reads parsed: {len(reads)}',
            f'Minimum trimmed read length: {read_min_length}',
            f'Minimum final sequence length: {final_min_length}',
            f'Reads kept after trimming: {kept_reads}',
            f'Passing final sequences written: {assembled_count}',
            f'Failed final sequences retained: {len(failed_final_records)}',
            f'Failed read sequences retained: {len(failed_read_records)}',
            f'Final QC accept: {accept_count}',
            f'Final QC manual review: {review_count}',
            f'Final QC resequence: {resequence_count}',
            f'Output groups failing final QC: {filtered_outputs}',
            'Processing modes: ' + ', '.join(f'{mode}={count}' for mode, count in sorted(mode_counts.items())),
            f'Final FASTA: {assembled_fasta}',
            f'Read QC: {read_qc_tsv}',
            f'Per-base error profile: {per_base_error_tsv}',
            f'Assembly report: {assembly_tsv}',
            f'Assembly read placements: {assembly_placements_tsv}',
            f'Resequencing recommendations: {recommendations_tsv}',
            f'Marker review template: {marker_review_template_tsv}',
            f'QC policy: {qc_policy_tsv}',
            f'Failed QC final FASTA: {failed_final_fasta}',
            f'Failed QC read FASTA: {failed_read_fasta}',
            f'Failed isolate manifest: {failed_manifest_tsv}',
            f'Failed read manifest: {failed_read_manifest_tsv}',
            f'Failed QC guide: {failed_qc_guide}',
            f'Visual report maximum height: {max_report_image_height} pixels',
            f'Read profile visual pages: {len(read_error_pngs)} in {visual_root / "read_error_profiles"}',
            f'Chromatogram visual pages: {len(chromatogram_pngs)} in {visual_root / "trace_chromatograms"}',
            f'Primer-read taxonomy concordance: {taxonomy_screen_path or "not requested"}',
            f'Assembly visual pages: {len(assembly_pngs)} in {visual_root / "assembly_overviews"}',
            f'Visual report manifest: {visual_manifest_tsv}',
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
        'assembly_placements_tsv': str(assembly_placements_tsv),
        'recommendations_tsv': str(recommendations_tsv),
        'marker_review_template_tsv': str(marker_review_template_tsv),
        'qc_policy_tsv': str(qc_policy_tsv),
        'failed_qc_dir': str(failed_qc_dir),
        'failed_final_fasta': str(failed_final_fasta),
        'failed_read_fasta': str(failed_read_fasta),
        'failed_manifest_tsv': str(failed_manifest_tsv),
        'failed_read_manifest_tsv': str(failed_read_manifest_tsv),
        'failed_qc_guide': str(failed_qc_guide),
        'read_error_pngs': read_error_pngs,
        'chromatogram_pngs': chromatogram_pngs,
        'taxonomy_screen_tsv': str(taxonomy_screen_path or ''),
        'assembly_pngs': assembly_pngs,
        'visual_manifest_tsv': str(visual_manifest_tsv),
        'summary': str(summary_txt),
    }
