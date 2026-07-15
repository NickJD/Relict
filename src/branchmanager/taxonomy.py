from __future__ import annotations

import re
from typing import Dict, Optional

RANK_ALIASES = {
    'domain': 'd',
    'superkingdom': 'd',
    'kingdom': 'k',
    'phylum': 'p',
    'class': 'c',
    'order': 'o',
    'family': 'f',
    'genus': 'g',
    'species': 's',
}

DOMAIN_QUERY_ALIASES = {
    'bacterial': 'bacteria',
    'bacterium': 'bacteria',
    'bacteria': 'bacteria',
    'archaeal': 'archaea',
    'archaea': 'archaea',
    'fungal': 'fungi',
    'fungus': 'fungi',
    'fungi': 'fungi',
    'eukaryotic': 'eukaryota',
    'eukarya': 'eukaryota',
    'eukaryota': 'eukaryota',
}


def _clean_taxon_value(value: str) -> str:
    if value is None:
        return ''
    value = str(value).strip().strip('"\'')
    value = re.sub(r'\s+', ' ', value)
    return value


def normalise_taxon_name(value: str) -> str:
    value = _clean_taxon_value(value)
    if not value:
        return ''
    value = value.lower()
    value = re.sub(r'[^a-z0-9]+', '', value)
    return value


def normalise_domain_query(value: str) -> str:
    norm = normalise_taxon_name(value)
    return DOMAIN_QUERY_ALIASES.get(norm, norm)


def parse_taxon_string(taxon: str) -> Dict[str, str]:
    """Parse lineage strings like d__Bacteria; p__Firmicutes into rank->name.

    Supports both short prefixes (d__/p__) and long prefixes (domain__/phylum__).
    Returns short one-letter rank keys where possible and also preserves unknown
    prefixes verbatim.
    """
    result: Dict[str, str] = {}
    if taxon is None:
        return result
    for raw_part in str(taxon).split(';'):
        part = raw_part.strip()
        if not part or '__' not in part:
            continue
        raw_rank, raw_value = part.split('__', 1)
        rank = raw_rank.strip().lower()
        value = _clean_taxon_value(raw_value)
        if not value:
            continue
        rank = RANK_ALIASES.get(rank, rank[:1] if rank else rank)
        result[rank] = value
    return result


def get_rank_value(taxon: str, rank: str) -> Optional[str]:
    parsed = parse_taxon_string(taxon)
    key = RANK_ALIASES.get(rank.lower(), rank.lower()[:1])
    return parsed.get(key)


def get_domain_or_kingdom(taxon: str) -> Optional[str]:
    parsed = parse_taxon_string(taxon)
    return parsed.get('d') or parsed.get('k')


def taxonomy_matches_kingdom(taxon: Optional[str], kingdom: str) -> bool:
    """Return True if lineage belongs to the requested kingdom/domain.

    Matching is rank-aware and prefers the domain/kingdom rank rather than a
    loose substring search across the whole lineage.
    """
    if not kingdom:
        return True
    wanted = normalise_domain_query(kingdom)
    if not wanted:
        return True
    parsed = parse_taxon_string(taxon)
    domain = normalise_domain_query(parsed.get('d', ''))
    kingdom_value = normalise_domain_query(parsed.get('k', ''))
    candidates = [value for value in (domain, kingdom_value) if value]
    if candidates:
        return wanted in candidates
    # fallback for malformed lineage strings with no explicit rank prefixes
    raw = normalise_domain_query(taxon)
    if not raw:
        return False
    return raw == wanted


def canonicalise_sequence_id(header: str) -> Optional[str]:
    """Conservative sequence-id canonicalisation.

    Keeps biologically meaningful suffixes such as `_1` / `_2` (ASV/OTU labels),
    while stripping transport/coordinate artefacts from FASTA headers.
    """
    if header is None:
        return None
    token = str(header).strip().split()[0]
    if not token:
        return None
    if '|' in token:
        token = token.split('|')[-1]
    if '#' in token:
        token = token.split('#', 1)[0]
    token = re.sub(r':\d+-\d+(?:\([^)]*\))?$', '', token)
    token = re.sub(r':\d+-\d+$', '', token)
    token = token.strip('()[]{}')
    return token or None


def reference_lookup_keys(header: Optional[str]) -> list[str]:
    """Return resilient lookup keys for reference/taxonomy ID matching."""
    if header is None:
        return []
    raw = str(header).strip()
    if not raw:
        return []
    keys = []
    first = raw.split()[0]
    for candidate in (raw, first, canonicalise_sequence_id(raw), canonicalise_sequence_id(first)):
        if candidate and candidate not in keys:
            keys.append(candidate)
    if '|' in first:
        last = first.split('|')[-1]
        if last and last not in keys:
            keys.append(last)
        canon_last = canonicalise_sequence_id(last)
        if canon_last and canon_last not in keys:
            keys.append(canon_last)
    return keys


def parse_reference_header_taxonomy(header: Optional[str]) -> tuple[Optional[str], Optional[str]]:
    """Extract (reference_id, taxonomy) from GTDB-style FASTA headers.

    Example header:
      RS_GCF_... d__Bacteria;p__... [locus_tag=...] ...
    """
    if header is None:
        return None, None
    parts = str(header).strip().split(None, 1)
    ref_id = parts[0] if parts else None
    remainder = parts[1] if len(parts) > 1 else ''
    if not remainder:
        return ref_id, None
    lineage = remainder.split('[', 1)[0].strip()
    if 'd__' not in lineage:
        return ref_id, None
    lineage = re.sub(r'\s*;\s*', ';', lineage)
    return ref_id, lineage or None
