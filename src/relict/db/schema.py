SCHEMA = """
CREATE TABLE IF NOT EXISTS sequences (
    id TEXT PRIMARY KEY,
    sequence TEXT,
    length INTEGER,
    dataset TEXT
);

CREATE TABLE IF NOT EXISTS taxonomy (
    id TEXT,
    taxonomy TEXT,
    confidence REAL,
    dataset TEXT
);

CREATE TABLE IF NOT EXISTS distances (
    id TEXT,
    dataset TEXT,
    nearest TEXT,
    identity REAL
);

CREATE TABLE IF NOT EXISTS colors (
    id TEXT PRIMARY KEY,
    color TEXT,
    source TEXT,
    dataset TEXT
);

CREATE TABLE IF NOT EXISTS seq_aliases (
    canonical_id TEXT PRIMARY KEY,
    original_header TEXT
);
"""

