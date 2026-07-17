SCHEMA = """
CREATE TABLE IF NOT EXISTS sequences (
    id TEXT PRIMARY KEY,
    sequence TEXT,
    length INTEGER,
    dataset TEXT,
    baseline_tier TEXT
);

CREATE TABLE IF NOT EXISTS taxonomy (
    id TEXT,
    taxonomy TEXT,
    confidence REAL,
    dataset TEXT
);

CREATE TABLE IF NOT EXISTS taxonomy_alt (
    id TEXT,
    ref_db TEXT,
    taxonomy TEXT,
    confidence REAL,
    best_hit TEXT,
    identity REAL
);

CREATE TABLE IF NOT EXISTS distances (
    id TEXT,
    dataset TEXT,
    nearest TEXT,
    identity REAL
);

CREATE TABLE IF NOT EXISTS classification_evidence (
    sequence_id TEXT NOT NULL,
    ref_db TEXT NOT NULL,
    best_hit TEXT,
    identity REAL,
    query_coverage REAL,
    target_coverage REAL,
    alignment_length INTEGER,
    query_length INTEGER,
    target_length INTEGER,
    mismatches INTEGER,
    gaps INTEGER,
    updated_at TEXT DEFAULT CURRENT_TIMESTAMP,
    PRIMARY KEY (sequence_id, ref_db)
);

CREATE TABLE IF NOT EXISTS colours (
    id TEXT PRIMARY KEY,
    colour TEXT,
    source TEXT,
    dataset TEXT
);

CREATE TABLE IF NOT EXISTS seq_aliases (
    canonical_id TEXT PRIMARY KEY,
    original_header TEXT
);

CREATE TABLE IF NOT EXISTS sequencing_metadata (
    id TEXT PRIMARY KEY,
    partner_id TEXT,
    dataset TEXT,
    selected_for_sequencing INTEGER DEFAULT 0,
    selected_for_wgs INTEGER DEFAULT 0,
    source_id TEXT,
    source_file TEXT,
    raw_selected_value TEXT,
    raw_commitment_value TEXT,
    operational_status TEXT DEFAULT 'RECEIVED',
    status_detail TEXT,
    manual_review_status TEXT DEFAULT 'NOT_REVIEWED'
);

CREATE TABLE IF NOT EXISTS dataset_roles (
    dataset TEXT PRIMARY KEY,
    role TEXT NOT NULL,
    genomes_available INTEGER DEFAULT 0,
    baseline_tier TEXT
);

CREATE TABLE IF NOT EXISTS assessment_snapshots (
    snapshot_id TEXT NOT NULL,
    sequence_id TEXT NOT NULL,
    dataset TEXT,
    source_path TEXT,
    assessment_json TEXT NOT NULL,
    created_at TEXT DEFAULT CURRENT_TIMESTAMP,
    PRIMARY KEY (snapshot_id, sequence_id)
);

CREATE TABLE IF NOT EXISTS selection_rounds (
    round_id TEXT PRIMARY KEY,
    mode TEXT NOT NULL,
    parameters_json TEXT NOT NULL,
    source_snapshot_count INTEGER DEFAULT 0,
    created_at TEXT DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS selection_round_members (
    round_id TEXT NOT NULL,
    sequence_id TEXT NOT NULL,
    role TEXT NOT NULL,
    round_rank INTEGER,
    recommendation_json TEXT NOT NULL,
    PRIMARY KEY (round_id, sequence_id)
);

CREATE TABLE IF NOT EXISTS schema_migrations (
    version INTEGER PRIMARY KEY,
    applied_at TEXT DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS project_runs (
    run_id TEXT PRIMARY KEY,
    workflow TEXT NOT NULL,
    dataset TEXT,
    status TEXT NOT NULL,
    manifest_path TEXT,
    started_at TEXT,
    completed_at TEXT,
    error TEXT
);

CREATE TABLE IF NOT EXISTS sequence_provenance (
    sequence_id TEXT PRIMARY KEY,
    source_manifest TEXT,
    source_sequence_file TEXT,
    source_sequence_sha256 TEXT,
    source_read_ids TEXT,
    source_read_files TEXT,
    marker_qc_class TEXT,
    marker_qc_recommendation TEXT,
    marker_qc_reasons TEXT,
    manual_review_status TEXT DEFAULT 'NOT_REVIEWED',
    chimera_call TEXT DEFAULT 'NOT_RUN',
    chimera_score REAL,
    updated_at TEXT DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS isolate_status (
    sequence_id TEXT PRIMARY KEY,
    status TEXT NOT NULL DEFAULT 'RECEIVED',
    status_detail TEXT,
    source_file TEXT,
    updated_at TEXT DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS isolate_status_events (
    event_id INTEGER PRIMARY KEY AUTOINCREMENT,
    sequence_id TEXT NOT NULL,
    old_status TEXT,
    new_status TEXT NOT NULL,
    detail TEXT,
    source_file TEXT,
    created_at TEXT DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS genome_records (
    genome_id TEXT PRIMARY KEY,
    sequence_id TEXT NOT NULL,
    accession TEXT,
    genome_status TEXT NOT NULL,
    genome_qc_pass INTEGER DEFAULT 0,
    completeness REAL,
    contamination REAL,
    gtdb_taxonomy TEXT,
    ani_cluster TEXT,
    genome_path TEXT,
    genome_sha256 TEXT,
    notes TEXT,
    source_file TEXT,
    updated_at TEXT DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS sequence_removals (
    removal_id INTEGER PRIMARY KEY AUTOINCREMENT,
    sequence_id TEXT NOT NULL,
    original_dataset TEXT,
    partner_id TEXT,
    sequence_length INTEGER,
    sequence_sha256 TEXT NOT NULL,
    taxonomy TEXT,
    reason TEXT NOT NULL,
    source_request TEXT,
    removed_records_json TEXT NOT NULL,
    removed_at TEXT DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX IF NOT EXISTS idx_genome_records_sequence ON genome_records(sequence_id);
CREATE INDEX IF NOT EXISTS idx_isolate_status_events_sequence ON isolate_status_events(sequence_id);
CREATE INDEX IF NOT EXISTS idx_sequence_removals_sequence ON sequence_removals(sequence_id);
"""
