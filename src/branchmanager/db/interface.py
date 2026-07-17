import sqlite3
import json
from branchmanager.db.schema import SCHEMA
from pathlib import Path
from branchmanager.utils.fasta import write_fasta
import hashlib
import logging
from branchmanager.taxonomy import canonicalise_sequence_id, parse_taxon_string


class Database:
    logger = logging.getLogger(__name__)

    def __init__(self, path):
        self.path = path

    def connect(self):
        if self.path != ':memory:':
            db_path = Path(self.path).expanduser()
            parent = db_path.parent
            if parent and str(parent) not in ('', '.'):
                parent.mkdir(parents=True, exist_ok=True)
            conn = sqlite3.connect(str(db_path), timeout=60.0)
        else:
            conn = sqlite3.connect(self.path, timeout=60.0)
        conn.execute('PRAGMA foreign_keys = ON')
        conn.execute('PRAGMA busy_timeout = 60000')
        return conn

    def initialise(self):
        with self.connect() as conn:
            conn.executescript(SCHEMA)
        # ensure older DBs are migrated to include dataset columns
        self._ensure_schema_up_to_date()

    def record_project_run(self, run_id, workflow, status, *, dataset=None,
                           manifest_path=None, started_at=None, completed_at=None, error=''):
        with self.connect() as conn:
            conn.execute(
                "INSERT INTO project_runs "
                "(run_id, workflow, dataset, status, manifest_path, started_at, completed_at, error) "
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?) "
                "ON CONFLICT(run_id) DO UPDATE SET workflow=excluded.workflow, "
                "dataset=excluded.dataset, status=excluded.status, manifest_path=excluded.manifest_path, "
                "started_at=excluded.started_at, completed_at=excluded.completed_at, error=excluded.error",
                (run_id, workflow, dataset, status, manifest_path, started_at, completed_at, str(error or '')),
            )

    def upsert_classification_evidence(self, rows):
        prepared = []
        for row in rows or []:
            sequence_id = str(row.get('sequence_id') or row.get('id') or '').strip()
            ref_db = str(row.get('ref_db') or '').strip()
            if not sequence_id or not ref_db:
                continue
            prepared.append((
                sequence_id, ref_db, row.get('best_hit'), row.get('identity'),
                row.get('query_coverage'), row.get('target_coverage'), row.get('alignment_length'),
                row.get('query_length'), row.get('target_length'), row.get('mismatches'), row.get('gaps'),
            ))
        if not prepared:
            return 0
        with self.connect() as conn:
            conn.executemany(
                'INSERT INTO classification_evidence '
                '(sequence_id, ref_db, best_hit, identity, query_coverage, target_coverage, '
                'alignment_length, query_length, target_length, mismatches, gaps) '
                'VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?) ON CONFLICT(sequence_id, ref_db) '
                'DO UPDATE SET best_hit=excluded.best_hit, identity=excluded.identity, '
                'query_coverage=excluded.query_coverage, target_coverage=excluded.target_coverage, '
                'alignment_length=excluded.alignment_length, query_length=excluded.query_length, '
                'target_length=excluded.target_length, mismatches=excluded.mismatches, '
                'gaps=excluded.gaps, updated_at=CURRENT_TIMESTAMP',
                prepared,
            )
        return len(prepared)

    def get_classification_evidence(self, sequence_ids=None, ref_db=None):
        clauses, params = [], []
        if sequence_ids:
            ids = [str(value) for value in sequence_ids]
            clauses.append('sequence_id IN (' + ','.join('?' for _ in ids) + ')')
            params.extend(ids)
        if ref_db:
            clauses.append('lower(ref_db) = lower(?)')
            params.append(str(ref_db))
        query = (
            'SELECT sequence_id, ref_db, best_hit, identity, query_coverage, target_coverage, '
            'alignment_length, query_length, target_length, mismatches, gaps '
            'FROM classification_evidence' + (' WHERE ' + ' AND '.join(clauses) if clauses else '')
        )
        with self.connect() as conn:
            rows = conn.execute(query, tuple(params)).fetchall()
        keys = ('sequence_id', 'ref_db', 'best_hit', 'identity', 'query_coverage', 'target_coverage',
                'alignment_length', 'query_length', 'target_length', 'mismatches', 'gaps')
        return {(str(row[0]), str(row[1])): dict(zip(keys, row)) for row in rows}

    def upsert_sequence_provenance(self, rows):
        prepared = []
        for row in rows or []:
            sequence_id = str(row.get('sequence_id') or row.get('id') or '').strip()
            if not sequence_id:
                continue
            prepared.append((
                sequence_id,
                str(row.get('source_manifest') or ''),
                str(row.get('source_sequence_file') or ''),
                str(row.get('source_sequence_sha256') or ''),
                str(row.get('source_read_ids') or ''),
                str(row.get('source_read_files') or ''),
                str(row.get('marker_qc_class') or ''),
                str(row.get('marker_qc_recommendation') or ''),
                str(row.get('marker_qc_reasons') or ''),
                str(row.get('manual_review_status') or 'NOT_REVIEWED'),
                str(row.get('chimera_call') or 'NOT_RUN'),
                row.get('chimera_score'),
            ))
        if not prepared:
            return 0
        with self.connect() as conn:
            conn.executemany(
                "INSERT INTO sequence_provenance "
                "(sequence_id, source_manifest, source_sequence_file, source_sequence_sha256, "
                "source_read_ids, source_read_files, marker_qc_class, marker_qc_recommendation, "
                "marker_qc_reasons, manual_review_status, chimera_call, chimera_score) "
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?) "
                "ON CONFLICT(sequence_id) DO UPDATE SET source_manifest=excluded.source_manifest, "
                "source_sequence_file=excluded.source_sequence_file, "
                "source_sequence_sha256=excluded.source_sequence_sha256, "
                "source_read_ids=excluded.source_read_ids, source_read_files=excluded.source_read_files, "
                "marker_qc_class=excluded.marker_qc_class, "
                "marker_qc_recommendation=excluded.marker_qc_recommendation, "
                "marker_qc_reasons=excluded.marker_qc_reasons, "
                "manual_review_status=excluded.manual_review_status, "
                "chimera_call=excluded.chimera_call, chimera_score=excluded.chimera_score, "
                "updated_at=CURRENT_TIMESTAMP",
                prepared,
            )
        return len(prepared)

    def get_sequence_provenance(self, sequence_ids=None):
        query = (
            "SELECT sequence_id, source_manifest, source_sequence_file, source_sequence_sha256, "
            "source_read_ids, source_read_files, marker_qc_class, marker_qc_recommendation, "
            "marker_qc_reasons, manual_review_status, chimera_call, chimera_score, "
            "updated_at FROM sequence_provenance"
        )
        params = ()
        if sequence_ids:
            ids = [str(value) for value in sequence_ids]
            query += ' WHERE sequence_id IN (' + ','.join('?' for _ in ids) + ')'
            params = tuple(ids)
        with self.connect() as conn:
            rows = conn.execute(query, params).fetchall()
        keys = (
            'sequence_id', 'source_manifest', 'source_sequence_file', 'source_sequence_sha256',
            'source_read_ids', 'source_read_files', 'marker_qc_class', 'marker_qc_recommendation',
            'marker_qc_reasons', 'manual_review_status', 'chimera_call', 'chimera_score', 'updated_at',
        )
        return {str(row[0]): dict(zip(keys, row)) for row in rows}

    def update_isolate_status(self, sequence_id, status, *, detail='', source_file=''):
        sequence_id = str(sequence_id).strip()
        status = str(status).strip().upper()
        if not sequence_id or not status:
            raise ValueError('sequence_id and status are required')
        with self.connect() as conn:
            previous = conn.execute(
                'SELECT status FROM isolate_status WHERE sequence_id = ?', (sequence_id,),
            ).fetchone()
            old_status = str(previous[0]) if previous else ''
            conn.execute(
                "INSERT INTO isolate_status (sequence_id, status, status_detail, source_file) "
                "VALUES (?, ?, ?, ?) ON CONFLICT(sequence_id) DO UPDATE SET "
                "status=excluded.status, status_detail=excluded.status_detail, "
                "source_file=excluded.source_file, updated_at=CURRENT_TIMESTAMP",
                (sequence_id, status, str(detail or ''), str(source_file or '')),
            )
            if old_status != status or detail:
                conn.execute(
                    "INSERT INTO isolate_status_events "
                    "(sequence_id, old_status, new_status, detail, source_file) VALUES (?, ?, ?, ?, ?)",
                    (sequence_id, old_status, status, str(detail or ''), str(source_file or '')),
                )
            conn.execute(
                "UPDATE sequencing_metadata SET operational_status = ?, status_detail = ? WHERE id = ?",
                (status, str(detail or ''), sequence_id),
            )

    def get_isolate_statuses(self, sequence_ids=None):
        query = 'SELECT sequence_id, status, status_detail, source_file, updated_at FROM isolate_status'
        params = ()
        if sequence_ids:
            ids = [str(value) for value in sequence_ids]
            query += ' WHERE sequence_id IN (' + ','.join('?' for _ in ids) + ')'
            params = tuple(ids)
        with self.connect() as conn:
            rows = conn.execute(query, params).fetchall()
        return {
            str(sequence_id): {
                'sequence_id': str(sequence_id), 'status': str(status),
                'status_detail': detail or '', 'source_file': source_file or '',
                'updated_at': updated_at or '',
            }
            for sequence_id, status, detail, source_file, updated_at in rows
        }

    def upsert_genome_records(self, rows):
        prepared = []
        for row in rows or []:
            genome_id = str(row.get('genome_id') or row.get('accession') or '').strip()
            sequence_id = str(row.get('sequence_id') or row.get('isolate_id') or '').strip()
            if not genome_id or not sequence_id:
                continue
            prepared.append((
                genome_id, sequence_id, str(row.get('accession') or ''),
                str(row.get('genome_status') or 'SEQUENCED').upper(),
                int(bool(row.get('genome_qc_pass'))), row.get('completeness'), row.get('contamination'),
                str(row.get('gtdb_taxonomy') or ''), str(row.get('ani_cluster') or ''),
                str(row.get('genome_path') or ''), str(row.get('genome_sha256') or ''),
                str(row.get('notes') or ''), str(row.get('source_file') or ''),
            ))
        if not prepared:
            return 0
        with self.connect() as conn:
            conn.executemany(
                "INSERT INTO genome_records "
                "(genome_id, sequence_id, accession, genome_status, genome_qc_pass, completeness, "
                "contamination, gtdb_taxonomy, ani_cluster, genome_path, genome_sha256, notes, source_file) "
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?) "
                "ON CONFLICT(genome_id) DO UPDATE SET sequence_id=excluded.sequence_id, "
                "accession=excluded.accession, genome_status=excluded.genome_status, "
                "genome_qc_pass=excluded.genome_qc_pass, completeness=excluded.completeness, "
                "contamination=excluded.contamination, gtdb_taxonomy=excluded.gtdb_taxonomy, "
                "ani_cluster=excluded.ani_cluster, genome_path=excluded.genome_path, "
                "genome_sha256=excluded.genome_sha256, notes=excluded.notes, "
                "source_file=excluded.source_file, updated_at=CURRENT_TIMESTAMP",
                prepared,
            )
        return len(prepared)

    def get_genome_records(self):
        with self.connect() as conn:
            rows = conn.execute(
                "SELECT genome_id, sequence_id, accession, genome_status, genome_qc_pass, "
                "completeness, contamination, gtdb_taxonomy, ani_cluster, genome_path, "
                "genome_sha256, notes, source_file, updated_at FROM genome_records "
                "ORDER BY sequence_id, genome_id"
            ).fetchall()
        keys = (
            'genome_id', 'sequence_id', 'accession', 'genome_status', 'genome_qc_pass',
            'completeness', 'contamination', 'gtdb_taxonomy', 'ani_cluster', 'genome_path',
            'genome_sha256', 'notes', 'source_file', 'updated_at',
        )
        return [dict(zip(keys, row)) for row in rows]

    def plan_sequence_removals(
        self, requests, *, allow_baseline=False, allow_genome_records=False,
    ):
        """Resolve removal requests and identify protected project records."""
        planned = []
        seen = set()
        with self.connect() as conn:
            cur = conn.cursor()
            for request in requests or []:
                requested_id = str(request.get('sequence_id') or '').strip()
                reason = str(request.get('reason') or '').strip()
                source_request = str(request.get('source_request') or '').strip()
                row = {
                    'requested_id': requested_id,
                    'sequence_id': '',
                    'dataset': '',
                    'dataset_role': '',
                    'partner_id': '',
                    'sequence_length': '',
                    'sequence_sha256': '',
                    'taxonomy': '',
                    'selected_for_genome_sequencing': False,
                    'genome_already_sequenced': False,
                    'genome_records': 0,
                    'status': 'NOT_FOUND',
                    'reason': reason,
                    'source_request': source_request,
                    '_sequence': '',
                }
                if not requested_id:
                    row['status'] = 'INVALID_REQUEST'
                    planned.append(row)
                    continue
                sequence_id = self._resolve_sequence_id_with_cursor(cur, requested_id)
                if not sequence_id:
                    planned.append(row)
                    continue
                if sequence_id in seen:
                    row['sequence_id'] = sequence_id
                    row['status'] = 'DUPLICATE_REQUEST'
                    planned.append(row)
                    continue
                seen.add(sequence_id)
                sequence_row = cur.execute(
                    'SELECT s.sequence, s.length, s.dataset, '
                    'COALESCE(m.partner_id, ""), COALESCE(m.selected_for_sequencing, 0), '
                    'COALESCE(m.selected_for_wgs, 0), COALESCE(r.role, "unregistered") '
                    'FROM sequences s LEFT JOIN sequencing_metadata m ON m.id=s.id '
                    'LEFT JOIN dataset_roles r ON r.dataset=s.dataset WHERE s.id = ?',
                    (sequence_id,),
                ).fetchone()
                if not sequence_row:
                    planned.append(row)
                    continue
                sequence, length, dataset, partner_id, selected, sequenced, dataset_role = sequence_row
                taxonomy_row = cur.execute(
                    'SELECT taxonomy FROM taxonomy WHERE id = ? '
                    'ORDER BY confidence DESC, rowid DESC LIMIT 1',
                    (sequence_id,),
                ).fetchone()
                genome_count = cur.execute(
                    'SELECT COUNT(*) FROM genome_records WHERE sequence_id = ?',
                    (sequence_id,),
                ).fetchone()[0]
                sequence_text = str(sequence or '')
                row.update({
                    'sequence_id': sequence_id,
                    'dataset': dataset or '',
                    'dataset_role': dataset_role or 'unregistered',
                    'partner_id': partner_id or '',
                    'sequence_length': int(length if length is not None else len(sequence_text)),
                    'sequence_sha256': hashlib.sha256(sequence_text.encode('utf-8')).hexdigest(),
                    'taxonomy': taxonomy_row[0] if taxonomy_row else '',
                    'selected_for_genome_sequencing': bool(selected),
                    'genome_already_sequenced': bool(sequenced or genome_count),
                    'genome_records': int(genome_count),
                    'status': 'READY',
                    '_sequence': sequence_text,
                })
                if dataset_role == 'baseline' and not allow_baseline:
                    row['status'] = 'BLOCKED_BASELINE'
                elif (sequenced or genome_count) and not allow_genome_records:
                    row['status'] = 'BLOCKED_GENOME_RECORD'
                planned.append(row)
        return planned

    def apply_sequence_removals(self, planned_rows, *, source_request=''):
        """Remove planned sequences atomically while retaining an audit tombstone."""
        rows = list(planned_rows or [])
        invalid = [row for row in rows if row.get('status') != 'READY']
        if invalid:
            statuses = ', '.join(
                f"{row.get('requested_id') or '<blank>'}:{row.get('status')}" for row in invalid
            )
            raise ValueError(f'Cannot apply an incomplete Exit Interview plan: {statuses}')
        deletion_specs = (
            ('classification_evidence', 'sequence_id = ?'),
            ('taxonomy_alt', 'id = ?'),
            ('taxonomy', 'id = ?'),
            ('colours', 'id = ?'),
            ('sequencing_metadata', 'id = ?'),
            ('assessment_snapshots', 'sequence_id = ?'),
            ('selection_round_members', 'sequence_id = ?'),
            ('sequence_provenance', 'sequence_id = ?'),
            ('isolate_status_events', 'sequence_id = ?'),
            ('isolate_status', 'sequence_id = ?'),
            ('genome_records', 'sequence_id = ?'),
            ('seq_aliases', 'canonical_id = ?'),
        )
        with self.connect() as conn:
            cur = conn.cursor()
            for row in rows:
                sequence_id = str(row['sequence_id'])
                current = cur.execute(
                    'SELECT sequence, dataset FROM sequences WHERE id = ?', (sequence_id,),
                ).fetchone()
                if not current:
                    raise ValueError(f'Sequence {sequence_id!r} disappeared after the plan was created')
                current_sha256 = hashlib.sha256(str(current[0] or '').encode('utf-8')).hexdigest()
                if current_sha256 != row['sequence_sha256']:
                    raise ValueError(f'Sequence {sequence_id!r} changed after the plan was created')

                removed_counts = {}
                cur.execute('DELETE FROM distances WHERE id = ? OR nearest = ?', (sequence_id, sequence_id))
                removed_counts['distances'] = cur.rowcount
                for table, predicate in deletion_specs:
                    cur.execute(f'DELETE FROM {table} WHERE {predicate}', (sequence_id,))
                    removed_counts[table] = cur.rowcount
                cur.execute('DELETE FROM sequences WHERE id = ?', (sequence_id,))
                removed_counts['sequences'] = cur.rowcount
                if removed_counts['sequences'] != 1:
                    raise RuntimeError(f'Exit Interview did not remove sequence {sequence_id!r}')

                cur.execute(
                    'INSERT INTO sequence_removals '
                    '(sequence_id, original_dataset, partner_id, sequence_length, sequence_sha256, '
                    'taxonomy, reason, source_request, removed_records_json) '
                    'VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)',
                    (
                        sequence_id, row.get('dataset', ''), row.get('partner_id', ''),
                        row.get('sequence_length'), row.get('sequence_sha256', ''),
                        row.get('taxonomy', ''), row.get('reason', ''),
                        row.get('source_request') or source_request or '',
                        json.dumps(removed_counts, sort_keys=True, separators=(',', ':')),
                    ),
                )
                row['status'] = 'REMOVED'
                row['removed_records'] = removed_counts
            conn.commit()
        return rows

    def _ensure_schema_up_to_date(self):
        """Apply idempotent migrations and fail if project state cannot be made current."""
        with self.connect() as conn:
            cur = conn.cursor()
            tables = {
                row[0] for row in cur.execute(
                    "SELECT name FROM sqlite_master WHERE type = 'table'"
                )
            }
            if 'colors' in tables:
                old_columns = {row[1] for row in cur.execute('PRAGMA table_info(colors)')}
                dataset_expr = 'dataset' if 'dataset' in old_columns else "''"
                cur.execute(
                    'INSERT OR REPLACE INTO colours (id, colour, source, dataset) '
                    f'SELECT id, color, source, {dataset_expr} FROM colors'
                )
                cur.execute('DROP TABLE colors')
            for table, column, definition in (
                ('sequences', 'dataset', "TEXT DEFAULT 'user'"),
                ('taxonomy', 'dataset', 'TEXT'),
                ('colours', 'dataset', 'TEXT'),
                ('distances', 'dataset', 'TEXT'),
            ):
                columns = {row[1] for row in cur.execute(f'PRAGMA table_info({table})')}
                if column not in columns:
                    cur.execute(f'ALTER TABLE {table} ADD COLUMN {column} {definition}')
            cur.execute("UPDATE distances SET dataset = 'gg2' WHERE dataset IS NULL")

            # Older prototypes had no uniqueness constraint. Retain the newest
            # row deterministically before adding the indexes used by upserts.
            cur.execute(
                'DELETE FROM taxonomy WHERE rowid NOT IN '
                '(SELECT MAX(rowid) FROM taxonomy GROUP BY id, dataset)'
            )
            cur.execute(
                'DELETE FROM distances WHERE rowid NOT IN '
                '(SELECT MAX(rowid) FROM distances GROUP BY id, dataset)'
            )
            cur.execute(
                'DELETE FROM taxonomy_alt WHERE rowid NOT IN '
                '(SELECT MAX(rowid) FROM taxonomy_alt GROUP BY id, ref_db)'
            )
            cur.execute('CREATE UNIQUE INDEX IF NOT EXISTS idx_taxonomy_id_dataset ON taxonomy (id, dataset)')
            cur.execute('CREATE UNIQUE INDEX IF NOT EXISTS idx_distances_id_dataset ON distances (id, dataset)')
            cur.execute('CREATE UNIQUE INDEX IF NOT EXISTS idx_taxonomy_alt_id_refdb ON taxonomy_alt (id, ref_db)')

            # sequencing_metadata: operational state and explicit marker-QC review.
            cur.execute("PRAGMA table_info(sequencing_metadata)")
            sequencing_columns = {row[1] for row in cur.fetchall()}
            for column, definition in (
                ('selected_for_sequencing', 'INTEGER DEFAULT 0'),
                ('raw_commitment_value', 'TEXT'),
                ('operational_status', "TEXT DEFAULT 'RECEIVED'"),
                ('status_detail', 'TEXT'),
                ('manual_review_status', "TEXT DEFAULT 'NOT_REVIEWED'"),
            ):
                if column not in sequencing_columns:
                    cur.execute(f'ALTER TABLE sequencing_metadata ADD COLUMN {column} {definition}')
            conn.commit()

            cur.execute("PRAGMA table_info(sequence_provenance)")
            provenance_columns = {row[1] for row in cur.fetchall()}
            for column, definition in (
                ('chimera_call', "TEXT DEFAULT 'NOT_RUN'"),
                ('chimera_score', 'REAL'),
            ):
                if column not in provenance_columns:
                    cur.execute(f'ALTER TABLE sequence_provenance ADD COLUMN {column} {definition}')
            conn.commit()

            # Record the first formal schema generation. CREATE TABLE statements
            # above remain idempotent for databases produced by earlier prototypes.
            cur.execute(
                'INSERT OR IGNORE INTO schema_migrations (version) VALUES (?)',
                (1,),
            )
            cur.execute(
                'INSERT OR IGNORE INTO schema_migrations (version) VALUES (?)',
                (2,),
            )
            cur.execute(
                'INSERT OR IGNORE INTO schema_migrations (version) VALUES (?)',
                (3,),
            )
            cur.execute(
                'INSERT OR IGNORE INTO schema_migrations (version) VALUES (?)',
                (4,),
            )
            conn.commit()


    def get_input_ids(self, fasta):
        from branchmanager.utils.fasta import read_fasta
        return [h for h, _ in read_fasta(fasta)]

    def register_filing_cabinet(self, fasta_path, taxa_tsv=None, colour_csv=None, source='filing_cabinet', dataset='Baseline', outdir=None, shorten_ids=False):
        """Register Filing Cabinet sequences and optional taxonomy/colours into the database.

        - fasta_path: FASTA file of reference sequences to add
        - taxa_tsv: optional TSV/CSV mapping FeatureID -> Taxon -> Confidence
        - colour_csv: optional CSV with id,colour to set explicit colours
        - source: provenance label for inserted colours
        """
        from branchmanager.utils.fasta import read_fasta
        self.logger.info("[DB][FILING CABINET] Reading fasta %s", fasta_path)
        records = [(h, s) for h, s in read_fasta(fasta_path)]
        self.logger.info("[DB][FILING CABINET] Read %d records from %s", len(records), fasta_path)
        alias_entries = []
        # prepare mapping containers so downstream mapping code can reference them
        used_ids = set(self.get_all_ids())
        orig_to_short = {}
        records_to_insert = []
        if records:
            # assign runtime IDs. By default these are the supplied headers;
            # --shorten-ids requests deterministic compact IDs.

            # When not shortening IDs, only check for collisions *within* the current
            # input batch — collisions against already-stored DB IDs are fine because
            # INSERT OR IGNORE will silently skip them.
            batch_ids: set = set()  # tracks IDs assigned during this Filing Cabinet operation only
            for orig_h, seq in records:
                if shorten_ids:
                    effective_id = self.choose_effective_sequence_id(orig_h, used_ids, shorten_ids=True)
                else:
                    provided_id = self._provided_id_from_header(orig_h)
                    if not provided_id:
                        raise ValueError(f"Could not derive an ID from header: {orig_h!r}")
                    if provided_id in batch_ids:
                        raise ValueError(
                            f"Duplicate ID in input while ID shortening is disabled: {provided_id!r}. "
                            "Ensure all input FASTA headers are unique."
                        )
                    batch_ids.add(provided_id)
                    used_ids.add(provided_id)
                    effective_id = provided_id
                records_to_insert.append((effective_id, seq))
                orig_to_short[orig_h] = effective_id

            # bulk-insert sequences and aliases in a single transaction to improve performance
            alias_entries = [(short, orig) for orig, short in orig_to_short.items()]
            try:
                with self.connect() as conn:
                    cur = conn.cursor()
                    self.logger.info("[DB][FILING CABINET] Inserting %d sequences into DB (dataset=%s)", len(records_to_insert), dataset)
                    seq_rows = []
                    for sid, seq in records_to_insert:
                        rid = self._provided_id_from_header(sid)
                        if rid is None:
                            continue
                        seq_rows.append((rid, seq, len(seq), dataset))
                    cur.executemany("INSERT OR IGNORE INTO sequences (id, sequence, length, dataset) VALUES (?, ?, ?, ?)", seq_rows)
                    if alias_entries:
                        cur.executemany("INSERT OR REPLACE INTO seq_aliases (canonical_id, original_header) VALUES (?, ?)", alias_entries)
                    conn.commit()
                    self.logger.info("[DB][FILING CABINET] Inserted %d sequences and %d alias mappings", len(seq_rows), len(alias_entries))
            except Exception as e:
                self.logger.warning("[DB][FILING CABINET] Bulk insert failed, falling back to per-record inserts: %s", e)
                self.insert_sequences(records_to_insert, dataset=dataset)
                if alias_entries:
                    self.insert_aliases(alias_entries)

        # load taxonomy TSV/CSV if provided
        if taxa_tsv:
            self.logger.info("[DB][FILING CABINET] Loading taxa table %s", taxa_tsv)
            tax_entries = []
            from branchmanager.taxonomy_io import iter_taxonomy_assignment_rows
            for row in iter_taxonomy_assignment_rows(taxa_tsv):
                tax_entries.append((row.get('id'), row.get('taxonomy'), row.get('confidence')))
            if tax_entries:
                # We only want taxonomy rows for sequences present in the provided FASTA
                # Build a canonical->short map for the records we just loaded so we only
                # insert taxonomy for those registered baseline sequences.
                canonical_to_short = {}
                for orig, short in orig_to_short.items():
                    try:
                        c = self._canonical_from_header(orig)
                    except Exception:
                        c = orig
                    canonical_to_short[c] = short

                mapped = []
                skipped = 0
                for fid, tax, conf in tax_entries:
                    sid = None
                    # direct match to original header (exact)
                    if fid in orig_to_short:
                        sid = orig_to_short[fid]
                    else:
                        # try exact canonical match against the registered baseline records
                        try:
                            cid = self._canonical_from_header(fid)
                        except Exception:
                            cid = fid
                        if cid in canonical_to_short:
                            sid = canonical_to_short[cid]

                    if sid is None:
                        # do NOT insert taxonomy rows for IDs not present in this Filing Cabinet operation
                        skipped += 1
                        continue

                    mapped.append((sid, tax, conf, dataset))

                self.logger.info("[DB][FILING CABINET] Taxa table contained %d rows; %d mapped to this dataset, %d skipped", len(tax_entries), len(mapped), skipped)

                if mapped:
                    try:
                        with self.connect() as conn:
                            cur = conn.cursor()
                            self.logger.info("[DB][FILING CABINET] Inserting %d taxonomy entries for dataset %s", len(mapped), dataset)
                            cur.executemany("INSERT OR REPLACE INTO taxonomy (id, taxonomy, confidence, dataset) VALUES (?, ?, ?, ?)", mapped)
                            conn.commit()
                            self.logger.info("[DB][FILING CABINET] Inserted taxonomy entries for dataset %s", dataset)
                    except Exception as e:
                        self.logger.warning("[DB][FILING CABINET] Bulk taxonomy insert failed: %s; falling back", e)
                        self.insert_taxonomy(mapped)

        # Generate deterministic colours when no explicit mapping is supplied.
        colour_entries = []
        if colour_csv:
            import csv
            with open(colour_csv) as cc:
                reader = csv.DictReader(cc)
                for row in reader:
                    rid = row.get('id') or row.get('ID')
                    col = row.get('colour') or row.get('Colour')
                    if rid and col:
                        # map rid to short id if possible
                        mapped_rid = orig_to_short.get(rid, None)
                        if mapped_rid:
                            colour_entries.append((mapped_rid, col, source, dataset))
                        else:
                            colour_entries.append((rid, col, source, dataset))

        if not colour_entries:
            try:
                from branchmanager.pipeline.itol import _name_to_colour
                with self.connect() as conn:
                    cur = conn.cursor()
                    cur.execute('SELECT id, taxonomy FROM taxonomy WHERE dataset = ?', (dataset,))
                    rows = cur.fetchall()
                for rid, tax in rows:
                    parsed = parse_taxon_string(tax or '')
                    ge = parsed.get('g') or parsed.get('f') or rid
                    # if we have a short id mapping for this rid, use dataset mapping
                    short = None
                    # rid here is expected to be canonical id stored in taxonomy; try to find short id
                    for orig, shortid in orig_to_short.items():
                        if self._canonical_from_header(orig) == rid or orig == rid:
                            short = shortid
                            break
                    if short:
                        colour_entries.append((short, _name_to_colour(ge), source, dataset))
                    else:
                        colour_entries.append((rid, _name_to_colour(ge), source, dataset))
            except Exception:
                # Preserve deterministic output even when taxonomy is unavailable.
                for rid, _ in records:
                    # map to short id if available
                    short = orig_to_short.get(rid, rid)
                    colour_entries.append((short, f"#{hash(rid) & 0xFFFFFF:06x}", source, dataset))

        if colour_entries:
            try:
                with self.connect() as conn:
                    cur = conn.cursor()
                    self.logger.info("[DB][FILING CABINET] Inserting %d colour entries", len(colour_entries))
                    cur.executemany("INSERT OR REPLACE INTO colours (id, colour, source, dataset) VALUES (?, ?, ?, ?)", colour_entries)
                    conn.commit()
                    self.logger.info("[DB][FILING CABINET] Inserted colours for dataset %s", dataset)
            except Exception as e:
                self.logger.warning("[DB][FILING CABINET] Bulk colour insert failed: %s; falling back", e)
                self.insert_colours(colour_entries)

        # write mapped FASTA with runtime IDs to outdir if requested so callers can avoid exporting from DB
        mapped_fasta_path = None
        if outdir and records_to_insert:
            try:
                p = Path(outdir)
                p.mkdir(parents=True, exist_ok=True)
                mapped_fasta_path = str(p / f"filing_cabinet_{dataset}_sequences.fasta")
                write_fasta([(short, seq) for short, seq in records_to_insert], mapped_fasta_path)
                self.logger.info("[DB][FILING CABINET] Wrote mapped fasta with runtime IDs to %s", mapped_fasta_path)
            except Exception as e:
                self.logger.warning("[DB][FILING CABINET] Failed to write mapped fasta to outdir: %s", e)

        # return alias entries and optional mapped fasta path
        self.logger.info("[DB][FILING CABINET] Baseline registration complete; %d alias entries created", len(alias_entries) if alias_entries else 0)
        return alias_entries, mapped_fasta_path

    def get_sequences_fasta(self, outpath, dataset=None):
        """Write sequences stored in the database to outpath FASTA.

        If `dataset` is provided, only sequences belonging to that dataset
        are exported. Returns True if sequences were written, False if none.
        """
        def _write_rows(rows):
            records = [(r[0], r[1]) for r in rows]
            write_fasta(records, outpath)

        p = Path(outpath)
        p.parent.mkdir(parents=True, exist_ok=True)

        with self.connect() as conn:
            cur = conn.cursor()
            if dataset:
                cur.execute("SELECT id, sequence FROM sequences WHERE dataset = ?", (dataset,))
            else:
                cur.execute("SELECT id, sequence FROM sequences")
            rows = cur.fetchall()

        if not rows:
            return False

        _write_rows(rows)
        return True

    def insert_aliases(self, alias_entries):
        """Insert alias mappings: iterable of (canonical_id, original_header)"""
        with self.connect() as conn:
            cur = conn.cursor()
            for cid, orig in alias_entries:
                cur.execute("INSERT OR REPLACE INTO seq_aliases (canonical_id, original_header) VALUES (?, ?)", (cid, orig))
            conn.commit()

    def get_datasets(self):
        """Return a list of dataset names present in the sequences table, excluding 'user' and NULL."""
        with self.connect() as conn:
            cur = conn.cursor()
            cur.execute("SELECT DISTINCT dataset FROM sequences WHERE dataset IS NOT NULL AND dataset != 'user' AND dataset != ''")
            rows = cur.fetchall()
        return [r[0] for r in rows]

    def upsert_dataset_role(self, dataset, role, genomes_available=False):
        """Record whether a dataset is a baseline or partner-candidate collection."""
        dataset = str(dataset or '').strip()
        role = str(role or '').strip().lower()
        if not dataset or role not in ('baseline', 'candidate'):
            raise ValueError('dataset role must be baseline or candidate')
        with self.connect() as conn:
            conn.execute(
                "INSERT OR REPLACE INTO dataset_roles (dataset, role, genomes_available) VALUES (?, ?, ?)",
                (dataset, role, 1 if genomes_available else 0),
            )
            conn.commit()

    def get_dataset_roles(self):
        """Return {dataset: {role, genomes_available}} for registered datasets."""
        with self.connect() as conn:
            rows = conn.execute(
                "SELECT dataset, role, genomes_available FROM dataset_roles"
            ).fetchall()
        return {
            str(dataset): {
                'role': str(role),
                'genomes_available': bool(genomes_available),
            }
            for dataset, role, genomes_available in rows
        }

    def get_dataset_names_by_role(self, role):
        wanted = str(role or '').strip().lower()
        return sorted(
            dataset for dataset, metadata in self.get_dataset_roles().items()
            if metadata.get('role') == wanted
        )

    def export_dataset_fasta(self, dataset, outpath):
        """Write sequences belonging to a given dataset to outpath FASTA. Returns True if any sequences written."""
        p = Path(outpath)
        p.parent.mkdir(parents=True, exist_ok=True)
        with self.connect() as conn:
            cur = conn.cursor()
            cur.execute("SELECT id, sequence FROM sequences WHERE dataset = ?", (dataset,))
            rows = cur.fetchall()
        if not rows:
            return False

        records = [(r[0], r[1]) for r in rows]
        write_fasta(records, outpath)
        return True

    def get_all_ids(self):
        """Return a set of all sequence IDs stored in the DB."""
        with self.connect() as conn:
            cur = conn.cursor()
            cur.execute("SELECT id FROM sequences")
            rows = cur.fetchall()
        return set(r[0] for r in rows)

    def assert_sequence_compatible(self, sid, sequence):
        """Reject reuse of a persistent ID for a different nucleotide sequence."""
        key = self.resolve_sequence_id(sid)
        if not key:
            return None
        with self.connect() as conn:
            row = conn.execute(
                'SELECT sequence FROM sequences WHERE id = ? LIMIT 1',
                (key,),
            ).fetchone()
        if not row:
            return None

        def normalise(value):
            return ''.join(str(value or '').split()).upper()

        if normalise(row[0]) != normalise(sequence):
            raise ValueError(
                f'Sequence ID {key!r} already exists in the project database with a '
                'different nucleotide sequence. Rolling sequence IDs are immutable; '
                'correct the input or assign the revised sequence a new unique ID.'
            )
        return key

    def generate_short_id(self, header: str, used_set: set):
        """Generate a deterministic short ID (3 letters + 2 digits) for a header.

        Uses md5(header|tries) to derive values, avoids collisions by consulting
        used_set (which is mutated to include the returned id).
        """
        i = 0
        while True:
            h = hashlib.md5((header + '|' + str(i)).encode('utf-8')).hexdigest()
            val = int(h[:8], 16)
            letters_val = val % (26 ** 3)
            dv = (val >> 12) % 100
            letters = ''
            vv = letters_val
            for _ in range(3):
                letters = chr(ord('A') + (vv % 26)) + letters
                vv //= 26
            sid = f"{letters}{dv:02d}"
            if sid not in used_set:
                used_set.add(sid)
                return sid
            i += 1

    def choose_effective_sequence_id(self, header: str, used_set: set, shorten_ids: bool = False):
        """Return the runtime ID to use for a sequence.

        When `shorten_ids` is True, generate a deterministic compact ID.
        Otherwise, use the source FASTA header exactly as supplied (trimmed)
        and fail fast on exact collisions so users do not silently merge
        distinct records.
        """
        if shorten_ids:
            return self.generate_short_id(header, used_set)
        provided_id = self._provided_id_from_header(header)
        if not provided_id:
            raise ValueError(f"Could not derive an ID from header: {header!r}")
        if provided_id in used_set:
            raise ValueError(
                f"ID collision while ID shortening is disabled: {provided_id!r}. "
                "Use unique input IDs or re-run with --shorten-ids."
            )
        used_set.add(provided_id)
        return provided_id

    def _provided_id_from_header(self, header: str) -> str:
        if header is None:
            return None
        value = str(header).strip()
        return value or None

    def _canonical_from_header(self, header: str) -> str:
        return canonicalise_sequence_id(header)

    def _candidate_ids_for_lookup(self, sid):
        candidates = []
        raw = self._provided_id_from_header(sid)
        if raw:
            candidates.append(raw)
        try:
            cid = self._canonical_from_header(sid)
        except Exception:
            cid = None
        if cid and cid not in candidates:
            candidates.append(cid)
        return candidates

    def resolve_sequence_id(self, sid):
        """Resolve a user/reference identifier to an existing sequence-table ID."""
        candidates = self._candidate_ids_for_lookup(sid)
        if not candidates:
            return None
        with self.connect() as conn:
            cur = conn.cursor()
            for cand in candidates:
                cur.execute("SELECT id FROM sequences WHERE id = ? LIMIT 1", (cand,))
                row = cur.fetchone()
                if row:
                    return row[0]
            for cand in candidates:
                cur.execute(
                    "SELECT canonical_id FROM seq_aliases WHERE original_header = ? OR canonical_id = ? LIMIT 1",
                    (cand, cand),
                )
                row = cur.fetchone()
                if row:
                    return row[0]
        return None

    def _resolve_sequence_id_with_cursor(self, cur, sid):
        candidates = self._candidate_ids_for_lookup(sid)
        if not candidates:
            return None
        for cand in candidates:
            cur.execute("SELECT id FROM sequences WHERE id = ? LIMIT 1", (cand,))
            row = cur.fetchone()
            if row:
                return row[0]
        for cand in candidates:
            cur.execute(
                "SELECT canonical_id FROM seq_aliases WHERE original_header = ? OR canonical_id = ? LIMIT 1",
                (cand, cand),
            )
            row = cur.fetchone()
            if row:
                return row[0]
        return None


    def insert_taxonomy_alt(self, alt_entries):
        """Insert alternative-database taxonomy entries.

        alt_entries: iterable of (id, ref_db, taxonomy, confidence, best_hit, identity)
        Only inserts rows for sequence IDs that exist in the sequences table.
        """
        with self.connect() as conn:
            cur = conn.cursor()
            to_insert = []
            for entry in alt_entries:
                sid, ref_db, tax, conf, best_hit, identity = entry
                key = self._resolve_sequence_id_with_cursor(cur, sid)
                if key:
                    to_insert.append((key, ref_db, tax, conf, best_hit, identity))
            if to_insert:
                cur.executemany(
                    "INSERT OR REPLACE INTO taxonomy_alt "
                    "(id, ref_db, taxonomy, confidence, best_hit, identity) VALUES (?, ?, ?, ?, ?, ?)",
                    to_insert,
                )
            conn.commit()

    def get_taxonomy_alt_for_ids(self, ids=None):
        """Return alt-db taxonomy rows for given IDs (or all if ids is None).

        Returns {id: {ref_db: (taxonomy, confidence, best_hit, identity)}}.
        """
        with self.connect() as conn:
            cur = conn.cursor()
            if ids:
                placeholders = ','.join('?' for _ in ids)
                cur.execute(
                    f"SELECT id, ref_db, taxonomy, confidence, best_hit, identity "
                    f"FROM taxonomy_alt WHERE id IN ({placeholders})",
                    tuple(ids),
                )
            else:
                cur.execute("SELECT id, ref_db, taxonomy, confidence, best_hit, identity FROM taxonomy_alt")
            rows = cur.fetchall()
        result: dict = {}
        for rid, ref_db, tax, conf, best_hit, identity in rows:
            result.setdefault(rid, {})[ref_db] = (tax, conf, best_hit, identity)
        return result

    def get_alt_ref_dbs(self):
        """Return a sorted list of distinct ref_db values stored in taxonomy_alt."""
        with self.connect() as conn:
            cur = conn.cursor()
            cur.execute("SELECT DISTINCT ref_db FROM taxonomy_alt ORDER BY ref_db")
            return [r[0] for r in cur.fetchall()]

    def insert_sequences(self, records, dataset='user'):
        """Insert sequence records into the DB.

        records: iterable of (id, seq)
        Uses INSERT OR IGNORE for sequences table (id is PK).
        """
        alias_entries = []
        with self.connect() as conn:
            cur = conn.cursor()
            for sid, seq in records:
                rid = self._provided_id_from_header(sid)
                if rid is None:
                    continue
                cur.execute("INSERT OR IGNORE INTO sequences (id, sequence, length, dataset) VALUES (?, ?, ?, ?)",
                            (rid, seq, len(seq), dataset))
                alias_entries.append((rid, sid))
            conn.commit()
        if alias_entries:
            # store alias mappings
            self.insert_aliases(alias_entries)

    def insert_taxonomy(self, tax_entries):
        """Insert taxonomy entries.

        tax_entries: iterable of (id, taxonomy, confidence)
        Uses INSERT OR REPLACE to update taxonomy info.
        """
        # Only persist taxonomy for ids that are present in the sequences table
        with self.connect() as conn:
            cur = conn.cursor()
            to_insert = []
            skipped = 0
            for entry in tax_entries:
                # accept either (sid, tax, conf) or (sid, tax, conf, dataset)
                if len(entry) == 3:
                    sid, tax, conf = entry
                    ds = ''
                else:
                    sid, tax, conf, ds = entry
                    ds = ds or ''
                key = self._resolve_sequence_id_with_cursor(cur, sid)
                # ensure this id was provided by the user (exists in sequences)
                if key:
                    to_insert.append((key, tax, conf, ds))
                else:
                    skipped += 1
            if to_insert:
                cur.executemany("INSERT OR REPLACE INTO taxonomy (id, taxonomy, confidence, dataset) VALUES (?, ?, ?, ?)", to_insert)
            conn.commit()
            if skipped:
                self.logger.info("[DB] Skipped %d taxonomy entries because ids are not present in sequences table", skipped)

    def insert_distances(self, dist_entries):
        """Insert nearest-hit distances.

        dist_entries: iterable of (id, nearest, identity)
        Uses INSERT OR REPLACE.
        """
        # Only persist distances for ids that are present in the sequences table
        with self.connect() as conn:
            cur = conn.cursor()
            to_insert = []
            skipped = 0
            for entry in dist_entries:
                # accept either (sid, nearest, identity) or (sid, dataset, nearest, identity)
                if len(entry) == 3:
                    sid, nearest, identity = entry
                    ds = ''
                else:
                    sid, ds, nearest, identity = entry
                    ds = ds or ''
                key = self._resolve_sequence_id_with_cursor(cur, sid)
                # Only insert if the id exists in sequences
                if not key:
                    skipped += 1
                    continue
                nearest_key = self._resolve_sequence_id_with_cursor(cur, nearest) if nearest is not None else None
                to_insert.append((key, ds, nearest_key or nearest, identity))
            if to_insert:
                cur.executemany("INSERT OR REPLACE INTO distances (id, dataset, nearest, identity) VALUES (?, ?, ?, ?)", to_insert)
            conn.commit()
            if skipped:
                self.logger.info("[DB] Skipped %d distance entries because ids are not present in sequences table", skipped)

    def insert_colours(self, colour_entries):
        """Insert or replace colour entries.

        ``colour_entries`` contains ``(id, colour, source[, dataset])`` tuples.
        """
        with self.connect() as conn:
            cur = conn.cursor()
            to_insert = []
            skipped = 0
            for entry in colour_entries:
                if len(entry) == 3:
                    sid, colour, source = entry
                    ds = ''
                else:
                    sid, colour, source, ds = entry
                    ds = ds or ''
                key = self._resolve_sequence_id_with_cursor(cur, sid)
                if key:
                    to_insert.append((key, colour, source, ds))
                else:
                    skipped += 1
            if to_insert:
                cur.executemany("INSERT OR REPLACE INTO colours (id, colour, source, dataset) VALUES (?, ?, ?, ?)", to_insert)
            conn.commit()
            if skipped:
                self.logger.info("[DB] Skipped %d colour entries because ids are not present in sequences table", skipped)

    def get_colour_map(self, ids=None):
        """Return the stored ``id -> colour`` mapping."""
        with self.connect() as conn:
            cur = conn.cursor()
            if ids:
                placeholders = ','.join('?' for _ in ids)
                cur.execute(f"SELECT id, colour FROM colours WHERE id IN ({placeholders})", tuple(ids))
            else:
                cur.execute("SELECT id, colour FROM colours")
            rows = cur.fetchall()
        return {r[0]: r[1] for r in rows}

    def upsert_sequencing_metadata(self, metadata_rows):
        """Insert/update rolling partner WGS-selection metadata.

        ``selected_for_sequencing`` is a project commitment. The legacy-named
        ``selected_for_wgs`` column means a genome is already available.
        """
        rows = []
        with self.connect() as conn:
            cur = conn.cursor()
            for row in metadata_rows:
                sid = row.get('id') if isinstance(row, dict) else None
                key = self._resolve_sequence_id_with_cursor(cur, sid) or self._provided_id_from_header(sid)
                if not key:
                    continue
                rows.append((
                    key,
                    row.get('partner_id') or key,
                    row.get('dataset') or '',
                    1 if bool(row.get('selected_for_sequencing')) else 0,
                    1 if bool(row.get('selected_for_wgs')) else 0,
                    row.get('source_id') or sid,
                    row.get('source_file') or '',
                    row.get('raw_selected_value') if row.get('raw_selected_value') is not None else '',
                    row.get('raw_commitment_value') if row.get('raw_commitment_value') is not None else '',
                ))
            if rows:
                cur.executemany(
                    "INSERT INTO sequencing_metadata "
                    "(id, partner_id, dataset, selected_for_sequencing, selected_for_wgs, "
                    "source_id, source_file, raw_selected_value, raw_commitment_value) "
                    "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?) ON CONFLICT(id) DO UPDATE SET "
                    "partner_id=excluded.partner_id, dataset=excluded.dataset, "
                    "selected_for_sequencing=excluded.selected_for_sequencing, "
                    "selected_for_wgs=excluded.selected_for_wgs, source_id=excluded.source_id, "
                    "source_file=excluded.source_file, raw_selected_value=excluded.raw_selected_value, "
                    "raw_commitment_value=excluded.raw_commitment_value",
                    rows,
                )
            conn.commit()
        return len(rows)

    def get_sequencing_metadata_for_ids(self, ids=None):
        with self.connect() as conn:
            cur = conn.cursor()
            if ids:
                resolved = []
                for sid in ids:
                    key = self._resolve_sequence_id_with_cursor(cur, sid)
                    if key and key not in resolved:
                        resolved.append(key)
                if not resolved:
                    return {}
                placeholders = ','.join('?' for _ in resolved)
                cur.execute(
                    f"SELECT id, partner_id, dataset, selected_for_sequencing, selected_for_wgs, "
                    f"source_id, source_file, raw_selected_value, raw_commitment_value "
                    f"FROM sequencing_metadata WHERE id IN ({placeholders})",
                    tuple(resolved),
                )
            else:
                cur.execute(
                    "SELECT id, partner_id, dataset, selected_for_sequencing, selected_for_wgs, "
                    "source_id, source_file, raw_selected_value, raw_commitment_value "
                    "FROM sequencing_metadata"
                )
            rows = cur.fetchall()
        return {
            rid: {
                'partner_id': partner_id or rid,
                'dataset': dataset or '',
                'selected_for_sequencing': bool(selected_for_sequencing),
                'selected_for_wgs': bool(selected_for_wgs),
                'genome_available': bool(selected_for_wgs),
                'source_id': source_id or rid,
                'source_file': source_file or '',
                'raw_selected_value': raw_selected_value or '',
                'raw_commitment_value': raw_commitment_value or '',
            }
            for rid, partner_id, dataset, selected_for_sequencing, selected_for_wgs,
            source_id, source_file, raw_selected_value, raw_commitment_value in rows
        }

    def save_assessment_snapshot(self, snapshot_id, assessment_rows, *, dataset='', source_path=''):
        """Persist one complete performance review without mutating genome status."""
        rows = []
        for row in assessment_rows:
            sequence_id = str(row.get('id') or row.get('ID') or row.get('SequenceID') or '').strip()
            if not sequence_id:
                continue
            rows.append((
                str(snapshot_id),
                sequence_id,
                str(dataset or ''),
                str(source_path or ''),
                json.dumps(dict(row), sort_keys=True, separators=(',', ':')),
            ))
        with self.connect() as conn:
            if rows:
                conn.executemany(
                    "INSERT OR REPLACE INTO assessment_snapshots "
                    "(snapshot_id, sequence_id, dataset, source_path, assessment_json) "
                    "VALUES (?, ?, ?, ?, ?)",
                    rows,
                )
            conn.commit()
        return len(rows)

    def get_latest_assessment_rows(self):
        """Return each sequence's latest assessment from its owning dataset.

        Explicit project-wide imports remain eligible. Candidate-run snapshots
        cannot supersede rows belonging to another partner dataset.
        """
        latest = {}
        with self.connect() as conn:
            candidate_datasets = {
                str(row[0]) for row in conn.execute(
                    "SELECT dataset FROM dataset_roles WHERE role = 'candidate'"
                ).fetchall()
            }
            sequence_datasets = {
                str(sequence_id): str(dataset or '')
                for sequence_id, dataset in conn.execute(
                    "SELECT id, dataset FROM sequences"
                ).fetchall()
            }
            rows = conn.execute(
                "SELECT snapshot_id, sequence_id, dataset, source_path, assessment_json, created_at "
                "FROM assessment_snapshots ORDER BY created_at, rowid"
            ).fetchall()
        for snapshot_id, sequence_id, dataset, source_path, payload, created_at in rows:
            snapshot_dataset = str(dataset or '')
            if (
                snapshot_dataset in candidate_datasets
                and sequence_datasets.get(str(sequence_id)) != snapshot_dataset
            ):
                continue
            try:
                assessment = json.loads(payload)
            except (TypeError, ValueError):
                continue
            assessment['_snapshot_id'] = snapshot_id
            assessment['_snapshot_dataset'] = snapshot_dataset
            assessment['_snapshot_source_path'] = source_path or ''
            assessment['_snapshot_created_at'] = created_at or ''
            latest[str(sequence_id)] = assessment
        return latest

    def save_selection_round(self, round_id, recommendations, *, mode='quarterly_review', parameters=None):
        """Persist one immutable project-wide recommendation round."""
        parameters = dict(parameters or {})
        source_snapshots = {
            str(row.get('source_snapshot') or '') for row in recommendations
            if row.get('source_snapshot')
        }
        with self.connect() as conn:
            conn.execute(
                "INSERT OR REPLACE INTO selection_rounds "
                "(round_id, mode, parameters_json, source_snapshot_count) VALUES (?, ?, ?, ?)",
                (
                    str(round_id), str(mode),
                    json.dumps(parameters, sort_keys=True, separators=(',', ':')),
                    len(source_snapshots),
                ),
            )
            conn.execute("DELETE FROM selection_round_members WHERE round_id = ?", (str(round_id),))
            member_rows = []
            for row in recommendations:
                sequence_id = str(row.get('sequence_id') or row.get('id') or '').strip()
                if not sequence_id:
                    continue
                rank = row.get('round_rank')
                try:
                    rank = int(rank) if rank not in (None, '', 'NA') else None
                except (TypeError, ValueError):
                    rank = None
                member_rows.append((
                    str(round_id), sequence_id, str(row.get('role') or 'NOT_SELECTED'), rank,
                    json.dumps(dict(row), sort_keys=True, separators=(',', ':')),
                ))
            if member_rows:
                conn.executemany(
                    "INSERT INTO selection_round_members "
                    "(round_id, sequence_id, role, round_rank, recommendation_json) "
                    "VALUES (?, ?, ?, ?, ?)",
                    member_rows,
                )
            conn.commit()
        return len(member_rows)
