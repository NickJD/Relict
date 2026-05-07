import sqlite3
from relict.db.schema import SCHEMA
from pathlib import Path
from relict.utils.fasta import write_fasta
import hashlib
import logging
from relict.taxonomy import canonicalize_sequence_id, parse_taxon_string


class Database:
    logger = logging.getLogger(__name__)

    def __init__(self, path):
        self.path = path

    def connect(self):
        return sqlite3.connect(self.path)

    def initialise(self):
        with self.connect() as conn:
            conn.executescript(SCHEMA)
        # ensure older DBs are migrated to include dataset columns
        self._ensure_schema_up_to_date()

    def _ensure_schema_up_to_date(self):
        """Perform lightweight migrations to add dataset columns to existing tables if missing.

        This function is safe to run repeatedly and will add columns or set defaults as needed.
        """
        with self.connect() as conn:
            cur = conn.cursor()
            # sequences: add dataset column if missing
            cur.execute("PRAGMA table_info(sequences)")
            cols = [r[1] for r in cur.fetchall()]
            if 'dataset' not in cols:
                try:
                    cur.execute("ALTER TABLE sequences ADD COLUMN dataset TEXT DEFAULT 'user'")
                    conn.commit()
                except Exception:
                    pass

            # taxonomy: add dataset column if missing
            cur.execute("PRAGMA table_info(taxonomy)")
            cols = [r[1] for r in cur.fetchall()]
            if 'dataset' not in cols:
                try:
                    cur.execute("ALTER TABLE taxonomy ADD COLUMN dataset TEXT")
                    conn.commit()
                except Exception:
                    pass

            # colors: add dataset column if missing
            cur.execute("PRAGMA table_info(colors)")
            cols = [r[1] for r in cur.fetchall()]
            if 'dataset' not in cols:
                try:
                    cur.execute("ALTER TABLE colors ADD COLUMN dataset TEXT")
                    conn.commit()
                except Exception:
                    pass

            # distances: add dataset column if missing and set default for existing rows
            cur.execute("PRAGMA table_info(distances)")
            cols = [r[1] for r in cur.fetchall()]
            if 'dataset' not in cols:
                try:
                    cur.execute("ALTER TABLE distances ADD COLUMN dataset TEXT")
                    cur.execute("UPDATE distances SET dataset = 'gg2' WHERE dataset IS NULL")
                    conn.commit()
                except Exception:
                    pass

            # ensure unique indexes exist so INSERT OR REPLACE behaves as intended
            try:
                cur.execute("CREATE UNIQUE INDEX IF NOT EXISTS idx_taxonomy_id_dataset ON taxonomy (id, dataset)")
                cur.execute("CREATE UNIQUE INDEX IF NOT EXISTS idx_distances_id_dataset ON distances (id, dataset)")
                conn.commit()
            except Exception:
                pass

            # taxonomy_alt: create table and unique index if missing
            cur.execute("PRAGMA table_info(taxonomy_alt)")
            if not cur.fetchall():
                try:
                    cur.execute(
                        "CREATE TABLE IF NOT EXISTS taxonomy_alt "
                        "(id TEXT, ref_db TEXT, taxonomy TEXT, confidence REAL, best_hit TEXT, identity REAL)"
                    )
                    conn.commit()
                except Exception:
                    pass
            try:
                cur.execute(
                    "CREATE UNIQUE INDEX IF NOT EXISTS idx_taxonomy_alt_id_refdb ON taxonomy_alt (id, ref_db)"
                )
                conn.commit()
            except Exception:
                pass


    def get_input_ids(self, fasta):
        from relict.utils.fasta import read_fasta
        return [h for h, _ in read_fasta(fasta)]

    def preload_from_files(self, fasta_path, taxa_tsv=None, color_csv=None, source='preload', dataset='preload', outdir=None, shorten_ids=True):
        """Preload sequences (and optionally taxonomy/colors) into the DB.

        - fasta_path: FASTA file of reference sequences to add
        - taxa_tsv: optional TSV mapping FeatureID -> Taxon -> Confidence
        - color_csv: optional CSV with id,color to set explicit colors
        - source: string to record as the source for inserted colors
        """
        from relict.utils.fasta import read_fasta
        self.logger.info("[DB][PRELOAD] Reading fasta %s", fasta_path)
        records = [(h, s) for h, s in read_fasta(fasta_path)]
        self.logger.info("[DB][PRELOAD] Read %d records from %s", len(records), fasta_path)
        alias_entries = []
        # prepare mapping containers so downstream mapping code can reference them
        used_ids = set(self.get_all_ids())
        orig_to_short = {}
        records_to_insert = []
        if records:
            # assign deterministic short IDs (3 letters + 2 digits) per original header

            # When not shortening IDs, only check for collisions *within* the current
            # input batch — collisions against already-stored DB IDs are fine because
            # INSERT OR IGNORE will silently skip them.
            batch_ids: set = set()  # tracks IDs assigned during this preload only
            for orig_h, seq in records:
                if shorten_ids:
                    effective_id = self.choose_effective_sequence_id(orig_h, used_ids, shorten_ids=True)
                else:
                    cid = self._canonical_from_header(orig_h)
                    if not cid:
                        raise ValueError(f"Could not derive a canonical ID from header: {orig_h!r}")
                    if cid in batch_ids:
                        raise ValueError(
                            f"Duplicate ID in input while --no-shorten-ids is active: {cid!r}. "
                            "Ensure all input FASTA headers are unique."
                        )
                    batch_ids.add(cid)
                    used_ids.add(cid)
                    effective_id = cid
                records_to_insert.append((effective_id, seq))
                orig_to_short[orig_h] = effective_id

            # bulk-insert sequences and aliases in a single transaction to improve performance
            alias_entries = [(short, orig) for orig, short in orig_to_short.items()]
            try:
                with self.connect() as conn:
                    cur = conn.cursor()
                    self.logger.info("[DB][PRELOAD] Inserting %d sequences into DB (dataset=%s)", len(records_to_insert), dataset)
                    seq_rows = []
                    for sid, seq in records_to_insert:
                        cid = self._canonical_from_header(sid)
                        if cid is None:
                            continue
                        seq_rows.append((cid, seq, len(seq), dataset))
                    cur.executemany("INSERT OR IGNORE INTO sequences (id, sequence, length, dataset) VALUES (?, ?, ?, ?)", seq_rows)
                    if alias_entries:
                        cur.executemany("INSERT OR REPLACE INTO seq_aliases (canonical_id, original_header) VALUES (?, ?)", alias_entries)
                    conn.commit()
                    self.logger.info("[DB][PRELOAD] Inserted %d sequences and %d alias mappings", len(seq_rows), len(alias_entries))
            except Exception as e:
                self.logger.warning("[DB][PRELOAD] Bulk insert failed, falling back to per-record inserts: %s", e)
                self.insert_sequences(records_to_insert, dataset=dataset)
                if alias_entries:
                    self.insert_aliases(alias_entries)

        # load taxonomy TSV if provided
        if taxa_tsv:
            self.logger.info("[DB][PRELOAD] Loading taxa TSV %s", taxa_tsv)
            tax_entries = []
            import gzip
            open_fn = gzip.open if str(taxa_tsv).endswith('.gz') else open
            with open_fn(taxa_tsv, 'rt') as t:
                # skip header if present
                first = t.readline()
                if 'Feature' in first or 'Taxon' in first:
                    pass
                else:
                    parts = first.strip().split('\t')
                    if len(parts) >= 2:
                        tax_entries.append((parts[0], parts[1], float(parts[2]) if len(parts) > 2 else None))
                for line in t:
                    parts = line.strip().split('\t')
                    if len(parts) < 2:
                        continue
                    conf = float(parts[2]) if len(parts) > 2 else None
                    fid = parts[0]
                    tax_entries.append((fid, parts[1], conf))
            if tax_entries:
                # We only want taxonomy rows for sequences present in the provided FASTA
                # Build a canonical->short map for the records we just loaded so we only
                # insert taxonomy for those preloaded sequences.
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
                        # try exact canonical match against the preloaded records
                        try:
                            cid = self._canonical_from_header(fid)
                        except Exception:
                            cid = fid
                        if cid in canonical_to_short:
                            sid = canonical_to_short[cid]

                    if sid is None:
                        # do NOT insert taxonomy rows for IDs not present in this preload
                        skipped += 1
                        continue

                    mapped.append((sid, tax, conf, dataset))

                self.logger.info("[DB][PRELOAD] Taxa TSV contained %d rows; %d mapped to this dataset, %d skipped", len(tax_entries), len(mapped), skipped)

                if mapped:
                    try:
                        with self.connect() as conn:
                            cur = conn.cursor()
                            self.logger.info("[DB][PRELOAD] Inserting %d taxonomy entries for dataset %s", len(mapped), dataset)
                            cur.executemany("INSERT OR REPLACE INTO taxonomy (id, taxonomy, confidence, dataset) VALUES (?, ?, ?, ?)", mapped)
                            conn.commit()
                            self.logger.info("[DB][PRELOAD] Inserted taxonomy entries for dataset %s", dataset)
                    except Exception as e:
                        self.logger.warning("[DB][PRELOAD] Bulk taxonomy insert failed: %s; falling back", e)
                        self.insert_taxonomy(mapped)

        # load colors if provided, else generate colors automatically for these records
        color_entries = []
        if color_csv:
            import csv
            with open(color_csv) as cc:
                reader = csv.DictReader(cc)
                for row in reader:
                    rid = row.get('id') or row.get('ID')
                    col = row.get('color') or row.get('colour') or row.get('Color')
                    if rid and col:
                        # map rid to short id if possible
                        mapped_rid = orig_to_short.get(rid, None)
                        if mapped_rid:
                            color_entries.append((mapped_rid, col, source, dataset))
                        else:
                            color_entries.append((rid, col, source, dataset))

        # if no color entries, derive colors from taxonomy/genus
        if not color_entries:
            try:
                from relict.pipeline.itol import _name_to_color
                with self.connect() as conn:
                    cur = conn.cursor()
                    # restrict color derivation to taxonomy rows for this dataset only
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
                        color_entries.append((short, _name_to_color(ge), source, dataset))
                    else:
                        color_entries.append((rid, _name_to_color(ge), source, dataset))
            except Exception:
                # fallback: color by id
                for rid, _ in records:
                    # map to short id if available
                    short = orig_to_short.get(rid, rid)
                    color_entries.append((short, f"#{hash(rid) & 0xFFFFFF:06x}", source, dataset))

        if color_entries:
            try:
                with self.connect() as conn:
                    cur = conn.cursor()
                    self.logger.info("[DB][PRELOAD] Inserting %d color entries", len(color_entries))
                    cur.executemany("INSERT OR REPLACE INTO colors (id, color, source, dataset) VALUES (?, ?, ?, ?)", color_entries)
                    conn.commit()
                    self.logger.info("[DB][PRELOAD] Inserted colors for dataset %s", dataset)
            except Exception as e:
                self.logger.warning("[DB][PRELOAD] Bulk color insert failed: %s; falling back", e)
                self.insert_colors(color_entries)

        # write mapped fasta (short IDs) to outdir if requested so callers can avoid exporting from DB
        mapped_fasta_path = None
        if outdir and records_to_insert:
            try:
                p = Path(outdir)
                p.mkdir(parents=True, exist_ok=True)
                mapped_fasta_path = str(p / f"preload_{dataset}_seqs.fasta")
                write_fasta([(short, seq) for short, seq in records_to_insert], mapped_fasta_path)
                self.logger.info("[DB][PRELOAD] Wrote mapped fasta with short IDs to %s", mapped_fasta_path)
            except Exception as e:
                self.logger.warning("[DB][PRELOAD] Failed to write mapped fasta to outdir: %s", e)

        # return alias entries and optional mapped fasta path
        self.logger.info("[DB][PRELOAD] Preload complete; %d alias entries created", len(alias_entries) if alias_entries else 0)
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

    def choose_effective_sequence_id(self, header: str, used_set: set, shorten_ids: bool = True):
        """Return the runtime ID to use for a sequence.

        When `shorten_ids` is True, generate a deterministic compact ID.
        Otherwise, use the canonicalized source ID directly and fail fast on
        collisions so users do not silently merge distinct records.
        """
        if shorten_ids:
            return self.generate_short_id(header, used_set)
        cid = self._canonical_from_header(header)
        if not cid:
            raise ValueError(f"Could not derive a canonical ID from header: {header!r}")
        if cid in used_set:
            raise ValueError(
                f"Canonical ID collision while --no-shorten-ids is active: {cid!r}. "
                "Use unique input IDs or re-run with ID shortening enabled."
            )
        used_set.add(cid)
        return cid

    def _canonical_from_header(self, header: str) -> str:
        return canonicalize_sequence_id(header)


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
                cid = self._canonical_from_header(sid) if sid else None
                key = cid or sid
                cur.execute("SELECT 1 FROM sequences WHERE id = ? LIMIT 1", (key,))
                if cur.fetchone():
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
                cid = self._canonical_from_header(sid)
                if cid is None:
                    continue
                cur.execute("INSERT OR IGNORE INTO sequences (id, sequence, length, dataset) VALUES (?, ?, ?, ?)",
                            (cid, seq, len(seq), dataset))
                alias_entries.append((cid, sid))
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
                cid = self._canonical_from_header(sid) if sid is not None else None
                key = cid or sid
                # ensure this id was provided by the user (exists in sequences)
                cur.execute("SELECT 1 FROM sequences WHERE id = ? LIMIT 1", (key,))
                if cur.fetchone():
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
                cid = self._canonical_from_header(sid) if sid is not None else None
                key = cid or sid
                # Only insert if the id exists in sequences
                cur.execute("SELECT 1 FROM sequences WHERE id = ? LIMIT 1", (key,))
                if not cur.fetchone():
                    skipped += 1
                    continue
                n_cid = None
                if nearest is not None:
                    n_cid = self._canonical_from_header(nearest)
                to_insert.append((key, ds, n_cid or nearest, identity))
            if to_insert:
                cur.executemany("INSERT OR REPLACE INTO distances (id, dataset, nearest, identity) VALUES (?, ?, ?, ?)", to_insert)
            conn.commit()
            if skipped:
                self.logger.info("[DB] Skipped %d distance entries because ids are not present in sequences table", skipped)

    def insert_colors(self, color_entries):
        """Insert or replace color entries.

        color_entries: iterable of (id, color, source)
        """
        # Only persist colors for ids that exist in the sequences table
        with self.connect() as conn:
            cur = conn.cursor()
            to_insert = []
            skipped = 0
            for entry in color_entries:
                # accept (id, color, source) or (id, color, source, dataset)
                if len(entry) == 3:
                    sid, color, source = entry
                    ds = ''
                else:
                    sid, color, source, ds = entry
                    ds = ds or ''
                cid = self._canonical_from_header(sid) if sid is not None else None
                key = cid or sid
                cur.execute("SELECT 1 FROM sequences WHERE id = ? LIMIT 1", (key,))
                if cur.fetchone():
                    to_insert.append((key, color, source, ds))
                else:
                    skipped += 1
            if to_insert:
                cur.executemany("INSERT OR REPLACE INTO colors (id, color, source, dataset) VALUES (?, ?, ?, ?)", to_insert)
            conn.commit()
            if skipped:
                self.logger.info("[DB] Skipped %d color entries because ids are not present in sequences table", skipped)

    def get_color_map(self, ids=None):
        """Return a dict id->color for provided ids or for all if ids is None."""
        with self.connect() as conn:
            cur = conn.cursor()
            if ids:
                placeholders = ','.join('?' for _ in ids)
                cur.execute(f"SELECT id, color FROM colors WHERE id IN ({placeholders})", tuple(ids))
            else:
                cur.execute("SELECT id, color FROM colors")
            rows = cur.fetchall()
        return {r[0]: r[1] for r in rows}


