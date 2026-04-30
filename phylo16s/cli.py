import argparse
from pathlib import Path
from phylo16s.pipeline import qc, derep, classify, novelty, report, tree
from phylo16s.db.interface import Database
import logging


def _configure_logging(outdir: str, level=logging.INFO):
    """Configure root logging to stream to console and write a debug log to outdir."""
    logger = logging.getLogger()
    logger.setLevel(level)
    # remove existing handlers to avoid duplicate logs in repeated runs
    for h in list(logger.handlers):
        logger.removeHandler(h)
    fmt = logging.Formatter("%(asctime)s %(levelname)s: %(message)s")
    ch = logging.StreamHandler()
    ch.setLevel(level)
    ch.setFormatter(fmt)
    logger.addHandler(ch)
    try:
        import os
        fh = logging.FileHandler(os.path.join(outdir, 'phylo16s.log'))
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(fmt)
        logger.addHandler(fh)
    except Exception:
        logger.warning("Could not set up file logging in %s", outdir)


logger = logging.getLogger(__name__)


def run_pipeline(args):
    """Run the full workflow.

    The pipeline expects a reference fasta (e.g. greengene2) passed with
    --ref. A sqlite database path is passed with --db to store metadata.

    The workflow:
    - initialise sqlite database
    - run QC and dereplication on user input
    - classify user sequences against the provided reference fasta
    - detect novelty against the reference set
    - add only the user sequences into a persistent tree/alignment kept in
      the output directory (the tree/alignment are initialised from the
      reference fasta on first run)
    - generate a simple report
    """
    db = Database(args.db)
    db.initialise()

    # preflight: report input counts

    def count_fasta(path):
        try:
            # lazy-import to avoid module-level dependency issues
            from utils.fasta import read_fasta
            return sum(1 for _ in read_fasta(path))
        except Exception:
            return 0

    input_count = count_fasta(args.input)
    logger.info("[PRE-FLIGHT] Input sequences: %d", input_count)

    # run QC with configurable thresholds
    qc_out = qc.run_qc(args.input, args.out, min_len=getattr(args, 'min_len', 1200), max_n=getattr(args, 'max_n', 5))
    qc_count = count_fasta(qc_out)
    logger.info("[PRE-FLIGHT] After QC: %d sequences", qc_count)

    derep_out = derep.run_derep(qc_out, args.out)
    derep_count = count_fasta(derep_out)
    logger.info("[PRE-FLIGHT] After dereplication: %d unique sequences", derep_count)

    # map user-provided dereplicated IDs to short IDs (3 letters + 2 digits)
    try:
        mapped_derep = Path(args.out) / 'derep_short.fasta'
        from utils.fasta import read_fasta, write_fasta
        used_ids = set(db.get_all_ids())
        orig_to_short = {}
        mapped_records = []
        for h, s in read_fasta(derep_out):
            short = db.generate_short_id(h, used_ids)
            mapped_records.append((short, s))
            orig_to_short[h] = short
        write_fasta(mapped_records, str(mapped_derep))
        logger.info("[DB] Mapped %d user sequence IDs to short IDs and wrote %s", len(mapped_records), mapped_derep)
        # insert mapped records into DB (INSERT OR IGNORE inside)
        db.insert_sequences(mapped_records, dataset='user')
        # write mapping file for user's reference
        try:
            map_path = Path(args.out) / 'user_id_map.tsv'
            with open(map_path, 'w') as m:
                m.write('short_id\toriginal_header\n')
                for orig, short in orig_to_short.items():
                    m.write(f"{short}\t{orig}\n")
            logger.info("[DB] Wrote user id mapping to %s", map_path)
        except Exception as e:
            logger.warning("[DB] Could not write user id mapping file: %s", e)
    except Exception as e:
        logger.warning("[DB] Failed to generate short IDs for user sequences: %s", e)
        mapped_derep = derep_out

    # classification and novelty use the reference fasta provided by --ref
    class_out = classify.run_classification(str(mapped_derep), args.out, ref_fasta=args.ref, taxa_tsv=getattr(args, 'taxa', None), threads=getattr(args, 'threads', None))
    # if classification produced only a header or is missing, try to recover
    # a simple taxonomy table from a matches.tsv produced in the outdir
    import os
    try:
        def _is_header_only(path):
            # helper to robustly map various id/header forms to short ids
            def _try_map_key(key):
                if not key:
                    return None
                # direct match
                if key in orig_to_short:
                    return orig_to_short[key]
                lk = key.lower()
                for k, v in orig_to_short.items():
                    if k.lower() == lk:
                        return v
                try:
                    cid = db._canonical_from_header(key)
                except Exception:
                    cid = key
                if cid in orig_to_short:
                    return orig_to_short[cid]
                if '|' in key:
                    last = key.split('|')[-1]
                    if last in orig_to_short:
                        return orig_to_short[last]
                # substring match
                for k, v in orig_to_short.items():
                    if k and key and (k in key or key in k):
                        return v
                return None

            try:
                with open(path) as f:
                    lines = [l for l in f.readlines() if l.strip()]
                return len(lines) <= 1
            except Exception:
                return True

        if _is_header_only(class_out):
            matches_path = os.path.join(args.out, 'matches.tsv')
            if os.path.exists(matches_path):
                logger.info("[CLASSIFY] taxonomy output %s empty; regenerating from %s", class_out, matches_path)
                with open(matches_path) as inp, open(class_out, 'w') as out:
                    out.write("ID\tBestHit\tIdentity\tTaxon\tConfidence\n")
                    for line in inp:
                        parts = line.strip().split("\t")
                        if len(parts) < 3:
                            continue
                        qid = parts[0]
                        sid = parts[1]
                        identity = parts[2]
                        out.write(f"{qid}\t{sid}\t{identity}\tNA\tNA\n")
                logger.info("[CLASSIFY] Wrote regenerated taxonomy to %s", class_out)
    except Exception as e:
        logger.warning("[CLASSIFY] Warning: failed to regenerate taxonomy from matches.tsv: %s", e)
    novelty_out = novelty.run_novelty(str(mapped_derep), args.ref, args.out, threads=getattr(args, 'threads', None))

    # (DB persistence already performed earlier)

    # run similarity computations against preloaded datasets (and gg2 entries from DB)
    sim_long = None
    sim_pivot = None
    try:
        from phylo16s.pipeline import similarity
        sim_long, sim_pivot = similarity.run_similarity(str(mapped_derep), args.out, args.ref, db, threads=getattr(args, 'threads', None))
        logger.info("[SIM] Wrote similarity files: %s, %s", sim_long, sim_pivot)
    except Exception as e:
        logger.warning("[SIM] Warning: similarity step failed: %s", e)

    # maintain an ever-growing tree/alignment in the output directory
    if getattr(args, 'skip_tree', False):
        logger.info("[TREE] Skipping tree/alignment update as requested (--skip-tree)")
    else:
        # optionally force rebuild backbone by removing current_alignment
        if getattr(args, 'rebuild_backbone', False):
            import os
            ca = os.path.join(args.out, 'current_alignment.fasta')
            ct = os.path.join(args.out, 'current_tree.nwk')
            for p in (ca, ct):
                if os.path.exists(p):
                    os.remove(p)
                    logger.info("[TREE] Removed existing %s to force rebuild", p)

        tree.initialise_or_update_tree(ref_fasta=args.ref,
                                       user_fasta=str(mapped_derep),
                                       outdir=args.out,
                                       db=db,
                                       threads=getattr(args, 'threads', None))

    report.generate_report(
        str(mapped_derep),
        class_out,
        novelty_out,
        args.out,
        similarity_pivot=sim_pivot
    )

    # generate iTOL color files based on taxonomy output
    try:
        from phylo16s.pipeline import itol

        # prefer to reuse combined taxonomy from a prior preload if provided
        combined_tax_path = None
        user_color_csv = getattr(args, 'user_colors', None)
        preload_dir = getattr(args, 'preload_dir', None)
        if preload_dir:
            p = Path(preload_dir)
            # Prefer explicit combined taxonomy files in the preload dir, but
            # allow searching recursively in case the user passed a parent or DB
            # directory. Try multiple candidate filenames.
            cand = p / 'preload_combined_taxonomy.tsv'
            if not cand.exists():
                cand = p / 'combined_taxonomy.tsv'
            # if still not found, search recursively inside the provided path
            if not cand.exists():
                try:
                    found = next(p.rglob('preload_combined_taxonomy.tsv'), None)
                    if not found:
                        found = next(p.rglob('combined_taxonomy.tsv'), None)
                    if found:
                        cand = found
                except Exception:
                    pass
            if cand.exists():
                combined_tax_path = cand
                logger.info("[ITOL] Using combined taxonomy from preload dir: %s", combined_tax_path)
                # try to find a user colors CSV saved in the preload dir if none supplied
                if not user_color_csv:
                    for cname in ('user_colors.csv', 'colors.csv', 'user_color_csv.csv'):
                        cc = p / cname
                        if cc.exists():
                            user_color_csv = str(cc)
                            logger.info("[ITOL] Using user colors CSV from preload dir: %s", user_color_csv)
                            break

        # if no preload combined taxonomy found, build combined_tax_path here
        if combined_tax_path is None:
            combined_tax_path = Path(args.out) / 'combined_taxonomy.tsv'
            try:
                with open(combined_tax_path, 'w') as out_tax:
                    out_tax.write('ID\tTaxon\tConfidence\n')
                    try:
                        with db.connect() as conn:
                            cur = conn.cursor()
                            cur.execute("SELECT id, taxonomy, confidence FROM taxonomy")
                            for rid, tax, conf in cur.fetchall():
                                out_tax.write(f"{rid}\t{tax if tax is not None else 'NA'}\t{conf if conf is not None else 'NA'}\n")
                    except Exception:
                        pass

                    try:
                        with open(class_out) as t:
                            next(t, None)
                            for line in t:
                                parts = line.strip().split('\t')
                                if len(parts) < 4:
                                    continue
                                qid = parts[0]
                                tax = parts[3] if len(parts) > 3 else 'NA'
                                conf = parts[4] if len(parts) > 4 else 'NA'
                                out_tax.write(f"{qid}\t{tax if tax is not None else 'NA'}\t{conf if conf is not None else 'NA'}\n")
                    except Exception:
                        pass
            except Exception:
                combined_tax_path = class_out

        # If preload_dir provided and contains a combined taxonomy, merge it with
        # this run's classification and regenerate iTOL files from the merged
        # combined taxonomy so legends and counts are correct and IDs (short ids)
        # remain synchronized.
        if preload_dir:
            p = Path(preload_dir)
            cand = p / 'preload_combined_taxonomy.tsv'
            if not cand.exists():
                cand = p / 'combined_taxonomy.tsv'

            merged = {}
            order = []
            # try to load any id_map produced during preload so we can map
            # original/long ids from the preload combined taxonomy to the
            # short IDs used throughout the DB/tree.
            orig_to_short = {}
            # If the preload combined taxonomy file is missing from the provided
            # preload_dir, fall back to using taxonomy rows from the DB so that
            # preloaded dataset taxonomies are still included in the merged
            # combined file. This covers cases where the user passed a DB path
            # or the preload outdir doesn't contain the pre-generated files.
            if not cand.exists():
                try:
                    # Only load taxonomy for sequences that belong to non-user datasets
                    # (i.e. preloaded datasets). This avoids pulling the entire gg2
                    # reference taxonomy into the merged file which would leave
                    # reference IDs (not short dataset IDs) in the combined output.
                    with db.connect() as conn:
                        cur = conn.cursor()
                        cur.execute("SELECT id FROM sequences WHERE dataset IS NOT NULL AND dataset != 'user'")
                        seq_rows = [r[0] for r in cur.fetchall()]
                        if seq_rows:
                            placeholders = ','.join('?' for _ in seq_rows)
                            cur.execute(f"SELECT id, taxonomy, confidence FROM taxonomy WHERE id IN ({placeholders})", tuple(seq_rows))
                            for rid, tax, conf in cur.fetchall():
                                merged[rid] = (tax if tax is not None else 'NA', conf if conf is not None else 'NA')
                                order.append(rid)
                    logger.info("[ITOL] No preload combined file found in %s; loaded %d dataset taxonomy rows from DB for merging", str(p), len(order))
                except Exception:
                    # leave merged empty and continue to attempt loading any
                    # classification entries below
                    merged = {}
                    order = []
            try:
                # look for a dataset id map file in the preload dir or nearby. Use
                # a recursive search so callers can pass a parent/DB path.
                id_map_candidates = []
                try:
                    id_map_candidates = list(p.glob('*_id_map.tsv'))
                except Exception:
                    id_map_candidates = []
                # if none directly in p, try recursive search (rglob) and parent dirs
                if not id_map_candidates:
                    try:
                        id_map_candidates = list(p.rglob('*_id_map.tsv'))
                    except Exception:
                        id_map_candidates = []
                if not id_map_candidates:
                    try:
                        id_map_candidates = list(p.parent.rglob('*_id_map.tsv'))
                    except Exception:
                        id_map_candidates = []
                if id_map_candidates:
                    # prefer a file that contains 'preload' or that matches the
                    # preload dir name; otherwise take the first candidate
                    id_map_path = None
                    for c in id_map_candidates:
                        name = str(c.name).lower()
                        if 'preload' in name or (p.name.lower() + '_id_map.tsv') in name:
                            id_map_path = c
                            break
                    if id_map_path is None:
                        id_map_path = id_map_candidates[0]
                    try:
                        with open(id_map_path) as im:
                            next(im, None)
                            for line in im:
                                parts = line.strip().split('\t')
                                if len(parts) < 2:
                                    continue
                                short, orig = parts[0], parts[1]
                                # store multiple lookup keys to make mapping robust
                                orig_to_short[orig] = short
                                try:
                                    cid = db._canonical_from_header(orig)
                                except Exception:
                                    cid = orig
                                if cid and cid not in orig_to_short:
                                    orig_to_short[cid] = short
                                # also map the last pipe-delimited token (common in headers)
                                if '|' in orig:
                                    last = orig.split('|')[-1]
                                    if last and last not in orig_to_short:
                                        orig_to_short[last] = short
                    except Exception:
                        orig_to_short = {}

                if cand.exists():
                    try:
                        with open(cand) as f:
                            next(f, None)
                            for line in f:
                                parts = line.rstrip('\n').split('\t')
                                if not parts:
                                    continue
                                iid = parts[0]
                                tax = parts[1] if len(parts) > 1 else 'NA'
                                conf = parts[2] if len(parts) > 2 else 'NA'
                                # map various id/header forms to short ids (inline to avoid scope issues)
                                mapped_iid = None
                                if iid in orig_to_short:
                                    mapped_iid = orig_to_short[iid]
                                else:
                                    lk = iid.lower()
                                    for k, v in orig_to_short.items():
                                        try:
                                            if k.lower() == lk:
                                                mapped_iid = v
                                                break
                                        except Exception:
                                            continue
                                if mapped_iid is None:
                                    try:
                                        cid = db._canonical_from_header(iid)
                                    except Exception:
                                        cid = iid
                                    if cid in orig_to_short:
                                        mapped_iid = orig_to_short[cid]
                                if mapped_iid is None and '|' in iid:
                                    last = iid.split('|')[-1]
                                    if last in orig_to_short:
                                        mapped_iid = orig_to_short[last]
                                if mapped_iid is None:
                                    for k, v in orig_to_short.items():
                                        if k and iid and (k in iid or iid in k):
                                            mapped_iid = v
                                            break
                                if not mapped_iid:
                                    mapped_iid = iid
                                # track mapping diagnostics
                                try:
                                    import re
                                    is_short = bool(re.match(r'^[A-Z]{3}[0-9]{2}$', str(mapped_iid)))
                                except Exception:
                                    is_short = False
                                if not hasattr(globals(), '_map_stats'):
                                    globals()['_map_stats'] = {'mapped': 0, 'unmapped': 0, 'examples': []}
                                if is_short:
                                    globals()['_map_stats']['mapped'] += 1
                                else:
                                    globals()['_map_stats']['unmapped'] += 1
                                    if len(globals()['_map_stats']['examples']) < 8:
                                        globals()['_map_stats']['examples'].append(iid)
                                merged[mapped_iid] = (tax, conf)
                                order.append(mapped_iid)
                    except Exception:
                        merged = {}
                        order = []

                # after attempting to read and map preload ids, log mapping stats
                try:
                    stats = globals().get('_map_stats', None)
                    if stats:
                        logger.info("[ITOL] Preload id mapping: mapped=%d unmapped=%d examples_unmapped=%s", stats.get('mapped', 0), stats.get('unmapped', 0), ','.join(stats.get('examples', [])))
                except Exception:
                    pass

            except Exception:
                # If any unexpected error occurs while attempting to read the
                # preload combined taxonomy or id_map, fall back to empty
                # merged structures so classification entries can still be
                # integrated below.
                orig_to_short = {}
                merged = {}
                order = []

            # integrate classification output: override/add entries for qids
            try:
                with open(class_out) as t:
                    next(t, None)
                    for line in t:
                        parts = line.strip().split('\t')
                        if len(parts) < 4:
                            continue
                        qid = parts[0]
                        tax = parts[3] if len(parts) > 3 else 'NA'
                        conf = parts[4] if len(parts) > 4 else 'NA'
                        if qid not in merged:
                            order.append(qid)
                        merged[qid] = (tax, conf)
            except Exception:
                pass

            # write merged combined taxonomy for run
            merged_path = Path(args.out) / 'combined_taxonomy.tsv'
            with open(merged_path, 'w') as out_tax:
                out_tax.write('ID\tTaxon\tConfidence\n')
                for iid in order:
                    tax, conf = merged.get(iid, ('NA', 'NA'))
                    out_tax.write(f"{iid}\t{tax if tax is not None else 'NA'}\t{conf if conf is not None else 'NA'}\n")

            # copy user color CSV from preload dir if present and user didn't supply one
            if not user_color_csv:
                import shutil
                for cname in ('user_colors.csv', 'colors.csv', 'user_color_csv.csv'):
                    cc = p / cname
                    if cc.exists():
                        try:
                            shutil.copy(str(cc), str(Path(args.out) / cc.name))
                            user_color_csv = str(Path(args.out) / cc.name)
                            logger.info("[ITOL] Copied user colors CSV from preload dir: %s", user_color_csv)
                            break
                        except Exception:
                            pass

            # If no explicit per-id user color CSV was found, try to extract per-id
            # user colors from the preload itol_combined_colors.csv (if present).
            # This file contains per-id columns: id,phylum_color,family_color,genus_color,user_color
            # We'll build a simple 'id,color' CSV from the user_color column (or fall
            # back to genus/family/phylum color if user_color is empty) and pass it
            # to the itol generator so preloaded sequence-specific colors are retained.
            if not user_color_csv:
                try:
                    import csv
                    preload_itol_comb = p / 'itol_combined_colors.csv'
                    if preload_itol_comb.exists():
                        out_uc = Path(args.out) / 'user_colors_from_preload.csv'
                        with open(preload_itol_comb) as pre, open(out_uc, 'w', newline='') as outc:
                            reader = csv.DictReader(pre)
                            writer = csv.writer(outc)
                            writer.writerow(['id', 'color'])
                            for row in reader:
                                uid = row.get('id') or row.get('ID') or row.get('Id')
                                if not uid:
                                    continue
                                # prefer explicit user_color column
                                user_col = row.get('user_color') or row.get('user_color'.upper()) or row.get('user_color'.title()) or ''
                                # fallback to genus/family/phylum colors if user_color missing
                                if not user_col:
                                    user_col = row.get('genus_color') or row.get('family_color') or row.get('phylum_color') or ''
                                if not user_col:
                                    continue
                                # map preload id to short id if id_map available
                                mapped_id = orig_to_short.get(uid, None)
                                if mapped_id is None:
                                    try:
                                        cid = db._canonical_from_header(uid)
                                    except Exception:
                                        cid = uid
                                    mapped_id = orig_to_short.get(cid, uid)
                                writer.writerow([mapped_id, user_col])
                        user_color_csv = str(out_uc)
                        logger.info("[ITOL] Built user color CSV from preload itol_combined_colors.csv -> %s", user_color_csv)
                except Exception:
                    # non-fatal; we'll proceed without per-id user colors
                    pass

            # regenerate iTOL files from merged taxonomy so legends/counts are correct
            itol.generate_itol_colors(str(merged_path), args.out, user_color_csv=user_color_csv, id_map=orig_to_short)
            logger.info("[ITOL] Regenerated iTOL color files in %s from merged preload taxonomy", args.out)

            # Produce an ITOL-format combined colorstrip file (itol_combined_colors.itol)
            # by reading the itol_combined_colors.csv produced by the generator and
            # converting to a DATASET_COLORSTRIP file mapping each id to a final
            # color (prefer user_color, then genus/family/phylum colors).
            try:
                import csv
                out_dir = Path(args.out)
                csv_path = out_dir / 'itol_combined_colors.csv'
                itol_path = out_dir / 'itol_combined_colors.itol'
                if csv_path.exists():
                    with open(csv_path) as cf:
                        reader = csv.DictReader(cf)
                        rows = [r for r in reader]

                    with open(itol_path, 'w') as of:
                        of.write('DATASET_COLORSTRIP\n')
                        of.write('SEPARATOR COMMA\n')
                        of.write('DATASET_LABEL,Combined colors\n')
                        # choose a neutral dataset color if needed
                        of.write('COLOR,#AAAAAA\n')
                        of.write('MARGIN,5\n')
                        of.write('SHOW_INTERNAL,0\n')
                        of.write('DATA\n')
                        for r in rows:
                            uid = r.get('id') or r.get('ID') or r.get('Id')
                            if not uid:
                                continue
                            # prefer explicit user color if present
                            user_col = (r.get('user_color') or '').strip()
                            final = user_col or r.get('genus_color') or r.get('family_color') or r.get('phylum_color') or ''
                            if not final:
                                continue
                            # map preload/original ids to short ids if mapping is available
                            try:
                                mapped_uid = orig_to_short.get(uid, uid)
                            except Exception:
                                mapped_uid = uid
                            of.write(f"{mapped_uid},{final}\n")
                    logger.info("[ITOL] Wrote combined ITOL colorstrip to %s", str(itol_path))
                else:
                    logger.warning("[ITOL] Expected %s produced by color generator but not found", str(csv_path))
            except Exception as e:
                logger.warning("[ITOL] Failed to write combined ITOL colorstrip: %s", e)
            # Build a single consolidated ITOL file (taxon + dataset membership)
            try:
                out_itol = Path(args.out) / 'itol_combined_colors.itol'
                if out_itol.exists():
                    taxon_text = out_itol.read_text()
                else:
                    taxon_text = ''

                # no preload directory available here, so only user ids vs others
                preload_ids = set()
                preload_color = None
                user_ids = set()
                try:
                    uid_map = Path(args.out) / 'user_id_map.tsv'
                    if uid_map.exists():
                        with open(uid_map) as uf:
                            next(uf, None)
                            for l in uf:
                                parts = l.strip().split('\t')
                                if parts and parts[0]:
                                    user_ids.add(parts[0])
                except Exception:
                    user_ids = set()

                run_color = '#ff3333'
                other_color = '#cccccc'

                # attempt to map any combined_taxonomy ids to short ids using
                # orig_to_short mapping if available so iTOL dataset ids match
                # the tree tip labels. orig_to_short may not exist in all code
                # paths, so guard against NameError.
                try:
                    id_map = orig_to_short
                except NameError:
                    id_map = {}

                def _map_iid(iid):
                    if not id_map:
                        return iid
                    if iid in id_map:
                        return id_map[iid]
                    lk = iid.lower()
                    for k, v in id_map.items():
                        try:
                            if k.lower() == lk:
                                return v
                        except Exception:
                            continue
                    if '|' in iid:
                        last = iid.split('|')[-1]
                        if last in id_map:
                            return id_map[last]
                    for k, v in id_map.items():
                        if k and iid and (k in iid or iid in k):
                            return v
                    try:
                        cid = db._canonical_from_header(iid)
                        if cid in id_map:
                            return id_map[cid]
                    except Exception:
                        pass
                    return iid

                combined_tax = Path(args.out) / 'combined_taxonomy.tsv'
                ids_in_order = []
                if combined_tax.exists():
                    with open(combined_tax) as ct:
                        next(ct, None)
                        for l in ct:
                            iid = l.strip().split('\t')[0]
                            if iid:
                                ids_in_order.append(_map_iid(iid))

                ds_lines = []
                ds_lines.append('')
                ds_lines.append('DATASET_COLORSTRIP')
                ds_lines.append('SEPARATOR COMMA')
                ds_lines.append('DATASET_LABEL,Dataset membership')
                ds_lines.append(f'COLOR,#AAAAAA')
                ds_lines.append('MARGIN,5')
                ds_lines.append('SHOW_INTERNAL,0')
                ds_lines.append('DATA')
                for iid in ids_in_order:
                    if iid in preload_ids:
                        ds_lines.append(f"{iid},{preload_color or '#1f78b4'}")
                    elif iid in user_ids:
                        ds_lines.append(f"{iid},{run_color}")
                    else:
                        ds_lines.append(f"{iid},{other_color}")

                # write consolidated file and cleanup other itol files
                out_itol.write_text('\n'.join([taxon_text.strip()] + ds_lines) + '\n')
                for f in Path(args.out).glob('itol_*.itol'):
                    try:
                        if f.name != out_itol.name:
                            f.unlink()
                    except Exception:
                        pass
                logger.info("[ITOL] Wrote single consolidated ITOL file to %s", str(out_itol))
                try:
                    # attempt to post-process with a robust rebuild script if present
                    import subprocess, sys
                    script = Path('/Users/nicholas/Git/rebuild_itol.py')
                    if script.exists():
                        subprocess.run([sys.executable, str(script)], check=False)
                        logger.info("[ITOL] Ran rebuild_itol.py to validate consolidated ITOL file")
                except Exception:
                    pass
            except Exception as e:
                logger.warning("[ITOL] Failed to build consolidated ITOL file: %s", e)
            # If a preload itol_combined_colors.csv exists, ensure its per-id
            # colours (short ids) are present in the run itol file. Append any
            # missing preload short ids and their colours so the run itol file
            # contains preload + user entries.
            try:
                preload_csv = p / 'itol_combined_colors.csv'
                out_itol = Path(args.out) / 'itol_combined_colors.itol'
                if preload_csv.exists() and out_itol.exists():
                    import csv
                    # read existing ids in the run itol file
                    existing_ids = set()
                    with open(out_itol) as of:
                        for ln in of:
                            ln = ln.strip()
                            if not ln or ln.startswith('#'):
                                continue
                            if ',' in ln:
                                a, b = ln.split(',', 1)
                                # only consider data lines (assume header lines earlier)
                                if a and any(c.isalnum() for c in a):
                                    existing_ids.add(a)

                    # iterate preload CSV and append missing ids
                    appended = 0
                    with open(preload_csv) as pf, open(out_itol, 'a') as of:
                        reader = csv.DictReader(pf)
                        for row in reader:
                            uid = (row.get('id') or row.get('ID') or row.get('Id') or '').strip()
                            if not uid:
                                continue
                            # prefer explicit user_color column, fallback to genus/family/phylum
                            user_col = (row.get('user_color') or '').strip()
                            if not user_col:
                                user_col = (row.get('genus_color') or row.get('family_color') or row.get('phylum_color') or '').strip()
                            if not user_col:
                                continue
                            # if preload id maps to a short id in orig_to_short, use that
                            mapped = orig_to_short.get(uid, uid)
                            if mapped in existing_ids:
                                continue
                            of.write(f"{mapped},{user_col}\n")
                            appended += 1
                    if appended:
                        logger.info("[ITOL] Appended %d preload ids to %s", appended, str(out_itol))
                    # Now build a single, consolidated ITOL file that contains both
                    # the taxon colorstrip (already in out_itol) and a dataset
                    # membership colorstrip (preload vs run vs others). Overwrite
                    # the existing itol_combined_colors.itol with the combined
                    # content so the run outdir ends up with one ITOL file only.
                    try:
                        # read current taxon colorstrip content
                        taxon_text = out_itol.read_text()
                        # collect preload ids and preload header color if available
                        preload_ids = set()
                        preload_color = None
                        try:
                            preload_itol = p / 'itol_dataset_preload.itol'
                            if preload_itol.exists():
                                with open(preload_itol) as pf:
                                    for ln in pf:
                                        ln = ln.strip()
                                        if not ln:
                                            continue
                                        if ln.startswith('COLOR,') and preload_color is None:
                                            try:
                                                preload_color = ln.split(',', 1)[1]
                                            except Exception:
                                                preload_color = None
                                        if ',' in ln and not ln.startswith('DATASET_') and not ln.startswith('SEPARATOR') and not ln.startswith('DATASET_LABEL') and not ln.startswith('COLOR') and not ln.startswith('MARGIN') and not ln.startswith('SHOW_INTERNAL') and not ln.startswith('DATA'):
                                            a, b = ln.split(',', 1)
                                            preload_ids.add(a)
                        except Exception:
                            preload_ids = set()
                        # collect user/run ids from user_id_map.tsv if present
                        user_ids = set()
                        try:
                            uid_map = Path(args.out) / 'user_id_map.tsv'
                            if uid_map.exists():
                                with open(uid_map) as uf:
                                    next(uf, None)
                                    for l in uf:
                                        parts = l.strip().split('\t')
                                        if parts and parts[0]:
                                            user_ids.add(parts[0])
                        except Exception:
                            user_ids = set()

                        # determine dataset colors: prefer preload header color, use a distinct run color, and a fallback for others
                        run_color = '#ff3333'
                        other_color = '#cccccc'
                        if preload_color is None:
                            preload_color = '#1f78b4'

                        # read combined taxonomy ordering so dataset strip matches tree ordering
                        try:
                            id_map = orig_to_short
                        except NameError:
                            id_map = {}

                        def _map_iid_local(iid):
                            if not id_map:
                                return iid
                            if iid in id_map:
                                return id_map[iid]
                            lk = iid.lower()
                            for k, v in id_map.items():
                                try:
                                    if k.lower() == lk:
                                        return v
                                except Exception:
                                    continue
                            if '|' in iid:
                                last = iid.split('|')[-1]
                                if last in id_map:
                                    return id_map[last]
                            for k, v in id_map.items():
                                if k and iid and (k in iid or iid in k):
                                    return v
                            try:
                                cid = db._canonical_from_header(iid)
                                if cid in id_map:
                                    return id_map[cid]
                            except Exception:
                                pass
                            return iid

                        combined_tax = Path(args.out) / 'combined_taxonomy.tsv'
                        ids_in_order = []
                        if combined_tax.exists():
                            with open(combined_tax) as ct:
                                next(ct, None)
                                for l in ct:
                                    iid = l.strip().split('\t')[0]
                                    if iid:
                                        ids_in_order.append(_map_iid_local(iid))

                        # build dataset membership section
                        ds_lines = []
                        ds_lines.append('')
                        ds_lines.append('DATASET_COLORSTRIP')
                        ds_lines.append('SEPARATOR COMMA')
                        ds_lines.append('DATASET_LABEL,Dataset membership')
                        ds_lines.append(f'COLOR,#AAAAAA')
                        ds_lines.append('MARGIN,5')
                        ds_lines.append('SHOW_INTERNAL,0')
                        ds_lines.append('DATA')
                        for iid in ids_in_order:
                            if iid in preload_ids:
                                ds_lines.append(f"{iid},{preload_color}")
                            elif iid in user_ids:
                                ds_lines.append(f"{iid},{run_color}")
                            else:
                                ds_lines.append(f"{iid},{other_color}")

                        # write combined file (overwrite existing itol file)
                        out_itol.write_text('\n'.join([taxon_text.strip()] + ds_lines) + '\n')
                        logger.info("[ITOL] Wrote single consolidated ITOL file to %s", str(out_itol))
                        try:
                            import subprocess, sys
                            script = Path('/Users/nicholas/Git/rebuild_itol.py')
                            if script.exists():
                                subprocess.run([sys.executable, str(script)], check=False)
                                logger.info("[ITOL] Ran rebuild_itol.py to validate consolidated ITOL file")
                        except Exception:
                            pass
                        # remove other per-track itol files so only the single file remains
                        for f in Path(args.out).glob('itol_*.itol'):
                            try:
                                if f.name != out_itol.name:
                                    f.unlink()
                            except Exception:
                                pass
                    except Exception as e:
                        logger.warning("[ITOL] Failed to build consolidated ITOL file: %s", e)
            except Exception as e:
                logger.warning("[ITOL] Failed to merge preload itol colours into run itol: %s", e)
        else:
            # no preload dir: generate fresh files
            try:
                id_map = orig_to_short
            except NameError:
                id_map = None
            itol.generate_itol_colors(str(combined_tax_path), args.out, user_color_csv=user_color_csv, id_map=id_map)
            logger.info("[ITOL] Generated iTOL color files in %s (source=%s)", args.out, combined_tax_path)
            # also produce an ITOL-format combined colorstrip for convenience
            try:
                import csv
                out_dir = Path(args.out)
                csv_path = out_dir / 'itol_combined_colors.csv'
                itol_path = out_dir / 'itol_combined_colors.itol'
                if csv_path.exists():
                    with open(csv_path) as cf:
                        reader = csv.DictReader(cf)
                        rows = [r for r in reader]
                    with open(itol_path, 'w') as of:
                        of.write('DATASET_COLORSTRIP\n')
                        of.write('SEPARATOR COMMA\n')
                        of.write('DATASET_LABEL,Combined colors\n')
                        of.write('COLOR,#AAAAAA\n')
                        of.write('MARGIN,5\n')
                        of.write('SHOW_INTERNAL,0\n')
                        of.write('DATA\n')
                        for r in rows:
                            uid = r.get('id') or r.get('ID') or r.get('Id')
                            if not uid:
                                continue
                            user_col = (r.get('user_color') or '').strip()
                            final = user_col or r.get('genus_color') or r.get('family_color') or r.get('phylum_color') or ''
                            if not final:
                                continue
                            of.write(f"{uid},{final}\n")
                    logger.info("[ITOL] Wrote combined ITOL colorstrip to %s", str(itol_path))
                else:
                    logger.warning("[ITOL] Expected %s produced by color generator but not found", str(csv_path))
            except Exception as e:
                logger.warning("[ITOL] Failed to write combined ITOL colorstrip: %s", e)
    except Exception as e:
        logger.warning("[ITOL] Warning: failed to generate iTOL files: %s", e)

    # Persist new sequences and classification into the database
    try:
        from utils.fasta import read_fasta

        # insert sequences from mapped derep fasta (IDs already short-mapped)
        seqs = [(h, s) for h, s in read_fasta(str(mapped_derep))]
        if seqs:
            db.insert_sequences(seqs)
            logger.info("[DB] Inserted %d sequences into DB", len(seqs))

        # insert taxonomy and distances
        tax_entries = []
        dist_entries = []
        with open(class_out) as t:
            next(t, None)
            for line in t:
                parts = line.strip().split("\t")
                if len(parts) < 3:
                    continue
                sid = parts[0]
                best = parts[1] if len(parts) > 1 else None
                identity = None
                try:
                    identity = float(parts[2])
                except Exception:
                    identity = None
                tax = parts[3] if len(parts) > 3 else None
                conf = None
                try:
                    conf = float(parts[4]) if len(parts) > 4 else None
                except Exception:
                    conf = None
                tax_entries.append((sid, tax, conf))
                # record distances with dataset label 'gg2' for classification hits
                dist_entries.append((sid, 'gg2', best, identity))

        if tax_entries:
            db.insert_taxonomy(tax_entries)
            logger.info("[DB] Inserted/updated taxonomy for %d ids", len(tax_entries))
        if dist_entries:
            db.insert_distances(dist_entries)
            logger.info("[DB] Inserted/updated distances for %d ids", len(dist_entries))
    except Exception as e:
        logger.warning("[DB] Warning: failed to persist to DB: %s", e)


def main():
    parser = argparse.ArgumentParser(prog="phylo16s")

    subparsers = parser.add_subparsers(dest="command")
    # run pipeline
    run = subparsers.add_parser("run")
    run.add_argument("--input", required=True)
    run.add_argument("--db", required=True)
    run.add_argument("-o", "--out", required=True, help="Output directory to write pipeline artifacts (itol, tree, matches)")
    run.add_argument("--ref", required=True, help="Reference fasta (e.g. greengene2) used for classification and initial tree")
    run.add_argument("--taxa", required=False, help="Optional taxa TSV mapping file (Feature ID<TAB>Taxon<TAB>Confidence) for annotating reference hits")
    run.add_argument("--min-len", dest='min_len', type=int, default=1200, help="Minimum sequence length to keep in QC (default: 1200)")
    run.add_argument("--max-n", dest='max_n', type=int, default=5, help="Maximum number of Ns allowed in sequence (default: 5)")
    run.add_argument("--threads", dest='threads', type=int, required=False, default=4, help="Number of threads to pass to vsearch (optional)")
    run.add_argument("--preload-dir", dest='preload_dir', required=False, help="Path to preload output directory to reuse combined taxonomy and colors")
    run.add_argument("--skip-tree", dest='skip_tree', action='store_true', help="Skip tree/alignment update step")
    run.add_argument("--rebuild-backbone", dest='rebuild_backbone', action='store_true', help="Force rebuild of backbone from DB or user sequences (remove current_alignment)")
    # preload DB with reference sequences (e.g., Hungate collection)
    preload = subparsers.add_parser("preload")
    preload.add_argument("--fasta", required=True, help="FASTA of reference sequences to preload into DB")
    preload.add_argument("--taxa", required=False, help="Optional taxa TSV matching the FASTA (FeatureID\tTaxon\tConfidence)")
    preload.add_argument("--colors", required=False, help="Optional CSV (id,color) to set colors for preloaded sequences")
    preload.add_argument("--db", required=True, help="Path to sqlite DB to preload into")
    preload.add_argument("--dataset", required=False, default='preload', help="Label/name for this preloaded dataset (e.g. Hungate)")
    preload.add_argument("-o", "--out", required=False, default='.', help="Output directory to write preload artifacts (itol, tree, matches)")
    preload.add_argument("--ref", required=False, help="Reference fasta (gg2) used for classification of preloaded dataset")
    preload.add_argument("--classify", dest='classify', action='store_true', help="Run classification of the preloaded fasta against --ref and insert distances into DB")
    preload.add_argument("--build-tree", dest='build_tree', action='store_true', help="Build a baseline tree/alignment from the preloaded dataset and DB sequences into --out")
    preload.add_argument("--dataset-color", dest='dataset_color', required=False, help="Hex color to use for the dataset iTOL colorstrip (e.g. #ff0000). If omitted a deterministic color is chosen")
    preload.add_argument("--threads", dest='threads', type=int, required=False, default=4, help="Number of threads to pass to vsearch for preload classification (optional)")
    preload.add_argument("--use-db-backbone", dest='use_db_backbone', action='store_true', help="When building a baseline tree during preload, include existing DB sequences as part of the backbone (default: False)")
    preload.add_argument("--db-backbone-dataset", dest='db_backbone_dataset', required=False, help="Optional dataset name in the DB to use when building the backbone (requires --use-db-backbone). If omitted and --use-db-backbone is set, all DB sequences will be used.")
    run.add_argument("--user-colors", dest='user_colors', required=False, help="Optional CSV of user-provided colors (header id,color) to override colors for specific sequences")

    args = parser.parse_args()

    if args.command == "run":
        # ensure output directory exists
        import os
        try:
            os.makedirs(args.out, exist_ok=True)
        except Exception:
            pass
        # configure logging now that outdir exists
        _configure_logging(args.out)
        run_pipeline(args)
    elif args.command == "preload":
        db = Database(args.db)
        db.initialise()
        outdir = getattr(args, 'out', '.')
        # ensure output directory exists for preload artifacts
        import os
        try:
            os.makedirs(outdir, exist_ok=True)
        except Exception:
            pass
        _configure_logging(outdir)
        logger.info("[PRELOAD] Starting preload of %s into %s", args.fasta, args.db)
        alias_entries, mapped_fasta = db.preload_from_files(args.fasta, taxa_tsv=getattr(args, 'taxa', None), color_csv=getattr(args, 'colors', None), source='preload', dataset=getattr(args, 'dataset', 'preload'), outdir=outdir)
        logger.info("[PRELOAD] Preloaded sequences from %s into %s", args.fasta, args.db)
        # write alias mapping to outdir for reference
        try:
            outdir = getattr(args, 'out', '.')
            map_path = Path(outdir) / f"{getattr(args,'dataset','preload')}_id_map.tsv"
            if alias_entries:
                with open(map_path, 'w') as m:
                    m.write('short_id\toriginal_header\n')
                    for short, orig in alias_entries:
                        m.write(f"{short}\t{orig}\n")
                logger.info("[PRELOAD] Wrote alias mapping to %s", map_path)
        except Exception as e:
            logger.warning("[PRELOAD] Could not write alias mapping file: %s", e)

        # determine dataset fasta to use for downstream steps. Prefer the mapped
        # fasta written by preload (short IDs) if available; otherwise export from DB.
        if mapped_fasta:
            ds_fasta = Path(mapped_fasta)
            exported = True
        else:
            ds_fasta = Path(outdir) / f"preload_{getattr(args,'dataset','preload')}_seqs.fasta"
            exported = False
            try:
                exported = db.export_dataset_fasta(getattr(args, 'dataset', 'preload'), str(ds_fasta))
            except Exception:
                exported = False

        # optionally classify the preloaded fasta against gg2 and insert distances
        if getattr(args, 'classify', False):
            if not getattr(args, 'ref', None):
                logger.info("[PRELOAD] --classify requested but --ref not provided; skipping classification")
            else:
                try:
                    input_for_classify = str(ds_fasta) if exported else args.fasta
                    logger.info("[PRELOAD] Classifying preloaded fasta %s against %s", input_for_classify, args.ref)
                    class_out = classify.run_classification(input_for_classify, outdir, ref_fasta=args.ref, taxa_tsv=getattr(args, 'taxa', None), threads=getattr(args, 'threads', None))
                    # parse classification and insert distances labeled with dataset
                    dist_entries = []
                    best_hits = set()
                    with open(class_out) as t:
                        next(t, None)
                        for line in t:
                            parts = line.strip().split("\t")
                            if len(parts) < 3:
                                continue
                            sid = parts[0]
                            best = parts[1] if len(parts) > 1 else None
                            identity = None
                            try:
                                identity = float(parts[2])
                            except Exception:
                                identity = None
                            # store as (sid, dataset, nearest, identity)
                            dist_entries.append((sid, getattr(args, 'dataset', 'preload'), best, identity))
                            if best and best != 'NA':
                                best_hits.add(best)

                    # Insert taxonomy rows for any gg2 best-hits found in this classification,
                    # but only those references — do not insert the entire gg2 taxa file.
                    if best_hits and getattr(args, 'taxa', None):
                        taxa_file = getattr(args, 'taxa')
                        import gzip
                        open_fn = gzip.open if str(taxa_file).endswith('.gz') else open
                        ref_tax_rows = []
                        # match by canonical id to avoid formatting mismatches
                        try:
                            needed_canon = set()
                            for b in best_hits:
                                try:
                                    cb = db._canonical_from_header(b)
                                except Exception:
                                    cb = b
                                if cb:
                                    needed_canon.add(cb)

                            with open_fn(taxa_file, 'rt') as tf:
                                header = tf.readline()
                                # if the header line is actually data, process it
                                if not ('Feature' in header or 'Taxon' in header):
                                    parts = header.strip().split('\t')
                                    if len(parts) >= 2:
                                        try:
                                            cid = db._canonical_from_header(parts[0])
                                        except Exception:
                                            cid = parts[0]
                                        if cid in needed_canon:
                                            conf = float(parts[2]) if len(parts) > 2 else None
                                            ref_tax_rows.append((cid, parts[1], conf, 'gg2'))
                                            needed_canon.discard(cid)
                                for line in tf:
                                    parts = line.strip().split('\t')
                                    if len(parts) < 2:
                                        continue
                                    try:
                                        cid = db._canonical_from_header(parts[0])
                                    except Exception:
                                        cid = parts[0]
                                    if cid in needed_canon:
                                        conf = float(parts[2]) if len(parts) > 2 else None
                                        ref_tax_rows.append((cid, parts[1], conf, 'gg2'))
                                        needed_canon.discard(cid)
                                    if not needed_canon:
                                        break
                        except Exception:
                            ref_tax_rows = []

                        if ref_tax_rows:
                            try:
                                db.insert_taxonomy(ref_tax_rows)
                                logger.info("[PRELOAD] Inserted %d gg2 taxonomy rows for matched references", len(ref_tax_rows))
                            except Exception as e:
                                logger.warning("[PRELOAD] Failed to insert gg2 taxonomy rows: %s", e)

                    if dist_entries:
                        # insert distances in chunks with periodic progress logs to avoid long silent pauses
                        total = len(dist_entries)
                        chunk = 1000
                        for i in range(0, total, chunk):
                            part = dist_entries[i:i+chunk]
                            db.insert_distances(part)
                            logger.info("[PRELOAD] Inserted distances %d-%d of %d for dataset %s", i+1, min(i+chunk, total), total, getattr(args, 'dataset', 'preload'))
                        logger.info("[PRELOAD] Inserted/updated distances for preloaded dataset %s (total %d)", getattr(args, 'dataset', 'preload'), total)
                except Exception as e:
                    logger.warning("[PRELOAD] Warning: classification during preload failed: %s", e)

        # write a simple iTOL colorstrip where all preloaded sequences get the same dataset color
        try:
            from phylo16s.pipeline import itol as _itol
            import csv
            dataset = getattr(args, 'dataset', 'preload')
            ds_color = getattr(args, 'dataset_color', None) or _itol._name_to_color(dataset)
            # export dataset ids
            ds_fasta = Path(outdir) / f"preload_{dataset}_seqs.fasta"
            exported = db.export_dataset_fasta(dataset, str(ds_fasta))
            if exported:
                from utils.fasta import read_fasta as _rf
                ids = [h for h, _ in _rf(str(ds_fasta))]
                itol_path = Path(outdir) / f"itol_dataset_{dataset}.itol"
                with open(itol_path, 'w') as f:
                    f.write('DATASET_COLORSTRIP\n')
                    f.write('SEPARATOR,COMMA\n')
                    f.write(f'DATASET_LABEL,{dataset} colors\n')
                    f.write(f'COLOR,{ds_color}\n')
                    f.write('MARGIN,5\n')
                    f.write('SHOW_INTERNAL,0\n')
                    f.write('DATA\n')
                    for id_ in ids:
                        f.write(f"{id_},{ds_color}\n")
                logger.info("[PRELOAD] Wrote iTOL dataset colorstrip to %s", itol_path)
        except Exception as e:
            logger.warning("[PRELOAD] Warning: failed to write dataset iTOL colorstrip: %s", e)

        # generate taxa-level iTOL color files for the preloaded DB so downstream
        # runs can reuse the same taxon palettes (phylum/family/genus)
        try:
            from phylo16s.pipeline import itol as _itol
            combined_tax_path = Path(outdir) / 'preload_combined_taxonomy.tsv'
            # Build combined taxonomy containing only dataset short IDs. For each dataset id
            # prefer taxonomy stored in DB; if missing, consult the distances table to find
            # the nearest reference and then look up taxonomy for that reference. If the
            # reference taxonomy is not present in the DB we will try to read it directly
            # from the provided taxa_tsv (without inserting the whole file into the DB).
            try:
                # determine dataset fasta (prefer mapped_fasta)
                if mapped_fasta:
                    ds_fasta_for_tax = Path(mapped_fasta)
                else:
                    ds_fasta_for_tax = Path(outdir) / f"preload_{getattr(args,'dataset','preload')}_seqs.fasta"

                # gather ids from dataset fasta if present
                ids = []
                try:
                    if ds_fasta_for_tax.exists():
                        from utils.fasta import read_fasta as _rf
                        ids = [h for h, _ in _rf(str(ds_fasta_for_tax))]
                except Exception:
                    ids = []

                # Prepare mappings by querying DB once for all ids to avoid per-id roundtrips
                db_tax = {}
                nearest_map = {}
                nearest_needed = set()
                with db.connect() as conn:
                    cur = conn.cursor()
                    if ids:
                        # fetch best nearest hits from distances (dataset-specific first, then any)
                        for iid in ids:
                            cur.execute("SELECT nearest FROM distances WHERE id = ? AND dataset = ? ORDER BY identity DESC LIMIT 1", (iid, getattr(args, 'dataset', 'preload')))
                            dr = cur.fetchone()
                            if not dr or not dr[0]:
                                cur.execute("SELECT nearest FROM distances WHERE id = ? ORDER BY identity DESC LIMIT 1", (iid,))
                                dr = cur.fetchone()
                            if dr and dr[0]:
                                nearest_map[iid] = dr[0]

                        # fetch taxonomy for dataset ids and any nearest references in one query
                        query_ids = list(set(ids) | set(nearest_map.values()))
                        if query_ids:
                            placeholders = ','.join('?' for _ in query_ids)
                            cur.execute(f"SELECT id, taxonomy, confidence FROM taxonomy WHERE id IN ({placeholders})", tuple(query_ids))
                            for rid, tax, conf in cur.fetchall():
                                db_tax[rid] = (tax, conf)

                        # track nearest references that still need file-based lookup
                        for n in set(nearest_map.values()):
                            if n not in db_tax:
                                nearest_needed.add(n)

                # If we need taxonomy for nearest references and a taxa_tsv was provided,
                # scan the taxa_tsv once and extract only the needed rows.
                file_tax_map = {}
                if nearest_needed and getattr(args, 'taxa', None):
                    taxa_file = getattr(args, 'taxa')
                    import gzip
                    open_fn = gzip.open if str(taxa_file).endswith('.gz') else open
                    needed = set(nearest_needed) | set(ids)
                    try:
                        with open_fn(taxa_file, 'rt') as tf:
                            header = tf.readline()
                            # determine whether header was read; if it's not a header and content
                            # looks like data we should process that line too
                            if 'Feature' in header or 'Taxon' in header:
                                pass
                            else:
                                parts = header.strip().split('\t')
                                if len(parts) >= 2:
                                    try:
                                        cid = db._canonical_from_header(parts[0])
                                    except Exception:
                                        cid = parts[0]
                                    if cid in needed:
                                        conf = float(parts[2]) if len(parts) > 2 else None
                                        file_tax_map[cid] = (parts[1], conf)
                            for line in tf:
                                parts = line.strip().split('\t')
                                if len(parts) < 2:
                                    continue
                                fid = parts[0]
                                try:
                                    cid = db._canonical_from_header(fid)
                                except Exception:
                                    cid = fid
                                if cid in needed:
                                    conf = float(parts[2]) if len(parts) > 2 else None
                                    file_tax_map[cid] = (parts[1], conf)
                    except Exception:
                        # if we can't read the taxa file, we simply won't have file-based fallbacks
                        file_tax_map = {}

                # write combined taxonomy file using DB taxonomy where present, otherwise
                # use nearest->DB taxonomy, otherwise nearest->taxa_tsv mapping, else NA
                # log some diagnostics about how many ids can be assigned
                assigned_count = 0
                if ids:
                    for iid in ids:
                        if iid in db_tax and db_tax[iid] and db_tax[iid][0]:
                            assigned_count += 1
                            continue
                        nearest = nearest_map.get(iid)
                        if nearest:
                            if nearest in db_tax and db_tax[nearest] and db_tax[nearest][0]:
                                assigned_count += 1
                                continue
                            if nearest in file_tax_map:
                                assigned_count += 1
                                continue
                logger.info("[PRELOAD] Combined taxonomy diagnostics: dataset_ids=%d, db_tax_rows_found=%d, nearest_refs=%d, file_tax_rows_found=%d, assigned=%d", len(ids), len(db_tax), len(nearest_map), len(file_tax_map), assigned_count)

                with open(combined_tax_path, 'w') as out_tax:
                    out_tax.write('ID\tTaxon\tConfidence\n')
                    if ids:
                        total = len(ids)
                        for idx, iid in enumerate(ids, start=1):
                            if idx == 1 or (idx % 500) == 0 or idx == total:
                                logger.info("[PRELOAD] Building combined taxonomy: processed %d/%d ids", idx, total)
                            tax = None
                            conf = None
                            if iid in db_tax and db_tax[iid] and db_tax[iid][0]:
                                tax, conf = db_tax[iid]
                            else:
                                nearest = nearest_map.get(iid)
                                if nearest:
                                    # try DB taxonomy for nearest
                                    if nearest in db_tax and db_tax[nearest] and db_tax[nearest][0]:
                                        tax, conf = db_tax[nearest]
                                    # then try taxa_tsv mapping we scanned
                                    elif nearest in file_tax_map:
                                        tax, conf = file_tax_map[nearest]
                            out_tax.write(f"{iid}\t{tax if tax is not None else 'NA'}\t{conf if conf is not None else 'NA'}\n")
                        logger.info("[PRELOAD] Finished writing combined taxonomy for %d ids to %s", total, combined_tax_path)
                    else:
                        logger.info("[PRELOAD] No dataset short IDs found in %s; writing header-only combined taxonomy for consistent palettes", ds_fasta_for_tax)

            except Exception:
                try:
                    with open(combined_tax_path, 'w') as out_tax:
                        out_tax.write('ID\tTaxon\tConfidence\n')
                except Exception:
                    pass

            # pass any explicit color CSV provided during preload as user overrides
            user_color_csv = getattr(args, 'colors', None)
            # Build a simple id_map from alias_entries (short, original) so the
            # generator can emit short IDs for the preload itol files.
            try:
                id_map = {}
                if alias_entries:
                    for short, orig in alias_entries:
                        id_map[orig] = short
                        id_map[short] = short
                _itol.generate_itol_colors(str(combined_tax_path), outdir, user_color_csv=user_color_csv, id_map=id_map)
            except TypeError:
                # fallback to previous call signature if the itol module doesn't
                # accept id_map for some reason
                _itol.generate_itol_colors(str(combined_tax_path), outdir, user_color_csv=user_color_csv)
            logger.info("[PRELOAD] Wrote taxa-level iTOL color files to %s", outdir)
        except Exception as e:
            logger.warning("[PRELOAD] Warning: failed to generate taxa-level iTOL files: %s", e)

        # optionally build a baseline tree/alignment using DB sequences
        if getattr(args, 'build_tree', False):
            if not getattr(args, 'ref', None):
                logger.info("[PRELOAD] --build-tree requested but --ref not provided; providing ref is recommended but not required")
            try:
                logger.info("[PRELOAD] Building baseline tree/alignment in %s from preloaded sequences", outdir)
                # tree.initialise_or_update_tree will try to use DB sequences when present
                user_fasta_for_tree = str(ds_fasta) if exported else args.fasta
                # By default do not include the entire DB as backbone when building a
                # baseline tree for a preload run (this can unintentionally pull in
                # thousands of reference sequences). If the caller explicitly
                # requested it via --use-db-backbone, pass `db` so DB sequences are
                # exported and used; otherwise pass None to build the backbone only
                # from the provided preloaded fasta.
                tree_db = db if getattr(args, 'use_db_backbone', False) else None
                tree_db_dataset = getattr(args, 'db_backbone_dataset', None)
                tree.initialise_or_update_tree(ref_fasta=getattr(args, 'ref', None), user_fasta=user_fasta_for_tree, outdir=outdir, db=tree_db, db_dataset=tree_db_dataset, threads=getattr(args, 'threads', None))
                logger.info("[PRELOAD] Baseline tree/alignment written to %s", outdir)
            except Exception as e:
                logger.warning("[PRELOAD] Warning: failed to build baseline tree: %s", e)


if __name__ == "__main__":
    main()
