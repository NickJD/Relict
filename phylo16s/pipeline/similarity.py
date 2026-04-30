from pathlib import Path
from utils.subprocess import run_cmd
from utils.fasta import read_fasta
import concurrent.futures
import os
import logging

logger = logging.getLogger(__name__)


def _safe_name(x: str) -> str:
    return ''.join(c if c.isalnum() or c in ('-', '_') else '_' for c in x)


def run_similarity(user_fasta: str, outdir: str, gg2_ref: str, db, threads=None):
    """Compute per-user-sequence similarity to gg2 (from DB) and to each preloaded dataset.

    - user_fasta: path to dereplicated user fasta
    - outdir: output directory for similarity files
    - gg2_ref: reference fasta path (not used directly here but kept for signature)
    - db: Database instance
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # collect user ids
    user_ids = [h.split()[0] for h, _ in read_fasta(user_fasta)]

    # start with gg2 distances as stored in DB (if any)
    gg2_map = {}
    try:
        with db.connect() as conn:
            cur = conn.cursor()
            cur.execute("SELECT id, nearest, identity FROM distances WHERE dataset = 'gg2'")
            for rid, nearest, identity in cur.fetchall():
                gg2_map[rid] = (nearest, identity)
    except Exception:
        gg2_map = {}

    similarity_long = []

    # add gg2 entries to similarity_long (ensure every user id has an entry, else NA)
    for uid in user_ids:
        if uid in gg2_map and gg2_map[uid][0] is not None:
            similarity_long.append((uid, 'gg2', gg2_map[uid][0], gg2_map[uid][1]))
        else:
            similarity_long.append((uid, 'gg2', 'NA', 'NA'))

    # iterate preloaded datasets and run vsearch in parallel
    datasets = db.get_datasets()
    jobs = []
    for ds in datasets:
        safe = _safe_name(ds)
        ds_fasta = outdir / f"dataset_{safe}.fasta"
        exported = db.export_dataset_fasta(ds, str(ds_fasta))
        if not exported:
            continue
        matches = outdir / f"matches_{safe}.tsv"
        thread_flag = f" --threads {int(threads)}" if threads and int(threads) > 0 else ""
        cmd = f"vsearch --usearch_global {user_fasta} --db {ds_fasta} --id 0.0 --blast6out {matches} --maxaccepts 1 --maxhits 1{thread_flag}"
        jobs.append((ds, safe, str(ds_fasta), str(matches), cmd))

    if jobs:
        max_workers = min(8, (os.cpu_count() or 1))
        logger.info("[SIM] Starting similarity searches for %d datasets (max_workers=%d)", len(jobs), max_workers)
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as exe:
            future_to_job = {exe.submit(run_cmd, job[4]): job for job in jobs}
            for fut in concurrent.futures.as_completed(future_to_job):
                job = future_to_job[fut]
                ds, safe, ds_fasta, matches, cmd = job
                try:
                    fut.result()
                    logger.info("[SIM] vsearch finished for dataset %s (matches=%s)", ds, matches)
                except Exception as e:
                    # vsearch failed for this dataset; record NA for all user ids and continue
                    for uid in user_ids:
                        similarity_long.append((uid, ds, 'NA', 'NA'))
                    print(f"[SIM] vsearch failed for dataset {ds}: {e}")
                    continue

                # parse matches file
                best_hits = {}
                try:
                    if Path(matches).exists():
                        with open(matches) as f:
                            for line in f:
                                parts = line.strip().split('\t')
                                if len(parts) < 3:
                                    continue
                                qid = parts[0]
                                sid = parts[1]
                                try:
                                    ident = float(parts[2])
                                except Exception:
                                    ident = None
                                best_hits[qid] = (sid, ident)
                except Exception:
                    best_hits = {}

                # append results for each user id
                for uid in user_ids:
                    if uid in best_hits:
                        sid, ident = best_hits[uid]
                        similarity_long.append((uid, ds, sid, f"{ident}" if ident is not None else 'NA'))
                    else:
                        similarity_long.append((uid, ds, 'NA', 'NA'))

                # write distances into DB for this dataset
                dist_entries = []
                for uid in user_ids:
                    if uid in best_hits:
                        sid, ident = best_hits[uid]
                        dist_entries.append((uid, ds, sid, ident))
                    else:
                        dist_entries.append((uid, ds, None, None))
                db.insert_distances(dist_entries)

    # write similarity_long.tsv
    long_path = outdir / 'similarity_long.tsv'
    with open(long_path, 'w') as out:
        out.write('id\tdataset\tbest_hit\tidentity\n')
        for row in similarity_long:
            out.write('\t'.join([str(x) if x is not None else 'NA' for x in row]) + '\n')

    # build pivot
    pivot = {}
    datasets_all = ['gg2'] + datasets
    for uid in user_ids:
        pivot[uid] = {}
    for uid, ds, best, ident in similarity_long:
        pivot[uid][ds] = (best, ident)

    pivot_path = outdir / 'similarity_pivot.tsv'
    with open(pivot_path, 'w') as out:
        headers = ['id']
        for ds in datasets_all:
            safe = _safe_name(ds)
            headers.extend([f"{safe}_best", f"{safe}_identity"])
        out.write('\t'.join(headers) + '\n')
        for uid in user_ids:
            row = [uid]
            for ds in datasets_all:
                val = pivot[uid].get(ds, ('NA', 'NA'))
                row.extend([str(val[0]), str(val[1])])
            out.write('\t'.join(row) + '\n')

    return str(long_path), str(pivot_path)

