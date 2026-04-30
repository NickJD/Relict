from pathlib import Path
from utils.subprocess import run_cmd
from utils.fasta import read_fasta, write_fasta
from collections import defaultdict
import re
import hashlib
import logging

logger = logging.getLogger(__name__)


def _make_unique_fasta(fasta_path, outdir):
    """If fasta has duplicate headers, write a new fasta with unique headers.

    Returns (path_to_fasta_to_use, mapping) where mapping maps unique_id -> original_id.
    If no duplicates, returns (fasta_path, {}).
    """
    records = list(read_fasta(fasta_path))
    seen = defaultdict(int)
    mapping = {}
    new_records = []
    duplicates = False

    import re

    def _sanitize_name(x: str) -> str:
        # take first token (drop whitespace metadata)
        token = str(x).split()[0]
        # prefer last pipe-field if present
        if '|' in token:
            token = token.split('|')[-1]
        # collapse runs of non-alnum to single underscore to produce Newick-safe names
        token = re.sub(r'[^0-9A-Za-z_.-]+', '_', token).strip('_')
        # fallback to a short hex if name becomes empty
        if not token:
            token = hashlib.md5(str(x).encode('utf-8')).hexdigest()[:8]
        return token

    for h, seq in records:
        base = _sanitize_name(h)
        if seen[base] == 0:
            new_h = base
        else:
            duplicates = True
            new_h = f"{base}__dup{seen[base]}"
        seen[base] += 1
        mapping[new_h] = h
        new_records.append((new_h, seq))

    if not duplicates:
        return fasta_path, {}

    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    unique_path = out / (Path(fasta_path).stem + "_unique.fasta")
    write_fasta(new_records, str(unique_path))

    # also write mapping file for debugging
    with open(out / 'id_map.tsv', 'w') as m:
        m.write('unique_id\toriginal_id\n')
        for new_id, orig in mapping.items():
            m.write(f"{new_id}\t{orig}\n")

    return str(unique_path), mapping


def initialise_or_update_tree(ref_fasta, user_fasta, outdir, db=None, db_dataset=None, threads=None):
    """Initialise or update a persistent tree/alignment in outdir.

    IMPORTANT: the backbone alignment and tree are built only from user
    sequences (and any sequences already stored in the DB if you later
    choose to provide them). The provided reference FASTA (e.g. gg2)
    is NOT used to build the alignment/tree — it is used only for
    classification/novelty steps elsewhere.

    Behaviour:
    - If an existing `current_alignment.fasta` is present it is used as
      the backbone.
    - Otherwise the backbone is created from the user sequences (the
      dereplicated input). The function will not align the full
      reference FASTA.
    - User sequences are added to the backbone via `mafft --addfragments`
      when the backbone already exists; if the backbone was created from
      the same user sequences no adding is needed.
    - A tree is built with FastTree from the final alignment.
    """
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    current_aln = out / "current_alignment.fasta"
    current_tree = out / "current_tree.nwk"

    # If no backbone alignment exists, try to create one from DB sequences
    # (if present), otherwise create from the user sequences. We do NOT
    # use the full reference fasta to create the backbone.
    if not current_aln.exists():
        base_aln = None

        # try to export DB sequences to FASTA and align them if db provided
        if db is not None:
            db_fasta = out / 'db_sequences.fasta'
            try:
                # allow exporting only a specific dataset from the DB to avoid
                # unintentionally pulling in the entire reference collection
                got = db.get_sequences_fasta(str(db_fasta), dataset=db_dataset)
            except Exception:
                got = False
            if got:
                logger.info("[TREE] Creating backbone alignment from DB sequences: %s", db_fasta)
                thread_flag = f" --thread {int(threads)}" if threads and int(threads) > 0 else ""
                cmd = f"mafft --auto{thread_flag} {db_fasta} > {out}/db_backbone_aln.fasta"
                logger.info("[TREE] Running mafft for DB backbone (threads=%s)", threads)
                run_cmd(cmd)
                logger.info("[TREE] mafft finished; db backbone written to %s/db_backbone_aln.fasta", out)
                base_aln = out / 'db_backbone_aln.fasta'

        if base_aln is None:
            logger.info("[TREE] No existing backbone alignment found; creating backbone from user sequences: %s", user_fasta)
            # align user sequences to create initial backbone
            thread_flag = f" --thread {int(threads)}" if threads and int(threads) > 0 else ""
            cmd = f"mafft --auto{thread_flag} {user_fasta} > {out}/user_backbone_aln.fasta"
            logger.info("[TREE] Running mafft for user backbone (threads=%s)", threads)
            run_cmd(cmd)
            logger.info("[TREE] mafft finished; user backbone written to %s/user_backbone_aln.fasta", out)
            base_aln = out / "user_backbone_aln.fasta"
    else:
        base_aln = current_aln

    # If there are user sequences, add them as fragments to the base alignment
    # If no user sequences (empty file), skip adding
    user_count = 0
    with open(user_fasta) as f:
        for line in f:
            if line.startswith(">"):
                user_count += 1

    if user_count == 0:
        # no new sequences; if no tree exists build from base_aln
        if not current_tree.exists():
            logger.info("[TREE] No user sequences and no tree present; building tree from backbone alignment")
            cmd = f"FastTree -nt {base_aln} > {current_tree}"
            try:
                logger.info("[TREE] Running FastTree to build tree from %s", base_aln)
                run_cmd(cmd)
                logger.info("[TREE] FastTree finished; tree written to %s", current_tree)
            except RuntimeError:
                logger.warning("[TREE] FastTree not available or failed to run. Skipping tree build.\n"
                               "Install FastTree (e.g. 'conda install -c bioconda fasttree') or provide its path.")
                # ensure current_aln exists
                if not current_aln.exists():
                    (out / "current_alignment.fasta").write_text(base_aln.read_text())
                return
        else:
            print("[TREE] No new user sequences and existing tree present; nothing to do")
        return

    # Add fragments (user sequences) to the backbone alignment
    logger.info("[TREE] Adding %d user sequences to backbone alignment", user_count)
    added_aln = out / "combined_aln.fasta"
    # If the backbone was just created from the user sequences, no need to
    # run --addfragments (the backbone already contains the user seqs).
    if str(base_aln).endswith('user_backbone_aln.fasta'):
        # backbone already contains the user sequences
        (out / "combined_aln.fasta").write_text(base_aln.read_text())
    else:
        # filter user_fasta to exclude sequences already present in base_aln or DB
        from utils.fasta import read_fasta, write_fasta
        base_ids = set()
        base_ids_norm = set()
        try:
            for h, _ in read_fasta(str(base_aln)):
                base_ids.add(h)
                base_ids_norm.add(_norm_id(h))
        except Exception:
            base_ids = set()
            base_ids_norm = set()

        db_ids = set()
        db_ids_norm = set()
        if db is not None:
            try:
                db_ids = db.get_all_ids()
                db_ids_norm = set(_norm_id(x) for x in db_ids)
            except Exception:
                db_ids = set()
                db_ids_norm = set()

        # write filtered user sequences
        filtered_user = out / 'user_new.fasta'
        new_records = []
        for h, s in read_fasta(user_fasta):
            nh = _norm_id(h)
            if nh in base_ids_norm or nh in db_ids_norm:
                continue
            new_records.append((h, s))

        if not new_records:
            logger.info('[TREE] No new sequences to add after filtering against DB/current alignment')
            # if no tree exists yet, build one from the existing backbone alignment
            if not current_tree.exists():
                logger.info('[TREE] No existing tree; building tree from backbone alignment')
                # ensure combined_aln points to base_aln so later code can build tree
                (out / 'combined_aln.fasta').write_text(base_aln.read_text())
                # ensure unique ids and run FastTree (reuse later logic)
                fasta_for_tree, id_map = _make_unique_fasta(str(out / 'combined_aln.fasta'), out)
                cmd = f"FastTree -nt {fasta_for_tree} > {current_tree}"
                try:
                    logger.info("[TREE] Running FastTree to build initial tree from backbone")
                    run_cmd(cmd)
                    logger.info("[TREE] FastTree initial tree written to %s", current_tree)
                except RuntimeError:
                    logger.warning("[TREE] FastTree not available or failed to run. Skipping tree build.\n"
                                   "Install FastTree (e.g. 'conda install -c bioconda fasttree') or provide its path.")
                    # still save backbone alignment as current_alignment
                    if not current_aln.exists():
                        (out / "current_alignment.fasta").write_text(base_aln.read_text())
                    return
                # remap ids if needed
                if id_map:
                    try:
                        newick = (out / 'current_tree.nwk').read_text()
                        for new_id, orig_id in id_map.items():
                            newick = newick.replace(new_id, orig_id)
                        (out / 'current_tree.nwk').write_text(newick)
                    except Exception as e:
                        print(f"[TREE] Warning: failed to remap IDs in tree: {e}")
                # save backbone as current_alignment for future runs
                if not current_aln.exists():
                    (out / "current_alignment.fasta").write_text(base_aln.read_text())
                print(f"[TREE] Updated tree written to {current_tree}")
                return
            else:
                print('[TREE] Existing tree present; nothing to do')
            return

        write_fasta(new_records, str(filtered_user))

        thread_flag = f" --thread {int(threads)}" if threads and int(threads) > 0 else ""
        cmd = f"mafft --addfragments{thread_flag} {filtered_user} {base_aln} > {added_aln}"
        logger.info("[TREE] Running mafft --addfragments (threads=%s) to add new user sequences", threads)
        run_cmd(cmd)
        logger.info("[TREE] mafft --addfragments finished; combined alignment at %s", added_aln)

    # Build tree from combined alignment
    logger.info("[TREE] Building tree from combined alignment (this may take a while)")
    # Ensure unique IDs in alignment (FastTree fails on duplicate names)
    fasta_for_tree, id_map = _make_unique_fasta(str(added_aln), out)

    cmd = f"FastTree -nt {fasta_for_tree} > {current_tree}"
    try:
        logger.info("[TREE] Running FastTree to build tree from combined alignment")
        run_cmd(cmd)
        logger.info("[TREE] FastTree finished; tree at %s", current_tree)
    except RuntimeError:
        logger.warning("[TREE] FastTree not available or failed to run. Skipping tree build.\n"
                       "Install FastTree (e.g. 'conda install -c bioconda fasttree') or provide its path.")
        # still save the combined alignment as current alignment so future runs have a backbone
        (out / "current_alignment.fasta").write_text(added_aln.read_text())
        return

    # If we used a modified fasta with unique ids, map names back in the tree
    if id_map:
        try:
            newick = (out / 'current_tree.nwk').read_text()
            for new_id, orig_id in id_map.items():
                # replace exact occurrences of new_id with orig_id
                newick = newick.replace(new_id, orig_id)
            (out / 'current_tree.nwk').write_text(newick)
        except Exception as e:
            print(f"[TREE] Warning: failed to remap IDs in tree: {e}")

    # Save the new alignment as current alignment for future runs
    (out / "current_alignment.fasta").write_text(added_aln.read_text())

    logger.info("[TREE] Updated tree written to %s", current_tree)


def _norm_id(x: str) -> str:
    """Normalize sequence/feature IDs for loose matching.

    Rules:
    - take first whitespace-separated token
    - if token contains pipes (|), take the last field
    - strip trailing underscore+digits (e.g. _1)
    """
    if x is None:
        return x
    x = str(x).split()[0]
    if '|' in x:
        x = x.split('|')[-1]
    x = re.sub(r'_[0-9]+$', '', x)
    return x
