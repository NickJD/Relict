"""Relict CLI — lightweight entrypoint for the relict package.

Implements the `preload`, `run`, and `regen-itol` commands for the src/ layout.
"""
import argparse
import logging
import os
import sys
from pathlib import Path

# When invoked directly (python src/relict/PhenGO-Predict.py) the package root (src)
# may not be on sys.path. Ensure the parent of this `relict` package is
# available so absolute imports like `relict.db.interface` work.
try:
    here = Path(__file__).resolve().parent
    src_root = str(here.parent)
    if src_root not in sys.path:
        sys.path.insert(0, src_root)
except Exception:
    pass

from relict.db.interface import Database
from relict.pipeline import classify, tree, itol, qc, derep, novelty
from relict.pipeline import cluster_report as _cluster_report
from relict.pipeline.classify import _derive_db_name as _classify_derive_db_name
from relict.pipeline.collapse import collapse_fasta_within_taxa
from relict.pipeline.workflow_helpers import (
    _assignment_source_is_fasta,
    build_orig_to_short_map as _build_orig_to_short_map_helper,
    build_placement_warning_rows,
    build_sequence_assessment_rows,
    classification_ids_matching_kingdom as _classification_ids_matching_kingdom_helper,
    collect_db_taxonomy_rows,
    iter_assignment_rows,
    load_classification_results_for_dataset,
    load_taxonomy_entries_from_assignments,
    merge_combined_taxonomy_rows,
    prune_dataset_by_kingdom as _prune_dataset_by_kingdom_helper,
    read_combined_taxonomy_ids,
    write_combined_taxonomy_tsv,
    write_placement_warning_tsv,
    write_sequence_assessment_tsv,
)
from relict.taxonomy import canonicalize_sequence_id, taxonomy_matches_kingdom
def _find_tree_file_in_dir(d: str):
    """Return path to a tree file in directory d if present, preferring current_tree.nwk."""
    p = Path(d)
    cand = p / 'current_tree.nwk'
    if cand.exists():
        return str(cand)
    # otherwise search for any .nwk or .tree file
    for ext in ('*.nwk', '*.tree', '*.tre'):
        found = next(p.glob(ext), None)
        if found:
            return str(found)
    return None


def _configure_logging(outdir: str):
    """Configure root logger: console INFO + optional file DEBUG in outdir.

    The root logger level must be set to DEBUG so that INFO/DEBUG records
    actually reach the handlers — handler-level filtering alone is not enough
    because Python's logging framework gates records at the logger level first.
    """
    logger = logging.getLogger()
    # remove existing handlers to avoid duplicate messages on repeated calls
    for h in list(logger.handlers):
        logger.removeHandler(h)
    # Set root logger to DEBUG so all records flow through to handlers;
    # each handler then applies its own level filter.
    logger.setLevel(logging.DEBUG)
    fmt = logging.Formatter("%(asctime)s %(levelname)s: %(message)s")
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(fmt)
    logger.addHandler(ch)
    try:
        fh = logging.FileHandler(os.path.join(outdir, 'relict.log'))
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(fmt)
        logger.addHandler(fh)
    except Exception:
        logger.warning("Could not write log file to %s", outdir)


def _build_orig_to_short(alias_entries):
    return _build_orig_to_short_map_helper(alias_entries, Database(':memory:'))


def _write_id_map_tsv(path: str | Path, entries, *, short_header: str = 'short_id', original_header: str = 'original_header'):
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, 'w') as handle:
        handle.write(f'{short_header}\t{original_header}\n')
        for short, original in entries:
            handle.write(f'{short}\t{original}\n')
    return str(p)


def _load_id_map_from_tsv(path: str | Path, db: Database | None = None):
    mapping = {}
    with open(path) as handle:
        next(handle, None)
        for line in handle:
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            short, orig = parts[0], parts[1]
            if not short or not orig:
                continue
            mapping[orig] = short
            mapping[short] = short
            if db is not None:
                try:
                    cid = db._canonical_from_header(orig)
                    if cid:
                        mapping[cid] = short
                except Exception:
                    pass
    return mapping


def _find_preferred_id_map(directory: str | Path):
    p = Path(directory)
    preferred = [
        p / 'preload_id_map.tsv',
        p / 'user_id_map.tsv',
        p / 'user_id_map.csv',
    ]
    for cand in preferred:
        if cand.exists():
            return cand
    try:
        for cand in p.glob('*_id_map.tsv'):
            if cand.name != 'id_map.tsv':
                return cand
    except Exception:
        pass
    return None


def _write_output_explanations(outdir: str):
    """Write short human-readable explanation files for common outputs under outdir.

    For each file discovered under outdir we create a companion file named
    <filename>.explain.txt containing a brief description of the file and how
    it is produced. This helps downstream users and CI to understand outputs.
    """
    p = Path(outdir)
    if not p.exists():
        return
    # mapping of filename substrings (lowercase) to explanation text
    patterns = [
        ('combined_taxonomy.tsv', 'Combined taxonomy TSV. Columns: ID\tTaxon\tConfidence. Used to generate iTOL color/legend files.'),
        ('preload_combined_taxonomy.tsv', 'Combined taxonomy TSV for preload dataset. Columns: ID\tTaxon\tConfidence.'),
        ('user_id_map.tsv', 'Mapping of short_id to original header produced when inserting user sequences into the DB.'),
        ('preload_id_map.tsv', 'Mapping of preload short IDs back to original FASTA headers. Use this to trace tree labels such as QUE06 back to source records.'),
        ('*_id_map.tsv', 'ID map mapping original headers to short (DB) ids. Useful for iTOL to display short ids.'),
        ('collapsed_map.tsv', 'Cluster map for collapsed sequences: rep_id\ttaxonomy\tcount.'),
        ('preload_collapsed_map.tsv', 'Cluster map for preload collapsed sequences: rep_id\ttaxonomy\tcount.'),
        ('collapsed_members.tsv', 'Member->representative mapping (member\trep) for collapsed clusters.'),
        ('preload_collapsed_members.tsv', 'Member->representative mapping for preload collapsed clusters.'),
        ('derep_short.fasta', 'Dereplicated FASTA where sequence headers are short ids assigned by the DB.'),
        ('derep_short_collapsed.fasta', 'Dereplicated FASTA after collapse; representatives for clusters kept with short ids.'),
        ('preload_short_collapsed.fasta', 'Collapsed preload FASTA; representatives retained for tree building.'),
        ('novelty_matches.tsv', 'vsearch BLAST-like output used to compute nearest-neighbour novelty identities.'),
        ('novelty_metrics.tsv', (
            'Per-sequence novelty metrics. Novelty and density are measured against YOUR previously '
            'submitted sequences (preload database + all prior run datasets), NOT against the global '
            'reference. Columns: ID, NearestIdentity, NearestHit, Novel, MatchesGE99, MatchesGE97, '
            'MatchesGE95, NoveltyScore, Crowding, SequencingPriority, DensitySource. '
            'DensitySource tells you which sequence pool was used for neighbourhood density.'
        )),
        ('sequence_assessment.tsv', (
            'Unified per-sequence assessment. '
            'COLUMN GROUPS: '
            '(1) TAXONOMY/CLASSIFICATION — Taxonomy, ClassificationHit, ClassificationIdentity, '
            'ClassificationConfidence: derived from the primary reference database (GTDB/SILVA). '
            'ClassificationHit is the reference accession vsearch matched. '
            'Repeated as Taxonomy_<DB>, ClassificationHit_<DB>, Identity_<DB>, Confidence_<DB> '
            'for each additional --alt-ref database. '
            '(2) NOVELTY — NearestHit, NearestIdentity, MatchesGE*, NoveltyScore, Crowding, '
            'SequencingPriority: derived from the PRELOAD sequences already in the DB. '
            'NearestHit is the closest sequence you have previously submitted. '
            'These columns are entirely independent of the reference database. '
            '(3) TREE/CLUSTER — InTree, ClusterRepresentative, ClusterSize, ClusteredMembers: '
            'records whether the sequence entered the phylogenetic tree directly or was '
            'represented by a cluster representative after --collapse.'
        )),
        ('taxonomy_input_warnings.tsv', 'Warnings about inconsistencies between the classifier reference FASTA and taxonomy TSV.'),
        ('tree_build_warnings.tsv', 'Warnings about weak phylogenetic signal, missing anchors, or poor alignment quality.'),
        ('tree_orientation_summary.tsv', 'Sequence-level audit of tree-input orientation checks. Reports which sequences were kept forward, reverse-complemented, or lacked orientation evidence before alignment.'),
        ('placement_warnings.tsv', 'Warnings about low-support placements, low identity matches, or potentially artefactual novelty assignments.'),
        ('taxa_assignments_classout.tsv', 'Synthetic classification-like TSV created when --taxa-assignments provided. Columns: id\tbest\tidentity\ttaxon\tconfidence.'),
        ('itol_dataset_membership.itol', 'iTOL DATASET_COLORSTRIP mapping sequence IDs to dataset colors (membership).'),
        ('itol_novelty.itol', 'iTOL colorstrip showing novelty (nearest identity) for run sequences.'),
        ('itol_dataset_preload.itol', 'iTOL colorstrip for the preload dataset; maps preload ids to the dataset color.'),
        ('.nwk', 'Newick tree file (phylogenetic tree). Commonly named current_tree.nwk.'),
        ('.itol', 'iTOL dataset file (text format) describing colors/strips/legends for visualization in iTOL.'),
        ('.log', 'Log file produced by the pipeline (relict.log) containing debug/info messages')
    ]

    rows = []
    # iterate files in outdir (non-recursive and recursive) and collect explanations
    for fp in sorted(p.rglob('*')):
        if fp.is_dir():
            continue
        if fp.name.endswith('.explain.txt'):
            try:
                fp.unlink()
            except Exception:
                pass
            continue
        name = fp.name
        lname = name.lower()
        expl = None
        for patt, text in patterns:
            # treat wildcard at start/end
            if patt.startswith('*') and patt.endswith('*'):
                key = patt.strip('*').lower()
                if key in lname:
                    expl = text
                    break
            elif patt.startswith('*'):
                key = patt.lstrip('*').lower()
                if lname.endswith(key):
                    expl = text
                    break
            elif patt.endswith('*'):
                key = patt.rstrip('*').lower()
                if lname.startswith(key):
                    expl = text
                    break
            elif patt.startswith('.'):
                if lname.endswith(patt):
                    expl = text
                    break
            else:
                if patt in lname:
                    expl = text
                    break
        if expl is None:
            # generic descriptions based on extension
            if lname.endswith('.fasta') or lname.endswith('.fa'):
                expl = 'FASTA file containing sequences. May be dereplicated, collapsed, or exported from the DB.'
            elif lname.endswith('.tsv') or lname.endswith('.csv'):
                expl = 'Tab/CSV-delimited table used by the pipeline.'
            elif lname.endswith('.uc'):
                expl = 'vsearch UC membership file describing clusters.'
            else:
                expl = 'Pipeline output file.'

        rows.append((fp.name, str(fp), expl))

    try:
        manifest = p / 'OUTPUT_EXPLANATIONS.tsv'
        with open(manifest, 'w') as ef:
            ef.write('File\tPath\tDescription\n')
            for name, path, expl in rows:
                ef.write(f"{name}\t{path}\t{expl}\n")
    except Exception:
        pass


def _classification_ids_matching_kingdom(classification_tsv: str, kingdom: str):
    return _classification_ids_matching_kingdom_helper(classification_tsv, kingdom)


def _prune_dataset_by_kingdom(db: Database, dataset: str, kingdom: str | None, log_prefix: str):
    deleted = _prune_dataset_by_kingdom_helper(db, dataset, kingdom)
    if deleted:
        logging.getLogger(__name__).info(
            "%s Removed %d sequences from dataset %s not matching kingdom %s",
            log_prefix,
            deleted,
            dataset,
            kingdom,
        )
    return deleted


def _same_path(path_a: str | None, path_b: str | None) -> bool:
    if not path_a or not path_b:
        return False
    try:
        return os.path.samefile(path_a, path_b)
    except Exception:
        try:
            return Path(path_a).resolve() == Path(path_b).resolve()
        except Exception:
            return str(path_a) == str(path_b)


def _build_alt_databases(args) -> list:
    """Return a list of (ref_fasta, taxa_tsv_or_None, db_name) for alt references."""
    alt_refs = getattr(args, 'alt_ref', None) or []
    alt_taxa = getattr(args, 'alt_taxa', None) or []
    alt_names = getattr(args, 'alt_ref_name', None) or []
    result = []
    for i, ref in enumerate(alt_refs):
        taxa = alt_taxa[i] if i < len(alt_taxa) else None
        raw_name = alt_names[i] if i < len(alt_names) else None
        name = raw_name if raw_name else _classify_derive_db_name(ref)
        result.append((ref, taxa, name))
    return result


def _store_alt_taxonomy_in_db(db: Database, all_results: dict, main_db_name: str):
    """Persist alt-db classification results to taxonomy_alt table.

    Iterates over *all_results* ``{db_name: {qid: (hit, pct, tax, conf)}}``,
    skipping the primary database (its results go into the standard taxonomy table).
    """
    log = logging.getLogger(__name__)
    for db_name, results in all_results.items():
        if db_name == main_db_name:
            continue
        alt_entries = []
        for qid, row in results.items():
            if not isinstance(row, (tuple, list)) or len(row) < 4:
                continue
            hit, pct, tax, conf = row[:4]
            if tax == 'NA' and hit == 'NA':
                continue  # skip unclassified rows to keep table lean
            alt_entries.append((qid, db_name, tax, conf, hit, pct))
        if alt_entries:
            try:
                db.insert_taxonomy_alt(alt_entries)
                log.info("[DB] Stored %d alt-db taxonomy entries for ref_db=%s", len(alt_entries), db_name)
            except Exception as e:
                log.warning("[DB] Failed to store alt-db taxonomy for %s: %s", db_name, e)


def _resolve_reference_inputs(
    ref_fasta: str | None,
    taxa_tsv: str | None,
    taxa_assignments: str | None,
    *,
    source_fasta_path: str | None,
    log_prefix: str,
):
    """Resolve whether --taxa-assignments is a TSV assignment file or a reference FASTA.

    External FASTA/FASTA.gz files provided via --taxa-assignments are treated as
    the effective classifier reference database. TSV inputs (and the special case
    where the assignments file is the same file as the source FASTA) keep the
    legacy direct-assignment behaviour.
    """
    effective_ref = ref_fasta
    effective_taxa = taxa_tsv
    assignments_tsv = None

    if not taxa_assignments:
        return effective_ref, effective_taxa, assignments_tsv

    use_as_reference = False
    try:
        use_as_reference = (
            not _same_path(taxa_assignments, source_fasta_path)
            and _assignment_source_is_fasta(taxa_assignments, source_fasta_path=source_fasta_path)
        )
    except Exception:
        use_as_reference = False

    if use_as_reference:
        log = logging.getLogger(__name__)
        if effective_ref and not _same_path(effective_ref, taxa_assignments):
            log.info(
                "%s Treating --taxa-assignments=%s as the GTDB/reference FASTA and using it instead of --ref=%s (supported for compatibility, but --ref is the preferred flag for reference databases)",
                log_prefix,
                taxa_assignments,
                effective_ref,
            )
        else:
            log.info(
                "%s Treating --taxa-assignments=%s as the GTDB/reference FASTA for classification (supported for compatibility, but --ref is the preferred flag for reference databases)",
                log_prefix,
                taxa_assignments,
            )
        if not effective_taxa:
            log.info(
                "%s No --taxa TSV provided; taxonomy will be parsed directly from reference FASTA headers",
                log_prefix,
            )
        effective_ref = taxa_assignments
        return effective_ref, effective_taxa, None

    assignments_tsv = taxa_assignments
    return effective_ref, effective_taxa, assignments_tsv


def cmd_preload(args):
    db = Database(args.db)
    db.initialise()
    outdir = args.out or '.'
    threads = int(getattr(args, 'threads', 4) or 4)
    requested_kingdom = getattr(args, 'kingdom', None)
    requested_kingdom = str(requested_kingdom) if requested_kingdom else None
    os.makedirs(outdir, exist_ok=True)
    _configure_logging(outdir)
    logging.getLogger(__name__).info("[PRELOAD] Starting preload into %s", args.db)

    alias_entries, mapped_fasta = db.preload_from_files(
        args.fasta,
        taxa_tsv=getattr(args, 'taxa', None),
        color_csv=getattr(args, 'colors', None),
        source='preload',
        dataset=getattr(args, 'dataset', 'preload'),
        outdir=outdir,
        shorten_ids=bool(getattr(args, 'shorten_ids', True)),
    )
    try:
        if alias_entries:
            preload_map_path = _write_id_map_tsv(Path(outdir) / 'preload_id_map.tsv', alias_entries)
            logging.getLogger(__name__).info('[PRELOAD] Wrote preload id mapping to %s', preload_map_path)
    except Exception as e:
        logging.getLogger(__name__).warning('[PRELOAD] Could not write preload id mapping file: %s', e)
    effective_ref, effective_taxa_tsv, assignment_tsv = _resolve_reference_inputs(
        getattr(args, 'ref', None),
        getattr(args, 'taxa', None),
        getattr(args, 'taxa_assignments', None),
        source_fasta_path=args.fasta,
        log_prefix='[PRELOAD]',
    )
    classification_requested = bool(getattr(args, 'classify', False) or (getattr(args, 'taxa_assignments', None) and not assignment_tsv))

    # If the user provided a TSV of predetermined taxa assignments, use it
    # instead of running the classifier. The TSV should be tab-separated with
    # at least two columns: sequence_header<TAB>taxonomy[<TAB>confidence]
    if assignment_tsv:
        taxa_file = assignment_tsv
        logging.getLogger(__name__).info("[PRELOAD] Using taxa assignments from %s (skipping classifier)", taxa_file)
        # build mapping orig->short from alias_entries
        orig_to_short = _build_orig_to_short(alias_entries)
        try:
            tax_entries = load_taxonomy_entries_from_assignments(
                taxa_file,
                orig_to_short,
                db,
                getattr(args, 'dataset', 'preload'),
                source_fasta_path=args.fasta,
            )
        except Exception as e:
            logging.getLogger(__name__).warning("[PRELOAD] Failed to read taxa assignments file %s: %s", taxa_file, e)
            tax_entries = []

        if tax_entries:
            db.insert_taxonomy(tax_entries)
            logging.getLogger(__name__).info("[PRELOAD] Inserted/updated taxonomy for %d preloaded ids from taxa_assignments", len(tax_entries))

    # If classification requested, run classifier on the mapped fasta (short ids)
    # unless the user supplied a taxa assignments TSV, in which case use that
    # instead and skip running the external classifier.
    if classification_requested and not assignment_tsv:
        if not effective_ref:
            logging.getLogger(__name__).info("[PRELOAD] Classification requested but no reference FASTA was provided via --ref or --taxa-assignments; skipping classification")
        else:
            input_for_classify = str(mapped_fasta) if mapped_fasta else args.fasta
            logging.getLogger(__name__).info("[PRELOAD] Classifying preloaded fasta %s against %s", input_for_classify, effective_ref)

            alt_databases = _build_alt_databases(args)
            ref_name = getattr(args, 'ref_name', None) or _classify_derive_db_name(effective_ref)
            main_db = getattr(args, 'main_ref', None) or ref_name
            all_results: dict = {}

            if alt_databases:
                logging.getLogger(__name__).info(
                    "[PRELOAD] Multi-database classification: primary=%s, alt=%s, main=%s",
                    ref_name, [n for _, _, n in alt_databases], main_db,
                )
                class_out, all_results = classify.run_all_classifications(
                    input_for_classify, outdir,
                    primary_ref=effective_ref,
                    primary_taxa=effective_taxa_tsv,
                    primary_name=ref_name,
                    alt_refs=alt_databases,
                    threads=threads,
                    main_db=main_db,
                )
                _store_alt_taxonomy_in_db(db, all_results, main_db)
            else:
                class_out = classify.run_classification(input_for_classify, outdir, ref_fasta=effective_ref, taxa_tsv=effective_taxa_tsv, threads=threads)

            # parse classification output and persist taxonomy/distances only for
            # the preloaded sequence ids (short ids present in the DB)
            orig_to_short = _build_orig_to_short(alias_entries)
            try:
                tax_entries, dist_entries = load_classification_results_for_dataset(
                    class_out,
                    orig_to_short,
                    db,
                    getattr(args, 'dataset', 'preload'),
                )

                if tax_entries:
                    db.insert_taxonomy(tax_entries)
                    logging.getLogger(__name__).info("[PRELOAD] Inserted/updated taxonomy for %d preloaded ids", len(tax_entries))
                if dist_entries:
                    db.insert_distances(dist_entries)
                    logging.getLogger(__name__).info("[PRELOAD] Inserted/updated distances for %d preloaded ids", len(dist_entries))
            except Exception as e:
                logging.getLogger(__name__).warning("[PRELOAD] Failed to parse classification output: %s", e)

    # Optionally build a baseline tree/alignment from the preloaded sequences
    if getattr(args, 'build_tree', False):
        try:
            # prefer mapped fasta (short ids) if present
            user_fasta = str(mapped_fasta) if mapped_fasta else args.fasta

            # optionally collapse preloaded sequences before building tree
            if getattr(args, 'collapse', False):
                # require classification to be run for safe taxon-based collapsing
                if not (classification_requested or assignment_tsv):
                    logging.getLogger(__name__).warning("[PRELOAD COLLAPSE] --collapse requires --classify to be set for safe taxon grouping; skipping collapse")
                else:
                    try:
                        threshold = float(getattr(args, 'collapse_threshold', 99.9))
                    except Exception:
                        threshold = 99.9

                    # build qid -> tax mapping from DB taxonomy for this dataset
                    qid_to_tax = {}
                    try:
                        with db.connect() as conn:
                            cur = conn.cursor()
                            cur.execute("SELECT s.id, t.taxonomy FROM sequences s LEFT JOIN taxonomy t ON s.id = t.id WHERE s.dataset = ?", (getattr(args, 'dataset', 'preload'),))
                            for rid, tax in cur.fetchall():
                                qid_to_tax[rid] = tax
                    except Exception:
                        qid_to_tax = {}

                    # group records by tax and cluster within groups
                    from relict.utils.fasta import read_fasta
                    taxa_groups = {}
                    try:
                        for h, s in read_fasta(user_fasta):
                            tax = qid_to_tax.get(h)
                            taxa_groups.setdefault(tax, []).append((h, s))
                    except Exception:
                        taxa_groups = {None: []}
                    try:
                        artifacts = collapse_fasta_within_taxa(
                            taxa_groups,
                            outdir,
                            'preload_short_collapsed.fasta',
                            'preload_collapsed_map.tsv',
                            'preload_collapsed_members.tsv',
                            threshold=threshold,
                            threads=threads,
                            log_prefix='[PRELOAD COLLAPSE]',
                        )
                        user_fasta = artifacts.collapsed_path
                        logging.getLogger(__name__).info("[PRELOAD COLLAPSE] Wrote collapsed preload fasta %s (reps=%d)", user_fasta, len(artifacts.collapsed_records))
                    except Exception as e:
                        logging.getLogger(__name__).warning("[PRELOAD COLLAPSE] Failed to write collapsed preload fasta: %s", e)

            logging.getLogger(__name__).info("[PRELOAD] Building baseline tree/alignment in %s from preloaded sequences", outdir)
            tree.initialise_or_update_tree(
                ref_fasta=effective_ref,
                user_fasta=user_fasta,
                outdir=outdir,
                db=None,
                threads=threads,
                anchor_file=getattr(args, 'anchors', None),
                tree_method=getattr(args, 'tree_method', 'fasttree'),
            )
            logging.getLogger(__name__).info("[PRELOAD] Baseline tree/alignment written to %s", outdir)
            try:
                warning_rows = tree.collect_tree_build_warnings(user_fasta=str(user_fasta), anchor_file=getattr(args, 'anchors', None), db=None)
                warning_rows.extend(tree.summarize_alignment_quality(str(Path(outdir) / 'current_alignment.fasta')))
                if warning_rows:
                    warn_path = tree.write_tree_warning_tsv(outdir, warning_rows)
                    logging.getLogger(__name__).warning("[PRELOAD] Tree/alignment warnings written to %s", warn_path)
            except Exception as e:
                logging.getLogger(__name__).warning("[PRELOAD] Failed to summarise tree/alignment quality: %s", e)
        except Exception as e:
            logging.getLogger(__name__).warning("[PRELOAD] Failed to build baseline tree: %s", e)

    # If kingdom filter requested for preload, remove any preloaded sequences
    # whose assigned taxonomy indicates they are not the requested kingdom.
    try:
        _prune_dataset_by_kingdom(
            db,
            getattr(args, 'dataset', 'preload'),
            requested_kingdom,
            '[PRELOAD]',
        )
    except Exception as e:
        logging.getLogger(__name__).warning("[PRELOAD] Kingdom-based pruning failed: %s", e)

    # Build combined taxonomy and generate iTOL color files for the preload dataset.
    try:
        out_p = Path(outdir)
        combined_tax = out_p / 'preload_combined_taxonomy.tsv'
        if not combined_tax.exists():
            try:
                with db.connect() as conn:
                    cur = conn.cursor()
                    cur.execute("SELECT s.id, t.taxonomy, t.confidence FROM sequences s LEFT JOIN taxonomy t ON s.id = t.id WHERE s.dataset = ?", (getattr(args, 'dataset', 'preload'),))
                    rows = cur.fetchall()
                write_combined_taxonomy_tsv(combined_tax, rows)
                logging.getLogger(__name__).info("[PRELOAD] Wrote combined taxonomy for %d ids to %s", len(rows), combined_tax)
            except Exception as e:
                logging.getLogger(__name__).warning("[PRELOAD] Failed to build combined taxonomy: %s", e)

        # Build id_map from alias_entries so the generator can emit short ids
        id_map = {}
        try:
            if alias_entries:
                for short, orig in alias_entries:
                    id_map[orig] = short
                    id_map[short] = short
                    try:
                        cid = db._canonical_from_header(orig)
                        if cid:
                            id_map[cid] = short
                    except Exception:
                        pass
        except Exception:
            id_map = {}

        # Call generator (fallback if id_map not accepted)
        try:
            tree_path = Path(outdir) / 'current_tree.nwk'
            tfile = str(tree_path) if tree_path.exists() else _find_tree_file_in_dir(outdir)
            itol.generate_itol_colors(str(combined_tax), outdir, user_color_csv=getattr(args, 'colors', None), id_map=id_map, tree_file=tfile, phylum_groups=getattr(args, 'group_phyla', None))
            logging.getLogger(__name__).info("[PRELOAD] Generated iTOL color files in %s", outdir)
        except TypeError:
            try:
                itol.generate_itol_colors(str(combined_tax), outdir, user_color_csv=getattr(args, 'colors', None), phylum_groups=getattr(args, 'group_phyla', None))
                logging.getLogger(__name__).info("[PRELOAD] Generated iTOL color files in %s (no id_map)", outdir)
            except Exception as e:
                logging.getLogger(__name__).warning("[PRELOAD] Failed to generate iTOL colors: %s", e)
        except Exception as e:
            logging.getLogger(__name__).warning("[PRELOAD] Failed to generate iTOL colors: %s", e)

        # ── Functional annotations (optional) ─────────────────────────────────
        try:
            func_tsv = getattr(args, 'functional', None)
            if func_tsv:
                try:
                    written = itol.write_functional_annotations(str(func_tsv), outdir, id_map=id_map)
                    logging.getLogger(__name__).info("[PRELOAD] Wrote functional annotation iTOL files: %s", ','.join(written) if written else '(none)')
                except Exception as e:
                    logging.getLogger(__name__).warning("[PRELOAD] Functional annotations generation failed: %s", e)
        except Exception:
            pass

        # ── Draft rumen functional groups (optional) ──────────────────────────
        if getattr(args, 'draft_rumen_functions', False):
            try:
                combined_tax_path = str(out_p / 'combined_taxonomy.tsv')
                tsv_out, itol_out = itol.generate_rumen_function_draft(
                    combined_tax_path, outdir, id_map=id_map
                )
                if tsv_out:
                    logging.getLogger(__name__).info(
                        "[PRELOAD] Draft rumen functional annotation: %s", tsv_out
                    )
                if itol_out:
                    logging.getLogger(__name__).info(
                        "[PRELOAD] Rumen functional iTOL file: %s", itol_out
                    )
            except Exception as e:
                logging.getLogger(__name__).warning(
                    "[PRELOAD] Draft rumen functions generation failed: %s", e
                )

        # Write a simple dataset colorstrip mapping preloaded ids to dataset color
        try:
            ds_color = getattr(args, 'dataset_color', None) if hasattr(args, 'dataset_color') else None
            if not ds_color:
                # deterministic dataset-level color (use distinct palette)
                ds_color = itol._name_to_dataset_color(getattr(args, 'dataset', 'preload'))
            itol_path = out_p / 'itol_dataset_preload.itol'
            dataset_label = getattr(args, 'dataset', 'preload')
            with db.connect() as conn:
                cur = conn.cursor()
                cur.execute("SELECT id FROM sequences WHERE dataset = ?", (getattr(args, 'dataset', 'preload'),))
                id_to_color = {iid: ds_color for (iid,) in cur.fetchall()}
            itol.write_dataset_colorstrip(str(itol_path), dataset_label, id_to_color, legend_title=f"{dataset_label} legend")
            logging.getLogger(__name__).info("[PRELOAD] Wrote dataset ITOL colorstrip to %s", itol_path)
        except Exception as e:
            logging.getLogger(__name__).warning("[PRELOAD] Failed to write dataset iTOL colorstrip: %s", e)
    except Exception:
        pass

    # write brief explanations for files produced by preload
    try:
        _write_output_explanations(outdir)
    except Exception:
        pass


def cmd_run(args):
    db = Database(args.db)
    db.initialise()
    outdir = args.out
    threads = int(getattr(args, 'threads', 4) or 4)
    os.makedirs(outdir, exist_ok=True)
    _configure_logging(outdir)
    logging.getLogger(__name__).info("[RUN] Starting run pipeline (input=%s)", args.input)
    effective_ref, effective_taxa_tsv, assignment_tsv = _resolve_reference_inputs(
        getattr(args, 'ref', None),
        getattr(args, 'taxa', None),
        getattr(args, 'taxa_assignments', None),
        source_fasta_path=args.input,
        log_prefix='[RUN]',
    )
    if not effective_ref:
        raise SystemExit("[RUN] A reference FASTA is required via --ref, or `--taxa-assignments` must point to a GTDB/reference FASTA rather than a TSV assignments file.")

    # QC
    qc_out = qc.run_qc(args.input, outdir, min_len=getattr(args, 'min_len', 1200), max_n=getattr(args, 'max_n', 5))

    # derep
    derep_out = derep.run_derep(qc_out, outdir)

    # map user-provided dereplicated IDs to short IDs and insert into DB
    from relict.utils.fasta import read_fasta, write_fasta
    mapped_derep = Path(outdir) / 'derep_short.fasta'
    used_ids = set(db.get_all_ids())
    orig_to_short = {}
    mapped_records = []
    skipped_existing = 0

    # If the user requested a kingdom filter run classification early on the
    # dereplicated fasta so we can keep only sequences assigned to the chosen
    # kingdom. This avoids inserting unwanted sequences into the DB.
    early_class_out = None
    allowed_qids = None
    kingdom = getattr(args, 'kingdom', None)
    if kingdom:
        kingdom_text = str(kingdom)
        if not effective_ref:
            logging.getLogger(__name__).warning("[RUN] --kingdom specified but no reference FASTA was available; cannot classify to filter; proceeding without kingdom filtering")
            allowed_qids = None
        else:
            try:
                logging.getLogger(__name__).info("[RUN] Running pre-insert classification on dereplicated fasta to filter by kingdom=%s", kingdom)
                early_class_out = classify.run_classification(str(derep_out), outdir, ref_fasta=effective_ref, taxa_tsv=effective_taxa_tsv, threads=threads)
                allowed_qids = _classification_ids_matching_kingdom(early_class_out, kingdom_text)
                logging.getLogger(__name__).info("[RUN] Kingdom filter: %d dereplicated sequences match %s", len(allowed_qids), kingdom)
            except Exception as e:
                logging.getLogger(__name__).warning("[RUN] Failed to run pre-insert classification for kingdom filtering: %s", e)
                allowed_qids = None

    for h, s in read_fasta(derep_out):
        # if kingdom filtering is active, skip sequences that did not match
        if allowed_qids is not None:
            # try a few candidate header forms for matching
            candidates = [h, h.split('|')[-1], h.split()[-1] if h.split() else h]
            hit = False
            for c in candidates:
                if c in allowed_qids:
                    hit = True
                    break
            if not hit:
                continue

        # If DB already contains this canonical sequence, we still want to include it
        # in the tree (for visualization), but we'll use its existing short ID.
        # The "INSERT OR IGNORE" during insert_sequences will prevent DB duplicates.
        try:
            cid = db._canonical_from_header(h)
        except Exception:
            cid = None
        if cid and cid in used_ids:
            skipped_existing += 1
            # Use the existing canonical ID as the short ID
            short = cid
            orig_to_short[h] = short
            orig_to_short[cid] = short
            # Still add to mapped_records so it appears in the tree
            mapped_records.append((short, s))
        else:
            try:
                short = db.choose_effective_sequence_id(
                    h,
                    used_ids,
                    shorten_ids=bool(getattr(args, 'shorten_ids', True)),
                )
            except ValueError as e:
                raise SystemExit(f"[RUN] {e}")
            mapped_records.append((short, s))
            orig_to_short[h] = short
            try:
                cid = canonicalize_sequence_id(h)
                if cid:
                    orig_to_short[cid] = short
            except Exception:
                pass
    if skipped_existing:
        logging.getLogger(__name__).info("[DB] Found %d sequences already in DB (by canonical id); keeping them for tree inclusion", skipped_existing)
    write_fasta(mapped_records, str(mapped_derep))
    logging.getLogger(__name__).info("[DB] Mapped %d user sequence IDs to short IDs and wrote %s", len(mapped_records), mapped_derep)
    if not mapped_records:
        logging.getLogger(__name__).warning("[DB] No sequences were mapped — check if input sequences were filtered")
    # insert mapped records into DB (dataset provided by user)
    run_dataset = getattr(args, 'dataset', 'user')
    db.insert_sequences(mapped_records, dataset=run_dataset)
    # write mapping file
    try:
        map_path = Path(outdir) / 'user_id_map.tsv'
        _write_id_map_tsv(map_path, ((short, orig) for orig, short in orig_to_short.items()))
        logging.getLogger(__name__).info("[DB] Wrote user id mapping to %s", map_path)
    except Exception as e:
        logging.getLogger(__name__).warning("[DB] Could not write user id mapping file: %s", e)

    for short, _seq in mapped_records:
        orig_to_short[short] = short

    # classify (or use external taxa assignments if provided)
    class_out = None
    all_class_results: dict = {}
    run_main_db_name: str = 'main'
    if assignment_tsv:
        taxa_file = assignment_tsv
        logging.getLogger(__name__).info("[RUN] Using taxa assignments from %s (skipping classifier)", taxa_file)
        # parse provided TSV and insert taxonomy for mapped short ids
        try:
            tax_entries_local = load_taxonomy_entries_from_assignments(
                taxa_file,
                orig_to_short,
                db,
                run_dataset,
                source_fasta_path=args.input,
            )
        except Exception as e:
            logging.getLogger(__name__).warning("[RUN] Failed to read taxa assignments file %s: %s", taxa_file, e)
            tax_entries_local = []

        # persist taxonomy for mapped ids
        if tax_entries_local:
            try:
                db.insert_taxonomy(tax_entries_local)
                logging.getLogger(__name__).info("[DB] Inserted/updated taxonomy for %d ids from taxa_assignments", len(tax_entries_local))
            except Exception as e:
                logging.getLogger(__name__).warning("[DB] Failed to insert taxonomy from taxa_assignments: %s", e)

        # create a synthetic classification-like file so downstream code that
        # expects `class_out` can run unchanged. Format: qid\tbest\tidentity\ttaxon\tconfidence
        try:
            class_out_path = Path(outdir) / 'taxa_assignments_classout.tsv'
            with open(class_out_path, 'w') as cf:
                cf.write('id\tbest\tidentity\ttaxon\tconfidence\n')
                try:
                    for row in iter_assignment_rows(taxa_file, source_fasta_path=args.input):
                        qid = row.get('qid', 'NA')
                        tax = row.get('tax') if row.get('tax') is not None else 'NA'
                        conf = row.get('confidence') if row.get('confidence') is not None else 'NA'
                        cf.write(f"{qid}\t\tNA\t{tax}\t{conf}\n")
                except Exception:
                    pass
            class_out = str(class_out_path)
        except Exception:
            class_out = None
    else:
        # Only run classification here if we did not already run an early
        # pre-insert classification for kingdom filtering and the user did
        # not supply taxa_assignments.
        if 'early_class_out' in locals() and early_class_out:
            class_out = early_class_out
        else:
            alt_databases = _build_alt_databases(args)
            ref_name = getattr(args, 'ref_name', None) or _classify_derive_db_name(effective_ref)
            run_main_db_name = getattr(args, 'main_ref', None) or ref_name

            if alt_databases:
                logging.getLogger(__name__).info(
                    "[RUN] Multi-database classification: primary=%s, alt=%s, main=%s",
                    ref_name, [n for _, _, n in alt_databases], run_main_db_name,
                )
                class_out, all_class_results = classify.run_all_classifications(
                    str(mapped_derep), outdir,
                    primary_ref=effective_ref,
                    primary_taxa=effective_taxa_tsv,
                    primary_name=ref_name,
                    alt_refs=alt_databases,
                    threads=threads,
                    main_db=run_main_db_name,
                )
                _store_alt_taxonomy_in_db(db, all_class_results, run_main_db_name)
            else:
                class_out = classify.run_classification(str(mapped_derep), outdir, ref_fasta=effective_ref, taxa_tsv=effective_taxa_tsv, threads=threads)

    # novelty
    target_fasta = getattr(args, 'target', None)
    novelty_out = novelty.run_novelty(str(mapped_derep), effective_ref, outdir, db=db, run_dataset=run_dataset, threads=threads, target_fasta=target_fasta)
    try:
        novelty_metrics_out = novelty.build_reference_novelty_metrics(
            str(mapped_derep),
            effective_ref,
            novelty_out,
            outdir,
            threads=threads,
            db=db,
            run_dataset=run_dataset,
            target_fasta=target_fasta,
        )
        logging.getLogger(__name__).info("[NOVELTY] Wrote novelty metrics to %s", novelty_metrics_out)
    except Exception as e:
        logging.getLogger(__name__).warning("[NOVELTY] Failed to build novelty metrics: %s", e)

    # persist taxonomy & distances for user sequences only
    try:
        tax_entries, dist_entries = load_classification_results_for_dataset(
            class_out or '',
            orig_to_short,
            db,
            run_dataset,
        )
    except Exception as e:
        logging.getLogger(__name__).warning("[RUN] Failed to parse classification output: %s", e)
        tax_entries, dist_entries = [], []

    if tax_entries:
        db.insert_taxonomy(tax_entries)
        logging.getLogger(__name__).info("[DB] Inserted/updated taxonomy for %d ids", len(tax_entries))
    if dist_entries:
        db.insert_distances(dist_entries)
        logging.getLogger(__name__).info("[DB] Inserted/updated distances for %d ids", len(dist_entries))
    # Safety-net: remove any run sequences with explicit non-matching kingdom assignments.
    try:
        _prune_dataset_by_kingdom(db, run_dataset, str(kingdom) if kingdom else None, '[RUN]')
    except Exception as e:
        logging.getLogger(__name__).warning("[RUN] Kingdom-based pruning failed: %s", e)

    # Initialise cluster-tracking variables used by both the tree section and
    # the assessment section below.  They will be populated if --collapse is
    # active; otherwise they stay empty and the assessment omits cluster columns.
    tree_fasta = mapped_derep
    run_member_to_rep: dict = {}
    run_rep_to_members: dict = {}

    # update tree/alignment
    try:
        # report DB sequence count for diagnostics so users can verify that
        # preloaded sequences exist and will be used to seed the backbone
        try:
            db_ids = db.get_all_ids()
            logging.getLogger(__name__).info("[TREE] DB contains %d sequences; these will be considered for backbone construction", len(db_ids))
        except Exception:
            logging.getLogger(__name__).info("[TREE] Could not determine DB sequence count before tree build")

        # Optionally collapse highly similar run sequences (for tree readability)
        if getattr(args, 'collapse', False):
            try:
                threshold = float(getattr(args, 'collapse_threshold', 99.9))
            except Exception:
                threshold = 99.9
            # build qid -> tax mapping from classification output so we only
            # cluster sequences that share the same taxonomic assignment
            qid_to_tax = {}
            try:
                if class_out:
                    with open(class_out) as cf:
                        next(cf, None)
                        for line in cf:
                            parts = line.strip().split('\t')
                            if not parts:
                                continue
                            q = parts[0]
                            tax = parts[3] if len(parts) > 3 else None
                            qid_to_tax[q] = tax
            except Exception:
                qid_to_tax = {}

            # group short ids by tax
            taxa_groups = {}
            try:
                for h, s in read_fasta(str(mapped_derep)):
                    tax = qid_to_tax.get(h)
                    taxa_groups.setdefault(tax, []).append((h, s))
            except Exception:
                taxa_groups = {None: []}
            try:
                artifacts = collapse_fasta_within_taxa(
                    taxa_groups,
                    outdir,
                    'derep_short_collapsed.fasta',
                    'collapsed_map.tsv',
                    'collapsed_members.tsv',
                    threshold=threshold,
                    threads=threads,
                    log_prefix='[COLLAPSE]',
                )
                tree_fasta = artifacts.collapsed_path
                run_member_to_rep = artifacts.member_to_rep or {}
                # Build rep -> members mapping
                for mem, rep in run_member_to_rep.items():
                    run_rep_to_members.setdefault(rep, []).append(mem)
                logging.getLogger(__name__).info("[COLLAPSE] Wrote collapsed fasta %s (reps=%d)", tree_fasta, len(artifacts.collapsed_records))
            except Exception as e:
                logging.getLogger(__name__).warning("[COLLAPSE] Failed to write collapsed fasta: %s", e)

        # Optionally filter by phylum before tree building
        phylum_filter = getattr(args, 'phylum', None)
        if phylum_filter:
            try:
                # Parse classification to extract phylum-level assignments
                phylum_ids = set()
                if class_out:
                    with open(class_out) as cf:
                        next(cf, None)  # skip header
                        for line in cf:
                            parts = line.strip().split('\t')
                            if len(parts) >= 4:
                                qid = parts[0]
                                tax = parts[3]  # Full taxonomy string
                                # Extract phylum from taxonomy (first level after kingdom)
                                # Format is usually: k__Kingdom; p__Phylum; c__Class; ...
                                try:
                                    tax_parts = [t.strip() for t in tax.split(';')]
                                    for part in tax_parts:
                                        if part.startswith('p__'):
                                            phylum_name = part[3:].strip()
                                            # Normalize for comparison (case-insensitive, strip underscores)
                                            if phylum_name.lower().replace('_', '') == phylum_filter.lower().replace('_', ''):
                                                phylum_ids.add(qid)
                                            break
                                except Exception:
                                    pass
                
                # Filter tree_fasta to only include sequences in phylum_ids
                if phylum_ids:
                    filtered_fasta = Path(outdir) / f'tree_sequences_phylum_{phylum_filter}.fasta'
                    filtered_records = []
                    total_checked = 0
                    try:
                        for h, s in read_fasta(str(tree_fasta)):
                            total_checked += 1
                            if h in phylum_ids:
                                filtered_records.append((h, s))
                    except Exception:
                        pass

                    if filtered_records:
                        write_fasta(filtered_records, str(filtered_fasta))
                        logging.getLogger(__name__).info(
                            "[PHYLUM] Filtered tree sequences from %d to %d for phylum: %s",
                            total_checked, len(filtered_records), phylum_filter,
                        )
                        tree_fasta = filtered_fasta
                    else:
                        logging.getLogger(__name__).warning("[PHYLUM] No sequences found for phylum '%s' in classification results", phylum_filter)
                else:
                    logging.getLogger(__name__).warning("[PHYLUM] Could not extract phylum information from classification output")
            except Exception as e:
                logging.getLogger(__name__).warning("[PHYLUM] Phylum filtering failed: %s", e)

        # pass the Database object so the tree builder can export existing
        # preloaded sequences from the DB to form the backbone alignment
        _force_rebuild = bool(getattr(args, 'force_rebuild', False))
        try:
            tree_seq_count = sum(1 for line in open(str(tree_fasta)) if line.startswith(">"))
            logging.getLogger(__name__).info("[TREE] About to build/update tree with %d sequences from: %s", tree_seq_count, tree_fasta)
        except Exception as e:
            logging.getLogger(__name__).warning("[TREE] Could not count sequences in tree_fasta: %s", e)
        if _force_rebuild:
            logging.getLogger(__name__).info("[TREE] Force rebuild requested (--force-rebuild/--rebuild-tree); tree will be rebuilt from scratch including all datasets")
        tree.initialise_or_update_tree(
            ref_fasta=effective_ref,
            user_fasta=str(tree_fasta),
            outdir=outdir,
            db=db,
            db_dataset=None,
            threads=threads,
            preload_dir=None if _force_rebuild else getattr(args, 'preload_dir', None),
            force_rebuild=_force_rebuild,
            anchor_file=getattr(args, 'anchors', None),
            tree_method=getattr(args, 'tree_method', 'fasttree'),
        )
        try:
            warning_rows = tree.collect_tree_build_warnings(user_fasta=str(tree_fasta), anchor_file=getattr(args, 'anchors', None), db=db, db_dataset=None)
            warning_rows.extend(tree.summarize_alignment_quality(str(Path(outdir) / 'current_alignment.fasta')))
            if warning_rows:
                warn_path = tree.write_tree_warning_tsv(outdir, warning_rows)
                logging.getLogger(__name__).warning("[TREE] Tree/alignment warnings written to %s", warn_path)
        except Exception as e:
            logging.getLogger(__name__).warning("[TREE] Failed to summarise tree/alignment quality: %s", e)
    except Exception as e:
        logging.getLogger(__name__).warning("[TREE] Tree update failed: %s", e)

    # Build combined taxonomy for iTOL from DB sequences and current run results.
    combined_path = Path(outdir) / 'combined_taxonomy.tsv'
    merged = {}
    order = []

    try:
        preload_dir = getattr(args, 'preload_dir', None)
        preload_ids = None
        if preload_dir:
            try:
                p = Path(str(preload_dir))
                cand = p / 'preload_combined_taxonomy.tsv'
                if not cand.exists():
                    cand = p / 'combined_taxonomy.tsv'
                if cand.exists():
                    preload_ids = read_combined_taxonomy_ids(cand)
            except Exception:
                preload_ids = None

        base_rows = collect_db_taxonomy_rows(db, preload_ids if preload_ids else None)
        merged_rows = merge_combined_taxonomy_rows(base_rows, class_out or '', orig_to_short, db)
        for rid, tax, conf in merged_rows:
            merged[rid] = (tax if tax is not None else 'NA', conf if conf is not None else 'NA')
            order.append(rid)
    except Exception:
        # fallback: leave merged empty and continue
        pass

    # write merged combined taxonomy preserving order (DB-order first, then new)
    try:
        write_combined_taxonomy_tsv(
            combined_path,
            [(iid, merged.get(iid, ('NA', 'NA'))[0], merged.get(iid, ('NA', 'NA'))[1]) for iid in order],
        )
        logging.getLogger(__name__).info("[ITOL] Wrote combined taxonomy to %s", combined_path)
    except Exception as e:
        logging.getLogger(__name__).warning("[ITOL] Failed to write combined taxonomy: %s", e)

    # Load preload id_map for iTOL when --preload-dir is provided.
    id_map_for_itol = None
    preload_dir = getattr(args, 'preload_dir', None)
    if preload_dir:
        try:
            id_map = {}
            cand_map = _find_preferred_id_map(str(preload_dir))
            if cand_map:
                id_map = _load_id_map_from_tsv(cand_map, db=db)
            if id_map:
                id_map_for_itol = id_map
        except Exception:
            id_map_for_itol = None

        try:
            tree_path = Path(outdir) / 'current_tree.nwk'
            tfile = str(tree_path) if tree_path.exists() else _find_tree_file_in_dir(outdir)
            itol.generate_itol_colors(str(combined_path), outdir, user_color_csv=getattr(args, 'user_colors', None), id_map=id_map_for_itol, tree_file=tfile, phylum_groups=getattr(args, 'group_phyla', None))
            logging.getLogger(__name__).info("[ITOL] Generated iTOL color files in %s", outdir)
        except Exception as e:
            logging.getLogger(__name__).warning("[ITOL] Failed to generate iTOL files: %s", e)
    # Generate iTOL color files for the run output directory.
    try:
        tree_path = Path(outdir) / 'current_tree.nwk'
        tfile = str(tree_path) if tree_path.exists() else _find_tree_file_in_dir(outdir)
        itol.generate_itol_colors(str(combined_path), outdir, user_color_csv=getattr(args, 'user_colors', None), id_map=id_map_for_itol, tree_file=tfile, phylum_groups=getattr(args, 'group_phyla', None))
        logging.getLogger(__name__).info("[ITOL] Generated iTOL color files in %s", outdir)
    except Exception as e:
        logging.getLogger(__name__).warning("[ITOL] Failed to generate iTOL files: %s", e)
    # Optional: write functional annotation iTOL datasets when provided
    try:
        func_tsv = getattr(args, 'functional', None)
        if func_tsv:
            try:
                written = itol.write_functional_annotations(str(func_tsv), outdir, id_map=id_map_for_itol)
                logging.getLogger(__name__).info("[ITOL] Wrote functional annotation iTOL files: %s", ','.join(written) if written else '(none)')
            except Exception as e:
                logging.getLogger(__name__).warning("[ITOL] Functional annotations generation failed: %s", e)
    except Exception:
        pass
    # Draft rumen functional groups (auto-generated from output taxonomy)
    if getattr(args, 'draft_rumen_functions', False):
        try:
            tsv_out, itol_out = itol.generate_rumen_function_draft(
                str(combined_path), outdir, id_map=id_map_for_itol
            )
            if tsv_out:
                logging.getLogger(__name__).info("[ITOL] Draft rumen functional annotation: %s", tsv_out)
            if itol_out:
                logging.getLogger(__name__).info("[ITOL] Rumen functional iTOL file: %s", itol_out)
        except Exception as e:
            logging.getLogger(__name__).warning("[ITOL] Draft rumen functions generation failed: %s", e)
    # produce dataset membership band (preload vs run)
    try:
        combined_tax = combined_path
        ids_in_order = []
        if combined_tax.exists():
            with open(combined_tax) as ct:
                next(ct, None)
                for l in ct:
                    iid = l.strip().split('\t')[0]
                    if iid:
                        ids_in_order.append(iid)
        # Query DB for exact dataset membership per id to support arbitrary dataset names
        ds_map = {}
        try:
            with db.connect() as conn:
                cur = conn.cursor()
                placeholders = ','.join('?' for _ in ids_in_order) if ids_in_order else ''
                if placeholders:
                    cur.execute(f"SELECT id, dataset FROM sequences WHERE id IN ({placeholders})", tuple(ids_in_order))
                    for iid_row, ds in cur.fetchall():
                        ds_map[iid_row] = ds or ''
        except Exception:
            ds_map = {}

        membership_path = Path(outdir) / 'itol_dataset_membership.itol'
        itol.write_dataset_membership_strip(str(membership_path), ids_in_order, ds_map)
        logging.getLogger(__name__).info("[ITOL] Wrote dataset membership ITOL to %s", membership_path)
    except Exception as e:
        logging.getLogger(__name__).warning("[ITOL] Failed to write dataset membership ITOL: %s", e)

    # Also produce a novelty gradient colorstrip for NEW (run) sequences
    try:
        try:
            run_ids = [r[0] for r in mapped_records]
        except Exception:
            run_ids = []
        if not run_ids:
            try:
                with db.connect() as conn:
                    cur = conn.cursor()
                    cur.execute("SELECT id FROM sequences WHERE dataset = ?", (run_dataset,))
                    run_ids = [r[0] for r in cur.fetchall()]
            except Exception:
                run_ids = []

        if run_ids:
            novelty.build_run_novelty_itol(
                outdir,
                run_ids,
                str(mapped_derep),
                db,
                run_dataset,
                orig_to_short,
                threads=getattr(args, 'threads', None),
            )
    except Exception as e:
        logging.getLogger(__name__).warning("[ITOL] Failed to write novelty ITOL: %s", e)

    # placement warning summary combining classification + novelty support
    try:
        if class_out:
            warning_rows = build_placement_warning_rows(class_out, novelty_out, orig_to_short, db)
            if warning_rows:
                warn_path = write_placement_warning_tsv(Path(outdir) / 'placement_warnings.tsv', warning_rows)
                logging.getLogger(__name__).warning("[RUN] Placement warnings written to %s", warn_path)
            try:
                run_ids_for_assessment = [r[0] for r in mapped_records]
                # Determine which sequences are actually in the tree
                try:
                    _tree_ids = set()
                    for _h, _s in read_fasta(str(tree_fasta)):
                        _tree_ids.add(_h)
                except Exception:
                    _tree_ids = None

                # Collect alt-db taxonomy from DB for the run sequences
                alt_ref_dbs: list = []
                alt_taxonomies: dict = {}
                try:
                    alt_ref_dbs = db.get_alt_ref_dbs()
                    if alt_ref_dbs:
                        alt_taxonomies = db.get_taxonomy_alt_for_ids(run_ids_for_assessment)
                        logging.getLogger(__name__).info(
                            "[RUN] Found alt-db taxonomy for ref_dbs: %s", alt_ref_dbs
                        )
                except Exception as e:
                    logging.getLogger(__name__).warning("[RUN] Could not load alt-db taxonomy: %s", e)

                assessment_rows = build_sequence_assessment_rows(
                    run_ids_for_assessment,
                    class_out,
                    str(Path(outdir) / 'novelty_metrics.tsv'),
                    warning_rows,
                    orig_to_short,
                    db,
                    member_to_rep=run_member_to_rep if run_member_to_rep else None,
                    rep_to_members=run_rep_to_members if run_rep_to_members else None,
                    tree_ids=_tree_ids,
                    alt_taxonomies=alt_taxonomies,
                    alt_ref_dbs=alt_ref_dbs,
                )
                assess_path = write_sequence_assessment_tsv(Path(outdir) / 'sequence_assessment.tsv', assessment_rows)
                logging.getLogger(__name__).info("[RUN] Wrote sequence assessment to %s", assess_path)

                # ── Cluster-level reports + phylogenetic isolation ───────────
                try:
                    _tree_nwk = str(Path(outdir) / 'current_tree.nwk')
                    _cluster_summary, _cluster_csvs, _backup_tsv = _cluster_report.generate_cluster_reports(
                        outdir=outdir,
                        assessment_rows=assessment_rows,
                        tree_path=_tree_nwk if Path(_tree_nwk).exists() else None,
                    )
                    # Re-write sequence_assessment.tsv now that phylo_isolation /
                    # investigation_score have been filled in by generate_cluster_reports
                    write_sequence_assessment_tsv(Path(outdir) / 'sequence_assessment.tsv', assessment_rows)
                    if _cluster_summary:
                        logging.getLogger(__name__).info(
                            "[CLUSTER] Wrote cluster summary → %s  (%d per-cluster CSVs in %s/clusters/)",
                            _cluster_summary, len(_cluster_csvs), outdir,
                        )
                    if _backup_tsv:
                        logging.getLogger(__name__).info(
                            "[CLUSTER] Wrote backup candidates table → %s", _backup_tsv,
                        )
                except Exception as _ce:
                    logging.getLogger(__name__).warning("[CLUSTER] Cluster report generation failed: %s", _ce)
                # Emit a user-friendly summary of HIGH priority candidates
                try:
                    high_priority = [r for r in assessment_rows if r.get('sequencing_priority') == 'HIGH']
                    medium_priority = [r for r in assessment_rows if r.get('sequencing_priority') == 'MEDIUM']
                    collapsed_away = [r for r in assessment_rows if r.get('in_tree') == 'No']
                    logging.getLogger(__name__).info(
                        "[ASSESSMENT SUMMARY] %d sequences assessed: %d HIGH priority for sequencing, "
                        "%d MEDIUM priority. %d sequences were clustered and excluded from tree "
                        "(their representatives are in the tree). "
                        "See sequence_assessment.tsv for full details.",
                        len(assessment_rows), len(high_priority), len(medium_priority), len(collapsed_away),
                    )
                except Exception:
                    pass
            except Exception as e:
                logging.getLogger(__name__).warning("[RUN] Failed to build sequence assessment: %s", e)
    except Exception as e:
        logging.getLogger(__name__).warning("[RUN] Failed to build placement warnings: %s", e)

    # At end of run, write plain-text explanation files for outputs in outdir
    try:
        _write_output_explanations(outdir)
    except Exception:
        pass


def _detect_taxon_rank(taxon_query: str, rank_arg: str) -> str:
    """Return the single-letter rank key to filter on.

    Priority:
      1. Explicit ``--rank`` argument (if not 'auto')
      2. GTDB-style prefix embedded in the query (``p__``, ``f__``, ``g__`` …)
      3. Known domain keywords: ``archaea`` / ``bacteria`` → 'd'
      4. Fallback: 'p' (phylum)
    """
    from relict.taxonomy import RANK_ALIASES
    if rank_arg and rank_arg.lower() != 'auto':
        return RANK_ALIASES.get(rank_arg.lower(), rank_arg.lower()[:1])
    if '__' in taxon_query:
        prefix = taxon_query.split('__')[0].lower().strip()
        return RANK_ALIASES.get(prefix, prefix[:1] if prefix else 'p')
    # Bare domain-level keywords — map to 'd' so we filter at domain rank
    if taxon_query.strip().lower() in ('archaea', 'bacteria', 'eukarya', 'eukaryota'):
        return 'd'
    return 'p'


def _strip_rank_prefix(taxon: str) -> str:
    """Strip a GTDB-style rank prefix: 'p__Bacillota' → 'Bacillota'."""
    if '__' in taxon:
        return taxon.split('__', 1)[1].strip()
    return taxon.strip()


def _taxon_name_matches(val: str, query: str) -> bool:
    """Case-insensitive, underscore/space-flexible exact match."""
    import re as _re
    def _n(s): return _re.sub(r'[\s_]+', '_', str(s).strip().lower())
    return bool(val) and _n(val) == _n(query)


def cmd_subtree(args):
    """Build a focused subtree for a specific taxon from an existing DB.

    Fast path: if a ``current_alignment.fasta`` already exists in
    ``--from-dir`` (or the output directory), the matching sequences are
    extracted from that *pre-built* alignment and FastTree is run on the
    subset — no re-alignment needed.

    Slow path (fallback): if no existing alignment is found, the matching
    sequences are exported from the DB and a full MAFFT + FastTree build
    is performed (same as a normal run).
    """
    from relict.pipeline import itol as itol_mod
    from relict.pipeline import tree as tree_mod  # noqa: F401 (used in _build_subtree)
    from relict.pipeline.workflow_helpers import write_combined_taxonomy_tsv
    from relict.taxonomy import parse_taxon_string
    from relict.utils.fasta import read_fasta, write_fasta  # noqa: F401

    outdir = args.out
    os.makedirs(outdir, exist_ok=True)
    _configure_logging(outdir)
    log = logging.getLogger(__name__)

    taxon_query = args.taxon
    rank_key = _detect_taxon_rank(taxon_query, getattr(args, 'rank', 'auto'))
    taxon_clean = _strip_rank_prefix(taxon_query)
    threads = getattr(args, 'threads', 4)
    min_seqs = getattr(args, 'min_seqs', 3)
    ref_fasta = getattr(args, 'ref', None)
    anchor_file = getattr(args, 'anchors', None)
    from_dir = getattr(args, 'from_dir', None) or outdir
    no_tree = getattr(args, 'no_tree', False)

    log.info(
        "[SUBTREE] Taxon query: '%s'  rank: '%s'  clean name: '%s'  from-dir: %s",
        taxon_query, rank_key, taxon_clean, from_dir,
    )

    # ── Query DB ──────────────────────────────────────────────────────────────
    db = Database(args.db)
    db.initialise()

    with db.connect() as conn:
        cur = conn.cursor()
        cur.execute(
            "SELECT s.id, s.sequence, t.taxonomy, t.confidence, s.dataset "
            "FROM sequences s LEFT JOIN taxonomy t ON s.id = t.id"
        )
        all_rows = cur.fetchall()

    # ── Filter to matching taxon ──────────────────────────────────────────────
    matched = []
    for sid, seq, tax, conf, dataset in all_rows:
        if not tax:
            continue
        parsed = parse_taxon_string(tax)
        val = parsed.get(rank_key, '')
        if _taxon_name_matches(val, taxon_clean):
            matched.append((sid, seq or '', tax, conf or 'NA', dataset or ''))

    log.info("[SUBTREE] DB has %d sequences total; %d match taxon '%s' at rank '%s'",
             len(all_rows), len(matched), taxon_clean, rank_key)

    if len(matched) < min_seqs:
        log.warning(
            "[SUBTREE] Only %d sequences matched (minimum required: %d). "
            "Check taxon spelling and rank. Available phyla in DB: %s",
            len(matched), min_seqs,
            sorted({parse_taxon_string(r[2]).get('p', 'unknown')
                    for r in all_rows if r[2]}),
        )
        return

    # Collect matched IDs as a set for fast lookup
    matched_ids = {sid for sid, *_ in matched}

    # ── Write taxonomy TSV ────────────────────────────────────────────────────
    combined_tax_path = Path(outdir) / 'subtree_combined_taxonomy.tsv'
    write_combined_taxonomy_tsv(
        combined_tax_path,
        [(sid, tax, conf) for sid, _, tax, conf, _ in matched],
    )
    log.info("[SUBTREE] Wrote taxonomy for %d sequences → %s", len(matched), combined_tax_path)

    # ── Write summary TSV ─────────────────────────────────────────────────────
    summary_path = Path(outdir) / 'subtree_sequence_list.tsv'
    with open(summary_path, 'w') as sf:
        sf.write('ID\tTaxonomy\tConfidence\tDataset\n')
        for sid, _, tax, conf, dataset in matched:
            sf.write(f"{sid}\t{tax}\t{conf}\t{dataset}\n")
    log.info("[SUBTREE] Wrote sequence list → %s", summary_path)

    # ── Tree building ─────────────────────────────────────────────────────────
    tree_path = Path(outdir) / 'subtree_tree.nwk'

    if not no_tree:
        _build_subtree(
            matched=matched,
            matched_ids=matched_ids,
            from_dir=from_dir,
            outdir=outdir,
            tree_path=tree_path,
            ref_fasta=ref_fasta,
            anchor_file=anchor_file,
            threads=threads,
            taxon_clean=taxon_clean,
            log=log,
        )
    else:
        log.info("[SUBTREE] Tree build skipped (--no-tree).")

    # ── iTOL files ────────────────────────────────────────────────────────────
    try:
        tfile = str(tree_path) if tree_path.exists() else None
        itol_mod.generate_itol_colors(
            str(combined_tax_path),
            outdir,
            tree_file=tfile,
            phylum_groups=getattr(args, 'group_phyla', None),
        )
        log.info("[SUBTREE] Generated iTOL colour files in %s", outdir)
    except Exception as e:
        log.warning("[SUBTREE] iTOL colour generation failed: %s", e)

    # Optional: write functional annotation datasets when provided
    try:
        func_tsv = getattr(args, 'functional', None)
        if func_tsv:
            try:
                written = itol_mod.write_functional_annotations(str(func_tsv), outdir, id_map=None)
                log.info("[SUBTREE] Wrote functional annotation iTOL files: %s", ','.join(written) if written else '(none)')
            except Exception as e:
                log.warning("[SUBTREE] Functional annotations generation failed: %s", e)
    except Exception:
        pass
    # Draft rumen functional groups
    if getattr(args, 'draft_rumen_functions', False):
        try:
            tsv_out, itol_out = itol_mod.generate_rumen_function_draft(
                str(combined_tax_path), outdir, id_map=None
            )
            if tsv_out:
                log.info("[SUBTREE] Draft rumen functional annotation: %s", tsv_out)
            if itol_out:
                log.info("[SUBTREE] Rumen functional iTOL file: %s", itol_out)
        except Exception as e:
            log.warning("[SUBTREE] Draft rumen functions generation failed: %s", e)

    # ── Dataset membership strip ──────────────────────────────────────────────
    # One colour per dataset label so users can see which sequences came from
    # which dataset (preload, run batch, etc.) in the same tree view.
    try:
        ids_in_order = [sid for sid, *_ in matched]
        ds_map = {sid: (dataset or 'unknown') for sid, _, _tax, _conf, dataset in matched}
        membership_path = Path(outdir) / 'itol_dataset_membership.itol'
        itol_mod.write_dataset_membership_strip(
            str(membership_path), ids_in_order, ds_map,
            dataset_label='Dataset membership',
        )
        log.info("[SUBTREE] Wrote dataset membership strip → %s", membership_path)
    except Exception as e:
        log.warning("[SUBTREE] Dataset membership strip failed: %s", e)

    # ── Done ──────────────────────────────────────────────────────────────────
    log.info(
        "[SUBTREE] Complete. %d sequences, taxon='%s', output=%s",
        len(matched), taxon_query, outdir,
    )
    print(
        f"[subtree] Done: {len(matched)} sequences for '{taxon_query}' "
        f"→ {outdir}"
    )


def _build_subtree(
    matched, matched_ids, from_dir, outdir, tree_path,
    ref_fasta, anchor_file, threads, taxon_clean, log,
):
    """Build a FastTree for the given matched sequences.

    Fast path: filter an existing ``current_alignment.fasta`` from
    ``from_dir`` to the matched IDs + any reference anchors already
    embedded in that alignment, then run FastTree directly — no
    re-alignment required.

    Slow path: write the raw sequences as a FASTA and run the full
    ``initialise_or_update_tree`` pipeline (MAFFT → FastTree).
    """
    from relict.pipeline.tree import (
        _run_fasttree, _make_unique_fasta,
        is_ref_anchor, get_anchor_file,
    )
    from relict.utils.fasta import read_fasta, write_fasta
    from relict.pipeline import tree as tree_mod
    import re

    out = Path(outdir)
    resolved_anchor = get_anchor_file(anchor_file)

    # ── Fast path ─────────────────────────────────────────────────────────────
    # Look for an existing alignment in from_dir (or outdir)
    aln_candidates = [
        Path(from_dir) / 'current_alignment.fasta',
        Path(outdir) / 'current_alignment.fasta',
    ]
    existing_aln = next((p for p in aln_candidates if p.exists()), None)

    if existing_aln:
        log.info("[SUBTREE] Fast path: filtering existing alignment %s", existing_aln)
        # Build a normalised ID lookup for the matched IDs
        def _norm_id(x):
            x = str(x).split()[0]
            if '|' in x:
                x = x.split('|')[-1]
            x = re.sub(r'_[0-9]+$', '', x)
            return x.lower()

        matched_norm = {_norm_id(sid) for sid in matched_ids}

        kept = []
        anchors_included = 0
        for header, seq in read_fasta(str(existing_aln)):
            hid = header.split()[0]
            if is_ref_anchor(hid):
                kept.append((header, seq))
                anchors_included += 1
            elif hid in matched_ids or _norm_id(hid) in matched_norm:
                kept.append((header, seq))

        n_data = len(kept) - anchors_included
        log.info(
            "[SUBTREE] Fast path: %d data sequences + %d anchors extracted from alignment",
            n_data, anchors_included,
        )

        if n_data < 3:
            log.warning(
                "[SUBTREE] Only %d data sequences found in existing alignment "
                "(IDs may not match). Falling back to slow path.", n_data,
            )
            existing_aln = None   # trigger slow path below
        else:
            subtree_aln = out / 'subtree_alignment.fasta'
            write_fasta(kept, str(subtree_aln))
            fasta_for_tree, id_map = _make_unique_fasta(str(subtree_aln), out)
            if _run_fasttree(Path(fasta_for_tree), tree_path):
                _finalise_tree_subtree(out, id_map, tree_path, log)
            else:
                log.warning("[SUBTREE] FastTree failed on fast path")
            return

    # ── Slow path ─────────────────────────────────────────────────────────────
    log.info("[SUBTREE] Slow path: full MAFFT + FastTree build")
    seqs_fasta = out / 'subtree_input_sequences.fasta'
    write_fasta([(sid, seq) for sid, seq, *_ in matched], str(seqs_fasta))

    # Empty FASTA as user_fasta; all sequences go in via the FASTA directly
    # by passing them as user_fasta and disabling DB pull (db=None)
    try:
        tree_mod.initialise_or_update_tree(
            ref_fasta=ref_fasta or '',
            user_fasta=str(seqs_fasta),
            outdir=outdir,
            db=None,
            threads=threads,
            anchor_file=anchor_file,
            force_rebuild=True,
        )
        # The standard pipeline writes to current_tree.nwk; copy/rename
        default_tree = out / 'current_tree.nwk'
        if default_tree.exists():
            import shutil
            shutil.copy2(str(default_tree), str(tree_path))
            log.info("[SUBTREE] Slow path tree written → %s", tree_path)
    except Exception as e:
        log.warning("[SUBTREE] Slow path tree build failed: %s", e)


def _finalise_tree_subtree(out: Path, id_map: dict, tree_path: Path, log) -> None:
    """Remap IDs, prune anchors, and label internal nodes for subtree output."""
    from relict.pipeline.tree import (
        _repair_legacy_internal_node_labels, _label_internal_nodes,
        _prune_anchor_leaves, REF_ANCHOR_PREFIX,
    )
    if not tree_path.exists():
        return
    newick = tree_path.read_text()
    if id_map:
        for new_id, orig_id in id_map.items():
            newick = newick.replace(new_id, orig_id)
    before = newick.count(REF_ANCHOR_PREFIX)
    newick = _prune_anchor_leaves(newick)
    after = newick.count(REF_ANCHOR_PREFIX)
    if before:
        log.info("[SUBTREE] Pruned %d anchor leaves from subtree newick", before - after)
    newick = _repair_legacy_internal_node_labels(newick)
    newick = _label_internal_nodes(newick)
    tree_path.write_text(newick)
    log.info("[SUBTREE] Subtree finalised → %s", tree_path)


def cmd_regen_itol(args):
    db = Database(args.db)
    db.initialise()
    outdir = args.out
    os.makedirs(outdir, exist_ok=True)
    _configure_logging(outdir)
    log = logging.getLogger(__name__)

    try:
        with db.connect() as conn:
            cur = conn.cursor()
            # If the outdir contains a preload combined taxonomy file, prefer
            # to use only those IDs so the regenerated iTOL matches preload
            # statistics. This mirrors the behaviour used during `run` when a
            # --preload-dir is supplied.
            try:
                p = Path(outdir)
                # Prefer an existing combined taxonomy file in the outdir, but
                # only if it contains data (header + >=1 data row). Fall back
                # to the DB-wide query otherwise to avoid regenerating empty
                # iTOL outputs when a stub file exists.
                preload_file = None
                for cand_name in ('preload_combined_taxonomy.tsv', 'combined_taxonomy.tsv'):
                    cand = p / cand_name
                    if not cand.exists():
                        continue
                    try:
                        # count lines cheaply
                        with open(cand) as pf:
                            cnt = sum(1 for _ in pf)
                        if cnt > 1:
                            preload_file = cand
                            break
                        else:
                            # file exists but only header or empty -> ignore
                            continue
                    except Exception:
                        continue
                preload_ids = None
                if preload_file is not None:
                    preload_ids = read_combined_taxonomy_ids(preload_file)
                else:
                    preload_ids = None
            except Exception:
                preload_ids = None

            preload_id_list = preload_ids if isinstance(preload_ids, list) else []

            # fetch id, taxonomy, confidence and dataset for all sequences (or filter)
            if getattr(args, 'include_datasets', None):
                ds_list = [d.strip() for d in args.include_datasets.split(',') if d.strip()]
                placeholders = ','.join('?' for _ in ds_list)
                cur.execute(f"SELECT s.id, t.taxonomy, t.confidence, s.dataset FROM sequences s LEFT JOIN taxonomy t ON s.id = t.id WHERE s.dataset IN ({placeholders})", tuple(ds_list))
            elif preload_id_list:
                placeholders = ','.join('?' for _ in preload_id_list)
                cur.execute(f"SELECT s.id, t.taxonomy, t.confidence, s.dataset FROM sequences s LEFT JOIN taxonomy t ON s.id = t.id WHERE s.id IN ({placeholders})", tuple(preload_id_list))
            else:
                cur.execute("SELECT s.id, t.taxonomy, t.confidence, s.dataset FROM sequences s LEFT JOIN taxonomy t ON s.id = t.id")
            rows = cur.fetchall()
    except Exception as e:
        log.warning("[REGEN-ITOL] Failed to query DB: %s", e)
        return

    try:
        kingdom = getattr(args, 'kingdom', None)
        if kingdom:
            kingdom_text = str(kingdom)
            rows = [row for row in rows if row[1] and taxonomy_matches_kingdom(str(row[1]), kingdom_text)]
            log.info("[REGEN-ITOL] Retained %d rows after kingdom filter=%s", len(rows), kingdom)
    except Exception as e:
        log.warning("[REGEN-ITOL] Kingdom filter failed: %s", e)

    combined_path = Path(outdir) / 'combined_taxonomy.tsv'
    try:
        write_combined_taxonomy_tsv(combined_path, [(rid, tax, conf) for rid, tax, conf, ds in rows])
        log.info("[REGEN-ITOL] Wrote combined taxonomy for %d ids to %s", len(rows), combined_path)
    except Exception as e:
        log.warning("[REGEN-ITOL] Failed to write combined taxonomy: %s", e)
        return

    # call itol generator
    try:
        tree_path = Path(outdir) / 'current_tree.nwk'
        tfile = str(tree_path) if tree_path.exists() else _find_tree_file_in_dir(outdir)
        itol.generate_itol_colors(str(combined_path), outdir, tree_file=tfile, phylum_groups=getattr(args, 'group_phyla', None))
        log.info("[REGEN-ITOL] Generated iTOL color files in %s", outdir)
    except Exception as e:
        log.warning("[REGEN-ITOL] Color generation failed: %s", e)

    # Optional: write functional annotation datasets when provided
    try:
        func_tsv = getattr(args, 'functional', None)
        if func_tsv:
            try:
                written = itol.write_functional_annotations(str(func_tsv), outdir, id_map=None)
                log.info("[REGEN-ITOL] Wrote functional annotation iTOL files: %s", ','.join(written) if written else '(none)')
            except Exception as e:
                log.warning("[REGEN-ITOL] Functional annotations generation failed: %s", e)
    except Exception:
        pass
    # Draft rumen functional groups
    if getattr(args, 'draft_rumen_functions', False):
        try:
            tsv_out, itol_out = itol.generate_rumen_function_draft(str(combined_path), outdir, id_map=None)
            if tsv_out:
                log.info("[REGEN-ITOL] Draft rumen functional annotation: %s", tsv_out)
            if itol_out:
                log.info("[REGEN-ITOL] Rumen functional iTOL file: %s", itol_out)
        except Exception as e:
            log.warning("[REGEN-ITOL] Draft rumen functions generation failed: %s", e)

    # build dataset membership strip
    try:
        ids_in_order = [r[0] for r in rows]
        ds_map = {r[0]: (r[3] or '') for r in rows}
        membership_path = Path(outdir) / 'itol_dataset_membership.itol'
        itol.write_dataset_membership_strip(str(membership_path), ids_in_order, ds_map)
        log.info("[REGEN-ITOL] Wrote dataset membership ITOL to %s", membership_path)
    except Exception as e:
        log.warning("[REGEN-ITOL] Failed to build/write dataset membership ITOL: %s", e)

    # write explanations for regenerated outputs
    try:
        _write_output_explanations(outdir)
    except Exception:
        pass


def build_parser():
    parser = argparse.ArgumentParser(
        prog='relict',
        description=(
            'Relict — reference-aware 16S novelty and phylogenetic context tool.\n\n'
            'Subcommands:\n'
            '  preload     Load a baseline dataset (e.g. Hungate) and build the backbone tree.\n'
            '  run         Process new sequences against the baseline; score novelty and update the tree.\n'
            '  subtree     Extract a focused tree and iTOL files for a specific taxon from an existing DB.\n'
            '  regen-itol  Regenerate iTOL colour files from an existing DB without re-running analysis.\n\n'
            'Typical workflow:\n'
            '  1. relict preload  --fasta baseline.fasta --db project.db --dataset Hungate --ref gtdb.fna --classify --build-tree -o preload_out\n'
            '  2. relict run      --input new_seqs.fasta --db project.db --dataset Batch1  --ref gtdb.fna --preload-dir preload_out -o run_out\n'
            '  3. relict subtree  --db project.db --taxon archaea --from-dir preload_out -o archaea_out\n'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest='command')

    # ── preload ───────────────────────────────────────────────────────────────
    preload = sub.add_parser(
        'preload',
        help='Load a baseline dataset and build the reference tree.',
        description=(
            'Load a baseline FASTA dataset (e.g. Hungate 16S) into the DB, optionally classify\n'
            'sequences against a reference (GTDB/SILVA), collapse near-identical sequences,\n'
            'and build the backbone phylogenetic tree.\n\n'
            'This is always the first step. All subsequent `run` calls measure novelty against\n'
            'sequences stored by this command.'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    preload.add_argument('--fasta', required=True,
        help='Input FASTA file containing the baseline sequences to load.')
    preload.add_argument('--db', required=True,
        help='Path to the Relict SQLite database (created if it does not exist).')
    preload.add_argument('-o', '--out', required=False, default='.',
        help='Output directory for tree, iTOL files, and reports (default: current directory).')
    preload.add_argument('--dataset', required=True,
        help='Label for this dataset stored in the DB (e.g. Hungate). Used to colour iTOL strips.')
    preload.add_argument('--shorten-ids', dest='shorten_ids',
        action=argparse.BooleanOptionalAction, default=False,
        help='Replace input headers with compact IDs (e.g. HUN001). Use --no-shorten-ids to replace source names.')
    preload.add_argument('--classify', action='store_true',
        help='Classify sequences against --ref and store taxonomy in the DB. Requires --ref.')
    preload.add_argument('--build-tree', action='store_true',
        help='Build the backbone MAFFT + FastTree phylogenetic tree after loading.')
    preload.add_argument('--ref', required=False,
        help='Reference FASTA (GTDB/SILVA reps) for classification and tree orientation. Preferred over --taxa-assignments for externally classified inputs.')
    preload.add_argument('--taxa', required=False,
        help='Tab-separated taxonomy file matching the IDs in --ref (id<TAB>lineage). Optional when --ref FASTA headers already contain GTDB lineages.')
    preload.add_argument('--ref-name', dest='ref_name', required=False, default=None,
        help='Display name for the primary reference database (default: derived from --ref filename). Used to label taxonomy columns.')
    preload.add_argument('--alt-ref', dest='alt_ref', action='append', default=None, metavar='FASTA',
        help='Additional reference FASTA to classify against (repeatable). Produces extra taxonomy columns in output files.')
    preload.add_argument('--alt-taxa', dest='alt_taxa', action='append', default=None, metavar='TSV',
        help='Taxa TSV for the corresponding --alt-ref (positionally paired; repeatable).')
    preload.add_argument('--alt-ref-name', dest='alt_ref_name', action='append', default=None, metavar='NAME',
        help='Display name for the corresponding --alt-ref (positionally paired; repeatable). Default: derived from filename.')
    preload.add_argument('--main-ref', dest='main_ref', required=False, default=None,
        help='Name of the reference database to use as the primary taxonomy source (default: primary --ref). Must match one of the --ref-name / --alt-ref-name values.')
    preload.add_argument('--taxa-assignments', '--taxa-aasignments',
        dest='taxa_assignments', required=False,
        help='Pre-computed taxonomy assignments for the INPUT sequences (TSV: query_id<TAB>lineage, or a FASTA with embedded lineages). Use this instead of --classify when you already have taxonomy.')
    preload.add_argument('--collapse', action='store_true',
        help='Collapse sequences that share ≥ --collapse-threshold identity AND the same taxonomy into a single representative for the tree. Saves time and reduces visual clutter.')
    preload.add_argument('--collapse-threshold', type=float, default=99.8,
        help='Identity threshold (percent) for collapsing duplicate-like sequences (default: 99.8).')
    preload.add_argument('--kingdom', required=False,
        help='Keep only sequences belonging to this kingdom/domain (case-insensitive). E.g. "bacteria" to drop archaeal sequences.')
    preload.add_argument('--anchors', required=False, default=None,
        help='Custom reference anchor FASTA for tree topology scaffolding. Defaults to the 26-sequence bundled anchor set (src/relict/data/reference_anchors.fasta).')
    preload.add_argument('--threads', type=int, required=False, default=4,
        help='Number of CPU threads for MAFFT and VSEARCH (default: 4).')
    preload.add_argument('--tree-method', dest='tree_method',
        choices=['fasttree', 'iqtree', 'iqtree-fast'], default='fasttree',
        help=(
            'Phylogenetic tree-building backend (default: fasttree). '
            'fasttree: approximate ML, GTR+CAT — fast. '
            'iqtree: full ML, GTR+G+I — more accurate and stable topology (slower; '
            'recommended for publication-quality trees or when FastTree produces unstable clades). '
            'iqtree-fast: IQ-TREE 2 with -fast flag — good compromise for exploratory runs. '
            'Requires iqtree2 in PATH when using iqtree/iqtree-fast.'
        ),
    )
    preload.add_argument('--colors', required=False,
        help='CSV file mapping sequence IDs to custom hex colours for iTOL (columns: id, color).')
    preload.add_argument(
        '--group-phyla', dest='group_phyla', action='append', default=None, metavar='SPEC',
        help=(
            'Collapse multiple phyla into a single colour in iTOL legends. Repeatable. '
            'Formats: "archaea" (all archaeal phyla), "bacteria" (all bacterial phyla), '
            '"Bacillota,Bacillota_I" (explicit list; label = first name), '
            '"Firmicutes:Bacillota,Bacillota_I" (named group).'
        ),
    )
    preload.add_argument('--functional', dest='functional', required=False, default=None,
        help=(
            'TSV file mapping sequence IDs to functional attributes (pathways, functions, '
            'traits, scores, etc.). Header row required; first column = sequence ID; '
            'subsequent columns = one functional attribute each. '
            'One iTOL file is generated per column: binary (0/1/yes/no) → DATASET_BINARY, '
            'numeric → DATASET_SIMPLEBAR, categorical → DATASET_COLORSTRIP.'
        ),
    )
    preload.add_argument('--draft-rumen-functions', dest='draft_rumen_functions',
        action='store_true', default=False,
        help=(
            'Auto-generate a draft rumen functional-group annotation from the output taxonomy. '
            'Maps each sequence to a broad ruminant microbiome functional category '
            '(e.g. Cellulolytic/Fibrolytic, Methanogenic Archaea, Butyrate Producers) '
            'and writes rumen_functions_draft.tsv + itol_func_Rumen_Functional_Group.itol. '
            'The draft TSV can be edited and re-supplied via --functional in future runs.'
        ),
    )

    # ── run ───────────────────────────────────────────────────────────────────
    run = sub.add_parser(
        'run',
        help='Process new sequences against the baseline; score novelty and update the tree.',
        description=(
            'Classify new sequences, score their novelty and neighbourhood density against the\n'
            'baseline (preload + all prior run datasets stored in the DB), and update the tree.\n\n'
            'Novelty is always relative to YOUR submitted data, not the full external reference.\n'
            'Each successive run extends the baseline, so scores become increasingly precise.'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    run.add_argument('--input', required=True,
        help='FASTA file of new sequences to analyse.')
    run.add_argument('--db', required=True,
        help='Path to the Relict SQLite database (must have been initialised with `relict preload`).')
    run.add_argument('-o', '--out', required=True,
        help='Output directory for this run (sequence_assessment.tsv, novelty_metrics.tsv, tree, iTOL files, etc.).')
    run.add_argument('--dataset', required=True,
        help='Label for this batch of sequences stored in the DB (e.g. Batch1). Used in iTOL dataset-membership strip.')
    run.add_argument('--ref', required=False,
        help='Reference FASTA (GTDB/SILVA reps) used for classification and tree orientation. Same file used in preload.')
    run.add_argument('--taxa', required=False,
        help='Tab-separated taxonomy file matching the IDs in --ref (id<TAB>lineage).')
    run.add_argument('--ref-name', dest='ref_name', required=False, default=None,
        help='Display name for the primary reference database (default: derived from --ref filename).')
    run.add_argument('--alt-ref', dest='alt_ref', action='append', default=None, metavar='FASTA',
        help='Additional reference FASTA to classify against (repeatable). Adds extra taxonomy columns to sequence_assessment.tsv and taxonomy_all_dbs.tsv.')
    run.add_argument('--alt-taxa', dest='alt_taxa', action='append', default=None, metavar='TSV',
        help='Taxa TSV for the corresponding --alt-ref (positionally paired; repeatable).')
    run.add_argument('--alt-ref-name', dest='alt_ref_name', action='append', default=None, metavar='NAME',
        help='Display name for the corresponding --alt-ref (positionally paired; repeatable).')
    run.add_argument('--main-ref', dest='main_ref', required=False, default=None,
        help='Name of the reference database to treat as primary (drives the main Taxonomy column). Default: primary --ref.')
    run.add_argument('--taxa-assignments', '--taxa-aasignments',
        dest='taxa_assignments', required=False,
        help='Pre-computed taxonomy for the INPUT sequences (TSV: query_id<TAB>lineage, or embedded-lineage FASTA).')
    run.add_argument('--preload-dir', dest='preload_dir', required=False,
        help='Path to the preload output directory. Used to seed the tree backbone alignment so only new sequences need aligning.')
    run.add_argument('--shorten-ids', dest='shorten_ids',
        action=argparse.BooleanOptionalAction, default=False,
        help='Replace input headers with compact IDs. Use --no-shorten-ids to replace source names.')
    run.add_argument('--min-len', dest='min_len', type=int, default=1200,
        help='Minimum sequence length to retain (bp, default: 1200). Shorter sequences are filtered out.')
    run.add_argument('--max-n', dest='max_n', type=int, default=5,
        help='Maximum number of ambiguous (N) bases allowed (default: 5).')
    run.add_argument('--collapse', action='store_true',
        help='Collapse near-identical same-taxonomy sequences into representatives for the tree.')
    run.add_argument('--collapse-threshold', type=float, default=99.8,
        help='Identity threshold (percent) for collapsing (default: 99.8).')
    run.add_argument('--kingdom', required=False,
        help='Keep only sequences belonging to this kingdom/domain (e.g. bacteria).')
    run.add_argument('--phylum', required=False,
        help='Filter iTOL output to sequences assigned to this phylum (e.g. Bacillota). Does not affect novelty scoring.')
    run.add_argument('--target', required=False, default=None,
        help=(
            'FASTA of sequences to measure novelty against instead of the DB. '
            'Leave unset to use all sequences previously stored in the DB (recommended).'
        ),
    )
    run.add_argument('--force-rebuild', '--rebuild-tree', dest='force_rebuild', action='store_true', default=False,
        help='Rebuild the entire tree from scratch even when an existing alignment is present. '
             'When combined with --preload-dir, ignores the preload backbone and jointly estimates '
             'tree topology across all datasets. (--rebuild-tree is an alias for this flag.)')
    run.add_argument('--anchors', required=False, default=None,
        help='Custom reference anchor FASTA for tree scaffolding. Defaults to bundled anchors.')
    run.add_argument('--threads', dest='threads', type=int, default=4,
        help='CPU threads for MAFFT and VSEARCH (default: 4).')
    run.add_argument('--tree-method', dest='tree_method',
        choices=['fasttree', 'iqtree', 'iqtree-fast'], default='fasttree',
        help=(
            'Phylogenetic tree-building backend (default: fasttree). '
            'fasttree: approximate ML, GTR+CAT. '
            'iqtree: full ML, GTR+G+I (recommended for production/publication runs). '
            'iqtree-fast: IQ-TREE 2 with -fast flag (good for exploratory incremental runs).'
        ),
    )
    run.add_argument('--user-colors', dest='user_colors', required=False,
        help='CSV file mapping sequence IDs to custom hex colours for iTOL (columns: id, color).')
    run.add_argument(
        '--group-phyla', dest='group_phyla', action='append', default=None, metavar='SPEC',
        help=(
            'Collapse multiple phyla into one colour in iTOL legends. Repeatable. '
            'Formats: "archaea", "bacteria", "Bacillota,Bacillota_I", "Firmicutes:Bacillota,Bacillota_I".'
        ),
    )
    run.add_argument('--functional', dest='functional', required=False, default=None,
        help=(
            'TSV file mapping sequence IDs to functional attributes. '
            'Header row required; first column = sequence ID; subsequent columns = functional attributes. '
            'Generates one iTOL file per column (DATASET_BINARY / DATASET_SIMPLEBAR / DATASET_COLORSTRIP).'
        ),
    )
    run.add_argument('--draft-rumen-functions', dest='draft_rumen_functions',
        action='store_true', default=False,
        help=(
            'Auto-generate a draft rumen functional-group annotation from the output taxonomy. '
            'Writes rumen_functions_draft.tsv and itol_func_Rumen_Functional_Group.itol. '
            'The draft TSV can be edited and re-supplied via --functional in later runs.'
        ),
    )

    # ── regen-itol ────────────────────────────────────────────────────────────
    regen = sub.add_parser(
        'regen-itol',
        help='Regenerate iTOL colour files from an existing DB without re-running analysis.',
        description=(
            'Re-generate all iTOL colour strips (phylum, family, genus, dataset membership)\n'
            'from the taxonomy already stored in the DB. Useful after changing --group-phyla\n'
            'options or after manually editing the database.'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    regen.add_argument('--db', required=True,
        help='Path to the Relict SQLite database.')
    regen.add_argument('-o', '--out', required=True,
        help='Output directory where iTOL files will be written (should be the preload or run output dir).')
    regen.add_argument('--include-datasets', required=False,
        help='Comma-separated list of dataset names to include (default: all datasets in the DB).')
    regen.add_argument('--kingdom', required=False,
        help='Include only sequences whose taxonomy contains this kingdom string (e.g. bacteria).')
    regen.add_argument(
        '--group-phyla', dest='group_phyla', action='append', default=None, metavar='SPEC',
        help=(
            'Collapse multiple phyla into one colour in iTOL legends. Repeatable. '
            'Formats: "archaea", "bacteria", "Bacillota,Bacillota_I", "Firmicutes:Bacillota,Bacillota_I".'
        ),
    )
    regen.add_argument('--functional', dest='functional', required=False, default=None,
        help=(
            'TSV file mapping sequence IDs to functional attributes. '
            'Generates one iTOL file per column (DATASET_BINARY / DATASET_SIMPLEBAR / DATASET_COLORSTRIP).'
        ),
    )
    regen.add_argument('--draft-rumen-functions', dest='draft_rumen_functions',
        action='store_true', default=False,
        help='Auto-generate rumen functional-group iTOL annotation from stored taxonomy.',
    )

    # ── subtree ───────────────────────────────────────────────────────────────
    subtree = sub.add_parser(
        'subtree',
        help='Build a focused tree and iTOL files for a specific taxon from an existing DB.',
        description=(
            'Extract all sequences matching a given taxon from the DB and build a focused\n'
            'phylogenetic tree for that group only.\n\n'
            'Fast path: if --from-dir points to a directory containing current_alignment.fasta\n'
            '(a preload or run output), sequences are sliced from the pre-built alignment and\n'
            'FastTree is run directly — no MAFFT re-alignment needed (~seconds for hundreds of seqs).\n\n'
            'Slow path: if no existing alignment is found, a full MAFFT + FastTree build is run.\n\n'
            'Taxon formats accepted:\n'
            '  archaea, bacteria          → all sequences at domain level\n'
            '  Bacillota, Bacteroidota    → phylum name (auto-detected)\n'
            '  p__Bacillota               → GTDB-prefixed phylum\n'
            '  f__Lachnospiraceae         → GTDB-prefixed family\n'
            '  g__Ruminococcus            → GTDB-prefixed genus\n'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    subtree.add_argument('--db', required=True,
        help='Path to the Relict SQLite database.')
    subtree.add_argument('-o', '--out', required=True,
        help='Output directory for the subtree results.')
    subtree.add_argument('--taxon', required=True,
        help='Taxon to extract (see description above for accepted formats).')
    subtree.add_argument(
        '--rank', required=False, default='auto',
        choices=['auto', 'domain', 'd', 'phylum', 'p', 'class', 'c',
                 'order', 'o', 'family', 'f', 'genus', 'g', 'species', 's'],
        help='Taxonomic rank to filter on. Default "auto" detects the rank from the taxon name or prefix.',
    )
    subtree.add_argument('--from-dir', dest='from_dir', required=False, default=None,
        help=(
            'Existing preload or run output directory containing current_alignment.fasta. '
            'Enables the fast path (sequences extracted from the existing MSA; no re-alignment).'
        ),
    )
    subtree.add_argument('--ref', required=False,
        help='Reference FASTA for orientation correction (slow-path full build only).')
    subtree.add_argument('--anchors', required=False, default=None,
        help='Custom reference anchor FASTA. Defaults to bundled anchors (26 type-strain sequences).')
    subtree.add_argument('--threads', type=int, default=4,
        help='CPU threads for FastTree / MAFFT (default: 4).')
    subtree.add_argument('--min-seqs', dest='min_seqs', type=int, default=3,
        help='Minimum sequences required to proceed with tree building (default: 3).')
    subtree.add_argument('--no-tree', dest='no_tree', action='store_true', default=False,
        help='Skip tree building; only write taxonomy TSV, sequence list, and iTOL colour files.')
    subtree.add_argument(
        '--group-phyla', dest='group_phyla', action='append', default=None, metavar='SPEC',
        help=(
            'Collapse multiple phyla into one colour in iTOL legends. Repeatable. '
            'Formats: "archaea", "bacteria", "Bacillota,Bacillota_I", "Firmicutes:Bacillota,Bacillota_I".'
        ),
    )
    subtree.add_argument('--functional', dest='functional', required=False, default=None,
        help=(
            'TSV file mapping sequence IDs to functional attributes. Header row required; first column = sequence ID; subsequent columns = functional attributes. '
            'Generates one iTOL file per column (DATASET_BINARY / DATASET_SIMPLEBAR / DATASET_COLORSTRIP).'
        ),
    )
    subtree.add_argument('--draft-rumen-functions', dest='draft_rumen_functions',
        action='store_true', default=False,
        help='Auto-generate rumen functional-group iTOL annotation from stored taxonomy.',
    )

    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == 'preload':
        cmd_preload(args)
    elif args.command == 'run':
        cmd_run(args)
    elif args.command == 'regen-itol':
        cmd_regen_itol(args)
    elif args.command == 'subtree':
        cmd_subtree(args)
    else:
        parser.print_help()


if __name__ == '__main__':
    main()

