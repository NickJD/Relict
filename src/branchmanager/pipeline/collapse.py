from __future__ import annotations

import logging
import shlex
import shutil
from dataclasses import dataclass
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Dict, List, Tuple

from branchmanager.utils.fasta import read_fasta, write_fasta
from branchmanager.utils.subprocess import run_cmd

logger = logging.getLogger(__name__)


@dataclass
class CollapseArtefacts:
    collapsed_path: Path
    map_path: Path
    members_path: Path
    collapsed_records: List[Tuple[str, str]]
    member_to_rep: Dict[str, str]


def collapse_fasta_within_taxa(
    taxa_groups,
    outdir: str,
    collapsed_fasta_name: str,
    map_name: str,
    members_name: str,
    threshold: float,
    threads: int = 1,
    log_prefix: str = '[COLLAPSE]',
    strict: bool = False,
):
    """Collapse records within taxonomic groups using vsearch cluster_fast.

    `taxa_groups` should be a mapping of taxon -> list[(header, sequence)].
    Output filenames are passed explicitly so Filing Cabinet/Performance Review can preserve their
    existing contracts.
    """
    outdir_p = Path(outdir)
    outdir_p.mkdir(parents=True, exist_ok=True)

    collapsed_records: List[Tuple[str, str]] = []
    cluster_map_lines: List[str] = []
    member_to_rep: Dict[str, str] = {}

    if shutil.which('vsearch') is None:
        if strict:
            raise RuntimeError('vsearch is required when sequence collapsing is requested')
        logger.warning('%s vsearch not found on PATH; skipping collapse step', log_prefix)
        for recs in taxa_groups.values():
            collapsed_records.extend(recs)
    else:
        for tax, recs in taxa_groups.items():
            if len(recs) <= 1:
                collapsed_records.extend(recs)
                continue
            member_map: Dict[str, List[str]] = {}
            with NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as tf:
                tmp_path = Path(tf.name)
                for h, seq in recs:
                    tf.write(f'>{h}\n{seq}\n')
            centroids = tmp_path.with_suffix('.centroids.fasta')
            ucfile = tmp_path.with_suffix('.uc')
            id_flag = threshold / 100.0
            try:
                cmd = (
                    f'vsearch --cluster_fast {shlex.quote(str(tmp_path))} --id {id_flag} '
                    f'--centroids {shlex.quote(str(centroids))} '
                    f'--uc {shlex.quote(str(ucfile))} --threads {int(threads)}'
                )
                logger.info('%s Running vsearch cluster for tax=%s (n=%d) threshold=%s', log_prefix, str(tax), len(recs), threshold)
                run_cmd(cmd)
                cent_recs = list(read_fasta(str(centroids))) if centroids.exists() else []
                if cent_recs:
                    try:
                        with open(ucfile) as uf:
                            for line in uf:
                                if not line.strip() or line.startswith('#'):
                                    continue
                                parts = line.split('\t')
                                # UC format columns (0-indexed):
                                #   0=type  1=cluster_id  8=seq_id  9=centroid_id
                                # Use parts[8] for the sequence ID (not parts[-1]
                                # which is the centroid ID for H rows, or '*' for S rows).
                                if parts[0] in ('H', 'S') and len(parts) > 8:
                                    seq_label = parts[8].strip()
                                    if seq_label:
                                        cluster_id = parts[1]
                                        member_map.setdefault(cluster_id, []).append(seq_label)
                    except Exception:
                        member_map = {}
                    for h, seq in cent_recs:
                        count = 1
                        found_cl = None
                        for cl, members in member_map.items():
                            if h in members:
                                count = len(members)
                                found_cl = cl
                                break
                        rep_h = h
                        collapsed_records.append((rep_h, seq))
                        cluster_map_lines.append(f'{rep_h}\t{tax}\t{count}\n')
                        if found_cl is not None:
                            cluster_members = member_map.get(found_cl, [])
                            for mem in cluster_members:
                                if mem != rep_h:  # only non-representative members
                                    member_to_rep[mem] = rep_h
                            if count > 1:
                                non_rep = [m for m in cluster_members if m != rep_h]
                                logger.info(
                                    '%s CLUSTER: representative=%s  size=%d  taxonomy=%s  '
                                    'collapsed_members=%s  '
                                    '→ representative enters the tree; members are EXCLUDED from tree '
                                    'but included in sequence_assessment.tsv',
                                    log_prefix, rep_h, count, str(tax),
                                    ';'.join(sorted(non_rep)) if non_rep else 'none',
                                )
                else:
                    if strict:
                        raise RuntimeError(
                            f'vsearch produced no centroid records for taxonomy group {tax!r}'
                        )
                    collapsed_records.extend(recs)
            except Exception as e:
                if strict:
                    raise RuntimeError(
                        f'vsearch clustering failed for taxonomy group {tax!r}: {e}'
                    ) from e
                logger.warning('%s vsearch clustering failed for tax=%s: %s', log_prefix, str(tax), e)
                collapsed_records.extend(recs)
            finally:
                try:
                    tmp_path.unlink()
                except Exception:
                    pass

    collapsed_path = outdir_p / collapsed_fasta_name
    map_path = outdir_p / map_name
    members_path = outdir_p / members_name

    write_fasta(collapsed_records, str(collapsed_path))
    with open(map_path, 'w') as m:
        m.write('rep_id\ttaxonomy\tcount\n')
        for line in cluster_map_lines:
            m.write(line)
    with open(members_path, 'w') as mm:
        mm.write('member\trep\n')
        for mem, rep in sorted(member_to_rep.items()):
            mm.write(f'{mem}\t{rep}\n')

    n_reps = len(collapsed_records)
    n_total = sum(len(recs) for recs in taxa_groups.values())
    n_clustered = len(member_to_rep)
    logger.info(
        '%s Collapse summary: %d input sequences → %d tree representatives '
        '(%d sequences clustered away, %.1f%% reduction). '
        'Cluster info written to %s and %s. '
        'All sequences (representatives AND their clustered members) appear in sequence_assessment.tsv.',
        log_prefix, n_total, n_reps, n_clustered,
        (100.0 * n_clustered / n_total) if n_total > 0 else 0.0,
        map_path.name, members_path.name,
    )

    return CollapseArtefacts(
        collapsed_path=collapsed_path,
        map_path=map_path,
        members_path=members_path,
        collapsed_records=collapsed_records,
        member_to_rep=member_to_rep,
    )
