"""Local phylogenetic-neighbourhood figures for assessed sequences."""

from __future__ import annotations

import csv
import logging
import math
from itertools import combinations_with_replacement
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from branchmanager.pipeline.cluster_report import _Node, _parse_node, _tokenise
from branchmanager.utils.fasta import read_fasta


logger = logging.getLogger(__name__)


def _walk(node: _Node) -> Iterable[_Node]:
    yield node
    for child in node.children:
        yield from _walk(child)


def _leaf_nodes(node: _Node) -> List[_Node]:
    if not node.children:
        return [node]
    leaves: List[_Node] = []
    for child in node.children:
        leaves.extend(_leaf_nodes(child))
    return leaves


def _parse_tree(path: str | Path) -> Tuple[_Node, Dict[int, _Node], Dict[str, _Node]]:
    text = Path(path).read_text().strip().rstrip(';')
    if not text:
        raise ValueError('tree is empty')
    root, _ = _parse_node(_tokenise(text), 0)
    parents: Dict[int, _Node] = {}
    leaves: Dict[str, _Node] = {}

    def visit(node: _Node) -> None:
        if not node.children and node.name:
            leaves[node.name.strip()] = node
        for child in node.children:
            parents[id(child)] = node
            visit(child)

    visit(root)
    return root, parents, leaves


def _leaf_count(node: _Node, cache: Dict[int, int]) -> int:
    key = id(node)
    if key not in cache:
        cache[key] = 1 if not node.children else sum(_leaf_count(c, cache) for c in node.children)
    return cache[key]


def _choose_context_root(
    leaf: _Node,
    parents: Dict[int, _Node],
    counts: Dict[int, int],
    min_context_leaves: int,
    max_context_leaves: int,
) -> _Node:
    current = parents.get(id(leaf), leaf)
    while _leaf_count(current, counts) < min_context_leaves:
        parent = parents.get(id(current))
        if parent is None or _leaf_count(parent, counts) > max_context_leaves:
            break
        current = parent
    return current


def _ancestor_chain(node: _Node, parents: Dict[int, _Node]) -> List[_Node]:
    chain = [node]
    while id(chain[-1]) in parents:
        chain.append(parents[id(chain[-1])])
    return chain


def _mrca(nodes: Iterable[_Node], parents: Dict[int, _Node]) -> _Node:
    selected = list(nodes)
    if not selected:
        raise ValueError('cannot find an MRCA for no nodes')
    common = {id(n) for n in _ancestor_chain(selected[0], parents)}
    for node in selected[1:]:
        common.intersection_update(id(n) for n in _ancestor_chain(node, parents))
    for node in _ancestor_chain(selected[0], parents):
        if id(node) in common:
            return node
    return selected[0]


def _group_targets(
    targets: List[dict],
    parents: Dict[int, _Node],
    counts: Dict[int, int],
    max_context_leaves: int,
) -> List[dict]:
    groups: List[dict] = []
    for target in targets:
        placed = False
        for group in groups:
            proposed = _mrca([group['root'], target['context_root']], parents)
            if id(group['root']) == id(target['context_root']) or _leaf_count(proposed, counts) <= max_context_leaves:
                group['targets'].append(target)
                group['root'] = proposed
                placed = True
                break
        if not placed:
            groups.append({'root': target['context_root'], 'targets': [target]})

    # A first-fit pass can leave mergeable sister groups. Merge to stability.
    changed = True
    while changed:
        changed = False
        for i in range(len(groups)):
            if changed:
                break
            for j in range(i + 1, len(groups)):
                proposed = _mrca([groups[i]['root'], groups[j]['root']], parents)
                if id(groups[i]['root']) == id(groups[j]['root']) or _leaf_count(proposed, counts) <= max_context_leaves:
                    groups[i]['targets'].extend(groups[j]['targets'])
                    groups[i]['root'] = proposed
                    groups.pop(j)
                    changed = True
                    break
    return groups


def _root_distances(root: _Node) -> Dict[int, float]:
    distances: Dict[int, float] = {}

    def visit(node: _Node, distance: float) -> None:
        distances[id(node)] = distance
        for child in node.children:
            visit(child, distance + max(0.0, float(child.branch_length or 0.0)))

    visit(root, 0.0)
    return distances


def _tree_distance(a: _Node, b: _Node, parents: Dict[int, _Node], distances: Dict[int, float]) -> float:
    ancestor_ids = {id(node) for node in _ancestor_chain(a, parents)}
    ancestor = next(node for node in _ancestor_chain(b, parents) if id(node) in ancestor_ids)
    return distances[id(a)] + distances[id(b)] - 2.0 * distances[id(ancestor)]


def pairwise_leaf_distances(tree_path: str | Path, leaf_ids: Optional[Iterable[str]] = None) -> Dict[tuple[str, str], float]:
    """Return patristic distances for requested tree leaves.

    Keys are stored in sorted order so callers can use ``tuple(sorted((a, b)))``.
    Missing leaves are omitted rather than assigned a synthetic distance.
    """
    root, parents, leaves = _parse_tree(tree_path)
    wanted = set(str(value) for value in leaf_ids) if leaf_ids is not None else set(leaves)
    names = sorted(name for name in leaves if name in wanted)
    root_distances = _root_distances(root)
    result: Dict[tuple[str, str], float] = {}
    for index, left in enumerate(names):
        result[(left, left)] = 0.0
        for right in names[index + 1:]:
            result[(left, right)] = _tree_distance(
                leaves[left], leaves[right], parents, root_distances,
            )
    return result


def _prune_copy(node: _Node, keep_names: set[str]) -> Optional[_Node]:
    if not node.children:
        if node.name.strip() not in keep_names:
            return None
        return _Node(node.name, node.branch_length)
    children = [copy for child in node.children if (copy := _prune_copy(child, keep_names)) is not None]
    if not children:
        return None
    if len(children) == 1:
        child = children[0]
        child.branch_length = float(child.branch_length or 0.0) + float(node.branch_length or 0.0)
        return child
    copy = _Node(node.name, node.branch_length)
    copy.children = children
    return copy


def _bounded_render_root(
    root: _Node,
    targets: List[dict],
    parents: Dict[int, _Node],
    distances: Dict[int, float],
    max_context_leaves: int,
    forced_leaf_names: Optional[set[str]] = None,
) -> _Node:
    leaves = _leaf_nodes(root)
    forced_leaf_names = set(forced_leaf_names or set())
    if len(leaves) <= max_context_leaves and forced_leaf_names.issubset(
        {leaf.name.strip() for leaf in leaves}
    ):
        return root
    target_leaves = [target['leaf'] for target in targets]
    ranked = sorted(
        leaves,
        key=lambda leaf: (
            min(_tree_distance(leaf, target, parents, distances) for target in target_leaves),
            leaf.name,
        ),
    )
    keep = {target.name.strip() for target in target_leaves}
    keep.update(forced_leaf_names)
    for leaf in ranked:
        if len(keep) >= max_context_leaves:
            break
        keep.add(leaf.name.strip())
    pruned = _prune_copy(root, keep)
    return pruned or root


def _terminal_taxon(taxonomy: object) -> str:
    parts = [p.strip() for p in str(taxonomy or '').split(';') if p.strip()]
    informative = [p for p in parts if not p.endswith('__') and p.lower() not in ('na', 'none')]
    return informative[-1] if informative else ''


def _baseline_dataset_names(rows: List[dict]) -> set[str]:
    names: set[str] = set()
    for row in rows:
        source = str(row.get('density_source') or '')
        if not source.startswith('baseline:'):
            continue
        for name in source.split(':', 1)[1].split(','):
            if name.strip():
                names.add(name.strip())
    return names


def _load_metadata(db, leaf_names: List[str], assessment_rows: List[dict]) -> Dict[str, dict]:
    metadata: Dict[str, dict] = {name: {'dataset': '', 'taxonomy': ''} for name in leaf_names}
    if db is not None and leaf_names:
        try:
            with db.connect() as conn:
                cur = conn.cursor()
                placeholders = ','.join('?' for _ in leaf_names)
                cur.execute(
                    "SELECT s.id, s.dataset, t.taxonomy "
                    "FROM sequences s LEFT JOIN taxonomy t ON s.id = t.id "
                    f"WHERE s.id IN ({placeholders})",
                    tuple(leaf_names),
                )
                for sid, dataset, taxonomy in cur.fetchall():
                    entry = metadata.setdefault(str(sid), {'dataset': '', 'taxonomy': ''})
                    entry['dataset'] = entry.get('dataset') or dataset or ''
                    entry['taxonomy'] = entry.get('taxonomy') or taxonomy or ''
        except Exception as exc:
            logger.warning('[NEIGHBOURHOOD] Could not load tree metadata: %s', exc)
        try:
            sequencing = db.get_sequencing_metadata_for_ids(leaf_names)
            for sid, values in sequencing.items():
                metadata.setdefault(str(sid), {}).update({
                    'partner_id': values.get('partner_id', ''),
                    'selected_for_wgs': bool(values.get('selected_for_wgs')),
                    'selected_for_sequencing': bool(values.get('selected_for_sequencing')),
                })
        except Exception as exc:
            logger.warning('[NEIGHBOURHOOD] Could not load genome-availability metadata: %s', exc)

    for row in assessment_rows:
        sid = str(row.get('id') or '')
        if not sid:
            continue
        entry = metadata.setdefault(sid, {})
        entry['taxonomy'] = entry.get('taxonomy') or row.get('taxonomy') or ''
        entry['partner_id'] = entry.get('partner_id') or row.get('partner_id') or ''
        if str(row.get('already_sequenced') or '').lower() == 'true':
            entry['selected_for_wgs'] = True
        if str(row.get('selected_for_genome_sequencing') or '').lower() == 'true':
            entry['selected_for_sequencing'] = True
        entry['sequencing_set_id'] = row.get('sequencing_set_id') or ''
        entry['sequencing_set_role'] = row.get('sequencing_set_role') or ''
        entry['sequencing_set_rank'] = row.get('sequencing_set_rank') or ''
        entry['sequencing_set_badge'] = row.get('sequencing_set_badge') or ''
    return metadata


def _load_alignment_sequences(path: Optional[str | Path]) -> Dict[str, str]:
    if not path or not Path(path).exists():
        return {}
    try:
        return {str(header).strip(): str(sequence).upper().replace('U', 'T') for header, sequence in read_fasta(str(path))}
    except Exception as exc:
        logger.warning('[NEIGHBOURHOOD] Could not read alignment %s for pident: %s', path, exc)
        return {}


def _msa_pident(left: str, right: str) -> tuple[Optional[float], int]:
    """Return identity across jointly unambiguous A/C/G/T MSA columns."""
    if not left or not right:
        return None, 0
    width = max(len(left), len(right))
    left = left.ljust(width, '-')
    right = right.ljust(width, '-')
    compared = 0
    matches = 0
    for a, b in zip(left, right):
        if a not in {'A', 'C', 'G', 'T'} or b not in {'A', 'C', 'G', 'T'}:
            continue
        compared += 1
        if a == b:
            matches += 1
    return ((100.0 * matches / compared) if compared else None), compared


def _selection_badge(metadata: dict) -> str:
    explicit = str(metadata.get('sequencing_set_badge') or '').strip()
    if explicit:
        return explicit
    role = str(metadata.get('sequencing_set_role') or '').strip().upper()
    rank = str(metadata.get('sequencing_set_rank') or '').strip()
    if role == 'PRIMARY':
        return f'P{rank}' if rank not in ('', 'NA', 'None') else 'PRIMARY'
    if role == 'BACKUP':
        return f'B{rank}' if rank not in ('', 'NA', 'None') else 'BACKUP'
    if role == 'DIVERSITY_CANDIDATE':
        return f'D{rank}' if rank not in ('', 'NA', 'None') else 'DIVERSITY'
    return {
        'COMMITTED': 'SELECTED',
        'SEQUENCED': 'SEQUENCED',
        'ALTERNATE': 'ALT',
        'REVIEW_EVIDENCE': 'REVIEW',
        'PANGENOME_BOUNDARY_REVIEW': 'BOUNDARY REVIEW',
        'BASELINE_REDUNDANT': 'BASELINE REDUNDANT',
        'TARGET_MET': 'TARGET MET',
    }.get(role, '')


def _star_points(cx: float, cy: float, outer: float = 7.0, inner: float = 3.2) -> List[tuple[float, float]]:
    points = []
    for index in range(10):
        angle = -math.pi / 2.0 + index * math.pi / 5.0
        radius = outer if index % 2 == 0 else inner
        points.append((cx + math.cos(angle) * radius, cy + math.sin(angle) * radius))
    return points


def _identity_anchor(targets: List[dict], aligned_sequences: Dict[str, str]) -> Optional[str]:
    role_order = {
        'PRIMARY': 0, 'BACKUP': 1, 'SEQUENCED': 2, 'COMMITTED': 2,
        'DIVERSITY_CANDIDATE': 3, 'ALTERNATE': 4,
    }

    def key(target: dict) -> tuple:
        row = target['row']
        role = str(row.get('sequencing_set_role') or '').upper()
        try:
            rank = int(float(row.get('sequencing_set_rank')))
        except (TypeError, ValueError):
            rank = 999
        return role_order.get(role, 99), rank, str(row.get('id') or '')

    for target in sorted(targets, key=key):
        leaf_id = str(target.get('tree_leaf_id') or '')
        if leaf_id in aligned_sequences:
            return leaf_id
    return None


def _baseline_context_for_group(group: dict, leaf_map: Dict[str, _Node]) -> Dict[str, dict]:
    context: Dict[str, dict] = {}
    for target in group['targets']:
        row = target['row']
        hit = str(row.get('nearest_hit') or '').strip()
        if not hit or hit in ('NA', 'None') or hit not in leaf_map:
            continue
        try:
            identity = float(row.get('nearest_identity'))
        except (TypeError, ValueError):
            identity = None
        entry = context.setdefault(hit, {'identity': identity, 'queries': []})
        if identity is not None and (entry.get('identity') is None or identity > entry['identity']):
            entry['identity'] = identity
        entry['queries'].append(str(row.get('id') or ''))
    return context


def _write_pairwise_pident_table(
    path: Path,
    leaf_names: List[str],
    aligned_sequences: Dict[str, str],
) -> bool:
    available = [name for name in leaf_names if name in aligned_sequences]
    if not available:
        return False
    with open(path, 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow([
            'SequenceA', 'SequenceB', 'MSAPercentIdentity',
            'ComparableACGTColumns', 'IdentityDefinition',
        ])
        for left, right in combinations_with_replacement(available, 2):
            identity, compared = _msa_pident(aligned_sequences[left], aligned_sequences[right])
            writer.writerow([
                left,
                right,
                f'{identity:.2f}' if identity is not None else 'NA',
                compared,
                'identical_ACGT_bases/jointly_unambiguous_ACGT_MSA_columns',
            ])
    return True


def _nice_scale(max_distance: float) -> float:
    target = max_distance / 5.0
    if target <= 0:
        return 0.0
    power = 10 ** math.floor(math.log10(target))
    for multiplier in (1, 2, 5, 10):
        value = multiplier * power
        if value >= target:
            return value
    return target


def _figure_geometry(
    root: _Node,
    targets: List[dict],
    metadata: Dict[str, dict],
    baseline_datasets: set[str],
    aligned_sequences: Optional[Dict[str, str]] = None,
    identity_anchor: Optional[str] = None,
    baseline_context: Optional[Dict[str, dict]] = None,
) -> dict:
    aligned_sequences = aligned_sequences or {}
    baseline_context = baseline_context or {}
    leaves = _leaf_nodes(root)
    assessed_by_leaf: Dict[str, List[str]] = {}
    for target in targets:
        assessed_by_leaf.setdefault(target['tree_leaf_id'], []).append(str(target['row'].get('id') or ''))

    y_start = 156.0
    row_height = 32.0
    leaf_y = {id(leaf): y_start + idx * row_height for idx, leaf in enumerate(leaves)}
    node_y: Dict[int, float] = {}

    def assign_y(node: _Node) -> float:
        if not node.children:
            node_y[id(node)] = leaf_y[id(node)]
        else:
            child_y = [assign_y(c) for c in node.children]
            node_y[id(node)] = sum(child_y) / len(child_y)
        return node_y[id(node)]

    assign_y(root)
    distances: Dict[int, float] = {}
    levels: Dict[int, int] = {}

    def assign_distance(node: _Node, distance: float, level: int) -> None:
        distances[id(node)] = distance
        levels[id(node)] = level
        for child in node.children:
            assign_distance(child, distance + max(0.0, float(child.branch_length or 0.0)), level + 1)

    assign_distance(root, 0.0, 0)
    max_distance = max((distances[id(leaf)] for leaf in leaves), default=0.0)
    max_level = max((levels[id(leaf)] for leaf in leaves), default=1) or 1
    plot_left = 54.0
    plot_width = 560.0
    node_x = {
        id(node): plot_left + (
            (distances[id(node)] / max_distance) if max_distance > 0
            else (levels[id(node)] / max_level)
        ) * plot_width
        for node in _walk(root)
    }
    label_x = plot_left + plot_width + 28.0
    display_labels = []
    for leaf in leaves:
        sid = leaf.name.strip()
        meta = metadata.get(sid, {})
        badge = _selection_badge(meta)
        context = ' | '.join(
            value for value in (
                str(meta.get('dataset') or ''),
                _terminal_taxon(meta.get('taxonomy')),
            ) if value
        )
        assessed = assessed_by_leaf.get(sid, [])
        assessed_note = f" | assessed: {', '.join(assessed)}" if assessed and assessed != [sid] else ''
        baseline_note = ''
        if sid in baseline_context:
            identity = baseline_context[sid].get('identity')
            baseline_note = (
                f' | nearest baseline context; max {identity:.2f}% vsearch'
                if identity is not None else ' | nearest baseline context'
            )
        pident_note = ''
        if identity_anchor and sid in aligned_sequences and identity_anchor in aligned_sequences:
            identity, compared = _msa_pident(aligned_sequences[sid], aligned_sequences[identity_anchor])
            if identity is not None:
                pident_note = f' | MSA pident {identity:.2f}% (n={compared})'
        prefix = f'[{badge}] ' if badge else ''
        display_labels.append(
            f'{prefix}{sid}{assessed_note}{(" | " + context) if context else ""}'
            f'{baseline_note}{pident_note}'
        )
    longest = max((len(label) for label in display_labels), default=20)
    width = max(1050.0, label_x + longest * 7.8 + 48.0)
    height = max(260.0, y_start + max(0, len(leaves) - 1) * row_height + 58.0)
    assessed_ids = sorted({sid for values in assessed_by_leaf.values() for sid in values})
    selected_count = sum(1 for leaf in leaves if metadata.get(leaf.name.strip(), {}).get('selected_for_wgs'))
    pending_count = sum(
        1 for leaf in leaves
        if metadata.get(leaf.name.strip(), {}).get('selected_for_sequencing')
        and not metadata.get(leaf.name.strip(), {}).get('selected_for_wgs')
    )
    baseline_count = sum(1 for leaf in leaves if metadata.get(leaf.name.strip(), {}).get('dataset') in baseline_datasets)
    return {
        'leaves': leaves,
        'assessed_by_leaf': assessed_by_leaf,
        'node_y': node_y,
        'node_x': node_x,
        'plot_left': plot_left,
        'plot_width': plot_width,
        'label_x': label_x,
        'display_labels': display_labels,
        'width': width,
        'height': height,
        'max_distance': max_distance,
        'assessed_ids': assessed_ids,
        'selected_count': selected_count,
        'pending_count': pending_count,
        'baseline_count': baseline_count,
        'identity_anchor': identity_anchor or '',
        'nearest_baseline_ids': sorted(baseline_context),
    }


def _load_png_font(size: int, *, bold: bool = False):
    try:
        from PIL import ImageFont
    except ImportError as exc:
        raise RuntimeError(
            'PNG neighbourhood figures require Pillow. Install it with `conda install pillow`.'
        ) from exc
    candidates = [
        '/System/Library/Fonts/Supplemental/Arial Bold.ttf' if bold else '/System/Library/Fonts/Supplemental/Arial.ttf',
        '/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf' if bold else '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',
        'DejaVuSans-Bold.ttf' if bold else 'DejaVuSans.ttf',
    ]
    for candidate in candidates:
        try:
            return ImageFont.truetype(candidate, size=size)
        except (OSError, ValueError):
            continue
    return ImageFont.load_default()


def _write_group_png(
    path: Path,
    group_number: int,
    root: _Node,
    targets: List[dict],
    metadata: Dict[str, dict],
    baseline_datasets: set[str],
    aligned_sequences: Optional[Dict[str, str]] = None,
    identity_anchor: Optional[str] = None,
    baseline_context: Optional[Dict[str, dict]] = None,
) -> dict:
    try:
        from PIL import Image, ImageDraw
    except ImportError as exc:
        raise RuntimeError(
            'PNG neighbourhood figures require Pillow. Install it with `conda install pillow`.'
        ) from exc

    geometry = _figure_geometry(
        root, targets, metadata, baseline_datasets,
        aligned_sequences=aligned_sequences,
        identity_anchor=identity_anchor,
        baseline_context=baseline_context,
    )
    leaves = geometry['leaves']
    assessed_by_leaf = geometry['assessed_by_leaf']
    node_y = geometry['node_y']
    node_x = geometry['node_x']
    plot_left = geometry['plot_left']
    plot_width = geometry['plot_width']
    label_x = geometry['label_x']
    display_labels = geometry['display_labels']
    width = geometry['width']
    height = geometry['height']
    max_distance = geometry['max_distance']
    assessed_ids = geometry['assessed_ids']
    selected_count = geometry['selected_count']
    pending_count = geometry['pending_count']
    baseline_count = geometry['baseline_count']
    identity_anchor = geometry['identity_anchor']
    nearest_baseline_ids = geometry['nearest_baseline_ids']

    resolution = 2
    def px(value):
        return int(round(float(value) * resolution))

    image = Image.new('RGB', (px(width), px(height)), '#ffffff')
    draw = ImageDraw.Draw(image)
    title_font = _load_png_font(px(18), bold=True)
    body_font = _load_png_font(px(13))
    label_font = _load_png_font(px(12))
    label_bold_font = _load_png_font(px(12), bold=True)
    small_font = _load_png_font(px(10))
    node_font = _load_png_font(px(9))

    draw.text((px(28), px(13)), f'Local phylogenetic neighbourhood {group_number:03d}', fill='#17202a', font=title_font)
    draw.text(
        (px(28), px(39)),
        f'{len(assessed_ids)} assessed isolate(s); {len(leaves)} tree leaves; '
        f'{baseline_count} baseline; {selected_count} sequenced; {pending_count} selected pending',
        fill='#4d5966',
        font=body_font,
    )
    draw.ellipse((px(28), px(76), px(40), px(88)), fill='#1565c0')
    draw.text((px(46), px(73)), 'assessed isolate', fill='#25313d', font=label_font)
    draw.rectangle((px(174), px(76), px(186), px(88)), fill='#2e7d32')
    draw.text((px(193), px(73)), 'baseline dataset', fill='#25313d', font=label_font)
    draw.polygon([(px(330), px(75)), (px(337), px(82)), (px(330), px(89)), (px(323), px(82))], fill='#ef6c00')
    draw.text((px(345), px(73)), 'already sequenced', fill='#25313d', font=label_font)
    draw.rectangle((px(468), px(76), px(480), px(88)), fill='#ffffff', outline='#b7791f', width=px(2))
    draw.text((px(487), px(73)), 'selected, genome pending', fill='#25313d', font=label_font)
    draw.polygon([(px(x), px(y)) for x, y in _star_points(680, 82)], fill='#8e24aa', outline='#5f146f')
    draw.text((px(693), px(73)), 'sequence next', fill='#25313d', font=label_font)
    draw.polygon([(px(x), px(y)) for x, y in _star_points(830, 82)], fill='#ffffff', outline='#8e24aa')
    draw.text((px(843), px(73)), 'backup', fill='#25313d', font=label_font)
    draw.text(
        (px(28), px(101)),
        'Selection labels: [P#] primary; [B#] backup; [D#] post-target diversity; boundary review/exclusion shown in labels',
        fill='#4d5966', font=small_font,
    )
    identity_text = (
        f'MSA pident labels are relative to {identity_anchor}; internal branches remain substitutions/site.'
        if identity_anchor else
        'Internal branches are modelled substitutions/site; no MSA pident anchor was available.'
    )
    draw.text((px(28), px(120)), identity_text, fill='#4d5966', font=small_font)

    for node in _walk(root):
        if not node.children:
            continue
        child_ys = [node_y[id(child)] for child in node.children]
        x = node_x[id(node)]
        draw.line((px(x), px(min(child_ys)), px(x), px(max(child_ys))), fill='#66717d', width=px(1.4))
        for child in node.children:
            cy = node_y[id(child)]
            draw.line((px(x), px(cy), px(node_x[id(child)]), px(cy)), fill='#66717d', width=px(1.4))
        if node.name and node.name.startswith('NODE'):
            draw.text((px(x + 4), px(node_y[id(node)] - 14)), node.name, fill='#7c8792', font=node_font)

    for leaf, label in zip(leaves, display_labels):
        sid = leaf.name.strip()
        y = node_y[id(leaf)]
        x = node_x[id(leaf)]
        meta = metadata.get(sid, {})
        assessed = sid in assessed_by_leaf
        baseline = meta.get('dataset') in baseline_datasets
        selected = bool(meta.get('selected_for_wgs'))
        pending = bool(meta.get('selected_for_sequencing')) and not selected
        selection_role = str(meta.get('sequencing_set_role') or '').upper()
        if assessed:
            draw.ellipse((px(x - 6), px(y - 6), px(x + 6), px(y + 6)), fill='#1565c0', outline='#0d3f78', width=px(1))
        elif baseline:
            draw.rectangle((px(x - 5), px(y - 5), px(x + 5), px(y + 5)), fill='#2e7d32')
        else:
            draw.ellipse((px(x - 4), px(y - 4), px(x + 4), px(y + 4)), fill='#8c98a4')
        if selected:
            sx = x + 12.0
            draw.polygon(
                [(px(sx), px(y - 6)), (px(sx + 6), px(y)), (px(sx), px(y + 6)), (px(sx - 6), px(y))],
                fill='#ef6c00',
                outline='#9a3f00',
            )
        elif pending:
            draw.rectangle(
                (px(x + 7), px(y - 5), px(x + 17), px(y + 5)),
                fill='#ffffff', outline='#b7791f', width=px(2),
            )
        elif selection_role in ('PRIMARY', 'DIVERSITY_CANDIDATE'):
            draw.polygon(
                [(px(px_x), px(px_y)) for px_x, px_y in _star_points(x + 13.0, y)],
                fill='#8e24aa', outline='#5f146f',
            )
        elif selection_role == 'BACKUP':
            draw.polygon(
                [(px(px_x), px(px_y)) for px_x, px_y in _star_points(x + 13.0, y)],
                fill='#ffffff', outline='#8e24aa',
            )
        draw.text(
            (px(label_x), px(y - 8)),
            label,
            fill='#0d3f78' if assessed else '#17202a',
            font=label_bold_font if assessed else label_font,
        )

    if max_distance > 0:
        scale_value = _nice_scale(max_distance)
        scale_width = (scale_value / max_distance) * plot_width
        if 18 <= scale_width <= plot_width:
            sy = height - 24.0
            draw.line((px(plot_left), px(sy), px(plot_left + scale_width), px(sy)), fill='#25313d', width=px(2))
            draw.line((px(plot_left), px(sy - 4), px(plot_left), px(sy + 4)), fill='#25313d', width=px(1))
            draw.line((px(plot_left + scale_width), px(sy - 4), px(plot_left + scale_width), px(sy + 4)), fill='#25313d', width=px(1))
            draw.text((px(plot_left), px(sy - 18)), f'{scale_value:g} substitutions/site', fill='#4d5966', font=small_font)
    else:
        draw.text((px(plot_left), px(height - 28)), 'Cladogram: branch lengths unavailable', fill='#4d5966', font=small_font)

    image.save(path, format='PNG', optimize=True, dpi=(144, 144))
    return {
        'leaf_count': len(leaves),
        'assessed_ids': assessed_ids,
        'baseline_count': baseline_count,
        'selected_count': selected_count,
        'pending_count': pending_count,
        'identity_anchor': identity_anchor,
        'nearest_baseline_ids': nearest_baseline_ids,
    }


def generate_local_neighbourhood_visuals(
    tree_path: str | Path,
    assessment_rows: List[dict],
    db,
    outdir: str | Path,
    *,
    alignment_path: Optional[str | Path] = None,
    min_context_leaves: int = 8,
    max_context_leaves: int = 30,
    image_format: str = 'png',
) -> dict:
    """Write grouped local-subtree images and link each assessment row to one.

    A local context expands from each assessed leaf until it contains at least
    ``min_context_leaves`` without exceeding ``max_context_leaves``. Contexts
    whose MRCA still fits within that limit are combined, so nearby assessed
    isolates share one figure rather than producing duplicate plots.
    """
    image_format = str(image_format or 'png').strip().lower()
    if image_format != 'png':
        raise ValueError(f'Unsupported neighbourhood image format: {image_format}; only png is supported')
    output_dir = Path(outdir)
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / 'neighbourhood_manifest.tsv'
    for row in assessment_rows:
        row['local_neighbourhood_figure'] = 'NA'
        row['local_pairwise_pident_table'] = 'NA'
        row['tree_context_leaf_count'] = 'NA'
        row['assessed_sequences_in_tree_context'] = 'NA'

    try:
        root, parents, leaf_map = _parse_tree(tree_path)
    except Exception as exc:
        logger.warning('[NEIGHBOURHOOD] Could not parse %s: %s', tree_path, exc)
        with open(manifest_path, 'w') as handle:
            handle.write(
                'SequenceID\tTreeLeafID\tFigure\tPairwisePidentTable\tIdentityAnchor\t'
                'NearestBaselineHitsShown\tAssessedSequencesInFigure\tTreeLeavesShown\t'
                'BaselineLeavesShown\tSequencedGenomeLeavesShown\tStatus\n'
            )
            for row in assessment_rows:
                handle.write(
                    f"{row.get('id', 'NA')}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\ttree_unavailable\n"
                )
        return {'manifest': str(manifest_path), 'figures': [], 'resolved': 0, 'unresolved': len(assessment_rows)}

    counts: Dict[int, int] = {}
    distances = _root_distances(root)
    global_order = {leaf.name.strip(): idx for idx, leaf in enumerate(_leaf_nodes(root))}
    targets: List[dict] = []
    unresolved: List[dict] = []
    for row in assessment_rows:
        sid = str(row.get('id') or '')
        representative = str(row.get('cluster_representative') or '')
        tree_leaf_id = sid if sid in leaf_map else representative
        if tree_leaf_id not in leaf_map:
            unresolved.append({'row': row, 'tree_leaf_id': 'NA'})
            continue
        leaf = leaf_map[tree_leaf_id]
        context_root = _choose_context_root(
            leaf,
            parents,
            counts,
            max(2, int(min_context_leaves)),
            max(2, int(max_context_leaves)),
        )
        targets.append({
            'row': row,
            'leaf': leaf,
            'tree_leaf_id': tree_leaf_id,
            'context_root': context_root,
        })

    targets.sort(key=lambda target: global_order.get(target['tree_leaf_id'], 10**9))
    groups = _group_targets(targets, parents, counts, max(2, int(max_context_leaves)))
    groups.sort(key=lambda group: min(global_order.get(leaf.name.strip(), 10**9) for leaf in _leaf_nodes(group['root'])))
    for group in groups:
        baseline_context = _baseline_context_for_group(group, leaf_map)
        forced_nodes = [leaf_map[sid] for sid in baseline_context if sid in leaf_map]
        group['baseline_context'] = baseline_context
        group['render_source_root'] = (
            _mrca([group['root'], *forced_nodes], parents) if forced_nodes else group['root']
        )
    all_leaf_names = [
        leaf.name.strip()
        for group in groups
        for leaf in _leaf_nodes(group['render_source_root'])
    ]
    metadata = _load_metadata(db, sorted(set(all_leaf_names)), assessment_rows)
    baseline_datasets = _baseline_dataset_names(assessment_rows)
    if alignment_path is None:
        tree_parent = Path(tree_path).parent
        candidates = [
            tree_parent / 'current_alignment.fasta',
            tree_parent.parent / 'current_alignment.fasta',
            tree_parent.parent / 'tree' / 'current_alignment.fasta',
        ]
        alignment_path = next((path for path in candidates if path.exists()), None)
    aligned_sequences = _load_alignment_sequences(alignment_path)

    manifest_rows: List[dict] = []
    figures: List[str] = []
    for index, group in enumerate(groups, start=1):
        filename = f'clade_{index:03d}.{image_format}'
        figure_path = output_dir / filename
        render_root = _bounded_render_root(
            group['render_source_root'],
            group['targets'],
            parents,
            distances,
            max(2, int(max_context_leaves)),
            forced_leaf_names=set(group['baseline_context']),
        )
        anchor = _identity_anchor(group['targets'], aligned_sequences)
        figure_meta = _write_group_png(
            figure_path,
            index,
            render_root,
            group['targets'],
            metadata,
            baseline_datasets,
            aligned_sequences=aligned_sequences,
            identity_anchor=anchor,
            baseline_context=group['baseline_context'],
        )
        figures.append(str(figure_path))
        pairwise_filename = f'clade_{index:03d}_pairwise_pident.tsv'
        pairwise_path = output_dir / pairwise_filename
        pairwise_written = _write_pairwise_pident_table(
            pairwise_path,
            [leaf.name.strip() for leaf in _leaf_nodes(render_root)],
            aligned_sequences,
        )
        relative_pairwise = f'neighbourhoods/{pairwise_filename}' if pairwise_written else 'NA'
        assessed_text = ';'.join(figure_meta['assessed_ids'])
        relative_figure = f'neighbourhoods/{filename}'
        baseline_text = ';'.join(figure_meta['nearest_baseline_ids']) or 'None'
        for target in group['targets']:
            row = target['row']
            row['local_neighbourhood_figure'] = relative_figure
            row['local_pairwise_pident_table'] = relative_pairwise
            row['tree_context_leaf_count'] = str(figure_meta['leaf_count'])
            row['assessed_sequences_in_tree_context'] = str(len(figure_meta['assessed_ids']))
            manifest_rows.append({
                'SequenceID': row.get('id', 'NA'),
                'TreeLeafID': target['tree_leaf_id'],
                'Figure': relative_figure,
                'PairwisePidentTable': relative_pairwise,
                'IdentityAnchor': figure_meta['identity_anchor'] or 'NA',
                'NearestBaselineHitsShown': baseline_text,
                'AssessedSequencesInFigure': assessed_text,
                'TreeLeavesShown': figure_meta['leaf_count'],
                'BaselineLeavesShown': figure_meta['baseline_count'],
                'SequencedGenomeLeavesShown': figure_meta['selected_count'],
                'Status': 'rendered',
            })

    for target in unresolved:
        manifest_rows.append({
            'SequenceID': target['row'].get('id', 'NA'),
            'TreeLeafID': 'NA',
            'Figure': 'NA',
            'PairwisePidentTable': 'NA',
            'IdentityAnchor': 'NA',
            'NearestBaselineHitsShown': 'NA',
            'AssessedSequencesInFigure': 'NA',
            'TreeLeavesShown': 'NA',
            'BaselineLeavesShown': 'NA',
            'SequencedGenomeLeavesShown': 'NA',
            'Status': 'not_present_in_tree',
        })

    fields = [
        'SequenceID', 'TreeLeafID', 'Figure', 'PairwisePidentTable',
        'IdentityAnchor', 'NearestBaselineHitsShown', 'AssessedSequencesInFigure',
        'TreeLeavesShown', 'BaselineLeavesShown', 'SequencedGenomeLeavesShown', 'Status',
    ]
    with open(manifest_path, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(sorted(manifest_rows, key=lambda row: str(row['SequenceID'])))

    logger.info(
        '[NEIGHBOURHOOD] Wrote %d grouped local-clade figure(s) for %d/%d assessed sequences to %s',
        len(figures), len(targets), len(assessment_rows), output_dir,
    )
    return {
        'manifest': str(manifest_path),
        'figures': figures,
        'resolved': len(targets),
        'unresolved': len(unresolved),
    }
