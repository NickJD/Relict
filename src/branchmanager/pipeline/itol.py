from pathlib import Path
import hashlib
from branchmanager.taxonomy import parse_taxon_string as _shared_parse_taxon_string


def _hash_to_hue(s: str) -> float:
    h = hashlib.md5(s.encode('utf-8')).hexdigest()
    val = int(h[:8], 16)
    return (val % 360) / 360.0


def _hsv_to_hex(h, s=0.65, v=0.95):
    import colorsys
    r, g, b = colorsys.hsv_to_rgb(h, s, v)
    return '#{:02x}{:02x}{:02x}'.format(int(r*255), int(g*255), int(b*255))


def _name_to_colour(name: str) -> str:
    h = _hash_to_hue(name)
    return _hsv_to_hex(h)


def _name_to_dataset_colour(name: str, bins: int = 8) -> str:
    """Generate a visually distinct colour for dataset-level labels.

    Hue quantisation keeps similarly hashed names visually distinct.
    """
    base_h = _hash_to_hue(name)
    idx = int(base_h * bins) % bins
    h = (idx + 0.5) / bins
    h_frac = _hash_to_hue(name + '_jitter')
    jitter = (h_frac - 0.5) * (0.5 / bins)
    h = (h + jitter) % 1.0
    return _hsv_to_hex(h, s=0.85, v=0.92)


def _identity_to_colour(identity: float, vmin: float = 0.0, vmax: float = 1.0) -> str:
    """Map an identity score to a red-to-green colour gradient.

    Values may be proportions or percentages; missing values use neutral grey.
    """
    if identity is None:
        return '#cccccc'
    try:
        val = float(identity)
    except Exception:
        return '#cccccc'
    if val > 1.0 and val <= 100.0:
        val = val / 100.0
    if val is None or val != val:  # NaN
        return '#cccccc'
    val = max(vmin, min(vmax, val))
    if vmax - vmin > 0:
        norm = (val - vmin) / (vmax - vmin)
    else:
        norm = 0.0
    h = norm * 0.33
    return _hsv_to_hex(h, s=0.95, v=0.95)


def _novelty_colour_for_pct(pct, min_pct=90, max_pct=100):
    """Return a distinct colour for a novelty category.

    pct may be '<90' or an int between min_pct and max_pct inclusive.
    Colours are evenly spaced between red (low) and green (high).
    """
    if pct == '<90' or pct is None:
        return _hsv_to_hex(0.0, s=0.9, v=0.9)
    try:
        p = int(pct)
    except Exception:
        return _hsv_to_hex(0.0, s=0.9, v=0.9)
    if p < min_pct:
        return _hsv_to_hex(0.0, s=0.9, v=0.9)
    if p > max_pct:
        p = max_pct
    norm = (p - min_pct) / float(max_pct - min_pct) if max_pct > min_pct else 1.0
    h = norm * 0.33
    return _hsv_to_hex(h, s=0.9, v=0.9)


def _name_to_colour_by_rank(name: str, rank: str = None) -> str:
    """Generate a deterministic colour adjusted by taxonomic rank.

    rank can be 'phylum', 'family', 'genus' or None. Different ranks receive
    a hue offset and slightly different saturation/value to increase contrast.
    """
    h = _hash_to_hue(name)
    rank_offsets = {
        'phylum': 0.0,
        'family': 0.33,
        'genus': 0.66,
    }
    offset = rank_offsets.get(rank, 0.0)
    h = (h + offset) % 1.0

    sv_map = {
        'phylum': (0.85, 0.97),
        'family': (0.65, 0.90),
        'genus': (0.45, 0.80),
    }
    s, v = sv_map.get(rank, (0.65, 0.90))
    return _hsv_to_hex(h, s=s, v=v)


def parse_taxon_string(taxon: str):
    return _shared_parse_taxon_string(taxon)


def _normalise_taxon_name(name: str, preserve_subgroup: bool = False) -> str:
    """Normalise taxon names produced by gg2/reference sources.

    Rules applied:
    - strip whitespace and surrounding quotes
    - collapse internal whitespace
    - drop terminal numeric accession-like suffixes (e.g. _368345)
    - keep short alphabetic subgroup suffixes (e.g. Bacillota_A) but remove
      trailing numeric ids that create many redundant variants
    """
    if not name:
        return 'unknown'
    try:
        s = str(name).strip()
    except Exception:
        return 'unknown'
    if not s:
        return 'unknown'
    if (s.startswith('"') and s.endswith('"')) or (s.startswith("'") and s.endswith("'")):
        s = s[1:-1].strip()
    s = ' '.join(s.split())
    import re
    parts = re.split(r'[\s_]+', s)
    if not parts:
        return s
    while parts and parts[-1].isdigit():
        parts.pop()
    # Preserve phylum subgroup suffixes such as Bacillota_A when requested.
    if not preserve_subgroup:
        while len(parts) > 1 and re.fullmatch(r'[A-Za-z]', parts[-1]):
            parts.pop()
    norm = '_'.join([p for p in parts if p])
    return norm or 'unknown'


def _parse_phylum_groups(specs, all_phyla, ph_to_domain):
    """Parse ``--group-phyla`` spec strings into a ``{phylum: group_label}`` map.

    Each *spec* in *specs* may be one of:

    ``archaea``
        Group **all** phyla whose domain contains "archaea" under the single
        label ``"Archaea"``.

    ``bacteria``
        Group all bacterial phyla under ``"Bacteria"``.

    ``"Bacillota,Bacillota_I,Bacillota_A"``
        Group these phyla together; the label defaults to the first name in the
        list (``"Bacillota"``).

    ``"Firmicutes:Bacillota,Bacillota_I"``
        Explicitly name the group ``"Firmicutes"`` and include the two phyla.

    Parameters
    ----------
    specs : list[str]
        Raw spec strings, typically from repeated ``--group-phyla`` CLI args.
    all_phyla : set[str]
        The normalised phylum names actually present in the current dataset.
    ph_to_domain : dict[str, str]
        Maps each phylum name to its domain string (e.g. ``"Archaea"``).

    Returns
    -------
    dict[str, str]
        ``{phylum_name: group_label}`` for every phylum that belongs to a
        group.  Phyla absent from the mapping retain their own name as their
        effective label.
    """
    if not specs:
        return {}

    import re as _re

    # Normalise a name to allow case/underscore-flexible matching
    def _norm(n):
        return _re.sub(r'[\s_]+', '_', str(n).strip().lower())

    # Build a lookup from normalised phylum -> actual phylum name
    norm_to_actual = {_norm(ph): ph for ph in all_phyla}

    phylum_to_group: dict = {}

    for spec in specs:
        spec = spec.strip()
        if not spec:
            continue

        lospec = spec.lower()

        if lospec in ('archaea', 'bacteria'):
            domain_key = lospec          # 'archaea' or 'bacteria'
            group_label = lospec.capitalize()
            for ph in all_phyla:
                dom = ph_to_domain.get(ph, '').lower()
                if domain_key in dom:
                    phylum_to_group[ph] = group_label
            continue

        if ':' in spec:
            label, rest = spec.split(':', 1)
            label = label.strip()
            phyla_part = rest.strip()
        else:
            phyla_part = spec
            label = None

        phyla_list = [p.strip() for p in phyla_part.split(',') if p.strip()]
        if not phyla_list:
            continue

        if label is None:
            label = phyla_list[0]   # use first phylum name as the group label

        for requested in phyla_list:
            # Try normalised matching against known phyla
            actual = norm_to_actual.get(_norm(requested))
            if actual:
                phylum_to_group[actual] = label
            else:
                # Phylum not yet in dataset — store anyway; may arrive later
                phylum_to_group[requested] = label

    return phylum_to_group


def write_dataset_colourstrip(output_path: str, dataset_label: str, id_to_colour: dict, legend_title: str = None):
    """Write a simple iTOL DATASET_COLORSTRIP file.

    Parameters
    ----------
    output_path:
        Destination path for the .itol file.
    dataset_label:
        DATASET_LABEL value.
    id_to_colour:
        Mapping of sequence ID to hex colour.
    legend_title:
        Optional legend title. When omitted, a simple legend is derived from
        the unique colours used in id_to_colour.
    """
    p = Path(output_path)
    p.parent.mkdir(parents=True, exist_ok=True)
    unique_items = []
    seen = set()
    for iid, colour in id_to_colour.items():
        key = (str(colour),)
        if key not in seen:
            seen.add(key)
            unique_items.append((iid, colour))

    lines = [
        'DATASET_COLORSTRIP',
        'SEPARATOR COMMA',
        f'DATASET_LABEL,{dataset_label}',
        'COLOR,#AAAAAA',
        'MARGIN,5',
        'SHOW_INTERNAL,0',
    ]
    if unique_items:
        lines.append(f"LEGEND_TITLE,{(legend_title or (dataset_label + ' legend')).replace(',', ';')}")
        lines.append('LEGEND_SHAPES,' + ','.join(['1'] * len(unique_items)))
        lines.append('LEGEND_COLORS,' + ','.join(colour for _, colour in unique_items))
        lines.append('LEGEND_LABELS,' + ','.join(str(iid).replace(',', ';') for iid, _ in unique_items))
    lines.append('DATA')
    for iid, colour in id_to_colour.items():
        lines.append(f'{iid},{colour}')
    p.write_text('\n'.join(lines) + '\n')
    return str(p)


def build_dataset_colour_map(dataset_names):
    dataset_names = sorted({d for d in dataset_names if d})
    QUAL_PALETTE = ['#1f78b4', '#e31a1c', '#33a02c', '#ff7f00', '#6a3d9a', '#b15928', '#a6cee3', '#fb9a99', '#b2df8a', '#fdbf6f', '#cab2d6', '#ffff99']
    ds_colour_map = {}
    if len(dataset_names) <= len(QUAL_PALETTE):
        for i, dn in enumerate(dataset_names):
            ds_colour_map[dn] = QUAL_PALETTE[i]
    else:
        for dn in dataset_names:
            ds_colour_map[dn] = _name_to_dataset_colour(dn, bins=max(12, len(dataset_names)))
    return ds_colour_map


def write_dataset_membership_strip(output_path: str, ids_in_order, ds_map, dataset_label: str = 'Dataset membership', other_colour: str = '#cccccc'):
    ds_colour_map = build_dataset_colour_map(ds_map.values())
    lines = [
        'DATASET_COLORSTRIP',
        'SEPARATOR COMMA',
        f'DATASET_LABEL,{dataset_label}',
        'COLOR,#AAAAAA',
        'MARGIN,5',
        'SHOW_INTERNAL,0',
    ]
    dataset_names = list(ds_colour_map.keys())
    if dataset_names:
        lines.append('LEGEND_TITLE,Dataset membership legend')
        lines.append('LEGEND_SHAPES,' + ','.join(['1'] * len(dataset_names)))
        lines.append('LEGEND_COLORS,' + ','.join(ds_colour_map[d] for d in dataset_names))
        lines.append('LEGEND_LABELS,' + ','.join(d.replace(',', ';') for d in dataset_names))
    lines.append('DATA')
    for iid in ids_in_order:
        ds_name = ds_map.get(iid, '')
        lines.append(f"{iid},{ds_colour_map.get(ds_name, other_colour)}")
    p = Path(output_path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text('\n'.join(lines) + '\n')
    return str(p)


def write_baseline_tier_strip(output_path: str, ids_in_order, tier_map, dataset_label: str = 'Baseline tier'):
    colours = {
        'priority': '#2e7d32',
        'secondary': '#00838f',
    }
    lines = [
        'DATASET_COLORSTRIP',
        'SEPARATOR COMMA',
        f'DATASET_LABEL,{dataset_label}',
        'COLOR,#AAAAAA',
        'MARGIN,5',
        'SHOW_INTERNAL,0',
        'LEGEND_TITLE,Baseline tier legend',
        'LEGEND_SHAPES,1,2',
        'LEGEND_COLORS,#2e7d32,#00838f',
        'LEGEND_LABELS,Hungate baseline (priority),Secondary rumen baseline',
        'DATA',
    ]
    for iid in ids_in_order:
        tier = str(tier_map.get(iid) or '').strip().lower()
        if tier in colours:
            lines.append(f'{iid},{colours[tier]}')
    p = Path(output_path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text('\n'.join(lines) + '\n')
    return str(p)


def generate_itol_colours(taxonomy_tsv: str, outdir: str, user_colour_csv: str = None, id_map: dict = None, tree_file: str = None, phylum_groups: list = None):
    """Generate one iTOL-compatible colour strip per taxonomy metadata type.

    Outputs DATASET_COLORSTRIP files in outdir for phylum, family, genus, and
    optional user colours. TREE_COLORS and DATASET_SYMBOL variants are not
    retained because they duplicate the same metadata in another visual
    encoding and clutter the workflow outputs.

    The user_colour_csv, if provided, must contain ``id`` and ``colour`` columns.

    phylum_groups, if provided, is a list of grouping spec strings passed from
    ``--group-phyla``.  See ``_parse_phylum_groups`` for the spec format.
    Examples::

        phylum_groups=['archaea']
        phylum_groups=['Bacillota,Bacillota_I,Bacillota_A']
        phylum_groups=['archaea', 'Firmicutes:Bacillota,Bacillota_I']
    """
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    # Discover nearby ID maps when callers do not provide one explicitly.
    if id_map is None:
        id_map = {}
        try:
            tdir = Path(taxonomy_tsv).parent
            candidates = []
            for cand_name in ("filing_cabinet_id_map.tsv", "user_id_map.tsv", "user_id_map.csv"):
                p = tdir / cand_name
                if p.exists():
                    candidates.append(p)
            try:
                for p in tdir.glob('*_id_map.tsv'):
                    candidates.append(p)
            except Exception:
                pass
            try:
                for cand_name in ("filing_cabinet_id_map.tsv", "user_id_map.tsv"):
                    p = out / cand_name
                    if p.exists():
                        candidates.append(p)
            except Exception:
                pass

            # Mapping files contain short and original IDs, with an optional header.
            for p in candidates:
                try:
                    with open(p) as fh:
                        first = fh.readline()
                        if first and ('short' in first.lower() or 'original' in first.lower() or '\t' in first):
                            try:
                                parts = first.strip().split('\t')
                                if len(parts) >= 2 and (not parts[0].isdigit() or not parts[1].isdigit()):
                                    pass
                                else:
                                    fh.seek(0)
                            except Exception:
                                fh.seek(0)
                        else:
                            fh.seek(0)
                        for line in fh:
                            parts = line.strip().split('\t')
                            if len(parts) < 2:
                                continue
                            short, orig = parts[0].strip(), parts[1].strip()
                            if not short or not orig:
                                continue
                            id_map[orig] = short
                            id_map[short] = short
                except Exception:
                    continue
        except Exception:
            id_map = {}

    def _map_id(qid: str):
        if not id_map:
            return qid
        if qid in id_map:
            return id_map[qid]
        lk = qid.lower()
        for k, v in id_map.items():
            try:
                if k.lower() == lk:
                    return v
            except Exception:
                continue
        if '|' in qid:
            last = qid.split('|')[-1]
            if last in id_map:
                return id_map[last]
        for k, v in id_map.items():
            if k and qid and (k in qid or qid in k):
                return v
        return qid

    taxa = {}
    with open(taxonomy_tsv) as f:
        next(f, None)
        for line in f:
            parts = line.strip().split('\t')
            if not parts:
                continue
            qid = parts[0]
            # Accept classification and combined taxonomy row layouts.
            tax = 'NA'
            if len(parts) > 3:
                tax = parts[3]
            elif len(parts) > 1:
                tax = parts[1]
            mapped_qid = _map_id(qid)
            if mapped_qid not in taxa:
                taxa[mapped_qid] = tax

    user_colours = {}
    if user_colour_csv:
        import csv
        with open(user_colour_csv) as uc:
            reader = csv.DictReader(uc)
            for row in reader:
                rid = row.get('id') or row.get('ID') or row.get('Id')
                col = row.get('colour') or row.get('Colour')
                if rid and col:
                    mapped = rid
                    if id_map:
                        mapped = _map_id(rid)
                    user_colours[mapped] = col

    combined_lines = []

    parsed_map = {}
    ph_names = {}
    fa_names = {}
    ge_names = {}
    ph_to_domain: dict = {}     # phylum -> domain, used for group-by-domain feature
    for qid, taxstr in taxa.items():
        parsed = parse_taxon_string(taxstr)
        ph = parsed.get('p', 'unknown')
        fa = parsed.get('f', 'unknown')
        ge = parsed.get('g', 'unknown')
        dom = parsed.get('d', '')
        # Keep meaningful phylum subgroups while collapsing numeric suffixes.
        ph = _normalise_taxon_name(ph, preserve_subgroup=True)
        fa = _normalise_taxon_name(fa)
        ge = _normalise_taxon_name(ge)
        parsed_map[qid] = (ph, fa, ge)
        ph_names[ph] = ph_names.get(ph, 0) + 1
        fa_names[fa] = fa_names.get(fa, 0) + 1
        ge_names[ge] = ge_names.get(ge, 0) + 1
        if ph and ph != 'unknown' and dom:
            ph_to_domain.setdefault(ph, dom)

    phylum_to_group = _parse_phylum_groups(
        phylum_groups or [], set(ph_names.keys()), ph_to_domain
    )

    if phylum_to_group:
        # Family and genus remain ungrouped.
        parsed_map = {
            qid: (phylum_to_group.get(ph, ph), fa, ge)
            for qid, (ph, fa, ge) in parsed_map.items()
        }
        group_ph_names: dict = {}
        for ph, count in ph_names.items():
            grp = phylum_to_group.get(ph, ph)
            group_ph_names[grp] = group_ph_names.get(grp, 0) + count
    else:
        group_ph_names = ph_names

    def _make_palette(names, rank):
        names_list = sorted(names.keys())
        n = len(names_list) if names_list else 1
        cmap = {}

        # ColorBrewer qualitative palettes for small category sets.
        BREWER = {
            3: ['#e41a1c', '#377eb8', '#4daf4a'],
            4: ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3'],
            5: ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00'],
            6: ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33'],
            7: ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628'],
            8: ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf'],
            9: ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf', '#999999'],
            10: ['#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd'],
            11: ['#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5'],
            12: ['#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5','#ffed6f']
        }

        if n in BREWER:
            palette = BREWER[n]
            for name, col in zip(names_list, palette):
                cmap[name] = col
            return cmap

        if rank == 'phylum':
            # Greedy Lab-space sampling preserves contrast across many phyla.
            import colorsys

            candidates = []
            for hdeg in range(0, 360, 6):
                h = (hdeg % 360) / 360.0
                for s in (0.95, 0.85, 0.75, 0.65):
                    for v in (0.96, 0.90, 0.80, 0.70):
                        r, g, b = colorsys.hsv_to_rgb(h, s, v)
                        rgb = (int(r*255), int(g*255), int(b*255))
                        hexc = '#{:02x}{:02x}{:02x}'.format(*rgb)
                        candidates.append((hexc, rgb))

            seen = set()
            uniq = []
            for hexc, rgb in candidates:
                if hexc not in seen:
                    seen.add(hexc)
                    uniq.append((hexc, rgb))

            def _rgb_to_lab(rgb):
                r,g,b=[x/255.0 for x in rgb]
                def lin(c):
                    if c <= 0.04045:
                        return c / 12.92
                    return ((c+0.055)/1.055)**2.4
                r,g,b=lin(r),lin(g),lin(b)
                X=r*0.4124564+g*0.3575761+b*0.1804375
                Y=r*0.2126729+g*0.7151522+b*0.0721750
                Z=r*0.0193339+g*0.1191920+b*0.9503041
                Xn,Yn,Zn=0.95047,1.0,1.08883
                x,y,z=X/Xn,Y/Yn,Z/Zn
                def f(t):
                    if t > 0.008856:
                        return t**(1/3)
                    return 7.787*t+16/116
                fx,fy,fz=f(x),f(y),f(z)
                L=116*fy-16
                a=500*(fx-fy)
                bval=200*(fy-fz)
                return (L,a,bval)

            cand_lab = [(_hex, _rgb_to_lab(rgb)) for (_hex, rgb) in uniq]

            chosen = []
            chosen_lab = []
            try:
                seed_idx = int((_hash_to_hue(''.join(names_list)) * 10000)) % len(cand_lab)
            except Exception:
                seed_idx = 0
            chosen.append(cand_lab.pop(seed_idx))
            chosen_lab.append(chosen[-1][1])

            def _sqdist_lab(a,b):
                return (a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2

            while len(chosen) < n and cand_lab:
                best_i = None
                best_min = -1
                for i,(hexc, lab) in enumerate(cand_lab):
                    mins = min(_sqdist_lab(lab, cl) for cl in chosen_lab)
                    if mins > best_min:
                        best_min = mins
                        best_i = i
                if best_i is None:
                    break
                chosen.append(cand_lab.pop(best_i))
                chosen_lab.append(chosen[-1][1])

            for name, (hexc, _) in zip(names_list, chosen):
                cmap[name] = hexc

            if len(cmap) < len(names_list):
                remaining = [nm for nm in names_list if nm not in cmap]
                m = len(remaining)
                for i, name in enumerate(remaining):
                    h = (i / float(m)) % 1.0
                    cmap[name] = _hsv_to_hex(h, s=0.85, v=0.95)
            return cmap

        # Offset evenly spaced hues by rank.
        rank_offsets = {
            'phylum': 0.0,
            'family': 0.15,
            'genus': 0.30,
        }
        offset = rank_offsets.get(rank, 0.0)
        sv_map = {
            'phylum': (0.85, 0.97),
            'family': (0.65, 0.90),
            'genus': (0.45, 0.80),
        }
        s, v = sv_map.get(rank, (0.65, 0.90))
        for i, name in enumerate(names_list):
            h = ((i / float(n)) + offset) % 1.0
            cmap[name] = _hsv_to_hex(h, s=s, v=v)
        return cmap

    phylum_map = _make_palette(group_ph_names, 'phylum')
    family_map = _make_palette(fa_names, 'family')
    genus_map = _make_palette(ge_names, 'genus')

    # Stable node labels let iTOL metadata target the final tree precisely.
    newick = None
    if tree_file:
        try:
            tf = Path(tree_file)
            if tf.exists():
                newick = tf.read_text()
                try:
                    from branchmanager.pipeline.tree import _repair_internal_node_label_delimiters
                    newick = _repair_internal_node_label_delimiters(newick)
                except Exception:
                    pass
                labelled_newick = newick
                try:
                    import re
                    existing = set(re.findall(r"\bNODE\d{4}\b", newick))
                    counter = 0
                    def _next_node():
                        nonlocal counter
                        while True:
                            name = f"NODE{counter:04d}"
                            counter += 1
                            if name not in existing:
                                existing.add(name)
                                return name

                    def _repl(m):
                        node = _next_node()
                        return f"){node}_{m.group(1)}"

                    labelled_newick = re.sub(r"\)\s*([0-9]+(?:\.[0-9]+)?)", _repl, labelled_newick)
                except Exception:
                    labelled_newick = newick
                try:
                    tf.write_text(labelled_newick)
                    newick = labelled_newick
                except Exception:
                    try:
                        out_labelled = out / 'current_tree_labelled.nwk'
                        out_labelled.write_text(labelled_newick)
                        newick = labelled_newick
                    except Exception:
                        pass
        except Exception:
            newick = None

    # User colours override rank palettes for individual sequences.
    for qid in list(taxa.keys()):
        ph, fa, ge = parsed_map[qid]
        user_col = user_colours.get(qid, '')
        ph_c = user_col or phylum_map.get(ph, _name_to_colour_by_rank(ph, 'phylum'))
        fa_c = user_col or family_map.get(fa, _name_to_colour_by_rank(fa, 'family'))
        ge_c = user_col or genus_map.get(ge, _name_to_colour_by_rank(ge, 'genus'))
        combined_lines.append((qid, ph_c, fa_c, ge_c, user_col))


    def write_colourstrip(path, title, id_colour_pairs, legend_pairs=None):
        """Write a minimal iTOL DATASET_COLORSTRIP file.

        legend_pairs contains optional ``(label, colour)`` entries.
        """
        with open(path, 'w') as f:
            f.write('DATASET_COLORSTRIP\n')
            f.write('SEPARATOR COMMA\n')
            f.write(f'DATASET_LABEL,{title}\n')
            if legend_pairs and len(legend_pairs) > 0:
                dataset_colour = legend_pairs[0][1]
            else:
                dataset_colour = id_colour_pairs[0][1] if id_colour_pairs else '#AAAAAA'
            f.write(f'COLOR,{dataset_colour}\n')
            f.write('MARGIN,5\n')
            f.write('SHOW_INTERNAL,0\n')

            if legend_pairs:
                labels = [lbl.replace(',', ';') for lbl, _ in legend_pairs]
                colours = [col for _, col in legend_pairs]
                shapes = ['1'] * len(legend_pairs)
                f.write(f"LEGEND_TITLE,{title} legend\n")
                f.write('LEGEND_SHAPES,' + ','.join(shapes) + '\n')
                f.write('LEGEND_COLORS,' + ','.join(colours) + '\n')
                f.write('LEGEND_LABELS,' + ','.join(labels) + '\n')

            f.write('DATA\n')
            for item in id_colour_pairs:
                if len(item) >= 2:
                    id_, col = item[0], item[1]
                else:
                    continue
                f.write(f"{id_},{col}\n")

    # Limit colour-strip rows to leaves present in the supplied tree.
    leaf_ids_in_tree = set()
    if tree_file:
        try:
            import re as _re
            if not newick:
                raise ValueError('no newick content')
            tokens = _re.findall(r"([^(),:;]+)(?=[:),;])", newick)
            for t in tokens:
                tt = t.strip()
                if not tt:
                    continue
                mapped = _map_id(tt)
                leaf_ids_in_tree.add(mapped)
        except Exception:
            leaf_ids_in_tree = set()

    ph_pairs = [(qid, ph_c, parsed_map[qid][0]) for qid, ph_c, _, _, _ in combined_lines]
    fa_pairs = [(qid, fa_c, parsed_map[qid][1]) for qid, _, fa_c, _, _ in combined_lines]
    ge_pairs = [(qid, ge_c, parsed_map[qid][2]) for qid, _, _, ge_c, _ in combined_lines]

    if leaf_ids_in_tree:
        # Two-column collapse maps resolve members to representatives in the tree.
        collapse_map = {}
        try:
            cand_paths = [
                out / 'filing_cabinet_collapsed_members.tsv',
                out / 'collapsed_members.tsv',
                out / 'filing_cabinet_collapsed_map.tsv',
                out / 'collapsed_map.tsv',
                Path(taxonomy_tsv).parent / 'filing_cabinet_collapsed_members.tsv',
                Path(taxonomy_tsv).parent / 'filing_cabinet_collapsed_map.tsv'
            ]
            for cp in cand_paths:
                try:
                    if not cp.exists():
                        continue
                    with open(cp) as cf:
                        next(cf, None)
                        for line in cf:
                            if not line.strip():
                                continue
                            parts = line.strip().split('\t')
                            if len(parts) == 2:
                                collapse_map[parts[0]] = parts[1]
                except Exception:
                    continue
        except Exception:
            collapse_map = {}

        def _candidate_forms(x):
            forms = []
            try:
                forms.append(_map_id(x))
            except Exception:
                forms.append(x)
            forms.append(x)
            if '|' in x:
                forms.append(x.split('|')[-1])
            forms.append(x.split()[-1] if x.split() else x)
            try:
                forms.append(forms[0].upper())
                forms.append(forms[0].lower())
            except Exception:
                pass
            import re as _re
            m = _re.sub(r'_n?\d+$','', forms[0] if forms else x)
            forms.append(m)
            seen = set()
            out = []
            for f in forms:
                if not f:
                    continue
                if f not in seen:
                    seen.add(f)
                    out.append(f)
            return out

        def _resolve_to_tree_id(qid):
            for cand in _candidate_forms(qid):
                if cand in leaf_ids_in_tree:
                    return cand
            return None

        def _build_tree_pairs(pairs):
            out = []
            used = set()
            missing = []
            for qid, col, tname in pairs:
                tree_id = _resolve_to_tree_id(qid)
                if not tree_id:
                    missing.append(qid)
                    continue
                if tree_id in used:
                    continue
                used.add(tree_id)
                out.append((tree_id, col, tname))
            try:
                if missing:
                    with open(outdir + '/itol_missing_ids.tsv', 'w') as mf:
                        mf.write('id\n')
                        for m in sorted(set(missing)):
                            mf.write(m + '\n')
            except Exception:
                pass
            return out

        ph_pairs_tree = _build_tree_pairs(ph_pairs)
        fa_pairs_tree = _build_tree_pairs(fa_pairs)
        ge_pairs_tree = _build_tree_pairs(ge_pairs)
    else:
        ph_pairs_tree = ph_pairs
        fa_pairs_tree = fa_pairs
        ge_pairs_tree = ge_pairs

    ph_legend = sorted([(name + f" ({group_ph_names.get(name,0)})", col) for name, col in phylum_map.items()], key=lambda x: x[0])
    fa_legend = sorted([(name + f" ({fa_names.get(name,0)})", col) for name, col in family_map.items()], key=lambda x: x[0])
    ge_legend = sorted([(name + f" ({ge_names.get(name,0)})", col) for name, col in genus_map.items()], key=lambda x: x[0])

    if leaf_ids_in_tree:
        ph_dataset_pairs = [(qid, col) for qid, col, _ in ph_pairs_tree]
        fa_dataset_pairs = [(qid, col) for qid, col, _ in fa_pairs_tree]
        ge_dataset_pairs = [(qid, col) for qid, col, _ in ge_pairs_tree]
    else:
        ph_dataset_pairs = [(qid, col) for qid, col, _ in ph_pairs]
        fa_dataset_pairs = [(qid, col) for qid, col, _ in fa_pairs]
        ge_dataset_pairs = [(qid, col) for qid, col, _ in ge_pairs]

    write_colourstrip(out / 'itol_phylum_colours.itol', 'Phylum colours', ph_dataset_pairs, legend_pairs=ph_legend)
    write_colourstrip(out / 'itol_family_colours.itol', 'Family colours', fa_dataset_pairs, legend_pairs=fa_legend)
    write_colourstrip(out / 'itol_genus_colours.itol', 'Genus colours', ge_dataset_pairs, legend_pairs=ge_legend)

    if user_colours:
        user_pairs = [(rid, col) for rid, col in user_colours.items()]
        write_colourstrip(out / 'itol_user_colours.itol', 'User colours', user_pairs)

    return str(out)


def write_functional_annotations(
    func_tsv: str,
    outdir: str,
    id_map: dict = None,
) -> list:
    """Generate iTOL annotation files from a functional attributes TSV.

    The TSV must have a header row.  The first column is the sequence ID; every
    subsequent column is a functional attribute (pathway, function, binary
    trait, numeric score, etc.).

    One iTOL file is written per column:

    * **Binary columns** (values limited to 0/1/yes/no/true/false/+/-) →
      ``DATASET_BINARY`` (presence/absence squares).
    * **Numeric columns** (all non-empty values parse as float) →
      ``DATASET_SIMPLEBAR`` (horizontal bar chart).
    * **Categorical columns** (everything else) →
      ``DATASET_COLORSTRIP`` (colour-coded strip with legend).

    Output files are named ``itol_func_<column_name>.itol`` inside *outdir*.

    Parameters
    ----------
    func_tsv:
        Path to the tab-separated annotation file.
    outdir:
        Directory where iTOL files will be written.
    id_map:
        Optional ``{original_id: short_id}`` mapping.  When provided,
        sequence IDs in the TSV are mapped to the short IDs used in the tree.

    Returns
    -------
    list of str
        Paths to the written iTOL files.
    """
    import csv
    import re as _re

    _BINARY_VALUES = {'0', '1', 'yes', 'no', 'true', 'false', '+', '-'}
    _EMPTY_VALUES  = {'', 'na', 'n/a', '-', 'nd', 'none', 'null', '.'}

    out_p = Path(outdir)
    out_p.mkdir(parents=True, exist_ok=True)

    def _map_id(qid):
        if not id_map:
            return qid
        if qid in id_map:
            return id_map[qid]
        lk = qid.lower()
        for k, v in id_map.items():
            try:
                if k.lower() == lk:
                    return v
            except Exception:
                continue
        return qid

    with open(func_tsv, newline='') as fh:
        reader = csv.reader(fh, delimiter='\t')
        header = next(reader, None)
        if not header or len(header) < 2:
            return []
        attr_cols = header[1:]
        rows = []
        for line in reader:
            if not line:
                continue
            row_id = _map_id(line[0].strip())
            vals = []
            for i in range(len(attr_cols)):
                vals.append(line[i + 1].strip() if i + 1 < len(line) else '')
            rows.append((row_id, vals))

    if not rows:
        return []

    n_attrs = len(attr_cols)
    written = []

    for col_idx in range(n_attrs):
        col_name  = attr_cols[col_idx].strip()
        col_values = [row[1][col_idx] for row in rows]

        safe_name = _re.sub(r'[^\w]+', '_', col_name).strip('_') or f'col{col_idx+1}'

        non_empty = [v.lower() for v in col_values if v.lower() not in _EMPTY_VALUES]

        if not non_empty:
            continue   # all values are empty — skip

        is_binary  = all(v in _BINARY_VALUES for v in non_empty)
        is_numeric = False
        if not is_binary:
            try:
                [float(v) for v in non_empty]
                is_numeric = True
            except ValueError:
                pass

        dest = out_p / f'itol_func_{safe_name}.itol'

        if is_binary:
            def _to_01(v):
                lv = v.lower()
                if lv in ('1', 'yes', 'true', '+'):
                    return '1'
                return '0'

            col_colour = _name_to_colour(col_name)
            lines_out = [
                'DATASET_BINARY',
                'SEPARATOR COMMA',
                f'DATASET_LABEL,{col_name}',
                f'COLOR,{col_colour}',
                f'FIELD_LABELS,{col_name}',
                f'FIELD_COLORS,{col_colour}',
                'FIELD_SHAPES,1',
                'SHOW_INTERNAL,0',
                'DATA',
            ]
            for row_id, vals in rows:
                v = vals[col_idx]
                lines_out.append(f'{row_id},{_to_01(v)}')
            dest.write_text('\n'.join(lines_out) + '\n')
            written.append(str(dest))

        elif is_numeric:
            col_colour = _name_to_colour(col_name)
            try:
                float_vals = [float(v) for v in non_empty]
                max_val = max(float_vals)
            except Exception:
                max_val = 1.0
            lines_out = [
                'DATASET_SIMPLEBAR',
                'SEPARATOR COMMA',
                f'DATASET_LABEL,{col_name}',
                f'COLOR,{col_colour}',
                'WIDTH,200',
                f'MAXIMUM_VALUE,{max_val}',
                'SHOW_INTERNAL,0',
                'DATA',
            ]
            for row_id, vals in rows:
                v = vals[col_idx]
                if v.lower() in _EMPTY_VALUES:
                    continue
                try:
                    lines_out.append(f'{row_id},{float(v)}')
                except ValueError:
                    pass
            dest.write_text('\n'.join(lines_out) + '\n')
            written.append(str(dest))

        else:
            categories = []
            seen_cats = set()
            for v in col_values:
                lv = v.lower()
                if lv not in _EMPTY_VALUES and v not in seen_cats:
                    seen_cats.add(v)
                    categories.append(v)

            cat_palette = _make_palette_simple(categories)

            legend_shapes  = ','.join(['1'] * len(categories))
            legend_colours  = ','.join(cat_palette[c] for c in categories)
            legend_labels  = ','.join(c.replace(',', ';') for c in categories)
            first_col = cat_palette[categories[0]] if categories else '#AAAAAA'

            lines_out = [
                'DATASET_COLORSTRIP',
                'SEPARATOR COMMA',
                f'DATASET_LABEL,{col_name}',
                f'COLOR,{first_col}',
                'MARGIN,5',
                'SHOW_INTERNAL,0',
                f'LEGEND_TITLE,{col_name} legend',
                f'LEGEND_SHAPES,{legend_shapes}',
                f'LEGEND_COLORS,{legend_colours}',
                f'LEGEND_LABELS,{legend_labels}',
                'DATA',
            ]
            for row_id, vals in rows:
                v = vals[col_idx]
                if v.lower() in _EMPTY_VALUES:
                    continue
                col_hex = cat_palette.get(v, '#cccccc')
                lines_out.append(f'{row_id},{col_hex},{v.replace(",", ";")}')
            dest.write_text('\n'.join(lines_out) + '\n')
            written.append(str(dest))

    return written


def _make_palette_simple(names: list) -> dict:
    """Return a ``{name: hex_colour}`` palette for category names.

    Uses the same greedy Lab-space sampling as ``_make_palette`` for phyla
    when there are many categories, falling back to ColorBrewer for small
    sets.
    """
    n = len(names)
    if n == 0:
        return {}
    BREWER = {
        1: ['#1f78b4'],
        2: ['#1f78b4', '#e31a1c'],
        3: ['#e41a1c', '#377eb8', '#4daf4a'],
        4: ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3'],
        5: ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00'],
        6: ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33'],
        7: ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628'],
        8: ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf'],
        9: ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf', '#999999'],
        10: ['#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd'],
        11: ['#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5'],
        12: ['#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5','#ffed6f'],
    }
    if n in BREWER:
        return {name: col for name, col in zip(names, BREWER[n])}
    return {name: _hsv_to_hex(i / float(n), s=0.80, v=0.92) for i, name in enumerate(names)}



# Broad-scale ruminant microbiome functional group assignments.
# Keys are lowercase genus, family, order, or phylum names (or partial matches).
# Values are human-readable functional group labels.
# Hierarchy: genus match takes priority over family, then order, then phylum.
# A sequence not matching any entry defaults to 'Other / Unclassified'.

# Curated at genus level first (most specific)
_RUMEN_GENUS: dict = {
    # Cellulolytic / Fibrolytic
    'ruminococcus': 'Cellulolytic/Fibrolytic',
    'hungateiclostridium': 'Cellulolytic/Fibrolytic',
    'clostridium': 'Cellulolytic/Fibrolytic',
    'caldicellulosiruptor': 'Cellulolytic/Fibrolytic',
    'cellulosilyticum': 'Cellulolytic/Fibrolytic',
    'pseudobutyrivibrio': 'Cellulolytic/Fibrolytic',
    # Proteolytic
    'prevotella': 'Proteolytic/Peptidolytic',
    'bacteroides': 'Proteolytic/Peptidolytic',
    'butyrivibrio': 'Proteolytic/Peptidolytic',
    'selenomonas': 'Proteolytic/Peptidolytic',
    'anaerovibrio': 'Lipolytic',
    # Methanogenic Archaea
    'methanobrevibacter': 'Methanogenic Archaea',
    'methanobacterium': 'Methanogenic Archaea',
    'methanomicrobium': 'Methanogenic Archaea',
    'methanosarcina': 'Methanogenic Archaea',
    'methanosphaera': 'Methanogenic Archaea',
    'methanomassiliicoccus': 'Methanogenic Archaea',
    'methanoregula': 'Methanogenic Archaea',
    'methanoplanus': 'Methanogenic Archaea',
    'methanoculleus': 'Methanogenic Archaea',
    'methanocorpusculum': 'Methanogenic Archaea',
    # Butyrate producers
    'roseburia': 'Butyrate Producers',
    'faecalibacterium': 'Butyrate Producers',
    'eubacterium': 'Butyrate Producers',
    'coprococcus': 'Butyrate Producers',
    'anaerobutyricum': 'Butyrate Producers',
    'anaerostipes': 'Butyrate Producers',
    'lachnospira': 'Butyrate Producers',
    'blautia': 'Butyrate Producers',
    # Succinate/Propionate producers
    'propionibacterium': 'Succinate/Propionate Producers',
    'propionigenium': 'Succinate/Propionate Producers',
    'veillonella': 'Succinate/Propionate Producers',
    'succiniclasticum': 'Succinate/Propionate Producers',
    'succinivibrio': 'Succinate/Propionate Producers',
    'succinimonas': 'Succinate/Propionate Producers',
    'ruminobacter': 'Succinate/Propionate Producers',
    'fibrobacter': 'Cellulolytic/Fibrolytic',
    # Lactate producers
    'streptococcus': 'Lactate Producers',
    'lactobacillus': 'Lactate Producers',
    'lactococcus': 'Lactate Producers',
    'ligilactobacillus': 'Lactate Producers',
    'enterococcus': 'Lactate Producers',
    'megasphaera': 'Lactate Producers',
    # Acetate producers / Acetogens
    'acetitomaculum': 'Acetogens',
    'acetoanaerobium': 'Acetogens',
    'moorella': 'Acetogens',
    'sporomusa': 'Acetogens',
    'ruminococcaceae': 'Acetogens',
    # Amylolytic / starch degraders
    'amylophilus': 'Amylolytic/Starch Degraders',
    'treponema': 'Amylolytic/Starch Degraders',
    # Sulfate-reducing bacteria
    'desulfovibrio': 'Sulfate-Reducing Bacteria',
    'desulfobulbus': 'Sulfate-Reducing Bacteria',
    'desulfobacter': 'Sulfate-Reducing Bacteria',
}

# Family-level fallback
_RUMEN_FAMILY: dict = {
    'ruminococcaceae': 'Cellulolytic/Fibrolytic',
    'lachnospiraceae': 'Butyrate Producers',
    'clostridiaceae': 'Cellulolytic/Fibrolytic',
    'bacteroidaceae': 'Proteolytic/Peptidolytic',
    'prevotellaceae': 'Proteolytic/Peptidolytic',
    'fibrobacteraceae': 'Cellulolytic/Fibrolytic',
    'methanobacteriaceae': 'Methanogenic Archaea',
    'methanomicrobiaceae': 'Methanogenic Archaea',
    'methanosaetaceae': 'Methanogenic Archaea',
    'methanosarcinaceae': 'Methanogenic Archaea',
    'veillonellaceae': 'Succinate/Propionate Producers',
    'selenomonadaceae': 'Succinate/Propionate Producers',
    'succinivibrionaceae': 'Succinate/Propionate Producers',
    'streptococcaceae': 'Lactate Producers',
    'lactobacillaceae': 'Lactate Producers',
    'desulfovibrionaceae': 'Sulfate-Reducing Bacteria',
    'spirochaetaceae': 'Amylolytic/Starch Degraders',
    'treponemataceae': 'Amylolytic/Starch Degraders',
    'butyrivibrionaceae': 'Cellulolytic/Fibrolytic',
}

# Phylum-level fallback
_RUMEN_PHYLUM: dict = {
    'fibrobacterota': 'Cellulolytic/Fibrolytic',
    'fibrobacteres': 'Cellulolytic/Fibrolytic',
    'methanobacteriota': 'Methanogenic Archaea',
    'euryarchaeota': 'Methanogenic Archaea',
    'thermoplasmota': 'Methanogenic Archaea',
    # Domain-level archaea fallback applied separately
    'bacteroidota': 'Proteolytic/Peptidolytic',
    'bacteroidetes': 'Proteolytic/Peptidolytic',
    'bacillota': 'Butyrate Producers',
    'firmicutes': 'Butyrate Producers',
    'spirochaetota': 'Amylolytic/Starch Degraders',
    'spirochaetes': 'Amylolytic/Starch Degraders',
    'proteobacteria': 'Succinate/Propionate Producers',
    'pseudomonadota': 'Succinate/Propionate Producers',
    'desulfobacterota': 'Sulfate-Reducing Bacteria',
    'actinobacteriota': 'Succinate/Propionate Producers',
    'actinobacteria': 'Succinate/Propionate Producers',
    'verrucomicrobiota': 'Glycan/Mucin Degraders',
    'verrucomicrobia': 'Glycan/Mucin Degraders',
}

_RUMEN_FUNC_PALETTE: dict = {
    'Cellulolytic/Fibrolytic':      '#1b7837',  # dark green
    'Proteolytic/Peptidolytic':     '#762a83',  # purple
    'Methanogenic Archaea':         '#4393c3',  # sky blue
    'Butyrate Producers':           '#d6604d',  # coral red
    'Succinate/Propionate Producers': '#f4a582', # light orange
    'Lactate Producers':            '#fdae61',  # amber
    'Acetogens':                    '#92c5de',  # pale blue
    'Amylolytic/Starch Degraders':  '#e7d4e8',  # lavender
    'Lipolytic':                    '#fee08b',  # yellow
    'Sulfate-Reducing Bacteria':    '#878787',  # grey
    'Glycan/Mucin Degraders':       '#a6d96a',  # lime
    'Other / Unclassified':         '#cccccc',  # neutral grey
}


def _assign_rumen_function(parsed: dict) -> str:
    """Return a broad ruminant functional group label for a parsed taxonomy dict.

    Parameters
    ----------
    parsed : dict
        Output of :func:`parse_taxon_string` with keys 'd','p','c','o','f','g','s'.

    Returns
    -------
    str
        Functional group label (one of the keys in :data:`_RUMEN_FUNC_PALETTE`).
    """
    def _lo(x):
        return (x or '').lower().strip()

    domain = _lo(parsed.get('d', ''))
    phylum = _lo(_normalise_taxon_name(parsed.get('p', ''), preserve_subgroup=True))
    family = _lo(_normalise_taxon_name(parsed.get('f', '')))
    genus  = _lo(_normalise_taxon_name(parsed.get('g', '')))

    if 'archaea' in domain:
        hit = _RUMEN_GENUS.get(genus)
        if hit:
            return hit
        return 'Methanogenic Archaea'

    for lookup, key in (
        (_RUMEN_GENUS,  genus),
        (_RUMEN_FAMILY, family),
        (_RUMEN_PHYLUM, phylum),
    ):
        hit = lookup.get(key)
        if hit:
            return hit
        # Prefix matching handles subgroup suffixes such as ``bacillota_a``.
        for k, v in lookup.items():
            if key and k and (key.startswith(k) or k.startswith(key)):
                return v

    return 'Other / Unclassified'


def generate_rumen_function_draft(
    taxonomy_tsv: str,
    outdir: str,
    id_map: dict = None,
) -> tuple:
    """Auto-generate a draft rumen functional annotation from a taxonomy TSV.

    Reads the combined taxonomy TSV produced by a Performance Review and maps each
    sequence to one of the broad ruminant microbiome functional groups defined
    in :data:`_RUMEN_FUNC_PALETTE`.

    Two output files are written in *outdir*:

    ``rumen_functions_draft.tsv``
        Tab-separated.  Columns: ``sequence_id``, ``Rumen_Functional_Group``.
        This file can be edited by the user and then supplied to future runs
        via ``--functional`` for a refined annotation.

    ``itol_func_Rumen_Functional_Group.itol``
        iTOL ``DATASET_COLORSTRIP`` file ready for direct upload to iTOL after
        the tree.

    Parameters
    ----------
    taxonomy_tsv : str
        Path to the taxonomy TSV (combined_taxonomy.tsv or similar).
        Supports the ``ID\\tTaxon\\tConfidence`` and
        ``ID\\tBestHit\\tIdentity\\tTaxon\\tConfidence`` formats.
    outdir : str
        Output directory.
    id_map : dict, optional
        ``{original_id: short_id}`` mapping to translate sequence IDs.

    Returns
    -------
    (tsv_path, itol_path) : tuple of str
        Paths to the written files.  Returns ``(None, None)`` on failure.
    """
    import csv as _csv

    out_p = Path(outdir)
    out_p.mkdir(parents=True, exist_ok=True)

    def _map_id(qid: str):
        if not id_map:
            return qid
        if qid in id_map:
            return id_map[qid]
        lk = qid.lower()
        for k, v in (id_map or {}).items():
            try:
                if k.lower() == lk:
                    return v
            except Exception:
                continue
        return qid

    rows: list = []
    try:
        import gzip as _gzip
        _opener = _gzip.open if str(taxonomy_tsv).endswith(('.gz', '.gzip')) else open
        with _opener(taxonomy_tsv, 'rt') as fh:
            reader = _csv.reader(fh, delimiter='\t')
            next(reader, None)
            for line in reader:
                if not line:
                    continue
                qid = line[0].strip()
                taxon = 'NA'
                if len(line) > 3:
                    taxon = line[3].strip()
                elif len(line) > 1:
                    taxon = line[1].strip()
                short_id = _map_id(qid)
                parsed   = parse_taxon_string(taxon)
                func_grp = _assign_rumen_function(parsed)
                rows.append((short_id, func_grp))
    except Exception as e:
        import logging as _log
        _log.getLogger(__name__).warning(
            "generate_rumen_function_draft: failed to read %s: %s", taxonomy_tsv, e
        )
        return None, None

    if not rows:
        return None, None

    tsv_path = out_p / 'rumen_functions_draft.tsv'
    try:
        with open(tsv_path, 'w', newline='') as fh:
            w = _csv.writer(fh, delimiter='\t')
            w.writerow(['sequence_id', 'Rumen_Functional_Group'])
            w.writerows(rows)
    except Exception as e:
        import logging as _log
        _log.getLogger(__name__).warning(
            "generate_rumen_function_draft: failed to write TSV: %s", e
        )
        return None, None

    # Prefer curated colours and generate any additional categories.
    cats_seen = []
    cats_set: set = set()
    for _, grp in rows:
        if grp not in cats_set:
            cats_set.add(grp)
            cats_seen.append(grp)

    cat_palette = {}
    unknown_cats = []
    for c in cats_seen:
        if c in _RUMEN_FUNC_PALETTE:
            cat_palette[c] = _RUMEN_FUNC_PALETTE[c]
        else:
            unknown_cats.append(c)
    if unknown_cats:
        extra = _make_palette_simple(unknown_cats)
        cat_palette.update(extra)

    legend_shapes = ','.join(['1'] * len(cats_seen))
    legend_colours = ','.join(cat_palette.get(c, '#cccccc') for c in cats_seen)
    legend_labels = ','.join(c.replace(',', ';') for c in cats_seen)
    first_col = cat_palette.get(cats_seen[0], '#cccccc') if cats_seen else '#cccccc'

    itol_lines = [
        'DATASET_COLORSTRIP',
        'SEPARATOR COMMA',
        'DATASET_LABEL,Rumen Functional Group',
        f'COLOR,{first_col}',
        'MARGIN,5',
        'SHOW_INTERNAL,0',
        'LEGEND_TITLE,Rumen Functional Group legend',
        f'LEGEND_SHAPES,{legend_shapes}',
        f'LEGEND_COLORS,{legend_colours}',
        f'LEGEND_LABELS,{legend_labels}',
        'DATA',
    ]
    for sid, grp in rows:
        col_hex = cat_palette.get(grp, '#cccccc')
        itol_lines.append(f'{sid},{col_hex},{grp.replace(",", ";")}')

    itol_path = out_p / 'itol_func_Rumen_Functional_Group.itol'
    try:
        itol_path.write_text('\n'.join(itol_lines) + '\n')
    except Exception as e:
        import logging as _log
        _log.getLogger(__name__).warning(
            "generate_rumen_function_draft: failed to write iTOL file: %s", e
        )
        return str(tsv_path), None

    return str(tsv_path), str(itol_path)
