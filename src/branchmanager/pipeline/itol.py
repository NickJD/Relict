from pathlib import Path
import hashlib
from branchmanager.taxonomy import parse_taxon_string as _shared_parse_taxon_string


def _hash_to_hue(s: str) -> float:
    # deterministic hash to [0,1)
    h = hashlib.md5(s.encode('utf-8')).hexdigest()
    val = int(h[:8], 16)
    return (val % 360) / 360.0


def _hsv_to_hex(h, s=0.65, v=0.95):
    # convert hsv to hex
    import colorsys
    r, g, b = colorsys.hsv_to_rgb(h, s, v)
    return '#{:02x}{:02x}{:02x}'.format(int(r*255), int(g*255), int(b*255))


def _name_to_color(name: str) -> str:
    h = _hash_to_hue(name)
    return _hsv_to_hex(h)


def _name_to_dataset_color(name: str, bins: int = 8) -> str:
    """Generate a visually distinct color for dataset-level labels.

    This quantizes the hue space into `bins` discrete buckets to ensure
    different dataset names map to visibly different hues (avoids two
    blue-ish colors for similar hashes). A small jitter is added for
    deterministic variation within a bucket.
    """
    # quantize hue into bins
    base_h = _hash_to_hue(name)
    idx = int(base_h * bins) % bins
    # center of bin
    h = (idx + 0.5) / bins
    # add a small deterministic jitter based on the hash to avoid exact
    # collisions across names that fall in the same bin center
    h_frac = _hash_to_hue(name + '_jitter')
    jitter = (h_frac - 0.5) * (0.5 / bins)
    h = (h + jitter) % 1.0
    # use higher saturation/value for dataset markers for strong contrast
    return _hsv_to_hex(h, s=0.85, v=0.92)


def _identity_to_color(identity: float, vmin: float = 0.0, vmax: float = 1.0) -> str:
    """Map an identity/similarity score to a color on a red->green gradient.

    - identity can be in 0..1 or 0..100. If None, returns a neutral gray.
    - vmin/vmax allow clipping/normalisation; values outside are clipped.
    """
    if identity is None:
        return '#cccccc'
    try:
        val = float(identity)
    except Exception:
        return '#cccccc'
    # handle percent values
    if val > 1.0 and val <= 100.0:
        val = val / 100.0
    # clip
    if val is None or val != val:  # NaN
        return '#cccccc'
    val = max(vmin, min(vmax, val))
    # normalize to 0..1
    if vmax - vmin > 0:
        norm = (val - vmin) / (vmax - vmin)
    else:
        norm = 0.0
    # map 0 -> red (h=0.0), 1 -> greenish (h=0.33)
    h = norm * 0.33
    # slightly stronger saturation/value for clearer distinction
    return _hsv_to_hex(h, s=0.95, v=0.95)


def _novelty_color_for_pct(pct, min_pct=90, max_pct=100):
    """Return a distinct color for a novelty category.

    pct may be '<90' or an int between min_pct and max_pct inclusive.
    Colors are evenly spaced between red (low) and green (high).
    """
    # red hue for novel (<90)
    if pct == '<90' or pct is None:
        return _hsv_to_hex(0.0, s=0.9, v=0.9)
    try:
        p = int(pct)
    except Exception:
        return _hsv_to_hex(0.0, s=0.9, v=0.9)
    # clamp
    if p < min_pct:
        return _hsv_to_hex(0.0, s=0.9, v=0.9)
    if p > max_pct:
        p = max_pct
    # normalized 0..1 across the range
    norm = (p - min_pct) / float(max_pct - min_pct) if max_pct > min_pct else 1.0
    h = norm * 0.33
    return _hsv_to_hex(h, s=0.9, v=0.9)


def _name_to_color_by_rank(name: str, rank: str = None) -> str:
    """Generate a deterministic color for a name, adjusted by rank to
    reduce color reuse between hierarchical ranks (phylum/family/genus).

    rank can be 'phylum', 'family', 'genus' or None. Different ranks receive
    a hue offset and slightly different saturation/value to increase contrast.
    """
    h = _hash_to_hue(name)
    # rank offsets chosen to separate palettes more strongly across the hue wheel
    # larger separation reduces color reuse between hierarchical ranks
    rank_offsets = {
        'phylum': 0.0,
        'family': 0.33,
        'genus': 0.66,
    }
    offset = rank_offsets.get(rank, 0.0)
    h = (h + offset) % 1.0

    # choose saturation/value per rank for additional separation; more distinct
    # settings create visually different palettes between ranks
    sv_map = {
        'phylum': (0.85, 0.97),
        'family': (0.65, 0.90),
        'genus': (0.45, 0.80),
    }
    s, v = sv_map.get(rank, (0.65, 0.90))
    return _hsv_to_hex(h, s=s, v=v)


def parse_taxon_string(taxon: str):
    return _shared_parse_taxon_string(taxon)


def _normalize_taxon_name(name: str, preserve_subgroup: bool = False) -> str:
    """Normalize taxon names produced by gg2/reference sources.

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
    # remove surrounding quotes
    if (s.startswith('"') and s.endswith('"')) or (s.startswith("'") and s.endswith("'")):
        s = s[1:-1].strip()
    # collapse whitespace
    s = ' '.join(s.split())
    # split on underscores and spaces to inspect final tokens
    import re
    parts = re.split(r'[\s_]+', s)
    if not parts:
        return s
    # drop trailing numeric-only tokens (accession-like ids)
    while parts and parts[-1].isdigit():
        parts.pop()
    # drop single-letter subgroup tokens (e.g. Bacillota_A -> Bacillota) only
    # when preserve_subgroup is False. When preserve_subgroup is True we keep
    # such subgroup suffixes (useful for phylum-level labels like Bacillota_A).
    # Only remove such tokens when there is more than one part; this preserves
    # legitimate short names that consist of a single token.
    if not preserve_subgroup:
        while len(parts) > 1 and re.fullmatch(r'[A-Za-z]', parts[-1]):
            parts.pop()
    # rejoin using underscore for stable formatting
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

        # ── Domain shortcut keywords ──────────────────────────────────────
        if lospec in ('archaea', 'bacteria'):
            domain_key = lospec          # 'archaea' or 'bacteria'
            group_label = lospec.capitalize()
            for ph in all_phyla:
                dom = ph_to_domain.get(ph, '').lower()
                if domain_key in dom:
                    phylum_to_group[ph] = group_label
            continue

        # ── Explicit spec: optional "Label:phyla,list" or just "phyla,list" ─
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
            # Try normalized matching against known phyla
            actual = norm_to_actual.get(_norm(requested))
            if actual:
                phylum_to_group[actual] = label
            else:
                # Phylum not yet in dataset — store anyway; may arrive later
                phylum_to_group[requested] = label

    return phylum_to_group


def write_dataset_colorstrip(output_path: str, dataset_label: str, id_to_color: dict, legend_title: str = None):
    """Write a simple iTOL DATASET_COLORSTRIP file.

    Parameters
    ----------
    output_path:
        Destination path for the .itol file.
    dataset_label:
        DATASET_LABEL value.
    id_to_color:
        Mapping of sequence id -> hex color.
    legend_title:
        Optional legend title. When omitted, a simple legend is derived from
        the unique colors used in id_to_color.
    """
    p = Path(output_path)
    p.parent.mkdir(parents=True, exist_ok=True)
    unique_items = []
    seen = set()
    for iid, color in id_to_color.items():
        key = (str(color),)
        if key not in seen:
            seen.add(key)
            unique_items.append((iid, color))

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
        lines.append('LEGEND_COLORS,' + ','.join(color for _, color in unique_items))
        lines.append('LEGEND_LABELS,' + ','.join(str(iid).replace(',', ';') for iid, _ in unique_items))
    lines.append('DATA')
    for iid, color in id_to_color.items():
        lines.append(f'{iid},{color}')
    p.write_text('\n'.join(lines) + '\n')
    return str(p)


def build_dataset_color_map(dataset_names):
    dataset_names = sorted({d for d in dataset_names if d})
    QUAL_PALETTE = ['#1f78b4', '#e31a1c', '#33a02c', '#ff7f00', '#6a3d9a', '#b15928', '#a6cee3', '#fb9a99', '#b2df8a', '#fdbf6f', '#cab2d6', '#ffff99']
    ds_color_map = {}
    if len(dataset_names) <= len(QUAL_PALETTE):
        for i, dn in enumerate(dataset_names):
            ds_color_map[dn] = QUAL_PALETTE[i]
    else:
        for dn in dataset_names:
            ds_color_map[dn] = _name_to_dataset_color(dn, bins=max(12, len(dataset_names)))
    return ds_color_map


def write_dataset_membership_strip(output_path: str, ids_in_order, ds_map, dataset_label: str = 'Dataset membership', other_color: str = '#cccccc'):
    ds_color_map = build_dataset_color_map(ds_map.values())
    lines = [
        'DATASET_COLORSTRIP',
        'SEPARATOR COMMA',
        f'DATASET_LABEL,{dataset_label}',
        'COLOR,#AAAAAA',
        'MARGIN,5',
        'SHOW_INTERNAL,0',
    ]
    dataset_names = list(ds_color_map.keys())
    if dataset_names:
        lines.append('LEGEND_TITLE,Dataset membership legend')
        lines.append('LEGEND_SHAPES,' + ','.join(['1'] * len(dataset_names)))
        lines.append('LEGEND_COLORS,' + ','.join(ds_color_map[d] for d in dataset_names))
        lines.append('LEGEND_LABELS,' + ','.join(d.replace(',', ';') for d in dataset_names))
    lines.append('DATA')
    for iid in ids_in_order:
        ds_name = ds_map.get(iid, '')
        lines.append(f"{iid},{ds_color_map.get(ds_name, other_color)}")
    p = Path(output_path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text('\n'.join(lines) + '\n')
    return str(p)


def generate_itol_colors(taxonomy_tsv: str, outdir: str, user_color_csv: str = None, id_map: dict = None, tree_file: str = None, phylum_groups: list = None):
    """Generate one iTOL-compatible colorstrip per taxonomy metadata type.

    Outputs DATASET_COLORSTRIP files in outdir for phylum, family, genus, and
    optional user colours. TREE_COLORS and DATASET_SYMBOL variants are not
    retained because they duplicate the same metadata in another visual
    encoding and clutter the workflow outputs.

    The user_color_csv, if provided, should be a CSV with header including 'id' and 'color'.

    phylum_groups, if provided, is a list of grouping spec strings passed from
    ``--group-phyla``.  See ``_parse_phylum_groups`` for the spec format.
    Examples::

        phylum_groups=['archaea']
        phylum_groups=['Bacillota,Bacillota_I,Bacillota_A']
        phylum_groups=['archaea', 'Firmicutes:Bacillota,Bacillota_I']
    """
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    # If no id_map provided, attempt to discover one near the taxonomy file
    # or in the output directory. This helps map long/original headers to the
    # short IDs used in the tree when callers forget to pass an id_map.
    if id_map is None:
        id_map = {}
        try:
            tdir = Path(taxonomy_tsv).parent
            # candidate files to look for
            candidates = []
            for cand_name in ("preload_id_map.tsv", "user_id_map.tsv", "user_id_map.csv"):
                p = tdir / cand_name
                if p.exists():
                    candidates.append(p)
            # any *_id_map.tsv files in the same dir
            try:
                for p in tdir.glob('*_id_map.tsv'):
                    candidates.append(p)
            except Exception:
                pass
            # also check output dir for user_id_map
            try:
                for cand_name in ("preload_id_map.tsv", "user_id_map.tsv"):
                    p = out / cand_name
                    if p.exists():
                        candidates.append(p)
            except Exception:
                pass

            # read mapping files (expect short\torig columns with optional header)
            for p in candidates:
                try:
                    with open(p) as fh:
                        # skip header if present
                        first = fh.readline()
                        if first and ('short' in first.lower() or 'original' in first.lower() or '\t' in first):
                            # check if first line looks like a header by presence of letters
                            # if it contains a tab and non-numeric text, treat as header
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
                        for l in fh:
                            parts = l.strip().split('\t')
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

    # Helper to map any input id to the short id used in the tree, when an
    # id_map dict is provided. Mapping attempts several reasonable fallbacks
    # (exact, case-insensitive, last pipe-field, substring matches) to be
    # resilient to common header transformations.
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
            # accept either classification-style rows (ID,BestHit,Identity,Taxon,Confidence)
            # or combined-style rows (ID,Taxon,Confidence). Support both to avoid losing
            # taxonomy when using the combined files produced by preload.
            tax = 'NA'
            if len(parts) > 3:
                tax = parts[3]
            elif len(parts) > 1:
                tax = parts[1]
            mapped_qid = _map_id(qid)
            # if multiple input ids map to the same short id, prefer the first
            # seen taxonomy (do not overwrite an existing assignment)
            if mapped_qid not in taxa:
                taxa[mapped_qid] = tax

    user_colors = {}
    if user_color_csv:
        import csv
        with open(user_color_csv) as uc:
            reader = csv.DictReader(uc)
            for row in reader:
                rid = row.get('id') or row.get('ID') or row.get('Id')
                col = row.get('color') or row.get('colour') or row.get('Color')
                if rid and col:
                    # map user-provided color ids to short ids when possible
                    mapped = rid
                    if id_map:
                        mapped = _map_id(rid)
                    user_colors[mapped] = col

    combined_lines = []

    # first pass: parse all taxon strings and collect unique names per rank
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
        # normalize names to remove accession-like numeric suffixes and
        # collapse redundant variants produced by gg2 (e.g. Bacillota_A_368345)
        # preserve single-letter subgroup tokens for phylum names (avoid
        # collapsing Bacillota_A -> Bacillota). Family/genus continue to be
        # normalized with subgroup collapsing.
        ph = _normalize_taxon_name(ph, preserve_subgroup=True)
        fa = _normalize_taxon_name(fa)
        ge = _normalize_taxon_name(ge)
        parsed_map[qid] = (ph, fa, ge)
        ph_names[ph] = ph_names.get(ph, 0) + 1
        fa_names[fa] = fa_names.get(fa, 0) + 1
        ge_names[ge] = ge_names.get(ge, 0) + 1
        if ph and ph != 'unknown' and dom:
            ph_to_domain.setdefault(ph, dom)

    # ── Apply phylum grouping ─────────────────────────────────────────────────
    # Build a {phylum: group_label} map from --group-phyla specs.  Phyla not
    # mentioned in any spec keep their own name as the effective label.
    phylum_to_group = _parse_phylum_groups(
        phylum_groups or [], set(ph_names.keys()), ph_to_domain
    )

    if phylum_to_group:
        # Remap parsed_map so grouped phyla share a single label.
        # Family and genus are intentionally left ungrouped.
        parsed_map = {
            qid: (phylum_to_group.get(ph, ph), fa, ge)
            for qid, (ph, fa, ge) in parsed_map.items()
        }
        # Merge counts: grouped phyla collapse into a single entry
        group_ph_names: dict = {}
        for ph, count in ph_names.items():
            grp = phylum_to_group.get(ph, ph)
            group_ph_names[grp] = group_ph_names.get(grp, 0) + count
    else:
        group_ph_names = ph_names

    # generate distinct palettes per rank by evenly spacing hues across unique names
    def _make_palette(names, rank):
        names_list = sorted(names.keys())
        n = len(names_list) if names_list else 1
        cmap = {}

        # ColorBrewer qualitative palettes for small numbers of categories
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

        # If we have a small number of categories, prefer a ColorBrewer palette
        if n in BREWER:
            palette = BREWER[n]
            for name, col in zip(names_list, palette):
                cmap[name] = col
            return cmap

        # For phylum-level palettes we want very distinct, high-contrast
        # colours even when there are >12 categories. Use a golden-ratio-based
        # hue stepping to avoid producing adjacent-similar hues when many
        # categories are present. Use stronger saturation/value for phylum to
        # increase perceptual separation.
        if rank == 'phylum':
            # Prefer a curated high-contrast palette (Kelly 22) when possible;
            # otherwise extend it by greedily sampling candidates maximizing
            # perceptual separation in Lab colour space (Delta E 1976 approx).
            # For robust distinctness across many phyla prefer to sample a
            # dense set of candidate colours and then greedily pick the set
            # that maximises minimal perceptual distance (Lab). Do not
            # pre-seed with Kelly since some Kelly entries are similar
            # (e.g. two yellows) and reduce minimal separation.
            import colorsys

            # build candidate pool across H, S, V with several levels
            candidates = []
            for hdeg in range(0, 360, 6):
                h = (hdeg % 360) / 360.0
                for s in (0.95, 0.85, 0.75, 0.65):
                    for v in (0.96, 0.90, 0.80, 0.70):
                        r, g, b = colorsys.hsv_to_rgb(h, s, v)
                        rgb = (int(r*255), int(g*255), int(b*255))
                        hexc = '#{:02x}{:02x}{:02x}'.format(*rgb)
                        candidates.append((hexc, rgb))

            # unique candidates preserving order
            seen = set()
            uniq = []
            for hexc, rgb in candidates:
                if hexc not in seen:
                    seen.add(hexc)
                    uniq.append((hexc, rgb))

            # helper: convert hex/rgb to Lab
            def _rgb_to_lab(rgb):
                r,g,b=[x/255.0 for x in rgb]
                def lin(c):
                    if c<=0.04045: return c/12.92
                    return ((c+0.055)/1.055)**2.4
                r,g,b=lin(r),lin(g),lin(b)
                X=r*0.4124564+g*0.3575761+b*0.1804375
                Y=r*0.2126729+g*0.7151522+b*0.0721750
                Z=r*0.0193339+g*0.1191920+b*0.9503041
                Xn,Yn,Zn=0.95047,1.0,1.08883
                x,y,z=X/Xn,Y/Yn,Z/Zn
                def f(t):
                    if t>0.008856: return t**(1/3)
                    return 7.787*t+16/116
                fx,fy,fz=f(x),f(y),f(z)
                L=116*fy-16
                a=500*(fx-fy)
                bval=200*(fy-fz)
                return (L,a,bval)

            cand_lab = [(_hex, _rgb_to_lab(rgb)) for (_hex, rgb) in uniq]

            # greedy farthest-point sampling in Lab space
            chosen = []
            chosen_lab = []
            # deterministic seed index based on concatenated names
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

            # assign chosen colours to names
            for name, (hexc, _) in zip(names_list, chosen):
                cmap[name] = hexc

            # fallback: evenly space if not enough candidates
            if len(cmap) < len(names_list):
                remaining = [nm for nm in names_list if nm not in cmap]
                m = len(remaining)
                for i, name in enumerate(remaining):
                    h = (i / float(m)) % 1.0
                    cmap[name] = _hsv_to_hex(h, s=0.85, v=0.95)
            return cmap

        # Otherwise evenly space hues across the wheel to ensure distinct colors
        # Apply a rank-specific offset so ranks have different base palettes
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

    # If a tree file is provided, repair/write stable internal node labels and
    # retain the newick text so colorstrip rows can be limited to actual leaves.
    newick = None
    if tree_file:
        try:
            tf = Path(tree_file)
            if tf.exists():
                newick = tf.read_text()
                try:
                    # Repair malformed NODE labels before further tree rewriting.
                    from branchmanager.pipeline.tree import _repair_internal_node_label_delimiters
                    newick = _repair_internal_node_label_delimiters(newick)
                except Exception:
                    pass
                # Inject deterministic textual node names for any internal node
                # that has a numeric-only support value so iTOL can match them.
                labeled_newick = newick
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

                    labeled_newick = re.sub(r"\)\s*([0-9]+(?:\.[0-9]+)?)", _repl, labeled_newick)
                except Exception:
                    labeled_newick = newick
                # Write the labeled tree back to the same file.
                try:
                    tf.write_text(labeled_newick)
                    newick = labeled_newick
                except Exception:
                    try:
                        out_labeled = out / 'current_tree_labeled.nwk'
                        out_labeled.write_text(labeled_newick)
                        newick = labeled_newick
                    except Exception:
                        pass
        except Exception:
            newick = None

    # second pass: build per-sequence combined lines using per-name palettes,
    # allowing a per-sequence user color to override if provided
    for qid in list(taxa.keys()):
        ph, fa, ge = parsed_map[qid]
        user_col = user_colors.get(qid, '')
        ph_c = user_col or phylum_map.get(ph, _name_to_color_by_rank(ph, 'phylum'))
        fa_c = user_col or family_map.get(fa, _name_to_color_by_rank(fa, 'family'))
        ge_c = user_col or genus_map.get(ge, _name_to_color_by_rank(ge, 'genus'))
        combined_lines.append((qid, ph_c, fa_c, ge_c, user_col))


    # helper to write iTOL DATASET_COLORSTRIP file
    def write_colorstrip(path, title, id_color_pairs, legend_pairs=None):
        """Write a minimal iTOL DATASET_COLORSTRIP file.

        legend_pairs: optional list of (label, color) to include in the legend.
        """
        with open(path, 'w') as f:
            f.write('DATASET_COLORSTRIP\n')
            f.write('SEPARATOR COMMA\n')
            f.write(f'DATASET_LABEL,{title}\n')
            # set a dataset-level COLOR (use first legend color if available,
            # otherwise fall back to the first id color or a neutral default)
            if legend_pairs and len(legend_pairs) > 0:
                dataset_color = legend_pairs[0][1]
            else:
                dataset_color = id_color_pairs[0][1] if id_color_pairs else '#AAAAAA'
            f.write(f'COLOR,{dataset_color}\n')
            f.write('MARGIN,5\n')
            f.write('SHOW_INTERNAL,0\n')

            # optional legend block
            if legend_pairs:
                # legend labels may contain commas; replace them to avoid CSV
                labels = [lbl.replace(',', ';') for lbl, _ in legend_pairs]
                colors = [col for _, col in legend_pairs]
                shapes = ['1'] * len(legend_pairs)
                f.write(f"LEGEND_TITLE,{title} legend\n")
                f.write('LEGEND_SHAPES,' + ','.join(shapes) + '\n')
                f.write('LEGEND_COLORS,' + ','.join(colors) + '\n')
                f.write('LEGEND_LABELS,' + ','.join(labels) + '\n')

            f.write('DATA\n')
            for item in id_color_pairs:
                # accept either (id, color) or (id, color, name)
                if len(item) >= 2:
                    id_, col = item[0], item[1]
                else:
                    continue
                f.write(f"{id_},{col}\n")

    # Determine set of leaf ids present in the provided tree (when available).
    # This lets us avoid writing colorstrip rows for ids that are not present
    # in the uploaded tree.
    leaf_ids_in_tree = set()
    if tree_file:
        try:
            import re as _re
            # find candidate leaf tokens (text preceding :, ,, ) or ;)
            if not newick:
                raise ValueError('no newick content')
            tokens = _re.findall(r"([^(),:;]+)(?=[:),;])", newick)
            for t in tokens:
                tt = t.strip()
                if not tt:
                    continue
                # map token to short id when possible (use same mapper used above)
                mapped = _map_id(tt)
                leaf_ids_in_tree.add(mapped)
        except Exception:
            leaf_ids_in_tree = set()

    # write per-rank color maps as iTOL colorstrip datasets
    # for per-rank we need per-sequence mappings: map each qid to its phylum/family/genus color
    # include the taxon name in pairs for stable tree-id resolution
    ph_pairs = [(qid, ph_c, parsed_map[qid][0]) for qid, ph_c, _, _, _ in combined_lines]
    fa_pairs = [(qid, fa_c, parsed_map[qid][1]) for qid, _, fa_c, _, _ in combined_lines]
    ge_pairs = [(qid, ge_c, parsed_map[qid][2]) for qid, _, _, ge_c, _ in combined_lines]

    # Only emit rows for ids that actually appear in the tree. Build filtered
    # versions of the pairs for colorstrip output.
    if leaf_ids_in_tree:
        # Resolve each qid to an id present in the tree. Load member->rep
        # collapse mappings so original members map to their representative
        # IDs in the (possibly collapsed) tree.
        # Only two-column files are treated as member->rep; three-column
        # rep_map files (rep\ttax\tcount) are ignored.
        collapse_map = {}
        try:
            cand_paths = [
                out / 'preload_collapsed_members.tsv',
                out / 'collapsed_members.tsv',
                out / 'preload_collapsed_map.tsv',
                out / 'collapsed_map.tsv',
                Path(taxonomy_tsv).parent / 'preload_collapsed_members.tsv',
                Path(taxonomy_tsv).parent / 'preload_collapsed_map.tsv'
            ]
            for cp in cand_paths:
                try:
                    if not cp.exists():
                        continue
                    with open(cp) as cf:
                        next(cf, None)
                        for l in cf:
                            if not l.strip():
                                continue
                            parts = l.strip().split('\t')
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
            # last pipe field
            if '|' in x:
                forms.append(x.split('|')[-1])
            # whitespace-separated last token
            forms.append(x.split()[-1] if x.split() else x)
            # uppercase/lowercase variations
            try:
                forms.append(forms[0].upper())
                forms.append(forms[0].lower())
            except Exception:
                pass
            # strip trailing underscore+digits (collapse suffixes)
            import re as _re
            m = _re.sub(r'_n?\d+$','', forms[0] if forms else x)
            forms.append(m)
            # unique preserving order
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
            # write missing ids summary for diagnostics
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
        # no tree or failed parse -> assume all ids may be present
        ph_pairs_tree = ph_pairs
        fa_pairs_tree = fa_pairs
        ge_pairs_tree = ge_pairs

    # prepare legend pairs (label,color) for each rank
    # prepare legend pairs including counts e.g. 'Firmicutes (12)'
    ph_legend = sorted([(name + f" ({group_ph_names.get(name,0)})", col) for name, col in phylum_map.items()], key=lambda x: x[0])
    fa_legend = sorted([(name + f" ({fa_names.get(name,0)})", col) for name, col in family_map.items()], key=lambda x: x[0])
    ge_legend = sorted([(name + f" ({ge_names.get(name,0)})", col) for name, col in genus_map.items()], key=lambda x: x[0])

    # When a tree is provided prefer to write dataset files using the actual
    # leaf ids present in the tree (ph_pairs_tree etc.) so iTOL datasets
    # apply correctly to the uploaded tree. Fall back to the original pairs
    # when no tree mapping was available.
    if leaf_ids_in_tree:
        ph_dataset_pairs = [(qid, col) for qid, col, _ in ph_pairs_tree]
        fa_dataset_pairs = [(qid, col) for qid, col, _ in fa_pairs_tree]
        ge_dataset_pairs = [(qid, col) for qid, col, _ in ge_pairs_tree]
    else:
        ph_dataset_pairs = [(qid, col) for qid, col, _ in ph_pairs]
        fa_dataset_pairs = [(qid, col) for qid, col, _ in fa_pairs]
        ge_dataset_pairs = [(qid, col) for qid, col, _ in ge_pairs]

    write_colorstrip(out / 'itol_phylum_colors.itol', 'Phylum colors', ph_dataset_pairs, legend_pairs=ph_legend)
    write_colorstrip(out / 'itol_family_colors.itol', 'Family colors', fa_dataset_pairs, legend_pairs=fa_legend)
    write_colorstrip(out / 'itol_genus_colors.itol', 'Genus colors', ge_dataset_pairs, legend_pairs=ge_legend)

    # write user-provided colors as iTOL dataset if any
    if user_colors:
        user_pairs = [(rid, col) for rid, col in user_colors.items()]
        write_colorstrip(out / 'itol_user_colors.itol', 'User colors', user_pairs)

    # Keep only one iTOL metadata file per metadata type: the colorstrip.
    # Branch TREE_COLORS files and DATASET_SYMBOL files are alternative
    # encodings of the same metadata and made workflow outputs noisy.
    for stale_name in (
        'itol_phylum_tree_colors.txt',
        'itol_family_tree_colors.txt',
        'itol_genus_tree_colors.txt',
        'itol_phylum_symbols.itol',
        'itol_family_symbols.itol',
        'itol_genus_symbols.itol',
        'itol_combined_colors.csv',
        'tree_colors_with_clades.txt',
    ):
        try:
            stale = out / stale_name
            if stale.exists():
                stale.unlink()
        except Exception:
            pass

    return str(out)


# ── Functional annotations ────────────────────────────────────────────────────

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

    # ── ID mapper ────────────────────────────────────────────────────────────
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

    # ── Read TSV ─────────────────────────────────────────────────────────────
    with open(func_tsv, newline='') as fh:
        reader = csv.reader(fh, delimiter='\t')
        header = next(reader, None)
        if not header or len(header) < 2:
            return []
        id_col  = header[0]
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

        # Sanitise column name for file name
        safe_name = _re.sub(r'[^\w]+', '_', col_name).strip('_') or f'col{col_idx+1}'

        # ── Detect column type ───────────────────────────────────────────────
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

        # ── DATASET_BINARY ────────────────────────────────────────────────────
        if is_binary:
            # Normalise to 0 / 1
            def _to_01(v):
                lv = v.lower()
                if lv in ('1', 'yes', 'true', '+'):
                    return '1'
                return '0'

            col_color = _name_to_color(col_name)
            lines_out = [
                'DATASET_BINARY',
                'SEPARATOR COMMA',
                f'DATASET_LABEL,{col_name}',
                f'COLOR,{col_color}',
                f'FIELD_LABELS,{col_name}',
                f'FIELD_COLORS,{col_color}',
                f'FIELD_SHAPES,1',
                'SHOW_INTERNAL,0',
                'DATA',
            ]
            for row_id, vals in rows:
                v = vals[col_idx]
                lines_out.append(f'{row_id},{_to_01(v)}')
            dest.write_text('\n'.join(lines_out) + '\n')
            written.append(str(dest))

        # ── DATASET_SIMPLEBAR ─────────────────────────────────────────────────
        elif is_numeric:
            col_color = _name_to_color(col_name)
            try:
                float_vals = [float(v) for v in non_empty]
                max_val = max(float_vals)
            except Exception:
                max_val = 1.0
            lines_out = [
                'DATASET_SIMPLEBAR',
                'SEPARATOR COMMA',
                f'DATASET_LABEL,{col_name}',
                f'COLOR,{col_color}',
                f'WIDTH,200',
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

        # ── DATASET_COLORSTRIP (categorical) ──────────────────────────────────
        else:
            # Build a per-category colour palette
            categories = []
            seen_cats = set()
            for v in col_values:
                lv = v.lower()
                if lv not in _EMPTY_VALUES and v not in seen_cats:
                    seen_cats.add(v)
                    categories.append(v)

            cat_palette = _make_palette_simple(categories)

            # Legend entries
            legend_shapes  = ','.join(['1'] * len(categories))
            legend_colors  = ','.join(cat_palette[c] for c in categories)
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
                f'LEGEND_COLORS,{legend_colors}',
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
    """Return a ``{name: hex_color}`` palette for a list of category names.

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
    # Fall back to evenly spaced hues
    return {name: _hsv_to_hex(i / float(n), s=0.80, v=0.92) for i, name in enumerate(names)}



# ── Rumen bacterial functional groups ────────────────────────────────────────

# Broad-scale ruminant microbiome functional group assignments.
# Keys are lowercase genus, family, order, or phylum names (or partial matches).
# Values are human-readable functional group labels.
# Hierarchy: genus match takes priority over family, then order, then phylum.
# A sequence not matching any entry defaults to 'Other / Unclassified'.

# Curated at genus level first (most specific)
_RUMEN_GENUS: dict = {
    # Cellulolytic / Fibrolytic
    'fibrobacter': 'Cellulolytic/Fibrolytic',
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

# Ordered palette for the known functional groups (for consistent, distinct colours)
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
    phylum = _lo(_normalize_taxon_name(parsed.get('p', ''), preserve_subgroup=True))
    family = _lo(_normalize_taxon_name(parsed.get('f', '')))
    genus  = _lo(_normalize_taxon_name(parsed.get('g', '')))

    # Archaea domain → Methanogenic Archaea (many rumen archaea are methanogens)
    if 'archaea' in domain:
        # Still allow a genus-level override first
        hit = _RUMEN_GENUS.get(genus)
        if hit:
            return hit
        return 'Methanogenic Archaea'

    # Genus → family → phylum hierarchy
    for lookup, key in (
        (_RUMEN_GENUS,  genus),
        (_RUMEN_FAMILY, family),
        (_RUMEN_PHYLUM, phylum),
    ):
        # exact match
        hit = lookup.get(key)
        if hit:
            return hit
        # partial / prefix match (handles subgroup suffixes like 'bacillota_a')
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

    Reads the combined taxonomy TSV produced by a BranchManager run and maps each
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

    # ── ID mapper ────────────────────────────────────────────────────────────
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

    rows: list = []   # [(short_id, func_group)]
    try:
        # Support gzipped taxonomy files
        import gzip as _gzip
        _opener = _gzip.open if str(taxonomy_tsv).endswith(('.gz', '.gzip')) else open
        with _opener(taxonomy_tsv, 'rt') as fh:
            reader = _csv.reader(fh, delimiter='\t')
            header = next(reader, None)  # skip header
            for line in reader:
                if not line:
                    continue
                qid = line[0].strip()
                # Accept both 5-column and 3-column taxonomy TSV formats
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

    # ── Write draft TSV ──────────────────────────────────────────────────────
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

    # ── Write iTOL DATASET_COLORSTRIP ────────────────────────────────────────
    # Build category palette: prefer curated colours, generate for any unknowns.
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
    # assign distinct colours to any unrecognised groups
    if unknown_cats:
        extra = _make_palette_simple(unknown_cats)
        cat_palette.update(extra)

    legend_shapes = ','.join(['1'] * len(cats_seen))
    legend_colors = ','.join(cat_palette.get(c, '#cccccc') for c in cats_seen)
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
        f'LEGEND_COLORS,{legend_colors}',
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
