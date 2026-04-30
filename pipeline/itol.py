from pathlib import Path
import hashlib


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
    # backward compatible simple mapping
    h = _hash_to_hue(name)
    return _hsv_to_hex(h)


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

    # choose saturation/value per rank for additional separation
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
    # taxon expected like: d__Bacteria; p__Firmicutes; c__Bacilli; ...
    parts = [p.strip() for p in taxon.split(';')]
    res = {}
    for p in parts:
        if '__' in p:
            k, v = p.split('__', 1)
            k = k.strip()
            v = v.strip()
            if v:
                res[k] = v
    return res


def generate_itol_colors(taxonomy_tsv: str, outdir: str, user_color_csv: str = None, id_map: dict = None):
    """Generate simple iTOL-compatible color mapping files for phylum, family, genus.

    Outputs files in outdir:
      - itol_phylum_colors.csv (id,color)
      - itol_family_colors.csv
      - itol_genus_colors.csv
      - itol_user_colors.csv (if provided)
      - itol_combined_colors.csv (id,phylum_color,family_color,genus_color,user_color)

    The user_color_csv, if provided, should be a CSV with header including 'id' and 'color'.
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
            for cand_name in ("user_id_map.tsv", "user_id_map.csv"):
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
                for cand_name in ("user_id_map.tsv",):
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
    for qid, taxstr in taxa.items():
        parsed = parse_taxon_string(taxstr)
        ph = parsed.get('p', 'unknown')
        fa = parsed.get('f', 'unknown')
        ge = parsed.get('g', 'unknown')
        parsed_map[qid] = (ph, fa, ge)
        ph_names[ph] = ph_names.get(ph, 0) + 1
        fa_names[fa] = fa_names.get(fa, 0) + 1
        ge_names[ge] = ge_names.get(ge, 0) + 1

    # generate distinct palettes per rank by evenly spacing hues across unique names
    def _make_palette(names, rank):
        names_list = sorted(names.keys())
        n = len(names_list) if names_list else 1
        cmap = {}

        # ColorBrewer qualitative palettes (hex) for small numbers of categories
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

        # Use deterministic per-name colors so the same taxon receives the same
        # color across different runs and datasets. Rely on _name_to_color_by_rank
        # which hashes the taxon name and applies a rank-specific offset.
        for name in names_list:
            cmap[name] = _name_to_color_by_rank(name, rank)
        return cmap

    phylum_map = _make_palette(ph_names, 'phylum')
    family_map = _make_palette(fa_names, 'family')
    genus_map = _make_palette(ge_names, 'genus')

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

    # write per-rank color maps as iTOL colorstrip datasets
    # for per-rank we need per-sequence mappings: map each qid to its phylum/family/genus color
    # include the taxon name in pairs to aid symbol assignment later
    ph_pairs = [(qid, ph_c, parsed_map[qid][0]) for qid, ph_c, _, _, _ in combined_lines]
    fa_pairs = [(qid, fa_c, parsed_map[qid][1]) for qid, _, fa_c, _, _ in combined_lines]
    ge_pairs = [(qid, ge_c, parsed_map[qid][2]) for qid, _, _, ge_c, _ in combined_lines]

    # prepare legend pairs (label,color) for each rank
    # prepare legend pairs including counts e.g. 'Firmicutes (12)'
    ph_legend = sorted([(name + f" ({ph_names.get(name,0)})", col) for name, col in phylum_map.items()], key=lambda x: x[0])
    fa_legend = sorted([(name + f" ({fa_names.get(name,0)})", col) for name, col in family_map.items()], key=lambda x: x[0])
    ge_legend = sorted([(name + f" ({ge_names.get(name,0)})", col) for name, col in genus_map.items()], key=lambda x: x[0])

    write_colorstrip(out / 'itol_phylum_colors.itol', 'Phylum colors', ph_pairs, legend_pairs=ph_legend)
    write_colorstrip(out / 'itol_family_colors.itol', 'Family colors', fa_pairs, legend_pairs=fa_legend)
    write_colorstrip(out / 'itol_genus_colors.itol', 'Genus colors', ge_pairs, legend_pairs=ge_legend)

    # write symbol strips to give an additional visual cue (shapes per name)
    def write_symbolstrip(path, title, id_shape_color_triplets, legend_pairs=None):
        with open(path, 'w') as f:
            f.write('DATASET_SYMBOL\n')
            f.write('SEPARATOR COMMA\n')
            f.write(f'DATASET_LABEL,{title}\n')
            # choose a neutral dataset color
            dataset_color = legend_pairs[0][1] if legend_pairs and len(legend_pairs) > 0 else (id_shape_color_triplets[0][2] if id_shape_color_triplets else '#AAAAAA')
            f.write(f'COLOR,{dataset_color}\n')
            f.write('MARGIN,5\n')
            f.write('SHOW_INTERNAL,0\n')
            if legend_pairs:
                labels = [lbl.replace(',', ';') for lbl, _ in legend_pairs]
                colors = [col for _, col in legend_pairs]
                # shapes: cycle through 1..12
                shapes = [str((i % 12) + 1) for i in range(len(legend_pairs))]
                f.write(f"LEGEND_TITLE,{title} legend\n")
                f.write('LEGEND_SHAPES,' + ','.join(shapes) + '\n')
                f.write('LEGEND_COLORS,' + ','.join(colors) + '\n')
                f.write('LEGEND_LABELS,' + ','.join(labels) + '\n')
            f.write('DATA\n')
            for item in id_shape_color_triplets:
                # accept triplet (id, shape, color)
                if len(item) >= 3:
                    id_, shape, col = item[0], item[1], item[2]
                else:
                    continue
                f.write(f"{id_},{shape},{col}\n")

    # prepare symbol triplets for each rank: id, shape, color
    # choose shapes by mapping unique names to shape numbers 1..12
    def _make_symbol_triplets(pairs, legend_map):
        # pairs: list of (qid, color, tax_name)
        # legend_map: list of (name_label, color)
        shapes = {}
        for i, (lbl, _) in enumerate(legend_map):
            shapes[lbl] = str((i % 12) + 1)
        triplets = []
        for qid, col, tax_name in pairs:
            # legend_map labels are like 'Name (count)'; extract base name
            matched_shape = None
            for lbl, _ in legend_map:
                name = lbl.split(' (')[0]
                if name == tax_name:
                    matched_shape = shapes.get(lbl)
                    break
            if not matched_shape:
                matched_shape = str((hash(qid) & 0x0f) + 1)
            triplets.append((qid, matched_shape, col))
        return triplets

    ph_symbol_triplets = _make_symbol_triplets(ph_pairs, ph_legend)
    fa_symbol_triplets = _make_symbol_triplets(fa_pairs, fa_legend)
    ge_symbol_triplets = _make_symbol_triplets(ge_pairs, ge_legend)

    write_symbolstrip(out / 'itol_phylum_symbols.itol', 'Phylum symbols', ph_symbol_triplets, legend_pairs=ph_legend)
    write_symbolstrip(out / 'itol_family_symbols.itol', 'Family symbols', fa_symbol_triplets, legend_pairs=fa_legend)
    write_symbolstrip(out / 'itol_genus_symbols.itol', 'Genus symbols', ge_symbol_triplets, legend_pairs=ge_legend)

    # write combined per-node CSV for convenience
    with open(out / 'itol_combined_colors.csv', 'w') as f:
        f.write('id,phylum_color,family_color,genus_color,user_color\n')
        for row in combined_lines:
            f.write(','.join(row) + '\n')

    # write user-provided colors as iTOL dataset if any
    if user_colors:
        user_pairs = [(rid, col) for rid, col in user_colors.items()]
        write_colorstrip(out / 'itol_user_colors.itol', 'User colors', user_pairs)

    return str(out)

