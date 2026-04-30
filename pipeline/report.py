def generate_report(derep, taxonomy, novelty, outdir, similarity_pivot=None):
    out = f"{outdir}/final_report.tsv"

    tax_map = {}
    with open(taxonomy) as f:
        for line in f:
            parts = line.strip().split("\t")
            tax_map[parts[0]] = parts[3] if len(parts) > 3 else "NA"

    nov_map = {}
    with open(novelty) as f:
        next(f)
        for line in f:
            i, n = line.strip().split("\t")
            nov_map[i] = n

    # load similarity pivot if provided
    pivot_headers = []
    pivot_map = {}
    if similarity_pivot:
        try:
            with open(similarity_pivot) as sp:
                hdr = sp.readline().strip().split("\t")
                # first header should be 'id'
                pivot_headers = hdr[1:]
                for line in sp:
                    parts = line.strip().split("\t")
                    if not parts:
                        continue
                    uid = parts[0]
                    pivot_map[uid] = parts[1:]
        except Exception:
            pivot_headers = []
            pivot_map = {}

    with open(out, "w") as f:
        # write base headers
        base_hdr = ["ID", "Taxonomy", "Novel"]
        # append pivot headers if present
        all_hdr = base_hdr + pivot_headers
        f.write("\t".join(all_hdr) + "\n")
        for i in tax_map:
            row = [i, tax_map[i], nov_map.get(i, 'NA')]
            if pivot_headers:
                vals = pivot_map.get(i, ['NA'] * len(pivot_headers))
                row.extend(vals)
            f.write("\t".join(str(x) for x in row) + "\n")
