"""
report.py — Final summary report generation for BranchManager.
"""
import logging

logger = logging.getLogger(__name__)


def generate_report(derep, taxonomy, novelty, outdir, similarity_pivot=None):
    """Generate a final summary TSV combining taxonomy, novelty and optional similarity data."""
    out = f"{outdir}/final_report.tsv"

    tax_map = {}
    try:
        with open(taxonomy) as f:
            next(f, None)  # skip header if present
            for line in f:
                parts = line.strip().split("\t")
                if not parts or not parts[0]:
                    continue
                # taxonomy.tsv: ID  BestHit  Identity  Taxon  Confidence
                tax_map[parts[0]] = parts[3] if len(parts) > 3 else "NA"
    except FileNotFoundError:
        logger.warning("[REPORT] Taxonomy file not found: %s", taxonomy)

    nov_map = {}
    try:
        with open(novelty) as f:
            next(f, None)  # skip header
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 2:
                    continue
                # novelty.tsv: ID  NearestIdentity  NearestHit  Novel
                nov_map[parts[0]] = parts[3] if len(parts) > 3 else parts[1]
    except FileNotFoundError:
        logger.warning("[REPORT] Novelty file not found: %s", novelty)

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
        except Exception as e:
            logger.warning("[REPORT] Failed to load similarity pivot %s: %s", similarity_pivot, e)
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

    logger.info("[REPORT] Wrote final report to %s", out)
    return out
