# Relict -  WARNING - BETA! - This tool is in early development. Expect bugs, breaking changes, and rough edges. Please report issues and contribute improvements!

**Reference-aware 16S novelty and phylogenetic context tool.**

Relict helps answer the question:

> Has this lineage already been seen and characterised, or is it still poorly represented enough to justify deeper follow-up such as whole-genome sequencing?

It combines taxonomic classification, nearest-neighbour novelty scoring, neighbourhood density (crowding), and phylogenetic tree visualisation into a single per-sequence assessment file.

---

## Quick start

```bash
# 1. Load a baseline dataset and build a reference tree
relict preload \
  --fasta hungate.fasta \
  --db project.db \
  --dataset Hungate \
  --ref gtdb_r232_bac_arch_ssu_reps.fna \
  --classify --build-tree \
  -o preload_out \
  --threads 10

# 2. Process new sequences — scores novelty against the baseline
relict run \
  --input new_sequences.fasta \
  --db project.db \
  --dataset Batch1 \
  --ref gtdb_r232_bac_arch_ssu_reps.fna \
  --preload-dir preload_out \
  -o batch1_out \
  --threads 10

# 3. Zoom into a specific taxon
relict subtree \
  --db project.db \
  --taxon archaea \
  --from-dir preload_out \
  -o archaea_out

# 4. Regenerate iTOL files with new grouping options
relict regen-itol \
  --db project.db \
  --out preload_out \
  --group-phyla archaea \
  --group-phyla "Bacillota,Bacillota_I,Bacillota_A"
```

---

## Core concepts

### Novelty is relative to YOUR submitted sequences

All novelty and neighbourhood-density calculations are made against the sequences submitted to Relict (the preload database plus all prior `run` datasets stored in the DB). They are **not** made against the full external reference (GTDB/SILVA) unless no user-submitted sequences are present.

This is intentional: Relict helps you judge whether a new sequence is worth investigating relative to what your lab or project has already characterised — not relative to all known biology.

Each successive `run` extends the reference pool, so novelty scores become increasingly precise as your project grows.

### Three-layer tree architecture

The phylogenetic tree is built from three independent layers:

| Layer | Source | Shown in iTOL? |
|---|---|---|
| 1 — Anchors | 26 NCBI RefSeq type strains (bundled) | **No** — invisible topology scaffolding |
| 2 — Preload | Your baseline dataset (e.g. Hungate) | Yes |
| 3 — Run sequences | Each new `relict run` batch | Yes |

Anchor sequences constrain phylum-level topology during MAFFT + FastTree inference, then are pruned from the stored newick. They are never shown in outputs, iTOL files, or novelty scoring.

---

## Subcommands

### `relict preload`

Load a baseline FASTA dataset, classify against a reference, and build the backbone tree. This is always the first step.

```bash
relict preload \
  --fasta baseline.fasta \
  --db project.db \
  --dataset Hungate \
  --ref gtdb_r232_bac_arch_ssu_reps.fna \
  --classify --build-tree \
  -o preload_out \
  --threads 10 \
  --collapse \
  --no-shorten-ids
```

| Parameter | Required | Default | Description |
|---|---|---|---|
| `--fasta` | ✓ | — | Input FASTA file containing the baseline sequences |
| `--db` | ✓ | — | Path to the Relict SQLite database (created if absent) |
| `-o / --out` | | `.` | Output directory for tree, iTOL files, and reports |
| `--dataset` | ✓ | — | Label stored in the DB (e.g. `Hungate`) |
| `--ref` | | — | Reference FASTA (GTDB/SILVA reps) for classification and tree orientation |
| `--classify` | | off | Classify sequences against `--ref` and store taxonomy |
| `--build-tree` | | off | Build the MAFFT + FastTree backbone tree |
| `--taxa` | | — | TSV mapping reference IDs to lineages (optional when `--ref` headers contain GTDB lineages) |
| `--taxa-assignments` | | — | Pre-computed taxonomy for the INPUT sequences (TSV or embedded-lineage FASTA) |
| `--shorten-ids / --no-shorten-ids` | | `--shorten-ids` | Replace input headers with compact IDs (e.g. `HUN001`) or keep source names |
| `--collapse` | | off | Collapse near-identical same-taxonomy sequences into one representative for the tree |
| `--collapse-threshold` | | `99.8` | Identity threshold (%) for collapsing |
| `--kingdom` | | — | Keep only sequences belonging to this domain (e.g. `bacteria`) |
| `--anchors` | | bundled | Custom anchor FASTA for tree topology scaffolding |
| `--threads` | | `4` | CPU threads for MAFFT and VSEARCH |
| `--colors` | | — | CSV mapping sequence IDs to custom hex colours for iTOL (columns: `id`, `color`) |
| `--group-phyla SPEC` | | — | Group phyla into a single colour in iTOL (repeatable; see *Phylum grouping*) |
| `--functional` | | — | TSV file mapping sequence IDs to functional attributes (first column = ID; subsequent columns = attributes). Generates one iTOL file per column (DATASET_BINARY / DATASET_SIMPLEBAR / DATASET_COLORSTRIP). |

---

### `relict run`

Process new sequences against the baseline; score novelty and update the phylogenetic tree.

```bash
relict run \
  --input new_sequences.fasta \
  --db project.db \
  --dataset Batch1 \
  --ref gtdb_r232_bac_arch_ssu_reps.fna \
  --preload-dir preload_out \
  -o batch1_out \
  --threads 10 \
  --collapse
```

| Parameter | Required | Default | Description |
|---|---|---|---|
| `--input` | ✓ | — | FASTA file of new sequences to analyse |
| `--db` | ✓ | — | Path to the Relict SQLite database |
| `-o / --out` | ✓ | — | Output directory (assessment TSV, novelty metrics, tree, iTOL files) |
| `--dataset` | ✓ | — | Label for this batch (used in iTOL dataset-membership strip) |
| `--ref` | | — | Reference FASTA for classification and tree orientation |
| `--preload-dir` | | — | Path to the preload output directory — enables fast incremental tree updates |
| `--taxa` | | — | TSV mapping reference IDs to lineages |
| `--taxa-assignments` | | — | Pre-computed taxonomy for the input sequences |
| `--shorten-ids / --no-shorten-ids` | | `--shorten-ids` | Replace headers with compact IDs or keep source names |
| `--min-len` | | `1200` | Minimum sequence length to retain (bp) |
| `--max-n` | | `5` | Maximum ambiguous (N) bases allowed |
| `--collapse` | | off | Collapse near-identical same-taxonomy sequences for the tree |
| `--collapse-threshold` | | `99.8` | Identity threshold (%) for collapsing |
| `--kingdom` | | — | Keep only sequences from this domain |
| `--phylum` | | — | Filter iTOL output to a specific phylum (does not affect novelty scoring) |
| `--target` | | — | FASTA to measure novelty against instead of the DB |
| `--force-rebuild` | | off | Force a full tree rebuild from scratch |
| `--anchors` | | bundled | Custom anchor FASTA for tree scaffolding |
| `--threads` | | `4` | CPU threads |
| `--user-colors` | | — | CSV mapping sequence IDs to custom hex colours for iTOL |
| `--group-phyla SPEC` | | — | Group phyla into a single colour in iTOL (repeatable) |
| `--functional` | | — | TSV file mapping sequence IDs to functional attributes (first column = ID; subsequent columns = attributes). Generates one iTOL file per column (DATASET_BINARY / DATASET_SIMPLEBAR / DATASET_COLORSTRIP). |

---

### `relict subtree`

Extract all sequences matching a given taxon from the DB and build a focused phylogenetic tree for that group only.

**Fast path** (recommended): if `--from-dir` points to a directory containing `current_alignment.fasta`, sequences are sliced from the pre-built MSA and only FastTree is run — seconds to minutes for any sized group.

**Slow path**: if no existing alignment is found, a full MAFFT + FastTree build is performed.

```bash
# Domain-level (auto-detected keyword)
relict subtree --db project.db --taxon archaea      --from-dir preload_out -o archaea_out
relict subtree --db project.db --taxon bacteria     --from-dir preload_out -o bacteria_out

# Phylum — plain name (rank auto-detected) or GTDB-prefixed
relict subtree --db project.db --taxon Bacteroidota       --from-dir preload_out -o bact_out
relict subtree --db project.db --taxon p__Bacillota       --from-dir preload_out -o firm_out

# Family
relict subtree --db project.db --taxon f__Lachnospiraceae --from-dir preload_out -o lachno_out

# Genus
relict subtree --db project.db --taxon g__Ruminococcus    --from-dir preload_out -o rumino_out
```

Accepted taxon formats:

| Input | Detected rank | Notes |
|---|---|---|
| `archaea` / `bacteria` | domain | Auto-matched case-insensitively |
| `Bacillota` | phylum (fallback) | No prefix → defaults to phylum |
| `p__Bacillota` | phylum | GTDB prefix auto-detected |
| `f__Lachnospiraceae` | family | |
| `g__Ruminococcus` | genus | |
| `d__Archaea` | domain | Explicit GTDB domain prefix |

| Parameter | Required | Default | Description |
|---|---|---|---|
| `--db` | ✓ | — | Path to the Relict SQLite database |
| `-o / --out` | ✓ | — | Output directory |
| `--taxon` | ✓ | — | Taxon to extract (see table above) |
| `--rank` | | `auto` | Override rank detection: `domain` / `phylum` / `family` / `genus` / `species` |
| `--from-dir` | | — | Existing preload/run output dir with `current_alignment.fasta` (fast path) |
| `--ref` | | — | Reference FASTA for orientation correction (slow path only) |
| `--anchors` | | bundled | Custom anchor FASTA |
| `--threads` | | `4` | CPU threads |
| `--min-seqs` | | `3` | Minimum matching sequences required to build a tree |
| `--no-tree` | | off | Skip tree building; only write taxonomy TSV and iTOL colour files |
| `--group-phyla SPEC` | | — | Group phyla into a single colour in iTOL (repeatable) |
| `--functional` | | — | TSV file mapping sequence IDs to functional attributes (first column = ID; subsequent columns = attributes). Generates one iTOL file per column (DATASET_BINARY / DATASET_SIMPLEBAR / DATASET_COLORSTRIP). |

Subtree outputs:

| File | Contents |
|---|---|
| `subtree_tree.nwk` | Focused newick tree (anchor-free, nodes labelled `NODE####`) |
| `subtree_alignment.fasta` | Filtered alignment slice used for the tree |
| `subtree_combined_taxonomy.tsv` | ID → taxonomy → confidence for matched sequences |
| `subtree_sequence_list.tsv` | ID, taxonomy, confidence, dataset for matched sequences |
| `itol_phylum_colors.itol` | iTOL colour strip by phylum |
| `itol_phylum_tree_colors.txt` | iTOL branch + clade range colours |
| `itol_dataset_membership.itol` | iTOL strip showing which dataset each sequence came from |

---

### `relict regen-itol`

Regenerate all iTOL colour files from taxonomy already stored in the DB — without re-classifying or rebuilding the tree. Useful when changing phylum groupings.

```bash
relict regen-itol \
  --db project.db \
  --out preload_out \
  --group-phyla archaea \
  --group-phyla "Bacillota,Bacillota_I,Bacillota_A"
```

| Parameter | Required | Default | Description |
|---|---|---|---|
| `--db` | ✓ | — | Path to the Relict SQLite database |
| `-o / --out` | ✓ | — | Output directory (use your preload or run output dir) |
| `--include-datasets` | | all | Comma-separated dataset names to include |
| `--kingdom` | | — | Include only sequences from this domain |
| `--group-phyla SPEC` | | — | Group phyla into a single colour (repeatable) |
| `--functional` | | — | TSV file mapping sequence IDs to functional attributes (first column = ID; subsequent columns = attributes). Generates one iTOL file per column (DATASET_BINARY / DATASET_SIMPLEBAR / DATASET_COLORSTRIP). |

---

## Phylum grouping for iTOL (`--group-phyla`)

Available on all four subcommands. Groups multiple phyla into a single colour entry in the iTOL legend.

| Spec format | Effect |
|---|---|
| `archaea` | All archaeal phyla (detected by `d__Archaea` domain) → one colour labelled **Archaea** |
| `bacteria` | All bacterial phyla → one colour labelled **Bacteria** |
| `"Bacillota,Bacillota_I,Bacillota_A"` | These three phyla → one colour; label = first name (`Bacillota`) |
| `"Firmicutes:Bacillota,Bacillota_I"` | These two phyla → one colour labelled **Firmicutes** |

Multiple `--group-phyla` arguments are independent and additive:

```bash
relict run ... \
  --group-phyla archaea \
  --group-phyla "Bacillota,Bacillota_I,Bacillota_A" \
  --group-phyla "Bacteroidota,Bacteroidota_A"
```

---

## Multiple dataset stacking

Datasets accumulate in the DB. Each `run` sees the growing baseline:

```bash
# Step 1 — baseline
relict preload --fasta hungate.fasta --db project.db --dataset Hungate ...

# Step 2 — first batch: novelty scored against Hungate
relict run --input batch1.fasta --db project.db --dataset Batch1 ...

# Step 3 — second batch: novelty scored against Hungate + Batch1
relict run --input batch2.fasta --db project.db --dataset Batch2 ...
```

---

## Clustering and tree redundancy reduction (`--collapse`)

When `--collapse` is enabled, Relict groups sequences sharing ≥ `--collapse-threshold` identity and the same taxonomy, keeping only one **cluster representative** per group for tree building.

#### Column Explanation

| Column | Meaning |
|---|---|
| ID | Short identifier (e.g. FLZ63) assigned by Relict for tree readability |
| Taxonomy | Full GTDB lineage assigned by the classifier |
| BestHit | Closest reference genome accession (VSEARCH best hit) |
| ClassificationIdentity | % identity to the BestHit reference (VSEARCH alignment) |
| ClassificationConfidence | Confidence of the taxonomy assignment (from the taxa TSV); NA when taxonomy is parsed directly from reference FASTA headers (no confidence column available) |
| NearestHit | Closest sequence among YOUR previously submitted sequences (preload + prior runs) |
| NearestIdentity | % identity to NearestHit |
| MatchesGE99 / GE97 / GE95 | Count of YOUR submitted sequences within 99% / 97% / 95% identity |
| NoveltyScore | 100 − NearestIdentity (higher = more novel vs. your collection) |
| Crowding | Neighbourhood density: crowded / moderate / sparse based on how many of your sequences are nearby |
| SequencingPriority | HIGH / MEDIUM / LOW — suggested priority for follow-up based on novelty + crowding |
| InTree | Yes = entered the phylogenetic tree; No = excluded (see ClusterRepresentative) |
| ClusterRepresentative | `self` = this sequence IS in the tree; an ID = collapsed into that representative; `duplicate` = exact duplicate removed during dereplication |
| ClusterSize | Total sequences in this cluster (1 = singleton) |
| ClusteredMembers | Semicolon-separated IDs of OTHER sequences collapsed under this representative |
| PlacementFlags | Warnings: LOW_CLASSIFICATION_IDENTITY, LOW_NEAREST_IDENTITY, NOVEL_BUT_ASSIGNED, etc. |

### `novelty_metrics.tsv`

---

## All output files

### `relict preload` outputs

| File | Description |
|---|---|
| `current_tree.nwk` | Backbone phylogenetic tree (newick, nodes labelled `NODE####`) |
| `current_alignment.fasta` | Full multiple sequence alignment (includes anchor sequences) |
| `preload_combined_taxonomy.tsv` | ID → taxonomy → confidence for all preloaded sequences |
| `taxonomy.tsv` | Raw classifier output (ID, best-hit, identity, taxonomy, confidence) |
| `itol_phylum_colors.itol` | iTOL colour strip by phylum |
| `itol_phylum_tree_colors.txt` | iTOL branch + clade range colours by phylum |
| `itol_phylum_symbols.itol` | iTOL symbol strip by phylum |
| `itol_family_*.itol / *.txt` | Family-level equivalents |
| `itol_genus_*.itol / *.txt` | Genus-level equivalents |
| `itol_dataset_preload.itol` | iTOL strip marking all sequences with the dataset colour |
| `preload_id_map.tsv` | Short ID → original FASTA header mapping |
| `preload_collapsed_map.tsv` | Cluster representative → taxonomy → count |
| `preload_collapsed_members.tsv` | Member → representative mapping |
| `tree_build_warnings.tsv` | Warnings about sequence quality or alignment issues |
| `tree_orientation_summary.tsv` | Per-sequence orientation audit (forward / RC / unknown) |
| `OUTPUT_EXPLANATIONS.tsv` | Manifest describing each output file |

### `relict run` outputs

| File | Description |
|---|---|
| `sequence_assessment.tsv` | **Main output.** Per-sequence novelty, taxonomy, crowding, priority, clustering, and placement flags |
| `novelty_metrics.tsv` | Per-sequence novelty and crowding summary for candidate ranking |
| `placement_warnings.tsv` | Sequences flagged for low identity or unusual novelty patterns |
| `combined_taxonomy.tsv` | ID → taxonomy → confidence for all sequences in this run |
| `current_tree.nwk` | Updated tree incorporating new sequences |
| `itol_*.itol / *.txt` | All iTOL colour/symbol/range files |
| `itol_dataset_membership.itol` | iTOL strip showing which dataset each sequence belongs to |
| `user_id_map.tsv` | Short ID → original FASTA header mapping for this run |
| `tree_build_warnings.tsv` | Tree-quality warnings |

### `sequence_assessment.tsv` columns

| Column | Meaning |
|---|---|
| `ID` | Short identifier assigned by Relict for tree readability |
| `Taxonomy` | Full GTDB lineage assigned by the classifier |
| `BestHit` | Closest reference genome accession |
| `ClassificationIdentity` | % identity to the BestHit reference |
| `ClassificationConfidence` | Confidence of the taxonomy assignment |
| `NearestHit` | Closest sequence among YOUR previously submitted sequences |
| `NearestIdentity` | % identity to NearestHit |
| `MatchesGE99 / GE97 / GE95` | Count of YOUR submitted sequences within 99% / 97% / 95% identity |
| `NoveltyScore` | 0–100 score; higher = more novel vs. your collection |
| `Crowding` | Neighbourhood density: `crowded` / `moderate` / `sparse` |
| `SequencingPriority` | `HIGH` / `MEDIUM` / `LOW` — suggested follow-up priority |
| `InTree` | `Yes` = in tree; `No` = excluded (see `ClusterRepresentative`) |
| `ClusterRepresentative` | `self`, a representative ID, or `duplicate` |
| `ClusterSize` | Sequences in this cluster |
| `ClusteredMembers` | Other sequences collapsed under this representative |
| `PlacementFlags` | Warning codes (e.g. `LOW_CLASSIFICATION_IDENTITY`) |

---

## Reference inputs — which flag to use

| Scenario | Recommended flags |
|---|---|
| Classify against GTDB reps FASTA | `--ref gtdb.fna --classify` |
| Input FASTA has embedded GTDB lineages in headers | `--taxa-assignments input.fasta` |
| Separate taxonomy TSV for the input sequences | `--taxa-assignments assignments.tsv` |
| GTDB FASTA + separate taxonomy TSV for the reference | `--ref gtdb.fna --taxa gtdb_taxonomy.tsv --classify` |

`--ref` = the reference database used for classification and tree orientation.  
`--taxa-assignments` = pre-computed taxonomy for the *input* sequences themselves.

---

## ID shortening

By default, Relict replaces input headers with compact IDs (e.g. `HUN001`, `QUE06`) so tree labels stay readable. Full original headers are always preserved in `preload_id_map.tsv` or `user_id_map.tsv`.

- `--shorten-ids` (default) — compact generated IDs
- `--no-shorten-ids` — keeps a canonicalized version of the source header

---

## Anchor sequences

Anchor sequences are **invisible scaffolding** constraining phylum-level topology during tree construction. They are pruned from the stored newick after each build and never appear in outputs, iTOL files, or novelty scoring.

The bundled anchor set (`src/relict/data/reference_anchors.fasta`) contains 26 NCBI RefSeq type-strain 16S sequences covering major gut/rumen phyla:

| Phylum | Type strain | Accession |
|---|---|---|
| Bacillota | *Ruminococcus albus* 7 | NR_025930 |
| Bacillota | *Clostridium butyricum* ATCC 19398 | NR_074545 |
| Bacillota | *Bacillus subtilis* 168 | NR_027552 |
| Bacillota | *Lactobacillus acidophilus* ATCC 4356 | NR_075051 |
| Bacillota | *Streptococcus equinus* NBRC 12553 | NR_113594 |
| Bacillota | *Butyrivibrio fibrisolvens* D1 | NR_044858 |
| Bacteroidota | *Bacteroides fragilis* ATCC 25285 | NR_041386 |
| Bacteroidota | *Prevotella bryantii* B14 | NR_044825 |
| Bacteroidota | *Porphyromonas gingivalis* ATCC 33277 | NR_040847 |
| Pseudomonadota | *Escherichia coli* K-12 | NR_102804 |
| Pseudomonadota | *Helicobacter pylori* 26695 | NR_073694 |
| Pseudomonadota | *Wolinella succinogenes* DSM 1740 | NR_043184 |
| Actinomycetota | *Bifidobacterium longum* NCC2705 | NR_040783 |
| Actinomycetota | *Streptomyces griseus* ATCC 23345 | NR_043823 |
| Spirochaetota | *Treponema pallidum* str. Nichols | NR_027243 |
| Spirochaetota | *Borrelia burgdorferi* B31 | NR_025890 |
| Fibrobacterota | *Fibrobacter succinogenes* S85 | NR_041558 |
| Fusobacteriota | *Fusobacterium nucleatum* ATCC 25586 | NR_026043 |
| Planctomycetota | *Planctomyces maris* DSM 8797 | NR_043399 |
| Verrucomicrobiota | *Akkermansia muciniphila* ATCC BAA-835 | NR_042817 |
| Thermotogota | *Thermotoga maritima* MSB8 | NR_043084 |
| Deinococcota | *Deinococcus radiodurans* R1 | NR_036779 |
| Cyanobacteriota | *Synechococcus elongatus* PCC 6301 | NR_043317 |
| Archaea (outgroup) | *Methanobrevibacter ruminantium* M1 | NR_044812 |
| Archaea (outgroup) | *Methanobacterium thermoautotrophicum* Marburg | NR_028243 |
| Archaea (outgroup) | *Sulfolobus acidocaldarius* DSM 639 | NR_043325 |

To use a custom anchor file:
```bash
relict preload --anchors /path/to/my_anchors.fasta ...
```
Custom anchor headers must begin with `RELICT_REF_`:
```
>RELICT_REF_MyPhylum accession=NR_XXXXXX source=SILVA138
```

To refresh the bundled anchors from NCBI:
```bash
python scripts/build_anchor_fasta.py --email your@institution.ac.uk
```

---

## Interpreting outputs

### Sequence priority

| Label | Meaning |
|---|---|
| `HIGH` | Novel, sparse neighbourhood — strong candidate for WGS follow-up |
| `MEDIUM` | Moderately novel or moderate crowding |
| `LOW` | Well-represented in your collection; limited novelty gain |

### Placement flags

| Flag | Meaning |
|---|---|
| `LOW_CLASSIFICATION_IDENTITY` | < 95% identity to best reference hit |
| `LOW_CONFIDENCE` | Classification confidence below threshold |
| `LOW_NEAREST_IDENTITY` | < 95% identity to nearest sequence in your collection |
| `NOVEL_BUT_ASSIGNED` | High novelty score despite a confident taxonomy assignment |

### Tree quality warnings

The tree is a context plot, not absolute proof of novelty. Always pair tree inspection with `sequence_assessment.tsv`. Warnings in `tree_build_warnings.tsv`:

| Warning | Meaning |
|---|---|
| `PARTIAL_16S_SEQUENCES` | Sequences shorter than 1200 bp reduce placement accuracy |
| `HIGH_N_CONTENT` | High ambiguous-base content can distort branch lengths |
| `VERY_SHORT_ALIGNED_FRAGMENTS` | Some aligned rows are mostly gaps after MAFFT |
| `MISSING_REFERENCE_ANCHORS` | Tree was built without topology scaffolding |

If a large `unknown` sector appears in iTOL, check:
- `combined_taxonomy.tsv` or `preload_combined_taxonomy.tsv`
- `taxonomy.tsv` (raw classifier output)
- `taxonomy_input_warnings.tsv`

---

## What "relict" means

Relict is aimed at finding lineages that may have been missed because of primer bias, sparse reference coverage, taxonomy lag, or conservative filtering. Many targets are not new in nature — they are often **new to the reference record** or **underrepresented in existing collections**.

That is why Relict combines reference-aware classification with novelty scoring and warning layers, rather than treating every long branch as a discovery.
