# BranchManager

**Working research prototype. Use with caution.**

**Marker-gene QC, taxonomy, novelty scoring, and isolate prioritisation for genome-sequencing candidate selection.**

BranchManager helps answer the question:

> Has this lineage already been seen and characterised, or is it still poorly represented enough to justify deeper follow-up such as whole-genome sequencing?

It combines Sanger/AB1 quality control, marker-gene taxonomic classification, nearest-neighbour novelty scoring, neighbourhood density (crowding), Most Wanted List matching, and phylogenetic tree visualisation into sequence assessment and selection reports.


---
## Core Ethos
* Absence of a hit is not automatically **novelty**. It may indicate poor sequence quality, a chimera, contamination, inadequate coverage, or reference incompleteness.
* Facts and recommendations remain separate. already_sequenced, genome QC, and DNA availability are factual states; PRIMARY and BACKUP are recommendations.
* Uncertainty must survive into the final report. Do not collapse everything into one opaque score.
* Selection-critical stages fail closed. An incomplete tree, novelty calculation, or candidate set must never produce a superficially successful report.
* Every decision is reproducible from raw AB1 files.
* Project-changing workflows use a staged SQLite copy and publish it only after required outputs and database integrity checks pass.
* 16S supports selection groups, not definitive species or strain boundaries. The 98.65% boundary is useful as a heuristic but is not equivalent to genomic species assignment; published work estimates meaningful exceptions even above that threshold.

## Workflow vocabulary

| Stage | Purpose |
|---|---|
| **Mailroom** | Inventory an AB1 delivery, reconcile supplier read IDs, and prepare the per-batch map |
| **Interview** | Run standalone AB1 conversion, assembly, QC, and resequencing triage before project evaluation |
| **Onboarding** | Validate partner IDs, metadata, and raw-file ownership before analysis |
| **Paper Trail** | Read AB1 base calls, Phred scores, peak positions, dye channels, and mixed-peak evidence |
| **Merge Meeting** | Trim primers and assemble multiple primer reads, or choose the best independent read |
| **Filing Cabinet** | Register cultured baseline isolates and establish the initial phylogenetic context |
| **Performance Review** | Classify markers, screen chimeras, score novelty, and update the MSA/tree |
| **Hiring Panel** | Propose primary and backup isolate sets for the current evidence state |
| **Quarterly Review** | Reconsider the complete collection for a later, budgeted genome tranche |
| **Status Meeting** | Import factual isolate lifecycle progress |
| **Records Update** | Import completed genome/QC/GTDB/ANI evidence |
| **Annual Report** | Produce the cumulative marker-to-genome close-out report |
| **Assistant to the Branch Manager** | Run Onboarding through the Hiring Panel for one submission |

The names are navigation aids. Manifests and reports always retain explicit scientific operation names and thresholds.

## Quick start

For a complete partner submission, use the guarded golden path:

```bash
branchmanager assistant \
  --sample-map ab1_mapping.tsv \
  --partner-metadata partner_metadata.tsv \
  --db project.sqlite --dataset QUB_01 \
  --ref gtdb_ssu_reps.fna.gz --ref-name GTDB \
  --chimera-ref silva_gold_ssu.fasta.gz \
  --baseline-fasta hungate.fasta --baseline-dataset Hungate \
  --mwl MWL.csv --threads 10 \
  -o runs/01_QUB_01
```

`assistant` runs IT Desk checks, Onboarding, Paper Trail/Merge Meeting, Performance Review, and Hiring Panel in numbered directories. Each stage remains independently runnable:

```bash
# 1. Establish the cultured baseline Filing Cabinet and build a reference tree
branchmanager filing-cabinet \
  --fasta hungate.fasta \
  --db project.db \
  --dataset Hungate \
  --ref gtdb_r232_bac_arch_ssu_reps.fna \
  --classify --build-tree \
  -o filing_cabinet_out \
  --threads 10

# 2. Process new sequences — scores novelty against the baseline
branchmanager performance-review \
  --input new_sequences.fasta \
  --db project.db \
  --dataset Batch1 \
  --ref gtdb_r232_bac_arch_ssu_reps.fna \
  --previous-review filing_cabinet_out \
  -o batch1_out \
  --threads 10

# 3. Zoom into a specific taxon
branchmanager org-chart \
  --db project.db \
  --taxon archaea \
  --from-dir filing_cabinet_out \
  -o archaea_out

# 4. Regenerate iTOL files with new grouping options
branchmanager label-maker \
  --db project.db \
  --out filing_cabinet_out \
  --group-phyla archaea \
  --group-phyla "Bacillota,Bacillota_I,Bacillota_A"
```

---

## Core concepts

### Novelty is relative to YOUR submitted sequences

All novelty and neighbourhood-density calculations are made against the sequences submitted to BranchManager (the Filing Cabinet baseline plus all prior Performance Review datasets stored in the DB). They are **not** made against the full external reference (GTDB/SILVA).

This is intentional: BranchManager helps you judge whether a new sequence is worth investigating relative to what your lab or project has already characterised — not relative to all known biology.

Each successive Performance Review extends the reference pool, so novelty scores become increasingly precise as your project grows.

### Three-layer tree architecture

The phylogenetic tree is built from three independent layers:

| Layer | Source | Shown in iTOL? |
|---|---|---|
| 1 — Anchors | 26 NCBI RefSeq anchor sequences (bundled) | **No** — invisible topology scaffolding |
| 2 — Filing Cabinet | Your cultured baseline dataset (e.g. Hungate) | Yes |
| 3 — Performance Review | Each new partner marker batch | Yes |

Anchor sequences constrain phylum-level topology during MAFFT + FastTree inference, then are pruned from the stored newick. They are never shown in outputs, iTOL files, or novelty scoring.

---

## Subcommands

### `branchmanager filing-cabinet` 

Load a baseline FASTA dataset, classify against a reference, and build the backbone tree. This is always the first step.

```bash
branchmanager filing-cabinet \
  --fasta baseline.fasta \
  --db project.db \
  --dataset Hungate \
  --ref gtdb_r232_bac_arch_ssu_reps.fna \
  --classify --build-tree \
  -o filing_cabinet_out \
  --threads 10 \
  --collapse
```

| Parameter | Required | Default | Description |
|---|---|---|---|
| `--fasta` | ✓ | — | Input FASTA file containing the baseline sequences |
| `--db` | ✓ | — | Path to the BranchManager SQLite database (created if absent) |
| `-o / --out` | | `.` | Output directory for tree, iTOL files, and reports |
| `--dataset` | ✓ | — | Label stored in the DB (e.g. `Hungate`) |
| `--ref` | | — | Reference FASTA (GTDB/SILVA reps) for classification and tree orientation |
| `--classify` | | off | Classify sequences against `--ref` and store taxonomy |
| `--build-tree` | | off | Build the MAFFT + FastTree backbone tree |
| `--taxa` | | — | TSV mapping reference IDs to lineages (optional when `--ref` headers contain GTDB lineages) |
| `--taxa-assignments` | | — | Pre-computed taxonomy for the INPUT sequences (TSV or embedded-lineage FASTA) |
| `--shorten-ids / --no-shorten-ids` | | `--no-shorten-ids` | Preserve input headers exactly by default; use `--shorten-ids` only when compact generated IDs are desired |
| `--collapse` | | off | Collapse near-identical same-taxonomy sequences into one representative for the tree |
| `--collapse-threshold` | | `99.8` | Identity threshold (%) for collapsing |
| `--sequence-domain` | | `bacteria` | Domain/profile to process: `bacteria`, `archaea`, `fungi`, or `mixed`/`all`/`none` to disable filtering. Use matching references, previous reviews, and anchors for non-bacterial runs |
| `--anchors` | | bundled | Custom anchor FASTA for tree topology scaffolding |
| `--threads` | | `4` | CPU threads for MAFFT and VSEARCH |
| `--colours` | | — | CSV mapping sequence IDs to custom hex colours for iTOL (columns: `id`, `colour`) |
| `--group-phyla SPEC` | | — | Group phyla into a single colour in iTOL (repeatable; see *Phylum grouping*) |
| `--functional` | | — | TSV file mapping sequence IDs to functional attributes (first column = ID; subsequent columns = attributes). Generates one iTOL file per column (DATASET_BINARY / DATASET_SIMPLEBAR / DATASET_COLORSTRIP). |

---

### `branchmanager mailroom`

Prepare the immutable per-batch AB1 map before Onboarding. Mailroom scans every physical chromatogram, extracts ABIF sample/run/instrument and base-call evidence, reconciles supplier sequencing IDs, and reports missing, duplicated, ambiguous, or unmapped reads.

```bash
branchmanager mailroom \
  --read-dir UoG/All_AB1 \
  --metadata UoG/supplier_metadata.csv \
  --dataset UoG_01 \
  --forward-primer 63F \
  --reverse-primer 1492R \
  -o UoG/UoG_01
```

Only use `--forward-primer` and `--reverse-primer` after the oligonucleotides have been confirmed from the partner protocol or sequencing submission. Supplier metadata may instead contain a `primer` column. Mailroom uses embedded ABIF primer names when present, but never treats sequence-motif similarity as proof of the exact primer.

Preferred supplier-metadata columns:

```tsv
sequencing_id	isolate_number	read	primer	processing_mode
KKX994	SW_0016	Forward	63F	assemble
KKY011	SW_0016	Reverse	1492R	assemble
```

`primer` and `processing_mode` are optional. With `--processing-mode auto`, paired/multi-read isolates use `assemble` and single-read isolates use `best_read`.

| Output | Purpose |
|---|---|
| `ab1_map.tsv` | Validated one-row-per-read input for Onboarding |
| `ab1_inventory.tsv` | Physical file and extracted ABIF metadata inventory |
| `mailroom_report.tsv` | Missing, duplicate, ambiguous, unreadable, unmapped, and unresolved-primer findings |
| `mailroom_summary.json` | Machine-readable counts and `PASS`, `REVIEW_REQUIRED`, or `FAIL` status |
| `run_manifest.json` | Checksums, inputs, outputs, and stage status |

Mailroom returns `REVIEW_REQUIRED` when files reconcile but primer names remain unresolved, and `FAIL` for structural/file errors. Correct the supplier metadata or confirm primers before Onboarding.

---

### `branchmanager interview`

Run AB1 conversion and QC without Onboarding, taxonomy, tree building, candidate selection, or any project-database access. Interview consumes the complete Mailroom output directory and uses the same versioned QC policy as `assistant`.

```bash
branchmanager interview \
  --mailroom UoG/UoG_01 \
  --screen-ref gtdb_r232_bac_arch_ssu_reps.fna.gz \
  --threads 10 \
  -o UoG/UoG_01_interview
```

`--screen-ref` is optional and only checks whether separate primer reads from the same isolate have discordant taxonomic hits. Without it, all chromatogram, trimming, assembly, quality, ambiguity, mixed-peak, and resequencing outputs are still produced.

Interview accepts either the Mailroom directory or its `ab1_map.tsv`. It requires the accompanying `mailroom_summary.json`, rejects a Mailroom `FAIL`, and stops on `REVIEW_REQUIRED` unless the batch has been reviewed and `--allow-mailroom-review` is supplied explicitly.

The principal hand-off files are `assembled.fasta`, `assembly_report.tsv`, `marker_review_template.tsv`, and `resequence_recommendations.tsv`. Once any manual-review decisions are completed, these can enter the full workflow without repeating chromatogram processing:

```bash
branchmanager assistant \
  --fasta UoG/UoG_01_interview/assembled.fasta \
  --marker-qc UoG/UoG_01_interview/assembly_report.tsv \
  --marker-review UoG/UoG_01_interview/marker_review_template.tsv \
  --partner-metadata Sequence_Metadata.csv \
  --db project.sqlite --dataset UoG_01 --ref gtdb_ssu_reps.fna.gz \
  -o runs/UoG_01
```

Do not supply the blank review template to Assistant. Add decisions only for the listed `PASS_WITH_WARNINGS` markers first; omit `--marker-review` when the template contains no rows.

---

### `branchmanager onboarding`

Every partner dataset is onboarded from exactly one marker source. AB1 submissions use a sample map; FASTA-only submissions use the FASTA directly. Both paths validate sequence IDs against the same cumulative project metadata ledger.

```bash
# AB1 or multiple primer reads
branchmanager onboarding \
  --sample-map UoG/UoG_01/ab1_map.tsv \
  --read-dir UoG_01_AB1/ \
  --partner-metadata project_partner_metadata.tsv \
  --partner-id UoG --dataset UoG_01 \
  -o UoG_01/01_onboarding

# Partner-supplied marker FASTA; .gz is accepted
branchmanager onboarding \
  --fasta QUB_01_16S.fasta.gz \
  --partner-metadata project_partner_metadata.tsv \
  --partner-id QUB --dataset QUB_01 \
  -o QUB_01/01_onboarding
```

FASTA Onboarding writes `normalised_input.fasta`. It does not invent chromatogram quality evidence. Supply a reviewed `--marker-qc` table to Performance Review where one exists, or use `--accept-unverified-marker-qc` as an explicit audited acceptance; unverified marker evidence remains visible in the reports.

Partner acronyms and submission labels are separate. For example, `UoG_01` and `UoG_02` are distinct `--dataset` values, while both sets of ledger rows use `partner_id=UoG`. Sequence IDs must be unique and stable across the entire project.

Keep exactly one cumulative metadata ledger with one row per durable isolate/marker ID. It owns `partner_id`, `selected_for_genome_sequencing`, and `already_sequenced`; chromatogram run IDs do not belong in it. Keep a separate read map under each partner/batch, for example `UoG/UoG_01/ab1_mapping.csv`, which maps every physical trace filename to its durable isolate ID, direction, and `processing_mode`. Duplicate ledger IDs, partner mismatches, missing traces, and traces assigned to more than one isolate are Onboarding errors.

Onboarding reports only actionable errors and warnings. Normalised submission/read-map tables omit optional columns that are empty throughout the submitted batch.

---

### `branchmanager paper-trail` (AB1 conversion and Merge Meeting)

Convert Sanger AB1 chromatograms or already-basecalled primer reads into a trimmed FASTA for the Performance Review. When multiple reads belong to the same isolate, for example `27F` and `907R`, BranchManager orients reverse-primer reads and builds a quality-aware overlap consensus.

```bash
branchmanager paper-trail \
  --input ab1_reads/ \
  -o paper_trail_out \
  --min-quality 20 \
  --min-mean-quality 20 \
  --min-length 800 \
  --min-overlap 40

branchmanager performance-review \
  --input paper_trail_out/assembled.fasta \
  --partner-metadata new_sequences_metadata.tsv \
  --db project.db \
  --dataset Batch1 \
  --ref gtdb_r232_bac_arch_ssu_reps.fna \
  -o batch1_out
```

If filenames are structured like `Iso001_27F.ab1` and `Iso001_907R.ab1`, BranchManager can infer the read information. Production partner submissions should nevertheless use one explicit `--sample-map` per batch, preferably with one row per physical read:

```tsv
sequence_id	dataset	read_file	primer	direction	processing_mode
Iso001	UoG_01	well_A01.ab1	27F	forward	assemble
Iso001	UoG_01	well_A02.ab1	907R	reverse	assemble
Iso002	UoG_01	well_B01.ab1	27F	forward	best_read
```

Multiple rows may share the same `sequence_id` when an isolate has multiple primer reads. Relative paths are resolved from the batch-map location. A wide one-row-per-isolate representation is also accepted when laboratories naturally supply separate primer columns:

```tsv
isolate_id	27F	907R	processing_mode
Iso001	well_A01.ab1	well_A02.ab1	assemble
Iso002	well_B01.ab1	well_B02.ab1	best_read
```

Run it with:

```bash
branchmanager paper-trail --sample-map sample_reads.tsv -o paper_trail_out
```

Onboarding writes `normalised_read_map.tsv` as the validated hand-off to Paper Trail. It is a generated run artefact, not another source ledger and not a file users should maintain.

`assemble`/`merge`/`consensus` orients the primer reads and tries to build one longer sequence. `best_read`/`highest_quality`/`select_best`/`independent` converts each read independently and writes only the highest-quality passing read for that isolate.

BranchManager applies the same versioned Sanger QC policy in `interview`, `paper-trail`, and `assistant`:

- `--min-quality` is a Phred cutoff used for Mott-style end trimming.
- `--mask-quality` masks internal bases below Q20 to `N` before assembly; masking is reported but does not by itself force manual review.
- Primer sequences are removed only after a confident IUPAC-aware leading match; customise them with repeatable `--primer-sequence NAME=SEQUENCE`.
- Mixed chromatogram peaks are measured from `PLOC` and four `DATA9`-`DATA12` dye channels. `--secondary-peak-ratio` and `--max-mixed-peak-percent` control the retained-region gate.
- Overlap consensus bases use posterior base probabilities from both Phred observations; qualities are not added as if they were independent certainty scores.
- `--screen-ref` optionally classifies each primer read independently and rejects family/genus discordance before consensus use.
- Final sequences must be at least 800 bp. Individual primer reads may be at least 300 bp so complementary reads can assemble into a valid final marker.
- Ambiguous bases use a tiered rule: more than 3% requires manual review and more than 5% fails QC by default.
- `--min-mean-quality`, expected-error limits, ambiguity, internal low-quality runs, mixed peaks, and overlap conflict density determine `PASS_HIGH_CONFIDENCE`, `PASS_WITH_WARNINGS`, or `FAIL_QC`.
- Final `FAIL_QC` sequences are withheld from `assembled.fasta` and listed as `RESEQUENCE`.
- Final `FAIL_QC` sequences and failed reads are retained under `failed_qc_sequences/` for manual review.
- `PASS_WITH_WARNINGS` sequences are included but listed as `MANUAL_REVIEW`.
- Visual reports are automatically split into numbered PNG pages before they exceed 2400 pixels high. Adjust the ceiling with `--max-report-image-height` (minimum 600 pixels).

The defaults follow Phred error probabilities and published Sanger guidance that Q20 is a reliable base call and that, as a general rule, trimmed sequence should contain less than 5% ambiguous bases. The 3% review boundary is deliberately more conservative for phylogenetic candidate selection. See [Ewing and Green 1998](https://doi.org/10.1101/gr.8.3.186) and [Crossley et al. 2020](https://doi.org/10.1177/1040638720905833).

Outputs:

| File | Description |
|---|---|
| `assembled.fasta` | One trimmed/assembled sequence per isolate, suitable for `branchmanager performance-review --input` |
| `failed_qc_sequences/failed_final_sequences.fasta` | Final isolate consensus/best-read sequences that failed final QC |
| `failed_qc_sequences/failed_read_sequences.fasta` | Individual read sequences that failed read-level QC, oriented when direction is known |
| `failed_qc_sequences/failed_qc_manifest.tsv` | One row per failed isolate, combining the final outcome with read IDs, source files, read-failure evidence, and the best recoverable read |
| `failed_qc_sequences/failed_read_manifest.tsv` | One row per failed physical read with trimming and QC evidence; kept separate to avoid duplicate isolate rows |
| `failed_qc_sequences/failed_qc_guide.txt` | Plain-language explanation of failure reason codes and the distinction between isolate- and read-level records |
| `trimmed_oriented_reads.fasta` | Individual reads after quality trimming and primer-direction orientation |
| `raw_reads.fasta` | Raw base calls extracted from AB1 or input sequence files |
| `read_qc.tsv` | Per-read trimming, quality, expected error, length, and filter status |
| `per_base_error.tsv` | Per-base quality/error probability table with left/right trim and retained-base status |
| `visual_reports/read_error_profiles/read_error_profiles_page_*.png` | Paginated per-read quality/error profiles showing retained trim windows |
| `visual_reports/trace_chromatograms/trace_chromatograms_page_*.png` | Paginated dye-channel chromatograms with retained windows and mixed-peak evidence |
| `read_taxonomy_concordance.tsv` | Independent primer-read classifications and within-isolate agreement when `--screen-ref` is used |
| `assembly_report.tsv` | Per-isolate assembly status, overlap identity, conflicts, unmerged reads, and contributing read IDs |
| `assembly_read_placements.tsv` | Per-read consensus coordinates, aligned blocks, and whether the read contributed to the selected consensus |
| `resequence_recommendations.tsv` | Per-isolate `ACCEPT`, `MANUAL_REVIEW`, or `RESEQUENCE` decision with reason codes and suggested action |
| `marker_review_template.tsv` | Pre-filled rows for manual-review markers; complete `decision`, `reviewer`, and `notes`, then provide it as `--marker-review` |
| `paper_trail_qc_policy.tsv` | Versioned thresholds actually used for the run, for reproducibility |
| `visual_reports/assembly_overviews/assembly_overview_page_*.png` | Paginated per-isolate overviews with reads placed on consensus coordinates, contribution status, and assembly diagnostics |
| `visual_reports/visual_report_manifest.tsv` | Page index with report type, record range, dimensions, and continuation notes |
| `paper_trail_summary.txt` | Short run summary |
| `run_manifest.json` / `.tsv` | Input checksums, software/tool versions, thresholds, stages, outputs, and final status |

Supported inputs are `.ab1`, `.abi`, `.fasta`, `.fa`, `.fna`, `.fastq`, and `.fq`, with `.gz` accepted for AB1/ABI, FASTA, and FASTQ files.

---

### `branchmanager performance-review`

Process new sequences against the baseline; score novelty and update the phylogenetic tree.

```bash
branchmanager performance-review \
  --input new_sequences.fasta \
  --db project.db \
  --dataset Batch1 \
  --ref gtdb_r232_bac_arch_ssu_reps.fna \
  --previous-review filing_cabinet_out \
  -o batch1_out \
  --threads 10 \
  --collapse
```

| Parameter | Required | Default | Description |
|---|---|---|---|
| `--input` | ✓ | — | FASTA file of new sequences to analyse |
| `--db` | ✓ | — | Path to the BranchManager SQLite database |
| `-o / --out` | ✓ | — | Output directory (assessment TSV, novelty metrics, tree, iTOL files) |
| `--dataset` | ✓ | — | Label for this batch (used in iTOL dataset-membership strip) |
| `--ref` | | — | Reference FASTA for classification and tree orientation |
| `--previous-review` | | — | Previous Filing Cabinet or Performance Review output directory used to seed an incremental MSA. When the same `-o` directory is reused, `tree/current_alignment.fasta` is detected automatically |
| `--taxa` | | — | TSV mapping reference IDs to lineages |
| `--taxa-assignments` | | — | Pre-computed taxonomy for the input sequences |
| `--partner-metadata / --sequencing-metadata` | Performance Review | — | Cumulative CSV/TSV ledger with sequence IDs, partner acronyms, optional `selected_for_genome_sequencing`, and required `already_sequenced` status |
| `--shorten-ids / --no-shorten-ids` | | `--no-shorten-ids` | Preserve input headers exactly by default; use `--shorten-ids` only when compact generated IDs are desired |
| `--min-len` | | `800` | Minimum sequence length to retain (bp) |
| `--max-n` | | `5` | Maximum ambiguous (N) bases allowed |
| `--marker-qc` | | auto | Paper Trail/Merge Meeting `assembly_report.tsv`; auto-discovered beside the FASTA |
| `--marker-review` | | — | Reviewed accept/reject decisions for warning or unverified marker evidence |
| `--accept-unverified-marker-qc` | | off | Explicit, audited acceptance of independently validated FASTA without Paper Trail provenance |
| `--chimera-ref` | | primary ref | Curated chimera-free marker reference used by reference UCHIME |
| `--skip-chimera-check` | | off | Explicitly omit UCHIME; affected isolates remain review-required |
| `--collapse` | | off | Collapse near-identical same-taxonomy sequences for the tree |
| `--collapse-threshold` | | `99.8` | Identity threshold (%) for collapsing |
| `--sequence-domain` | | `bacteria` | Domain/profile to process: `bacteria`, `archaea`, `fungi`, or `mixed`/`all`/`none` to disable filtering. Run archaea/fungi separately with matching `--ref`, `--alt-ref`, `--baseline-fasta`, `--previous-review`, and `--anchors` |
| `--phylum` | | — | Filter iTOL output to a specific phylum (does not affect novelty scoring) |
| `--target` | | — | FASTA to measure novelty against instead of the DB |
| `--baseline-fasta` | | — | Baseline/context FASTA to load before evaluating, e.g. Hungate |
| `--baseline-dataset` | | `Baseline` | Dataset label for `--baseline-fasta` and default cultured-baseline novelty pool |
| `--novelty-baseline-dataset` | | — | Existing DB dataset to include in the baseline novelty pool; repeatable |
| `--baseline-shorten-ids / --no-baseline-shorten-ids` | | `--no-baseline-shorten-ids` | Preserve baseline IDs exactly by default |
| `--force-rebuild` | | off | Force a full tree rebuild from scratch |
| `--anchors` | | bundled | Custom anchor FASTA for tree scaffolding |
| `--threads` | | `4` | CPU threads |
| `--neighbourhood-format` | | `png` | Local-clade figure format; PNG is the only supported output |
| `--user-colours` | | — | CSV mapping sequence IDs to custom hex colours for iTOL |
| `--group-phyla SPEC` | | — | Group phyla into a single colour in iTOL (repeatable) |
| `--functional` | | — | TSV file mapping sequence IDs to functional attributes (first column = ID; subsequent columns = attributes). Generates one iTOL file per column (DATASET_BINARY / DATASET_SIMPLEBAR / DATASET_COLORSTRIP). |

### Rolling project runs

Use one persistent SQLite database for the project and give each partner submission a unique `--dataset` label. The input FASTA may contain only the genuinely new isolates; `--partner-metadata` may be a cumulative table containing both new and previously stored partner IDs, so changed genome-selection flags are applied before the collection is reassessed.

For an auditable history, use a new output directory for every Performance Review and pass the preceding review with `--previous-review`. Reusing the same output directory is also supported and automatically resumes from `tree/current_alignment.fasta`, but replaces that directory's previous reports. Keep reference-database versions, domain, QC settings, collapse settings, anchors, and tree method constant when comparing reviews. Do not use `--force-rebuild` unless a deliberate full reanalysis is required.

Sequence IDs are persistent keys. Resubmitting an existing ID with the same sequence is allowed; resubmitting it with a different sequence is rejected rather than silently changing the project record. External-reference and cultured-baseline evidence should remain stable when inputs and settings are unchanged. Project crowding, nearest-project hit, genome-collection coverage, pangenome gap, selection-set role, and recommendation are deliberately recalculated because new isolates and updated WGS commitments change those questions.

Incremental runs retain the old MSA columns and add new sequences with MAFFT. The tree builder still estimates a new tree from the expanded alignment, so small changes in old topology or branch lengths remain scientifically possible. Use the preceding Performance Review directory as `--previous-review`, keep parameters/reference versions fixed, inspect `assessment/decision_changes.tsv`, and reserve `--force-rebuild` for deliberate recalibration. Strictly fixed old placements require a validated fixed-backbone phylogenetic-placement workflow and are not claimed by BranchManager.

Domain/profile handling:

BranchManager defaults to `--sequence-domain bacteria` for Filing Cabinet and Performance Review. This keeps bacterial 16S runs clean by filtering non-bacterial assignments before DB insertion and tree building where taxonomy is available. Process archaea and fungi as separate reviews with their own DB, outputs, baselines, references, and anchors:

```bash
branchmanager performance-review ... --sequence-domain archaea --ref archaeal_16s_refs.fasta --previous-review archaea_filing_cabinet -o archaea_review
branchmanager performance-review ... --sequence-domain fungi --ref fungal_its_or_18s_refs.fasta --anchors fungal_anchors.fasta -o fungi_review
```

Use `--sequence-domain mixed` when you intentionally want to keep all domains in one run.

---

### `branchmanager org-chart`

Extract all sequences matching a given taxon from the DB and build a focused phylogenetic tree for that group only.

**Fast path** (recommended): if `--from-dir` points to a directory containing `current_alignment.fasta`, sequences are sliced from the pre-built MSA and only FastTree is run — seconds to minutes for any sized group.

**Slow path**: if no existing alignment is found, a full MAFFT + FastTree build is performed.

```bash
# Domain-level (auto-detected keyword)
branchmanager org-chart --db project.db --taxon archaea      --from-dir filing_cabinet_out -o archaea_out
branchmanager org-chart --db project.db --taxon bacteria     --from-dir filing_cabinet_out -o bacteria_out

# Phylum — plain name (rank auto-detected) or GTDB-prefixed
branchmanager org-chart --db project.db --taxon Bacteroidota       --from-dir filing_cabinet_out -o bact_out
branchmanager org-chart --db project.db --taxon p__Bacillota       --from-dir filing_cabinet_out -o firm_out

# Family
branchmanager org-chart --db project.db --taxon f__Lachnospiraceae --from-dir filing_cabinet_out -o lachno_out

# Genus
branchmanager org-chart --db project.db --taxon g__Ruminococcus    --from-dir filing_cabinet_out -o rumino_out
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
| `--db` | ✓ | — | Path to the BranchManager SQLite database |
| `-o / --out` | ✓ | — | Output directory |
| `--taxon` | ✓ | — | Taxon to extract (see table above) |
| `--rank` | | `auto` | Override rank detection: `domain` / `phylum` / `family` / `genus` / `species` |
| `--from-dir` | | — | Existing Filing Cabinet or Performance Review output directory with `current_alignment.fasta` or `tree/current_alignment.fasta` (fast path) |
| `--ref` | | — | Reference FASTA for orientation correction (slow path only) |
| `--anchors` | | bundled | Custom anchor FASTA |
| `--threads` | | `4` | CPU threads |
| `--min-seqs` | | `3` | Minimum matching sequences required to build a tree |
| `--no-tree` | | off | Skip tree building; only write taxonomy TSV and iTOL colour files |
| `--group-phyla SPEC` | | — | Group phyla into a single colour in iTOL (repeatable) |
| `--functional` | | — | TSV file mapping sequence IDs to functional attributes (first column = ID; subsequent columns = attributes). Generates one iTOL file per column (DATASET_BINARY / DATASET_SIMPLEBAR / DATASET_COLORSTRIP). |

Org Chart outputs:

| File | Contents |
|---|---|
| `org_chart_tree.nwk` | Focused newick tree (anchor-free, nodes labelled `NODE####`) |
| `org_chart_alignment.fasta` | Filtered alignment slice used for the tree |
| `org_chart_combined_taxonomy.tsv` | ID → taxonomy → confidence for matched sequences |
| `org_chart_sequence_list.tsv` | ID, taxonomy, confidence, dataset for matched sequences |
| `itol_phylum_colours.itol` | iTOL colour strip by phylum |
| `itol_dataset_membership.itol` | iTOL strip showing which dataset each sequence came from |

---

### `branchmanager label-maker`

Regenerate all iTOL colour files from taxonomy already stored in the DB — without re-classifying or rebuilding the tree. Useful when changing phylum groupings.

```bash
branchmanager label-maker \
  --db project.db \
  --out filing_cabinet_out \
  --group-phyla archaea \
  --group-phyla "Bacillota,Bacillota_I,Bacillota_A"
```

| Parameter | Required | Default | Description |
|---|---|---|---|
| `--db` | ✓ | — | Path to the BranchManager SQLite database |
| `-o / --out` | ✓ | — | Output directory (use your Filing Cabinet or Performance Review output directory) |
| `--include-datasets` | | all | Comma-separated dataset names to include |
| `--sequence-domain` | | — | Include only sequences from this domain/profile |
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
branchmanager performance-review ... \
  --group-phyla archaea \
  --group-phyla "Bacillota,Bacillota_I,Bacillota_A" \
  --group-phyla "Bacteroidota,Bacteroidota_A"
```

---

## Multiple dataset stacking

Datasets accumulate in the DB. Each Performance Review sees the growing baseline and project collection:

```bash
# Step 1 — baseline
branchmanager filing-cabinet --fasta hungate.fasta --db project.db --dataset Hungate ...

# Step 2 — first batch: novelty scored against Hungate
branchmanager performance-review --input batch1.fasta --db project.db --dataset Batch1 ...

# Step 3 — second batch: novelty scored against Hungate + Batch1
branchmanager performance-review --input batch2.fasta --db project.db --dataset Batch2 ...
```

---

## Clustering and tree redundancy reduction (`--collapse`)

When `--collapse` is enabled, BranchManager groups sequences sharing ≥ `--collapse-threshold` identity and the same taxonomy, keeping only one **cluster representative** per group for tree building.

#### Column Explanation

| Column | Meaning |
|---|---|
| ID | User-supplied sequence ID (preserved unless ID shortening is explicitly requested) |
| Taxonomy | Full GTDB lineage assigned by the classifier |
| BestHit | Closest reference genome accession (VSEARCH best hit) |
| ClassificationIdentity | % identity to the BestHit reference (VSEARCH alignment) |
| ClassificationConfidence | Confidence of the taxonomy assignment (from the taxa TSV); NA when taxonomy is parsed directly from reference FASTA headers (no confidence column available) |
| NearestHit | Closest sequence among YOUR previously submitted sequences (Filing Cabinet + prior Performance Reviews) |
| NearestIdentity | % identity to NearestHit |
| MatchesGE99 / GE97 / GE95 | Count of YOUR submitted sequences within 99% / 97% / 95% identity |
| NoveltyScore | 100 − NearestIdentity (higher = more novel vs. your collection) |
| Crowding | Neighbourhood density: crowded / moderate / sparse based on how many of your sequences are nearby |
| SequencingPriority | HIGH / MEDIUM / LOW — suggested priority for follow-up based on novelty + crowding |
| Reference* | Parallel nearest hit, identity, score, crowding, and priority against the selected external reference FASTA, usually GTDB |
| InTree | Yes = entered the phylogenetic tree; No = excluded (see ClusterRepresentative) |
| ClusterRepresentative | `self` = this sequence IS in the tree; an ID = collapsed into that representative; `duplicate` = exact duplicate removed during dereplication |
| ClusterSize | Total sequences in this cluster (1 = singleton) |
| ClusteredMembers | Semicolon-separated IDs of OTHER sequences collapsed under this representative |
| PlacementFlags | Warnings: LOW_CLASSIFICATION_IDENTITY, LOW_NEAREST_IDENTITY, NOVEL_BUT_ASSIGNED, etc. |

### `novelty_metrics.tsv`

This table keeps the novelty lenses separate:

- `Baseline*`: nearest hit and density against explicit cultured/baseline datasets such as Hungate.
- `Project*`: nearest hit and density against all partner candidates accumulated across runs, including the current collection but excluding each query's self-hit.
- `Reference*`: nearest hit and density against the selected external reference FASTA supplied with `--ref`, usually GTDB.
- `GenomeCollection*` / `Pangenome*`: rolling genome coverage. Every baseline isolate counts because its genome is available; partner isolates count once selected. The default target is nine committed genomes per exact GTDB species or unresolved local clade.
- `SelectionGroupType`: `BASELINE_PANGENOME_EXTENSION` (`BMEXT_*`) means at least one baseline genome anchors the exact GTDB species; `CANDIDATE_PANGENOME_GROUP` (`BMSET_*`) has no baseline-genome anchor.
- `BaselineExtension*`: membership gate for a baseline-pangenome extension. The candidate and nearest baseline must have the same exact GTDB species, with >=98.65% 16S identity across >=95% of the query by default. A failed gate is `PANGENOME_BOUNDARY_REVIEW`, not an automatic rejection: it may represent a separate candidate lineage.
- `BaselineRedundancy*`: hard eligibility gate for uncommitted candidates. By default, a nearest cultured-baseline hit at >=99.8% identity across >=95% of the query is reported as `BASELINE_REDUNDANT` and receives no panel rank.
- `SelectionDiversityDistance`: marginal patristic distance from baseline and already committed genome markers when tree distances are available. This drives farthest-first diversity capture; cultured-baseline, project, and GTDB-reference divergence provide the fallback, with MWL evidence used as an additional tie-break signal.

A species represented by one baseline genome therefore starts at `1/9`, leaving eight primary genome slots. Eligible partner isolates are ranked to extend that baseline pangenome; additional ranked isolates can be retained as extraction-failure/diversity backups. MWL evidence can order eligible isolates but cannot override the species/coverage boundary or the near-identical baseline exclusion.

---

## All output files

### Filing Cabinet outputs

| File | Description |
|---|---|
| `current_tree.nwk` | Backbone phylogenetic tree (newick, nodes labelled `NODE####`) |
| `current_alignment.fasta` | Full multiple sequence alignment (includes anchor sequences) |
| `filing_cabinet_combined_taxonomy.tsv` | ID → taxonomy → confidence for all registered baseline sequences |
| `taxonomy.tsv` | Raw classifier output (ID, best-hit, identity, taxonomy, confidence) |
| `itol_phylum_colours.itol` | iTOL colour strip by phylum |
| `itol_family_colours.itol` | iTOL colour strip by family |
| `itol_genus_colours.itol` | iTOL colour strip by genus |
| `filing_cabinet_dataset.itol` | iTOL strip marking all sequences with the dataset colour |
| `filing_cabinet_id_map.tsv` | Short ID → original FASTA header mapping |
| `filing_cabinet_collapsed_map.tsv` | Cluster representative → taxonomy → count |
| `filing_cabinet_collapsed_members.tsv` | Member → representative mapping |
| `tree_build_warnings.tsv` | Warnings about sequence quality or alignment issues |
| `tree_orientation_summary.tsv` | Per-sequence orientation audit (forward / RC / unknown) |
| `OUTPUT_EXPLANATIONS.tsv` | Manifest describing each output file |

### Performance Review and Hiring Panel outputs

| File | Description |
|---|---|
| `performance_review_dashboard.html` | Compact linked view of the current Hiring Panel recommendations |
| `assessment/sequence_assessment.tsv` | **Full audit table.** Per-sequence novelty, taxonomy, crowding, priority, clustering, MWL matches, and placement flags |
| `assessment/selection_summary.tsv` | **SAB decision table.** Recommendation, evidence quality, component novelty evidence, same-species genome coverage, candidate-set role, rationale, and local-tree link |
| `assessment/sequencing_sets.tsv` | **Rolling nine-member diversity plan.** One row per candidate in a `BMEXT_*` baseline-pangenome extension or `BMSET_*` candidate-only group, with rank, marginal diversity, MWL context, panel completeness, and `PRIMARY`, `BACKUP`, `DIVERSITY_CANDIDATE`, `ALTERNATE`, `PANGENOME_BOUNDARY_REVIEW`, `BASELINE_REDUNDANT`, `SEQUENCED`, or `REVIEW_EVIDENCE` role |
| `assessment/novelty_metrics.tsv` | Per-sequence novelty and crowding summary for candidate ranking |
| `assessment/neighbourhoods/clade_*.png` | Labelled local phylogenetic neighbourhoods (default). Nearby assessed isolates are grouped into one figure rather than duplicated across figures |
| `assessment/neighbourhoods/clade_*_pairwise_pident.tsv` | Complete long-form MSA percent-identity table for every pair of displayed tree leaves, including compared-column counts and the identity definition |
| `assessment/neighbourhoods/neighbourhood_manifest.tsv` | Maps every assessed sequence to its image and pident table, including the P1 identity anchor, forced nearest-baseline hits, assessed isolates, baseline leaves, and already-sequenced genomes shown |
| `assessment/mwl_matches.tsv` | Most Wanted List hits when `--mwl` is supplied |
| `assessment/cluster_summary.tsv` | Cluster-level prioritisation report when clustering/tree reports are available |
| `assessment/clusters.csv` | One consolidated membership table containing all clusters and all member isolates; replaces the former one-file-per-cluster directory |
| `assessment/backup_candidates.tsv` | Ranked alternative isolates within clusters when a primary candidate cannot be sequenced |
| `baseline/baseline_hits.tsv` | Nearest-hit report against Hungate or other provided baseline datasets |
| `baseline/nearest_baseline_hits_raw.tsv` | Raw nearest-hit output used for baseline novelty scoring |
| `taxonomy/<DB>.tsv` | Taxonomic assignment report for each reference database, e.g. `taxonomy/GTDB.tsv`, `taxonomy/GG2.tsv` |
| `taxonomy/all_databases.tsv` | Multi-database assignment summary |
| `taxonomy/tree_taxonomy.tsv` | ID → taxonomy → confidence table used for tree/iTOL metadata |
| `tree/current_tree.nwk` | Updated tree incorporating new sequences |
| `tree/current_alignment.fasta` | MSA used to build the tree |
| `tree/tree_build_warnings.tsv` | Tree-quality warnings |
| `itol/phylum.itol` | iTOL colour strip by phylum |
| `itol/family.itol` | iTOL colour strip by family |
| `itol/genus.itol` | iTOL colour strip by genus |
| `itol/dataset_membership.itol` | iTOL strip showing which dataset each sequence belongs to |
| `itol/novelty.itol` | iTOL strip showing nearest-hit novelty |
| `ids/user_id_map.tsv` | Short ID → original FASTA header mapping for this run |
| `intermediate/` | QC, dereplication, collapse, and classifier scratch outputs retained for debugging |
| `logs/branchmanager.log` | Run log |

BranchManager keeps iTOL `DATASET_COLORSTRIP` files, one per metadata type. Older `TREE_COLORS` branch/range files and symbol-strip variants are removed because they encoded the same metadata in additional visual styles.

### `selection_summary.tsv`

This is the board-facing table for scientific advisory board discussions. It avoids an opaque adjusted score and keeps the evidence components visible: cultured-baseline novelty, rolling project coverage, GTDB-reference divergence, rank-aware MWL evidence, exact same-species genome counts, marker evidence quality, candidate-set role, and local-tree figure.

`SelectedForGenomeSequencing` and `GenomeAlreadySequenced` keep pending commitments distinct from completed genomes. `SelectedPendingGenomesSameSpecies` contributes to `CommittedGenomesSameSpecies` but not to nearest-available-genome identity.

The decision labels are:

| Decision | Meaning |
|---|---|
| `PRIORITISE - SET PRIMARY` | Fills one of the missing genomes needed to reach the per-species pangenome target |
| `RESERVE - SET BACKUP` | Additional phylogenetically spread isolate retained for extraction failure and strain-diversity capture |
| `STRONG CANDIDATE` / `SECONDARY CANDIDATE` | Useful evidence, but not currently assigned a primary/backup place in the working set |
| `SECONDARY - STRAIN DIVERSITY` | The numerical target is met, but the isolate adds a ranked diversity direction not represented by baseline/committed markers |
| `REVIEW - PANGENOME BOUNDARY` | The species, identity, or coverage evidence does not justify placing the isolate in the baseline-pangenome extension; review it as a possible separate lineage |
| `EXCLUDE - BASELINE REDUNDANT` | >=99.8% identity to a cultured baseline marker across >=95% query coverage; retained in reports but excluded from ranking |
| `LOWER PRIORITY - TARGET MET` | Available genomes plus committed pending selections meet the requested exact-GTDB-species target |
| `REVIEW BEFORE SELECTION` | Marker classification, QC, disagreement, or required comparison evidence needs review |
| `ALREADY SELECTED - GENOME PENDING` | The isolate is already committed to sequencing but no usable genome is available yet |
| `ALREADY SEQUENCED` | The current isolate already has an available genome and is not a new sequencing recommendation |

### Project-wide expansion rounds

The default Performance Review targets nine committed genomes per exact GTDB species, or per local phylogenetic clade when species assignment is unresolved. It writes up to nine ranked eligible partner candidates per group. Baseline and completed/pending partner genomes determine the remaining target gap, while candidates outside that gap remain ranked as backups or post-target diversity options.

Near-identical cultured-baseline matches are excluded only when both identity and alignment extent support the claim: >=99.8% marker identity and >=95% query coverage by default. This prevents a short high-identity fragment from being treated as redundant. Diversity is ranked first using marginal tree distance from baseline and committed genomes where available, then baseline/project/reference divergence and phylogenetic isolation. MWL rank and score add priority among otherwise comparable diversity choices but never override the baseline-redundancy gate.

After partner submissions have accumulated and sequencing statuses have been updated, run a Quarterly Review to choose the next tranche across the whole project:

```bash
branchmanager quarterly-review \
  --db project.sqlite \
  --genome-budget 24 \
  --backups-per-primary 1 \
  --tree runs/final/tree/current_tree.nwk \
  --alignment runs/final/tree/current_alignment.fasta \
  --partner-metadata partner_status.tsv \
  -o quarterly_review_01
```

`--genome-budget` is the number of new primary nominations, not including backups. Quarterly Review fills remaining nine-genome coverage gaps first, filters baseline-redundant candidates, then balances diversity across species. Within those tiers it uses marginal tree distance where available, nearest-genome marker identity, cultured-baseline novelty, GTDB-reference context, phylogenetic isolation, and MWL evidence. When `--alignment` is supplied, nearest available-genome identity is recomputed against the current genome collection, so newly completed genomes affect the next round immediately. It does not change `already_sequenced`; recommendations become available genomes only after metadata is updated in a later round.

Every current Performance Review stores its normalised assessment rows in the project database. For an older review made before assessment snapshots were introduced, import its full audit table once:

```bash
branchmanager quarterly-review \
  --db project.sqlite \
  --assessment runs/01_Batch1/assessment/sequence_assessment.tsv \
  --genome-budget 24 \
  --tree runs/01_Batch1/tree/current_tree.nwk \
  -o quarterly_review_01
```

Quarterly Review outputs are:

| File | Description |
|---|---|
| `next_genome_set.tsv` | Current `PRIMARY` nominations and their `BACKUP` isolates |
| `quarterly_review_summary.tsv` | Full audit containing all candidates, already-sequenced isolates, review exclusions, tiers, marginal distances, and reasons |
| `quarterly_review_manifest.tsv` | Round parameters, recommendation counts, and scientific scope |
| `neighbourhoods/clade_*.png` | Round-specific local-clade figures whose stars show this Quarterly Review's primary and backup nominations |

Quarterly Review is marker-gene decision support. Once genomes are available, strain-level expansion should be checked with ANI, phylogenomics, and pangenome structure rather than treating 16S distance as a strain boundary.

### `sequence_assessment.tsv` columns

| Column | Meaning |
|---|---|
| `ID` | Sequence ID supplied by the user; IDs are preserved unless shortening is explicitly requested |
| `GTDBTaxonomy` | Full authoritative GTDB lineage assigned by the classifier |
| `GTDBClassificationHit` | Closest GTDB reference accession used for classification |
| `GTDBClassificationIdentity` | % identity to the classification hit |
| `GTDBClassificationConfidence` | Confidence of the GTDB assignment |
| `BaselineNearestHit` | Closest sequence in the cultured baseline, e.g. Hungate |
| `BaselineNearestHitDataset` | Dataset label for that baseline hit |
| `BaselineNearestIdentity` | % identity to the cultured-baseline hit |
| `BaselineMatchesGE99 / GE97 / GE95` | Count of baseline sequences within 99% / 97% / 95% identity |
| `BaselineNoveltyScore` | 0–100 score; higher = more divergent and less crowded versus cultured isolates |
| `BaselineCrowding` | Baseline neighbourhood density: `crowded` / `moderate` / `sparse` / `isolated` |
| `BaselinePriority` | `HIGH` / `MEDIUM` / `LOW` heuristic from baseline novelty and density |
| `Project*` | Parallel novelty hit, identity, score, crowding, and priority against the rolling partner-candidate collection, excluding self |
| `GTDBReference*` | Parallel nearest hit, identity, score, and crowding against the authoritative GTDB reference |
| `PartnerID` | Partner acronym loaded from `--partner-metadata`, e.g. `QUB` or `UoG` |
| `SelectedForGenomeSequencing` | Whether the project has committed this isolate to sequencing; the genome may still be pending |
| `GenomeAlreadySequenced` | Whether the current partner isolate already has an available genome; a recommendation alone does not set this field |
| `GenomeCollection*` | Nearest hit and density across all baseline genomes plus already-sequenced partner genomes |
| `BaselineGenomesSameAssessmentSpecies` / `SequencedPartnerGenomesSameAssessmentSpecies` | Existing exact-species genome coverage by source; assessment species means GTDB for bacterial/archaeal Performance Reviews |
| `CommittedGenomesSameAssessmentSpecies` / `SelectedPendingGenomesSameAssessmentSpecies` / `PangenomeTarget` / `PangenomeGap` | Rolling completed-genome coverage, pending commitments, and the number of additional selections required |
| `SequencingSet*` | Stable group ID and proposed primary/backup/alternate action for the isolate |
| `LocalNeighbourhoodFigure` | Relative path to the local-clade PNG containing this assessed sequence |
| `TreeContextLeafCount` | Number of tree leaves displayed in that figure; this is display context, not a crowding score |
| `AssessedSequencesInTreeContext` | Number of current assessed isolates sharing the figure |
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

### Partner sequencing metadata

`branchmanager performance-review` requires `--partner-metadata` / `--sequencing-metadata`. Provide this as a simple sidecar `.csv`, `.tsv`, or gzipped CSV/TSV alongside the FASTA. It must contain:

- A sequence ID column such as `sequence_id`, `isolate_id`, `sample_id`, or `ID`. Values may match IDs in the current FASTA or partner isolates already stored in the same project SQLite database.
- A partner acronym column such as `partner_id`, `partner`, or `partner_acronym`, with values like `QUB` or `UoG`.
- An optional selection-commitment column such as `selected_for_genome_sequencing`. Use `yes` after the project has actually committed that isolate to sequencing.
- An already-sequenced/genome-available column such as `already_sequenced`, `genome_available`, or `genome_sequenced`. This column is required and must be explicit `yes` or `no`.

Example:

```tsv
sequence_id	partner_id	selected_for_genome_sequencing	already_sequenced
Iso001	QUB	yes	no
Iso002	UoG	yes	yes
Iso003	UoG	no	no
```

Set `selected_for_genome_sequencing=yes` only after the advisory/project decision has been made. It reserves that isolate in rolling pangenome commitments but does not put it in the nearest-genome comparison pool. Set `already_sequenced=yes` only when a usable genome is genuinely available; that state does enter WGS coverage and nearest-genome calculations. `PRIMARY` and `BACKUP` are BranchManager recommendations and do not update either factual state automatically.

Keep one cumulative ledger for the project. It may contain rows from earlier submissions even when the current FASTA contains only new sequences. Each Performance Review refreshes matched stored isolates before recalculating commitments, available-genome coverage, and recommendations.

---

## ID shortening

By default, BranchManager preserves the FASTA IDs exactly as supplied. This is especially important for Hungate/baseline datasets and partner-provided isolate IDs.

- `--no-shorten-ids` (default) — keep supplied IDs.
- `--shorten-ids` — generate compact IDs when explicitly requested.
- `--no-baseline-shorten-ids` (default for Performance Review baselines) — keep baseline IDs exactly.
- `--baseline-shorten-ids` — generate compact baseline IDs when explicitly requested.

When IDs are shortened, BranchManager writes an ID map under `ids/`. When IDs are preserved, the map is still useful as an audit trail but should normally be identity-to-identity.

---

## Anchor sequences

Anchor sequences are **invisible scaffolding** constraining phylum-level topology during tree construction. They are pruned from the stored newick after each build and never appear in outputs, iTOL files, or novelty scoring.

The bundled anchor set (`src/branchmanager/data/reference_anchors.fasta`) contains 26 NCBI RefSeq 16S sequences covering major gut/rumen phyla. A companion metadata table (`src/branchmanager/data/reference_anchors.tsv`) explains what each anchor represents, its rumen relevance, and whether it is a core rumen anchor or broader topology scaffold.

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

This is a pragmatic rumen/gut scaffold rather than a final rumen-only reference panel. The TSV flags which anchors are core rumen representatives and which are broader topology anchors that can be replaced as better vetted rumen-specific sequences become available.

To use a custom anchor file:
```bash
branchmanager filing-cabinet --anchors /path/to/my_anchors.fasta ...
```
Custom anchor headers must begin with `BRANCHMANAGER_REF_`:
```
>BRANCHMANAGER_REF_MyPhylum accession=NR_XXXXXX source=SILVA138
```

To refresh the bundled anchors from NCBI:
```bash
python scripts/build_anchor_fasta.py --email your@institution.ac.uk
```

---

## Interpreting outputs

### Selection evidence

| Metric | Meaning |
|---|---|
| `CulturedGap` | `LARGE` below 97%, `MODERATE` from 97% to <98.65%, and `SMALL` at >=98.65% against the cultured baseline |
| `ProjectCoverage` | Number of other rolling partner candidates at >=97% identity: uncovered, sparse (1-3), moderate (4-10), or dense (>10) |
| `ReferenceContext` | `HIGH DIVERGENCE` below 94.5%, `MODERATE DIVERGENCE` from 94.5% to <98.65%, or `CLOSE REFERENCE` at >=98.65% |
| `GenomeCoverage` | Exact same-GTDB-species available-genome count and remaining target gap; nearest-genome identity is shown separately as context |
| `EvidenceQuality` | Whether classification identity and placement warnings support using the marker result for selection |

The 94.5% and 98.65% values are commonly used full-length 16S heuristics for genus- and species-range divergence. BranchManager uses them as decision-support categories, not as declarations of a new genus or species. Sequence quality, alignment coverage, marker copy variation, phylogenetic placement, and ultimately genome-based analysis still matter.

Threshold basis: [Yarza et al. (2014), genus and higher-rank 16S boundaries](https://doi.org/10.1038/nrmicro3330) and [Kim et al. (2014), 98.65% 16S/ANI species-demarcation correspondence](https://pubmed.ncbi.nlm.nih.gov/24505072/). Their use in rumen-focused taxonomy is discussed by [Henderson et al. (2019)](https://doi.org/10.7717/peerj.6496), including limitations for unusually diverse genera.

The local-clade images expand around each assessed tree leaf to show a compact context (normally 8-30 leaves). Overlapping contexts are merged when their common clade remains compact. The nearest cultured-baseline hit for each grouped context is forced into the image even when reaching it crosses the normal leaf limit; less informative non-target leaves are pruned first. An orange diamond and `[SEQUENCED]` mark a confirmed available genome. A filled purple star marks the current sequence-next recommendation (`[P1]`, or another primary if more than one genome is missing); outlined purple stars mark recommended backups (`[B2]` and later ranks). `[ALT]` isolates are contextual alternatives, not current recommendations.

Each leaf label reports MSA pident relative to the figure's P1 sequence, and `clade_*_pairwise_pident.tsv` contains all displayed leaf-to-leaf comparisons. MSA pident is defined as identical A/C/G/T bases divided by alignment columns where both sequences contain an unambiguous A/C/G/T base; `ComparableACGTColumns` reports the overlap supporting the value. Terminal overhangs, gaps, and ambiguous bases are excluded rather than being misreported as substitutions. Internal nodes are inferred ancestors without observed marker sequences, so BranchManager does not assign pident to them; branch geometry and the scale bar remain modelled substitutions per site. These plots support visual inspection, while quantitative crowding remains the vsearch `MatchesGE*` count against explicitly named pools.

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
- `combined_taxonomy.tsv` or `filing_cabinet_combined_taxonomy.tsv`
- `taxonomy.tsv` (raw classifier output)
- `taxonomy_input_warnings.tsv`

---

## Operational updates and project close-out

Recommendations are not genome facts. Use `status-meeting` for laboratory progress and `records-update` only when genome evidence exists.

```bash
# Laboratory/SAB lifecycle changes
branchmanager status-meeting --db project.sqlite --input isolate_status.tsv -o updates/status_01

# Completed assemblies, with automatic 90% completeness / 5% contamination defaults
branchmanager records-update --db project.sqlite --input genome_results.tsv -o updates/genomes_01

# Re-open the complete collection after the initial target-nine selection phase
branchmanager quarterly-review --db project.sqlite --genome-budget 24 \
  --from-dir runs/latest/03_performance_review_hiring_panel -o quarterly_reviews/round_02

# Final cumulative report after sequencing and genome QC
branchmanager annual-report --db project.sqlite -o final_annual_report
```

`isolate_status.tsv` requires `sequence_id` and `status`; supported states include `RECEIVED`, `TRACE_REVIEW`, `MARKER_QC_PASSED`, `PROPOSED`, `SAB_APPROVED`, `DNA_EXTRACTION_FAILED`, `SEQUENCED`, `GENOME_QC_FAILED`, and `GENOME_QC_PASSED`.

`genome_results.tsv` requires `sequence_id` and `genome_id`/`accession`. Recommended fields are `genome_status`, `completeness`, `contamination`, `gtdb_taxonomy`, `ani_cluster`, and `genome_path`. A `SEQUENCED` row with both quality metrics is automatically assigned genome-QC pass/fail using the configured thresholds. Explicit QC pass values that conflict with supplied metrics are rejected. Only QC-passed genomes update the available genome collection.

Every Performance Review, Quarterly Review, and Records Update works on a locked staged database and atomically publishes the result after validation. Annual Report writes an HTML dashboard, isolate and genome ledgers, workflow history, and the latest decision-change audit.

Office-stage reports are named `onboarding_report.tsv`, `onboarding_summary.json`, `status_meeting_report.tsv`, `records_update_report.tsv`, `annual_report.html`, `it_desk_report.tsv`, and `it_desk_summary.json`.

Before a production run:

```bash
branchmanager it-desk --db project.sqlite --ref gtdb_ssu_reps.fna.gz \
  --tree-method fasttree --strict -o it_desk
```

## What BranchManager does

BranchManager is aimed at ranking marker-gene sequence evidence for isolate follow-up. It helps identify lineages that may have been missed because of primer bias, sparse reference coverage, taxonomy lag, or conservative filtering. Many targets are not new in nature — they are often **new to the reference record** or **underrepresented in existing collections**.

That is why BranchManager combines reference-aware classification with novelty scoring, Sanger QC, baseline/reference comparisons, partner metadata, and warning layers, rather than treating every long branch as a discovery.
