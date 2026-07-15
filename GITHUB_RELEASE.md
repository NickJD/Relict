# BranchManager First Working Prototype: Marker-Gene QC, Novelty Scoring, and Isolate Prioritisation

BranchManager is a substantially expanded marker-gene based candidate-selection system for isolate collections: Sanger/AB1 QC, taxonomy against multiple 16S/marker databases, novelty scoring against baseline and project collections, MWL matching, tree/iTOL output, and concise candidate-selection summaries.

## Breaking Changes

- The toolkit is now branded **BranchManager** throughout.
- Python imports use the `branchmanager.*` namespace.
- CLI examples and help now use `branchmanager`.
- Reference anchor headers use the `BRANCHMANAGER_REF_` prefix.
- Run logs are now written as `branchmanager.log`.
- IDs are preserved by default; generated shortened IDs are only used when requested.
- Filing Cabinet and Performance Review default to `--sequence-domain bacteria`; use `--sequence-domain archaea`, `fungi`, or `mixed` as needed.
- The Sanger/AB1 workflow now applies stricter QC by default and withholds `FAIL_QC` sequences from `assembled.fasta`.

## Highlights

- Local-clade PNG figures now distinguish confirmed sequenced isolates from primary and backup recommendations; force nearest cultured-baseline hits into bounded contexts; annotate leaf pident relative to P1; and include complete per-figure pairwise MSA-pident tables.

- Added the core **Performance Review** (`performance-review`) for rolling partner-sequence assessment and candidate prioritisation.
- Added baseline loading for cultured/reference collections such as Hungate via `--baseline-fasta` and `--baseline-dataset`.
- Added multi-database taxonomic reporting for GTDB plus alternate references such as GG2, SILVA, and NCBI.
- Added MWL matching using GTDB assignments as the authoritative taxonomy layer.
- Added a cumulative project metadata ledger with separate partner acronym, selected-for-sequencing commitment, and confirmed already-sequenced/genome-available status.
- Onboarding now accepts either an AB1/sample map or a partner-supplied FASTA, including gzipped FASTA. FASTA submissions retain the same ID and metadata checks while bypassing chromatogram-only stages.
- Added explicit rolling collections: cultured baseline genomes, all partner candidates, and the genome collection formed from every baseline plus partner isolates with an available genome.
- Added separate novelty metrics against the cultured baseline, rolling partner collection, and GTDB external reference.
- Added exact same-GTDB-species genome counts with a configurable pangenome target (`--pangenome-target`, default 3).
- Added `sequencing_sets.tsv`: group-level primary and backup proposals, normally up to four candidates per GTDB species or unresolved local clade.
- Candidate sets use phylogenetic spread where branch lengths are available; primary rows fill the pangenome gap and backup rows provide extraction-failure resilience plus strain-level diversity.
- Added concise `selection_summary.tsv` output for scientific advisory board review.
- Added reproducible **Quarterly Review** rounds for selecting a later project-wide genome tranche after the default three-genome species coverage phase, with an explicit genome budget, backups, assessment snapshots, and immutable factual genome status.
- Added detailed `sequence_assessment.tsv` audit output for traceability.
- Added grouped local-clade images for every tree-resolved assessed isolate, with full leaf labels, baseline/project context, orange already-sequenced markers, filled primary-recommendation stars, outlined backup stars, and a sequence-to-figure manifest.
- All visual reports now use PNG, including chromatograms, read-error profiles, assembly overviews, and local-clade figures.
- Paper Trail visual reports are height-bounded and automatically paginated, with a TSV page manifest linking each image to its read or isolate range.
- Consolidated the former per-cluster CSV directory into one `assessment/clusters.csv` containing all clusters and members.
- Reworked `selection_summary.tsv` around transparent component evidence, exact genome coverage, and candidate-set role instead of a one-hit adjusted priority.
- Broad MWL matches at phylum/class level remain context; family/genus/species matches can contribute directly to selection priority.

## Paper Trail / Merge Meeting

- Added `branchmanager paper-trail` for AB1 conversion, error trimming, review, and multi-primer assembly.
- Supports AB1/ABI, FASTA, and FASTQ inputs, including gzipped files.
- Supports read-level metadata and one-row-per-isolate sample maps.
- Supports per-isolate handling modes:
  - `assemble` / `merge` / `consensus`
  - `best_read` / `highest_quality` / `select_best` / `independent`
- Added Phred/Mott-style end trimming.
- Added internal low-quality masking to `N`.
- Added expected-error metrics for reads and final outputs.
- Added quality-weighted consensus for overlapping primer reads.
- Added final QC classes: `PASS_HIGH_CONFIDENCE`, `PASS_WITH_WARNINGS`, and `FAIL_QC`.
- Added `resequence_recommendations.tsv` with `ACCEPT`, `MANUAL_REVIEW`, and `RESEQUENCE` decisions.
- Added `failed_qc_sequences/` so final failed outputs and failed reads are retained as FASTA plus a manifest for manual review.
- Added reproducible `paper_trail_qc_policy.tsv`.
- Added PNG visual reports for read quality/trimming and isolate assembly.

## Taxonomy and Background Checks

- Added reusable **Background Check** (`background-check`) for reference/baseline datasets.
- Added high-quality taxonomic disagreement reports where one query has multiple plausible conflicting hits.
- Taxonomy assignments can now be read from FASTA headers or external CSV/TSV files.
- FASTA, CSV, and TSV taxonomy inputs can be gzipped.
- Added stronger reference/taxonomy consistency checks.
- Added domain-aware filtering for bacteria, archaea, fungi, and mixed/all runs.

## Novelty and Candidate Scoring

- Added baseline novelty scoring, usually against cultured rumen isolates such as Hungate.
- Added project novelty scoring against all partner candidates accumulated across runs, excluding each query's self-hit.
- Added external reference novelty scoring against the chosen primary reference, usually GTDB.
- Every supplied baseline isolate now counts as genome available.
- Selected partner isolates reserve planned pangenome slots; only completed/available partner genomes enter nearest-genome and WGS-coverage comparisons.
- Added pangenome coverage and candidate-set roles without collapsing the evidence into an adjusted novelty score.
- Added crowding/density metrics across novelty lenses.

## Outputs and Reporting

- Reorganised outputs into clearer directories for assessment, taxonomy, baseline hits, tree/MSA files, iTOL files, IDs, logs, and intermediates.
- Replaced legacy stage-branded filenames with the office workflow vocabulary, including `it_desk_report.tsv`, `onboarding_report.tsv`, `performance_review_dashboard.html`, `quarterly_review_summary.tsv`, and `annual_report.html`.
- Removed redundant iTOL symbol and tree-colour files; BranchManager now emits one colour-strip style output per metadata type.
- Added `OUTPUT_EXPLANATIONS.tsv` and `OUTPUT_GUIDE.md`.
- Added per-sequence QC rejection reasons.
- Added `selection_summary.tsv` as the board-facing candidate-decision report.
- Preserved baseline and partner FASTA IDs by default.

## Tree and iTOL

- Maintained anchored MAFFT/FastTree workflows while pruning anchor leaves from final trees.
- Added domain/profile-aware Filing Cabinet and Performance Review filtering.
- Improved tree metadata outputs for iTOL viewing.
- Kept current tree and alignment files available for external inspection.

## Validation

- Full test suite passes: `103 passed` plus office-command parser subtests.
- CLI smoke test passes: `python src/branchmanager/cli.py --help`.
- Sanger workflow was exercised on a 120-read AB1 batch and produced reproducible QC, assembly, and resequencing recommendation outputs.

## Upgrade Notes

- Use `branchmanager ...` as the only CLI spelling.
- Use `branchmanager.*` as the only Python import namespace.
- Use `BRANCHMANAGER_REF_...` for custom anchor headers.
- Review Paper Trail outputs carefully after upgrading: stricter QC may classify previously retained reads as `MANUAL_REVIEW` or `RESEQUENCE`.
