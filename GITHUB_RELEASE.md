# BranchManager First Working Prototype: Marker-Gene QC, Novelty Scoring, and Isolate Prioritisation

BranchManager is a substantially expanded marker-gene based candidate-selection system for isolate collections: Sanger/AB1 QC, taxonomy against multiple 16S/marker databases, novelty scoring against baseline and project collections, MWL matching, tree/iTOL output, and concise candidate-selection summaries.

## Breaking Changes

- The toolkit is now branded **BranchManager** throughout.
- Python imports use the `branchmanager.*` namespace.
- CLI examples and help now use `branchmanager`.
- Reference anchor headers use the `BRANCHMANAGER_REF_` prefix.
- Run logs are now written as `branchmanager.log`.
- IDs are preserved by default; generated shortened IDs are only used when requested.
- `preload` and `evaluate` default to `--sequence-domain bacteria`; use `--sequence-domain archaea`, `fungi`, or `mixed` as needed.
- The Sanger/AB1 workflow now applies stricter QC by default and withholds `FAIL_QC` sequences from `assembled.fasta`.

## Highlights

- Added a core `evaluate` workflow for rolling partner-sequence assessment and candidate prioritisation.
- Added baseline loading for cultured/reference collections such as Hungate via `--baseline-fasta` and `--baseline-dataset`.
- Added multi-database taxonomic reporting for GTDB plus alternate references such as GG2, SILVA, and NCBI.
- Added MWL matching using GTDB assignments as the authoritative taxonomy layer.
- Added partner metadata support with partner acronym and selected-for-genome-sequencing status.
- Added rolling WGS-selection context so close clades already represented by selected genomes can be deprioritised.
- Added separate novelty metrics against baseline datasets, all known project datasets, and the selected external reference.
- Added concise `selection_summary.tsv` output for scientific advisory board review.
- Added detailed `sequence_assessment.tsv` audit output for traceability.

## Sanger / AB1 Processing

- Added `branchmanager sanger` / `branchmanager ab1` / `branchmanager ab1-to-fasta`.
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
- Added reproducible `sanger_qc_policy.tsv`.
- Added SVG visual reports for read quality/trimming and isolate assembly.

## Taxonomy and Preclassification

- Added reusable `preclassify` workflow for reference/baseline datasets.
- Added high-quality taxonomic disagreement reports where one query has multiple plausible conflicting hits.
- Taxonomy assignments can now be read from FASTA headers or external CSV/TSV files.
- FASTA, CSV, and TSV taxonomy inputs can be gzipped.
- Added stronger reference/taxonomy consistency checks.
- Added domain-aware filtering for bacteria, archaea, fungi, and mixed/all runs.

## Novelty and Candidate Scoring

- Added baseline novelty scoring, usually against cultured rumen isolates such as Hungate.
- Added all-known novelty scoring against all non-current datasets in the project database.
- Added external reference novelty scoring against the chosen primary reference, usually GTDB.
- Added selected-genome clade context using partner metadata.
- Added adjusted genome-sequencing priority when a close 16S clade already has a selected genome.
- Added crowding/density metrics across novelty lenses.

## Outputs and Reporting

- Reorganised outputs into clearer directories for assessment, taxonomy, baseline hits, tree/MSA files, iTOL files, IDs, logs, and intermediates.
- Removed redundant iTOL symbol and tree-colour files; BranchManager now emits one colour-strip style output per metadata type.
- Added `OUTPUT_EXPLANATIONS.tsv` and `OUTPUT_GUIDE.md`.
- Added per-sequence QC rejection reasons.
- Added `selection_summary.tsv` as the board-facing candidate-decision report.
- Preserved baseline and partner FASTA IDs by default.

## Tree and iTOL

- Maintained anchored MAFFT/FastTree workflows while pruning anchor leaves from final trees.
- Added domain/profile-aware preload/evaluate filtering.
- Improved tree metadata outputs for iTOL viewing.
- Kept current tree and alignment files available for external inspection.

## Validation

- Focused test suite passes: `67 passed`.
- CLI smoke test passes: `python src/branchmanager/cli.py --help`.
- Sanger workflow was exercised on a 120-read AB1 batch and produced reproducible QC, assembly, and resequencing recommendation outputs.

## Upgrade Notes

- Use `branchmanager ...` as the only CLI spelling.
- Use `branchmanager.*` as the only Python import namespace.
- Use `BRANCHMANAGER_REF_...` for custom anchor headers.
- Review Sanger outputs carefully after upgrading: stricter QC may classify previously retained reads as `MANUAL_REVIEW` or `RESEQUENCE`.
