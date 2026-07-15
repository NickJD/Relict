#!/usr/bin/env bash
set -euo pipefail

source /Users/nicholas/anaconda3/bin/activate /Users/nicholas/anaconda3/envs/MWL

BM="/Users/nicholas/Git/BranchManager/src/branchmanager/cli.py"
RUNS="/Users/nicholas/Nextcloud/Current_Work/RGW/16s/workflow/BranchManager_Bac_runs"
DB="$RUNS/project.sqlite"
PROJECT_META="$RUNS/project_partner_metadata.tsv"

GTDB="/Users/nicholas/Nextcloud/Current_Work/RGW/16s/workflow/Taxa_Databases/GTDB/gtdb_r232_bac_arch_ssu_reps.fna.gz"
GG2="/Users/nicholas/Nextcloud/Current_Work/RGW/16s/workflow/Taxa_Databases/GG2/2024.09.backbone.full-length.fna.gz"
GG2_TAXA="/Users/nicholas/Nextcloud/Current_Work/RGW/16s/workflow/Taxa_Databases/GG2/2024.09.taxonomy.id.tsv.gz"
HUNGATE="/Users/nicholas/Nextcloud/Current_Work/RGW/16s/workflow/Taxa_Databases/Hungate/Sanger_Hun.fasta"
MWL="/Users/nicholas/Nextcloud/Current_Work/RGW/16s/MWL/MWL.csv"

export MPLCONFIGDIR="$RUNS/.matplotlib"
mkdir -p "$RUNS" "$MPLCONFIGDIR"

if [[ ! -s "$PROJECT_META" ]]; then
  printf '%s\n' \
    "Create $PROJECT_META before running." \
    "Required header:" \
    $'sequence_id\tpartner_id\tselected_for_genome_sequencing\talready_sequenced' >&2
  exit 2
fi

python "$BM" it-desk --ref "$GTDB" --ref "$GG2" \
  --tree-method fasttree --strict -o "$RUNS/00_it_desk"

if [[ ! -s "$DB" ]]; then
  python "$BM" filing-cabinet \
    --fasta "$HUNGATE" --db "$DB" --dataset Hungate \
    --ref "$GTDB" --ref-name GTDB \
    --alt-ref "$GG2" --alt-taxa "$GG2_TAXA" --alt-ref-name GG2 \
    --main-ref GTDB --classify --build-tree \
    --sequence-domain bacteria --threads 10 \
    -o "$RUNS/00_filing_cabinet"
fi

performance_review() {
  local dataset="$1"
  local marker_fasta="$2"
  local previous="$3"
  local marker_qc="$4"
  local accept_unverified="$5"
  local out="$RUNS/$dataset"
  local qc_args=()

  if [[ -n "$marker_qc" ]]; then
    qc_args+=(--marker-qc "$marker_qc")
  elif [[ "$accept_unverified" == "yes" ]]; then
    qc_args+=(--accept-unverified-marker-qc)
  else
    printf 'No marker-QC evidence or explicit acceptance supplied for %s\n' "$dataset" >&2
    return 2
  fi

  python "$BM" performance-review \
    --input "$marker_fasta" "${qc_args[@]}" \
    --partner-metadata "$PROJECT_META" \
    --db "$DB" --dataset "$dataset" \
    --ref "$GTDB" --ref-name GTDB \
    --alt-ref "$GG2" --alt-taxa "$GG2_TAXA" --alt-ref-name GG2 \
    --main-ref GTDB --chimera-ref "$GTDB" --mwl "$MWL" \
    --previous-review "$previous" --sequence-domain bacteria \
    --pangenome-target 3 --candidate-set-size 4 \
    --neighbourhood-format png --threads 10 \
    -o "$out/03_performance_review_hiring_panel"
}

process_ab1_partner() {
  local dataset="$1"
  local map="$2"
  local reads="$3"
  local previous="$4"
  local out="$RUNS/$dataset"

  python "$BM" onboarding \
    --sample-map "$map" --partner-metadata "$PROJECT_META" \
    --read-dir "$reads" -o "$out/01_onboarding"

  python "$BM" paper-trail \
    --sample-map "$out/01_onboarding/normalised_read_map.tsv" \
    --screen-ref "$GTDB" --min-quality 20 --min-mean-quality 25 \
    --min-read-length 300 --min-length 800 --min-overlap 40 \
    --threads 10 -o "$out/02_paper_trail_merge_meeting"

  performance_review "$dataset" \
    "$out/02_paper_trail_merge_meeting/assembled.fasta" "$previous" \
    "$out/02_paper_trail_merge_meeting/assembly_report.tsv" no
}

process_fasta_partner() {
  local dataset="$1"
  local fasta="$2"
  local previous="$3"
  local marker_qc="${4:-}"
  local out="$RUNS/$dataset"

  python "$BM" onboarding \
    --fasta "$fasta" --partner-metadata "$PROJECT_META" \
    -o "$out/01_onboarding"

  # With no marker-QC sidecar this is an explicit audited acceptance. The
  # assessment still records the marker as quality-unverified.
  performance_review "$dataset" \
    "$out/01_onboarding/normalised_input.fasta" "$previous" \
    "$marker_qc" yes
}

# Example AB1 submission:
# process_ab1_partner "UoG_01" \
#   "/path/to/UoG/ab1_mapping.csv" "/path/to/UoG/All_AB1" \
#   "$RUNS/00_filing_cabinet"

# Example FASTA-only submission after UoG_01:
# process_fasta_partner "QUB_01" "/path/to/QUB_01_16S.fasta.gz" \
#   "$RUNS/UoG_01/03_performance_review_hiring_panel"

# When all desired submissions have been processed:
# python "$BM" annual-report --db "$DB" -o "$RUNS/current_annual_report"
