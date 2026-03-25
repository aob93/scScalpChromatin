#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_FILE="${ROOT_DIR}/biowulf/scscalp.env"

if [[ -f "${ENV_FILE}" ]]; then
  # shellcheck disable=SC1090
  source "${ENV_FILE}"
fi

: "${SCSCALP_PROJECT_ROOT:=${ROOT_DIR}}"
: "${SCSCALP_RESULTS_ROOT:=${SCSCALP_PROJECT_ROOT}/results}"
: "${R_BIN:=Rscript}"

export SCSCALP_PROJECT_ROOT
export SCSCALP_RESULTS_ROOT

RNA_STEPS=(
  "rna_preprocessing/01_scRNA_preprocess.R"
  "rna_preprocessing/02_preclustering.R"
  "rna_preprocessing/03_precluster_qc.R"
  "rna_preprocessing/04_precluster_markers.R"
  "rna_preprocessing/05_reclustering.R"
  "rna_preprocessing/06_cluster_qc.R"
  "rna_preprocessing/07_cluster_markers.R"
  "rna_preprocessing/08_assign_RNA_clusters.R"
  "rna_preprocessing/09a_subgroup_clustering.R"
  "rna_preprocessing/09b_Fibroblast_cluster_markers.R"
  "rna_preprocessing/09b_HF_Kc_cluster_markers.R"
  "rna_preprocessing/09b_endothelial_cluster_markers.R"
  "rna_preprocessing/09b_keratinocyte_cluster_markers.R"
  "rna_preprocessing/09b_lymphoid_cluster_markers.R"
  "rna_preprocessing/09b_myeloid_cluster_markers.R"
  "rna_preprocessing/10_remerged_RNA.R"
)

ATAC_STEPS=(
  "atac_preprocessing/01_prepare_archr_proj.R"
  "atac_preprocessing/02a_subgroup_archr.R"
  "atac_preprocessing/02b_recluster_Endothelial.R"
  "atac_preprocessing/02b_recluster_Fibroblasts.R"
  "atac_preprocessing/02b_recluster_Keratinocytes.R"
  "atac_preprocessing/02b_recluster_Lymphoid.R"
  "atac_preprocessing/02b_recluster_Myeloid.R"
  "atac_preprocessing/03_call_peaks_archr.R"
  "atac_preprocessing/03b_p2gLinks_Endothelial.R"
  "atac_preprocessing/03b_p2gLinks_Fibroblasts.R"
  "atac_preprocessing/03b_p2gLinks_Keratinocytes.R"
  "atac_preprocessing/03b_p2gLinks_Lymphoid.R"
  "atac_preprocessing/03b_p2gLinks_Myeloid.R"
  "atac_preprocessing/04_full_RNA_integration.R"
  "atac_preprocessing/04b_Keratinocytes_analyses.R"
  "atac_preprocessing/04c_HF_Kc_analysis.R"
)

run_steps() {
  local stage="$1"
  shift
  local stage_dir
  stage_dir="$(printf '%s' "${stage}" | tr '[:upper:]' '[:lower:]')"
  local log_dir="${SCSCALP_RESULTS_ROOT}/pipeline_logs/${stage_dir}"
  mkdir -p "${log_dir}"

  for script in "$@"; do
    local script_name
    script_name="$(basename "${script}" .R)"
    local script_log="${log_dir}/${script_name}_$(date +%Y%m%d-%H%M%S).log"

    printf '[%s] %s\n' "${stage}" "${script}"
    printf '[%s] live log: %s\n' "${stage}" "${script_log}"
    "${R_BIN}" "${SCSCALP_PROJECT_ROOT}/${script}" 2>&1 | tee "${script_log}"
  done
}

case "${1:-all}" in
  rna)
    run_steps "RNA" "${RNA_STEPS[@]}"
    ;;
  atac)
    run_steps "ATAC" "${ATAC_STEPS[@]}"
    ;;
  all)
    run_steps "RNA" "${RNA_STEPS[@]}"
    run_steps "ATAC" "${ATAC_STEPS[@]}"
    ;;
  *)
    echo "Usage: $0 [rna|atac|all]" >&2
    exit 1
    ;;
esac
