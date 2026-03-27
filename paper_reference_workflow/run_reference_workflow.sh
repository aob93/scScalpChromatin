#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_FILE="${ROOT_DIR}/paper_reference_workflow/reference.env"

if [[ -f "${ENV_FILE}" ]]; then
  # shellcheck disable=SC1090
  source "${ENV_FILE}"
fi

: "${SCSCALP_PROJECT_ROOT:=${ROOT_DIR}}"
: "${SCSCALP_RESULTS_ROOT:=${SCSCALP_PROJECT_ROOT}/results}"
: "${R_BIN:=Rscript}"

export SCSCALP_PROJECT_ROOT
export SCSCALP_RESULTS_ROOT

BASELINE_STEPS=(
  "atac_preprocessing/01_prepare_archr_proj.R"
)

SUBCLUSTER_STEPS=(
  "atac_preprocessing/02a_subgroup_archr.R"
  "atac_preprocessing/02b_recluster_Endothelial.R"
  "atac_preprocessing/02b_recluster_Fibroblasts.R"
  "atac_preprocessing/02b_recluster_Keratinocytes.R"
  "atac_preprocessing/02b_recluster_Lymphoid.R"
  "atac_preprocessing/02b_recluster_Myeloid.R"
)

PEAK_STEPS=(
  "atac_preprocessing/03_call_peaks_archr.R"
  "atac_preprocessing/03b_p2gLinks_Endothelial.R"
  "atac_preprocessing/03b_p2gLinks_Fibroblasts.R"
  "atac_preprocessing/03b_p2gLinks_Keratinocytes.R"
  "atac_preprocessing/03b_p2gLinks_Lymphoid.R"
  "atac_preprocessing/03b_p2gLinks_Myeloid.R"
)

INTEGRATION_STEPS=(
  "atac_preprocessing/04_full_RNA_integration.R"
)

ADVANCED_STEPS=(
  "atac_preprocessing/04b_Keratinocytes_analyses.R"
  "atac_preprocessing/04c_HF_Kc_analysis.R"
)

run_validate() {
  local stage="${1:-all}"
  "${R_BIN}" "${ROOT_DIR}/paper_reference_workflow/validate_inputs.R" "${stage}"
}

run_steps() {
  local stage="$1"
  shift
  local stage_dir
  stage_dir="$(printf '%s' "${stage}" | tr '[:upper:]' '[:lower:]')"
  local log_dir="${SCSCALP_RESULTS_ROOT}/pipeline_logs/reference_workflow/${stage_dir}"
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
  validate)
    run_validate "${2:-all}"
    ;;
  baseline)
    run_validate "baseline"
    run_steps "BASELINE" "${BASELINE_STEPS[@]}"
    ;;
  subcluster)
    run_validate "subcluster"
    run_steps "SUBCLUSTER" "${SUBCLUSTER_STEPS[@]}"
    ;;
  peaks)
    run_validate "peaks"
    run_steps "PEAKS" "${PEAK_STEPS[@]}"
    ;;
  integrate)
    run_validate "integrate"
    run_steps "INTEGRATE" "${INTEGRATION_STEPS[@]}"
    ;;
  advanced)
    run_validate "advanced"
    run_steps "ADVANCED" "${ADVANCED_STEPS[@]}"
    ;;
  all)
    run_validate "all"
    run_steps "BASELINE" "${BASELINE_STEPS[@]}"
    run_steps "SUBCLUSTER" "${SUBCLUSTER_STEPS[@]}"
    run_steps "PEAKS" "${PEAK_STEPS[@]}"
    run_steps "INTEGRATE" "${INTEGRATION_STEPS[@]}"
    ;;
  *)
    echo "Usage: $0 [validate|baseline|subcluster|peaks|integrate|advanced|all]" >&2
    exit 1
    ;;
esac
