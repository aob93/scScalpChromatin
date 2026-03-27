#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_FILE="${ROOT_DIR}/scmore_melanoma/scmore.env"

if [[ -f "${ENV_FILE}" ]]; then
  # shellcheck disable=SC1090
  source "${ENV_FILE}"
fi

: "${SCMORE_PROJECT_ROOT:=${ROOT_DIR}}"
: "${SCMORE_RESULTS_ROOT:=${SCMORE_PROJECT_ROOT}/scmore_melanoma/results}"
: "${SCMORE_TRAIT_NAME:=melanoma}"
: "${R_BIN:=Rscript}"

export SCMORE_PROJECT_ROOT
export SCMORE_RESULTS_ROOT
export SCMORE_TRAIT_NAME

mkdir -p "${SCMORE_RESULTS_ROOT}/logs" "${SCMORE_RESULTS_ROOT}/inputs" "${SCMORE_RESULTS_ROOT}/outputs"

run_stage() {
  local stage="$1"
  local script="$2"
  local script_log="${SCMORE_RESULTS_ROOT}/logs/${stage}_$(date +%Y%m%d-%H%M%S).log"

  printf '[%s] %s\n' "${stage}" "${script}"
  printf '[%s] live log: %s\n' "${stage}" "${script_log}"
  "${R_BIN}" "${SCMORE_PROJECT_ROOT}/scmore_melanoma/${script}" 2>&1 | tee "${script_log}"
}

case "${1:-validate}" in
  validate)
    run_stage "VALIDATE" "validate_inputs.R"
    ;;
  prep-fuma)
    run_stage "PREP_FUMA" "00_prepare_fuma_sumstats.R"
    ;;
  prep-gwas)
    run_stage "PREP_GWAS" "01_prepare_sumstats.R"
    ;;
  prep-gene)
    run_stage "PREP_GENE" "02_prepare_gene_results.R"
    ;;
  build-object)
    run_stage "BUILD_OBJECT" "03_export_multimodal_object.R"
    ;;
  run)
    run_stage "RUN_SCMORE" "04_run_scmore.R"
    ;;
  all)
    run_stage "VALIDATE" "validate_inputs.R"
    run_stage "PREP_FUMA" "00_prepare_fuma_sumstats.R"
    run_stage "PREP_GWAS" "01_prepare_sumstats.R"
    run_stage "PREP_GENE" "02_prepare_gene_results.R"
    run_stage "BUILD_OBJECT" "03_export_multimodal_object.R"
    run_stage "RUN_SCMORE" "04_run_scmore.R"
    ;;
  *)
    echo "Usage: $0 [validate|prep-fuma|prep-gwas|prep-gene|build-object|run|all]" >&2
    exit 1
    ;;
esac
