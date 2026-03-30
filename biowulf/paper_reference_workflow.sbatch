#!/usr/bin/env bash
#SBATCH --job-name=scscalp_paper
#SBATCH --cpus-per-task=16
#SBATCH --mem=150g
#SBATCH --time=12:00:00
#SBATCH --gres=lscratch:200
#SBATCH --output=%x-%j.out

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if [[ -f "${PROJECT_ROOT}/biowulf/scscalp.env" ]]; then
  # shellcheck disable=SC1091
  source "${PROJECT_ROOT}/biowulf/scscalp.env"
fi

if [[ -f "${PROJECT_ROOT}/paper_reference_workflow/reference.env" ]]; then
  # shellcheck disable=SC1090
  source "${PROJECT_ROOT}/paper_reference_workflow/reference.env"
fi

module purge
module load "${SCSCALP_BIOWULF_R_MODULE:-R/4.5.2}"

# Optional hook for MACS2 or other environment setup needed by peak calling.
if [[ -n "${SCSCALP_BIOWULF_POST_MODULE_HOOK:-}" ]]; then
  # shellcheck disable=SC1090
  source "${SCSCALP_BIOWULF_POST_MODULE_HOOK}"
fi

export SCSCALP_PROJECT_ROOT="${SCSCALP_PROJECT_ROOT:-${PROJECT_ROOT}}"
export SCSCALP_RESULTS_ROOT="${SCSCALP_RESULTS_ROOT:-${SCSCALP_PROJECT_ROOT}/results}"
export TMPDIR="${TMPDIR:-/lscratch/${SLURM_JOB_ID}}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

mkdir -p "${TMPDIR}"
cd "${SCSCALP_PROJECT_ROOT}"

bash "${SCSCALP_PROJECT_ROOT}/paper_reference_workflow/run_reference_workflow.sh" all
