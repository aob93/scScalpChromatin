#!/usr/bin/env bash
#SBATCH --job-name=scscalp
#SBATCH --cpus-per-task=16
#SBATCH --mem=192g
#SBATCH --time=48:00:00
#SBATCH --output=%x-%j.out

set -euo pipefail

module purge
module load R/4.3.3

if [[ -f "$PWD/biowulf/scscalp.env" ]]; then
  # shellcheck disable=SC1091
  source "$PWD/biowulf/scscalp.env"
fi

export SCSCALP_PROJECT_ROOT="${SCSCALP_PROJECT_ROOT:-$PWD}"
export SCSCALP_RESULTS_ROOT="${SCSCALP_RESULTS_ROOT:-${SCSCALP_PROJECT_ROOT}/results}"

bash "${SCSCALP_PROJECT_ROOT}/biowulf/run_pipeline.sh" all
