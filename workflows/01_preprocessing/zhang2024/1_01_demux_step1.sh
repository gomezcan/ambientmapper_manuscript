#!/usr/bin/env bash
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=25
#SBATCH --mem=10G
#SBATCH --job-name=scifi_step1_01
#SBATCH --partition=standard
#SBATCH --output=_logs/%x_%A_%a.log
#SBATCH --array=0
#
# scifi-demux step 1, plan: split one sequencing run into 100 chunks for the array job 1_02_demux_step1.sh.
# Run once per SRA run (library id scifi_<dataset>_rep<N>_<run>; nine runs for the Zhang et al. 2024 libraries).
# Usage (from ${PROJECT_ROOT}/2_CleanReads):
#   LIB=scifi_B73Mo17_rep1_1 DESIGN=<repo>/config/PlateDesign_scifi_B73Mo17_rep1.txt sbatch 1_01_demux_step1.sh

set -euo pipefail
: "${PROJECT_ROOT:?set PROJECT_ROOT to the analysis tree that holds 1_RawData/ and 2_CleanReads/}"
REPO_ROOT="${REPO_ROOT:-$(git rev-parse --show-toplevel 2>/dev/null || pwd)}"
# conda activate <env from environment.yml>
cd "${PROJECT_ROOT}/2_CleanReads"
mkdir -p _logs

THREADS="${SLURM_CPUS_PER_TASK:-1}"
INPUT_DIR="${PROJECT_ROOT}/1_RawData"
LIB="${LIB:?e.g. scifi_B73Mo17_rep1_1 (one sequencing run of one replicate)}"
DESIGN="${DESIGN:?plate design, e.g. ${REPO_ROOT}/config/PlateDesign_scifi_B73Mo17_rep1.txt}"

# work root of this run
WORKDIR="01_${LIB}"
PLAN="${WORKDIR}/run_plan.step1.chunks.tsv"

[[ -d "$INPUT_DIR" ]] || { echo "Missing INPUT_DIR: $INPUT_DIR" >&2; exit 2; }
[[ -f "$DESIGN"   ]]  || { echo "Missing DESIGN file: $DESIGN" >&2; exit 3; }
mkdir -p "$WORKDIR"

if [[ ! -f "$PLAN" ]]; then
  echo "[$(date)] Plan not found, generating at: $PLAN"
  scifi-demux step1 plan \
    --library "$LIB" \
    --raw-dir "$INPUT_DIR" \
    --chunks 100 \
    --work-root "$WORKDIR"
  [[ -f "$PLAN" ]] || { echo "Plan was not created at $PLAN" >&2; exit 4; }
else
  echo "[$(date)] Using existing plan: $PLAN"
fi
