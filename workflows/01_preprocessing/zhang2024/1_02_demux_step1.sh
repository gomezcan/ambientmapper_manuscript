#!/usr/bin/env bash
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=10G
#SBATCH --job-name=scifi_step1_run
#SBATCH --partition=standard
#SBATCH --output=_logs/%x_%A_%a.log
#SBATCH --array=1-100
#
# scifi-demux step 1, run: demultiplex and barcode-correct one chunk of one sequencing run per array task
# (16 bp 10x barcode from the index read into the read names; 5 bp + 5 bp Tn5 well barcodes corrected to the
# closest well with at most one mismatch per half; reads written per well according to the plate design).
# Usage (from ${PROJECT_ROOT}/2_CleanReads, after 1_01_demux_step1.sh with the same LIB and DESIGN):
#   LIB=scifi_B73Mo17_rep1_1 DESIGN=<repo>/config/PlateDesign_scifi_B73Mo17_rep1.txt sbatch 1_02_demux_step1.sh

set -euo pipefail
: "${PROJECT_ROOT:?set PROJECT_ROOT to the analysis tree that holds 1_RawData/ and 2_CleanReads/}"
REPO_ROOT="${REPO_ROOT:-$(git rev-parse --show-toplevel 2>/dev/null || pwd)}"
# conda activate <env from environment.yml>
cd "${PROJECT_ROOT}/2_CleanReads"
mkdir -p _logs

# --- params ---
THREADS="${SLURM_CPUS_PER_TASK:-1}"
INPUT_DIR="${PROJECT_ROOT}/1_RawData"
LIB="${LIB:?e.g. scifi_B73Mo17_rep1_1 (one sequencing run of one replicate)}"
WORKDIR="01_${LIB}"                    # must match step 1
PLAN="${WORKDIR}/run_plan.step1.chunks.tsv"
DESIGN="${DESIGN:?plate design, e.g. ${REPO_ROOT}/config/PlateDesign_scifi_B73Mo17_rep1.txt}"

[[ -d "$INPUT_DIR" ]] || { echo "Missing INPUT_DIR: $INPUT_DIR" >&2; exit 2; }
[[ -f "$DESIGN"   ]]  || { echo "Missing DESIGN file: $DESIGN" >&2; exit 3; }
[[ -f "$PLAN"     ]]  || { echo "Missing PLAN file: $PLAN" >&2; exit 4; }

# derive chunk count from plan (assumes header row)
CHUNKS=$(tail -n +2 "$PLAN" | wc -l | awk '{print $1}')
[[ "$CHUNKS" -ge 1 ]] || { echo "No chunks in $PLAN" >&2; exit 5; }

# --- array index mapping: plan is 1-based; SLURM is already set to 1-CHUNKS above ---
IDX="${SLURM_ARRAY_TASK_ID}"

echo "[$(date)] step1.run idx=${IDX}/${CHUNKS} threads=${THREADS}"

# --- execute ---
scifi-demux step1 run \
  --mode hpc \
  --library "$LIB" \
  --raw-dir "$INPUT_DIR" \
  --design "$DESIGN" \
  --threads "$THREADS" \
  --work-root "$WORKDIR" \
  --chunks 100
