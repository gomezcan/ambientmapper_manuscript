#!/usr/bin/env bash
#SBATCH --job-name=syn_decontam_nod
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:45:00
#SBATCH --array=0-29
#SBATCH --output=_logs/03_22_syn_decontam_nod_%A_%a.log
#
# 03_22_decontam_synthetic_nodesign.sh — without-design decontam for both
# synthetic tracks (Track B = synthetic/, Track B-disc = synthetic_disc/).
#
# Goal: characterize decontam behavior on synthetic data with known ground
# truth (alpha titration, known contaminant). 15 datasets x 2 tracks = 30 tasks.
#
# Source cells_calls: the Phase 2 factorial C0 baseline (run tag factorial_phase2_2026-04-09).
#   genotyping_runs/factorial_phase2_2026-04-09/C0/<sample>_cells_calls.tsv.gz
# C0 chosen because Phase 2 evidence (eval/phase2_2026-04-09/) shows it gives
# sens_strict=100% on Track B-disc at alpha=0; xmap-on configs broke things
# on contaminated multi-genotype data. C0 = clean baseline an external user
# would get out-of-the-box.
#
# Decontam params mirror the SM2v2 (01_10b) and B73Mo17 (04_06) without-design
# blocks: alpha=0.05, top1_rescue, doublet=top1, indist=top1, safe-keep-delta=3,
# min_reads_post_clean=100, min_allowed_frac_post_clean=0.90.
#
# Output: <track>/<dataset>/decontam_without_design_alpha05_C0/
#         (the pre/post barcode_genome_counts tables are read by Fig. 4C)
#
# Array index decode:
#   TRACK_IDX = TASK_ID / N_DATASETS    [0..1]
#   DS_IDX    = TASK_ID % N_DATASETS    [0..14]
#
# Walltime budget: 3 chunks/dataset, 208 BCs/dataset; per-task <5 min real,
# 45 min generous to absorb queue/IO jitter.
# =============================================================================

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"   # root of the data tree (3_Mapping/, 5_AmbientDetection/, ...)

cd "${PROJECT_ROOT}/5_AmbientDetection"

# conda activate <env from environment.yml>

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export GOTO_NUM_THREADS=1

# -----------------------------------------------------------------------------
# Dispatch
# -----------------------------------------------------------------------------
DATASETS=(
    alpha_000
    alpha_002_Il14H alpha_005_Il14H alpha_010_Il14H alpha_020_Il14H
    alpha_030_Il14H alpha_040_Il14H alpha_050_Il14H
    alpha_002_Ki11  alpha_005_Ki11  alpha_010_Ki11  alpha_020_Ki11
    alpha_030_Ki11  alpha_040_Ki11  alpha_050_Ki11
)
TRACKS=(synthetic synthetic_disc)

N_DATASETS=${#DATASETS[@]}   # 15
N_TRACKS=${#TRACKS[@]}       # 2
N_TOTAL=$(( N_DATASETS * N_TRACKS ))   # 30

TASK_ID=${SLURM_ARRAY_TASK_ID:-0}
if [[ ${TASK_ID} -ge ${N_TOTAL} ]]; then
    echo "ERROR: TASK_ID=${TASK_ID} out of range (max ${N_TOTAL})" >&2
    exit 1
fi

TRACK_IDX=$(( TASK_ID / N_DATASETS ))
DS_IDX=$(( TASK_ID % N_DATASETS ))

TRACK="${TRACKS[$TRACK_IDX]}"
DATASET="${DATASETS[$DS_IDX]}"

# -----------------------------------------------------------------------------
# Source genotyping config + paths
# -----------------------------------------------------------------------------
CONFIG="C0"
FACTORIAL_TAG="factorial_phase2_2026-04-09"

DS_DIR="${TRACK}/${DATASET}"
GENO_DIR="${DS_DIR}/genotyping_runs/${FACTORIAL_TAG}/${CONFIG}"
CELLS_CALLS="${GENO_DIR}/${DATASET}_cells_calls.tsv.gz"
CHUNKS_DIR="${DS_DIR}/cell_map_ref_chunks"
ASSIGN_GLOB="${CHUNKS_DIR}/${DATASET}_chunk*_filtered.tsv.gz"

OUTDIR="${DS_DIR}/decontam_without_design_alpha05_C0"

threads="${SLURM_CPUS_PER_TASK:-8}"

# Decontam params (mirror SM2v2 / B73Mo17 without-design blocks)
alpha="0.05"
chunksize="1000000"
ambiguous_policy="top1_rescue"
doublet_policy="top1"
indist_policy="top1"
safe_keep_delta_as="3"
min_reads_post_clean="100"
min_allowed_frac_post_clean="0.90"
# Both post-clean thresholds are reporting flags in the AmbientMapper build used
# here, not enforced gates (see the stage README).

echo "============================================================"
echo "[$(date)] 03_22 synthetic decontam (no design) — task ${TASK_ID}"
echo "  track        : ${TRACK}"
echo "  dataset      : ${DATASET}"
echo "  cells_calls  : ${CELLS_CALLS}"
echo "  assign_glob  : ${ASSIGN_GLOB}"
echo "  outdir       : ${OUTDIR}"
echo "  alpha        : ${alpha}"
echo "  policies     : ambiguous=${ambiguous_policy} doublet=${doublet_policy} indist=${indist_policy}"
echo "============================================================"

# -----------------------------------------------------------------------------
# Pre-flight checks. Use `find` (not `ls glob`) for E2BIG safety on dirs that
# may grow large (here they are small, but the convention is uniform).
# -----------------------------------------------------------------------------
[[ -f "${CELLS_CALLS}" ]] || { echo "ERROR: missing cells_calls: ${CELLS_CALLS}" >&2; exit 1; }
[[ -d "${CHUNKS_DIR}" ]]  || { echo "ERROR: missing chunks dir: ${CHUNKS_DIR}" >&2; exit 1; }

N_FILTERED=$(find "${CHUNKS_DIR}" -maxdepth 1 -name "${DATASET}_chunk*_filtered.tsv.gz" | wc -l)
if [[ ${N_FILTERED} -eq 0 ]]; then
    echo "ERROR: no filtered chunks matching ${DATASET}_chunk*_filtered.tsv.gz in ${CHUNKS_DIR}" >&2
    exit 1
fi
echo "  found ${N_FILTERED} filtered chunk files"

mkdir -p "${OUTDIR}"

# -----------------------------------------------------------------------------
# Run decontam (without design)
# -----------------------------------------------------------------------------
ambientmapper decontam \
  --cells-calls "${CELLS_CALLS}" \
  --out-dir "${OUTDIR}" \
  --threads "${threads}" \
  --assign-glob "${ASSIGN_GLOB}" \
  --ambiguous-policy "${ambiguous_policy}" \
  --doublet-policy "${doublet_policy}" \
  --indist-policy "${indist_policy}" \
  --decontam-alpha "${alpha}" \
  --chunksize "${chunksize}" \
  --min-reads-post-clean "${min_reads_post_clean}" \
  --min-allowed-frac-post-clean "${min_allowed_frac_post_clean}" \
  --safe-keep-delta-as "${safe_keep_delta_as}"

echo ""
echo "============================================================"
echo "[$(date)] DONE: ${TRACK}/${DATASET} (no-design decontam)"
echo "  outputs in ${OUTDIR}"
echo "============================================================"
