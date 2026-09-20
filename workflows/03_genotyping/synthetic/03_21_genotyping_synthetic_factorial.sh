#!/usr/bin/env bash
#SBATCH --job-name=syn_phase2
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=_logs/03_21_syn_phase2_%A_%a.log
#SBATCH --array=0-299%60
#
# 03_21_genotyping_synthetic_factorial.sh — Phase 2 synthetic factorial
#
# Runs 10 single-knob ambientmapper genotyping configs across 15 datasets ×
# 2 tracks (Track B = synthetic/, Track B-disc = synthetic_disc/) using
# C0-aligned flag values (--w-as 1.0 etc.). Output feeds eval_phase_factorial.R
# phase2 -> Fig. S5, and Fig. 4B reads the C0 cells_calls directly.
#
# Phase 2 question: which Phase 1 single-knob findings generalize to data
# with real contamination? Especially: does friend rescue help on contaminated
# data even though it hurts on the 26-genome Root1 reference panel?
#
# Total tasks: 10 configs × 15 datasets × 2 tracks = 300.
#
# Array index decode:
#   TRACK_IDX  = TASK / (N_DATASETS * N_CONFIGS)   [0..1]
#   DS_IDX     = (TASK / N_CONFIGS) % N_DATASETS   [0..14]
#   CONFIG_IDX = TASK % N_CONFIGS                  [0..9]
#
# Prerequisites:
#   - 03_21a_synthetic_friend_relabel.sh has completed (the
#     `cell_map_ref_chunks_friendwithout/` siblings exist for both tracks)
#   - synthetic[_disc]/<dataset>/config.json exists for every dataset
#   - synthetic[_disc]/<dataset>/cell_map_ref_chunks/*_filtered.tsv.gz exist
#   - ambientmapper installed with --winner-discount and --topk-reclass flags
#
# =============================================================================

set -euo pipefail

REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"   # this repository (sourced helpers, sibling scripts)
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
# Source the shared genotyping helper. Provides:
#   init_c0_defaults, apply_config_overrides, friend_mode_for_config,
#   print_effective_flags, check_ambientmapper_cli_features,
#   run_ambientmapper_genotyping
# -----------------------------------------------------------------------------
source "${REPO_ROOT}/workflows/03_genotyping/_genotyping_configs.sh"

# Verify ambientmapper version + key flags
ambientmapper --version 2>/dev/null || echo "(ambientmapper --version not available)"
python -c "import ambientmapper; print('ambientmapper:', getattr(ambientmapper, '__version__', '?'))" 2>/dev/null || true
check_ambientmapper_cli_features

# -----------------------------------------------------------------------------
# Phase 2 design — 10 single-knob configs
# -----------------------------------------------------------------------------
CONFIGS=(
    C0                       # baseline
    C1c_nofriend             # does friend rescue help with real contamination?
    C2a_xmap                 # does xmap help with real contamination?
    C2c_xmap_eta0            # xmap × eta_iters=0 cross-tier interaction
    C3b_mq50                 # does MAPQ≥50 generalize?
    C4a_bic3                 # bic_margin sensitivity on synthetic
    C4c_eta0                 # eta-trap test (Track B-disc has ground-truth doublets)
    C4c_eta5                 # does MORE eta iters help with real contamination?
    C4d_wamb005_wdon         # boundary point
    C4d_wamb05_wdon          # does w_amb=0.5 generalize?
)

DATASETS=(
    alpha_000
    alpha_002_Il14H alpha_005_Il14H alpha_010_Il14H alpha_020_Il14H
    alpha_030_Il14H alpha_040_Il14H alpha_050_Il14H
    alpha_002_Ki11  alpha_005_Ki11  alpha_010_Ki11  alpha_020_Ki11
    alpha_030_Ki11  alpha_040_Ki11  alpha_050_Ki11
)

TRACKS=(synthetic synthetic_disc)

N_CONFIGS=${#CONFIGS[@]}     # 10
N_DATASETS=${#DATASETS[@]}   # 15
N_TRACKS=${#TRACKS[@]}       # 2
N_TOTAL=$(( N_CONFIGS * N_DATASETS * N_TRACKS ))   # 300

if [[ ${N_TOTAL} -ne 300 ]]; then
    echo "ERROR: expected 300 total tasks, got ${N_TOTAL} (configs=${N_CONFIGS}, ds=${N_DATASETS}, tracks=${N_TRACKS})" >&2
    exit 1
fi

TASK_ID=${SLURM_ARRAY_TASK_ID:-0}
if [[ ${TASK_ID} -ge ${N_TOTAL} ]]; then
    echo "ERROR: TASK_ID=${TASK_ID} out of range (max ${N_TOTAL})" >&2
    exit 1
fi

TRACK_IDX=$(( TASK_ID / (N_DATASETS * N_CONFIGS) ))
DS_IDX=$(( (TASK_ID / N_CONFIGS) % N_DATASETS ))
CONFIG_IDX=$(( TASK_ID % N_CONFIGS ))

TRACK="${TRACKS[$TRACK_IDX]}"
DATASET="${DATASETS[$DS_IDX]}"
CONFIG="${CONFIGS[$CONFIG_IDX]}"

# -----------------------------------------------------------------------------
# Path construction
# -----------------------------------------------------------------------------
WORKDIR="${PROJECT_ROOT}/5_AmbientDetection"
DS_DIR="${WORKDIR}/${TRACK}/${DATASET}"
CONFIG_JSON="${DS_DIR}/config.json"

# Friend rescue routing — synthetic uses sibling dirs, NOT alpha-suffixed
# (sub1k convention: cell_map_ref_chunks_alpha005_friend{with,without})
FRIEND=$(friend_mode_for_config "${CONFIG}")
if [[ "${FRIEND}" == "with" ]]; then
    CHUNKS_DIR="${DS_DIR}/cell_map_ref_chunks"
else
    CHUNKS_DIR="${DS_DIR}/cell_map_ref_chunks_friendwithout"
fi
ASSIGN_GLOB="${CHUNKS_DIR}/*_filtered.tsv.gz"

# Output dir tag includes the date so re-runs don't collide
FACTORIAL_TAG="factorial_phase2_$(date +%Y-%m-%d)"
OUTDIR="${DS_DIR}/genotyping_runs/${FACTORIAL_TAG}/${CONFIG}"

echo "============================================================"
echo "[$(date)] Phase 2 synthetic factorial — task ${TASK_ID}"
echo "  track       = ${TRACK}"
echo "  dataset     = ${DATASET}"
echo "  config      = ${CONFIG}"
echo "  friend mode = ${FRIEND}"
echo "  chunks dir  = ${CHUNKS_DIR}"
echo "  outdir      = ${OUTDIR}"
echo "============================================================"

if [[ ! -f "${CONFIG_JSON}" ]]; then
    echo "ERROR: config JSON not found: ${CONFIG_JSON}" >&2
    exit 1
fi
if [[ ! -d "${CHUNKS_DIR}" ]]; then
    echo "ERROR: chunks dir not found: ${CHUNKS_DIR}" >&2
    if [[ "${FRIEND}" == "without" ]]; then
        echo "       Did you run 03_21a_synthetic_friend_relabel.sh first?" >&2
    fi
    exit 1
fi
N_FILTERED=$(ls "${CHUNKS_DIR}"/*_filtered.tsv.gz 2>/dev/null | wc -l)
if [[ ${N_FILTERED} -eq 0 ]]; then
    echo "ERROR: no *_filtered.tsv.gz files in ${CHUNKS_DIR}" >&2
    exit 1
fi
echo "  found ${N_FILTERED} filtered files"

mkdir -p "${OUTDIR}"

# -----------------------------------------------------------------------------
# Apply C0 defaults + the per-config overrides, then run genotyping
# -----------------------------------------------------------------------------
init_c0_defaults
apply_config_overrides "${CONFIG}"

echo ""
print_effective_flags
echo ""

run_ambientmapper_genotyping "${CONFIG_JSON}" "${ASSIGN_GLOB}" "${OUTDIR}"

echo ""
echo "============================================================"
echo "[$(date)] Done: ${TRACK}/${DATASET}/${CONFIG}"
echo "  outputs in ${OUTDIR}"
echo "============================================================"
