#!/usr/bin/env bash
#SBATCH --job-name=sub1k_phase3
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=_logs/02_11_sub1k_phase3_%A_%a.log
#SBATCH --array=0-37
#
# 02_11_genotyping_stacked_phase3_sub1k.sh — Phase 3 stacked factorial
#
# Tests 13 stacked configurations (S01..S13) plus 6 reference configs on
# the existing Root1 sub1k panels A and B at α=0.05 only. The α=1e-6 axis
# is dropped after Phase 1 determined that α=0.05 is the right assign bundle.
#
# Stacked configs combine winning single-knob changes from Phase 1.
# Reference configs (C0, C1c, C1g, C3b, C4a, C4d_wamb005) provide:
#   - C0: baseline (same as Phase 1 C0, for in-figure comparison)
#   - C1c/C1g: lower bounds (friend OFF / all rescue OFF)
#   - C3b/C4a/C4d_wamb005: Phase 2 synthetic winners, cross-validated on Root1
#
# The S01..S13 stack table lives in _genotyping_configs.sh. Output feeds
# eval_phase_factorial.R phase3 -> Fig. S6.
#
# Total tasks: 19 configs × 2 panels × 1 alpha = 38.
#
# Array index decode (configs rotate fastest, then panels):
#   PANEL_IDX  = TASK / N_CONFIGS    [0..1]
#   CONFIG_IDX = TASK % N_CONFIGS    [0..12]
#
# Prerequisites:
#   - Sub1k panels A and B exist (Root1_rep1/sub1k_{A,B}/)
#   - cell_map_ref_chunks_alpha005_friend{with,without}/ exist for each panel
#   - configs/Root1_rep1_sub1k_{A,B}.ambientmapper.json exist
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
# Source the shared genotyping helper
# -----------------------------------------------------------------------------
source "${REPO_ROOT}/workflows/03_genotyping/_genotyping_configs.sh"

ambientmapper --version 2>/dev/null || echo "(ambientmapper --version not available)"
python -c "import ambientmapper; print('ambientmapper:', getattr(ambientmapper, '__version__', '?'))" 2>/dev/null || true
check_ambientmapper_cli_features

# -----------------------------------------------------------------------------
# Phase 3 design — 13 stacked configurations
# -----------------------------------------------------------------------------
CONFIGS=(
    # --- Reference configs (Phase 1 single-knob, for in-figure comparison) ---
    C0                         # baseline (same as Phase 1 C0)
    C1c_nofriend               # friend OFF only (Phase 2 showed opposite on synthetic)
    C1g_naked                  # all rescue OFF (lower bound)
    C3b_mq50                   # best single-knob generalizer (Phase 1 + Phase 2)
    C4a_bic3                   # bic_margin=3 (strong on Track B-disc)
    C4d_wamb005_wdon           # w_amb=0.05 (synthetic Track B winner — cross-validate)
    # --- Phase 3 stacked configs ---
    S01_nofr_mq50              # minimal stack: friend off + mq50
    S02_nofr_mq50_xaUL         # + XA unlimited
    S03_nofr_mq50_xmap         # + xmap (with eta0)
    S04_nofr_mq50_xmap_xaUL    # + xmap + XA unlim
    S05_nofr_mq50_wamb         # + w_amb=0.5 (no xmap)
    S06_nofr_mq50_wamb_xaUL    # + w_amb=0.5 + XA unlim
    S07_nofr_mq50_xmap_wamb    # + xmap + w_amb=0.5
    S08_full_xa0               # FULL STACK with XA=0
    S09_full_xaUL              # FULL STACK with XA unlim
    S10_full_bic3              # FULL STACK + bic_margin=3 bonus
    S11_full_friendon          # sanity: does friend ON break the full stack?
    S12_full_wamb01            # sanity: w_amb=0.1 inside full stack
    S13_full_eta2              # sanity: eta_iters=2 inside full stack
)

PANELS=(A B)
ALPHA_TAG=alpha005   # Phase 3 uses α=0.05 only

N_CONFIGS=${#CONFIGS[@]}   # 19
N_PANELS=${#PANELS[@]}     # 2
N_TOTAL=$(( N_CONFIGS * N_PANELS ))   # 38

if [[ ${N_TOTAL} -ne 38 ]]; then
    echo "ERROR: expected 38 total tasks, got ${N_TOTAL}" >&2
    exit 1
fi

TASK_ID=${SLURM_ARRAY_TASK_ID:-0}
if [[ ${TASK_ID} -ge ${N_TOTAL} ]]; then
    echo "ERROR: TASK_ID=${TASK_ID} out of range (max ${N_TOTAL})" >&2
    exit 1
fi

PANEL_IDX=$(( TASK_ID / N_CONFIGS ))
CONFIG_IDX=$(( TASK_ID % N_CONFIGS ))

PANEL="${PANELS[$PANEL_IDX]}"
CONFIG="${CONFIGS[$CONFIG_IDX]}"

# -----------------------------------------------------------------------------
# Path construction (mirrors 02_10:137-160)
# -----------------------------------------------------------------------------
SAMPLE=sub1k_${PANEL}
WORKDIR=${PROJECT_ROOT}/5_AmbientDetection/Root1_rep1
CONFIG_JSON=configs/Root1_rep1_sub1k_${PANEL}.ambientmapper.json
SAMPLE_DIR=${WORKDIR}/${SAMPLE}

# Friend rescue routing — uses the existing sub1k _alpha005_friend{with,without}/
# chunks dirs created by 02_09_pipeline_sub1k.sh.
FRIEND=$(friend_mode_for_config "${CONFIG}")
CHUNKS_DIR=${SAMPLE_DIR}/cell_map_ref_chunks_${ALPHA_TAG}_friend${FRIEND}
ASSIGN_GLOB="${CHUNKS_DIR}/*_filtered.tsv.gz"

FACTORIAL_TAG="factorial_phase3_$(date +%Y-%m-%d)"
OUTDIR=${SAMPLE_DIR}/genotyping_runs/${FACTORIAL_TAG}/${CONFIG}

echo "============================================================"
echo "[$(date)] Phase 3 sub1k stacked factorial — task ${TASK_ID}"
echo "  panel       = ${PANEL}"
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
    echo "       Did you run 02_09_pipeline_sub1k.sh first?" >&2
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
echo "[$(date)] Done: panel=${PANEL} / ${CONFIG}"
echo "  outputs in ${OUTDIR}"
echo "============================================================"
