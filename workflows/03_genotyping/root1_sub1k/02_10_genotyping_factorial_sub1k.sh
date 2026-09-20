#!/usr/bin/env bash
#SBATCH --job-name=sub1k_fact
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --array=0-167
#SBATCH --output=_logs/02_10_sub1k_fact_%A_%a.log
#
# 02_10_genotyping_factorial_sub1k.sh — Sub1k factorial genotyping runs
#
# Runs all (panel × alpha × config) combinations of the Phase 1 C0..C4 factorial
# (42 single-knob configurations, listed below and mirrored in _genotyping_configs.sh).
# Output feeds eval_phase_factorial.R (phase 3 consistency check) and Fig. S6.
#
# Array layout: 2 panels × 2 alphas × 42 configs = 168 tasks
#   task_id = panel_idx * (N_ALPHA * N_CONFIG) + alpha_idx * N_CONFIG + config_idx
#
# Prerequisites:
#   1. 02_09_pipeline_sub1k.sh has completed for both panels (8 chunks dirs)
#   2. ambientmapper installed with topk-reclass + winner-discount support
#
# Output naming:
#   Root1_rep1/sub1k_<panel>/genotyping_runs/factorial_2026-04-07/<config>_alpha<tag>/
#
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

FACTORIAL_TAG=factorial_2026-04-07

# Pin ambientmapper version for reproducibility
# AMBIENTMAPPER_REPO: optional local checkout of the ambientmapper tool, used only to log its commit
AM_VERSION=$(pip show ambientmapper 2>/dev/null | awk '/^Version/ {print $2}')
AM_COMMIT=$(git -C "${AMBIENTMAPPER_REPO:-.}" \
    rev-parse --short HEAD 2>/dev/null || echo "unknown")
echo "ambientmapper: v${AM_VERSION} @ ${AM_COMMIT}"

# -----------------------------------------------------------------------------
# Task → (panel, alpha, config) mapping
# -----------------------------------------------------------------------------

PANELS=(A B)
ALPHA_TAGS=(alpha005 alpha1e6)

# 42 configurations — same names as the case statement in _genotyping_configs.sh
CONFIGS=(
    # C0 (1)
    C0
    # C1: rescue knockout 2^3 (7 new; C0 is the 8th cell)
    C1a_noreclass
    C1b_nowdisc
    C1c_nofriend
    C1d_noreclass_nowdisc
    C1e_noreclass_nofriend
    C1f_nowdisc_nofriend
    C1g_naked
    # C2: xmap layered + v2-fixeta-ma08 disentanglement (5)
    C2a_xmap
    C2b_xmap_noreclass
    C2c_xmap_eta0
    C2d_xmap_ma08
    C2e_xmap_eta0_ma08
    # C3: read filter ablations (5)
    C3a_mq10
    C3b_mq50
    C3c_xa_unlim
    C3d_mq10_xa_unlim
    C3e_mq50_xa_unlim
    # C4a: BIC singlet thresholds (10)
    C4a_mass04
    C4a_mass05
    C4a_mass07
    C4a_ratio12
    C4a_ratio18
    C4a_ratio20
    C4a_ma03
    C4a_ma08
    C4a_bic3
    C4a_bic10
    # C4b: scoring weights (6)
    C4b_was05
    C4b_was15
    C4b_wmq05
    C4b_wmq20
    C4b_wnm05
    C4b_wnm20
    # C4c: eta iters (2)
    C4c_eta0
    C4c_eta5
    # C4d: w_ambiguous × wdisc structural cross (6)
    C4d_wamb005_wdon
    C4d_wamb005_wdoff
    C4d_wamb02_wdon
    C4d_wamb02_wdoff
    C4d_wamb05_wdon
    C4d_wamb05_wdoff
)

N_PANELS=${#PANELS[@]}
N_ALPHAS=${#ALPHA_TAGS[@]}
N_CONFIGS=${#CONFIGS[@]}
N_TOTAL=$(( N_PANELS * N_ALPHAS * N_CONFIGS ))
if [[ ${N_TOTAL} -ne 168 ]]; then
    echo "ERROR: config count mismatch: expected 168, got ${N_TOTAL}" >&2
    echo "  (${N_PANELS} panels × ${N_ALPHAS} alphas × ${N_CONFIGS} configs)" >&2
    exit 1
fi

TASK_ID=${SLURM_ARRAY_TASK_ID:-0}
if [[ ${TASK_ID} -ge ${N_TOTAL} ]]; then
    echo "ERROR: task id ${TASK_ID} out of range [0, ${N_TOTAL})" >&2
    exit 1
fi

PANEL_IDX=$(( TASK_ID / (N_ALPHAS * N_CONFIGS) ))
ALPHA_IDX=$(( (TASK_ID / N_CONFIGS) % N_ALPHAS ))
CONFIG_IDX=$(( TASK_ID % N_CONFIGS ))

PANEL=${PANELS[$PANEL_IDX]}
ALPHA_TAG=${ALPHA_TAGS[$ALPHA_IDX]}
CONFIG=${CONFIGS[$CONFIG_IDX]}

SAMPLE=sub1k_${PANEL}
WORKDIR=${PROJECT_ROOT}/5_AmbientDetection/Root1_rep1
CONFIG_JSON=configs/Root1_rep1_sub1k_${PANEL}.ambientmapper.json
SAMPLE_DIR=${WORKDIR}/${SAMPLE}

# -----------------------------------------------------------------------------
# Friend-mode routing
# -----------------------------------------------------------------------------
#
# C1c/C1e/C1f/C1g knock out friend rescue. All other configs use the friendwith
# chunks dir. The knockout is the rescued -> ambiguous relabel performed by
# 02_09_pipeline_sub1k.sh (ambientmapper has no assign-step flag for it).
case "${CONFIG}" in
    C1c_nofriend|C1e_noreclass_nofriend|C1f_nowdisc_nofriend|C1g_naked)
        FRIEND=without
        ;;
    *)
        FRIEND=with
        ;;
esac

CHUNKS_DIR=${SAMPLE_DIR}/cell_map_ref_chunks_${ALPHA_TAG}_friend${FRIEND}
ASSIGN_GLOB="${CHUNKS_DIR}/*_filtered.tsv.gz"
OUTDIR=${SAMPLE_DIR}/genotyping_runs/${FACTORIAL_TAG}/${CONFIG}_${ALPHA_TAG}

echo "============================================================"
echo "[$(date)] sub1k factorial — task ${TASK_ID}"
echo "  panel       = ${PANEL}"
echo "  alpha       = ${ALPHA_TAG}"
echo "  config      = ${CONFIG}"
echo "  friend mode = ${FRIEND}"
echo "  chunks dir  = ${CHUNKS_DIR}"
echo "  outdir      = ${OUTDIR}"
echo "============================================================"

if [[ ! -d "${CHUNKS_DIR}" ]]; then
    echo "ERROR: chunks dir not found: ${CHUNKS_DIR}" >&2
    echo "       Run 02_09_pipeline_sub1k.sh first." >&2
    exit 1
fi
N_FILTERED=$(ls "${CHUNKS_DIR}"/*_filtered.tsv.gz 2>/dev/null | wc -l)
if [[ ${N_FILTERED} -eq 0 ]]; then
    echo "ERROR: no *_filtered.tsv.gz files in ${CHUNKS_DIR}" >&2
    exit 1
fi
echo "  found ${N_FILTERED} filtered files"

# -----------------------------------------------------------------------------
# C0 baseline: initialize all overridable flags to the armored-defaults values
# -----------------------------------------------------------------------------
#
# Every config below only overrides what it needs to; everything else stays at
# C0. This is what makes the factorial decomposable.

# Read-quality filters
MAPQ_MIN=20
XA_MAX_FLAG="--max-xa-genotyping 0"     # use empty string to disable the XA filter

# Topk + scoring weights
TOPK=2
W_CONFIDENT=1.0
W_AMBIG=0.1
W_AS=1.0
W_MAPQ=1.0
W_NM=1.0

# BIC thresholds
MASS_MIN=0.6
RATIO_MIN=1.5
MAX_ALPHA=0.5
BIC_MARGIN=6
DOUBLET_MINOR_MIN=0.20

# Eta
ETA_ITERS=2

# Empty-cell thresholds (frozen across the factorial)
EMPTY_BIC_MARGIN=6
EMPTY_TOP1_MAX=0.6
EMPTY_RATIO12_MAX=2
EMPTY_SEED_BIC_MIN=10
EMPTY_TAU_QUANTILE=0.95

# Rescue features
RECLASS_FLAG="--topk-reclass"
WDISC_FLAG="--winner-discount"
WDISC_MODE_FLAG="--winner-discount-mode winner_ratio"
XMAP_FLAG="--no-xmap"

# -----------------------------------------------------------------------------
# Config-specific overrides
# -----------------------------------------------------------------------------
case "${CONFIG}" in
    C0)
        : # pure baseline — no overrides
        ;;

    # ----- C1: rescue knockout factorial ----------------------------------
    C1a_noreclass)
        RECLASS_FLAG="--no-topk-reclass"
        ;;
    C1b_nowdisc)
        WDISC_FLAG="--no-winner-discount"
        ;;
    C1c_nofriend)
        : # friend rescue is knocked out via the friendwithout chunks dir
        ;;
    C1d_noreclass_nowdisc)
        RECLASS_FLAG="--no-topk-reclass"
        WDISC_FLAG="--no-winner-discount"
        ;;
    C1e_noreclass_nofriend)
        RECLASS_FLAG="--no-topk-reclass"
        ;;
    C1f_nowdisc_nofriend)
        WDISC_FLAG="--no-winner-discount"
        ;;
    C1g_naked)
        RECLASS_FLAG="--no-topk-reclass"
        WDISC_FLAG="--no-winner-discount"
        ;;

    # ----- C2: xmap layered + fixeta-ma08 disentanglement -----------------
    C2a_xmap)
        XMAP_FLAG=""   # omit --no-xmap → xmap enabled
        ;;
    C2b_xmap_noreclass)
        XMAP_FLAG=""
        RECLASS_FLAG="--no-topk-reclass"
        ;;
    C2c_xmap_eta0)
        XMAP_FLAG=""
        ETA_ITERS=0
        ;;
    C2d_xmap_ma08)
        XMAP_FLAG=""
        MAX_ALPHA=0.8
        ;;
    C2e_xmap_eta0_ma08)
        XMAP_FLAG=""
        ETA_ITERS=0
        MAX_ALPHA=0.8
        ;;

    # ----- C3: read filter ablations --------------------------------------
    C3a_mq10)
        MAPQ_MIN=10
        ;;
    C3b_mq50)
        MAPQ_MIN=50
        ;;
    C3c_xa_unlim)
        XA_MAX_FLAG=""   # omit the XA filter entirely
        ;;
    C3d_mq10_xa_unlim)
        MAPQ_MIN=10
        XA_MAX_FLAG=""
        ;;
    C3e_mq50_xa_unlim)
        MAPQ_MIN=50
        XA_MAX_FLAG=""
        ;;

    # ----- C4a: BIC singlet thresholds ------------------------------------
    C4a_mass04)  MASS_MIN=0.4 ;;
    C4a_mass05)  MASS_MIN=0.5 ;;
    C4a_mass07)  MASS_MIN=0.7 ;;
    C4a_ratio12) RATIO_MIN=1.2 ;;
    C4a_ratio18) RATIO_MIN=1.8 ;;
    C4a_ratio20) RATIO_MIN=2.0 ;;
    C4a_ma03)    MAX_ALPHA=0.3 ;;
    C4a_ma08)    MAX_ALPHA=0.8 ;;
    C4a_bic3)    BIC_MARGIN=3 ;;
    C4a_bic10)   BIC_MARGIN=10 ;;

    # ----- C4b: scoring weights -------------------------------------------
    C4b_was05) W_AS=0.5 ;;
    C4b_was15) W_AS=1.5 ;;
    C4b_wmq05) W_MAPQ=0.5 ;;
    C4b_wmq20) W_MAPQ=2.0 ;;
    C4b_wnm05) W_NM=0.5 ;;
    C4b_wnm20) W_NM=2.0 ;;

    # ----- C4c: eta iters -------------------------------------------------
    C4c_eta0) ETA_ITERS=0 ;;
    C4c_eta5) ETA_ITERS=5 ;;

    # ----- C4d: w_ambiguous × wdisc structural cross ----------------------
    C4d_wamb005_wdon)  W_AMBIG=0.05 ;;
    C4d_wamb005_wdoff) W_AMBIG=0.05; WDISC_FLAG="--no-winner-discount" ;;
    C4d_wamb02_wdon)   W_AMBIG=0.2 ;;
    C4d_wamb02_wdoff)  W_AMBIG=0.2;  WDISC_FLAG="--no-winner-discount" ;;
    C4d_wamb05_wdon)   W_AMBIG=0.5 ;;
    C4d_wamb05_wdoff)  W_AMBIG=0.5;  WDISC_FLAG="--no-winner-discount" ;;

    *)
        echo "ERROR: unknown config '${CONFIG}'" >&2
        exit 1
        ;;
esac

# -----------------------------------------------------------------------------
# Sanity checks
# -----------------------------------------------------------------------------
# Use Typer/Click introspection on the genotyping subcommand's parameter list.
# The previous `--help | grep` approach was flaky on compute nodes because
# Click's rich-help renderer wraps flag names based on the detected terminal
# width — narrow widths truncate `--winner-discount` to `--winner-discoun…`
# and grep misses it. ambientmapper uses Typer, so we resolve the underlying
# Click command via `typer.main.get_command(app).commands['genotyping']`.
_GENO_FLAGS=$(python -c "import ambientmapper.cli as _c, typer; g = typer.main.get_command(_c.app).commands['genotyping']; names = set(); [names.update(p.opts) for p in g.params]; print(' '.join(sorted(names)))" 2>/dev/null)
if [[ -z "${_GENO_FLAGS}" ]]; then
    echo "ERROR: failed to introspect ambientmapper.cli.app — is the package installed?" >&2
    exit 1
fi
if [[ ! " ${_GENO_FLAGS} " == *" --winner-discount "* ]]; then
    echo "ERROR: --winner-discount flag not found in installed ambientmapper" >&2
    exit 1
fi
if [[ ! " ${_GENO_FLAGS} " == *" --topk-reclass "* ]]; then
    echo "ERROR: --topk-reclass flag not found in installed ambientmapper" >&2
    exit 1
fi
unset _GENO_FLAGS

mkdir -p "${OUTDIR}"

# -----------------------------------------------------------------------------
# Run ambientmapper genotyping
# -----------------------------------------------------------------------------
# Flags are grouped in the order of the design doc C0 section.

echo ""
echo "[$(date)] Running ambientmapper genotyping ..."
echo "  effective flags:"
echo "    --min-mapq-genotyping ${MAPQ_MIN}"
echo "    ${XA_MAX_FLAG:-(no XA filter)}"
echo "    --topk-genomes ${TOPK}"
echo "    --w-confident ${W_CONFIDENT} --w-ambiguous ${W_AMBIG}"
echo "    --w-as ${W_AS} --w-mapq ${W_MAPQ} --w-nm ${W_NM}"
echo "    --single-mass-min ${MASS_MIN} --ratio-top1-top2-min ${RATIO_MIN}"
echo "    --max-alpha ${MAX_ALPHA} --bic-margin ${BIC_MARGIN}"
echo "    --eta-iters ${ETA_ITERS}"
echo "    ${RECLASS_FLAG}  ${WDISC_FLAG}  ${WDISC_MODE_FLAG}  ${XMAP_FLAG}"
echo ""

# shellcheck disable=SC2086
ambientmapper genotyping \
    --config "${CONFIG_JSON}" \
    --assign "${ASSIGN_GLOB}" \
    --outdir "${OUTDIR}" \
    --no-resume \
    --threads 8 \
    --pass1-workers 8 \
    --no-winner-only \
    --beta 10 \
    --min-reads 5 \
    --chunk-rows 100000 \
    --w-confident ${W_CONFIDENT} \
    --w-ambiguous ${W_AMBIG} \
    --w-as ${W_AS} \
    --w-mapq ${W_MAPQ} \
    --w-nm ${W_NM} \
    --topk-genomes ${TOPK} \
    --doublet-minor-min ${DOUBLET_MINOR_MIN} \
    --single-mass-min ${MASS_MIN} \
    --ratio-top1-top2-min ${RATIO_MIN} \
    --max-alpha ${MAX_ALPHA} \
    --bic-margin ${BIC_MARGIN} \
    --eta-iters ${ETA_ITERS} \
    --empty-bic-margin ${EMPTY_BIC_MARGIN} \
    --empty-top1-max ${EMPTY_TOP1_MAX} \
    --empty-ratio12-max ${EMPTY_RATIO12_MAX} \
    --empty-seed-bic-min ${EMPTY_SEED_BIC_MIN} \
    --empty-tau-quantile ${EMPTY_TAU_QUANTILE} \
    --min-mapq-genotyping ${MAPQ_MIN} \
    ${XA_MAX_FLAG} \
    ${RECLASS_FLAG} \
    ${WDISC_FLAG} \
    ${WDISC_MODE_FLAG} \
    ${XMAP_FLAG}

echo ""
echo "============================================================"
echo "[$(date)] Done: ${CONFIG} / panel=${PANEL} / ${ALPHA_TAG}"
echo "  outputs in ${OUTDIR}"
echo "============================================================"
