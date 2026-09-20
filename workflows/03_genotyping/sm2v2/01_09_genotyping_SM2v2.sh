#!/usr/bin/env bash
#SBATCH --job-name=sm2v2_geno
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --output=_logs/01_09_sm2v2_geno_%A.log
#
# 01_09_genotyping_SM2v2.sh — SM2v2 step 6/7: genotyping (C0 baseline).
#
# Uses the shared _genotyping_configs.sh helper to apply the C0 sub1k-validated
# baseline (mq20/xa=0/topk=2 + topk-reclass + winner-discount + friend rescue ON,
# xmap OFF). Friend rescue is ON by default (no _friendwithout sibling needed
# for SM2v2 since we are only running one config, not the factorial).
#
# Output:
#   SM2v2/genotyping_runs/SM2v2_C0/SM2v2_cells_calls.tsv.gz
#   SM2v2/genotyping_runs/SM2v2_C0/SM2v2_C_all.tsv.gz
#   SM2v2/genotyping_runs/SM2v2_C0/SM2v2_eta.tsv.gz
#   ...
#
# Also writes SM2v2/final/ symlink for downstream decontam discovery.
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

source "${REPO_ROOT}/workflows/03_genotyping/_genotyping_configs.sh"

ambientmapper --version 2>/dev/null || echo "(ambientmapper --version not available)"
check_ambientmapper_cli_features

SAMPLE="SM2v2"
WORKDIR="${PROJECT_ROOT}/5_AmbientDetection/SM2v2"
CONFIG_JSON="configs/SM2v2.ambientmapper.json"
CONFIG="C0"

CHUNKS_DIR="${WORKDIR}/cell_map_ref_chunks"
ASSIGN_GLOB="${CHUNKS_DIR}/*_filtered.tsv.gz"
OUTDIR="${WORKDIR}/genotyping_runs/${SAMPLE}_${CONFIG}"

echo "============================================================"
echo "[$(date)] SM2v2 step 6/7: genotyping (${CONFIG})"
echo "  sample        = ${SAMPLE}"
echo "  config        = ${CONFIG}"
echo "  chunks dir    = ${CHUNKS_DIR}"
echo "  outdir        = ${OUTDIR}"
echo "============================================================"

[[ -f "${CONFIG_JSON}" ]] || { echo "ERROR: config JSON not found: ${CONFIG_JSON}" >&2; exit 1; }
[[ -d "${CHUNKS_DIR}" ]] || { echo "ERROR: chunks dir not found: ${CHUNKS_DIR}" >&2; exit 1; }
N_FILTERED=$(ls "${CHUNKS_DIR}"/*_filtered.tsv.gz 2>/dev/null | wc -l)
[[ ${N_FILTERED} -gt 0 ]] || { echo "ERROR: no *_filtered.tsv.gz in ${CHUNKS_DIR}" >&2; exit 1; }
echo "  found ${N_FILTERED} filtered files"

mkdir -p "${OUTDIR}"

init_c0_defaults
apply_config_overrides "${CONFIG}"

echo ""
print_effective_flags
echo ""

run_ambientmapper_genotyping "${CONFIG_JSON}" "${ASSIGN_GLOB}" "${OUTDIR}"

# Promote cells_calls into final/ for downstream decontam (mirrors existing SM2 layout)
mkdir -p "${WORKDIR}/final"
ln -sf "../genotyping_runs/${SAMPLE}_${CONFIG}/${SAMPLE}_cells_calls.tsv.gz" \
       "${WORKDIR}/final/${SAMPLE}_cells_calls.tsv.gz"

echo ""
echo "============================================================"
echo "[$(date)] genotyping DONE — outputs in ${OUTDIR}"
echo "[$(date)] symlink: ${WORKDIR}/final/${SAMPLE}_cells_calls.tsv.gz → ../genotyping_runs/${SAMPLE}_${CONFIG}/"
echo "============================================================"
