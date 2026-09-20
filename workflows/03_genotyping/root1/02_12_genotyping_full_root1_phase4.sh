#!/usr/bin/env bash
#SBATCH --job-name=root1_phase4
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=_logs/02_12_root1_phase4_%A_%a.log
#SBATCH --array=0-4
#
# 02_12_genotyping_full_root1_phase4.sh — Phase 4 full Root1 genotyping
#
# Runs 5 configurations on the FULL Root1_rep1 dataset (~27K barcodes ≥250
# reads, ~15.9K ≥500 reads). The outputs become the manuscript-grade numbers
# and are used for the sub1k-vs-full consistency check in the Phase 4 eval.
#
# Selected from Phase 3 sub1k F1_relaxed ranking (mean across panels × b2-b6):
#   C0                       — baseline anchor (F1=0.523)
#   C3b_mq50                 — mq50 only (isolates filter contribution; F1=0.590)
#   C4a_bic3                 — bic_margin=3 only (isolates BIC effect; F1=0.573)
#   S04_nofr_mq50_xmap_xaUL  — runner-up, prec=1.000 alternative (F1=0.674)
#   S10_full_bic3            — Phase 3 winner, 89/996 strict singlets (F1=0.681)
#
# Total tasks: 5 (one per config).
#
# Prerequisites:
#   - 02_12a_pipeline_full_root1_alpha005.sh has completed
#       → Root1_rep1/cell_map_ref_chunks_alpha005_friendwith/ exists
#       → Root1_rep1/cell_map_ref_chunks_alpha005_friendwithout/ exists
#   - configs/Root1_rep1.ambientmapper.json exists (full 26-genome config)
#   - Phase 3 eval has run and the WINNER_CONFIG name is known
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
# Phase 4 design — 5 configs
# -----------------------------------------------------------------------------
CONFIGS=(
    C0                         # baseline anchor (F1=0.523 on sub1k)
    C3b_mq50                   # mq50 only ablation (F1=0.590)
    C4a_bic3                   # bic_margin=3 only ablation (F1=0.573)
    S04_nofr_mq50_xmap_xaUL    # xmap-based runner-up, prec=1.000 (F1=0.674)
    S10_full_bic3              # Phase 3 winner, 8.9% strict singlets (F1=0.681)
)

N_CONFIGS=${#CONFIGS[@]}   # 5
TASK_ID=${SLURM_ARRAY_TASK_ID:-0}
if [[ ${TASK_ID} -ge ${N_CONFIGS} ]]; then
    echo "ERROR: TASK_ID=${TASK_ID} out of range (max ${N_CONFIGS})" >&2
    exit 1
fi

CONFIG="${CONFIGS[$TASK_ID]}"

# -----------------------------------------------------------------------------
# Path construction (full Root1, NOT sub1k)
# -----------------------------------------------------------------------------
SAMPLE=Root1_rep1
WORKDIR=${PROJECT_ROOT}/5_AmbientDetection/Root1_rep1
CONFIG_JSON=configs/Root1_rep1.ambientmapper.json

# Friend rescue routing.
#   - friend ON  : the canonical α=0.05 output written by 02_12a directly
#                  into cell_map_ref_chunks/ (no _alpha005_friendwith sibling
#                  exists — the original IS the friend-with version).
#   - friend OFF : the post-hoc rescued→ambiguous relabel sibling, written
#                  by 02_12a step 2 / repaired by 02_12b.
FRIEND=$(friend_mode_for_config "${CONFIG}")
if [[ "${FRIEND}" == "with" ]]; then
    CHUNKS_DIR=${WORKDIR}/cell_map_ref_chunks
else
    CHUNKS_DIR=${WORKDIR}/cell_map_ref_chunks_alpha005_friendwithout
fi
ASSIGN_GLOB="${CHUNKS_DIR}/*_filtered.tsv.gz"

FACTORIAL_TAG="factorial_phase4_$(date +%Y-%m-%d)"
OUTDIR=${WORKDIR}/genotyping_runs/${FACTORIAL_TAG}/${CONFIG}

echo "============================================================"
echo "[$(date)] Phase 4 full Root1 genotyping — task ${TASK_ID}"
echo "  sample        = ${SAMPLE} (FULL data, not sub1k)"
echo "  config        = ${CONFIG}"
echo "  friend mode   = ${FRIEND}"
echo "  chunks dir    = ${CHUNKS_DIR}"
echo "  outdir        = ${OUTDIR}"
echo "============================================================"

if [[ ! -f "${CONFIG_JSON}" ]]; then
    echo "ERROR: config JSON not found: ${CONFIG_JSON}" >&2
    exit 1
fi
if [[ ! -d "${CHUNKS_DIR}" ]]; then
    echo "ERROR: chunks dir not found: ${CHUNKS_DIR}" >&2
    echo "       Did you run 02_12a_pipeline_full_root1_alpha005.sh first?" >&2
    exit 1
fi
N_FILTERED=$(ls "${CHUNKS_DIR}"/*_filtered.tsv.gz 2>/dev/null | wc -l)
if [[ ${N_FILTERED} -eq 0 ]]; then
    echo "ERROR: no *_filtered.tsv.gz files in ${CHUNKS_DIR}" >&2
    exit 1
fi
echo "  found ${N_FILTERED} filtered files"

# Sanity-check the chunks dir: header has assigned_class + at least one winner row.
# Note: filtered.tsv.gz is one row per (read × genome-hit), so per-row winner% is
# structurally ~2% even on healthy α=0.05 data (per-read winner% would be ~27%).
# This dilution makes per-row thresholds useless for α=0.05 vs α=1e-6 detection
# (legacy ~24% per-read vs current ~27% per-read — signal too weak), so we just
# verify the dir is non-empty and well-formed. Alpha-bundle correctness is
# enforced by 02_12a's post-conditions and the chunks-dir naming convention.
echo "  sanity-checking ${CHUNKS_DIR} (10-chunk sample)..."
python - <<PYEOF
import glob
import gzip
import os
import sys
from collections import Counter

src = "${CHUNKS_DIR}"
files = sorted(glob.glob(os.path.join(src, "*_filtered.tsv.gz")))
if not files:
    sys.exit("ERROR: no filtered files in " + src)

n_sample = min(10, len(files))
step = max(1, len(files) // n_sample)
sampled = files[::step][:n_sample]

c = Counter()
for f in sampled:
    with gzip.open(f, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if "assigned_class" not in header:
            sys.exit(f"ERROR: 'assigned_class' column missing in {f}")
        i = header.index("assigned_class")
        for line in fh:
            c[line.rstrip("\n").split("\t")[i]] += 1
total = sum(c.values())
print(f"  sampled {len(sampled)} chunks: winner={c.get('winner',0):,} "
      f"rescued={c.get('rescued',0):,} ambiguous={c.get('ambiguous',0):,} total={total:,}")
if total == 0:
    sys.exit("ERROR: sampled chunks contain zero rows")
if c.get("winner", 0) == 0:
    sys.exit("ERROR: zero winner rows in sampled chunks — assign step likely didn't run")
PYEOF

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
echo "[$(date)] Done: ${SAMPLE} / ${CONFIG}"
echo "  outputs in ${OUTDIR}"
echo "============================================================"
