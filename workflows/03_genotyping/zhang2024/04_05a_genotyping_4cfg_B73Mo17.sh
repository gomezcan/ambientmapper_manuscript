#!/usr/bin/env bash
#SBATCH --job-name=4cfg_b73mo17
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=36:00:00
#SBATCH --output=_logs/04_05a_4cfg_b73mo17_%A_%a.log
#SBATCH --array=0-7
#
# 04_05a_genotyping_4cfg_B73Mo17.sh
#
# B73Mo17 half of the 4-configuration genotyping grid (the first combined array hit
# the 22 h walltime on rep1 at the Pass 2.75 merge; rep2 finished in about 20 h).
# multiGenotypes_rep1 is handled separately in 04_05b (needs 128 G mem + 48 h time).
# Only the C0 output is displayed (Fig. 4H to K read <sample>_cells_calls.tsv.gz).
#
# Output dir FACTORIAL_TAG="4cfg_2026-05-01" kept intact so downstream Fig 4
# scripts that reference that name continue to resolve.
#
# 8-task array: 2 datasets × 4 configs.
#   DATASET_IDX = TASK_ID / 4
#     0  B73Mo17_rep1  (2 genomes, 64,111 chunks)
#     1  B73Mo17_rep2  (2 genomes, 55,539 chunks)
#   CONFIG_IDX = TASK_ID % 4
#     0  C0           baseline (mq20, xmap OFF)
#     1  C2a_xmap     C0 + xmap ON
#     2  C3b_mq50     C0 + mq50
#     3  Cxmap_mq50   C0 + xmap ON + mq50
#
# Resource budget rationale:
#   - rep2 / C0 finished in ~20 h on 64 G in the first run.
#   - rep1 has 1.15× more chunks (64 K vs 55 K) → projects to ~23 h.
#   - 36 h gives 13 h headroom for slower Pass 2.75 merge with more shards,
#     cluster noise, or marginal config drift.
#   - 64 G held for all 4 rep2 configs in May; no OOMs in any B73Mo17 task.
#
# topk-genomes = 2 (C0 default). Friend rescue ON (canonical chunks dir).
#
# Prerequisites (unchanged from original 04_05):
#   - {sample}/cell_map_ref_chunks/*_filtered.tsv.gz exists
#   - configs/{B73Mo17_rep1,B73Mo17_rep2}.ambientmapper.json present
#   - ambientmapper installed with --winner-discount and --topk-reclass flags
#   - No stale 4cfg_2026-05-01/<config>/ output under the sample (ambientmapper
#     runs with --no-resume, so leftover sentinels are ignored, but stale shard
#     tmpfiles waste disk).
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

source "${REPO_ROOT}/workflows/03_genotyping/_genotyping_configs.sh"

ambientmapper --version 2>/dev/null || echo "(ambientmapper --version not available)"
python -c "import ambientmapper; print('ambientmapper:', getattr(ambientmapper, '__version__', '?'))" 2>/dev/null || true
check_ambientmapper_cli_features

# -----------------------------------------------------------------------------
# Task → (dataset, config) dispatch
# -----------------------------------------------------------------------------
DATASETS=(B73Mo17_rep1 B73Mo17_rep2)
CONFIGS=(C0 C2a_xmap C3b_mq50 Cxmap_mq50)

N_DATASETS=${#DATASETS[@]}    # 2
N_CONFIGS=${#CONFIGS[@]}      # 4
N_TOTAL=$(( N_DATASETS * N_CONFIGS ))   # 8

TASK_ID=${SLURM_ARRAY_TASK_ID:-0}
if [[ ${TASK_ID} -ge ${N_TOTAL} ]]; then
    echo "ERROR: TASK_ID=${TASK_ID} out of range (max ${N_TOTAL})" >&2
    exit 1
fi

DATASET_IDX=$(( TASK_ID / N_CONFIGS ))
CONFIG_IDX=$(( TASK_ID % N_CONFIGS ))

SAMPLE=${DATASETS[$DATASET_IDX]}
CONFIG=${CONFIGS[$CONFIG_IDX]}

WORKDIR=${PROJECT_ROOT}/5_AmbientDetection/${SAMPLE}
CONFIG_JSON=configs/${SAMPLE}.ambientmapper.json

CHUNKS_DIR=${WORKDIR}/cell_map_ref_chunks
ASSIGN_GLOB="${CHUNKS_DIR}/*_filtered.tsv.gz"

FACTORIAL_TAG="4cfg_2026-05-01"
OUTDIR=${WORKDIR}/genotyping_runs/${FACTORIAL_TAG}/${CONFIG}

echo "============================================================"
echo "[$(date)] 4cfg B73Mo17 rerun — task ${TASK_ID}"
echo "  dataset       = ${SAMPLE}     (DATASET_IDX=${DATASET_IDX})"
echo "  config        = ${CONFIG}     (CONFIG_IDX=${CONFIG_IDX})"
echo "  config_json   = ${CONFIG_JSON}"
echo "  chunks dir    = ${CHUNKS_DIR}"
echo "  outdir        = ${OUTDIR}"
echo "  walltime      = 36:00:00, mem = 64G"
echo "============================================================"

# Pre-flight checks
if [[ ! -f "${CONFIG_JSON}" ]]; then
    echo "ERROR: config JSON not found: ${CONFIG_JSON}" >&2
    exit 1
fi
if [[ ! -d "${CHUNKS_DIR}" ]]; then
    echo "ERROR: chunks dir not found: ${CHUNKS_DIR}" >&2
    exit 1
fi
# `find`, not a glob: 64 K chunk files exceed ARG_MAX (E2BIG)
N_FILTERED=$(find "${CHUNKS_DIR}" -maxdepth 1 -name '*_filtered.tsv.gz' | wc -l)
if [[ ${N_FILTERED} -eq 0 ]]; then
    echo "ERROR: no *_filtered.tsv.gz files in ${CHUNKS_DIR}" >&2
    exit 1
fi
echo "  found ${N_FILTERED} filtered files"

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
if [[ -f "${OUTDIR}/${SAMPLE}_cells_calls.tsv.gz" ]]; then
    SIZE=$(stat -c '%s' "${OUTDIR}/${SAMPLE}_cells_calls.tsv.gz" 2>/dev/null || \
           stat -f '%z' "${OUTDIR}/${SAMPLE}_cells_calls.tsv.gz" 2>/dev/null)
    echo "  ${SAMPLE}_cells_calls.tsv.gz: ${SIZE} bytes"
else
    echo "  WARNING: ${SAMPLE}_cells_calls.tsv.gz not found at expected path"
fi
echo "============================================================"
