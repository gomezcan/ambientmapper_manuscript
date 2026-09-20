#!/usr/bin/env bash
#SBATCH --job-name=syn_friend_relabel
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=_logs/03_21a_friend_relabel_%A_%a.log
#SBATCH --array=0-29
#
# 03_21a_synthetic_friend_relabel.sh — Phase 2 friend rescue knockout
#
# For each (track, dataset) pair, copies the existing
# `cell_map_ref_chunks/` directory to `cell_map_ref_chunks_friendwithout/`
# and rewrites the `assigned_class` column from "rescued" to "ambiguous" in
# every `*_filtered.tsv.gz`. This is the synthetic equivalent of the
# friend-rescue knockout that 02_09_pipeline_sub1k.sh produces for sub1k
# panels A and B.
#
# Reuses the python relabel block from 02_09:127-164 verbatim. Column
# ordering and gzipping are verified identical between sub1k and synthetic
# chunks files.
#
# Why post-hoc relabel: ambientmapper has no `--no-friend-rescue` flag at
# the assign step (the call at assign_streaming.py:1584 is unconditional).
# Since genotyping.py:719 weights are purely `cls == "ambiguous"`, rewriting
# rescued labels to ambiguous is functionally equivalent to disabling friend
# rescue from the genotyping step's perspective.
#
# Output: a sibling dir next to the original — `cell_map_ref_chunks/` is
# preserved unchanged (it acts as the "friendwith" variant for Phase 2).
#
# Array layout: 30 tasks = 2 tracks × 15 datasets.
#   tasks  0-14:  synthetic/      (Track B)
#   tasks 15-29:  synthetic_disc/ (Track B-disc)
#
# Prerequisites:
#   - synthetic/<dataset>/cell_map_ref_chunks/*_filtered.tsv.gz exist
#   - synthetic_disc/<dataset>/cell_map_ref_chunks/*_filtered.tsv.gz exist
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

# --- Datasets (15) ---
DATASETS=(
    alpha_000
    alpha_002_Il14H alpha_005_Il14H alpha_010_Il14H alpha_020_Il14H
    alpha_030_Il14H alpha_040_Il14H alpha_050_Il14H
    alpha_002_Ki11  alpha_005_Ki11  alpha_010_Ki11  alpha_020_Ki11
    alpha_030_Ki11  alpha_040_Ki11  alpha_050_Ki11
)

# --- Tracks (2) ---
TRACKS=(synthetic synthetic_disc)
N_DATASETS=${#DATASETS[@]}   # 15

TASK_ID=${SLURM_ARRAY_TASK_ID:-0}
TRACK_IDX=$(( TASK_ID / N_DATASETS ))
DS_IDX=$(( TASK_ID % N_DATASETS ))

if [[ ${TRACK_IDX} -ge ${#TRACKS[@]} ]]; then
    echo "ERROR: TRACK_IDX=${TRACK_IDX} out of range" >&2
    exit 1
fi

TRACK="${TRACKS[$TRACK_IDX]}"
DATASET="${DATASETS[$DS_IDX]}"

WORKDIR="${PROJECT_ROOT}/5_AmbientDetection"
DS_DIR="${WORKDIR}/${TRACK}/${DATASET}"
SRC_DIR="${DS_DIR}/cell_map_ref_chunks"
DST_DIR="${DS_DIR}/cell_map_ref_chunks_friendwithout"

echo "============================================================"
echo "[$(date)] 03_21a friend relabel — task ${TASK_ID}"
echo "  track   = ${TRACK}"
echo "  dataset = ${DATASET}"
echo "  src     = ${SRC_DIR}"
echo "  dst     = ${DST_DIR}"
echo "============================================================"

if [[ ! -d "${SRC_DIR}" ]]; then
    echo "ERROR: source chunks dir not found: ${SRC_DIR}" >&2
    exit 1
fi

N_FILTERED=$(ls "${SRC_DIR}"/*_filtered.tsv.gz 2>/dev/null | wc -l)
if [[ ${N_FILTERED} -eq 0 ]]; then
    echo "ERROR: no *_filtered.tsv.gz files in ${SRC_DIR}" >&2
    exit 1
fi
echo "  found ${N_FILTERED} filtered files in source"

# Wipe and recreate destination
rm -rf "${DST_DIR}"
mkdir -p "${DST_DIR}"

# 1) Copy chunk .txt files to the friendwithout variant (for discoverability;
#    they're identical between with/without — only the filtered tsvs differ)
shopt -s nullglob
TXT_FILES=( "${SRC_DIR}"/*_cell_map_ref_chunk_*.txt )
shopt -u nullglob
if [[ ${#TXT_FILES[@]} -gt 0 ]]; then
    cp "${TXT_FILES[@]}" "${DST_DIR}/"
    echo "  copied ${#TXT_FILES[@]} chunk .txt files"
fi

# 2) Relabel rescued→ambiguous into the _friendwithout variant
#    (Reuses the python helper from 02_09_pipeline_sub1k.sh:127-164 verbatim)
python - <<PYEOF
import glob
import gzip
import os
import sys

src = "${SRC_DIR}"
dst = "${DST_DIR}"

files = sorted(glob.glob(os.path.join(src, "*_filtered.tsv.gz")))
if not files:
    sys.exit("ERROR: no filtered files in " + src)

n_files = 0
n_reads = 0
n_relabeled = 0
for f in files:
    base = os.path.basename(f)
    out = os.path.join(dst, base)
    with gzip.open(f, "rt") as fin, gzip.open(out, "wt") as fout:
        header = fin.readline()
        fout.write(header)
        cols = header.rstrip("\n").split("\t")
        try:
            i_cls = cols.index("assigned_class")
        except ValueError:
            sys.exit(f"ERROR: no assigned_class col in {f}, got {cols}")
        for line in fin:
            parts = line.rstrip("\n").split("\t")
            n_reads += 1
            if parts[i_cls] == "rescued":
                parts[i_cls] = "ambiguous"
                n_relabeled += 1
            fout.write("\t".join(parts) + "\n")
    n_files += 1
print(f"  relabeled {n_relabeled:,} rescued→ambiguous across {n_files} files ({n_reads:,} reads)")
if n_relabeled == 0:
    print("  WARNING: no rescued reads found — verify upstream assign step ran with friend rescue enabled")
PYEOF

echo ""
echo "============================================================"
echo "[$(date)] Done: ${TRACK}/${DATASET}"
echo "  output dir: ${DST_DIR}"
echo "============================================================"
