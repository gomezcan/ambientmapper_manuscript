#!/usr/bin/env bash
#SBATCH --job-name=sub1k_pipeline
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --array=0-1
#SBATCH --output=_logs/02_09_pipeline_sub1k_%A_%a.log
#
# 02_09_pipeline_sub1k.sh — AmbientMapper pipeline for Root1 sub1k panels
#
# Runs extract → filter → chunks → assign (×2 alpha bundles) per panel, then
# creates "friend rescue knockout" variants by post-hoc relabeling the
# assigned_class column (rescued → ambiguous) in the filtered output.
#
# Why post-hoc relabel: ambientmapper's _friend_rescue is hardcoded in
# assign_streaming.py:1584,1873 (no --no-friend-rescue flag). Since
# genotyping.py:719 classifies weights by `cls == "ambiguous"`, rewriting
# rescued labels to ambiguous is functionally equivalent to disabling friend
# rescue from the genotyping step's perspective.
#
# Array layout: 2 tasks (one per panel).
#   task 0 = sub1k_A
#   task 1 = sub1k_B
#
# Each task produces 4 chunks-style directories under
# Root1_rep1/sub1k_<panel>/:
#
#   cell_map_ref_chunks_alpha005_friendwith/      (Track B reconciled bundle, FR on)
#   cell_map_ref_chunks_alpha005_friendwithout/   (relabeled rescued→ambiguous)
#   cell_map_ref_chunks_alpha1e6_friendwith/      (legacy Root1 bundle,      FR on)
#   cell_map_ref_chunks_alpha1e6_friendwithout/   (relabeled rescued→ambiguous)
#
# Prerequisites:
#   1. 02_07_subsample_root1_balanced.py has been run (panel barcode lists)
#   2. 02_08_filter_bams_subsample.sh has completed (panel BAMs under sub1k_<panel>/bams/)
#   3. configs/Root1_rep1_sub1k_{A,B}.ambientmapper.json exist
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

# --- Task → panel mapping ---
PANELS=(A B)
PANEL=${PANELS[$SLURM_ARRAY_TASK_ID]}
SAMPLE=sub1k_${PANEL}
WORKDIR=${PROJECT_ROOT}/5_AmbientDetection/Root1_rep1
CONFIG=configs/Root1_rep1_sub1k_${PANEL}.ambientmapper.json

SAMPLE_DIR=${WORKDIR}/${SAMPLE}
CHUNKS_LIVE=${SAMPLE_DIR}/cell_map_ref_chunks
THREADS=16

# Pin ambientmapper version for reproducibility
# AMBIENTMAPPER_REPO: optional local checkout of the ambientmapper tool, used only to log its commit
AM_VERSION=$(pip show ambientmapper 2>/dev/null | awk '/^Version/ {print $2}')
AM_COMMIT=$(git -C "${AMBIENTMAPPER_REPO:-.}" \
    rev-parse --short HEAD 2>/dev/null || echo "unknown")

echo "============================================================"
echo "[$(date)] sub1k pipeline — panel ${PANEL}"
echo "  sample       = ${SAMPLE}"
echo "  workdir      = ${WORKDIR}"
echo "  sample dir   = ${SAMPLE_DIR}"
echo "  config       = ${CONFIG}"
echo "  ambientmapper: v${AM_VERSION} @ ${AM_COMMIT}"
echo "============================================================"

if [[ ! -f "${CONFIG}" ]]; then
    echo "ERROR: config not found: ${CONFIG}" >&2
    exit 1
fi
if [[ ! -d "${SAMPLE_DIR}/bams" ]]; then
    echo "ERROR: panel BAM dir not found: ${SAMPLE_DIR}/bams" >&2
    echo "       Run 02_08_filter_bams_subsample.sh first." >&2
    exit 1
fi

# -----------------------------------------------------------------------------
# Helper: archive & relabel
# -----------------------------------------------------------------------------
#
# After an assign run, this takes the live chunks dir's filtered output and:
#   1) copies it to <archive>_friendwith/ along with the chunk .txt files
#   2) creates <archive>_friendwithout/ by relabeling rescued→ambiguous
#
# $1 = label suffix (e.g. "alpha005", "alpha1e6")
# -----------------------------------------------------------------------------
archive_and_relabel() {
    local label="$1"
    local with_dir="${SAMPLE_DIR}/cell_map_ref_chunks_${label}_friendwith"
    local without_dir="${SAMPLE_DIR}/cell_map_ref_chunks_${label}_friendwithout"

    echo "[$(date)] archive_and_relabel: ${label}"

    rm -rf "${with_dir}" "${without_dir}"
    mkdir -p "${with_dir}" "${without_dir}"

    # 1) Copy chunk .txt files to both variants (for discoverability)
    cp "${CHUNKS_LIVE}"/*_cell_map_ref_chunk_*.txt "${with_dir}/"
    cp "${CHUNKS_LIVE}"/*_cell_map_ref_chunk_*.txt "${without_dir}/"

    # 2) Move filtered files from live dir to the _friendwith archive
    local n_filtered
    n_filtered=$(ls "${CHUNKS_LIVE}"/*_filtered.tsv.gz 2>/dev/null | wc -l)
    if [[ ${n_filtered} -eq 0 ]]; then
        echo "  ERROR: no *_filtered.tsv.gz found in ${CHUNKS_LIVE}" >&2
        exit 1
    fi
    mv "${CHUNKS_LIVE}"/*_filtered.tsv.gz "${with_dir}/"
    echo "  moved ${n_filtered} filtered files to ${with_dir}"

    # 3) Relabel rescued→ambiguous into the _friendwithout variant
    python - <<PYEOF
import glob
import gzip
import os
import sys

src = "${with_dir}"
dst = "${without_dir}"

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
PYEOF
}

# -----------------------------------------------------------------------------
# Step 1: extract
# -----------------------------------------------------------------------------
echo ""
echo "[$(date)] Step 1/4: extract"
ambientmapper extract --config "${CONFIG}" --threads "${THREADS}" --no-resume

# -----------------------------------------------------------------------------
# Step 2: filter
# -----------------------------------------------------------------------------
echo ""
echo "[$(date)] Step 2/4: filter"
ambientmapper filter --config "${CONFIG}" --threads 2 --min-barcode-freq 5 --no-resume

# -----------------------------------------------------------------------------
# Step 3: chunks
# -----------------------------------------------------------------------------
echo ""
echo "[$(date)] Step 3/4: chunks"
ambientmapper chunks --config "${CONFIG}" --chunk-size-cells 100 --no-resume

N_CHUNKS=$(ls "${CHUNKS_LIVE}"/*_cell_map_ref_chunk_*.txt 2>/dev/null | wc -l)
echo "  created ${N_CHUNKS} chunk files"

# -----------------------------------------------------------------------------
# Step 4a: assign with Track B reconciled bundle (alpha=0.05)
# -----------------------------------------------------------------------------
echo ""
echo "[$(date)] Step 4a: assign alpha=0.05 (Track B reconciled bundle)"
ambientmapper assign \
    --config "${CONFIG}" \
    --threads "${THREADS}" \
    --alpha 0.05 \
    --k 10 \
    --mapq-min 10 \
    --xa-max 2 \
    --chunksize 500000 \
    --batch-size 6 \
    --no-resume

archive_and_relabel "alpha005"

# -----------------------------------------------------------------------------
# Step 4b: assign with legacy Root1 bundle (alpha=1e-6)
# -----------------------------------------------------------------------------
echo ""
echo "[$(date)] Step 4b: assign alpha=1e-6 (legacy Root1 bundle)"
ambientmapper assign \
    --config "${CONFIG}" \
    --threads "${THREADS}" \
    --alpha 0.000001 \
    --k 5 \
    --mapq-min 10 \
    --xa-max 0 \
    --chunksize 500000 \
    --batch-size 6 \
    --skip-edges \
    --skip-ecdf \
    --no-resume

archive_and_relabel "alpha1e6"

# -----------------------------------------------------------------------------
# Summary + quick assigned_class distribution per chunks dir
# -----------------------------------------------------------------------------
echo ""
echo "============================================================"
echo "[$(date)] sub1k pipeline COMPLETE — panel ${PANEL}"
echo ""
echo "Chunks directories produced under ${SAMPLE_DIR}:"
for d in cell_map_ref_chunks_alpha005_friendwith \
         cell_map_ref_chunks_alpha005_friendwithout \
         cell_map_ref_chunks_alpha1e6_friendwith \
         cell_map_ref_chunks_alpha1e6_friendwithout; do
    target="${SAMPLE_DIR}/${d}"
    n=$(ls "${target}"/*_filtered.tsv.gz 2>/dev/null | wc -l)
    echo "  ${d}: ${n} filtered files"
done

echo ""
echo "assigned_class distribution (sampling first _filtered.tsv.gz per dir):"
python - <<PYEOF
import glob
import gzip
import os
from collections import Counter

root = "${SAMPLE_DIR}"
for d in ("cell_map_ref_chunks_alpha005_friendwith",
          "cell_map_ref_chunks_alpha005_friendwithout",
          "cell_map_ref_chunks_alpha1e6_friendwith",
          "cell_map_ref_chunks_alpha1e6_friendwithout"):
    files = sorted(glob.glob(os.path.join(root, d, "*_filtered.tsv.gz")))
    if not files:
        print(f"  {d}: (empty)")
        continue
    f = files[0]
    c = Counter()
    with gzip.open(f, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        i = header.index("assigned_class")
        for line in fh:
            c[line.rstrip("\n").split("\t")[i]] += 1
    total = sum(c.values())
    parts = []
    for k in ("winner", "rescued", "ambiguous"):
        n = c.get(k, 0)
        parts.append(f"{k}={n:,} ({100*n/total:.1f}%)" if total else f"{k}=0")
    print(f"  {d}/{os.path.basename(f)}: {' '.join(parts)}")
PYEOF

echo "============================================================"
