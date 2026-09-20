#!/usr/bin/env bash
#SBATCH --job-name=root1_assign_a005
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=_logs/02_12a_root1_assign_a005_%A.log
#
# 02_12a_pipeline_full_root1_alpha005.sh — Phase 4 prerequisite
#
# Re-runs the ambientmapper assign step on the FULL Root1_rep1 dataset using
# the α=0.05/k=10 reconciled bundle (the same bundle used for sub1k C0 in
# Phase 1).
#
# The assign step writes its α=0.05 output (*_filtered.tsv.gz) directly into
# the existing cell_map_ref_chunks/ dir. The legacy α=1e-6 output was backed
# up to cell_map_ref_chunks_alpha1e6_legacy_backup/ on a prior run and is NOT
# restored — Phase 4 numbers come from α=0.05, so cell_map_ref_chunks/
# becomes the canonical α=0.05 friendwith dir going forward.
#
# One new sibling dir is produced:
#
#   Root1_rep1/cell_map_ref_chunks_alpha005_friendwithout/
#       (post-hoc relabeled rescued → ambiguous, friend rescue = "OFF")
#
# Phase 4 genotyping reads directly from cell_map_ref_chunks/ for friend-ON
# configs, and from _alpha005_friendwithout/ for friend-OFF configs.
#
# Why a separate "_friendwithout" sibling: ambientmapper has no
# `--no-friend-rescue` flag at the assign step (the call at
# assign_streaming.py:1584 is unconditional). Since genotyping.py:719 weights
# are purely `cls == "ambiguous"`, rewriting rescued → ambiguous is
# functionally equivalent to disabling friend rescue from genotyping's
# perspective. The same trick is used in 02_09_pipeline_sub1k.sh:103-165.
#
# Single-task SLURM job (no array). Expected wall: ~120 hours for 9023 chunks.
# Resume-safe: resubmit after timeout and completed chunks are skipped.
#
# Prerequisites:
#   - configs/Root1_rep1.ambientmapper.json exists (full 26-genome config)
#   - The full Root1_rep1 extract+filter+chunks steps were already done in
#     a prior pipeline run (the existing Root1_rep1/cell_map_ref_chunks/ proves
#     this). We DO NOT re-run extract/filter/chunks here — only assign.
#   - ambientmapper is installed
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

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
SAMPLE=Root1_rep1
WORKDIR=${PROJECT_ROOT}/5_AmbientDetection/Root1_rep1
CONFIG=configs/Root1_rep1.ambientmapper.json
THREADS=16

# Canonical chunks dir — holds α=0.05 output after this script
CHUNKS_DIR=${WORKDIR}/cell_map_ref_chunks

# Friend rescue knockout sibling (rescued → ambiguous relabel)
WITHOUT_DIR=${WORKDIR}/cell_map_ref_chunks_alpha005_friendwithout

# Pin ambientmapper version for reproducibility
# AMBIENTMAPPER_REPO: optional local checkout of the ambientmapper tool, used only to log its commit
AM_VERSION=$(pip show ambientmapper 2>/dev/null | awk '/^Version/ {print $2}')
AM_COMMIT=$(git -C "${AMBIENTMAPPER_REPO:-.}" \
    rev-parse --short HEAD 2>/dev/null || echo "unknown")

echo "============================================================"
echo "[$(date)] Phase 4 pre-step: re-assign full Root1 at α=0.05"
echo "  sample        = ${SAMPLE}"
echo "  workdir       = ${WORKDIR}"
echo "  config        = ${CONFIG}"
echo "  chunks dir    = ${CHUNKS_DIR} (α=0.05 canonical)"
echo "  without dir   = ${WITHOUT_DIR}"
echo "  ambientmapper : v${AM_VERSION} @ ${AM_COMMIT}"
echo "============================================================"

if [[ ! -f "${CONFIG}" ]]; then
    echo "ERROR: config not found: ${CONFIG}" >&2
    exit 1
fi
if [[ ! -d "${CHUNKS_DIR}" ]]; then
    echo "ERROR: chunks dir not found: ${CHUNKS_DIR}" >&2
    echo "       The full Root1 extract/filter/chunks pipeline must already be done." >&2
    exit 1
fi

# Check that the chunks .txt files exist — these are the barcode chunks that
# the assign step consumes (alpha-independent, define the barcode partitioning).
N_CHUNK_TXT=$(find "${CHUNKS_DIR}" -maxdepth 1 -name '*_cell_map_ref_chunk_*.txt' | wc -l)
if [[ ${N_CHUNK_TXT} -eq 0 ]]; then
    echo "ERROR: no chunk .txt files in ${CHUNKS_DIR}" >&2
    echo "       Did the chunks step actually run? Check ambientmapper extract+filter+chunks output." >&2
    exit 1
fi
echo "  found ${N_CHUNK_TXT} existing chunk .txt files"

# -----------------------------------------------------------------------------
# Step 1 — run assign at α=0.05 against cell_map_ref_chunks/
# -----------------------------------------------------------------------------
#
# The assign step writes *_filtered.tsv.gz into cell_map_ref_chunks/.
# Legacy α=1e-6 output was backed up to _alpha1e6_legacy_backup/ on a prior
# run. We do NOT restore it — cell_map_ref_chunks/ becomes the canonical
# α=0.05 dir. Resume-safe: existing *_filtered.tsv.gz are skipped.
# -----------------------------------------------------------------------------

echo ""
echo "[$(date)] Step 1/2: assign full Root1 at α=0.05 (reconciled bundle)"
N_EXISTING=$(find "${CHUNKS_DIR}" -maxdepth 1 -name '*_filtered.tsv.gz' | wc -l)
echo "  ${N_EXISTING} filtered files already present (will be skipped on resume)"
ambientmapper assign \
    --config "${CONFIG}" \
    --threads "${THREADS}" \
    --alpha 0.05 \
    --k 10 \
    --mapq-min 10 \
    --xa-max 2 \
    --chunksize 500000 \
    --batch-size 6 \
    --score-batch-size 50 \
    --score-workers 2

# -----------------------------------------------------------------------------
# Step 2 — create friendwithout sibling via rescued→ambiguous relabel
# -----------------------------------------------------------------------------
N_NEW_FILTERED=$(find "${CHUNKS_DIR}" -maxdepth 1 -name '*_filtered.tsv.gz' | wc -l)
if [[ ${N_NEW_FILTERED} -eq 0 ]]; then
    echo "ERROR: assign step produced no *_filtered.tsv.gz in ${CHUNKS_DIR}" >&2
    exit 1
fi
echo ""
echo "[$(date)] Step 2/2: create friendwithout sibling (${N_NEW_FILTERED} files)"

rm -rf "${WITHOUT_DIR}"
mkdir -p "${WITHOUT_DIR}"

# Copy chunk .txt files (for discoverability — alpha-independent).
find "${CHUNKS_DIR}" -maxdepth 1 -name '*_cell_map_ref_chunk_*.txt' \
    -exec cp -t "${WITHOUT_DIR}/" {} +

# Relabel rescued → ambiguous into the friendwithout sibling
# (Reuses the python helper from 02_09_pipeline_sub1k.sh:127-164 verbatim)
echo "[$(date)] relabel rescued → ambiguous for friendwithout sibling"
python - <<PYEOF
import glob
import gzip
import os
import sys

src = "${CHUNKS_DIR}"
dst = "${WITHOUT_DIR}"

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
    print("  WARNING: no rescued reads found — friend rescue may not have run")
PYEOF

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
echo ""
echo "============================================================"
echo "[$(date)] Phase 4 pre-step COMPLETE — full Root1 α=0.05 chunks dirs ready"
echo ""
echo "Chunks directories under ${WORKDIR}:"
for d in cell_map_ref_chunks \
         cell_map_ref_chunks_alpha005_friendwithout; do
    target="${WORKDIR}/${d}"
    n_filt=$(find "${target}" -maxdepth 1 -name '*_filtered.tsv.gz' | wc -l)
    n_txt=$(find "${target}" -maxdepth 1 -name '*_cell_map_ref_chunk_*.txt' | wc -l)
    echo "  ${d}: ${n_filt} filtered files, ${n_txt} chunk .txt files"
done

echo ""
echo "assigned_class distribution (sampling first _filtered.tsv.gz per dir):"
python - <<PYEOF
import glob
import gzip
import os
from collections import Counter

root = "${WORKDIR}"
for d in ("cell_map_ref_chunks",
          "cell_map_ref_chunks_alpha005_friendwithout"):
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
echo "Next: submit 02_12_genotyping_full_root1_phase4.sh after the Phase 3 winner"
echo "      stacked config has been picked from Phase 3 eval."
echo "============================================================"
