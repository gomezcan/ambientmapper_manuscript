#!/bin/bash
#SBATCH --job-name=assign_barcodes
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=_logs/03_05_assign_barcodes_%j.log

# =============================================================================
# 03_05_assign_barcodes.sh — Phase 3: Assign reads to synthetic barcodes
#
# Titration design: 208 template barcodes × 15 datasets (1 baseline +
# 7 alpha levels × 2 contaminants). Cell reads are shared across datasets;
# only contamination reads change.
#
# Requires: Phase 2 output (synthetic/reads/{B73,Il14H,Ki11}_{1,2}.fq.gz)
# Output:   synthetic/barcoded/templates.tsv + 15 dataset directories
# =============================================================================

REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"   # this repository (sourced helpers, sibling scripts)
PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"   # root of the data tree (3_Mapping/, 5_AmbientDetection/, ...)
cd "${PROJECT_ROOT}/5_AmbientDetection"

# conda activate <env from environment.yml>

echo "============================================"
echo "[$(date)] Phase 3: Barcode assignment (titration)"
echo "  Host: $(hostname)"
echo "  Memory: ${SLURM_MEM_PER_NODE}MB"
echo "============================================"

python "${REPO_ROOT}/workflows/03_genotyping/synthetic/03_05_assign_barcodes.py" \
  --reads-dir synthetic/reads \
  --outdir synthetic/barcoded \
  --seed 42 \
  --depth-jitter 0.10 \
  --n-singlet-per-bin 25 \
  --n-doublet-reps 2 \
  --n-empty 10

echo ""
echo "============================================"
echo "[$(date)] Phase 3 complete (exit code: $?)"
echo "============================================"
