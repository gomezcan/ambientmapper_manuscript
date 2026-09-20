#!/bin/bash
#SBATCH --job-name=extract_qc
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=_logs/03_07_extract_metrics_%A_%a.log
#SBATCH --array=0-14

# =============================================================================
# 03_07_extract_metrics.sh — Phase 5: Extract QCMapping from synthetic BAMs
#
# Parses BWA-mapped BAMs with pysam, writes per-genome QCMapping files.
# SLURM array 0-14: one task per dataset.
#
# Requires: Phase 4 output (synthetic/mapping/{dataset}/syn_{genome}.sorted.bam)
# Output:   synthetic/{dataset}/qc/{genome}_QCMapping.txt
# =============================================================================

REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"   # this repository (sourced helpers, sibling scripts)
PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"   # root of the data tree (3_Mapping/, 5_AmbientDetection/, ...)
cd "${PROJECT_ROOT}/5_AmbientDetection"

# conda activate <env from environment.yml>

# --- Dataset list (must match Phase 3/4) ---
DATASETS=(
    alpha_000
    alpha_002_Il14H
    alpha_005_Il14H
    alpha_010_Il14H
    alpha_020_Il14H
    alpha_030_Il14H
    alpha_040_Il14H
    alpha_050_Il14H
    alpha_002_Ki11
    alpha_005_Ki11
    alpha_010_Ki11
    alpha_020_Ki11
    alpha_030_Ki11
    alpha_040_Ki11
    alpha_050_Ki11
)

DATASET="${DATASETS[$SLURM_ARRAY_TASK_ID]}"

if [ -z "${DATASET}" ]; then
    echo "Error: Invalid task ID ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

BAM_DIR="synthetic/mapping/${DATASET}"
OUTDIR="synthetic/${DATASET}/qc"

echo "============================================"
echo "[$(date)] Phase 5: Extract QCMapping"
echo "  Task:    ${SLURM_ARRAY_TASK_ID}"
echo "  Dataset: ${DATASET}"
echo "  BAMs:    ${BAM_DIR}/syn_{B73v5,Il14H,Ki11}.sorted.bam"
echo "  Output:  ${OUTDIR}/"
echo "============================================"

python "${REPO_ROOT}/workflows/03_genotyping/synthetic/03_07_extract_metrics.py" \
    --bam-dir "${BAM_DIR}" \
    --outdir "${OUTDIR}" \
    --genome-map B73v5:B73 Il14H:Il14H Ki11:Ki11

EXIT_CODE=$?

echo ""
echo "============================================"
echo "[$(date)] Phase 5 complete for ${DATASET} (exit code: ${EXIT_CODE})"
echo "  Output files:"
for G in B73 Il14H Ki11; do
    F="${OUTDIR}/${G}_QCMapping.txt"
    if [ -f "${F}" ]; then
        N=$(wc -l < "${F}")
        echo "    ${F}: $((N-1)) reads"
    fi
done
echo "============================================"

exit ${EXIT_CODE}
