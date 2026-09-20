#!/usr/bin/env bash
#SBATCH --job-name=sm2v2_merge
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=_logs/01_00_sm2v2_merge_%A.log
#
# 01_00_merge_SM2v2_inputs.sh: merge the per-library SM2 BAMs into one BAM per reference.
#
# The legacy chain writes one cleaned BAM per plate half (SM2_B73, SM2_At) and reference. AmbientMapper
# reads one BAM per reference for the whole library, so the two halves are merged here into
# 3_Mapping/ambientmapper_input/, the layout also used by the Zhang et al. 2024 libraries.
#
# Inputs (per reference, merge the two plate halves SM2_B73 + SM2_At):
#   B73v5 reference:
#     3_Mapping/SM2_B73/SM2_B73_Zm_B73v5_scifiATAC.mq10.BC.rmdup.mm.bam
#     3_Mapping/SM2_At/SM2_At_Zm_B73v5_scifiATAC.mq10.BC.rmdup.mm.bam
#   TAIR10 reference:
#     3_Mapping/SM2_B73/SM2_B73_AraTAIR10_scifiATAC.mq10.BC.rmdup.mm.bam
#     3_Mapping/SM2_At/SM2_At_AraTAIR10_scifiATAC.mq10.BC.rmdup.mm.bam
#
# Outputs (merged + indexed):
#   3_Mapping/ambientmapper_input/SM2_B73v5_scifiATAC.mq10.BC.rmdup.mm.bam
#   3_Mapping/ambientmapper_input/SM2_TAIR10_scifiATAC.mq10.BC.rmdup.mm.bam
#
# Resume-safe: skips a merge if the output BAM + .bai already exist.
# Usage (from ${PROJECT_ROOT}/5_AmbientDetection): sbatch 01_00_merge_SM2v2_inputs.sh
# =============================================================================

set -euo pipefail
: "${PROJECT_ROOT:?set PROJECT_ROOT to the analysis tree that holds 3_Mapping/ and 5_AmbientDetection/}"

cd "${PROJECT_ROOT}/5_AmbientDetection"
mkdir -p _logs

# conda activate <env from environment.yml>   (samtools)

THREADS="${SLURM_CPUS_PER_TASK:-8}"

MAP="${PROJECT_ROOT}/3_Mapping"
OUT=${MAP}/ambientmapper_input
mkdir -p "${OUT}"

merge_one() {
    local LABEL=$1
    local OUTBAM=$2
    shift 2
    local INPUTS=("$@")

    if [[ -f "${OUTBAM}" && -f "${OUTBAM}.bai" ]]; then
        echo "[skip] ${LABEL}: output already exists"
        ls -lh "${OUTBAM}" "${OUTBAM}.bai"
        return 0
    fi

    for f in "${INPUTS[@]}"; do
        [[ -f "${f}" ]] || { echo "ERROR: missing input ${f}" >&2; exit 1; }
    done

    echo ""
    echo "[$(date)] merging ${LABEL}"
    echo "  inputs:"
    for f in "${INPUTS[@]}"; do echo "    ${f}"; done
    echo "  output: ${OUTBAM}"

    samtools merge -@ "${THREADS}" -f "${OUTBAM}" "${INPUTS[@]}"
    samtools index -@ "${THREADS}" "${OUTBAM}"

    echo "[$(date)] ${LABEL} done"
    ls -lh "${OUTBAM}" "${OUTBAM}.bai"
}

echo "============================================================"
echo "[$(date)] merge SM2_B73 + SM2_At per-genome BAMs"
echo "  threads = ${THREADS}"
echo "  out dir = ${OUT}"
echo "============================================================"

merge_one "B73v5" \
    "${OUT}/SM2_B73v5_scifiATAC.mq10.BC.rmdup.mm.bam" \
    "${MAP}/SM2_B73/SM2_B73_Zm_B73v5_scifiATAC.mq10.BC.rmdup.mm.bam" \
    "${MAP}/SM2_At/SM2_At_Zm_B73v5_scifiATAC.mq10.BC.rmdup.mm.bam"

merge_one "TAIR10" \
    "${OUT}/SM2_TAIR10_scifiATAC.mq10.BC.rmdup.mm.bam" \
    "${MAP}/SM2_B73/SM2_B73_AraTAIR10_scifiATAC.mq10.BC.rmdup.mm.bam" \
    "${MAP}/SM2_At/SM2_At_AraTAIR10_scifiATAC.mq10.BC.rmdup.mm.bam"

echo ""
echo "============================================================"
echo "[$(date)] merge DONE"
ls -lh "${OUT}/SM2_B73v5_scifiATAC.mq10.BC.rmdup.mm.bam" \
       "${OUT}/SM2_TAIR10_scifiATAC.mq10.BC.rmdup.mm.bam"
echo "============================================================"
