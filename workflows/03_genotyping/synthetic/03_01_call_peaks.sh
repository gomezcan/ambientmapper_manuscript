#!/bin/bash
#SBATCH --job-name=macs3_peaks
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=_logs/03_01_call_peaks_%j.log

# =============================================================================
# 03_01_call_peaks.sh — Phase 0: Peak calling on Root1_rep1 B73 BAM
#
# Calls accessible chromatin regions (ACRs) using MACS3 on the cleaned
# Root1_rep1 B73 BAM. Post-processes to filter for autosomes, minimum
# peak width, and top peaks by q-value.
#
# Output: synthetic/peaks/B73_peaks_filtered.bed
# =============================================================================

PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"   # root of the data tree (3_Mapping/, 5_AmbientDetection/, ...)
cd "${PROJECT_ROOT}/5_AmbientDetection"

# conda activate <env from environment.yml>

# Input BAM of the run of record. 4_MappingCleaning/ was later consolidated into
# 3_Mapping/ (the same file now sits under 3_Mapping/Root1_rep1/).
BAM="${PROJECT_ROOT}/4_MappingCleaning/Root1_rep1/Root1_rep1_B73v5.mq10.BC.rmdup.mm.bam"
OUTDIR="synthetic/peaks"
PREFIX="Root1_rep1_B73"
TOP_N=10000
MIN_WIDTH=150

mkdir -p "${OUTDIR}"

echo "============================================"
echo "[$(date)] Phase 0: Peak calling"
echo "  BAM:    ${BAM}"
echo "  Output: ${OUTDIR}"
echo "============================================"

# --- Step 1: MACS3 peak calling ---
echo "[$(date)] Running MACS3 callpeak..."

macs3 callpeak \
  -t "${BAM}" \
  -f BAMPE \
  --nomodel \
  --shift -75 \
  --extsize 150 \
  -g 1.6e9 \
  --keep-dup all \
  -q 0.05 \
  -n "${PREFIX}" \
  --outdir "${OUTDIR}" \
  2>&1

echo "[$(date)] MACS3 done."
echo "  Raw peaks: $(wc -l < "${OUTDIR}/${PREFIX}_peaks.narrowPeak")"

# --- Step 2: Post-processing ---
# Filter to:
#   1. Autosomes only (chr1-chr10 for maize B73v5, named "chr1"..."chr10" or "1"..."10")
#   2. Peak width >= MIN_WIDTH bp
#   3. Top TOP_N peaks by q-value (column 9 in narrowPeak, -log10 scale, higher = better)

echo "[$(date)] Filtering peaks..."

NARROWPEAK="${OUTDIR}/${PREFIX}_peaks.narrowPeak"
FILTERED="${OUTDIR}/B73_peaks_filtered.bed"

awk -v min_w="${MIN_WIDTH}" '
BEGIN { OFS="\t" }
{
    # Accept maize autosome names: chr1-chr10, Zm_chr1-Zm_chr10, 1-10
    chr = $1
    if (chr ~ /^(chr|Zm_chr)?[0-9]+$/) {
        # Extract the chromosome number
        gsub(/^(chr|Zm_chr)/, "", chr)
        num = int(chr)
        if (num >= 1 && num <= 10) {
            width = $3 - $2
            if (width >= min_w) {
                print $0
            }
        }
    }
}' "${NARROWPEAK}" \
| sort -k9,9gr \
| head -n "${TOP_N}" \
| sort -k1,1 -k2,2n \
> "${FILTERED}"

N_FILTERED=$(wc -l < "${FILTERED}")
echo "[$(date)] Filtered peaks: ${N_FILTERED} (top ${TOP_N}, autosomes, width >= ${MIN_WIDTH}bp)"

# --- Summary stats ---
echo ""
echo "=== Peak summary ==="
echo "  Raw narrowPeak: $(wc -l < "${NARROWPEAK}")"
echo "  Filtered BED:   ${N_FILTERED}"
if [ "${N_FILTERED}" -gt 0 ]; then
    echo "  Width range:    $(awk '{print $3-$2}' "${FILTERED}" | sort -n | head -1) - $(awk '{print $3-$2}' "${FILTERED}" | sort -n | tail -1) bp"
    echo "  Median width:   $(awk '{print $3-$2}' "${FILTERED}" | sort -n | awk 'NR==int(NR/2){print}')"
    echo "  Q-value range:  $(awk '{print $9}' "${FILTERED}" | sort -n | head -1) - $(awk '{print $9}' "${FILTERED}" | sort -n | tail -1)"
fi

echo ""
echo "============================================"
echo "[$(date)] Phase 0 complete."
echo "  Output: ${FILTERED}"
echo "============================================"
