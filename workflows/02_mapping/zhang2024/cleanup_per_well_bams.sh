#!/usr/bin/env bash
# cleanup_per_well_bams.sh
#
# Remove per-well BAM intermediates after the merge is complete.
# Keeps: bc_counts, metrics, proper_pairs, BEDs, sentinels.
# Removes: raw.bam, mq10.BC.bam, mq10.BC.rmdup.bam, mq10.BC.rmdup.mm.bam + .bai
#
# Run on the HPC for speed. Verify that the merged BAMs exist before running.
#
# Usage (from ${PROJECT_ROOT}/3_Mapping):
#   bash cleanup_per_well_bams.sh B73Mo17_rep1
#   bash cleanup_per_well_bams.sh B73Mo17_rep2
#   bash cleanup_per_well_bams.sh multiGenotypes_rep1  # after its merge completes

set -euo pipefail
: "${PROJECT_ROOT:?set PROJECT_ROOT to the analysis tree that holds 3_Mapping/}"

DATASET_REP="${1:?Usage: bash cleanup_per_well_bams.sh <DATASET_REP>}"
BASEDIR="${PROJECT_ROOT}/3_Mapping"
WELL_DIR="${BASEDIR}/${DATASET_REP}"
MERGED_DIR="${BASEDIR}/ambientmapper_input"

# Safety check: verify merged BAMs exist
n_merged=$(find "$MERGED_DIR" -name "${DATASET_REP}_*_scifiATAC.mq10.BC.rmdup.mm.bam" 2>/dev/null | wc -l)
if [[ "$n_merged" -eq 0 ]]; then
  echo "ERROR: No merged BAMs found in ${MERGED_DIR} for ${DATASET_REP}"
  echo "Run the merge job first before cleaning up per-well BAMs."
  exit 1
fi
echo "[$(date)] Found $n_merged merged BAMs for ${DATASET_REP}: safe to clean"

# Count before
before=$(du -sh "$WELL_DIR" | cut -f1)

# Remove per-well BAMs (all intermediate + final per-well)
echo "[$(date)] Removing per-well BAMs..."
find "$WELL_DIR" -name "*.raw.bam" -delete -print | wc -l | xargs -I{} echo "  raw.bam: {} files"
find "$WELL_DIR" -name "*.rawSort.bam" -delete -print | wc -l | xargs -I{} echo "  rawSort.bam: {} files"
find "$WELL_DIR" -name "*.mq10.BC.bam" -delete -print | wc -l | xargs -I{} echo "  mq10.BC.bam: {} files"
find "$WELL_DIR" -name "*.mq10.BC.rmdup.bam" -delete -print | wc -l | xargs -I{} echo "  mq10.BC.rmdup.bam: {} files"
find "$WELL_DIR" -name "*.mq10.BC.rmdup.mm.bam" -delete -print | wc -l | xargs -I{} echo "  mq10.BC.rmdup.mm.bam: {} files"
find "$WELL_DIR" -name "*.bam.bai" -delete -print | wc -l | xargs -I{} echo "  bam.bai: {} files"
find "$WELL_DIR" -name "*.raw.sam" -delete -print | wc -l | xargs -I{} echo "  raw.sam: {} files"

# Clean sentinels (no longer needed)
find "$WELL_DIR" -name "*.ok.json" -delete -print | wc -l | xargs -I{} echo "  sentinels: {} files"
find "$WELL_DIR" -type d -name "_sentinels" -empty -delete 2>/dev/null

after=$(du -sh "$WELL_DIR" | cut -f1)
echo ""
echo "[$(date)] ${DATASET_REP}: ${before} -> ${after}"
echo "Kept: bc_counts, metrics, proper_pairs, BEDs"
