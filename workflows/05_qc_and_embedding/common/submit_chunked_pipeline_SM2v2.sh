#!/usr/bin/env bash
#
# submit_chunked_pipeline_SM2v2.sh — per-sample 4-stage afterok chain.
#
# Replaces the monolithic 1_QC_scifiATAC_SM2v2.sh. Submits:
#   1) 1_1a_chunk_bed.sh         hash-partition BED into N chunks   (1 task,  8G,  6 cpu)
#   2) 1_1b_per_chunk.sh         per-chunk Socrates build (array)   (N tasks, 130G)
#   3) 1_1c_merge_and_qc.sh      mergeSocratesRDS + isCellv2        (1 task,  180G)
#   4) 1_2_1_3_run.sh            filter + meta-QC cascade            (1 task,  20G)
#
# Each stage depends on the previous via --dependency=afterok. Sentinels in
# each stage make the chain idempotent — re-submitting after a partial failure
# skips already-complete stages.
#
# Prerequisites:
#   - MACS3 peaks for <POOL> exist at _data/_PeakFiles/<POOL>/<POOL>_peaks.narrowPeak
#     (produced by sbatch 0_scripts/setup/0_07_macs3_SM2v2.sh).
#   - tn5 BED for <POOL> exists at
#     _data/_BED_files/<POOL>_ZmATcombined_scifiATAC.mq10.BC.rmdup.mm.tn5.bed.gz
#
# Usage (per sample):
#   POOL=SM2_At          OUTDIR=SM2          bash 0_scripts/common/submit_chunked_pipeline_SM2v2.sh
#   POOL=SM2_B73         OUTDIR=SM2          bash 0_scripts/common/submit_chunked_pipeline_SM2v2.sh
#   POOL=Clean.SM2v2_At   OUTDIR=SM2v2_clean bash 0_scripts/common/submit_chunked_pipeline_SM2v2.sh
#   POOL=Clean.SM2v2_B73  OUTDIR=SM2v2_clean bash 0_scripts/common/submit_chunked_pipeline_SM2v2.sh
#
# Optional knobs:
#   N=5             number of barcode-hash chunks (default 5)
#   DEPTH=200       depth floor (unique insertions) for 1_3 meta-QC (default 200)
#   KEEP_CHUNKS=1   skip cleanup of per-chunk intermediates after stage 3

set -euo pipefail

POOL="${POOL:?POOL not set; usage: POOL=<sample> OUTDIR=<dir> [N=5] bash submit_chunked_pipeline_SM2v2.sh}"
OUTDIR="${OUTDIR:?OUTDIR not set; SM2 for PreClean, SM2v2_clean for PostClean}"
N="${N:-5}"
DEPTH="${DEPTH:-200}"   # depth floor for 1_3 meta-QC; 200 matches the established SM2 analysis

PROJECT_ROOT="${PROJECT_ROOT:-.}"
BASE="${PROJECT_ROOT}/6_socrates"
SCRIPTS="${SCRIPTS:-${BASE}/0_scripts}"
LOGDIR="${BASE}/_logs"
mkdir -p "$LOGDIR"

# Pre-flight: MACS3 peaks must exist (0_07 run independently before this).
PEAKS="${BASE}/_data/_PeakFiles/${POOL}/${POOL}_peaks.narrowPeak"
if [ ! -s "$PEAKS" ]; then
  echo "ERROR: MACS3 peaks not found for $POOL:" >&2
  echo "       $PEAKS" >&2
  echo "       Run 'sbatch 0_scripts/setup/0_07_macs3_SM2v2.sh' first." >&2
  exit 1
fi

# Pre-flight: tn5 BED must exist.
BED="${BASE}/_data/_BED_files/${POOL}_ZmATcombined_scifiATAC.mq10.BC.rmdup.mm.tn5.bed.gz"
if [ ! -s "$BED" ]; then
  echo "ERROR: tn5 BED not found for $POOL:" >&2
  echo "       $BED" >&2
  echo "       For Clean.SM2v2_* samples run 0_06_make_tn5bed_SM2v2_combined.sh first." >&2
  exit 1
fi

echo " - chunked Socrates QC chain for $POOL"
echo "   OUTDIR : $OUTDIR"
echo "   chunks : $N"
echo "   peaks  : $PEAKS"
echo "   BED    : $BED"
echo "   logs   : $LOGDIR/QCs_*<jobid>*.log"
echo

# ---- Stage 1: split the BED ------------------------------------------------
echo " - stage 1: 1_1a_chunk_bed.sh (split)"
SPLIT_JID=$(sbatch --parsable \
  --export=ALL,POOL=${POOL},OUTDIR=${OUTDIR},N=${N} \
  --job-name="QCs_${POOL}_split" \
  "$SCRIPTS/common/1_1a_chunk_bed.sh")
echo "   submitted: $SPLIT_JID"
echo

# ---- Stage 2: per-chunk array (depends on split) ---------------------------
echo " - stage 2: 1_1b_per_chunk.sh (array of $N chunks, depends on $SPLIT_JID)"
ARRAY_JID=$(sbatch --parsable \
  --dependency=afterok:${SPLIT_JID} \
  --export=ALL,POOL=${POOL},OUTDIR=${OUTDIR},N=${N} \
  --array=0-$((N-1))%${N} \
  --job-name="QCs_${POOL}_chunk" \
  "$SCRIPTS/common/1_1b_per_chunk.sh")
echo "   submitted: $ARRAY_JID  (array 0-$((N-1)))"
echo

# ---- Stage 3: merge + isCellv2 (depends on array all-ok) -------------------
echo " - stage 3: 1_1c_merge_and_qc.sh (merge, depends on $ARRAY_JID)"
MERGE_JID=$(sbatch --parsable \
  --dependency=afterok:${ARRAY_JID} \
  --export=ALL,POOL=${POOL},OUTDIR=${OUTDIR},N=${N}${KEEP_CHUNKS:+,KEEP_CHUNKS=${KEEP_CHUNKS}} \
  --job-name="QCs_${POOL}_merge" \
  "$SCRIPTS/common/1_1c_merge_and_qc.sh")
echo "   submitted: $MERGE_JID"
echo

# ---- Stage 4: 1_2 + 1_3 (depends on merge) ---------------------------------
echo " - stage 4: 1_2_1_3_run.sh (filter + meta-QC, depends on $MERGE_JID)"
DOWN_JID=$(sbatch --parsable \
  --dependency=afterok:${MERGE_JID} \
  --export=ALL,POOL=${POOL},OUTDIR=${OUTDIR},DEPTH=${DEPTH} \
  --job-name="QCs_${POOL}_down" \
  "$SCRIPTS/common/1_2_1_3_run.sh")
echo "   submitted: $DOWN_JID"
echo

echo " - chain submitted:"
echo "     $SPLIT_JID  (split)"
echo "  -> $ARRAY_JID  (chunk x $N)"
echo "  -> $MERGE_JID  (merge + isCellv2)"
echo "  -> $DOWN_JID   (1_2 + 1_3)"
echo
echo " - monitor:"
echo "     squeue -u \$USER -o \"%.10i %.20j %.8T %.10M %.6D %R\""
echo "     tail -f $LOGDIR/QCs_${POOL}_split-*.log"
echo "     tail -f $LOGDIR/QCs_${POOL}_chunk-${ARRAY_JID}_0.log"
echo "     tail -f $LOGDIR/QCs_${POOL}_merge-${MERGE_JID}.log"
echo "     tail -f $LOGDIR/QCs_${POOL}_down-${DOWN_JID}.log"
