#!/bin/bash
#SBATCH --job-name=macs3_SM2v2
#SBATCH --partition=standard
#SBATCH --array=0-3
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=_logs/macs3_SM2v2_%A_%a.log
#SBATCH --error=_logs/macs3_SM2v2_%A_%a.log
#
# 0_07_macs3_SM2v2.sh
#
# MACS3 ACR pre-pass for the SM2v2 chunked Socrates pipeline. One peak set per
# sample (raw At, raw B73, Clean.SM2v2_At, Clean.SM2v2_B73) so 1_1b_per_chunk
# can load pre-computed peaks instead of calling MACS2 inline.
#
# Parameters match the inline MACS2 invocation that the monolithic 1_1 used
# (1_1_QC_scifiATAC_data.R:429): genomesize=1.6e9+0.089e9, shift=-75,
# extsize=150, fdr=0.1, --keep-dup all. ENCODE Tn5 single-end recipe (-f BAM).
#
# Outputs (per sample, in _data/_PeakFiles/<sample>/):
#   <sample>_peaks.narrowPeak
#   <sample>_summits.bed
#   <sample>_peaks.xls
#
# Submission:
#   sbatch 0_scripts/setup/0_07_macs3_SM2v2.sh
#
# Array task layout:
#   0: SM2_At              (raw,  PreClean)
#   1: SM2_B73             (raw,  PreClean)
#   2: Clean.SM2v2_At      (Clean BAM from 0_05, PostClean)
#   3: Clean.SM2v2_B73     (Clean BAM from 0_05, PostClean)

set -euo pipefail

# Paths are anchored on PROJECT_ROOT ($BASH_SOURCE self-location does not survive sbatch).
PROJECT_ROOT="${PROJECT_ROOT:-.}"
BASE="${PROJECT_ROOT}/6_socrates"
SCRIPTS="${SCRIPTS:-${BASE}/0_scripts}"
cd "${BASE}"

# conda activate ambientmapper-manuscript   (environment.yml at the repo root)

# --- sample table -----------------------------------------------------------
samples=(SM2_At SM2_B73 Clean.SM2v2_At Clean.SM2v2_B73)
# Raw BAMs (tasks 0-1): the canonical combined-genome mapping lives in 3_Mapping/.
bams=(
  "${PROJECT_ROOT}/3_Mapping/SM2_At/SM2_At_ZmATcombined_scifiATAC.mq10.BC.rmdup.mm.bam"
  "${PROJECT_ROOT}/3_Mapping/SM2_B73/SM2_B73_ZmATcombined_scifiATAC.mq10.BC.rmdup.mm.bam"
  "${BASE}/_data/SM2v2_clean_bams_combined/Clean.SM2v2_At_ZmATcombined_scifiATAC.mq10.BC.rmdup.mm.bam"
  "${BASE}/_data/SM2v2_clean_bams_combined/Clean.SM2v2_B73_ZmATcombined_scifiATAC.mq10.BC.rmdup.mm.bam"
)

idx="${SLURM_ARRAY_TASK_ID:?Array index required (sbatch --array=0-3)}"
sample="${samples[$idx]}"
bam="${bams[$idx]}"

OUT_DIR="${BASE}/_data/_PeakFiles/${sample}"
mkdir -p "${OUT_DIR}"

peaks="${OUT_DIR}/${sample}_peaks.narrowPeak"
summits="${OUT_DIR}/${sample}_summits.bed"

echo "============================================================"
echo "[$(date)] MACS3 task ${idx}: ${sample}"
echo "  bam     : ${bam}"
echo "  outdir  : ${OUT_DIR}"
echo "============================================================"

[[ -f "${bam}"     ]] || { echo "ERROR: missing BAM ${bam}"; exit 1; }
[[ -f "${bam}.bai" ]] || { echo "ERROR: missing index ${bam}.bai"; exit 1; }

# Resume: skip if outputs already exist.
if [[ -s "${peaks}" && -s "${summits}" ]]; then
  echo "[$(date)] outputs already present, skipping"
  echo "  peaks  : $(wc -l < "${peaks}") lines"
  echo "  summits: $(wc -l < "${summits}") lines"
  exit 0
fi

# ENCODE Tn5 single-end recipe; genome size = 1.6e9 (Zm B73v5 mappable)
# + 0.089e9 (At TAIR10 mappable) for the combined-genome BAM = 1.689e9.
# -q 0.1 matches the existing Socrates inline MACS2 fdr=0.1.
echo "[$(date)] running MACS3..."
macs3 callpeak \
  -t "${bam}" \
  -f BAM \
  -g 1.689e9 \
  --nomodel \
  --shift -75 --extsize 150 \
  --keep-dup all \
  -q 0.1 \
  --call-summits \
  -n "${sample}" \
  --outdir "${OUT_DIR}"

echo "[$(date)] done"
echo "  peaks  : $(wc -l < "${peaks}") lines"
echo "  summits: $(wc -l < "${summits}") lines"
ls -lh "${OUT_DIR}"/${sample}*
