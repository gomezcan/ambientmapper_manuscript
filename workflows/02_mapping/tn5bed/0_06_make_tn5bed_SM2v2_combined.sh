#!/usr/bin/env bash
#SBATCH --job-name=SM2v2_tn5bed_combined
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --array=0-1
#SBATCH --output=_logs/0_06_tn5bed_SM2v2_combined_%A_%a.log
#
# 0_06_make_tn5bed_SM2v2_combined.sh
#
# PostClean Tn5 BEDs of the cleaned concatenated-reference BAMs Clean.SM2v2_{At,B73}_ZmATcombined_*.bam
# (written by workflows/04_decontamination/combined_genome/0_05_clean_bams_SM2v2_combined.sh) for the
# combined-genome Socrates track. Output naming places them alongside the legacy SM2 BEDs in
# 6_socrates/_data/_BED_files/.
#
# PreClean BEDs are NOT regenerated here: the raw ZmATcombined BAM is identical between the legacy SM2
# and the SM2v2 cleaning runs (only the drop list differs), so the existing
#   _data/_BED_files/SM2_{At,B73}_ZmATcombined_*.tn5.bed.gz
# serve as the PreClean input.
#
# Uses the per-chromosome parallel extractor 00_bam_to_tn5bed_parallel.py (the same one the other datasets use).
#
# Array task layout:
#   0: Clean.SM2v2_At_ZmATcombined_*
#   1: Clean.SM2v2_B73_ZmATcombined_*
# Usage (from ${PROJECT_ROOT}/6_socrates): sbatch 0_06_make_tn5bed_SM2v2_combined.sh

set -euo pipefail
: "${PROJECT_ROOT:?set PROJECT_ROOT to the analysis tree that holds 6_socrates/}"
REPO_ROOT="${REPO_ROOT:-$(git rev-parse --show-toplevel 2>/dev/null || pwd)}"

BASE="${PROJECT_ROOT}/6_socrates"
cd "${BASE}"
mkdir -p _logs

# conda activate <env from environment.yml>   (python with pysam, pigz)

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export GOTO_NUM_THREADS=1

CLEAN_BAM_DIR="${BASE}/_data/SM2v2_clean_bams_combined"
OUT_DIR="${BASE}/_data/_BED_files"
mkdir -p "${OUT_DIR}"

idx="${SLURM_ARRAY_TASK_ID:?Array index required (sbatch --array=0-1)}"
species_arr=(At B73)
species="${species_arr[$idx]}"

clean_bam="${CLEAN_BAM_DIR}/Clean.SM2v2_${species}_ZmATcombined_scifiATAC.mq10.BC.rmdup.mm.bam"
out_bed_gz="${OUT_DIR}/Clean.SM2v2_${species}_ZmATcombined_scifiATAC.mq10.BC.rmdup.mm.tn5.bed.gz"

[[ -f "${clean_bam}"      ]] || { echo "ERROR: missing clean BAM ${clean_bam}  (run 0_05 first)"; exit 1; }
[[ -f "${clean_bam}.bai"  ]] || { echo "ERROR: missing index ${clean_bam}.bai"; exit 1; }

threads="${SLURM_CPUS_PER_TASK:-12}"
pigz_threads=2
worker_threads=$(( threads - pigz_threads ))
(( worker_threads >= 1 )) || worker_threads=1
scratch="${TMPDIR:-/tmp}"

SIDECAR="${REPO_ROOT}/workflows/02_mapping/tn5bed/00_bam_to_tn5bed_parallel.py"
[[ -f "${SIDECAR}" ]] || { echo "ERROR: missing extractor ${SIDECAR}"; exit 1; }

echo "============================================================"
echo "[$(date)] SM2v2 tn5-bed task ${idx} (species=${species})"
echo "  clean_bam : ${clean_bam} ($(du -h "${clean_bam}" | cut -f1))"
echo "  out_bed   : ${out_bed_gz}"
echo "  threads   : ${threads} (workers=${worker_threads}, pigz=${pigz_threads})"
echo "  scratch   : ${scratch}"
echo "============================================================"

python "${SIDECAR}" \
  --bam "${clean_bam}" \
  --out "${out_bed_gz}" \
  --threads "${worker_threads}" \
  --pigz-threads "${pigz_threads}" \
  --scratch-dir "${scratch}" \
  --skip-existing

[[ -s "${out_bed_gz}" ]] || { echo "ERROR: tn5 bed missing/empty: ${out_bed_gz}"; exit 1; }

echo "[$(date)] DONE species=${species}"
ls -la "${out_bed_gz}"
