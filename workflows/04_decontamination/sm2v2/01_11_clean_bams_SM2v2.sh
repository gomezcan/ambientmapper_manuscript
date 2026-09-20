#!/usr/bin/env bash
#SBATCH --job-name=sm2v2_clean
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --array=0-3
#SBATCH --output=_logs/01_11_clean_bams_SM2v2_%A_%a.log
#
# 01_11_clean_bams_SM2v2.sh — post-decontam BAM cleaning for SM2v2.
#
# Runs BOTH decontam passes (with_design = WD, without_design = ND) since each
# produces a different reads_to_drop set and both feed the downstream QC and
# embedding stage (05_qc_and_embedding, WD and ND arms). Uses the per-chrom
# parallel sidecar workflows/02_mapping/tn5bed/00_bam_to_tn5bed_parallel.py
# (shipped by stage 02) for the tn5-bed step.
#
# Input BAMs (created by 01_00_merge_SM2v2_inputs.sh):
#   3_Mapping/ambientmapper_input/SM2_B73v5_scifiATAC.mq10.BC.rmdup.mm.bam   (~3.2 GB)
#   3_Mapping/ambientmapper_input/SM2_TAIR10_scifiATAC.mq10.BC.rmdup.mm.bam  (~525 MB)
# (BAM stems are SM2_*, not SM2v2_* — the inputs are shared between v1
# and v2; only the AmbientMapper processing differs.)
#
# Output stems mirror inputs:
#   SM2v2/clean_bams_alpha05_C0_wd/SM2_{B73v5,TAIR10}_*.Clean.bam{,.bai,.tn5.bed.gz}
#   SM2v2/clean_bams_alpha05_C0_nd/SM2_{B73v5,TAIR10}_*.Clean.bam{,.bai,.tn5.bed.gz}
#
# Array task layout:
#   0: B73v5  / with_design     2: B73v5  / without_design
#   1: TAIR10 / with_design     3: TAIR10 / without_design
#
# Resource sizing: BAMs are small (~3.7 GB combined per pass), so 2 h
# wall + 16 GB mem + 8 cpu (6 chrom workers + 2 pigz) is generous.
# =============================================================================

set -euo pipefail

REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"   # this repository (sourced helpers, sibling scripts)
PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"   # root of the data tree (3_Mapping/, 5_AmbientDetection/, ...)

cd "${PROJECT_ROOT}/5_AmbientDetection"

# conda activate <env from environment.yml>

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export GOTO_NUM_THREADS=1

SAMPLE="SM2v2"
INPUT_DIR="${PROJECT_ROOT}/3_Mapping/ambientmapper_input"

genomes=(B73v5  TAIR10 B73v5         TAIR10)
passes=(wd     wd     nd            nd)

idx="${SLURM_ARRAY_TASK_ID:?This script must be submitted as an array job (sbatch --array=0-3 ...)}"
genome="${genomes[$idx]}"
pass="${passes[$idx]}"

case "$pass" in
  wd) decontam_subdir="decontam_with_design_alpha05_v2";    out_subdir="clean_bams_alpha05_C0_wd" ;;
  nd) decontam_subdir="decontam_without_design_alpha05_v2"; out_subdir="clean_bams_alpha05_C0_nd" ;;
  *)  echo "ERROR: unknown pass=$pass"; exit 1 ;;
esac

bam_in="${INPUT_DIR}/SM2_${genome}_scifiATAC.mq10.BC.rmdup.mm.bam"
drop_file="${SAMPLE}/${decontam_subdir}/${SAMPLE}_reads_to_drop.tsv.gz"
out_dir="${SAMPLE}/${out_subdir}"

[[ -f "$bam_in"       ]] || { echo "Missing input BAM: $bam_in"; exit 1; }
[[ -f "${bam_in}.bai" ]] || { echo "Missing input BAI: ${bam_in}.bai"; exit 1; }
[[ -f "$drop_file"    ]] || { echo "Missing reads_to_drop: $drop_file"; exit 1; }

mkdir -p "$out_dir"

threads="${SLURM_CPUS_PER_TASK:-8}"
pigz_threads=2
worker_threads=$(( threads - pigz_threads ))
(( worker_threads >= 1 )) || worker_threads=1

bam_base="$(basename "$bam_in" .bam)"
bam_clean="${out_dir}/${bam_base}.Clean.bam"
out_bed_gz="${bam_clean%.bam}.tn5.bed.gz"
scratch="${TMPDIR:-/tmp}"

echo "============================================================"
echo "[$(date)] SM2v2 BAM-clean task ${idx}: genome=${genome} pass=${pass}"
echo "  bam_in     : ${bam_in} ($(du -h "$bam_in" | cut -f1))"
echo "  drop_file  : ${drop_file}"
echo "  out_dir    : ${out_dir}"
echo "  bam_clean  : ${bam_clean}"
echo "  out_bed_gz : ${out_bed_gz}"
echo "  threads    : ${threads} (workers=${worker_threads}, pigz=${pigz_threads})"
echo "============================================================"

# 1. clean-bams (drops reads listed in $drop_file from $bam_in, indexes output)
if [[ -s "$bam_clean" && -s "${bam_clean}.bai" ]]; then
  echo "[clean-bams] SKIP: ${bam_clean} already exists with index."
else
  ambientmapper clean-bams \
    --reads-to-drop "$drop_file" \
    --bam "$bam_in" \
    --out-dir "$out_dir" \
    --out-suffix .Clean.bam
  [[ -s "$bam_clean" ]] || { echo "ERROR: Clean BAM missing: $bam_clean"; exit 1; }
  echo "[clean-bams] -> ${bam_clean}"
fi

# 2. tn5.bed.gz via per-chrom parallel sidecar
python "${REPO_ROOT}/workflows/02_mapping/tn5bed/00_bam_to_tn5bed_parallel.py" \
  --bam "$bam_clean" \
  --out "$out_bed_gz" \
  --threads "$worker_threads" \
  --pigz-threads "$pigz_threads" \
  --scratch-dir "$scratch" \
  --skip-existing

[[ -s "$out_bed_gz" ]] || { echo "ERROR: tn5.bed.gz missing or empty: $out_bed_gz"; exit 1; }

echo "============================================================"
echo "[$(date)] task ${idx} DONE for genome=${genome} pass=${pass}"
ls -la "$out_dir"/"${bam_base}".*
echo "============================================================"
