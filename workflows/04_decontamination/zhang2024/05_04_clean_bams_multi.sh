#!/usr/bin/env bash
#SBATCH --job-name=multi_clean
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --array=0-6
#SBATCH --output=_logs/05_04_clean_bams_multi_%A_%a.log
#
# 05_04_clean_bams_multi.sh — post-decontam BAM cleaning for multiGenotypes_rep1.
#
# One genome BAM per array task (7 BAMs, ~80 GB each = ~560 GB total). Reads
# to drop come from the C0/without_design decontam (05_03). The Clean BAMs are
# the "clean" arm of the WASP and purity chain (03b_variant_based_comparison).
#
# Array task layout (genome token = file infix):
#   0: B73v5      1: B97       2: Ky21      3: M162W
#   4: Mo18W      5: Oh7B      6: Tzi8
#
# Per task:
#   1. `ambientmapper clean-bams` drops contaminated reads from the per-genome
#      BAM using multiGenotypes_rep1_reads_to_drop.tsv.gz. Indexed BAM written
#      by default.
#   2. Tn5 BED: pysam-based site extraction → sort → uniq → pigz, for downstream
#      Socrates peak calling.
#
# Output dir: multiGenotypes_rep1/clean_bams_alpha05_C0_nd/
#
# Memory bumped to 48G vs 04_07's 32G — drop set is larger (73K chunks × deeper
# contamination) and BAMs are ~1.6× larger than B73Mo17_rep1.
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

INPUT_DIR="${PROJECT_ROOT}/3_Mapping/ambientmapper_input"
BAM_TO_BED="${PROJECT_ROOT}/3_Mapping/_archive/1_6_scifi_makeTn5bed.py"

sample="multiGenotypes_rep1"
genomes=(B73v5 B97 Ky21 M162W Mo18W Oh7B Tzi8)

idx="${SLURM_ARRAY_TASK_ID:?This script must be submitted as an array job (sbatch --array=0-6 ...)}"
genome="${genomes[$idx]}"

bam_in="${INPUT_DIR}/${sample}_${genome}_scifiATAC.mq10.BC.rmdup.mm.bam"
drop_file="${sample}/decontam_without_design_alpha05_C0/${sample}_reads_to_drop.tsv.gz"
out_dir="${sample}/clean_bams_alpha05_C0_nd"

[[ -f "$bam_in"     ]] || { echo "Missing BAM: $bam_in"; exit 1; }
[[ -f "$drop_file"  ]] || { echo "Missing reads_to_drop: $drop_file"; exit 1; }
[[ -f "$BAM_TO_BED" ]] || { echo "Missing Tn5 BED script: $BAM_TO_BED"; exit 1; }

mkdir -p "$out_dir"

threads="${SLURM_CPUS_PER_TASK:-5}"
bam_base="$(basename "$bam_in" .bam)"
bam_clean="${out_dir}/${bam_base}.Clean.bam"

echo "============================================================"
echo "[$(date)] multi BAM-clean task ${idx}: ${sample} / ${genome}"
echo "  bam_in    : ${bam_in} ($(du -h "$bam_in" | cut -f1))"
echo "  drop_file : ${drop_file}"
echo "  out_dir   : ${out_dir}"
echo "  threads   : ${threads}"
echo "============================================================"

ambientmapper clean-bams \
  --reads-to-drop "$drop_file" \
  --bam "$bam_in" \
  --out-dir "$out_dir" \
  --out-suffix .Clean.bam

[[ -f "$bam_clean" ]] || { echo "ERROR: Clean BAM missing: $bam_clean"; exit 1; }
echo "[clean-bams] -> ${bam_clean}"

tn5_bed="${bam_clean%.bam}.tn5.bed"
python "$BAM_TO_BED" "$bam_clean" \
  | sort -k1,1 -k2,2n \
  | uniq > "$tn5_bed"
pigz -p "$threads" "$tn5_bed"
echo "[tn5] -> ${tn5_bed}.gz"

echo "============================================================"
echo "[$(date)] task ${idx} DONE for ${sample} / ${genome}"
ls -la "${out_dir}/${bam_base}".*
echo "============================================================"
