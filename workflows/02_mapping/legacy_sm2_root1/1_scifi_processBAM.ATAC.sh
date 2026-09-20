#!/usr/bin/env bash
########## BATCH Lines for Resource Request ##########
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100G
#SBATCH --job-name=bamQC_clean
#SBATCH --partition=standard
#SBATCH --output=_logs/Clean_bams_%A_%a.log
#SBATCH --array=1-3
#
# Legacy BAM cleaning chain (SM2, Root1): one raw BAM (<lib>_<genome>_scifiATAC.raw.bam) per array task.
#   sort -> BC tag from the read name + MAPQ >= 10 proper pairs -> barcode counting (more than 5 reads)
#   -> barcode correction (10x segment within two substitutions to the closest whitelist entry, Tn5 halves
#   within one) -> BC tags rewritten -> Picard MarkDuplicates (BARCODE_TAG=BC) -> multi-mapping filter and
#   library tag (1_5_scifi_fixBC.pl) -> Tn5 insertion BED (1_6_scifi_makeTn5bed.py), sorted and unique.
# Output: ${PROJECT_ROOT}/3_Mapping/<lib>/<base>.mq10.BC.rmdup.mm.bam (+ intermediates, metrics, barcode counts)
#         ${PROJECT_ROOT}/3_Mapping/<lib>/bed/<base>.mq10.tn5.bed.gz
# Usage (from ${PROJECT_ROOT}/3_Mapping): sbatch --array=1-<n raw BAMs> 1_scifi_processBAM.ATAC.sh <lib>
# 1_3_scifi_correctBCs.10x.v2.pl reads tn5_bcs.txt and 737K-cratac-v1.txt (10x Genomics scATAC barcode
# whitelist, not shipped) from the working directory, which this script sets to its own directory.

: "${PROJECT_ROOT:?set PROJECT_ROOT to the analysis tree that holds 3_Mapping/}"
REPO_ROOT="${REPO_ROOT:-$(git rev-parse --show-toplevel 2>/dev/null || pwd)}"
SCRIPTS="${REPO_ROOT}/workflows/02_mapping/legacy_sm2_root1"

###################################
#######       Modules     #########
###################################
# conda activate <env from environment.yml>   (samtools, picard, GNU parallel, perl with String::Approx
#                                              and List::MoreUtils, python with pysam, pigz)

###################################
#######  Input & Output Dirs  #####
###################################

SAMPLE=$1

INPUT_DIC="${PROJECT_ROOT}/3_Mapping/${SAMPLE}"
OUT_DIC="${PROJECT_ROOT}/3_Mapping/${SAMPLE}"
OUT_DICBED="${PROJECT_ROOT}/3_Mapping/${SAMPLE}/bed"

# collect all raw BAMs into an array
mapfile -t bam_files < <(ls "${INPUT_DIC}"/*_scifiATAC.raw.bam | sort)
if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > ${#bam_files[@]} )); then
  echo "ERROR: SLURM_ARRAY_TASK_ID out of range" >&2
  exit 1
fi

# select this job's BAM
bam_file="${bam_files[$SLURM_ARRAY_TASK_ID-1]}"
base="$(basename "$bam_file" .raw.bam)"

echo "[$(date)] Processing sample: $base"

# make sample-specific dirs
mkdir -p "$OUT_DIC"
mkdir -p "$OUT_DICBED"

# the Perl barcode corrector reads its whitelists from the working directory
cd "${SCRIPTS}"

# helper function
doCall(){
  local base="$1"
  local threads="$SLURM_CPUS_PER_TASK"

  # 1) sort
  samtools sort -@ "$threads" -o "$OUT_DIC/$base.rawSort.bam"  "$bam_file"

  # 2) tag BC from the read name & filter MAPQ >= 10, keep proper pairs (-f 3)
  perl "$SCRIPTS/1_1_scifi_modufy_BC_flag.pl" "$OUT_DIC/$base.rawSort.bam" \
  | samtools view -@ "$threads" -hb -q 10 -f 3 - \
    > "$OUT_DIC/$base.mq10.bam"
  echo ".. done BC tagging & MQ filter .."

  # 3) count barcodes, keep those seen on more than 5 reads
  perl "$SCRIPTS/1_2_scifi_countBCs.BAM.pl" "$OUT_DIC/$base.mq10.bam" \
  | awk '$2>5' \
    > "$OUT_DIC/$base.mq10.barcodes.txt"
  echo ".. done barcode counting .."

  # 4) barcode correction
  cat "$OUT_DIC/$base.mq10.barcodes.txt" \
  | parallel --pipe -k -j "$threads" -N 1000 \
      perl "$SCRIPTS/1_3_scifi_correctBCs.10x.v2.pl" \
    > "$OUT_DIC/$base.mq10.barcodes.corrected.txt"
  echo ".. done barcode correction .."

  # 5) update BAM with corrected BC tag
  perl "$SCRIPTS/1_4_scifi_correctBAM.pl" \
      "$OUT_DIC/$base.mq10.barcodes.corrected.txt" \
      "$OUT_DIC/$base.mq10.bam" \
  | samtools view -@ "$threads" -bhS -f 3 - \
    > "$OUT_DIC/$base.mq10.BC.bam"
  echo ".. done updating BC tags .."

  # 6) remove duplicates
  echo "removing duplicates for $base ..."
  picard MarkDuplicates \
    I="$OUT_DIC/$base.mq10.BC.bam" \
    O="$OUT_DIC/$base.mq10.BC.rmdup.bam" \
    METRICS_FILE="$OUT_DIC/$base.metrics" \
    REMOVE_DUPLICATES=true \
    BARCODE_TAG=BC \
    ASSUME_SORT_ORDER=coordinate \
    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000
  echo ".. done deduplication .."

  # 7) fix multi-mapping & BC
  perl "$SCRIPTS/1_5_scifi_fixBC.pl" "$threads" \
  	"$OUT_DIC/${base}.mq10.BC.rmdup.bam" \
  	"$OUT_DIC/${base}.mq10.BC.rmdup.mm.bam" \
  	"$OUT_DIC/${base}_bc_counts.txt" \
  	"$base"

  echo ".. done fixing BC & multi-mapping .."

  # 8) make Tn5 BED
  samtools index -@ "$threads" "$OUT_DIC/$base.mq10.BC.rmdup.mm.bam"
  python "$SCRIPTS/1_6_scifi_makeTn5bed.py" \
    "$OUT_DIC/$base.mq10.BC.rmdup.mm.bam" \
  | sort -k1,1 -k2,2n \
  | uniq \
    > "$OUT_DICBED/$base.mq10.tn5.bed"
  pigz -p "$threads" "$OUT_DICBED/$base.mq10.tn5.bed"
  echo ".. done BED file .."

  # read counts before and after the multi-mapping filter
  samtools view -@ "$threads" -c "$OUT_DIC/$base.mq10.BC.rmdup.bam" > "$OUT_DIC/$base.mq10.BC.rmdup.proper_pairs.txt"
  samtools view -@ "$threads" -c "$OUT_DIC/$base.mq10.BC.rmdup.mm.bam" > "$OUT_DIC/$base.mq10.BC.rmdup.mm.proper_pairs.txt"
}

# run it
doCall "$base"
