#!/usr/bin/env bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=80G
#SBATCH --job-name=umi_tools_R3
#SBATCH --partition=standard
#SBATCH --output=_logs/%x_%A_%a.log
#SBATCH --array=0-1
#
# Legacy chain, step 1 (R3): move the 16 bp 10x cell barcode from the index read (R2) into the read
# names of R3 with umi_tools extract; the input is split 20 ways with seqkit and processed with GNU parallel.
# Libraries: SM2_ATAC (raw FASTQs under 1_RawData/SMs/) and Root1_rep1, one per array task.
# Output: <name>_R3.bc1.fastq.gz in ${PROJECT_ROOT}/2_CleanReads. Companion: 1_1_UMItools.ATAC.R1_parallel.sh.
# Usage (from ${PROJECT_ROOT}/2_CleanReads): sbatch 1_1_UMItools.ATAC.R3_parallel.sh

: "${PROJECT_ROOT:?set PROJECT_ROOT to the analysis tree that holds 1_RawData/ and 2_CleanReads/}"
# conda activate <env from environment.yml>   (umi_tools, seqkit, GNU parallel, pigz)
cd "${PROJECT_ROOT}/2_CleanReads"
mkdir -p _logs

# library prefixes of the raw FASTQs <name>_R{1,2,3}.fastq.gz, one array task per library
SAMPLES=(SM2_ATAC Root1_rep1)

# pick sample
name="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
echo "[INFO] Sample: $name (Job $SLURM_JOB_ID)"

# input paths
case "$name" in
  SM2_ATAC) INPUT_DIR="${PROJECT_ROOT}/1_RawData/SMs" ;;
  *)        INPUT_DIR="${PROJECT_ROOT}/1_RawData" ;;
esac

R3_IN="${INPUT_DIR}/${name}_R3.fastq.gz"
R2_IN="${INPUT_DIR}/${name}_R2.fastq.gz"

# output dirs and files
export OUT_ROOT="${name}_bc1_R3"
export CHUNK_DIR="${OUT_ROOT}/chunks"
export TEMP_OUT="${OUT_ROOT}/chunks_out"
export FINAL_R3="${name}_R3.bc1.fastq.gz"

mkdir -p "$CHUNK_DIR" "$TEMP_OUT"
export TMPDIR=./tem_file
mkdir -p "$TMPDIR"

# verify inputs
for f in "$R3_IN" "$R2_IN"; do
  [[ -f "$f" ]] || { echo "[ERROR] Missing $f"; exit 1; }
done

# split into 20 parts on all threads
echo "[INFO] Splitting into 20 parts..."
seqkit split2 \
 --by-part 20 \
 -j "$SLURM_CPUS_PER_TASK" \
 -O "$CHUNK_DIR" \
 -1 "$R3_IN" \
 -2 "$R2_IN"

echo "[INFO] Processing chunks with umi_tools extract..."
process_chunk(){
  local r3_chunk="$1"
  local r2_chunk="${r3_chunk//_R3./_R2.}"
  local base="$(basename "$r3_chunk" .fastq.gz)"
  local out="${TEMP_OUT}/${base}.bc1.fastq"

  umi_tools extract \
    --bc-pattern=NNNNNNNNNNNNNNNN \
    --stdin="$r2_chunk" \
    --read2-in="$r3_chunk" \
    --read2-out="$out"
}
export -f process_chunk

find "$CHUNK_DIR" -name "*_R3*.fastq.gz" | \
  parallel --tmpdir "$TMPDIR" --compress --halt soon,fail=1 -j "$SLURM_CPUS_PER_TASK" process_chunk {}

# concatenate & compress
echo "[INFO] Concatenating & compressing R3-only..."
cat "$TEMP_OUT"/*.bc1.fastq > "${OUT_ROOT}/all_R3.bc1.fastq"
pigz --processes "$SLURM_CPUS_PER_TASK" --force -c "${OUT_ROOT}/all_R3.bc1.fastq" > "$FINAL_R3"

echo "[DONE] Final file: $FINAL_R3"

# cleanup of the chunk intermediates
rm -rf "$CHUNK_DIR" "$TEMP_OUT" "${OUT_ROOT}/all_R3.bc1.fastq"
