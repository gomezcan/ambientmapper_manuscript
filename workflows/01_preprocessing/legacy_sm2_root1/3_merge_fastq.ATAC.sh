#!/bin/bash
########## BATCH Lines for Resource Request ##########
#SBATCH --time=0:40:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=5G
#SBATCH --job-name=MergeFastq
#SBATCH --partition=standard
#SBATCH --output=_logs/combine_results_%j.log
#
# Legacy chain, step 5 (SM2 only): concatenate the corrected chunks per sample into
# <sample>_R{1,3}.bc1.bc2.fastq.gz (SM2_B73, SM2_At, SM2_MUDR), the input of the alignment step.
# Usage (from ${PROJECT_ROOT}/2_CleanReads): sbatch 3_merge_fastq.ATAC.sh SM2_ATAC

: "${PROJECT_ROOT:?set PROJECT_ROOT to the analysis tree that holds 2_CleanReads/}"
# conda activate <env from environment.yml>   (GNU parallel)
cd "${PROJECT_ROOT}/2_CleanReads"
mkdir -p _logs

########## Define Variables #########
# Check that SAMPLE is provided as a command-line argument.
if [ $# -lt 1 ]; then
    echo "Usage: $0 SAMPLE"
    exit 1
fi

export SAMPLE="$1"

# Set directories.
INPUT_DIC="${SAMPLE}_chunks_Corrected"

export INPUT_DIC
export OUTPUT_DIC

########## Step 1: Get Pools to Merge #########
echo "[$(date)] Getting list of all pools for ${SAMPLE}..."

# List files in the input directory, extract the pool (sample) name from the file name,
# and write the unique pool names to the Pools file.
ls "${INPUT_DIC}/${SAMPLE}"_R* 2>/dev/null | \
    cut -d '/' -f2 | \
    cut -d'.' -f4 | \
    awk -F'_' '{print $(NF-2)"_"$(NF-1)"_"$NF}' | \
    sort -u > "Pools_${SAMPLE}.txt"

if [ ! -s "Pools_${SAMPLE}.txt" ]; then
    echo "[$(date)] Error: No pools found in ${INPUT_DIC}."
    exit 1
fi

########## Step 2: Merge FASTQ Files for Each Pool #########
Merge_fq(){
    local line="$1"
    # Remove any potential whitespace.
    pool="${line//[[:space:]]/}"

    # Define output filenames.
    OUTPUT_R1="${pool}_R1.bc1.bc2.fastq.gz"
    OUTPUT_R3="${pool}_R3.bc1.bc2.fastq.gz"

    # Combine R1 FASTQ files for the sample pool.
    echo "[$(date)] Combining R1 FASTQ files for pool ${pool}..."
    cat "${INPUT_DIC}/${SAMPLE}"_R1.bc1.bc2.part_*_"${pool}".fastq.gz > "${OUTPUT_R1}"
    echo "[$(date)] Created ${OUTPUT_R1}"

    # Combine R3 FASTQ files for the sample pool.
    echo "[$(date)] Combining R3 FASTQ files for pool ${pool}..."
    cat "${INPUT_DIC}/${SAMPLE}"_R3.bc1.bc2.part_*_"${pool}".fastq.gz > "${OUTPUT_R3}"
    echo "[$(date)] Created ${OUTPUT_R3}"
}

export -f Merge_fq

echo "[$(date)] Combining FASTQ files per sample..."
# Count the number of pools.
num_jobs=$(wc -l < "Pools_${SAMPLE}.txt")
# Process each pool concurrently using GNU parallel.
parallel -j "$num_jobs" Merge_fq :::: "Pools_${SAMPLE}.txt"

echo "[$(date)] Combining results completed."
