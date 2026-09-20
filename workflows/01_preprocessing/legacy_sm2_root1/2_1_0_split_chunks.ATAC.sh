#!/bin/bash
########## BATCH Lines for Resource Request ##########
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=20G
#SBATCH --job-name=FastqSplit
#SBATCH --partition=standard
#SBATCH --output=_logs/split_fastqs_ATACs_%j.log
#
# Legacy chain, step 3 (SM2 only): split <lib>_R{1,3}.bc1.bc2.fastq.gz into 60 chunk pairs with seqkit
# for the per-chunk demultiplexing array 2_1_1_fastq_sample_assign_and_fixTn5.ATAC.sh.
# Usage (from ${PROJECT_ROOT}/2_CleanReads): sbatch 2_1_0_split_chunks.ATAC.sh SM2_ATAC

: "${PROJECT_ROOT:?set PROJECT_ROOT to the analysis tree that holds 2_CleanReads/}"
# conda activate <env from environment.yml>   (seqkit)
cd "${PROJECT_ROOT}/2_CleanReads"
mkdir -p _logs

########## Define Variables #########
# Read SAMPLE from command-line argument
SAMPLE=$1

# Define input files
READ1=${SAMPLE}_R1.bc1.bc2.fastq.gz
READ2=${SAMPLE}_R3.bc1.bc2.fastq.gz

# Define split output directory
SPLIT_OUTPUT_DIR="${SAMPLE}_chunks"

########## Step 1: Split FASTQ Files #########
echo "[$(date)] Starting FASTQ splitting into 60 chunks using seqkit..."

# Create output directory if it doesn't exist
mkdir -p "${SPLIT_OUTPUT_DIR}"

echo "Splitting ${SAMPLE} into 60 chunks..."
seqkit split2 --by-part 60 -j 10 -O "${SPLIT_OUTPUT_DIR}" -1 "${READ1}" -2 "${READ2}"

echo "Splitting ${SAMPLE} completed."
