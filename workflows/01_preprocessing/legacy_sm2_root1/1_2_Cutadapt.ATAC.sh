#!/bin/bash
########## BATCH Lines for Resource Request ##########
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --job-name=cutadap_tools_
#SBATCH --partition=standard
#SBATCH --output=_logs/%x-%A_%a.log
#SBATCH --array=0
#
# Legacy chain, step 2 (SM2 only; Root1 is a 10x scATAC library without Tn5 barcodes and skips this step):
# cutadapt trims the 5 bp Tn5 well barcodes from R1 and R3 (error rate 0.2 against the mosaic end
# AGATGTGTATAAGAGACAG) and appends them to the read name: <name>_R{1,3}.bc1.fastq.gz -> <name>_R{1,3}.bc1.bc2.fastq.gz
# Usage (from ${PROJECT_ROOT}/2_CleanReads): sbatch 1_2_Cutadapt.ATAC.sh

: "${PROJECT_ROOT:?set PROJECT_ROOT to the analysis tree that holds 2_CleanReads/}"
# conda activate <env from environment.yml>   (cutadapt)
cd "${PROJECT_ROOT}/2_CleanReads"
mkdir -p _logs

########## Input File #########
# library prefixes (one per array task); SM2 = SM2_ATAC
SAMPLES=(SM2_ATAC)

# Get the sample name based on the array task ID
name=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

if [ -z "$name" ]; then
    echo "Error: No sample found for SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID"
    exit 1
fi

echo "Processing sample: $name"

########## Define a Single Function #########

process_sample() {
    local name=$1

    # Define file paths
    local tenxBC=${name}_R2.fastq.gz
    local R1=${name}_R1.fastq.gz
    local R2=${name}_R3.fastq.gz
    local R1tenx=${name}_R1.bc1.fastq.gz
    local R2tenx=${name}_R3.bc1.fastq.gz
    local R1tenxtn5=${name}_R1.bc1.bc2.fastq.gz
    local R2tenxtn5=${name}_R3.bc1.bc2.fastq.gz

    # Append Tn5 bc from R1 and R2 to read-name
    cutadapt -e 0.2 \
        --pair-filter=any \
        -j 4 \
        --rename='{id}_{r1.cut_prefix}_{r2.cut_prefix} {comment}' \
        -u 5 \
        -U 5 \
        -g AGATGTGTATAAGAGACAG \
        -G AGATGTGTATAAGAGACAG \
        -o $R1tenxtn5 \
        -p $R2tenxtn5 \
        $R1tenx $R2tenx
    if [ $? -ne 0 ]; then
        echo "Error: cutadapt failed for sample $name"
        exit 1
    fi

    # Log completion
    echo "Job completed successfully for sample $name on $(date)"
}

########## Call the Function for the Current Sample #########
process_sample "$name"
