#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=50G
#SBATCH --job-name=BWA_SMs
#SBATCH --partition=standard
#SBATCH --output=_logs/BWA_multigenome_%A_%a.log
#SBATCH --array=1-10
#
# Legacy multi-reference alignment: bwa mem -M of <lib>_R{1,3}.bc1.fastq.gz against one genome per array
# task, read from GENOME_LIST (line = array index; Genome_list_NAM_{1,2,3}.txt cover the 26 NAM genomes of
# the Root1 panel, Genome_list_SMs.txt the SM2 references). Output <lib>_<genome>_scifiATAC.raw.bam under
# ${PROJECT_ROOT}/3_Mapping/<lib>/, the input of 1_scifi_processBAM.ATAC.sh.
# Usage (from ${PROJECT_ROOT}/3_Mapping): GENOME_LIST=<list> sbatch --array=1-<n genomes> 1_1_align.scifi.ATAC.NAM.sh <lib>

: "${PROJECT_ROOT:?set PROJECT_ROOT to the analysis tree that holds 2_CleanReads/ and 3_Mapping/}"
: "${GENOMES_DIR:?set GENOMES_DIR to the directory holding Zea/NAN_Indexes/Index_Zm_<genome>_bwa}"
REPO_ROOT="${REPO_ROOT:-$(git rev-parse --show-toplevel 2>/dev/null || pwd)}"

########## Load Modules #########
# conda activate <env from environment.yml>   (samtools)
module load Bioinformatics bwa/0.7.17-mil4ns7   # HPC module for BWA 0.7.17; adapt to the local installation

########## Define Variables #########
THREADS=30
SAMPLE_ID="$1"  # library name, e.g. Root1_rep1

# Check if SAMPLE_ID is provided
if [ -z "$SAMPLE_ID" ]; then
    echo "Error: No sample name provided."
    echo "Usage: sbatch 1_1_align.scifi.ATAC.NAM.sh <SampleName>"
    exit 1
fi

# Paths to input files and the BWA index directory. The archived copy also carried a second,
# commented-out index directory (GenomesIndex/Mixed_index); the NAM index directory is the default here.
INPUT_DIR="${PROJECT_ROOT}/2_CleanReads"
BWA_INDEX_DIR="${BWA_INDEX_DIR:-${GENOMES_DIR}/Zea/NAN_Indexes}"
GENOME_LIST="${GENOME_LIST:-${REPO_ROOT}/workflows/02_mapping/legacy_sm2_root1/Genome_list_NAM_2.txt}"

# Input FASTQs carry the 10x barcode in the read name (.bc1); the scifi chain emits .bc1.bc2 files
FQ1_FILE="${INPUT_DIR}/${SAMPLE_ID}_R1.bc1.fastq.gz"
FQ2_FILE="${INPUT_DIR}/${SAMPLE_ID}_R3.bc1.fastq.gz"

echo "[$(date)] Verifying input FASTQ files: ${FQ1_FILE} and ${FQ2_FILE}"
if [ ! -f "${FQ1_FILE}" ] || [ ! -f "${FQ2_FILE}" ]; then
    echo "Error: FASTQ files for sample ${SAMPLE_ID} not found at ${INPUT_DIR}"
    exit 1
fi

########## Function: BWA Mapping #########
MappingStart() {
    # Accepts a genome target as its first argument
    GENOME_TARGET="$1"

    # Define local variables: BWA index and FASTQ files
    local REF="${BWA_INDEX_DIR}/Index_Zm_${GENOME_TARGET}_bwa"
    local FASTQ_R1="${INPUT_DIR}/${SAMPLE_ID}_R1.bc1.fastq.gz"
    local FASTQ_R2="${INPUT_DIR}/${SAMPLE_ID}_R3.bc1.fastq.gz"

    # Output directory for BWA alignment results
    local BWA_OUTPUT_DIR="${PROJECT_ROOT}/3_Mapping/${SAMPLE_ID}"

    # Create the output directory if it doesn't exist
    mkdir -p "${BWA_OUTPUT_DIR}"

    # Run BWA alignment
    echo "[$(date)] Starting BWA alignment for sample ${SAMPLE_ID} using genome target ${GENOME_TARGET}..."
    bwa mem -M -t 40 "${REF}" "${FASTQ_R1}" "${FASTQ_R2}" > "${BWA_OUTPUT_DIR}/${SAMPLE_ID}_${GENOME_TARGET}_scifiATAC.raw.sam"

    # Compress the SAM file to BAM format using samtools
    samtools view -@ "${THREADS}" -bS "${BWA_OUTPUT_DIR}/${SAMPLE_ID}_${GENOME_TARGET}_scifiATAC.raw.sam" > "${BWA_OUTPUT_DIR}/${SAMPLE_ID}_${GENOME_TARGET}_scifiATAC.raw.bam"
    echo "[$(date)] BWA alignment and SAM compression complete for sample ${SAMPLE_ID} with genome target ${GENOME_TARGET}."
    # Clean up
    rm "${BWA_OUTPUT_DIR}/${SAMPLE_ID}_${GENOME_TARGET}_scifiATAC.raw.sam"

}
export -f MappingStart

# Retrieve the genome target from the genome list using the SLURM array index.
GENOME_TARGET=$(sed -n "$((SLURM_ARRAY_TASK_ID))p" "${GENOME_LIST}")
if [ -z "${GENOME_TARGET}" ]; then
    echo "Error: No genome target found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} in ${GENOME_LIST}"
    exit 1
fi
echo "[$(date)] Genome target selected: ${GENOME_TARGET}"

# Execute the mapping function with the selected genome target
MappingStart "${GENOME_TARGET}"
