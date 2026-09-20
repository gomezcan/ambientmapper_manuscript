#!/bin/bash
########## BATCH Lines for Resource Request ##########
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=60G
#SBATCH --job-name=BWA_at
#SBATCH --partition=standard
#SBATCH --output=_logs/%x-%j.log
#
# Legacy single-reference alignment wrapper (Arabidopsis TAIR10 index): bwa mem -M of the
# <lib>_R{1,3}.bc1.fastq.gz pair (10x barcode in the read name; the scifi chain emits .bc1.bc2), SAM to BAM.
# Archived as used; see the README note on which wrapper produced the SM2 raw BAMs.
# Usage (from ${PROJECT_ROOT}/3_Mapping): sbatch 1_1_align.scifi.at.sh <lib>

: "${PROJECT_ROOT:?set PROJECT_ROOT to the analysis tree that holds 2_CleanReads/ and 3_Mapping/}"
: "${GENOMES_DIR:?set GENOMES_DIR to the directory holding the BWA indexes}"

########## Load Modules #########
ml Bioinformatics bwa/0.7.17-mil4ns7 samtools/1.13-fwwss5n   # HPC modules; adapt to the local installation

# INPUT
REF="${GENOMES_DIR}/Arabidopsis/Index_AraTAIR10_bwa"
INPUT="${PROJECT_ROOT}/2_CleanReads"
SAMPLE=$1;

# VARIABLES
r1="$INPUT/${SAMPLE}_R1.bc1.fastq.gz"
r2="$INPUT/${SAMPLE}_R3.bc1.fastq.gz"

bwa mem -M -t 30 $REF $r1 $r2 > ${SAMPLE}_At_scifiATAC.raw.sam;

samtools view -@ 30 -bS ${SAMPLE}_At_scifiATAC.raw.sam > ${SAMPLE}_At_scifiATAC.raw.bam
