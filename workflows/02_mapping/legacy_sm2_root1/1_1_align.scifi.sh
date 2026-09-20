#!/bin/bash
########## BATCH Lines for Resource Request ##########
#SBATCH --time=24:00:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=80G
#SBATCH --job-name=BWA_multiome
#SBATCH --partition=standard
#SBATCH --output=_logs/%x-%j.log
#
# Legacy single-reference alignment wrapper (maize B73v5 index): bwa mem -M of the demultiplexed
# <lib>_R{1,3}_001.bc1.bc2.fastq.gz pair, SAM to BAM. Archived as used; see the README note on which
# wrapper produced the SM2 raw BAMs.
# Usage (from ${PROJECT_ROOT}/3_Mapping): sbatch 1_1_align.scifi.sh <lib>

: "${PROJECT_ROOT:?set PROJECT_ROOT to the analysis tree that holds 2_CleanReads/ and 3_Mapping/}"
: "${GENOMES_DIR:?set GENOMES_DIR to the directory holding the BWA indexes}"

########## Load Modules #########
ml Bioinformatics bwa/0.7.17-mil4ns7 samtools/1.13-fwwss5n   # HPC modules; adapt to the local installation

# INPUT
REF="${GENOMES_DIR}/Zea/Index_B73v5_btw"   # BWA index of Zm-B73-REFERENCE-NAM-5.0 (chromosomes, Mt, Pt)
INPUT="${PROJECT_ROOT}/2_CleanReads"
SAMPLE=$1;

# VARIABLES
r1="$INPUT/${SAMPLE}_R1_001.bc1.bc2.fastq.gz"
r2="$INPUT/${SAMPLE}_R3_001.bc1.bc2.fastq.gz"

bwa mem -M -t 72 $REF $r1 $r2 > ${SAMPLE}_scifiATAC.raw.sam;

samtools view -@ 72 -bS ${SAMPLE}_scifiATAC.raw.sam > ${SAMPLE}_scifiATAC.raw.bam
