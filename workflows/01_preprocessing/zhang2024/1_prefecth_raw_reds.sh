#!/usr/bin/env bash
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=10G
#SBATCH --job-name=srr_array
#SBATCH --partition=standard
#SBATCH --output=_logs/srr_array_%A_%a.log
#SBATCH --array=0-8
#
# Download the nine SRA runs of the Zhang et al. 2024 libraries (SRR25320539 to SRR25320547: B73Mo17_rep1,
# B73Mo17_rep2 and multiGenotypes_rep1, three sequencing runs each) with prefetch + fasterq-dump, one run
# per array task, into ${PROJECT_ROOT}/1_RawData/<SRR>/. Follow with rename_sample.sh.
# Usage (from ${PROJECT_ROOT}/1_RawData): sbatch <repo>/workflows/01_preprocessing/zhang2024/1_prefecth_raw_reds.sh

set -euo pipefail
: "${PROJECT_ROOT:?set PROJECT_ROOT to the analysis tree that holds 1_RawData/}"
REPO_ROOT="${REPO_ROOT:-$(git rev-parse --show-toplevel 2>/dev/null || pwd)}"
cd "${PROJECT_ROOT}/1_RawData"
mkdir -p _logs

# HPC module for the SRA toolkit (sratoolkit 3.1.1 was used); adapt to the local installation
ml Bioinformatics sratoolkit/3.1.1

LIST="${REPO_ROOT}/config/scifi_Metadata_sra.only.txt"

# Read the SRR ID based on array index (1-indexed for sed)
SRR_ID="$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$LIST" | tr -d '\r' | xargs)"

if [[ -z "${SRR_ID}" ]]; then
  echo "ERROR: No SRR ID found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} in ${LIST}" >&2
  exit 1
fi

echo "Processing ${SRR_ID} on task ${SLURM_ARRAY_TASK_ID} (job ${SLURM_JOB_ID})"

# one directory per SRR
OUTDIR="${SRR_ID}"
mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

# Download to the SRA cache (prefetch) then convert to FASTQ
prefetch --max-size 100G "${SRR_ID}"

fasterq-dump \
  --split-files \
  --include-technical \
  --threads "${SLURM_CPUS_PER_TASK}" \
  "${SRR_ID}"

# Compress
pigz -p "${SLURM_CPUS_PER_TASK}" "${SRR_ID}"_*.fastq

echo "Done: ${SRR_ID}"
