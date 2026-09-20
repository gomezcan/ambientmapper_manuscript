#!/usr/bin/env bash
# Rename the fasterq-dump outputs <SRR>/<SRR>_*.fastq.gz to the library ids used downstream
# (e.g. SRR25320547/SRR25320547_R3.fastq.gz -> SRR25320547/scifi_B73Mo17_rep1_1_R3.fastq.gz)
# following config/scifi_Metadata_sra.clean.txt (columns SampleID, Run, LibraryLayout).
# Usage (from ${PROJECT_ROOT}/1_RawData): bash <repo>/workflows/01_preprocessing/zhang2024/rename_sample.sh

set -euo pipefail
: "${PROJECT_ROOT:?set PROJECT_ROOT to the analysis tree that holds 1_RawData/}"
REPO_ROOT="${REPO_ROOT:-$(git rev-parse --show-toplevel 2>/dev/null || pwd)}"
cd "${PROJECT_ROOT}/1_RawData"
META="${REPO_ROOT}/config/scifi_Metadata_sra.clean.txt"

tail -n +2 "${META}" | while read -r sample srr layout; do
  for f in "${srr}"/"${srr}"_*.gz; do
    [ -e "$f" ] || continue  # skip if the file does not exist
    name="$(basename "$f")"
    newname="${name/${srr}/${sample}}"
    mv "$f" "${srr}/${newname}"
  done
done
