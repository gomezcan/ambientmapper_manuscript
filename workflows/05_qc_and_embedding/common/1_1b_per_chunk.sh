#!/usr/bin/env bash
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=130G
#SBATCH --job-name=QCs_chunkSOC
#SBATCH --partition=standard
#SBATCH --output=_logs/%x-%A_%a.log
#
# 1_1b_per_chunk.sh — SLURM array task: run 1_1b_per_chunk.R on one BED chunk.
#
# Submitted by submit_chunked_pipeline_SM2v2.sh with --array=0-(N-1) and
# --export=ALL,POOL=...,OUTDIR=...,N=...
#
# Sentinel resume: skipped per-chunk if .soc.rds already exists.

set -euo pipefail

# Paths are anchored on PROJECT_ROOT ($BASH_SOURCE self-location does not survive sbatch).
PROJECT_ROOT="${PROJECT_ROOT:-.}"
BASE="${PROJECT_ROOT}/6_socrates"
SCRIPTS="${SCRIPTS:-${BASE}/0_scripts}"
cd "${BASE}"

POOL="${POOL:?POOL not set; submit via submit_chunked_pipeline_SM2v2.sh}"
OUTDIR="${OUTDIR:?OUTDIR not set}"
N="${N:-5}"
IDX="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID not set}"

# conda activate ambientmapper-manuscript   (environment.yml at the repo root)
# Socrates' loadBEDandGenomeData() runs .preRunChecks() which requires macs2
# on PATH even when we never call it (peaks come from 0_07's MACS3 output).
# module load macs2   (site-specific; Socrates' .preRunChecks() needs macs2 on PATH even with precomputed peaks)

CHUNKDIR="${BASE}/${OUTDIR}/chunks/${POOL}"
PEAKROOT="${BASE}/_data/_PeakFiles/${POOL}"

ann="${BASE}/_data/_GenomeInfo/ZmATcombined.gtf"
chr="${BASE}/_data/_GenomeInfo/ZmATcombined.chrs.size.txt"
peakpath="${PEAKROOT}/${POOL}_peaks.narrowPeak"

chunk_bed="${CHUNKDIR}/${POOL}.chunk${IDX}.bed.gz"
out_prefix="${POOL}.chunk${IDX}"
out_soc="${CHUNKDIR}/${out_prefix}.soc.rds"

[[ -e "$chunk_bed" ]] || { echo "ERROR: chunk BED not found: $chunk_bed (run 1_1a first)"; exit 1; }
[[ -e "$peakpath"  ]] || { echo "ERROR: MACS3 peaks not found: $peakpath (run 0_07 first)"; exit 1; }
[[ -e "$ann"       ]] || { echo "ERROR: GTF not found: $ann"; exit 1; }
[[ -e "$chr"       ]] || { echo "ERROR: chr_sizes not found: $chr"; exit 1; }

if [ -s "$out_soc" ]; then
  echo " - SKIP chunk $IDX: $out_soc already present ($(du -h "$out_soc" | cut -f1))"
  exit 0
fi

# R script outputs are prefix-relative; cd into the chunk dir.
cd "$CHUNKDIR"

echo " - running 1_1b for $POOL chunk $IDX of $N"
echo "   BED   : $chunk_bed"
echo "   peaks : $peakpath"
echo "   ann   : $ann"
echo "   chr   : $chr"
Rscript "$SCRIPTS/common/1_1b_per_chunk.R" "$chunk_bed" "$out_prefix" "$ann" "$chr" "$peakpath"

echo " - chunk $IDX done: $out_soc"
ls -lh "$out_soc"
