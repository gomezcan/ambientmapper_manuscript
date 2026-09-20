#!/usr/bin/env bash
#SBATCH --job-name=sm2v2_tn5regen
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --array=0-3
#SBATCH --output=_logs/01_11c_tn5bed_regen_SM2v2_%A_%a.log
#
# 01_11c_tn5bed_regen_SM2v2.sh: PostClean per-genome Tn5 BEDs for SM2v2 from the cleaned BAMs of both
# decontamination passes (ND = without plate design, WD = with plate design; B73v5 and TAIR10).
#
# --skip-existing is deliberately omitted, so an existing BED is recomputed and replaced atomically
# (the extractor writes <out>.tmp and renames it; a failed run leaves the old BED untouched).
# These BEDs are consumed by the independent, per-genome Socrates track; the combined-genome track uses
# the ZmATcombined BEDs of 0_06_make_tn5bed_SM2v2_combined.sh.
#
# Input:  5_AmbientDetection/SM2v2/clean_bams_alpha05_C0_{nd,wd}/SM2_{B73v5,TAIR10}_scifiATAC.mq10.BC.rmdup.mm.Clean.bam (+ .bai)
# Output: same directory, <bam without .bam>.tn5.bed.gz
# Filename convention: the working directory is SM2v2/ but the cleaned BAMs keep the SM2_ file prefix.
#
# Array task layout (pass x genome):
#   0: nd / B73v5
#   1: nd / TAIR10
#   2: wd / B73v5
#   3: wd / TAIR10
# Usage (from ${PROJECT_ROOT}/5_AmbientDetection): sbatch 01_11c_tn5bed_regen_SM2v2.sh
# =============================================================================

set -euo pipefail
: "${PROJECT_ROOT:?set PROJECT_ROOT to the analysis tree that holds 5_AmbientDetection/}"
REPO_ROOT="${REPO_ROOT:-$(git rev-parse --show-toplevel 2>/dev/null || pwd)}"

cd "${PROJECT_ROOT}/5_AmbientDetection"
mkdir -p _logs

# conda activate <env from environment.yml>   (python with pysam, pigz)

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export GOTO_NUM_THREADS=1

WORKDIR="SM2v2"                 # working directory
FILE_PREFIX="SM2"              # Clean BAM filename prefix (NOT the dir name)
SIDECAR="${REPO_ROOT}/workflows/02_mapping/tn5bed/00_bam_to_tn5bed_parallel.py"
passes=(nd     nd     wd     wd)
genomes=(B73v5 TAIR10 B73v5  TAIR10)

idx="${SLURM_ARRAY_TASK_ID:?This script must be submitted as an array job (sbatch --array=0-3 ...)}"
pass="${passes[$idx]}"
genome="${genomes[$idx]}"

bam_in="${WORKDIR}/clean_bams_alpha05_C0_${pass}/${FILE_PREFIX}_${genome}_scifiATAC.mq10.BC.rmdup.mm.Clean.bam"
out_bed_gz="${bam_in%.bam}.tn5.bed.gz"

[[ -f "$bam_in"       ]] || { echo "Missing Clean BAM: $bam_in"; exit 1; }
[[ -f "${bam_in}.bai" ]] || { echo "Missing index:     ${bam_in}.bai"; exit 1; }
[[ -f "$SIDECAR"      ]] || { echo "Missing extractor: $SIDECAR"; exit 1; }

threads="${SLURM_CPUS_PER_TASK:-12}"
# Reserve a couple of cpus for pigz; the rest fan out across chromosomes.
pigz_threads=2
worker_threads=$(( threads - pigz_threads ))
(( worker_threads >= 1 )) || worker_threads=1

scratch="${TMPDIR:-/tmp}"

echo "============================================================"
echo "[$(date)] SM2v2 tn5-bed REGEN task ${idx}: ${pass} / ${genome}"
echo "  bam_in         : ${bam_in} ($(du -h "$bam_in" | cut -f1))"
echo "  out_bed_gz     : ${out_bed_gz} (will be OVERWRITTEN)"
if [[ -s "$out_bed_gz" ]]; then
  echo "  prior bed      : $(du -h "$out_bed_gz" | cut -f1) (existing version, replacing)"
fi
echo "  threads        : ${threads} (workers=${worker_threads}, pigz=${pigz_threads})"
echo "  scratch        : ${scratch}"
echo "============================================================"

# --skip-existing intentionally OMITTED so an existing bed is overwritten.
python "${SIDECAR}" \
  --bam "$bam_in" \
  --out "$out_bed_gz" \
  --threads "$worker_threads" \
  --pigz-threads "$pigz_threads" \
  --scratch-dir "$scratch"

[[ -s "$out_bed_gz" ]] || { echo "ERROR: tn5.bed.gz missing or empty: $out_bed_gz"; exit 1; }

# ---- correctness gate: confirm the output is NOT collapsed --------------------
# A correct tn5 bed has many cells sharing insertion sites, so total records >> distinct (chrom,start).
# Fail loudly if total == distinct.
echo "[verify] collapse check (first 1,000,000 records)..."
stats=$(gzip -dc "$out_bed_gz" 2>/dev/null | head -1000000 \
        | awk '{k=$1"\t"$2; if(k!=p){d++; p=k} t++} END{printf "%d %d", t, d}') || true
tot=${stats% *}; pos=${stats#* }
echo "[verify] sampled records=${tot:-0} distinct(chrom,start)=${pos:-0}"
if [[ -z "${tot:-}" || -z "${pos:-}" || "$tot" -eq "$pos" ]]; then
  echo "ERROR: collapse signature present (records == distinct) or check failed. Aborting." >&2
  exit 1
fi
echo "[verify] OK, not collapsed (ratio ~$(awk "BEGIN{printf \"%.2f\", $tot/$pos}")x records:positions)."

echo "============================================================"
echo "[$(date)] task ${idx} DONE for ${pass} / ${genome}"
ls -la "$(dirname "$out_bed_gz")"/"${FILE_PREFIX}_${genome}_scifiATAC.mq10.BC.rmdup.mm.Clean".*
echo "============================================================"
