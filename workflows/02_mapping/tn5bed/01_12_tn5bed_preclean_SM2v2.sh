#!/usr/bin/env bash
#SBATCH --job-name=sm2v2_tn5preclean
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --array=0-1
#SBATCH --output=_logs/01_12_tn5bed_preclean_SM2v2_%A_%a.log
#
# 01_12_tn5bed_preclean_SM2v2.sh: PreClean per-genome Tn5 BEDs for SM2v2, the counterpart of the PostClean
# per-genome BEDs written by 01_11c_tn5bed_regen_SM2v2.sh.
#
# They are derived from the merged AmbientMapper input BAMs (01_00_merge_SM2v2_inputs.sh) rather than from
# the older per-library BEDs, so PreClean and PostClean share one BAM lineage: the PreClean BED is the same
# BAM before the reads_to_drop list is applied (hence no ND/WD split), merged across both plate halves and
# carrying both library barcode suffixes.
#
# Input BAMs:
#   3_Mapping/ambientmapper_input/SM2_B73v5_scifiATAC.mq10.BC.rmdup.mm.bam
#   3_Mapping/ambientmapper_input/SM2_TAIR10_scifiATAC.mq10.BC.rmdup.mm.bam
# Outputs (consumed by the independent, per-genome Socrates track as the PreClean BEDs):
#   5_AmbientDetection/SM2v2/preclean_tn5beds/SM2_B73v5_scifiATAC.mq10.BC.rmdup.mm.tn5.bed.gz
#   5_AmbientDetection/SM2v2/preclean_tn5beds/SM2_TAIR10_scifiATAC.mq10.BC.rmdup.mm.tn5.bed.gz
#
# Array task layout:  0: B73v5   1: TAIR10
# Usage (from ${PROJECT_ROOT}/5_AmbientDetection): sbatch 01_12_tn5bed_preclean_SM2v2.sh
# =============================================================================

set -euo pipefail
: "${PROJECT_ROOT:?set PROJECT_ROOT to the analysis tree that holds 3_Mapping/ and 5_AmbientDetection/}"
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

INPUT_DIR="${PROJECT_ROOT}/3_Mapping/ambientmapper_input"
OUT_DIR="SM2v2/preclean_tn5beds"
SIDECAR="${REPO_ROOT}/workflows/02_mapping/tn5bed/00_bam_to_tn5bed_parallel.py"
genomes=(B73v5 TAIR10)

idx="${SLURM_ARRAY_TASK_ID:?This script must be submitted as an array job (sbatch --array=0-1 ...)}"
genome="${genomes[$idx]}"

bam_in="${INPUT_DIR}/SM2_${genome}_scifiATAC.mq10.BC.rmdup.mm.bam"
bam_base="$(basename "$bam_in" .bam)"
out_bed_gz="${OUT_DIR}/${bam_base}.tn5.bed.gz"

[[ -f "$bam_in"       ]] || { echo "Missing input BAM: $bam_in"; exit 1; }
[[ -f "${bam_in}.bai" ]] || { echo "Missing input BAI: ${bam_in}.bai"; exit 1; }
[[ -f "$SIDECAR"      ]] || { echo "Missing extractor: $SIDECAR"; exit 1; }

mkdir -p "$OUT_DIR"

threads="${SLURM_CPUS_PER_TASK:-12}"
pigz_threads=2
worker_threads=$(( threads - pigz_threads ))
(( worker_threads >= 1 )) || worker_threads=1

scratch="${TMPDIR:-/tmp}"

echo "============================================================"
echo "[$(date)] SM2v2 tn5-bed PRECLEAN task ${idx}: ${genome}"
echo "  bam_in     : ${bam_in} ($(du -h "$bam_in" | cut -f1))"
echo "  out_bed_gz : ${out_bed_gz}"
echo "  threads    : ${threads} (workers=${worker_threads}, pigz=${pigz_threads})"
echo "  scratch    : ${scratch}"
echo "============================================================"

python "${SIDECAR}" \
  --bam "$bam_in" \
  --out "$out_bed_gz" \
  --threads "$worker_threads" \
  --pigz-threads "$pigz_threads" \
  --scratch-dir "$scratch" \
  --skip-existing

[[ -s "$out_bed_gz" ]] || { echo "ERROR: tn5.bed.gz missing or empty: $out_bed_gz"; exit 1; }

# ---- correctness gate: confirm the output is NOT collapsed --------------------
# A correct tn5 bed has many cells sharing insertion sites, so total records >> distinct (chrom,start)
# (the dedup key of the extractor is the full record: chrom, cut, barcode, strand).
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
echo "[$(date)] task ${idx} DONE for ${genome}"
ls -la "$OUT_DIR"/"${bam_base}".*
echo "============================================================"
