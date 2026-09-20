#!/usr/bin/env bash
#SBATCH --job-name=b73mo17r1_tn5regen
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --array=0-1
#SBATCH --output=_logs/04_07c_tn5bed_regen_B73Mo17_rep1_%A_%a.log
#
# 04_07c_tn5bed_regen_B73Mo17_rep1.sh — forced regeneration of the tn5-bed step
# for B73Mo17_rep1 (B73v5, Mo17) after a `sort -u` collapse bug fix in
# workflows/02_mapping/tn5bed/00_bam_to_tn5bed_parallel.py (the earlier sidecar
# deduplicated on (chrom,start) only, dropping ~90% of records).
#
# The sidecar's --skip-existing flag no-ops when --out exists and is non-empty, so
# the collapsed beds would be skipped by a plain re-run. This script omits
# --skip-existing so the sidecar recomputes and ATOMICALLY OVERWRITES the bed
# (it writes <out>.tmp then os.replace; failure leaves the old bed untouched).
#
# INPUTS unaffected by the bug: the .Clean.bam + .bai from 04_07 are correct (the
# bug only hit the tn5 derivative). This step re-derives tn5 from those BAMs.
#
# Array task layout:
#   0: B73Mo17_rep1 / B73v5
#   1: B73Mo17_rep1 / Mo17
# =============================================================================

set -euo pipefail

REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"   # this repository (sourced helpers, sibling scripts)
PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"   # root of the data tree (3_Mapping/, 5_AmbientDetection/, ...)

cd "${PROJECT_ROOT}/5_AmbientDetection"

# conda activate <env from environment.yml>

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export GOTO_NUM_THREADS=1

SAMPLE="B73Mo17_rep1"
genomes=(B73v5 Mo17)

idx="${SLURM_ARRAY_TASK_ID:?This script must be submitted as an array job (sbatch --array=0-1 ...)}"
genome="${genomes[$idx]}"

bam_in="${SAMPLE}/clean_bams_alpha05_C0_nd/${SAMPLE}_${genome}_scifiATAC.mq10.BC.rmdup.mm.Clean.bam"
out_bed_gz="${bam_in%.bam}.tn5.bed.gz"

[[ -f "$bam_in"       ]] || { echo "Missing Clean BAM: $bam_in"; exit 1; }
[[ -f "${bam_in}.bai" ]] || { echo "Missing index:     ${bam_in}.bai"; exit 1; }

threads="${SLURM_CPUS_PER_TASK:-12}"
# Reserve a couple of cpus for pigz; rest fan out across chroms.
pigz_threads=2
worker_threads=$(( threads - pigz_threads ))
(( worker_threads >= 1 )) || worker_threads=1

scratch="${TMPDIR:-/tmp}"

echo "============================================================"
echo "[$(date)] B73Mo17_rep1 tn5-bed REGEN task ${idx}: ${SAMPLE} / ${genome}"
echo "  bam_in         : ${bam_in} ($(du -h "$bam_in" | cut -f1))"
echo "  out_bed_gz     : ${out_bed_gz} (will be OVERWRITTEN)"
if [[ -s "$out_bed_gz" ]]; then
  echo "  prior bed      : $(du -h "$out_bed_gz" | cut -f1) (prior version, replacing)"
fi
echo "  threads        : ${threads} (workers=${worker_threads}, pigz=${pigz_threads})"
echo "  scratch        : ${scratch}"
echo "============================================================"

# NOTE: --skip-existing intentionally OMITTED so the corrupt bed is overwritten.
python "${REPO_ROOT}/workflows/02_mapping/tn5bed/00_bam_to_tn5bed_parallel.py" \
  --bam "$bam_in" \
  --out "$out_bed_gz" \
  --threads "$worker_threads" \
  --pigz-threads "$pigz_threads" \
  --scratch-dir "$scratch"

[[ -s "$out_bed_gz" ]] || { echo "ERROR: tn5.bed.gz missing or empty: $out_bed_gz"; exit 1; }

# ---- correctness gate: confirm the collapse signature is GONE ----------------
# A correct tn5 bed has many cells sharing insertion sites, so
# total records >> distinct (chrom,start). The collapse bug produced beds where
# total == distinct. Fail loudly if that signature is still present.
echo "[verify] collapse check (first 1,000,000 records)..."
stats=$(gzip -dc "$out_bed_gz" 2>/dev/null | head -1000000 \
        | awk '{k=$1"\t"$2; if(k!=p){d++; p=k} t++} END{printf "%d %d", t, d}') || true
tot=${stats% *}; pos=${stats#* }
echo "[verify] sampled records=${tot:-0} distinct(chrom,start)=${pos:-0}"
if [[ -z "${tot:-}" || -z "${pos:-}" || "$tot" -eq "$pos" ]]; then
  echo "ERROR: collapse signature present (records == distinct) or check failed. Aborting." >&2
  exit 1
fi
echo "[verify] OK — not collapsed (ratio ~$(awk "BEGIN{printf \"%.2f\", $tot/$pos}")x records:positions)."

echo "============================================================"
echo "[$(date)] task ${idx} DONE for ${SAMPLE} / ${genome}"
ls -la "$(dirname "$out_bed_gz")"/"${SAMPLE}_${genome}_scifiATAC.mq10.BC.rmdup.mm.Clean".*
echo "============================================================"
