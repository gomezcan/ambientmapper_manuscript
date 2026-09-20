#!/usr/bin/env bash
#SBATCH --job-name=06_47_prep
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=_logs/06_47_%x_%j.log
#
# Chunked WASP — STEP 1/4: split the (raw|clean) BAM into 10 per-chromosome slices.
# =============================================================================
# Region slicing (samtools view -b $BAM $chrom) uses the .bai to seek — it does NOT
# scan the whole BAM, so 10 slices cost ~one read of the file. Slices are named
# ${BASE}.${chrom}.bam and are coordinate-sorted (region subset of a sorted BAM), so
# WASP stage 1 can stream them with --is_sorted (no per-slice .bai needed).
# Idempotent: skips any slice that already exists non-empty.
#
#   sbatch --export=ALL,SAMPLE=B73Mo17_rep1,MODE=raw workflows/03b_variant_based_comparison/06_47_wasp_chunked_prep.sh
# =============================================================================
set -euo pipefail

REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"   # this repository (sourced helpers, sibling scripts)
PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"   # root of the data tree (3_Mapping/, 5_AmbientDetection/, ...)
: "${MODE:?need MODE=raw|clean}"
# SLURM copies the batch script to a spool dir, so $BASH_SOURCE can't find siblings;
# the shared helper is located through REPO_ROOT instead.
source "${REPO_ROOT}/workflows/03b_variant_based_comparison/06_47_wasp_chunked_common.sh"

# samtools: base has it; fall back to a module. unset PYTHONPATH afterward (harmless).
command -v samtools >/dev/null 2>&1 || module load Bioinformatics samtools 2>/dev/null || true
unset PYTHONPATH PYTHONHOME
command -v samtools >/dev/null 2>&1 || { echo "samtools not found"; exit 1; }

echo "[06_47/prep] sample=$SAMPLE mode=$MODE"
echo "[06_47/prep]   BAM=$BAM"
[[ -f "$BAM" ]]     || { echo "Missing BAM: $BAM"; exit 1; }
[[ -f "$BAM.bai" ]] || { echo "Missing BAM index (.bai): $BAM.bai (region slicing needs it)"; exit 1; }
[[ -e "$WASP/mapping/find_intersecting_snps.py" ]] || { echo "WASP not set up: $WASP — run 06_46 STEP=setup"; exit 1; }
[[ -s "$SNPDIR/chr1.snps.txt.gz" ]] || { echo "snp_dir empty: $SNPDIR — run 06_46 STEP=setup,SAMPLE=$SAMPLE first"; exit 1; }
[[ -f "$BWA_IDX.bwt" ]] || { echo "bwa index missing: $BWA_IDX.bwt"; exit 1; }

mkdir -p "$SLICES"

run_slice() {  # $1 = chrom ; uses exported SLICES/BASE/BAM
  local chrom="$1"
  local out="$SLICES/${BASE}.${chrom}.bam"
  if [[ -s "$out" ]]; then echo "[06_47/prep]   $chrom: slice exists, skip"; return 0; fi
  local tmp="$out.tmp.$$"
  if samtools view -b -o "$tmp" "$BAM" "$chrom"; then
    mv -f "$tmp" "$out"
    echo "[06_47/prep]   $chrom: $(samtools view -c "$out") reads -> $out"
  else
    rm -f "$tmp"; echo "[06_47/prep]   $chrom: FAILED"; return 1
  fi
}
export -f run_slice
export SLICES BASE BAM

echo "[06_47/prep] splitting into 10 per-chrom slices ($NCPU-way) -> $SLICES"
printf '%s\n' $CHROMS_WS | xargs -P "$NCPU" -I{} bash -c 'run_slice "$@"' _ {}

echo "[06_47/prep] done. slices:"
ls -lh "$SLICES"/*.bam | awk '{print "  "$5"\t"$9}'
