#!/usr/bin/env bash
#SBATCH --job-name=06_47_wasp
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --array=1-10
#SBATCH --output=_logs/06_47_%x_%A_%a.log
#
# Chunked WASP — STEP 2/4: run WASP on ONE chromosome slice (array task -> chrom).
# =============================================================================
# Per array task (chrom = CHROMS_ARR[$SLURM_ARRAY_TASK_ID-1]):
#   1. find_intersecting_snps.py  (--is_paired_end --is_sorted)  -> keep + to.remap + remap.fq{1,2}
#   2. bwa mem -M  (ORIGINAL aligner/params, full genome-wide index) | samtools sort
#   3. filter_remapped_reads.py   -> remap.keep.bam
# then delete this chrom's big temps, keeping only canonical keep.bam + remap.keep.bam.
#
# Behavior-preserving vs monolithic 06_46 (see plan doc): WASP already flushes its
# read-pair cache at each chromosome boundary, remap uses the full bwa index, and
# inter-chromosomal pairs are dropped identically. Two WASP-source watch-outs are
# honored: (A) do NOT coordinate-sort to.remap.bam (filter relies on FIS write order);
# (D) single-end remap reads are not remapped (== 06_46) — their count is surfaced.
#
# Idempotency guard: skips a chrom whose keep.bam + remap.keep.bam already exist, so a
# single failed chrom re-fires cheaply:  sbatch --array=3 --export=ALL,... this_script
# 96 G is conservative (monolithic OOM'd >72 G on the WHOLE genome; one chrom ≈ 16%).
#
#   sbatch --export=ALL,SAMPLE=B73Mo17_rep1,MODE=raw workflows/03b_variant_based_comparison/06_47_wasp_chunked_array.sh
# =============================================================================
set -euo pipefail

REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"   # this repository (sourced helpers, sibling scripts)
PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"   # root of the data tree (3_Mapping/, 5_AmbientDetection/, ...)
: "${MODE:?need MODE=raw|clean}"
: "${SLURM_ARRAY_TASK_ID:?run as a SLURM array (--array=1-10)}"
# SLURM copies the batch script to a spool dir, so $BASH_SOURCE can't find siblings;
# the shared helper is located through REPO_ROOT instead.
source "${REPO_ROOT}/workflows/03b_variant_based_comparison/06_47_wasp_chunked_common.sh"

chrom="${CHROMS_ARR[$((SLURM_ARRAY_TASK_ID-1))]:-}"
[[ -n "$chrom" ]] || { echo "bad array id $SLURM_ARRAY_TASK_ID (need 1-10)"; exit 1; }

# bwa via module (0.7.17 that mapped these BAMs); scrub the py3.9/PYTHONPATH injection
# so $PY runs on clean base py3.12 (see the env note in 06_47_wasp_chunked_common.sh).
module load Bioinformatics "$BWA_MODULE" 2>/dev/null || true
command -v samtools >/dev/null 2>&1 || module load samtools 2>/dev/null || true
unset PYTHONPATH PYTHONHOME
BWA="${BWA:-bwa}"; SAMTOOLS="${SAMTOOLS:-samtools}"

slice="$(slice_path "$chrom")"
cdir="$(chrom_dir "$chrom")"
keep="$(chrom_keep "$chrom")"      # $cdir/keep.bam        (canonical; WASP output renamed)
rkeep="$(chrom_rkeep "$chrom")"    # $cdir/remap.keep.bam  (filter_remapped_reads output)

echo "[06_47/wasp] sample=$SAMPLE mode=$MODE chrom=$chrom (array $SLURM_ARRAY_TASK_ID)"

# idempotency: both canonical outputs present & non-empty -> already done
if [[ -s "$keep" && -s "$rkeep" ]]; then
  echo "[06_47/wasp]   $chrom already done (keep + remap.keep present) — skip"; exit 0
fi
[[ -s "$slice" ]]        || { echo "Missing slice: $slice — run prep first"; exit 1; }
[[ -f "$BWA_IDX.bwt" ]]  || { echo "bwa index missing: $BWA_IDX.bwt"; exit 1; }
command -v "$BWA" >/dev/null || { echo "bwa not found ($BWA) — module load Bioinformatics $BWA_MODULE failed; set BWA=/abs/path/bwa"; exit 1; }
echo "[06_47/wasp]   bwa=$($BWA 2>&1 | awk '/^Version/{print $2; exit}')  samtools=$($SAMTOOLS --version 2>/dev/null | awk 'NR==1{print $2}')  py=$($PY --version 2>&1)"

# fresh per-chrom workspace (we are (re)doing this chrom; wipe any partial state)
rm -rf "$cdir"; mkdir -p "$cdir"

# 1. find reads overlapping SNPs -> allele-swapped remap FASTQs (this chrom only)
echo "[06_47/wasp] (1) find_intersecting_snps ($chrom)"
$PY "$WASP/mapping/find_intersecting_snps.py" \
  --is_paired_end --is_sorted \
  --output_dir "$cdir" --snp_dir "$SNPDIR" \
  "$slice"
# resolve WASP outputs by glob (independent of WASP's basename scheme)
kout=$(ls "$cdir"/*.keep.bam            2>/dev/null | head -1 || true)
toremap=$(ls "$cdir"/*.to.remap.bam     2>/dev/null | head -1 || true)
fq1=$(ls "$cdir"/*.remap.fq1.gz         2>/dev/null | head -1 || true)
fq2=$(ls "$cdir"/*.remap.fq2.gz         2>/dev/null | head -1 || true)
single=$(ls "$cdir"/*.remap.single.fq.gz 2>/dev/null | head -1 || true)
[[ -s "$kout" && -s "$toremap" && -s "$fq1" && -s "$fq2" ]] \
  || { echo "[06_47/wasp] ERROR: expected WASP outputs missing in $cdir"; ls -la "$cdir"; exit 1; }
mv -f "$kout" "$keep"

# caveat D: single-end remap reads are NOT remapped (matches monolithic 06_46, which
# only remapped fq1/fq2). Surface the count loudly instead of dropping it silently.
if [[ -n "$single" && -s "$single" ]]; then
  nsingle=$(( $(zcat "$single" | wc -l) / 4 ))
  [[ "$nsingle" -eq 0 ]] || echo "[06_47/wasp]   WARNING: $chrom has $nsingle single-end remap reads — dropped (as in 06_46). Nonzero for ATAC warrants a look."
fi

# 2. remap allele-swapped reads with the ORIGINAL aligner/params (bwa mem -M).
#    Uses the FULL genome-wide index — slicing groups input reads, never restricts remap.
echo "[06_47/wasp] (2) bwa mem remap ($chrom)"
"$BWA" mem -M -t "$NCPU" "$BWA_IDX" "$fq1" "$fq2" \
  | "$SAMTOOLS" sort -@ "$NCPU" -o "$cdir/remap.sorted.bam" -
"$SAMTOOLS" index -@ "$NCPU" "$cdir/remap.sorted.bam"

# 3. keep only reads that remapped back to the same locus.
#    NB (caveat A): do NOT coordinate-sort $toremap — the pair-cache relies on FIS order.
echo "[06_47/wasp] (3) filter_remapped_reads ($chrom)"
$PY "$WASP/mapping/filter_remapped_reads.py" "$toremap" "$cdir/remap.sorted.bam" "$rkeep"

# per-chrom read accounting
echo "[06_47/wasp]   $chrom: slice=$($SAMTOOLS view -c "$slice") keep=$($SAMTOOLS view -c "$keep") to.remap=$($SAMTOOLS view -c "$toremap") remap.kept=$($SAMTOOLS view -c "$rkeep")"

# drop this chrom's big temps; keep only keep.bam + remap.keep.bam for combine
rm -f "$toremap" "$cdir/remap.sorted.bam" "$cdir/remap.sorted.bam.bai" "$fq1" "$fq2" "$single" "$slice"
echo "[06_47/wasp] $chrom done -> $keep + $rkeep"
