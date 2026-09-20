#!/usr/bin/env bash
#SBATCH --job-name=06_47_combine
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=_logs/06_47_%x_%j.log
#
# Chunked WASP — STEP 3/4: merge the 10 per-chrom (keep + remap.keep) bams into one
# WASP-passed BAM, reconcile read accounting, then re-count alleles with the SAME
# 06_44 extractor (identical flags to monolithic 06_46 -> directly comparable numbers).
# =============================================================================
# Merge assumes ~coordinate-sorted inputs and the explicit sort afterward guarantees a
# correctly-sorted final BAM regardless; the read-accounting check (wasp.pass == Σkeep +
# Σremap.kept) proves the merge is lossless. Same merge->sort->extract as 06_46.
# Optional: CLEAN_PARTS=1 deletes the per-chrom work dirs after a successful run.
#
#   sbatch --export=ALL,SAMPLE=B73Mo17_rep1,MODE=raw workflows/03b_variant_based_comparison/06_47_wasp_chunked_combine.sh
# =============================================================================
set -euo pipefail

REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"   # this repository (sourced helpers, sibling scripts)
PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"   # root of the data tree (3_Mapping/, 5_AmbientDetection/, ...)
: "${MODE:?need MODE=raw|clean}"
# SLURM copies the batch script to a spool dir, so $BASH_SOURCE can't find siblings;
# the shared helper is located through REPO_ROOT instead.
source "${REPO_ROOT}/workflows/03b_variant_based_comparison/06_47_wasp_chunked_common.sh"

# samtools: base has it; fall back to a module. $PY stays absolute-base after unset.
command -v samtools >/dev/null 2>&1 || module load Bioinformatics samtools 2>/dev/null || true
unset PYTHONPATH PYTHONHOME
SAMTOOLS="${SAMTOOLS:-samtools}"
command -v "$SAMTOOLS" >/dev/null || { echo "samtools not found"; exit 1; }

echo "[06_47/combine] sample=$SAMPLE mode=$MODE"

# guard: all 10 chroms have both canonical bams
missing=0; parts=()
for c in $CHROMS_WS; do
  k="$(chrom_keep "$c")"; r="$(chrom_rkeep "$c")"
  if [[ -s "$k" && -s "$r" ]]; then
    parts+=("$k" "$r")
  else
    echo "[06_47/combine]   MISSING $c: keep=$([[ -s "$k" ]] && echo ok || echo NO) remap.keep=$([[ -s "$r" ]] && echo ok || echo NO)"
    missing=1
  fi
done
[[ "$missing" -eq 0 ]] || { echo "[06_47/combine] ERROR: not all 10 chroms complete — re-fire the missing array task(s)"; exit 1; }

WASP_BAM=$OUT/${SAMPLE}_${MODE}.wasp.bam
UNSORTED=$OUT/${SAMPLE}_${MODE}.wasp.unsorted.bam
echo "[06_47/combine] merging ${#parts[@]} per-chrom bams -> $WASP_BAM"
"$SAMTOOLS" merge -f -@ "$NCPU" "$UNSORTED" "${parts[@]}"
"$SAMTOOLS" sort -@ "$NCPU" -o "$WASP_BAM" "$UNSORTED"
"$SAMTOOLS" index -@ "$NCPU" "$WASP_BAM"
rm -f "$UNSORTED"

# BC-tag sanity — 06_44 groups reads by the BC tag
nbc=$("$SAMTOOLS" view "$WASP_BAM" 2>/dev/null | head -100000 | grep -c 'BC:Z:' || true)
echo "[06_47/combine]   BC-tagged reads in first 100k of wasp.bam: $nbc"
[[ "$nbc" -gt 0 ]] || { echo "[06_47/combine] ERROR: BC tag missing — 06_44 grouping would fail"; exit 1; }

# read-accounting reconciliation: wasp.pass == Σ keep + Σ remap.kept
sum_keep=0; sum_rkeep=0
for c in $CHROMS_WS; do
  sum_keep=$((  sum_keep  + $("$SAMTOOLS" view -c "$(chrom_keep  "$c")") ))
  sum_rkeep=$(( sum_rkeep + $("$SAMTOOLS" view -c "$(chrom_rkeep "$c")") ))
done
wasp_pass=$("$SAMTOOLS" view -c "$WASP_BAM")
echo "[06_47/combine]   Σkeep=$sum_keep  Σremap.kept=$sum_rkeep  sum=$(( sum_keep + sum_rkeep ))  wasp.pass=$wasp_pass"
[[ "$wasp_pass" -eq "$(( sum_keep + sum_rkeep ))" ]] \
  || { echo "[06_47/combine] ERROR: read-accounting mismatch — merge lost/gained reads"; exit 1; }

# re-count alleles with the SAME 06_44 extractor (identical flags to 06_46 step 5)
set_group_args
echo "[06_47/combine] 06_44 extractor on wasp.bam (group_mode=$GROUP_MODE)"
$PY "${REPO_ROOT}/workflows/03b_variant_based_comparison/06_44_extract_allele_counts.py" \
  --bam "$WASP_BAM" --vcf "$EXTRACT_VCF" \
  "${GROUP_ARGS[@]}" \
  --bc-tag BC --min-mapq 10 --threads "$NCPU" --chroms "$CHROMS_CSV" \
  --out "$OUT/${SAMPLE}_${MODE}_wasp_allele_counts.tsv.gz"

if [[ "${CLEAN_PARTS:-0}" == 1 ]]; then
  echo "[06_47/combine] CLEAN_PARTS=1 -> removing per-chrom work dirs"
  for c in $CHROMS_WS; do rm -rf "$(chrom_dir "$c")"; done
fi
echo "[06_47/combine] done mode=$MODE -> $WASP_BAM + $OUT/${SAMPLE}_${MODE}_wasp_allele_counts.tsv.gz"
