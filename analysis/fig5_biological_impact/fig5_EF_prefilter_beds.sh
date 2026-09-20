#!/usr/bin/env bash
# =============================================================================
# fig5_EF_prefilter_beds.sh  -  heavy prep for fig5_EF_qc_fripfair.R (Fig 5 panel F)
# Filters the combined-genome tn5 BEDs down to the paired fixed-set cells, per species/stage,
# so the FRiP recompute is fast. Writes <scratch>/{at,b73}_cells.txt and
# <scratch>/{at,b73}_{pre,post}_reads.tsv (chrom start end bc). Barcodes are stripped to the
# leading cell barcode (sub -.*) so pre/post match. Scratch is transient.
# Inputs  $DATA/socrates/compare/SM2_{At,B73}.fixedSet.minDepth200.pre_post.txt
#         $DATA/socrates/_data/_BED_files/{SM2,Clean.SM2v2}_{At,B73}_ZmATcombined_scifiATAC.mq10.BC.rmdup.mm.tn5.bed.gz
# Usage   bash analysis/fig5_biological_impact/fig5_EF_prefilter_beds.sh <scratch_dir>  (~5 min)
# Then    Rscript analysis/fig5_biological_impact/fig5_EF_qc_fripfair.R <scratch_dir>
# =============================================================================
set -euo pipefail
SC="${1:?usage: fig5_EF_prefilter_beds.sh <scratch_dir>}"
REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
cd "$REPO_ROOT"
DATA=data/processed/scifiATAC_B73_Arabidopsis
SOC="$DATA/socrates"
BEDDIR="$SOC/_data/_BED_files"
mkdir -p "$SC"

# species -> (compare prefix, pre BED base, post BED base)
for row in "at:SM2_At:SM2_At_ZmATcombined:Clean.SM2v2_At_ZmATcombined" \
           "b73:SM2_B73:SM2_B73_ZmATcombined:Clean.SM2v2_B73_ZmATcombined"; do
  IFS=: read -r sp cmpp preb postb <<< "$row"
  cmp="$SOC/compare/${cmpp}.fixedSet.minDepth200.pre_post.txt"
  [ -f "$cmp" ] || { echo "MISSING: $cmp" >&2; exit 1; }
  # cells present in BOTH stages -> stripped leading barcode
  awk -F'\t' 'NR==1{next} $3==1{c[$1]++} END{for(k in c) if(c[k]==2){b=k; sub(/-.*/,"",b); print b}}' \
      "$cmp" > "$SC/${sp}_cells.txt"
  echo "$sp paired cells: $(wc -l < "$SC/${sp}_cells.txt")"
  for stg in pre post; do
    if [ "$stg" = pre ]; then base="$preb"; else base="$postb"; fi
    bed="$BEDDIR/${base}_scifiATAC.mq10.BC.rmdup.mm.tn5.bed.gz"
    [ -f "$bed" ] || { echo "MISSING: $bed" >&2; exit 1; }
    gunzip -c "$bed" | awk 'NR==FNR{k[$1]=1;next}{b=$4; sub(/-.*/,"",b); if(b in k) print $1"\t"$2"\t"$3"\t"b}' \
        "$SC/${sp}_cells.txt" - > "$SC/${sp}_${stg}_reads.tsv"
    echo "  $sp $stg reads: $(wc -l < "$SC/${sp}_${stg}_reads.tsv")"
  done
done
echo "[done] prefiltered reads in $SC"
