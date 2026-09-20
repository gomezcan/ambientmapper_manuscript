#!/usr/bin/env bash
# =============================================================================
# figS7_prep_beds.sh  -  heavy prep for Fig S7 panel D (per-genome QC across cleaning modes).
# Filters the six PLATE-SPLIT tn5 BEDs (2 genomes x 3 stages) down to the cells SHARED by all
# three stages, so the FRiP recompute + peak rarefaction in figS7_fripnorm.R is fast and exactly
# paired. Writes <scratch>/{TAIR10,B73v5}_cells.txt and <scratch>/{TAIR10,B73v5}_{pre,wd,nd}_reads.tsv.
# Inputs  $DATA/socrates/SM2v2_plate/step0_qc/<stage>_<genome>.minDepth200.updated_metadata_v1.txt
#         $DATA/socrates/_data/_BED_files/<stage>_<genome>_scifiATAC.mq10.BC.rmdup.mm.tn5.bed.gz
# Usage   bash analysis/supplementary/figS7_prep_beds.sh <scratch_dir>     (~8 min, ~6 GB scratch, transient)
# Then    Rscript analysis/supplementary/figS7_fripnorm.R <scratch_dir>
# =============================================================================
#
#   genome        Pre                   wd                        nd
#   At->TAIR10    SM2_At_TAIR10         Clean.SM2v2wd_At_TAIR10   Clean.SM2v2_At_TAIR10
#   B73->B73v5    SM2_B73_B73v5         Clean.SM2v2wd_B73_B73v5   Clean.SM2v2_B73_B73v5
# wd and nd are INDEPENDENT treatments of the same raw input, NOT a chain.
# The step0_qc metadata are R-exported with row.names: data rows carry ONE MORE field than the
#   header (field 1 is a duplicated cellID rowname). field1 == field2 on every row of every
#   file; this script reads field 1 and ASSERTS that equality rather than trusting it.
# The BED barcode (column 4) is the FULL cellID and matches the metadata cellID verbatim -- no
#   barcode stripping (the combined-genome prep needs `sub -.*` because both libraries share one
#   BED; here each BED is already one library).
# =============================================================================
set -euo pipefail
SC="${1:?usage: bash figS7_prep_beds.sh <scratch_dir>}"
REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
cd "$REPO_ROOT"
DATA=data/processed/scifiATAC_B73_Arabidopsis
SOC="$DATA/socrates"
QC="$SOC/SM2v2_plate/step0_qc"
BEDDIR="$SOC/_data/_BED_files"
mkdir -p "$SC"

# METADATA VERSION -- v1 IS DELIBERATE, DO NOT "UPGRADE" IT TO v6.
# step0_qc ships updated_metadata_v1..v6. The metric values (total, pTSS, FRiP)
# are BYTE-IDENTICAL across all six -- the ladder only removes cells. v6 keeps
# ONLY cells with qc_check == 1, i.e. it has already applied a gate defined on
# the very metrics this figure plots, which would make "cleaning preserves QC"
# true by construction. v1 is the pre-gate set (all barcodes >= 200 reads,
# qc_check both values), so the comparison is honest.
#   At/Pre:  v1 1525 (1374 pass)  ->  v6 1090        B73/Pre: v1 18482 -> v6 15179
MDVER="${MDVER:-v1}"

# genome : metadata/BED suffix
for row in "TAIR10:At_TAIR10" "B73v5:B73_B73v5"; do
  IFS=: read -r g suf <<< "$row"

  pre="$QC/SM2_${suf}.minDepth200.updated_metadata_${MDVER}.txt"
  wd="$QC/Clean.SM2v2wd_${suf}.minDepth200.updated_metadata_${MDVER}.txt"
  nd="$QC/Clean.SM2v2_${suf}.minDepth200.updated_metadata_${MDVER}.txt"
  for f in "$pre" "$wd" "$nd"; do
    [ -f "$f" ] || { echo "MISSING metadata: $f" >&2; exit 1; }
  done

  # cells present in ALL THREE stages. The rowname assertion is deliberate:
  # if upstream ever writes these files without row.names, field 1 becomes
  # `total` (an integer) and this would silently build a garbage cell list.
  awk -F'\t' '
    FNR==1 { next }
    { if ($1 != $2) { print "ROWNAME ASSERT FAILED in " FILENAME " line " FNR > "/dev/stderr"; exit 1 }
      c[$1]++ }
    END { for (k in c) if (c[k] == 3) print k }
  ' "$pre" "$wd" "$nd" > "$SC/${g}_cells.txt"

  n_pre=$(( $(wc -l < "$pre") - 1 ))
  n_sh=$(wc -l < "$SC/${g}_cells.txt")
  echo "[$g] metadata $MDVER | shared cells: $n_sh  (Pre had $n_pre)"
  [ "$n_sh" -gt 0 ] || { echo "empty shared cell set for $g" >&2; exit 1; }

  # stage tag : BED basename prefix
  for srow in "pre:SM2_${suf}" "wd:Clean.SM2v2wd_${suf}" "nd:Clean.SM2v2_${suf}"; do
    IFS=: read -r stg base <<< "$srow"
    bed="$BEDDIR/${base}_scifiATAC.mq10.BC.rmdup.mm.tn5.bed.gz"
    [ -f "$bed" ] || { echo "MISSING bed: $bed" >&2; exit 1; }
    out="$SC/${g}_${stg}_reads.tsv"
    gunzip -c "$bed" \
      | awk 'NR==FNR { k[$1]=1; next } ($4 in k) { print $1"\t"$2"\t"$3"\t"$4 }' \
            "$SC/${g}_cells.txt" - > "$out"
    echo "  [$g/$stg] reads kept: $(wc -l < "$out")"
  done
done

echo "[done] prefiltered reads in $SC"
