#!/bin/bash
########## BATCH Lines for Resource Request ##########
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --job-name=Split_BED_by_plate
#SBATCH --partition=standard
#SBATCH --output=_logs/Split_BED_by_plate_%A_%a.log
#SBATCH --array=0-5
#
# 0_09_split_beds_by_plate.sh
# ---------------------------------------------------------------------------
# PLATE-SPLIT ("normal user") arm for SM2v2 — BED preparation.
#
# WHY THIS EXISTS
#   The SM2v2_indep objects map EVERY barcode to BOTH references (the substrate
#   AmbientMapper needs to detect ambient). That is NOT what a scifi-ATAC user
#   does downstream: they split samples per plate and map each to its EXPECTED
#   genome. This arm reproduces that normal workflow, Pre vs Post-AM, so the
#   practical value of cleaning can be measured in the deployment scenario.
#
#     At  plate  ->  TAIR10   (its expected genome)
#     B73 plate  ->  B73v5    (its expected genome)
#
#   Cross-plate barcodes (At-plate reads on B73v5, B73-plate reads on TAIR10)
#   are DISCARDED here by design — a normal user never sees them. They remain
#   the subject of the multi-reference arm (SM2v2_indep), which this does NOT
#   replace: the Step-4 cl5 contamination finding lives there and CANNOT be
#   reproduced here, because cl5 *is* the At-plate-in-maize population.
#
# SCOPE NOTE
#   This is a read-level filter only. Nothing is re-mapped; the plate tag is
#   already carried in the BED barcode column, e.g.
#     CGCGCAACAAGTCTGTAGTAAGTTTC-SM2_At_AraTAIR10_scifiATAC
#     CATTGGACACGACGAACGTATACTAA-SM2_B73_Zm_B73v5_scifiATAC
#   so the split is an exact substring match on field 4 -- no ambiguity.
#
# OUTPUT (into _data/_BED_files/, alongside the existing indep BEDs)
#   SM2_At_TAIR10_*            SM2_B73_B73v5_*             (PreClean)
#   Clean.SM2v2_At_TAIR10_*    Clean.SM2v2_B73_B73v5_*     (Post nd)
#   Clean.SM2v2wd_At_TAIR10_*  Clean.SM2v2wd_B73_B73v5_*   (Post wd, headline)
#
#   The names are chosen so the downstream driver's species parser
#   (species="${name##*_}") still resolves to B73v5 / TAIR10 unchanged.
#
# Usage:  sbatch 0_scripts/part2_split/step0_qc/0_09_split_beds_by_plate.sh
#         Submit from the 6_socrates dir so _logs/ resolves.
#         Then: sbatch 0_scripts/part2_split/step0_qc/1_QC_scifiATAC_SM2v2_plate.sh
# ---------------------------------------------------------------------------

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-.}"
BASE="${PROJECT_ROOT}/6_socrates"
BEDDIR="${BASE}/_data/_BED_files"
SUF="_scifiATAC.mq10.BC.rmdup.mm.tn5.bed.gz"

# src_prefix | plate_tag (literal substring of field 4) | out_prefix
TASKS=(
  "SM2_TAIR10|-SM2_At_|SM2_At_TAIR10"
  "SM2_B73v5|-SM2_B73_|SM2_B73_B73v5"
  "Clean.SM2v2_TAIR10|-SM2_At_|Clean.SM2v2_At_TAIR10"
  "Clean.SM2v2_B73v5|-SM2_B73_|Clean.SM2v2_B73_B73v5"
  "Clean.SM2v2wd_TAIR10|-SM2_At_|Clean.SM2v2wd_At_TAIR10"
  "Clean.SM2v2wd_B73v5|-SM2_B73_|Clean.SM2v2wd_B73_B73v5"
)

idx="${SLURM_ARRAY_TASK_ID:-${1:-}}"
[[ -n "$idx" ]] || { echo "ERROR: no array index (sbatch --array, or pass 0-5 as \$1)"; exit 1; }
entry="${TASKS[$idx]:-}"
[[ -n "$entry" ]] || { echo "ERROR: no task for index $idx"; exit 1; }

IFS='|' read -r src_pfx tag out_pfx <<< "$entry"
src="${BEDDIR}/${src_pfx}${SUF}"
dst="${BEDDIR}/${out_pfx}${SUF}"

[[ -f "$src" ]] || { echo "ERROR: missing source BED $src"; exit 2; }

# pigz if available (these BEDs run to ~433 MB gz); gzip otherwise.
if command -v pigz >/dev/null 2>&1; then
  DECOMP=(pigz -dc); COMP=(pigz -c -p "${SLURM_CPUS_PER_TASK:-4}")
else
  DECOMP=(gzip -cd);  COMP=(gzip -c)
fi

echo "$(date): plate-split  ${src_pfx}  --[keep ${tag}]-->  ${out_pfx}"

# index() = literal substring match on the barcode field; no regex escaping games.
"${DECOMP[@]}" "$src" \
  | awk -F'\t' -v t="$tag" 'index($4, t)' \
  | "${COMP[@]}" > "${dst}.tmp"
mv -f "${dst}.tmp" "$dst"

# ---- report: reads + unique barcodes kept, and what was dropped --------------
kept_reads=$("${DECOMP[@]}" "$dst" | wc -l)
kept_bcs=$("${DECOMP[@]}" "$dst" | cut -f4 | sort -u | wc -l)
src_reads=$("${DECOMP[@]}" "$src" | wc -l)

echo "     source reads : ${src_reads}"
echo "     kept   reads : ${kept_reads}  ($(awk -v a="$kept_reads" -v b="$src_reads" 'BEGIN{printf "%.1f", (b?100*a/b:0)}')% of source)"
echo "     kept barcodes: ${kept_bcs}"
echo "     -> ${dst}"
echo "$(date): DONE ${out_pfx}"
