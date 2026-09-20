#!/usr/bin/env bash
# =============================================================================
# browser/split.sh  -  STEP 2 of the Fig 5 panel L genome-browser chain: split each plate tn5 BED
#   into per-group BED.gz (group = consensus type call), producing the inputs bw.sh consumes.
# Inputs  $DATA/socrates/_data/_BED_files/<prefix>_scifiATAC.mq10.BC.rmdup.mm.tn5.bed.gz  (6 arms)
#         $DATA/socrates/_data/_GenomeInfo/SM2v2_{B73v5,TAIR10}.chrs.size.txt
#         figures/main/fig5/browser/groups/<prefix>.cell2group.tsv        (browser/groups.R)
# Output  figures/main/fig5/browser/bed_by_type/<prefix>.<group>.bed.gz (+ ALLcells), browser/split_counts.tsv
# Usage   bash analysis/fig5_biological_impact/browser/split.sh                    (all six arms)
#         bash analysis/fig5_biological_impact/browser/split.sh SM2_At_TAIR10      (one arm)
#   Needs only gzip/awk/sort (pigz used if present). At arms ~1 min each; maize arms ~5-10 min.
# =============================================================================
#
# One streaming pass per arm (awk hash on BED column 4 = the FULL cellID, verbatim, as in
# figS7_prep_beds.sh). Reads are kept only if (a) the barcode belongs to a grouped cell and
# (b) the chrom is in the genome's size file -- the size-file whitelist guarantees bedtools
# genomecov never sees an unknown contig. Every kept read is written TWICE: to its group's file
# and to the arm's ALLcells aggregate.
# Outputs raw (unslopped, input-ordered) tn5 sites: the slop / sort / RPM happen in bw.sh.
# =============================================================================
set -euo pipefail

REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
cd "$REPO_ROOT"
DATA=data/processed/scifiATAC_B73_Arabidopsis
BEDDIR="$DATA/socrates/_data/_BED_files"
GINFO="$DATA/socrates/_data/_GenomeInfo"
BROWSER=figures/main/fig5/browser
# NOT `GROUPS` -- that is a reserved bash array (the user's group IDs);
#   assignments to it are silently ignored and $GROUPS expanded to "20".
GRPDIR="$BROWSER/groups"
OUTDIR="$BROWSER/bed_by_type"
LOG="$BROWSER/split_counts.tsv"
mkdir -p "$OUTDIR"

GZ="gzip"; command -v pigz >/dev/null 2>&1 && GZ="pigz"

# arm prefix -> BED basename + chrom sizes.
# No associative arrays: macOS /bin/bash is 3.2 (`declare -A` is unavailable).
#   case-functions are portable to both machines.
KNOWN_ARMS="SM2_B73_B73v5 Clean.SM2v2wd_B73_B73v5 Clean.SM2v2_B73_B73v5 SM2_At_TAIR10 Clean.SM2v2wd_At_TAIR10 Clean.SM2v2_At_TAIR10"
bed_of()  {
  case " $KNOWN_ARMS " in
    *" $1 "*) echo "$BEDDIR/${1}_scifiATAC.mq10.BC.rmdup.mm.tn5.bed.gz" ;;
    *) echo "unknown arm: $1" >&2; return 1 ;;
  esac
}
size_of() {
  case "$1" in
    *At_TAIR10) echo "$GINFO/SM2v2_TAIR10.chrs.size.txt" ;;
    *B73_B73v5) echo "$GINFO/SM2v2_B73v5.chrs.size.txt"  ;;
    *) echo "unknown genome in arm: $1" >&2; return 1 ;;
  esac
}

ARMS=("$@")
if [ ${#ARMS[@]} -eq 0 ]; then
  ARMS=(SM2_B73_B73v5 Clean.SM2v2wd_B73_B73v5 Clean.SM2v2_B73_B73v5
        SM2_At_TAIR10 Clean.SM2v2wd_At_TAIR10 Clean.SM2v2_At_TAIR10)
fi

[ -f "$LOG" ] || printf "prefix\tgroup\tn_reads\n" > "$LOG"

for pf in "${ARMS[@]}"; do
  bed=$(bed_of "$pf")
  siz=$(size_of "$pf")
  c2g="$GRPDIR/${pf}.cell2group.tsv"
  for f in "$bed" "$siz" "$c2g"; do
    [ -f "$f" ] || { echo "MISSING: $f (run browser/groups.R first?)" >&2; exit 1; }
  done

  tmp="$OUTDIR/_tmp_${pf}"
  rm -rf "$tmp"; mkdir -p "$tmp"
  echo "$(date '+%H:%M:%S')  [$pf] splitting $(basename "$bed")"

  # one pass: barcode -> group hash, chrom whitelist, per-group writers
  "$GZ" -dc "$bed" | awk -F'\t' -v OFS='\t' \
      -v C2G="$c2g" -v SIZ="$siz" -v TMP="$tmp" -v PF="$pf" -v LOG="$LOG" '
    BEGIN {
      while ((getline line < C2G) > 0) {
        split(line, a, "\t"); grp[a[1]] = a[2]
      }
      close(C2G)
      while ((getline line < SIZ) > 0) {
        split(line, a, "\t"); okchr[a[1]] = 1
      }
      close(SIZ)
    }
    {
      if (!($1 in okchr)) { drop_chr++; next }
      if (!($4 in grp))   { drop_bc++;  next }
      g = grp[$4]
      print $0 > (TMP "/" g ".bed")
      print $0 > (TMP "/ALLcells.bed")
      n[g]++; n["ALLcells"]++
    }
    END {
      for (g in n) printf "%s\t%s\t%d\n", PF, g, n[g] >> LOG
      printf "  kept-to-ALL %d | dropped: %d off-whitelist chrom, %d non-member barcode\n",
             n["ALLcells"] + 0, drop_chr + 0, drop_bc + 0 > "/dev/stderr"
    }'

  # compress into place, explicit filenames, verify each
  for f in "$tmp"/*.bed; do
    g=$(basename "$f" .bed)
    out="$OUTDIR/${pf}.${g}.bed.gz"
    "$GZ" -c "$f" > "$out"
    [ -s "$out" ] || { echo "FAILED to write $out" >&2; exit 1; }
  done
  nfiles=$(ls "$tmp"/*.bed | wc -l | tr -d ' ')
  rm -rf "$tmp"
  echo "$(date '+%H:%M:%S')  [$pf] wrote $nfiles group BED.gz -> $OUTDIR/"
done

echo
echo "split complete. Per-group read counts: $LOG"
echo "next (HPC, from the repo root): sbatch analysis/fig5_biological_impact/browser/bw.sh"
