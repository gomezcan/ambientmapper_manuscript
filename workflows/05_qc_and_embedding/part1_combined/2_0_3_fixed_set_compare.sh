#!/usr/bin/env bash
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --job-name=fixedSetCompare
#SBATCH --partition=standard
#SBATCH --output=_logs/%x-%j.log
#
# 2_0_3_fixed_set_compare.sh — fixed-barcode pre/post comparison for SM2v2.
# Runs 2_0_3_fixed_set_compare.R for both species (At, B73). Metadata-only:
# reads the 1_2 outputs (<pool>.updated_metadata.txt), which already exist and
# are independent of the 1_3 depth cascade, so this can run any time after 1_2.
#
#   pre  : SM2/step0_qc/<pool>.updated_metadata.txt
#   post : SM2v2_clean/step0_qc/Clean.SM2v2_<species>.updated_metadata.txt
#   out  : compare/<pool>.fixedSet.minDepth<DEPTH>.pre_post.txt (+ .thresholds.txt)
#
# Optional: DEPTH=200 (raw-depth cell-call floor; matches the 1_3 default).

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-.}"
BASE="${PROJECT_ROOT}/6_socrates"
SCRIPTS="${SCRIPTS:-${BASE}/0_scripts}"
DEPTH="${DEPTH:-200}"
OUTDIR="${BASE}/compare"

# conda activate ambientmapper-manuscript   (environment.yml at the repo root)

mkdir -p "$OUTDIR" "${BASE}/_logs"
cd "$BASE"

# pool : pre_dir : post_prefix
samples=(
  "SM2_At:SM2:Clean.SM2v2_At"
  "SM2_B73:SM2:Clean.SM2v2_B73"
)

for s in "${samples[@]}"; do
  pool="${s%%:*}"; rest="${s#*:}"; predir="${rest%%:*}"; postpref="${rest##*:}"
  pre="${BASE}/${predir}/step0_qc/${pool}.updated_metadata.txt"
  post="${BASE}/SM2v2_clean/step0_qc/${postpref}.updated_metadata.txt"
  [[ -s "$pre"  ]] || { echo "ERROR: missing pre metadata: $pre";  exit 1; }
  [[ -s "$post" ]] || { echo "ERROR: missing post metadata: $post"; exit 1; }
  echo " - comparing ${pool}  (DEPTH=${DEPTH})"
  echo "     pre  = $pre"
  echo "     post = $post"
  Rscript "${SCRIPTS}/part1_combined/2_0_3_fixed_set_compare.R" "$pre" "$post" "$pool" "$DEPTH" "$OUTDIR"
done

echo " - done; outputs in ${OUTDIR}/"
ls -lh "${OUTDIR}"/*.fixedSet.* 2>/dev/null || true
