#!/bin/bash
########## BATCH Lines for Resource Request ##########
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --job-name=P1indep_3_1_coembed
#SBATCH --partition=standard
#SBATCH --output=_logs/P1indep_3_1_coembed_%A_%a.log
#SBATCH --array=0-1
#
# 3_1_coembed_build.sh -- STAGE 1: build the multi-reference co-embedding object.
# ---------------------------------------------------------------------------
#   Task 0: Pre  (raw)  SM2_B73v5           + SM2_TAIR10           -> SM2v2_coembed_Pre
#   Task 1: Post (wd)   Clean.SM2v2wd_B73v5 + Clean.SM2v2wd_TAIR10 -> SM2v2_coembed_Post_wd
#
# Stacks the two per-genome step2 tile matrices (union of v7 cells, zero-fill) into ONE
# joint object; Genome = plate-of-origin ground truth. Pre alone replaces Fig 1F/G/H;
# Post is the cleaning-effect counterpart (Fig 5 territory). Needs only the Matrix package.
# ---------------------------------------------------------------------------
set -euo pipefail
# conda activate ambientmapper-manuscript   (environment.yml at the repo root)

PROJECT_ROOT="${PROJECT_ROOT:-.}"
BASE="${PROJECT_ROOT}/6_socrates"
cd "$BASE"
S2="${BASE}/SM2v2_indep/step2_cluster"
OUT="${BASE}/SM2v2_indep/coembed"
SCR="${SCRIPTS:-${BASE}/0_scripts}/part1_indep"
mkdir -p "$OUT"

B73_TAG="pcs_20.k_near_20.min_dis_0.05.minc_50.res_0.5"
AT_TAG="pcs_8.k_near_30.min_dis_0.3.minc_50.res_1"

i="${SLURM_ARRAY_TASK_ID:-${1:?set via --array or pass 0/1}}"
case "$i" in
  0) STAGE="Pre";     B73_RDS="${S2}/SM2_B73v5.full.SocObj_v7.${B73_TAG}.rds";
                      AT_RDS="${S2}/SM2_TAIR10.full.SocObj_v7.${AT_TAG}.rds";
                      PREFIX="${OUT}/SM2v2_coembed_Pre" ;;
  1) STAGE="Post_wd"; B73_RDS="${S2}/Clean.SM2v2wd_B73v5.full.SocObj_v7.${B73_TAG}.rds";
                      AT_RDS="${S2}/Clean.SM2v2wd_TAIR10.full.SocObj_v7.${AT_TAG}.rds";
                      PREFIX="${OUT}/SM2v2_coembed_Post_wd" ;;
  *) echo "bad array index $i"; exit 2 ;;
esac
for f in "$B73_RDS" "$AT_RDS"; do [[ -s "$f" ]] || { echo "MISSING: $f"; exit 2; }; done

echo "$(date): [$STAGE] co-embed build"
if [[ -s "${PREFIX}.coembed.soc.rds" ]]; then
  echo "  exists, skipping (rm to rebuild): ${PREFIX}.coembed.soc.rds"
else
  Rscript "${SCR}/3_1_coembed_build.R" "$B73_RDS" "$AT_RDS" "$PREFIX" "$STAGE"
fi
echo "$(date): DONE [$STAGE] -> ${PREFIX}.coembed.soc.rds + .coembed.meta.tsv"
