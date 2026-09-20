#!/bin/bash
########## BATCH Lines for Resource Request ##########
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=120G
#SBATCH --job-name=Step2_cluster_plate_v5v6
#SBATCH --partition=standard
#SBATCH --output=_logs/Step2_cluster_plate_v5v6_%A_%a.log
#SBATCH --array=0-11
#
# 3_0_eval_plate_v5v6.sh  --  Step 2 (Part 2, PLATE arm), clustering on BOTH metaQC versions.
# ---------------------------------------------------------------------------
# WHY v5 AND v6
#   The v6 qc_check ("last isCell") step is size-disproportionate: it removes ~7% of the
#   small At-plate configs but only ~0.6% of the large B73-plate configs (population-normalized
#   z-cutoff, verified). We cluster BOTH the v6 (qc_check applied) and v5 (qc_check
#   dropped) cell sets so the comparison can decide whether the qc_check-removed At cells carry
#   real structure (=> drop qc_check, adopt v5) or scatter as edge noise (=> keep v6).
#   The SAME frozen config is applied to v5 and v6 within a config, so the v5-vs-v6 contrast is
#   the filter alone, not tuning.
#
# CONFIG SOURCE (borrowed, documented)
#   Uses the indep frozen configs SM2v2_indep/{B73v5,TAIR10}.cluster_config.tsv as the starting
#   point (B73v5: pcs20/k20/md0.05/minc50/res0.5 ; TAIR10: pcs8/k30/md0.3/minc50/res1). These were
#   optimized on the indep (mixed, larger) objects, so they are a REASONABLE starting point for the
#   plate arm but NOT a re-optimization. If v5 is adopted as canonical, re-run Step-1 optimization
#   for the plate cell sets (esp. At: ~1090 pure cells vs the 2937 mixed the config was tuned on).
#   The v5-vs-v6 DIFF is unaffected by config choice (identical config both sides).
#
# OUTPUT (SM2v2_plate/step2_cluster/)
#   <sample>.mQCv{5,6}.{full.SocObj_v7,updated_metadata_v7,reduced_dimensions_v7}.<tag>.{rds,txt}
#   The .mQCv5/.mQCv6 tag on the prefix is what keeps the two runs from colliding (same param tag).
#
# Array layout (At first = fast + decisive; B73 after = heavy):
#   0  SM2_At_TAIR10             mQCv5   TAIR10        6  SM2_B73_B73v5             mQCv5  B73v5
#   1  SM2_At_TAIR10             mQCv6   TAIR10        7  SM2_B73_B73v5             mQCv6  B73v5
#   2  Clean.SM2v2wd_At_TAIR10   mQCv5   TAIR10        8  Clean.SM2v2wd_B73_B73v5   mQCv5  B73v5
#   3  Clean.SM2v2wd_At_TAIR10   mQCv6   TAIR10        9  Clean.SM2v2wd_B73_B73v5   mQCv6  B73v5
#   4  Clean.SM2v2_At_TAIR10     mQCv5   TAIR10       10  Clean.SM2v2_B73_B73v5     mQCv5  B73v5
#   5  Clean.SM2v2_At_TAIR10     mQCv6   TAIR10       11  Clean.SM2v2_B73_B73v5     mQCv6  B73v5
#
#   sbatch 0_scripts/part2_split/step2_cluster/3_0_eval_plate_v5v6.sh                # all 12
#   sbatch --array=0-5 0_scripts/part2_split/step2_cluster/3_0_eval_plate_v5v6.sh    # At decision arm only (fast)
#
# PREREQ: 1_QC_scifiATAC_SM2v2_plate.sh done (raw.soc.rds + v5 + v6 for all 6 configs).
# ---------------------------------------------------------------------------

set -euo pipefail
PROJECT_ROOT="${PROJECT_ROOT:-.}"
BASE="${PROJECT_ROOT}/6_socrates"
SCRIPTS="${SCRIPTS:-${BASE}/0_scripts}"
cd "${BASE}"
# conda activate ambientmapper-manuscript   (environment.yml at the repo root)

OUTDIR="SM2v2_plate"

# sample | metaQC version (5 or 6)
TASKS=(
  "SM2_At_TAIR10|5"             "SM2_At_TAIR10|6"
  "Clean.SM2v2wd_At_TAIR10|5"   "Clean.SM2v2wd_At_TAIR10|6"
  "Clean.SM2v2_At_TAIR10|5"     "Clean.SM2v2_At_TAIR10|6"
  "SM2_B73_B73v5|5"             "SM2_B73_B73v5|6"
  "Clean.SM2v2wd_B73_B73v5|5"   "Clean.SM2v2wd_B73_B73v5|6"
  "Clean.SM2v2_B73_B73v5|5"     "Clean.SM2v2_B73_B73v5|6"
)

entry="${TASKS[$SLURM_ARRAY_TASK_ID]:-}"
[[ -n "$entry" ]] || { echo "ERROR: no task for array index ${SLURM_ARRAY_TASK_ID}"; exit 1; }
IFS='|' read -r SAMPLE VER <<< "$entry"

GENOME="${SAMPLE##*_}"                                  # trailing token -> B73v5 | TAIR10
OUTPREFIX="${SAMPLE}.mQCv${VER}"

SOC="${OUTDIR}/step0_qc/${SAMPLE}.raw.soc.rds"
META="${OUTDIR}/step0_qc/${SAMPLE}.minDepth200.updated_metadata_v${VER}.txt"
CFG="SM2v2_indep/${GENOME}.cluster_config.tsv"

[[ -s "$SOC"  ]] || { echo "ERROR: missing SocObj: $SOC";       exit 1; }
[[ -s "$META" ]] || { echo "ERROR: missing v${VER} metadata: $META"; exit 1; }
[[ -s "$CFG"  ]] || { echo "ERROR: missing borrowed config: $CFG";   exit 1; }

# pull frozen params by column name (robust to column order)
read -r PCS K MD MINC RES < <(awk -F'\t' 'NR==1{for(x=1;x<=NF;x++)h[$x]=x}
  NR==2{print $h["pcs"], $h["k_near"], $h["min_dist"], $h["min_c"], $h["resolution"]}' "$CFG")
[[ -n "$PCS" && -n "$K" && -n "$MD" && -n "$MINC" && -n "$RES" ]] || { echo "ERROR: could not parse config $CFG"; exit 1; }

echo "$(date): Step-2 plate cluster ${OUTPREFIX} (${GENOME}) | v${VER} | pcs=${PCS} k=${K} min_dist=${MD} min.c=${MINC} res=${RES}"
Rscript "${SCRIPTS}/part2_split/step2_cluster/3_0_0b_cluster_pergenome.R" \
  "$SOC" "$META" "$OUTDIR" "$OUTPREFIX" "$PCS" "$K" "$MD" "$MINC" "$RES" 1 "$GENOME"

echo "$(date): DONE ${OUTPREFIX} -> ${OUTDIR}/step2_cluster/${OUTPREFIX}.updated_metadata_v7.*"
