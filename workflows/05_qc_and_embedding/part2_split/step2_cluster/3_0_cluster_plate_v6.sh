#!/bin/bash
########## BATCH Lines for Resource Request ##########
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=120G
#SBATCH --job-name=Step2_cluster_plate_v6
#SBATCH --partition=standard
#SBATCH --output=_logs/Step2_cluster_plate_v6_%A_%a.log
#SBATCH --array=0-5
#
# 3_0_cluster_plate_v6.sh  --  Step 2 (Part 2, PLATE arm): FINAL clustering on the FROZEN
#                              plate configs, v6 only. Supersedes 3_0_eval_plate_v5v6.sh.
# ---------------------------------------------------------------------------
# WHY THIS (vs 3_0_eval_plate_v5v6.sh)
#   The v5-vs-v6 question is settled (v6 canonical) and the borrowed indep configs
#   were wrong for the plate cell sets (At shattered to 94 clusters). This driver applies the
#   RE-OPTIMIZED, FROZEN plate configs and runs v6 only:
#     SM2v2_plate/{TAIR10,B73v5}.plate.cluster_config.tsv
#     At  (TAIR10): pcs5/k20/md0.3/minc50/res0.3/m.clst40  -> 5 clusters (chosen from the 3_0_1b/1c panels;
#                   continuum bins, markers adjudicate). m.clst=40 is NON-default -> Step-2 R now
#                   reads m_clst from the config and tags the output .mclst_40 (no collision).
#     B73 (B73v5):  pcs20/k20/md0.05/minc50/res0.3/m.clst50 -> ~5 clusters (scan peak-stability).
#   Same frozen config applied to Pre / wd / nd within a genome (plan decision B) so any
#   Pre-vs-Post difference is cleaning, not tuning.
#
# Array layout (At first = fast; B73 after = heavy), all v6:
#   0  SM2_At_TAIR10              (Pre)      3  SM2_B73_B73v5              (Pre)
#   1  Clean.SM2v2wd_At_TAIR10    (wd)       4  Clean.SM2v2wd_B73_B73v5    (wd)
#   2  Clean.SM2v2_At_TAIR10      (nd)       5  Clean.SM2v2_B73_B73v5      (nd)
#
#   sbatch --array=0-2 0_scripts/part2_split/step2_cluster/3_0_cluster_plate_v6.sh   # At (frozen)  -- fast
#   sbatch --array=3-5 0_scripts/part2_split/step2_cluster/3_0_cluster_plate_v6.sh   # B73 (confirm config first)
#   sbatch            0_scripts/part2_split/step2_cluster/3_0_cluster_plate_v6.sh   # all 6
#
# OUTPUT: SM2v2_plate/step2_cluster/<sample>.mQCv6.{full.SocObj_v7,updated_metadata_v7,reduced_dimensions_v7}<tag>.{rds,txt}
#   At tag carries .mclst_40 -> supersedes (does not overwrite) the broken borrowed-config .mQCv6 objects.
# PREREQ: 1_QC_scifiATAC_SM2v2_plate.sh done; plate configs frozen.
# ---------------------------------------------------------------------------

set -euo pipefail
PROJECT_ROOT="${PROJECT_ROOT:-.}"
BASE="${PROJECT_ROOT}/6_socrates"
SCRIPTS="${SCRIPTS:-${BASE}/0_scripts}"
cd "${BASE}"
# conda activate ambientmapper-manuscript   (environment.yml at the repo root)

OUTDIR="SM2v2_plate"
VER=6

# v6-only samples: At Pre/wd/nd, then B73 Pre/wd/nd
SAMPLES=(
  "SM2_At_TAIR10"            "Clean.SM2v2wd_At_TAIR10"   "Clean.SM2v2_At_TAIR10"
  "SM2_B73_B73v5"            "Clean.SM2v2wd_B73_B73v5"   "Clean.SM2v2_B73_B73v5"
)
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]:-}"
[[ -n "$SAMPLE" ]] || { echo "ERROR: no sample for array index ${SLURM_ARRAY_TASK_ID}"; exit 1; }

GENOME="${SAMPLE##*_}"                                  # trailing token -> B73v5 | TAIR10
OUTPREFIX="${SAMPLE}.mQCv${VER}"

SOC="${OUTDIR}/step0_qc/${SAMPLE}.raw.soc.rds"
META="${OUTDIR}/step0_qc/${SAMPLE}.minDepth200.updated_metadata_v${VER}.txt"
CFG="${OUTDIR}/${GENOME}.plate.cluster_config.tsv"

[[ -s "$SOC"  ]] || { echo "ERROR: missing SocObj: $SOC";        exit 1; }
[[ -s "$META" ]] || { echo "ERROR: missing v${VER} metadata: $META"; exit 1; }
[[ -s "$CFG"  ]] || { echo "ERROR: missing plate config: $CFG";   exit 1; }

# pull frozen params by column name (robust to order); m_clst defaults to 50 if the column is absent
read -r PCS K MD MINC RES MCLST < <(awk -F'\t' 'NR==1{for(x=1;x<=NF;x++)h[$x]=x}
  NR==2{mc=(("m_clst" in h)?$h["m_clst"]:"50"); print $h["pcs"], $h["k_near"], $h["min_dist"], $h["min_c"], $h["resolution"], mc}' "$CFG")
[[ -n "$PCS" && -n "$K" && -n "$MD" && -n "$MINC" && -n "$RES" && -n "$MCLST" ]] || { echo "ERROR: could not parse config $CFG"; exit 1; }

echo "$(date): Step-2 plate FINAL cluster ${OUTPREFIX} (${GENOME}) | v${VER} | pcs=${PCS} k=${K} min_dist=${MD} min.c=${MINC} res=${RES} m.clst=${MCLST}"
Rscript "${SCRIPTS}/part2_split/step2_cluster/3_0_0b_cluster_pergenome.R" \
  "$SOC" "$META" "$OUTDIR" "$OUTPREFIX" "$PCS" "$K" "$MD" "$MINC" "$RES" 1 "$GENOME" "$MCLST"

echo "$(date): DONE ${OUTPREFIX} -> ${OUTDIR}/step2_cluster/${OUTPREFIX}.updated_metadata_v7.*"