#!/bin/bash
########## BATCH Lines for Resource Request ##########
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --job-name=Step4_smooth_SM2v2_plate
#SBATCH --partition=standard
#SBATCH --output=_logs/Step4_smooth_SM2v2_plate_%A_%a.log
#SBATCH --array=0-2
#
# 4_3f_smooth_plate.sh -- Step-4 STEP 2 (plate arm): kNN-smooth the gene activity.
# ---------------------------------------------------------------------------
# Modular split (mirrors clustering: opt -> cluster -> compare):
#   [1] gene activity  = 4_0 (done, *.plate.genes.sparse.rds)
#   [2] SMOOTH         = THIS (Markov diffusion in v7 PC space) -> *.plate.genes.smoothed.rds  (reusable)
#   [3] annotate       = 4_3g_annotate_plate.sh (on the smoothed artifact; re-tunable w/o re-smoothing)
#
#   0: SM2_At_TAIR10   Pre    3: SM2_B73_B73v5           Pre
#   1: Clean.SM2v2wd_At_TAIR10 wd   4: Clean.SM2v2wd_B73_B73v5 wd
#   2: Clean.SM2v2_At_TAIR10   nd   5: Clean.SM2v2_B73_B73v5   nd
#   sbatch --array=0-2 ...  # At (decisive) -- default ;  sbatch --array=3-5 ...  # B73 (heavier: 15k x 32k dense ~3.8 GB)
#
# Prereq: 4_0 gene matrices built (4_3_eval_plate.sh step [1/2]) + Step-2 v7 reduced_dimensions present.
# ---------------------------------------------------------------------------

# conda activate ambientmapper-manuscript   (environment.yml at the repo root)
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-.}"
BASE="${PROJECT_ROOT}/6_socrates"
SCRIPTS="${SCRIPTS:-${BASE}/0_scripts}"
STEP2="${BASE}/SM2v2_plate/step2_cluster"
OUT="${BASE}/SM2v2_plate/step3_compare"
s4_3f="${SCRIPTS}/common/4_3f_smooth_gene_activity.R"
K=25; STEP=3
NORM="${GENE_LENGTH_NORM:-raw}"        # raw (default) | perkb -> consume the length-normalized matrix
case "$NORM" in raw) VAR="" ;; perkb) VAR=".perkb" ;; *) echo "ERROR: GENE_LENGTH_NORM must be raw|perkb"; exit 1 ;; esac

declare -a NAMES=( SM2_At_TAIR10 Clean.SM2v2wd_At_TAIR10 Clean.SM2v2_At_TAIR10 SM2_B73_B73v5 Clean.SM2v2wd_B73_B73v5 Clean.SM2v2_B73_B73v5 )
declare -a STAGES=( Pre wd nd Pre wd nd )
i="${SLURM_ARRAY_TASK_ID:?set via --array}"
name="${NAMES[$i]:?bad array index $i}"; stage="${STAGES[$i]}"

mat="${OUT}/${name}.plate${VAR}.genes.sparse.rds"
rd="$(ls ${STEP2}/${name}.mQCv6.reduced_dimensions_v7.*.txt 2>/dev/null | head -1)"
[[ -f "$mat" ]]                || { echo "ERROR: missing gene matrix $mat (run 4_3_eval_plate.sh first)"; exit 2; }
[[ -n "${rd:-}" && -f "$rd" ]] || { echo "ERROR: no v7 reduced_dimensions for $name in $STEP2"; exit 2; }

echo "$(date): Step-4 SMOOTH | $name ($stage norm=$NORM) | matrix=$(basename "$mat") rd=$(basename "$rd") k=$K step=$STEP"
Rscript "$s4_3f" "$OUT" "${name}.plate${VAR}" "$mat" "$rd" "$stage" "$K" "$STEP"
echo "$(date): DONE $name -> ${OUT}/${name}.plate${VAR}.genes.smoothed.rds"
