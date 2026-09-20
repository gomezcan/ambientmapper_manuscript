#!/bin/bash
########## BATCH Lines for Resource Request ##########
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --job-name=Step5_1_acr_matrix_plate
#SBATCH --partition=standard
#SBATCH --output=_logs/Step5_1_acr_matrix_plate_%A_%a.log
#SBATCH --array=0-1
#
# 5_1_build_acr_matrix_plate.sh -- Step 5.1 (plate arm): consensus-ACR x cell matrices + the
#                                  portable SEACells bundle.
# ---------------------------------------------------------------------------
#   0: TAIR10 (At)      1: B73v5 (B73)
#
# ONE TASK PER GENOME, NOT PER OBJECT -- deliberate. The feature set must be FROZEN across stages,
# so Pre has to run FIRST and emit .acrs.features.txt, which wd and nd then consume as
# FEATURE_WHITELIST. A 6-way array cannot guarantee that ordering; three sequential calls inside one
# task can. (Measured: filtering per object gave At 7105/6917/8684 features -- a 25% swing
# driven purely by nd's ~65% deeper cells, which would make the meta-cell distance metric a function
# of stage. Same confound the frozen peak set exists to prevent.)
#
# Prereq: Step 5.0 consensus peaks; Step-2 v7 metadata + reduced dims; plate tn5 BEDs.
# Emits per object:  .acrs.sparse.rds  .acrs.features.txt  .acrs.stats.tsv
#                    .seacells.{mtx.gz,svd.tsv,obs.tsv,var.tsv}
#
# Env knobs:
#   MIN_CELLS_FRAC (0.005)  MIN_CELLS_ABS (10)
#     Fractional BY DESIGN. Socrates' cleanData(min.c=50) is absolute: measured, t=50 retains
#     157/39,923 At ACRs (0.4%) vs 42,706/284,132 maize (15.0%) -- a 272x difference at the SAME
#     threshold. NOTE the absolute floor BINDS for At (10 > 0.005*1061=5.3), so At sits at 0.94% of
#     cells vs maize 0.51% -- i.e. ~1.8x STRICTER on At. MIN_CELLS_ABS=5 would give At ~20,500 ACRs
#     instead of 7,105. This is a bandwidth-like parameter: DECLARE IT AND SWEEP IT.
#   PEAKSET (consensus_Pre)  -> set to union_PreWdNd for the reverse-direction supplementary arm.
#
# Local run (no SLURM; ~1 min At, ~10 min maize):
#   bash 0_scripts/part2_split/step5_metacell/5_1_build_acr_matrix_plate.sh 0
# ---------------------------------------------------------------------------

# conda activate ambientmapper-manuscript   (environment.yml at the repo root)
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-.}"
BASE="${PROJECT_ROOT}/6_socrates"
SCRIPTS="${SCRIPTS:-${BASE}/0_scripts}"
PEAKS="${BASE}/SM2v2_plate/_consensus_peaks"
STEP2="${BASE}/SM2v2_plate/step2_cluster"
BEDDIR="${BASE}/_data/_BED_files"
OUT="${BASE}/SM2v2_plate/step5_metacell"
s5_1="${SCRIPTS}/common/5_1_build_acr_matrix.R"
s5_1b="${SCRIPTS}/common/5_1b_export_for_seacells.R"

SUF="_scifiATAC.mq10.BC.rmdup.mm.tn5.bed.gz"
TAG="mQCv6"
PEAKSET="${PEAKSET:-consensus_Pre}"
# SWEPT AND SET (feature-filter sweep, not shipped): t=5 is optimal for BOTH genomes.
# Adjusted kNN purity in an independent ACR-LSI space, vs the frozen Leiden labels:
#   At    t=1 .331 | t=5 .362 | t=10 .361 | t=20 .359 | t=50 .312     -> flat 5-20, peak at 5
#   maize t=1 .660 | t=5 .663 | t=10 .662 | t=20 .645 | t=50 .607 | t=100 .535  -> declines past 10
# The earlier fractional default was the WRONG FIX: a feature's information content depends on how
# many times it was OBSERVED, not on what fraction of the dataset that is, so a low ABSOLUTE floor
# serves both genomes (t=5 = 0.47% of At cells but 0.035% of maize cells -- wildly different
# fractions, same optimum). cleanData(min.c=50) was bad because it was TOO HIGH, not because it was
# absolute. MIN_CELLS_FRAC=0 disables the fractional term.
export MIN_CELLS_FRAC="${MIN_CELLS_FRAC:-0}"
export MIN_CELLS_ABS="${MIN_CELLS_ABS:-5}"

# frozen Step-2 configs (see SM2v2_plate/{TAIR10,B73v5}.plate.cluster_config.tsv)
CFG_TAIR10="pcs_5.k_near_20.min_dis_0.3.minc_50.res_0.3.mclst_40"
CFG_B73v5="pcs_20.k_near_20.min_dis_0.05.minc_50.res_0.3"

declare -a GENOMES=( TAIR10 B73v5 )
declare -a SPECIES=( At     B73   )
i="${SLURM_ARRAY_TASK_ID:-${1:?set via --array or pass 0/1}}"
g="${GENOMES[$i]:?bad array index $i}"; sp="${SPECIES[$i]}"
case "$g" in TAIR10) cfg="$CFG_TAIR10" ;; B73v5) cfg="$CFG_B73v5" ;; esac

PEAKBED="${PEAKS}/${g}.${PEAKSET}.bed"
[[ -s "$PEAKBED" ]] || { echo "ERROR: missing peak set $PEAKBED -- run 5_0 first"; exit 2; }
mkdir -p "$OUT"

# Pre MUST be first: it emits the feature whitelist the other two consume.
declare -a NAMES=( "SM2_${sp}_${g}" "Clean.SM2v2wd_${sp}_${g}" "Clean.SM2v2_${sp}_${g}" )
PRE="${NAMES[0]}"
WL="${OUT}/${PRE}.plate.acrs.features.txt"

echo "$(date): Step 5.1 ${g} | peakset=${PEAKSET} | MIN_CELLS_ABS=${MIN_CELLS_ABS} MIN_CELLS_FRAC=${MIN_CELLS_FRAC}"

for k in 0 1 2; do
  name="${NAMES[$k]}"
  bed="${BEDDIR}/${name}${SUF}"
  meta="${STEP2}/${name}.${TAG}.updated_metadata_v7.${cfg}.txt"
  rdim="${STEP2}/${name}.${TAG}.reduced_dimensions_v7.${cfg}.txt"
  for f in "$bed" "$meta" "$rdim"; do
    [[ -s "$f" ]] || { echo "ERROR: missing input $f"; exit 3; }
  done

  echo "--- [${k}] ${name}"
  if [[ $k -eq 0 ]]; then
    unset FEATURE_WHITELIST || true          # Pre DEFINES the feature space
  else
    [[ -s "$WL" ]] || { echo "ERROR: whitelist $WL not produced by the Pre run"; exit 4; }
    export FEATURE_WHITELIST="$WL"
  fi
  Rscript "$s5_1"  "${OUT}/${name}.plate" "$PEAKBED" "$bed" "$meta"
  Rscript "$s5_1b" "${OUT}/${name}.plate.acrs.sparse.rds" "$rdim" "$meta" "${OUT}/${name}.seacells"
done

echo "$(date): done ${g}"
