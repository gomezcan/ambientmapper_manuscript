#!/bin/bash
########## BATCH Lines for Resource Request ##########
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --job-name=Step5_3_seacells_plate
#SBATCH --partition=standard
#SBATCH --output=_logs/Step5_3_seacells_plate_%A_%a.log
#SBATCH --array=0-5
#
# 5_3_run_seacells_plate.sh -- Step 5.3 (plate arm): meta-cells via SEACells.
# ---------------------------------------------------------------------------
#   0: SM2_At_TAIR10            Pre    3: SM2_B73_B73v5            Pre
#   1: Clean.SM2v2wd_At_TAIR10  wd     4: Clean.SM2v2wd_B73_B73v5  wd
#   2: Clean.SM2v2_At_TAIR10    nd     5: Clean.SM2v2_B73_B73v5    nd
#
# PREREQ (once, INTERACTIVELY -- compute nodes usually have no network):
#     bash 0_scripts/setup/0_10_install_seacells_env.sh
# PREREQ (data): Step 5.1 bundles in SM2v2_plate/step5_metacell/<name>.seacells.*
#
# THE ONE THING SEACELLS CANNOT DO FOR US: it takes ONE AnnData and cannot know Pre/wd/nd are meant
# to be compared. Measured: At nd cells carry +37% median counts vs Pre, so at fixed
# cells-per-meta-cell nd's meta-cells come out deeper FOR FREE -- a stage-dependent bandwidth
# pointing the same way as the cleaning effect under study. RAREFY_TO fixes it, and it MUST be the
# SAME value across the three stages of a genome, set from the SHALLOWEST stage.
#
#   Measured F_common (10th pct of pooled depth, all stages, greedy-kNN proxy -- 5_2s_sweep_k_F):
#      At    k=25 -> 3738   k=50 -> 8058   k=75 -> 12977   k=100 -> 16986
#      maize k=25 -> 5485   k=50 -> 12872  k=75 -> 19449   k=100 -> 27256
#   SEACells' assignment differs from that proxy, so RE-READ the "pooled fragments/meta-cell" line
#   the runner prints for all three stages, then re-run with RAREFY_TO set from the shallowest.
#   Two passes is the intended workflow -- do NOT guess it in one.
#
# CELLS_PER_SEACELL: SEACells' published default is ~75. Our sweep found the choice barely matters
# for maize (adjusted purity 0.876 -> 0.797 across k=10..100). For At it is a documented limit case:
# 1061 cells / 75 = 14 meta-cells, and <=20 meta-cells gave a between-seed spread that swamped any
# stage effect. Run At anyway -- a measured floor beats an omission -- but read it as a limit, and
# use SEEDS>1 so the instability is visible rather than implied.
#
# Env knobs: CELLS_PER_SEACELL(75) RAREFY_TO(0=off) SEEDS(1) NWAYPOINT(10) NNEIGHBORS(15)
#            SEEDS_PER_TASK(=SEEDS) OUT_TAG("") MEMBERSHIP_ONLY(0) SKIP_EXISTING(1)
#
# 100-SEED REPLICATE SWEEP -- the run that puts an error bar on co-assignment/ARI.
#   Array layout, as in 5_7:   sample = task % 6 ,  chunk = task / 6
#   Write to a NEW dir via OUT_TAG. Re-running seeds 1-5 inside seacells_cps50_F14235_N14/
#     would overwrite the At deliverable that svdnull_pooled_n1000/ was scored against.
#   Use MEMBERSHIP_ONLY=1. The readout needs only cell_to_seacell.tsv; the .metacells.h5ad is
#     134 MB per maize stage per seed, i.e. ~40 GB over 100 seeds x 3 stages (~300 MB with it on).
#
#   At    (~25 s/fit  => ~45 min/task, all 100 seeds in one task):
#     SEEDS=100 SEEDS_PER_TASK=100 MEMBERSHIP_ONLY=1 OUT_TAG=_S100 \
#       CELLS_PER_SEACELL=50 RAREFY_TO=14235 N_SEACELLS=14 sbatch --array=0-2 5_3_run_seacells_plate.sh
#
#   maize (~11 min/fit => ~2 h/task at 10 seeds; 10 chunks x 3 stages = 30 tasks):
#     SEEDS=100 SEEDS_PER_TASK=10 MEMBERSHIP_ONLY=1 OUT_TAG=_S100 \
#       CELLS_PER_SEACELL=50 RAREFY_TO=44658 N_SEACELLS=286 sbatch --array=3-59:6 ...   # Pre
#                                                            sbatch --array=4-59:6 ...   # wd
#                                                            sbatch --array=5-59:6 ...   # nd
#   SKIP_EXISTING=1 makes this resume-safe -- just resubmit the same array if a task is killed.
# ---------------------------------------------------------------------------

set -euo pipefail

# conda activate seacells   (created once by setup/0_10_install_seacells_env.sh)

PROJECT_ROOT="${PROJECT_ROOT:-.}"
BASE="${PROJECT_ROOT}/6_socrates"
S5="${BASE}/SM2v2_plate/step5_metacell"
RUNNER="${SCRIPTS:-${BASE}/0_scripts}/common/5_3_run_seacells.py"

CPS="${CELLS_PER_SEACELL:-75}"
RARE="${RAREFY_TO:-0}"
SEEDS="${SEEDS:-1}"
NWP="${NWAYPOINT:-10}"
NNB="${NNEIGHBORS:-15}"
# N_SEACELLS: explicit meta-cell COUNT, overriding the cells-per-meta-cell heuristic. 0 = off.
# Why this exists: with CPS alone the count is derived from n_cells, so a stage that lost cells to
# cleaning gets FEWER meta-cells for free -- pass 1 gave At 14 (Pre) / 14 (wd) / 9 (nd). The
# The count-matched controls showed unmatched counts MANUFACTURE the stage effect (that is what killed the
# "maize nd lowest 17/17" claim). RAREFY_TO matches DEPTH; this matches COUNT. For a stage comparison
# you generally need BOTH, and both must take the SAME value across the stages of one genome.
NSC="${N_SEACELLS:-0}"

# --- seed chunking + output isolation (for the 100-seed replicate sweep) -----
# SEEDS_PER_TASK lets ONE array task cover a RANGE of seeds. Without it, SEEDS=100 runs 100 fits
# serially in a single task: fine for At (~25 s/fit => ~45 min) but ~18 h for maize (~11 min/fit),
# which walls out against --time. Array layout follows 5_7's convention:
#       sample = task % 6 ,  chunk = task / 6
# Unset, SEEDS_PER_TASK defaults to SEEDS, so chunk 0 is the only chunk and any --array=0-5
# submission behaves exactly as before.
SPT="${SEEDS_PER_TASK:-$SEEDS}"
# OUT_TAG appends to the output directory name. Use it to keep a large replicate sweep OUT of the
# frozen deliverable dirs (seacells_cps50_F14235_N14/, seacells_cps50_F44658_N286/) that the
# pooled svdnull scoring in svdnull_pooled_n1000/ was run against. Re-running seeds 1-5 in place
# would overwrite that deliverable.
OUT_TAG="${OUT_TAG:-}"
# MEMBERSHIP_ONLY skips the aggregated matrix and .metacells.h5ad. The co-assignment / ARI readout
# needs ONLY cell_to_seacell.tsv, and the h5ad is 134 MB per maize stage per seed -- 100 seeds x 3
# stages would be ~40 GB. Leave at 0 if you need aggregated profiles (the rZ / cell-type track).
MEMONLY="${MEMBERSHIP_ONLY:-0}"
# SKIP_EXISTING makes the sweep resume-safe: a seed whose cell_to_seacell.tsv is already present
# and non-empty is skipped, so a task killed on walltime can simply be resubmitted.
SKIP="${SKIP_EXISTING:-1}"

declare -a NAMES=( SM2_At_TAIR10 Clean.SM2v2wd_At_TAIR10 Clean.SM2v2_At_TAIR10 \
                   SM2_B73_B73v5 Clean.SM2v2wd_B73_B73v5 Clean.SM2v2_B73_B73v5 )
declare -a STAGES=( Pre wd nd Pre wd nd )

task="${SLURM_ARRAY_TASK_ID:-${1:?set via --array or pass 0-5}}"
i=$(( task % 6 )); chunk=$(( task / 6 ))
name="${NAMES[$i]:?bad array index $i}"; stage="${STAGES[$i]}"
bundle="${S5}/${name}.seacells"

for ext in mtx.gz svd.tsv obs.tsv var.tsv; do
  [[ -s "${bundle}.${ext}" ]] || { echo "ERROR: missing ${bundle}.${ext} -- run Step 5.1 first"; exit 2; }
done

# Build the suffix with an `if`, NOT `VAR="...$( cond && echo x )"`: with RARE=0 the `[[ ]]` is
# false, `&&` short-circuits, the command substitution exits 1, and a simple assignment inherits
# that status -- so `set -e` aborts the script before the first echo (0-byte logs). It only bites
# RAREFY_TO=0, i.e. pass 1. `cond && VAR=x` on its own line has the same trap; the `if` form is safe.
SUFFIX=""
if [[ "$RARE" != 0 ]]; then SUFFIX="${SUFFIX}_F${RARE}"; fi
if [[ "$NSC"  != 0 ]]; then SUFFIX="${SUFFIX}_N${NSC}"; fi
OUT="${S5}/seacells_cps${CPS}${SUFFIX}${OUT_TAG}"
mkdir -p "$OUT"

# Seed range for THIS task. Same `if` discipline as the SUFFIX block above -- `(( ))` returns 1
# when its expression evaluates to 0, so `(( a > b )) && x=y` as a standalone statement is the
# same set -e trap.
s0=$(( chunk * SPT + 1 ))
s1=$(( s0 + SPT - 1 ))
if (( s1 > SEEDS )); then s1="$SEEDS"; fi
if (( s0 > SEEDS )); then
  echo "$(date): task ${task} -> seeds ${s0}..${s1} lies beyond SEEDS=${SEEDS}; nothing to do"
  exit 0
fi

MEMFLAG=()
if [[ "$MEMONLY" != 0 ]]; then MEMFLAG=(--membership-only); fi

echo "$(date): Step 5.3 | ${name} (${stage}) | cells_per_seacell=${CPS} rarefy_to=${RARE}" \
     "n_seacells=${NSC} | task=${task} chunk=${chunk} seeds ${s0}..${s1} of ${SEEDS}" \
     "| membership_only=${MEMONLY} | out=${OUT}"
for sd in $(seq "$s0" "$s1"); do
  odir="$OUT"
  if [[ "$SEEDS" -gt 1 ]]; then odir="${OUT}/seed${sd}"; mkdir -p "$odir"; fi
  marker="${odir}/${name}.seacells.cell_to_seacell.tsv"
  if [[ "$SKIP" != 0 && -s "$marker" ]]; then
    echo "--- seed ${sd}: SKIP, already present -- ${marker}"
    continue
  fi
  echo "--- seed ${sd} -> ${odir}"
  python "$RUNNER" --bundle "$bundle" --out "$odir" \
         --cells-per-seacell "$CPS" --n-seacells "$NSC" --rarefy-to "$RARE" --seed "$sd" \
         --n-waypoint-eigs "$NWP" --n-neighbors "$NNB" ${MEMFLAG[@]+"${MEMFLAG[@]}"}
done
echo "$(date): done ${name} (task ${task}, seeds ${s0}..${s1})"
