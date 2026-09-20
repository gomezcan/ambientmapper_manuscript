#!/bin/bash
########## BATCH Lines for Resource Request ##########
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --job-name=Step5_9_consensus
#SBATCH --partition=standard
#SBATCH --output=_logs/Step5_9_consensus_%A_%a.log
#SBATCH --array=0-5
#
# 5_9_consensus_plate.sh -- Step 5.9 (plate arm): CONSENSUS meta-cells from the 100-seed sweep.
# ---------------------------------------------------------------------------------------------
#   0: SM2_At_TAIR10            Pre    3: SM2_B73_B73v5            Pre
#   1: Clean.SM2v2wd_At_TAIR10  wd     4: Clean.SM2v2wd_B73_B73v5  wd
#   2: Clean.SM2v2_At_TAIR10    nd     5: Clean.SM2v2_B73_B73v5    nd
#
# SLURM IS NOT REQUIRED. Measured: F builds in ~19 s for maize's 14,427 cells and the
#   hclust is ~4 s (dist 0.78 GB, dense coercion 1.6 GB => ~2.4 GB peak, one stage at a time). The
#   whole consensus is well under a minute per stage; only the optional split-half sweep costs
#   real time (~1 min At, ~12 min maize at 25 replicates). The header is here for parity with the
#   rest of step5_metacell -- running `bash 5_9_consensus_plate.sh <0-5>` locally is fine.
#
# PREREQ: the 100-seed sweeps from 5_3
#     SEEDS=100 MEMBERSHIP_ONLY=1 OUT_TAG=_S100 ...  ->  seacells_cps50_F{14235,44658}_N{14,286}_S100/
#
# k IS PINNED to the deliverable's n_SEACells (At 14, maize 286) and is IDENTICAL across the three
#   stages of a genome. Unmatched group counts manufacture stage effects (the count-matched controls).
#   Pinning k does NOT equalise meta-cell SIZE (At 75.8/73.7/49.7 cells per meta-cell), so condition
#   on `n_cells` in any cross-stage per-meta-cell comparison.
#
# THE ENGINE HARD-STOPS IF THE ACCEPTANCE GATE FAILS -- mean ARI(consensus, seed) must exceed
#   mean ARI(seed, seed). Set ALLOW_GATE_FAIL=1 only to inspect a known-bad object.
#
# NO THRESHOLD IS APPLIED. `loyalty` / `conf_meanF` are emitted continuous. If a figure draws a
#   line on them it must REPORT THE RETAINED FRACTION PER STAGE -- a confidence filter removes far
#   more of Pre than of wd and is therefore a stage-dependent treatment (`min.c=50` failure mode).
#
# Env knobs: TAG("") N_SEEDS(0=auto) N_SEEDPAIRS(500) KSWEEP_REPS(25) KSWEEP_MAXN(20000) LINKAGE(average)
# ---------------------------------------------------------------------------------------------

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-.}"
BASE="${PROJECT_ROOT}/6_socrates"
S5="${BASE}/SM2v2_plate/step5_metacell"
ENGINE="${SCRIPTS:-${BASE}/0_scripts}/common/4_3u_metacell_consensus.R"

# TAG separates runs that differ only in a knob (LINKAGE, N_SEEDS). Unset => the deliverable path.
OUT="${S5}/consensus${TAG:+_${TAG}}"

NSEEDS="${N_SEEDS:-0}"
NPAIRS="${N_SEEDPAIRS:-500}"
KREPS="${KSWEEP_REPS:-25}"
KMAXN="${KSWEEP_MAXN:-20000}"
LINK="${LINKAGE:-average}"

declare -a NAMES=( SM2_At_TAIR10 Clean.SM2v2wd_At_TAIR10 Clean.SM2v2_At_TAIR10 \
                   SM2_B73_B73v5 Clean.SM2v2wd_B73_B73v5 Clean.SM2v2_B73_B73v5 )
declare -a STAGES=( Pre wd nd Pre wd nd )
declare -a KPIN=( 14 14 14 286 286 286 )
SWEEP_At="${S5}/seacells_cps50_F14235_N14_S100"
SWEEP_maize="${S5}/seacells_cps50_F44658_N286_S100"

mkdir -p "$OUT"

INDICES=( "${@:-0 1 2 3 4 5}" )
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then INDICES=( "$SLURM_ARRAY_TASK_ID" ); fi
read -r -a INDICES <<< "${INDICES[*]}"

for i in "${INDICES[@]}"; do
  name="${NAMES[$i]:?bad index $i}"; stage="${STAGES[$i]}"; k="${KPIN[$i]}"
  if [[ "$i" -lt 3 ]]; then sweep="$SWEEP_At"; else sweep="$SWEEP_maize"; fi

  [[ -d "$sweep" ]] || { echo "ERROR: no sweep dir $sweep -- run 5_3 with SEEDS=100 OUT_TAG=_S100"; exit 2; }
  nseed_found=$(ls -1d "${sweep}"/seed* 2>/dev/null | wc -l)
  [[ "$nseed_found" -gt 1 ]] || { echo "ERROR: ${sweep} has ${nseed_found} seed dir(s)"; exit 2; }
  [[ -s "${sweep}/seed1/${name}.seacells.cell_to_seacell.tsv" ]] || {
    echo "ERROR: missing ${sweep}/seed1/${name}.seacells.cell_to_seacell.tsv"; exit 2; }

  echo "--- [${i}] ${name} (${stage}) | k=${k} | ${nseed_found} seeds | ${sweep##*/} -> ${OUT}"
  Rscript "$ENGINE" "$OUT" "$name" "$sweep" "$k" "$stage" \
          "$NSEEDS" "$NPAIRS" "$KREPS" "$KMAXN" "$LINK"
done

echo "$(date): done -> ${OUT}"
