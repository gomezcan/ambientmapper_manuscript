#!/bin/bash
########## BATCH Lines for Resource Request ##########
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --job-name=Step5_0_consensus_peaks
#SBATCH --partition=standard
#SBATCH --output=_logs/Step5_0_consensus_peaks_%A_%a.log
#SBATCH --array=0-1
#
# 5_0_make_consensus_peaks.sh -- Step 5.0 (plate arm): FROZEN consensus ACR set.
# ---------------------------------------------------------------------------
# Emits the feature space that Step 5 (meta-cells / SEACells) and Step 6 (DAR) BOTH consume,
# from the existing whole-object bulk MACS2 calls in step0_qc/*_macs2_temp/.
#
#   0: TAIR10 (At)      1: B73v5 (B73)
#
# THE RULE: one peak set, frozen across Pre/wd/nd. Stage-specific sets confound a "lost" peak with
# a peak-calling threshold shift (same logic as the frozen cluster config, plan decision B).
# PRIMARY = Pre (defines the measurement space; wd/nd are measured in it). The union is emitted for
# the reverse-direction supplementary arm.
# NOT pseudobulk-per-cluster: Step 5 feeds the boundary test, so cluster-derived features would make
# that test circular by construction.
#
# Measured blind spot: wd-unique vs Pre 0.04-0.64%; nd-unique 4.19-8.45% but median
# width 150 bp = the MACS floor (--extsize 150), i.e. the weakest possible calls.
#
# Local run (no SLURM needed, ~1 min):
#   bash 0_scripts/part2_split/step5_metacell/5_0_make_consensus_peaks.sh 0
# ---------------------------------------------------------------------------

# conda activate ambientmapper-manuscript   (environment.yml at the repo root)
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-.}"
BASE="${PROJECT_ROOT}/6_socrates"
SCRIPTS="${SCRIPTS:-${BASE}/0_scripts}"
QC="${BASE}/SM2v2_plate/step0_qc"
OUT="${BASE}/SM2v2_plate/_consensus_peaks"
s5_0="${SCRIPTS}/common/5_0_make_consensus_peaks.R"

declare -a GENOMES=( TAIR10 B73v5 )
declare -a SPECIES=( At     B73   )

i="${SLURM_ARRAY_TASK_ID:-${1:?set via --array or pass 0/1}}"
g="${GENOMES[$i]:?bad array index $i}"; sp="${SPECIES[$i]}"

np () {  # $1 = object prefix  ->  its genome-wide combined narrowPeak
  echo "${QC}/$1_macs2_temp/$1_peaks_combined_peaks.narrowPeak"
}
PRE=$(np "SM2_${sp}_${g}")
WD=$(np  "Clean.SM2v2wd_${sp}_${g}")
ND=$(np  "Clean.SM2v2_${sp}_${g}")

for f in "$PRE" "$WD" "$ND"; do
  [[ -s "$f" ]] || { echo "ERROR: missing or empty narrowPeak: $f"; exit 2; }
done

mkdir -p "$OUT"
echo "$(date): Step 5.0 consensus peaks -> ${g}"
Rscript "$s5_0" "$OUT" "$g" "$PRE" "$WD" "$ND"
echo "$(date): done ${g}"
