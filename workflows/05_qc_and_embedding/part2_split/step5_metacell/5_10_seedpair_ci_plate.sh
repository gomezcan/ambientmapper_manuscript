#!/bin/bash
#SBATCH --job-name=Step5_10_seedpair_ci
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --partition=standard
#SBATCH --output=_logs/Step5_10_seedpair_ci_%j.log
## =================================================================================================
## 5_10_seedpair_ci_plate.sh -- driver for 4_3y (AMI + bootstrap CIs on the 100-seed sweep).
## SLURM OPTIONAL: ~1 min for At, ~5-15 min for maize, <2 GB -- runs fine on either computer.
##   bash 5_10_seedpair_ci_plate.sh At          # or: maize
##   bash 5_10_seedpair_ci_plate.sh maize 10000 # optional bootstrap B (default 10000)
## Reads  : step5_metacell/consensus/<prefix>.consensus.{seedpairs,summary}.tsv
##          step5_metacell/seacells_cps50_F*_N*_S100/seed*/<prefix>.seacells.cell_to_seacell.tsv
## Writes : step5_metacell/consensus/<genome>.{seedpair_ami.tsv,ari_ami_ci.tsv,ari_ami_ci.md}
## No pre-existing file is modified.
## =================================================================================================
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-.}"
AMB_BASE="${AMB_BASE:-${PROJECT_ROOT}/6_socrates}"
SCRIPTS="${SCRIPTS:-${AMB_BASE}/0_scripts}"
[[ -d "${AMB_BASE}" ]] || { echo "ERROR: 6_socrates not found at ${AMB_BASE} (set PROJECT_ROOT)" >&2; exit 1; }
cd "${AMB_BASE}"

G="${1:?usage: bash 5_10_seedpair_ci_plate.sh At|maize [B]}"
B="${2:-10000}"
M=SM2v2_plate/step5_metacell
case "$G" in
  At)    SWEEP="$M/seacells_cps50_F14235_N14_S100" ;;
  maize) SWEEP="$M/seacells_cps50_F44658_N286_S100" ;;
  *) echo "ERROR: genome must be At or maize, got '$G'" >&2; exit 1 ;;
esac
[[ -d "$SWEEP" ]] || { echo "ERROR: sweep dir missing: $SWEEP" >&2; exit 1; }

Rscript "${SCRIPTS}/common/4_3y_seedpair_ami_ci.R" "$M/consensus" "$SWEEP" "$M/consensus" "$G" "$B"
