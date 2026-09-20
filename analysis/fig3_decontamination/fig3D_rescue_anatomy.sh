#!/bin/bash
#SBATCH --job-name=fig3D_rescue
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=_logs/fig3D_rescue_%j.out
#SBATCH --error=_logs/fig3D_rescue_%j.err
#
# fig3D_rescue_anatomy.sh -- SLURM wrapper for analysis/fig3_decontamination/fig3D_rescue_anatomy.R
# (Fig 3 panel D). Run from the repository root:
#   mkdir -p _logs && sbatch analysis/fig3_decontamination/fig3D_rescue_anatomy.sh
# FIG3D_EXPLORATORY=1 also builds the exploratory panels that are not in the manuscript.

set -euo pipefail
REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
cd "$REPO_ROOT"
mkdir -p _logs

# conda activate <env from environment.yml>

echo "=== fig3D_rescue_anatomy.R (SM2v2) ==="
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "PWD:  $(pwd)"
echo "R:    $(which Rscript)"
echo ""

Rscript analysis/fig3_decontamination/fig3D_rescue_anatomy.R

echo ""
echo "Finished: $(date)"
