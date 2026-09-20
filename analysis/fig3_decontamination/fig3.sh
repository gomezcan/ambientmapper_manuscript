#!/bin/bash
#SBATCH --job-name=fig3
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=_logs/fig3_%j.out
#SBATCH --error=_logs/fig3_%j.err
#
# fig3.sh -- SLURM wrapper for analysis/fig3_decontamination/fig3.R (Fig 3 panels A, B, C, E).
# Run from the repository root:  mkdir -p _logs && sbatch analysis/fig3_decontamination/fig3.sh
# (or `bash analysis/fig3_decontamination/fig3.sh` in the foreground).

set -euo pipefail
REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
cd "$REPO_ROOT"
mkdir -p _logs

# conda activate <env from environment.yml>

echo "=== fig3.R (SM2v2) ==="
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "PWD:  $(pwd)"
echo "R:    $(which Rscript)"
echo ""

Rscript analysis/fig3_decontamination/fig3.R

echo ""
echo "Finished: $(date)"
