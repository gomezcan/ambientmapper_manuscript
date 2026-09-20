#!/bin/bash
#SBATCH --job-name=figS5_synthetic
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:45:00
#SBATCH --output=_logs/%x_%j.log

# Supplementary Fig S5 (synthetic benchmark validation). Submit from the repo root:
#   mkdir -p _logs && sbatch analysis/supplementary/figS5.sh
# conda activate <env from environment.yml>

set -euo pipefail
REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
cd "$REPO_ROOT"

echo "=== figS5.R (synthetic validation) ==="
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "PWD:  $(pwd)"
echo "R:    $(which Rscript)"
echo ""

Rscript analysis/supplementary/figS5.R

echo ""
echo "Finished: $(date)"
