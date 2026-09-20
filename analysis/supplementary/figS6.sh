#!/bin/bash
#SBATCH --job-name=figS6_root1
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=_logs/%x_%j.log

# Supplementary Fig S6 (Root1 sub1k_B validation). Submit from the repo root:
#   mkdir -p _logs && sbatch analysis/supplementary/figS6.sh
# conda activate <env from environment.yml>

set -euo pipefail
REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
cd "$REPO_ROOT"

echo "=== figS6.R (Root1 validation) ==="
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "PWD:  $(pwd)"
echo "R:    $(which Rscript)"
echo ""

Rscript analysis/supplementary/figS6.R

echo ""
echo "Finished: $(date)"
