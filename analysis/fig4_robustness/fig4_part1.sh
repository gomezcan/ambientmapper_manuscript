#!/bin/bash
#SBATCH --job-name=fig4_part1
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=_logs/%x_%j.log

# Fig 4 panels A to C (synthetic benchmark). Submit from the repo root:
#   mkdir -p _logs && sbatch analysis/fig4_robustness/fig4_part1.sh
# conda activate <env from environment.yml>

set -euo pipefail
REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
cd "$REPO_ROOT"

echo "=== fig4_part1.R ==="
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "PWD:  $(pwd)"
echo "R:    $(which Rscript)"
echo ""

Rscript analysis/fig4_robustness/fig4_part1.R

echo ""
echo "Finished: $(date)"
