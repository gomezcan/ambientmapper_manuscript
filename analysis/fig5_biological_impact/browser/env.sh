#!/usr/bin/env bash
# =============================================================================
# browser/env.sh  -  ONE-TIME environment setup for the Fig 5 panel L browser chain (HPC only):
#   a dedicated conda env for pyGenomeTracks (github.com/deeptools/pyGenomeTracks).
# A SEPARATE env, not the mapping/QC env: pyGenomeTracks pins matplotlib/intervaltree versions
#   and upgrading inside the main env could disturb the bed -> bw pipeline that env serves.
# Run inside an interactive job, not on a login node: the conda-forge repodata solve needs
#   several GB and a login-node cgroup can kill it. Compute nodes need outbound internet.
#     salloc --partition=standard --time=01:00:00 --cpus-per-task=4 --mem=16G
#     bash analysis/fig5_biological_impact/browser/env.sh
# Fallback if the conda solve still dies:  PIP=1 bash analysis/fig5_biological_impact/browser/env.sh
#   (tiny python-only conda env + pip wheels for pyGenomeTracks/pyBigWig, no giant repodata solve)
# =============================================================================
set -euo pipefail

# source <conda install>/etc/profile.d/conda.sh   # initialise conda for this shell first

if conda env list | grep -q "^pygenometracks "; then
  echo "env 'pygenometracks' already exists - nothing to do"
elif [ -n "${PIP:-}" ]; then
  echo "creating minimal env + pip install (low-memory path) ..."
  conda create -y -n pygenometracks -c conda-forge python=3.11 pip
  conda activate pygenometracks
  pip install pyGenomeTracks
  conda deactivate
else
  # mamba if available (much faster solver), conda otherwise
  SOLVER=conda; command -v mamba >/dev/null 2>&1 && SOLVER=mamba
  echo "creating env with $SOLVER ..."
  "$SOLVER" create -y -n pygenometracks -c conda-forge -c bioconda pygenometracks
fi

conda activate pygenometracks
echo
echo "=== verification ==="
pyGenomeTracks --version
python -c "import pyBigWig; print('pyBigWig', pyBigWig.__version__)"
echo
echo "OK. Next (after the bigwigs exist, from the repo root):"
echo "  mkdir -p figures/main/fig5/browser/_logs_plots && sbatch analysis/fig5_biological_impact/browser/plot.sh"
