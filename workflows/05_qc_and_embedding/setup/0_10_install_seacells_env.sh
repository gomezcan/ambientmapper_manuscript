#!/bin/bash
# =============================================================================
# 0_10_install_seacells_env.sh -- one-shot: create the `seacells` conda env on the HPC.
#
# RUN ONCE, INTERACTIVELY (not via sbatch -- compute nodes usually have no network):
#     bash 0_scripts/setup/0_10_install_seacells_env.sh
#
# If the login node forbids heavy builds, grab a short interactive job first:
#     salloc --partition=standard --cpus-per-task=4 --mem=16G --time=1:00:00
#
# WHY A DEDICATED ENV:
#   SEACells pins scanpy / anndata / numpy / jax. Installing those into the R + Socrates environment
#   that every other script in this stage uses would put that stack at risk. A separate env is
#   isolated and reversible (`conda env remove -n seacells`).
#
# On macOS arm64 this recipe resolved to:
#   SEACells 0.3.3 | scanpy 1.11.5 | anndata 0.11.4 | numpy 2.2.6 | scipy 1.15.3 | jax 0.6.2
# Python is pinned to 3.10 deliberately: SEACells is 2023-era and is not reliable on 3.12+.
# =============================================================================
set -euo pipefail

CONDA_SH="${CONDA_SH:-$(conda info --base 2>/dev/null)/etc/profile.d/conda.sh}"
ENV_NAME="${ENV_NAME:-seacells}"

[[ -f "$CONDA_SH" ]] || { echo "ERROR: conda.sh not found at $CONDA_SH  (override with CONDA_SH=...)"; exit 2; }
source "$CONDA_SH"

if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  echo "env '$ENV_NAME' already exists -- skipping create (remove with: conda env remove -n $ENV_NAME)"
else
  echo "== creating env '$ENV_NAME' (python 3.10)"
  conda create -n "$ENV_NAME" python=3.10 -y
fi

conda activate "$ENV_NAME"
echo "== python: $(python -V)  at $(which python)"

echo "== installing cmake (SEACells build dependency)"
python -m pip install --quiet cmake

echo "== installing SEACells"
python -m pip install --quiet SEACells

# REQUIRED, not optional. SEACells' build_graph imports `tqdm.notebook`, which raises
#   ImportError: IProgress not found. Please update jupyter and ipywidgets
# on the FIRST kernel-construction call in any headless/script context -- i.e. every SLURM job.
echo "== installing ipywidgets (SEACells imports tqdm.notebook -- fails headless without it)"
python -m pip install --quiet ipywidgets

echo "== verifying"
python - <<'PY'
import SEACells, scanpy as sc, anndata as ad, numpy as np, scipy, ipywidgets
from SEACells.core import SEACells as SC
import SEACells.evaluate as ev
print("  SEACells 0.3.3 ok | scanpy", sc.__version__, "| anndata", ad.__version__,
      "| numpy", np.__version__, "| scipy", scipy.__version__,
      "| ipywidgets", ipywidgets.__version__)
assert all(hasattr(SC, m) for m in
           ["construct_kernel_matrix", "initialize_archetypes", "fit",
            "get_hard_assignments", "get_soft_assignments"])
assert all(hasattr(ev, m) for m in ["compactness", "separation", "compute_celltype_purity"])

# End-to-end headless check on a tiny synthetic AnnData. This is the call that fails without
# ipywidgets, so a plain `import SEACells` is NOT sufficient verification.
import scipy.sparse as sp
rng = np.random.default_rng(0)
A = ad.AnnData(X=sp.csr_matrix((rng.random((300, 60)) < 0.1).astype(float)))
A.obsm["X_svd"] = rng.standard_normal((300, 6))
m = SC(A, build_kernel_on="X_svd", n_SEACells=6, verbose=False)
m.construct_kernel_matrix(); m.initialize_archetypes()
# SEACells RAISES (does not warn) when it hits max_iter without converging. Random synthetic data
# will not converge in a few iterations, and convergence is irrelevant to "does the install work",
# so tolerate it here -- what we are verifying is that the code path runs headless.
try:
    m.fit(max_iter=5, min_iter=2)
except RuntimeWarning:
    pass
assert m.get_hard_assignments().shape[0] == 300
print("  headless end-to-end check passed")
PY

echo
echo "DONE. Use it from SLURM with:"
echo "    source $CONDA_SH && conda activate $ENV_NAME"
