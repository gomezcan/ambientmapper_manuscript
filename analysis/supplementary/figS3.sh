#!/usr/bin/env bash
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=30G
#SBATCH --job-name=figS3_gridscan
#SBATCH --partition=standard
#SBATCH --output=_logs/figS3_gridscan_%j.log
#
# Fig S3: UMAP parameter grid scan (figS3.R, panels B and C plus the metrics table) followed by the
# shipped panel A (figS3A_replot.R). Needs the Socrates R environment; about 30 GB of memory.
# Usage (from the repository root): sbatch analysis/supplementary/figS3.sh

set -euo pipefail
# conda activate <env from environment.yml>

REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
cd "${REPO_ROOT}"
mkdir -p _logs

DATA="data/processed/scifiATAC_B73_Arabidopsis/socrates"
soc_rds="${DATA}/SM2/step1_integrate/SM2.full.SocObj.rds"
meta_tsv="${DATA}/SM2/step2_metaqc/SM2.FRiP0.2.FULL.minDepth200.stagePreClean.updated_metadata_v4.txt"
outdir="figures/supplementary/figS3"; mkdir -p "${outdir}/plots"

Rscript analysis/supplementary/figS3.R "${soc_rds}" "${meta_tsv}" "${outdir}"
Rscript analysis/supplementary/figS3A_replot.R "${outdir}/UMAP_grid_scan.metrics.tsv" "${outdir}/FigS3A_knn_preservation_grid_scan"

echo ' .. done ..'
