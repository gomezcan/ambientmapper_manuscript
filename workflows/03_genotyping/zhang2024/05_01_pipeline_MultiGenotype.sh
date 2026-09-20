#!/usr/bin/env bash
# 05_01_pipeline_MultiGenotype.sh — single-shot `ambientmapper run` on multiGenotypes_rep1
# (7 NAM genomes), the first pass over this library.
# Input : a configs TSV (sample, genome, bam, workdir), passed as $1. That TSV is not
#         preserved in configs/; configs/multiGenotypes_rep1.ambientmapper.json is the
#         JSON of record for the staged pipeline (04_01 to 04_05b) whose C0 run is the
#         one displayed (Fig. 4H to N).
# Run   : cd ${PROJECT_ROOT}/5_AmbientDetection && sbatch <repo>/workflows/03_genotyping/zhang2024/05_01_pipeline_MultiGenotype.sh <config_list.tsv>
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=35
#SBATCH --mem=80G
#SBATCH --job-name=ambientmapper_MultiGenotype
#SBATCH --partition=standard
#SBATCH --output=_logs/05_01_pipeline_MultiGenotype_%A.log

# conda activate <env from environment.yml>

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export GOTO_NUM_THREADS=1

# Usage: sbatch 05_01_pipeline_MultiGenotype.sh <config_list.tsv>
CONFIGS_TSV="$1"
if [[ -z "${CONFIGS_TSV}" ]]; then
  echo "Usage: sbatch $0 <config_list.tsv>"
  exit 1
fi

PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"   # root of the data tree (3_Mapping/, 5_AmbientDetection/, ...)
cd "${PROJECT_ROOT}/5_AmbientDetection"

echo "============================================"
echo "[$(date)] Starting AmbientMapper run: MultiGenotype"
echo "  Config: ${CONFIGS_TSV}"
echo "  7 maize genomes, independent wells"
echo "============================================"

# topk-genomes = 7 (all genomes); no winner-only to allow decontam across wells
ambientmapper run \
  --configs "${CONFIGS_TSV}" \
  --threads 35 \
  --min-barcode-freq 3 \
  --chunk-size-cells 50 \
  --assign-alpha 0.05 \
  --assign-k 10 \
  --assign-mapq-min 10 \
  --assign-xa-max 2 \
  --assign-chunksize 500000 \
  --assign-batch-size 6 \
  --no-genotyping-winner-only \
  --genotyping-pass1-workers 35 \
  --genotyping-beta 1 \
  --genotyping-min-reads 5 \
  --genotyping-chunk-rows 1000000 \
  --genotyping-bic-margin 6 \
  --genotyping-topk-genomes 7 \
  --genotyping-doublet-minor-min 0.20 \
  --genotyping-single-mass-min 0.6 \
  --genotyping-ratio-top1-top2-min 2.0 \
  --genotyping-empty-bic-margin 10 \
  --genotyping-empty-top1-max 0.6 \
  --genotyping-empty-ratio12-max 1.5 \
  --genotyping-eta-iters 2 \
  --genotyping-empty-seed-bic-min 10 \
  --genotyping-empty-tau-quantile 0.95

echo "[$(date)] Done: MultiGenotype steps 1-5"
