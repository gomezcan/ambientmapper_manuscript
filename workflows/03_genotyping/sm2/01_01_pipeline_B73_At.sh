#!/usr/bin/env bash
# 01_01_pipeline_B73_At.sh — SM2 (combined-reference arm) AmbientMapper run (single `ambientmapper run`:
# extract -> filter -> chunks -> assign -> genotyping) on the maize B73 + Arabidopsis
# library, each genome mapped independently.
#
# Input : configs/SM2_AtB73.list.tsv (sample, genome, bam, workdir), passed as $1
# Output: SM2/final/SM2_cells_calls.tsv.gz and SM2/cell_map_ref_chunks/*_filtered.tsv.gz,
#         consumed by workflows/04_decontamination/sm2/01_02_decontam_B73_At.sh, whose
#         tables are the inputs of Fig. 1B to E.
# Run   : cd ${PROJECT_ROOT}/5_AmbientDetection && sbatch <repo>/workflows/03_genotyping/sm2/01_01_pipeline_B73_At.sh configs/SM2_AtB73.list.tsv
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=35
#SBATCH --mem=40G
#SBATCH --job-name=ambientmapper_SM2
#SBATCH --partition=standard
#SBATCH --output=_logs/01_01_pipeline_B73_At_%A.log

# conda activate <env from environment.yml>

PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"   # root of the data tree (3_Mapping/, 5_AmbientDetection/, ...)
cd "${PROJECT_ROOT}/5_AmbientDetection"

# Prevent thread oversubscription inside numpy/pandas/etc.
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export GOTO_NUM_THREADS=1

CONFIGS_TSV="$1"

ambientmapper run \
  --configs "${CONFIGS_TSV}" \
  --threads 36 \
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
  --genotyping-topk-genomes 3 \
  --genotyping-doublet-minor-min 0.30 \
  --genotyping-single-mass-min 0.7 \
  --genotyping-ratio-top1-top2-min 2.0 \
  --genotyping-empty-bic-margin 10 \
  --genotyping-empty-top1-max 0.6 \
  --genotyping-empty-ratio12-max 1.5 \
  --genotyping-eta-iters 2 \
  --genotyping-empty-seed-bic-min 10 \
  --genotyping-empty-tau-quantile 0.95
