#!/usr/bin/env bash
# 02_01_pipeline_Root1.sh — Root1_rep1 (26 NAM genomes): extract, filter, chunks
# and the first assign pass.
#
# Input : configs/Root1_rep1.ambientmapper.json (26 per-genome BAMs from 3_Mapping/Root1_rep1/)
# Output: Root1_rep1/{qc,filtered_QCFiles,cell_map_ref_chunks}/
# Note  : this assign pass used the legacy bundle (alpha 1e-6, k 5). The assign
#         that feeds Fig. 4D to G was re-run at alpha 0.05 by 02_12a, then genotyped
#         by 02_12. Only extract, filter and chunks from this script are reused.
# Run   : cd ${PROJECT_ROOT}/5_AmbientDetection && sbatch <repo>/workflows/03_genotyping/root1/02_01_pipeline_Root1.sh

########## BATCH Lines for Resource Request ##########
#SBATCH --time=20:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=160G
#SBATCH --job-name=ambientmapper_root
#SBATCH --partition=standard
#SBATCH --output=_logs/02_01_pipeline_Root1_%A.log


###################################
#######   Conda / Modules   #######
###################################
# conda activate <env from environment.yml>

# Prevent thread oversubscription inside numpy/pandas/etc.
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export GOTO_NUM_THREADS=1

PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"   # root of the data tree (3_Mapping/, 5_AmbientDetection/, ...)
cd "${PROJECT_ROOT}/5_AmbientDetection"

CONFIGS_json="configs/Root1_rep1.ambientmapper.json"

if [[ -z "${CONFIGS_json}" ]]; then
  echo "Usage: sbatch $0 <config.json>"
  exit 1
fi

# Resolve to absolute path (needed for the dependent job submission)
CONFIGS_json="$(realpath "${CONFIGS_json}")"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

echo "============================================"
echo "[$(date)] Starting ambientmapper pipeline"
echo "  Config: ${CONFIGS_json}"
echo "============================================"

#################################################
#######   Steps 1-4: extract → assign      #####
#################################################

# 1) extract
echo "[$(date)] Step 1: extract ..."
ambientmapper extract --config "${CONFIGS_json}" --threads 30
echo "[$(date)] Done: step 1 (extract)"

# 2) filter
echo "[$(date)] Step 2: filter ..."
ambientmapper filter -c "${CONFIGS_json}" -t 2 --min-barcode-freq 5
echo "[$(date)] Done: step 2 (filter)"

# 3) chunks
echo "[$(date)] Step 3: chunks ..."
ambientmapper chunks --config "${CONFIGS_json}" --chunk-size-cells 100
echo "[$(date)] Done: step 3 (chunks)"

# 4) assign
echo "[$(date)] Step 4: assign ..."
ambientmapper assign --config "${CONFIGS_json}" \
    --threads 60 \
    --alpha 0.000001 \
    --k 5 \
    --mapq-min 10 \
    --xa-max 0 \
    --chunksize 1000000 \
    --edges-subsample 100000 \
    --ecdf-subsample  100000 \
    --batch-size 4 \
    --score-workers 4 \
    --score-duckdb \
    --duckdb-threads 3 \
    --ecdf-duckdb-threads 4 \
    --skip-ecdf \
    --skip-edges

echo "[$(date)] Done: step 4 (assign)"

echo "============================================"
echo "[$(date)] Pipeline steps 1-4 complete."
echo "  Next: 02_12a_pipeline_full_root1_alpha005.sh (assign at alpha 0.05), then 02_12_genotyping_full_root1_phase4.sh."
echo "============================================"
