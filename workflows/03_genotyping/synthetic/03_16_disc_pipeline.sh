#!/bin/bash
#SBATCH --job-name=disc_pipeline
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=_logs/03_16_disc_pipeline_%j.log

# =============================================================================
# 03_16_disc_pipeline.sh — Discriminative Peak Benchmark (Phases 2-6)
#
# Runs the full pipeline (ART → barcoding → BWA → extract → chunks+assign)
# on the discriminative peak set (843 peaks with ≥1 SNP per 75bp read).
#
# Single SLURM job because the smaller peak set makes each phase fast.
# Prerequisite: run 03_03_filter_disc_peaks.py first to create
#               synthetic_disc/orthologs/
#
# Output: synthetic_disc/{dataset}/cell_map_ref_chunks/*_filtered.tsv.gz
#         for all 15 datasets (same design as original benchmark)
# =============================================================================

set -euo pipefail

REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"   # this repository (sourced helpers, sibling scripts)
PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"   # root of the data tree (3_Mapping/, 5_AmbientDetection/, ...)

cd "${PROJECT_ROOT}/5_AmbientDetection"


export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export GOTO_NUM_THREADS=1

# --- Configuration ---
BASEDIR="synthetic_disc"
ORTHO_DIR="${BASEDIR}/orthologs"
READS_DIR="${BASEDIR}/reads"
BARCODED_DIR="${BASEDIR}/barcoded"
MAPPING_DIR="${BASEDIR}/mapping"
GENOMES=(B73 Il14H Ki11)
GENOME_KEYS=(B73v5 Il14H Ki11)  # BWA index naming
# BWA indexes of the NAM assemblies (Index_Zm_<genome>_bwa). BWA_INDEX_ROOT holds the Zea/ tree.
BWA_INDEX_ROOT="${BWA_INDEX_ROOT:?set BWA_INDEX_ROOT to the directory holding Zea/NAN_Indexes/}"
BWA_INDEX_BASE="${BWA_INDEX_ROOT}/Zea/NAN_Indexes"
THREADS=16

DATASETS=(
    alpha_000
    alpha_002_Il14H alpha_005_Il14H alpha_010_Il14H alpha_020_Il14H
    alpha_030_Il14H alpha_040_Il14H alpha_050_Il14H
    alpha_002_Ki11  alpha_005_Ki11  alpha_010_Ki11  alpha_020_Ki11
    alpha_030_Ki11  alpha_040_Ki11  alpha_050_Ki11
)

echo "============================================"
echo "[$(date)] Discriminative Peak Benchmark Pipeline"
echo "  Base dir:  ${BASEDIR}"
echo "  Orthologs: ${ORTHO_DIR}"
echo "============================================"

# --- Verify Phase 1b output ---
if [ ! -f "${ORTHO_DIR}/ortholog_map.tsv" ]; then
    echo "ERROR: ${ORTHO_DIR}/ortholog_map.tsv not found."
    echo "  Run: python 03_03_filter_disc_peaks.py first."
    exit 1
fi
N_PEAKS=$(tail -n +2 "${ORTHO_DIR}/ortholog_map.tsv" | wc -l)
echo "  Peaks: ${N_PEAKS}"

# =====================================================================
# PHASE 2: ART Read Simulation (250x to compensate for fewer peaks)
# =====================================================================
echo ""
echo "========== PHASE 2: ART Read Simulation =========="

# conda activate <env providing art_illumina (ART)>

FOLD_COV=450          # higher than original (100x) to compensate for fewer peaks
READ_LEN=75
INSERT_MEAN=166
INSERT_SD=30
SEED=42

mkdir -p "${READS_DIR}"

for G in "${GENOMES[@]}"; do
    FA="${ORTHO_DIR}/${G}_peak_seqs.fa"
    PREFIX="${READS_DIR}/${G}_"

    echo "[$(date)] ART: ${G} (${FOLD_COV}x)..."
    art_illumina \
        -ss HS25 \
        -i "${FA}" \
        -p \
        -l "${READ_LEN}" \
        -f "${FOLD_COV}" \
        -m "${INSERT_MEAN}" \
        -s "${INSERT_SD}" \
        -o "${PREFIX}" \
        --rndSeed "${SEED}" \
        --noALN \
        -q 2>&1

    N_READS=$(grep -c "^@" "${PREFIX}1.fq" 2>/dev/null || echo 0)
    echo "  ${G}: ${N_READS} read pairs"
    gzip -f "${PREFIX}1.fq" &
    gzip -f "${PREFIX}2.fq" &
    wait
done

# conda activate <env from environment.yml>

# =====================================================================
# PHASE 3: Barcode Assignment (reuses existing script)
# =====================================================================
echo ""
echo "========== PHASE 3: Barcode Assignment =========="

python "${REPO_ROOT}/workflows/03_genotyping/synthetic/03_05_assign_barcodes.py" \
    --reads-dir "${READS_DIR}" \
    --outdir "${BARCODED_DIR}" \
    --seed 42

echo "[$(date)] Phase 3 complete"

module load Bioinformatics bwa/0.7.17-mil4ns7

# =====================================================================
# PHASE 4: BWA Mapping (45 jobs in parallel batches)
# =====================================================================
echo ""
echo "========== PHASE 4: BWA Mapping =========="

PARALLEL=4  # concurrent BWA jobs

bwa_map_one() {
    local DS="$1"
    local GIDX="$2"
    local G="${GENOMES[$GIDX]}"
    local GK="${GENOME_KEYS[$GIDX]}"

    local R1="${BARCODED_DIR}/${DS}/all_reads_R1.fq.gz"
    local R2="${BARCODED_DIR}/${DS}/all_reads_R2.fq.gz"
    local OUTBAM="${MAPPING_DIR}/${DS}/syn_${GK}.sorted.bam"
    local BWA_IDX="${BWA_INDEX_BASE}/Index_Zm_${GK}_bwa"

    mkdir -p "${MAPPING_DIR}/${DS}"

    bwa mem -M -t 4 "${BWA_IDX}" "${R1}" "${R2}" 2>/dev/null \
        | samtools sort -@ 2 -o "${OUTBAM}" -
    samtools index "${OUTBAM}"
    if [ ! -s "${OUTBAM}" ]; then
        echo "ERROR: BWA failed for ${DS} x ${GK}" >&2
        return 1
    fi
}

JOB_COUNT=0
for DS in "${DATASETS[@]}"; do
    for GIDX in 0 1 2; do
        bwa_map_one "${DS}" "${GIDX}" &
        JOB_COUNT=$((JOB_COUNT + 1))
        if [ $((JOB_COUNT % PARALLEL)) -eq 0 ]; then
            wait
            echo "[$(date)] BWA: ${JOB_COUNT}/45 jobs done"
        fi
    done
done
wait
echo "[$(date)] BWA mapping complete (${JOB_COUNT} jobs)"

# =====================================================================
# PHASE 5: QCMapping Extraction (reuses existing script)
# =====================================================================
echo ""
echo "========== PHASE 5: QCMapping Extraction =========="

for DS in "${DATASETS[@]}"; do
    QC_DIR="${BASEDIR}/${DS}/qc"
    mkdir -p "${QC_DIR}"

    python "${REPO_ROOT}/workflows/03_genotyping/synthetic/03_07_extract_metrics.py" \
        --bam-dir "${MAPPING_DIR}/${DS}" \
        --outdir "${QC_DIR}" \
        --genome-map B73v5:B73 Il14H:Il14H Ki11:Ki11
done
echo "[$(date)] Phase 5 complete"

# =====================================================================
# PHASE 6: AmbientMapper Chunks + Assign
# =====================================================================
echo ""
echo "========== PHASE 6: Chunks + Assign =========="

WORKDIR_ABS="${PROJECT_ROOT}/5_AmbientDetection/${BASEDIR}"

for DS in "${DATASETS[@]}"; do
    echo "[$(date)] Phase 6: ${DS}"
    DS_DIR="${BASEDIR}/${DS}"

    # Copy QCMapping to filtered_QCFiles/
    FILT_DIR="${DS_DIR}/filtered_QCFiles"
    mkdir -p "${FILT_DIR}"
    for G in "${GENOMES[@]}"; do
        cp "${DS_DIR}/qc/${G}_QCMapping.txt" "${FILT_DIR}/filtered_${G}_QCMapping.txt"
    done

    # Auto-generate config.json
    CONFIG="${DS_DIR}/config.json"
    BAM_DIR="${WORKDIR_ABS}/../${MAPPING_DIR}/${DS}"
    cat > "${CONFIG}" <<CONFIGEOF
{
    "sample": "${DS}",
    "workdir": "${WORKDIR_ABS}",
    "min_barcode_freq": 3,
    "chunk_size_cells": 100,
    "genomes": {
        "B73": "${BAM_DIR}/syn_B73v5.sorted.bam",
        "Il14H": "${BAM_DIR}/syn_Il14H.sorted.bam",
        "Ki11": "${BAM_DIR}/syn_Ki11.sorted.bam"
    }
}
CONFIGEOF

    # Chunks
    ambientmapper chunks --config "${CONFIG}" --chunk-size-cells 100

    # Assign
    ambientmapper assign \
        --config "${CONFIG}" \
        --threads "${THREADS}" \
        --alpha 0.05 \
        --k 10 \
        --mapq-min 10 \
        --xa-max 2 \
        --chunksize 500000 \
        --batch-size 6
done

# =====================================================================
# SUMMARY
# =====================================================================
echo ""
echo "============================================"
echo "[$(date)] Discriminative pipeline complete"
echo "============================================"
echo ""
echo "Datasets processed:"
for DS in "${DATASETS[@]}"; do
    N=$(ls "${BASEDIR}/${DS}/cell_map_ref_chunks/"*_filtered.tsv.gz 2>/dev/null | wc -l)
    echo "  ${DS}: ${N} filtered chunks"
done
echo ""
echo "Next: 03_21a_synthetic_friend_relabel.sh, then 03_21_genotyping_synthetic_factorial.sh"
