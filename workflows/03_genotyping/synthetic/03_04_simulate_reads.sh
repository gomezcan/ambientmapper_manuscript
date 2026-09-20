#!/bin/bash
#SBATCH --job-name=art_simulate
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=_logs/03_04_simulate_reads_%j.log

# =============================================================================
# 03_04_simulate_reads.sh — Phase 2: Simulate reads from orthologous peaks
#
# Uses ART to generate paired-end Illumina reads from B73, Il14H, and Ki11
# peak sequences. Reads carry realistic quality scores and error profiles.
#
# Requires: Phase 1 output (synthetic/orthologs/{B73,Il14H,Ki11}_peak_seqs.fa)
# Output:   synthetic/reads/{genome}_{1,2}.fq.gz per genome
# =============================================================================

PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"   # root of the data tree (3_Mapping/, 5_AmbientDetection/, ...)
cd "${PROJECT_ROOT}/5_AmbientDetection"

# conda activate <env providing art_illumina (ART)>

# --- Parameters ---
ORTHO_DIR="synthetic/orthologs"
OUTDIR="synthetic/reads"
READ_LEN=75           # paired-end read length (bp)
FOLD_COV=100          # fold coverage per genome (~3.8M reads each)
INSERT_MEAN=166       # mean fragment size (from MACS3, Root1_rep1)
INSERT_SD=30          # insert size std dev
SEED=42               # reproducibility

GENOMES=(B73 Il14H Ki11)

# --- Verify inputs ---
for G in "${GENOMES[@]}"; do
    FA="${ORTHO_DIR}/${G}_peak_seqs.fa"
    if [ ! -f "${FA}" ]; then
        echo "Error: ${FA} not found. Run 03_02_find_orthologs.sh first."
        exit 1
    fi
done

mkdir -p "${OUTDIR}"

echo "============================================"
echo "[$(date)] Phase 2: ART read simulation"
echo "  Input:       ${ORTHO_DIR}/{B73,Il14H,Ki11}_peak_seqs.fa"
echo "  Read length: ${READ_LEN}bp PE"
echo "  Coverage:    ${FOLD_COV}x per genome"
echo "  Insert:      ${INSERT_MEAN} +/- ${INSERT_SD} bp"
echo "  Output:      ${OUTDIR}"
echo "============================================"

for G in "${GENOMES[@]}"; do
    FA="${ORTHO_DIR}/${G}_peak_seqs.fa"
    PREFIX="${OUTDIR}/${G}_"

    echo ""
    echo "[$(date)] Simulating reads for ${G}..."
    N_SEQS=$(grep -c "^>" "${FA}")
    TOTAL_BP=$(grep -v "^>" "${FA}" | tr -d '\n' | wc -c)
    echo "  Input: ${N_SEQS} sequences, ${TOTAL_BP} bp total"

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

    if [ $? -ne 0 ]; then
        echo "ERROR: ART failed for ${G}"
        exit 1
    fi

    # Count reads and compress
    N_READS=$(grep -c "^@" "${PREFIX}1.fq" 2>/dev/null || echo 0)
    echo "  Generated: ${N_READS} read pairs"

    echo "  Compressing..."
    gzip -f "${PREFIX}1.fq" &
    gzip -f "${PREFIX}2.fq" &
    wait

    echo "  Output: ${PREFIX}{1,2}.fq.gz"
done

# --- Summary ---
echo ""
echo "=== Read simulation summary ==="
for G in "${GENOMES[@]}"; do
    R1="${OUTDIR}/${G}_1.fq.gz"
    if [ -f "${R1}" ]; then
        SIZE=$(ls -lh "${R1}" | awk '{print $5}')
        echo "  ${G}: ${R1} (${SIZE})"
    fi
done

echo ""
echo "============================================"
echo "[$(date)] Phase 2 complete."
echo "  Output: ${OUTDIR}/{B73,Il14H,Ki11}_{1,2}.fq.gz"
echo "============================================"
