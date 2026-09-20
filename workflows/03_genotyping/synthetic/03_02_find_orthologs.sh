#!/bin/bash
#SBATCH --job-name=find_orthologs
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=_logs/03_02_find_orthologs_%j.log

# =============================================================================
# 03_02_find_orthologs.sh — Phase 1: Find orthologous peak regions in Il14H/Ki11
#
# 1. Extract B73 peak sequences (bedtools getfasta)
# 2. Align to Il14H and Ki11 (minimap2 -x asm5)
# 3. Parse alignments, filter, extract orthologous sequences (Python)
#
# Requires: Phase 0 output (synthetic/peaks/B73_peaks_filtered.bed)
# Output:   synthetic/orthologs/{B73,Il14H,Ki11}_peak_seqs.fa + ortholog_map.tsv
# =============================================================================

REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"   # this repository (sourced helpers, sibling scripts)
PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"   # root of the data tree (3_Mapping/, 5_AmbientDetection/, ...)
cd "${PROJECT_ROOT}/5_AmbientDetection"

# conda activate <env from environment.yml>
ml Bioinformatics minimap2/2.14-jvscyil

# --- Paths ---
PEAKS="synthetic/peaks/B73_peaks_filtered.bed"
OUTDIR="synthetic/orthologs"
ALN_DIR="${OUTDIR}/alignments"
THREADS=8

# Reference FASTAs (MaizeGDB NAM assemblies). GENOMES_ROOT holds the Zea/ tree.
GENOMES_ROOT="${GENOMES_ROOT:?set GENOMES_ROOT to the directory holding Zea/Zm_<NAM>_REFERENCE_NAM_*/}"
B73_FASTA="${GENOMES_ROOT}/Zea/Zm_B73_REFERENCE_NAM_5.0/Zm-B73-REFERENCE-NAM-5.0.chrs.fa"
IL14H_FASTA="${GENOMES_ROOT}/Zea/Zm_Il14H_REFERENCE_NAM_1.0/Zm-Il14H-REFERENCE-NAM-1.0.fa"
KI11_FASTA="${GENOMES_ROOT}/Zea/Zm_Ki11_REFERENCE_NAM_1.0/Zm-Ki11-REFERENCE-NAM-1.0.fa"

# --- Verify inputs ---
if [ ! -f "${PEAKS}" ]; then
    echo "Error: Peaks file not found: ${PEAKS}"
    echo "Run 03_01_call_peaks.sh first."
    exit 1
fi

mkdir -p "${OUTDIR}" "${ALN_DIR}"

echo "============================================"
echo "[$(date)] Phase 1: Ortholog finding"
echo "  Peaks:   ${PEAKS} ($(wc -l < "${PEAKS}") peaks)"
echo "  Genomes: B73, Il14H, Ki11"
echo "  Output:  ${OUTDIR}"
echo "============================================"

# --- Step 1: Extract B73 peak sequences ---
echo "[$(date)] Extracting B73 peak sequences..."

# bedtools getfasta needs a FASTA index
if [ ! -f "${B73_FASTA}.fai" ]; then
    echo "  Creating FASTA index for B73..."
    samtools faidx "${B73_FASTA}"
fi

B73_PEAK_FA="${ALN_DIR}/B73_peaks_query.fa"
bedtools getfasta \
    -fi "${B73_FASTA}" \
    -bed "${PEAKS}" \
    -fo "${B73_PEAK_FA}" \
    -name

N_SEQS=$(grep -c "^>" "${B73_PEAK_FA}")
echo "  Extracted ${N_SEQS} peak sequences -> ${B73_PEAK_FA}"

# --- Step 2: Align B73 peaks to Il14H and Ki11 ---
for GENOME in Il14H Ki11; do
    if [ "${GENOME}" = "Il14H" ]; then
        TARGET_FASTA="${IL14H_FASTA}"
    else
        TARGET_FASTA="${KI11_FASTA}"
    fi

    # Index target if needed
    if [ ! -f "${TARGET_FASTA}.fai" ]; then
        echo "[$(date)] Creating FASTA index for ${GENOME}..."
        samtools faidx "${TARGET_FASTA}"
    fi

    SAM_OUT="${ALN_DIR}/B73_to_${GENOME}.sam"
    echo "[$(date)] Aligning B73 peaks to ${GENOME}..."

    minimap2 \
        -a \
        -x asm5 \
        -t "${THREADS}" \
        --secondary=no \
        "${TARGET_FASTA}" \
        "${B73_PEAK_FA}" \
        > "${SAM_OUT}" \
        2> "${ALN_DIR}/minimap2_${GENOME}.log"

    N_ALN=$(grep -cv "^@" "${SAM_OUT}")
    echo "  ${GENOME}: ${N_ALN} alignments -> ${SAM_OUT}"
done

# --- Step 3: Parse alignments and extract orthologous sequences ---
echo "[$(date)] Parsing alignments and extracting orthologs..."

python3 "${REPO_ROOT}/workflows/03_genotyping/synthetic/03_02_find_orthologs.py" \
    --peaks "${PEAKS}" \
    --b73-fasta "${B73_FASTA}" \
    --target-genomes Il14H Ki11 \
    --target-fastas "${IL14H_FASTA}" "${KI11_FASTA}" \
    --sam-dir "${ALN_DIR}" \
    --outdir "${OUTDIR}" \
    --min-coverage 0.90 \
    --min-identity 0.90

EXIT_CODE=$?

echo ""
echo "============================================"
echo "[$(date)] Phase 1 complete (exit code: ${EXIT_CODE})"
echo "  Ortholog map: ${OUTDIR}/ortholog_map.tsv"
echo "  B73 seqs:     ${OUTDIR}/B73_peak_seqs.fa"
echo "  Il14H seqs:   ${OUTDIR}/Il14H_peak_seqs.fa"
echo "  Ki11 seqs:    ${OUTDIR}/Ki11_peak_seqs.fa"
echo "============================================"

exit ${EXIT_CODE}
