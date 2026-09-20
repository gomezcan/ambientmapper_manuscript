#!/bin/bash
########## SBATCH Resource Request ##########
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=14
#SBATCH --mem=160G
#SBATCH --job-name=Geno.soup
#SBATCH --partition=standard
#SBATCH --output=_logs/%x-%j.log

# ===========================================================================
# Souporcell (v2.1) supervised genotyping for the Zhang 2024 scifi-ATAC libraries,
# followed by cluster-to-genotype assignment (06_03). Its clusters.tsv and
# Genotype_ID_key.v2.txt are the SNP-based calls compared with AmbientMapper in
# Fig. 4H to K.
#
# Usage: sbatch 06_01_souporcell.sh <SAMPLE> [min_frags]
#   SAMPLE:    B73Mo17_rep1, B73Mo17_rep2, or multiGenotypes_rep1
#   min_frags: Minimum fragments per barcode (default: 500)
#
# VCF routing:
#   B73Mo17*         -> Final_Mo17_relative_to_B73.vcf (B73 + Mo17)
#   multiGenotypes*  -> 25NAM_full.vcf.gz (filtered to 7 genotypes)
#
# Output: <SAMPLE>/souporcell/supervised/<SAMPLE>.min<N>/
# ===========================================================================

set -euo pipefail

REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"   # this repository (sourced helpers, sibling scripts)
PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"   # root of the data tree (3_Mapping/, 5_AmbientDetection/, ...)

# --- Arguments ---
SAMPLE="${1:-}"
MIN_FRAGS="${2:-500}"

if [[ -z "$SAMPLE" ]]; then
    echo "Usage: $0 <SAMPLE> [min_frags]" >&2
    echo "  SAMPLE:    B73Mo17_rep1, B73Mo17_rep2, or multiGenotypes_rep1" >&2
    echo "  min_frags: Min fragments per barcode (default: 500)" >&2
    exit 1
fi

# --- Paths ---
PROJ="${PROJECT_ROOT}/5_Genotyping"
# Souporcell v2.1 Singularity image and the B73 NAM 5.0 reference (MaizeGDB).
LOCALINSTALL="${LOCALINSTALL:-${HOME}/LocalInstall}"
SOUPORCELL_SIF="${SOUPORCELL_SIF:-${LOCALINSTALL}/souporcell_release.sif}"
GENOMES_ROOT="${GENOMES_ROOT:?set GENOMES_ROOT to the directory holding Zea/Zm_B73_REFERENCE_NAM_5.0/}"
FASTA="${GENOMES_ROOT}/Zea/Zm_B73_REFERENCE_NAM_5.0/Zm-B73-REFERENCE-NAM-5.0.chrs.mt.pt.fa"
SAMPLE_DB="${PROJ}/configs/Pools_DB.by_Sample.txt"
INPUT_DIR="${PROJECT_ROOT}/3_Mapping/ambientmapper_input"
THREADS="${SLURM_CPUS_PER_TASK:-14}"

# Output directory
OUT_DIR="${PROJ}/${SAMPLE}/souporcell/supervised/${SAMPLE}.min${MIN_FRAGS}"

# --- VCF routing ---
if [[ "$SAMPLE" == B73Mo17* ]]; then
    VCF_FILE="${PROJ}/configs/Final_Mo17_relative_to_B73.vcf"
elif [[ "$SAMPLE" == multiGenotypes* ]]; then
    VCF_FILE="${PROJ}/configs/25NAM_full.vcf.gz"
else
    echo "Error: Unknown sample '$SAMPLE'. Expected B73Mo17* or multiGenotypes*" >&2
    exit 1
fi

# --- Load modules ---
ml Bioinformatics vcftools/0.1.15 bcftools/1.12-g4b275e singularity

# --- Validate inputs ---
BAM="${INPUT_DIR}/${SAMPLE}_B73v5_scifiATAC.mq10.BC.rmdup.mm.bam"

if [[ ! -f "$BAM" ]]; then
    echo "Error: BAM file not found: $BAM" >&2
    exit 1
fi

if [[ ! -f "$SOUPORCELL_SIF" ]]; then
    echo "Error: SouporCell image not found: $SOUPORCELL_SIF" >&2
    exit 1
fi

if [[ ! -f "$SAMPLE_DB" ]]; then
    echo "Error: Sample database not found: $SAMPLE_DB" >&2
    exit 1
fi

if [[ ! -f "$VCF_FILE" ]]; then
    echo "Error: VCF not found: $VCF_FILE" >&2
    exit 1
fi

# --- Setup ---
mkdir -p "$OUT_DIR" "${PROJ}/_logs"

echo "============================================"
echo "SouporCell genotyping — supervised"
echo "Sample:    $SAMPLE"
echo "Min frags: $MIN_FRAGS"
echo "BAM:       $BAM"
echo "VCF:       $VCF_FILE"
echo "Output:    $OUT_DIR"
echo "============================================"

# --- BAM: filter to chr1-10 + convert BC->CB for SouporCell ---
# Filter out Mt/Pt/scaffold contigs that cause naming conflicts.
# Shared across thresholds — store at sample/souporcell level
CHROMS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10"
BAM_CB="${PROJ}/${SAMPLE}/souporcell/supervised/${SAMPLE}.CB.bam"

if [[ ! -f "$BAM_CB" ]]; then
    echo "Filtering to chr1-10, converting BC:Z -> CB:Z, and sorting..."
    mkdir -p "$(dirname "$BAM_CB")"
    samtools view -@ "$THREADS" -h "$BAM" $CHROMS | sed 's/\tBC:Z:/\tCB:Z:/g' | \
        samtools sort -@ "$THREADS" -o "$BAM_CB"
    samtools index -@ "$THREADS" "$BAM_CB"
else
    echo "CB BAM already exists: $BAM_CB"
fi

# --- Barcodes: extract passing barcodes at requested threshold ---
BC_FILE="${OUT_DIR}/${SAMPLE}.min${MIN_FRAGS}.barcodes.txt"

if [[ ! -f "$BC_FILE" ]]; then
    echo "Extracting barcodes (min $MIN_FRAGS fragments)..."
    samtools view -@ 4 "$BAM_CB" | \
        grep -o 'CB:Z:[^[:space:]]*' | cut -d':' -f3- | \
        sort | uniq -c | awk -v min="$MIN_FRAGS" '{if($1>=min) print $2}' > "$BC_FILE"
    echo "  Found $(wc -l < "$BC_FILE") barcodes"
else
    echo "Barcode file already exists: $BC_FILE ($(wc -l < "$BC_FILE") barcodes)"
fi

# --- Pool genotypes: extract expected genotypes for this sample ---
# Pools_DB.by_Sample.txt format: genotype<TAB>sample (no header)
GENO_FILE="${PROJ}/${SAMPLE}/souporcell/supervised/${SAMPLE}.genotypes.txt"

if [[ ! -f "$GENO_FILE" ]]; then
    mkdir -p "$(dirname "$GENO_FILE")"
    awk -F'\t' -v s="$SAMPLE" '$2 == s {print $1}' "$SAMPLE_DB" | sort -u > "$GENO_FILE"
fi

NUM_K=$(wc -l < "$GENO_FILE")
GENOS=$(tr '\n' ' ' < "$GENO_FILE" | sed 's/ $//')
echo "Sample: k=$NUM_K | genotypes: $GENOS"

if [[ "$NUM_K" -eq 0 ]]; then
    echo "Error: No genotypes found for sample $SAMPLE in $SAMPLE_DB" >&2
    exit 1
fi

# --- Filter VCF to sample genotypes ---
FILTERED_VCF="${PROJ}/${SAMPLE}/souporcell/supervised/${SAMPLE}.filtered.vcf"

if [[ ! -f "$FILTERED_VCF" ]]; then
    echo "Filtering VCF to sample genotypes..."
    bcftools view --threads "$THREADS" -S "$GENO_FILE" "$VCF_FILE" -o "$FILTERED_VCF"
else
    echo "Filtered VCF already exists: $FILTERED_VCF"
fi

# --- Run SouporCell (supervised mode) ---
# Resume marker: consensus.done is the atomic marker SouporCell writes LAST.
if [[ ! -f "$OUT_DIR/consensus.done" || ! -s "$OUT_DIR/cluster_genotypes.vcf" ]]; then
    echo "Running SouporCell (supervised mode, min_frags=$MIN_FRAGS)..."

    # --cleanenv prevents conda compiler conflict inside container
    singularity exec --cleanenv "$SOUPORCELL_SIF" souporcell_pipeline.py \
        --bam "$BAM_CB" \
        --barcodes "$BC_FILE" \
        --fasta "$FASTA" \
        --threads "$THREADS" \
        --no_umi True \
        --cell_tag CB \
        --out_dir "$OUT_DIR" \
        -k "$NUM_K" \
        --known_genotypes "$FILTERED_VCF" \
        --known_genotypes_sample_names $GENOS \
        --skip_remap SKIP_REMAP
else
    echo "SouporCell already completed: $OUT_DIR/cluster_genotypes.vcf"
fi

# --- Genotype assignment via Pearson correlation ---
if [[ ! -f "$OUT_DIR/ref_clust_pearson_correlations.v2.tsv" ]]; then
    echo "Running genotype assignment (Pearson correlation + B73 rescue)..."
    Rscript "${REPO_ROOT}/workflows/03b_variant_based_comparison/06_03_assign_genotype_by_pearson.R" \
        "$FILTERED_VCF" \
        "$OUT_DIR/cluster_genotypes.vcf" \
        "$OUT_DIR"
else
    echo "Genotype assignment already completed"
fi

echo "============================================"
echo "Done: $SAMPLE (supervised, min_frags=$MIN_FRAGS)"
echo "Output: $OUT_DIR"
echo "============================================"
