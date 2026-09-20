#!/bin/bash
########## SBATCH Resource Request ##########
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=PrepVCF
#SBATCH --partition=standard
#SBATCH --output=_logs/%x-%j.log

# ===========================================================================
# One-time VCF prep: bgzip + tabix the Mo17-vs-B73 VCF.
#
# Usage: cd ${PROJECT_ROOT}/5_Genotyping && sbatch <repo>/workflows/03b_variant_based_comparison/06_00_prep_mo17_vcf.sh
#
# Souporcell (06_01) reads the plain VCF; the WASP and purity steps (06_46 to
# 06_48) read the indexed .vcf.gz produced here.
# Input : configs/Final_Mo17_relative_to_B73.vcf (data, not shipped, see README)
# Output: configs/Final_Mo17_relative_to_B73.vcf.gz + .tbi
# ===========================================================================

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"   # root of the data tree (3_Mapping/, 5_AmbientDetection/, ...)

PROJ="${PROJECT_ROOT}/5_Genotyping"
VCF_IN="${PROJ}/configs/Final_Mo17_relative_to_B73.vcf"
VCF_OUT="${PROJ}/configs/Final_Mo17_relative_to_B73.vcf.gz"

ml Bioinformatics bcftools/1.12-g4b275e

if [[ ! -f "$VCF_IN" ]]; then
    echo "Error: VCF not found: $VCF_IN" >&2
    exit 1
fi

if [[ -f "$VCF_OUT" && -f "${VCF_OUT}.tbi" ]]; then
    echo "Already done: $VCF_OUT and ${VCF_OUT}.tbi exist"
    exit 0
fi

mkdir -p "${PROJ}/_logs"

echo "bgzip + tabix: $VCF_IN"
bgzip -c "$VCF_IN" > "$VCF_OUT"
tabix -p vcf "$VCF_OUT"

echo "Done: $(ls -lh "$VCF_OUT" | awk '{print $5}') compressed"
