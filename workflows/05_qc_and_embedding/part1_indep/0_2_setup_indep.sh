#!/usr/bin/env bash
# 0_2_setup_indep.sh
#
# One-shot setup for the SM2v2 SPLIT-genome (independent-mapping) Socrates run.
# Stages the three things 1_0_qc_run.sh needs that don't exist yet:
#   (1) SM2-specific chr_sizes, derived from the AM-input BAM @SQ headers so the
#       contig names match the BEDs exactly (B73v5 keeps Mt:AGPv4/Pt:AGPv4; the
#       generic _GenomeInfo/B73v5.chrs.size.txt from B73Mo17 would NOT match).
#   (2) TAIR10 annotation symlink -> Ensembl TAIR10 r60 gff3 (chroms 1-5/Mt/Pt,
#       an exact match to the independent TAIR10 BED; no fixChrNames needed).
#   (3) PostClean (nd, design-free) per-genome BED symlinks into _data/_BED_files/,
#       reconciling the file's Clean-suffix name to the Clean.-prefix convention.
#   (3b) PostClean (wd, design-guided) BED symlinks, named Clean.SM2v2wd_* so they
#       coexist with the nd links (sensitivity fast-follow).
#
# Run interactively on HPC (needs samtools; takes seconds). No SBATCH.
#
#   bash 0_scripts/part1_indep/0_2_setup_indep.sh

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-.}"
BASE="${PROJECT_ROOT}/6_socrates"
GENOME_INFO="${BASE}/_data/_GenomeInfo"
BEDDIR="${BASE}/_data/_BED_files"

# Upstream inputs. AMIN = the AmbientMapper input BAMs (mapping-stage output). The two GFF3
# paths point at genome annotations outside this project tree and are taken from the environment:
#   AT_GFF3    Arabidopsis_thaliana.TAIR10.60.clean.gff3 (Ensembl Plants release 60)
#   B73V5_GFF3 Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 (MaizeGDB)
AMIN="${AMIN:-${PROJECT_ROOT}/3_Mapping/ambientmapper_input}"
AT_GFF3="${AT_GFF3:-}"

# conda activate ambientmapper-manuscript   (environment.yml at the repo root; provides samtools)
command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools not on PATH"; exit 2; }

mkdir -p "${GENOME_INFO}"

# ---- (1) chr_sizes from AM-input BAM @SQ headers ----------------------------
declare -A BAM
BAM[B73v5]="${AMIN}/SM2_B73v5_scifiATAC.mq10.BC.rmdup.mm.bam"
BAM[TAIR10]="${AMIN}/SM2_TAIR10_scifiATAC.mq10.BC.rmdup.mm.bam"
for g in B73v5 TAIR10; do
  bam="${BAM[$g]}"; out="${GENOME_INFO}/SM2v2_${g}.chrs.size.txt"
  [[ -f "$bam" ]] || { echo "WARN: missing AM-input BAM for ${g}: ${bam}"; continue; }
  echo " - chr_sizes ${g}  <-  ${bam##*/}"
  samtools view -H "$bam" \
    | awk -F'\t' '$1=="@SQ"{sn="";ln="";for(i=2;i<=NF;i++){if($i~/^SN:/)sn=substr($i,4);else if($i~/^LN:/)ln=substr($i,4)} if(sn!=""&&ln!="")printf "%s\t%s\n",sn,ln}' > "$out"
  echo "     $(wc -l < "$out") contigs -> ${out}"
done

# ---- (2) TAIR10 annotation symlink (Ensembl r60 gff3; chroms already 1-5/Mt/Pt)
if [[ -f "$AT_GFF3" ]]; then
  ln -sfn "${AT_GFF3}" "${GENOME_INFO}/TAIR10.gff3"
  echo " - symlink TAIR10.gff3 -> ${AT_GFF3}"
else
  echo "WARN: TAIR10 gff3 not found at ${AT_GFF3}"
fi

# ---- (2b) B73v5 annotation symlink (original NAM5 gff3 — proper ID/Parent).
#           The *.gtf variants (fixed_gene_id.gtf, Zm.gtf) carry GFF3 keys in GTF
#           syntax with NO gene_id -> makeTxDbFromGFF(format=gtf) fails at TxDb build.
B73V5_GFF3="${B73V5_GFF3:-}"
if [[ -f "$B73V5_GFF3" ]]; then
  ln -sfn "${B73V5_GFF3}" "${GENOME_INFO}/B73v5.gff3"
  echo " - symlink B73v5.gff3 -> ${B73V5_GFF3}"
else
  echo "WARN: B73v5 gff3 not found at ${B73V5_GFF3}"
fi

# ---- (3) PostClean (nd) BED symlinks (relative -> portable across mounts) ----
for g in B73v5 TAIR10; do
  rel="../../../5_AmbientDetection/SM2v2/clean_bams_alpha05_C0_nd/SM2_${g}_scifiATAC.mq10.BC.rmdup.mm.Clean.tn5.bed.gz"
  dst="${BEDDIR}/Clean.SM2v2_${g}_scifiATAC.mq10.BC.rmdup.mm.tn5.bed.gz"
  [[ -f "${BEDDIR}/${rel}" ]] || { echo "WARN: clean nd BED missing for ${g}: ${BEDDIR}/${rel}"; continue; }
  ln -sfn "$rel" "$dst"
  echo " - symlink ${dst##*/}  ->  ${rel}"
done

# ---- (3b) PostClean (wd, design-guided) BED symlinks — sensitivity fast-follow ----
#           Named Clean.SM2v2wd_* so they coexist with the nd links above.
for g in B73v5 TAIR10; do
  rel="../../../5_AmbientDetection/SM2v2/clean_bams_alpha05_C0_wd/SM2_${g}_scifiATAC.mq10.BC.rmdup.mm.Clean.tn5.bed.gz"
  dst="${BEDDIR}/Clean.SM2v2wd_${g}_scifiATAC.mq10.BC.rmdup.mm.tn5.bed.gz"
  [[ -f "${BEDDIR}/${rel}" ]] || { echo "WARN: clean wd BED missing for ${g}: ${BEDDIR}/${rel}"; continue; }
  ln -sfn "$rel" "$dst"
  echo " - symlink ${dst##*/}  ->  ${rel}"
done

echo
echo "Done. Verify:"
echo "  ls -l ${GENOME_INFO}/SM2v2_*.chrs.size.txt ${GENOME_INFO}/TAIR10.gff3"
echo "  ls -lL ${BEDDIR}/SM2_{B73v5,TAIR10}_*.tn5.bed.gz ${BEDDIR}/Clean.SM2v2{,wd}_{B73v5,TAIR10}_*.tn5.bed.gz"
