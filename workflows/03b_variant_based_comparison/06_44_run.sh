#!/usr/bin/env bash
#SBATCH --job-name=06_44_conc
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=36:00:00
#SBATCH --output=_logs/06_44_%x_%j.log
#
# 2d (count-based) PRIMARY — reference-anchored allele concordance, RAW vs CLEAN.
#
# Counts reads carrying each founder allele at KNOWN reference-VCF sites
# (genome-wide, all reads, no SNP discovery / MAF / LD / missingness filters),
# then estimates per-(group x founder) concordance + GRM-style correlation with
# block-bootstrap CIs.
#
# Manuscript role: the concordance panel produced by STEP=analyze is NOT in the
# manuscript. The prep() step is load-bearing: it builds the marker-site panel
# (multi_7geno_sites.vcf.gz for multiGenotypes_rep1, genotype_panel.tsv) and the
# barcode-to-genome table (bc_to_genome1.tsv for the B73Mo17 replicates) that the
# WASP and purity chain (06_46 to 06_48, Fig. 4L to N, Table S4) reads.
#
# Pseudo-bulk labels (the only per-sample difference):
#   B73Mo17_rep1/rep2   -> AM genome_1 'top-1' call   (bc_key_map mode)
#   multiGenotypes_rep1 -> plate-of-origin, ground truth (plate mode)
#
# Usage (extract is per SAMPLE+MODE; analyze is per SAMPLE, reads both modes):
#   sbatch --export=ALL,STEP=extract,SAMPLE=B73Mo17_rep1,MODE=raw   workflows/03b_variant_based_comparison/06_44_run.sh
#   sbatch --export=ALL,STEP=extract,SAMPLE=B73Mo17_rep1,MODE=clean workflows/03b_variant_based_comparison/06_44_run.sh
#   sbatch --export=ALL,STEP=analyze,SAMPLE=B73Mo17_rep1           workflows/03b_variant_based_comparison/06_44_run.sh
#
# Output: 5_Genotyping/<SAMPLE>/diagnostics/06_44_concordance/
# =============================================================================
set -euo pipefail

REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"   # this repository (sourced helpers, sibling scripts)
PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"   # root of the data tree (3_Mapping/, 5_AmbientDetection/, ...)

STEP="${STEP:?need STEP=extract|analyze}"
SAMPLE="${SAMPLE:?need SAMPLE}"

PROJ=${PROJECT_ROOT}
GENO=$PROJ/5_Genotyping
AM=$PROJ/5_AmbientDetection
MAP=$PROJ/3_Mapping
cd "$GENO"

# conda activate <env from environment.yml>
# NOTE: do NOT `ml Bioinformatics` — it shadows base python + injects a py3.9
# numpy via PYTHONPATH that crashes base py3.12.

NCPU="${SLURM_CPUS_PER_TASK:-16}"
CHROMS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10"   # exclude Mt/Pt
MIN_CELL_READS="${MIN_CELL_READS:-200}"   # B73Mo17 cell def: founder genome_1 & n_reads>=this
OUT=$GENO/${SAMPLE}/diagnostics/06_44_concordance
mkdir -p "$OUT"

# --- per-sample dispatch ----------------------------------------------------
case "$SAMPLE" in
  B73Mo17_rep1|B73Mo17_rep2)
    GROUP_MODE=bc_key_map
    GENOTYPES="B73,Mo17"
    EXTRACT_VCF=$GENO/configs/Final_Mo17_relative_to_B73.vcf.gz   # already B73-vs-Mo17 only
    AM_CALLS=$AM/${SAMPLE}/genotyping_runs/4cfg_2026-05-01/C0/${SAMPLE}_cells_calls.tsv.gz
    BC_MAP=$OUT/bc_to_genome1.tsv
    ;;
  multiGenotypes_rep1)
    GROUP_MODE=plate
    GENOTYPES="B73,B97,Ky21,M162W,Mo18W,Oh7B,Tzi8"
    WELL_MAP=$GENO/configs/Well_to_Genotype_multiGenotypes_rep1.txt
    EXTRACT_VCF=$OUT/multi_7geno_sites.vcf.gz                     # built below (polymorphic among 7)
    ;;
  *) echo "Unknown SAMPLE: $SAMPLE"; exit 1 ;;
esac

# --- shared prep: extraction VCF + group map --------------------------------
# Race-safe: raw+clean for one sample share these per-sample files. Exactly one
# job builds them (atomic temp+mv) under a lockdir; concurrent jobs wait for the
# sentinel. Re-runs skip if the sentinel exists.
prep() {
  local sentinel="$OUT/.prep.done"
  local lock="$OUT/.prep.lockdir"
  [[ -e "$sentinel" ]] && return 0   # -e (exists), NOT -s: sentinel is touch'd (0 bytes)
  if mkdir "$lock" 2>/dev/null; then
    trap 'rmdir "'"$lock"'" 2>/dev/null || true' EXIT
    if [[ "$SAMPLE" == multiGenotypes_rep1 && ! -s "$EXTRACT_VCF" ]]; then
      echo "[06_44] building 7-genotype polymorphic sites VCF"
      bcftools view -s "$GENOTYPES" -m2 -M2 -v snps --min-af 0.01:minor \
        -Oz -o "$EXTRACT_VCF.tmp" "$GENO/configs/25NAM_full.vcf.gz"
      tabix -p vcf "$EXTRACT_VCF.tmp"
      mv -f "$EXTRACT_VCF.tmp.tbi" "$EXTRACT_VCF.tbi"
      mv -f "$EXTRACT_VCF.tmp" "$EXTRACT_VCF"
      echo "[06_44]   sites: $(zcat "$EXTRACT_VCF" | grep -vc '^#')"
    fi
    if [[ "$GROUP_MODE" == bc_key_map && ! -s "$BC_MAP" ]]; then
      echo "[06_44] building bc_key -> AM genome_1 map (founder g1 & n_reads>=$MIN_CELL_READS)"
      { echo -e "bc_key\tgroup\tcall";
        zcat "$AM_CALLS" | tail -n +2 | awk -F'\t' -v mr="$MIN_CELL_READS" 'BEGIN{OFS="\t"}
          ($4=="B73"||$4=="Mo17") && ($8+0)>=mr { split($1,a,"-"); print a[1],$4,$2 }'
      } > "$BC_MAP.tmp"
      mv -f "$BC_MAP.tmp" "$BC_MAP"
      echo "[06_44]   barcodes per group:"; tail -n +2 "$BC_MAP" | cut -f2 | sort | uniq -c
    fi
    touch "$sentinel"
    rmdir "$lock" 2>/dev/null || true
    trap - EXIT
  else
    echo "[06_44] another job is building prep; waiting for $sentinel"
    local i=0
    while [[ ! -e "$sentinel" ]]; do
      sleep 5; i=$((i + 1))
      [[ $i -gt 360 ]] && { echo "[06_44] ERROR: prep wait timed out (30 min)"; exit 1; }
    done
    echo "[06_44] prep ready"
  fi
}

# =============================================================== STEP: extract
if [[ "$STEP" == "extract" ]]; then
  MODE="${MODE:?extract needs MODE=raw|clean}"
  case "$MODE" in
    raw)   BAM=$MAP/ambientmapper_input/${SAMPLE}_B73v5_scifiATAC.mq10.BC.rmdup.mm.bam ;;
    clean) BAM=$AM/${SAMPLE}/clean_bams_alpha05_C0_nd/${SAMPLE}_B73v5_scifiATAC.mq10.BC.rmdup.mm.Clean.bam ;;
    *) echo "MODE must be raw|clean"; exit 1 ;;
  esac
  [[ -f $BAM ]] || { echo "Missing BAM: $BAM"; exit 1; }
  prep
  OUTCOUNTS=$OUT/${SAMPLE}_${MODE}_allele_counts.tsv.gz
  echo "[06_44] EXTRACT sample=$SAMPLE mode=$MODE group_mode=$GROUP_MODE"
  echo "[06_44]   BAM=$BAM"
  if [[ "$GROUP_MODE" == bc_key_map ]]; then
    python "${REPO_ROOT}/workflows/03b_variant_based_comparison/06_44_extract_allele_counts.py" \
      --bam "$BAM" --vcf "$EXTRACT_VCF" --group-mode bc_key_map --bc-map "$BC_MAP" \
      --bc-tag BC --min-mapq 10 --threads "$NCPU" --chroms "$CHROMS" --out "$OUTCOUNTS"
  else
    python "${REPO_ROOT}/workflows/03b_variant_based_comparison/06_44_extract_allele_counts.py" \
      --bam "$BAM" --vcf "$EXTRACT_VCF" --group-mode plate --well-map "$WELL_MAP" \
      --bc-tag BC --min-mapq 10 --threads "$NCPU" --chroms "$CHROMS" --out "$OUTCOUNTS"
  fi
  echo "[06_44] EXTRACT done -> $OUTCOUNTS"; exit 0
fi

# =============================================================== STEP: analyze
if [[ "$STEP" == "analyze" ]]; then
  prep
  RAW=$OUT/${SAMPLE}_raw_allele_counts.tsv.gz
  CLN=$OUT/${SAMPLE}_clean_allele_counts.tsv.gz
  for f in "$RAW" "$CLN"; do [[ -s $f ]] || { echo "Missing counts (run extract first): $f"; exit 1; }; done
  PANEL=$OUT/genotype_panel.tsv
  if [[ ! -s $PANEL ]]; then
    echo "[06_44] building genotype panel ($GENOTYPES) from extraction VCF"
    bcftools query -s "$GENOTYPES" -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$EXTRACT_VCF" > "$PANEL"
  fi
  echo "[06_44] ANALYZE sample=$SAMPLE"
  Rscript "${REPO_ROOT}/workflows/03b_variant_based_comparison/06_44_concordance.R" \
    --counts-raw "$RAW" --counts-clean "$CLN" \
    --genotype-panel "$PANEL" --genotypes "$GENOTYPES" \
    --out-prefix "$OUT/${SAMPLE}" \
    --min-reads 2 --n-boot 200 --block-mb 1
  echo "[06_44] ANALYZE done -> $OUT/${SAMPLE}_*"; exit 0
fi

echo "Unknown STEP: $STEP (need extract|analyze)"; exit 1
