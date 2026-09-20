#!/usr/bin/env bash
#SBATCH --job-name=06_48_purity
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G   # 120G (not 180): 180G ~= a full standard node -> queue-stalls waiting for an empty
                     # node. 06_44 ran this exact scan (same BAMs, same 438M/461M panel, 16 workers) at
                     # 120G; the (barcode x block) accumulator peaks ~15-20G on rep2. Per-sample override
                     # at submit if needed (multi is the largest; size from rep1/rep2 MaxRSS).
#SBATCH --time=4:00:00
#SBATCH --output=_logs/06_48_%x_%j.log
#
# 06_48 — barcode-resolved, depth- & call-class-stratified purity on the
# WASP-corrected BAMs (raw vs clean). Extends 06_44 from pooled (group x site)
# to (barcode x 1-Mb block); anchor = AM call (top1=genome_1), C1-C4 locked.
#
#   extract : scan a WASP BAM -> per-(barcode,block) match counts (heavy).
#   analyze : per sample, both modes -> stratified table + decomposition +
#             beta-binomial (glmmTMB, if installed) + ambiguous margins.
#
# Usage:
#   sbatch --export=ALL,STEP=extract,SAMPLE=B73Mo17_rep1,MODE=raw   workflows/03b_variant_based_comparison/06_48_run.sh
#   sbatch --export=ALL,STEP=extract,SAMPLE=B73Mo17_rep1,MODE=clean workflows/03b_variant_based_comparison/06_48_run.sh
#   sbatch --export=ALL,STEP=analyze,SAMPLE=B73Mo17_rep1           workflows/03b_variant_based_comparison/06_48_run.sh
#
# Output: 5_Genotyping/<SAMPLE>/diagnostics/06_48_barcode_purity/
#   <SAMPLE>_barcode_purity.tsv.gz -> Fig. 4L to N;  <SAMPLE>_weak_doublet_diag.tsv -> Table S4
# =============================================================================
set -euo pipefail

REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"   # this repository (sourced helpers, sibling scripts)
PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"   # root of the data tree (3_Mapping/, 5_AmbientDetection/, ...)

STEP="${STEP:?need STEP=extract|analyze}"
SAMPLE="${SAMPLE:?need SAMPLE}"

PROJ=${PROJECT_ROOT}
GENO=$PROJ/5_Genotyping
AM=$PROJ/5_AmbientDetection
cd "$GENO"

# conda activate <env from environment.yml>
# (do not `ml Bioinformatics` here: it shadows the base python 3.12.)

NCPU="${SLURM_CPUS_PER_TASK:-16}"
CHROMS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10"
MIN_CELL_READS="${MIN_CELL_READS:-200}"
CONC44=$GENO/${SAMPLE}/diagnostics/06_44_concordance   # reuse its panel + VCF
OUT=$GENO/${SAMPLE}/diagnostics/06_48_barcode_purity
mkdir -p "$OUT"

# --- per-sample dispatch ----------------------------------------------------
# WASP BAMs: rep1/multi chunked (06_47), rep2 monolithic (06_46).
case "$SAMPLE" in
  B73Mo17_rep1)
    GENOTYPES="B73,Mo17";  WASP_DIR=$GENO/$SAMPLE/diagnostics/06_47_wasp_chunked
    EXTRACT_VCF=$GENO/configs/Final_Mo17_relative_to_B73.vcf.gz ;;
  B73Mo17_rep2)
    GENOTYPES="B73,Mo17";  WASP_DIR=$GENO/$SAMPLE/diagnostics/06_46_wasp
    EXTRACT_VCF=$GENO/configs/Final_Mo17_relative_to_B73.vcf.gz ;;
  multiGenotypes_rep1)
    GENOTYPES="B73,B97,Ky21,M162W,Mo18W,Oh7B,Tzi8"
    WASP_DIR=$GENO/$SAMPLE/diagnostics/06_47_wasp_chunked
    EXTRACT_VCF=$CONC44/multi_7geno_sites.vcf.gz
    WELL_MAP=$GENO/configs/Well_to_Genotype_multiGenotypes_rep1.txt ;;
  *) echo "Unknown SAMPLE: $SAMPLE"; exit 1 ;;
esac
AM_CALLS=$AM/${SAMPLE}/genotyping_runs/4cfg_2026-05-01/C0/${SAMPLE}_cells_calls.tsv.gz
PANEL=$CONC44/genotype_panel.tsv          # GT-bearing, reused from 06_44
ATTR_MAP=$OUT/bc_attr_map.tsv

# --- shared prep: GT panel + barcode attr-map (race-safe temp+mv) -----------
prep() {
  local sentinel="$OUT/.prep.done" lock="$OUT/.prep.lockdir"
  [[ -e "$sentinel" ]] && return 0
  if mkdir "$lock" 2>/dev/null; then
    trap 'rmdir "'"$lock"'" 2>/dev/null || true' EXIT
    if [[ ! -s "$PANEL" ]]; then
      echo "[06_48] 06_44 panel absent — building GT panel ($GENOTYPES) from $EXTRACT_VCF"
      mkdir -p "$CONC44"
      bcftools query -s "$GENOTYPES" -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' \
        "$EXTRACT_VCF" > "$PANEL.tmp" && mv -f "$PANEL.tmp" "$PANEL"
    fi
    if [[ ! -s "$ATTR_MAP" ]]; then
      echo "[06_48] building bc attr-map (genome_1 in {$GENOTYPES} & n_reads>=$MIN_CELL_READS)"
      { echo -e "bc_key\tcall\ttop1\ttop2\traw_depth";
        zcat "$AM_CALLS" | tail -n +2 | awk -F'\t' -v mr="$MIN_CELL_READS" -v gl="$GENOTYPES" \
          'BEGIN{OFS="\t"; n=split(gl,g,","); for(i=1;i<=n;i++) F[g[i]]=1}
           ($4 in F) && ($8+0)>=mr { split($1,a,"-"); print a[1],$2,$4,$5,$8 }'
      } > "$ATTR_MAP.tmp"
      mv -f "$ATTR_MAP.tmp" "$ATTR_MAP"
      echo "[06_48]   barcodes per call class:"; tail -n +2 "$ATTR_MAP" | cut -f2 | sort | uniq -c
    fi
    touch "$sentinel"; rmdir "$lock" 2>/dev/null || true; trap - EXIT
  else
    echo "[06_48] waiting on concurrent prep ($sentinel)"; local i=0
    while [[ ! -e "$sentinel" ]]; do sleep 5; i=$((i+1));
      [[ $i -gt 360 ]] && { echo "[06_48] ERROR: prep wait timed out"; exit 1; }; done
  fi
}

# =============================================================== STEP: extract
if [[ "$STEP" == "extract" ]]; then
  MODE="${MODE:?extract needs MODE=raw|clean}"
  BAM=$WASP_DIR/${SAMPLE}_${MODE}.wasp.bam
  [[ -f $BAM ]]     || { echo "Missing WASP BAM: $BAM"; exit 1; }
  [[ -f $BAM.bai || -f ${BAM%.bam}.bai ]] || { echo "Missing .bai for $BAM"; exit 1; }
  prep
  WELL_ARG=(); [[ -n "${WELL_MAP:-}" ]] && WELL_ARG=(--well-map "$WELL_MAP")
  echo "[06_48] EXTRACT sample=$SAMPLE mode=$MODE (WASP) BAM=$BAM"
  python "${REPO_ROOT}/workflows/03b_variant_based_comparison/06_48_extract_barcode_block_counts.py" \
    --bam "$BAM" --genotype-panel "$PANEL" --genotypes "$GENOTYPES" \
    --bc-attr-map "$ATTR_MAP" "${WELL_ARG[@]}" \
    --bc-tag BC --min-mapq 10 --threads "$NCPU" --chroms "$CHROMS" \
    --out "$OUT/${SAMPLE}_${MODE}"
  echo "[06_48] EXTRACT done -> $OUT/${SAMPLE}_${MODE}.*"; exit 0
fi

# =============================================================== STEP: analyze
if [[ "$STEP" == "analyze" ]]; then
  prep
  for M in raw clean; do
    [[ -s "$OUT/${SAMPLE}_${M}.counts.tsv.gz" ]] || { echo "Missing ${M} counts (run extract)"; exit 1; }
  done
  echo "[06_48] ANALYZE sample=$SAMPLE"
  Rscript "${REPO_ROOT}/workflows/03b_variant_based_comparison/06_48_purity_model.R" \
    --raw-prefix "$OUT/${SAMPLE}_raw" --clean-prefix "$OUT/${SAMPLE}_clean" \
    --genotypes "$GENOTYPES" --sample "$SAMPLE" \
    --out-prefix "$OUT/${SAMPLE}" --n-boot 200
  echo "[06_48] ANALYZE done -> $OUT/${SAMPLE}_*"; exit 0
fi

echo "Unknown STEP: $STEP (need extract|analyze)"; exit 1
