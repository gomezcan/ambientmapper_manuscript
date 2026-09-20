#!/usr/bin/env bash
#SBATCH --job-name=06_46_wasp
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=72G
#SBATCH --time=48:00:00
#SBATCH --output=_logs/06_46_%x_%j.log
#
# WASP reference-mapping-bias correction — PILOT on B73Mo17_rep1 (raw + clean).
# =============================================================================
# WHY. 06_44 refbias.tsv shows a strong single-reference (B73) mapping bias:
# reads carrying the ALT (non-B73) allele carry ~4-4.8 mismatches vs ~0.4-1.6 for
# REF reads, so ALT coverage is depressed and every non-B73 genome's apparent
# B73(=REF) fraction is inflated. That is exactly the bias WASP removes: for each
# read overlapping a known SNP it flips the allele, remaps, and keeps the read
# ONLY if it returns to the same locus (van de Geijn et al. 2015, Nat Methods).
#
# GOAL of the pilot. Run WASP on BOTH the raw and the AM-clean BAM, THEN re-count
# with the SAME 06_44 extractor. Go/no-go read:
#   * refbias.tsv: Mo17 ALT-vs-REF mean-NM asymmetry (4.77 vs 1.66) collapses,
#     ALT observations rise (fewer dropped);
#   * Mo17 absolute self-concordance climbs from ~0.67 toward truth;
#   * the raw<->clean Delta stays ~ +0.08 (WASP and AM are orthogonal).
# If it works here, roll to rep2 + multi (multi is where the Mo18W argmax payoff
# lands). WASP is orthogonal to AM cleaning, so applying it to both arms and
# counting after cleanly isolates mapping-bias from the ambient-removal effect.
#
# ORTHOGONAL to 06_44: writes to <SAMPLE>/diagnostics/06_46_wasp/, reuses
# 06_44_extract_allele_counts.py + 06_44_concordance.R unchanged, and does NOT
# touch the validated 06_44_concordance/ outputs.
#
# Defaults (documented, change here if undesired):
#   * SKIP WASP's own rmdup    -- the BAMs are already Picard BC-deduped.
#   * Reuse 06_44 counting     -- min-mapq 10, chr1-10, BC tag -> directly
#                                 comparable wasp vs non-wasp numbers.
#
# Usage (setup once on a LOGIN node -- it git-clones WASP, needs internet):
#   STEP=setup SAMPLE=B73Mo17_rep1 bash workflows/03b_variant_based_comparison/06_46_wasp_pilot.sh
# then the two heavy jobs + analyze (sbatch, compute nodes):
#   sbatch --export=ALL,STEP=wasp,SAMPLE=B73Mo17_rep1,MODE=raw   workflows/03b_variant_based_comparison/06_46_wasp_pilot.sh
#   sbatch --export=ALL,STEP=wasp,SAMPLE=B73Mo17_rep1,MODE=clean workflows/03b_variant_based_comparison/06_46_wasp_pilot.sh
#   sbatch --export=ALL,STEP=analyze,SAMPLE=B73Mo17_rep1         workflows/03b_variant_based_comparison/06_46_wasp_pilot.sh
# =============================================================================
set -euo pipefail

REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"   # this repository (sourced helpers, sibling scripts)
PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"   # root of the data tree (3_Mapping/, 5_AmbientDetection/, ...)

STEP="${STEP:?need STEP=setup|wasp|analyze}"
SAMPLE="${SAMPLE:?need SAMPLE (pilot: B73Mo17_rep1)}"

PROJ=${PROJECT_ROOT}
GENO=$PROJ/5_Genotyping
AM=$PROJ/5_AmbientDetection
MAP=$PROJ/3_Mapping
cd "$GENO"

# --- env. base has pysam/numpy + bcftools/bgzip/tabix + R. bwa is NOT in base;
#     it lives in the Bioinformatics module (bwa/0.7.17-mil4ns7 = the 0.7.17 that
#     mapped these BAMs). `ml Bioinformatics` breaks base python TWO ways (shadows
#     `python` + injects a py3.9 numpy via PYTHONPATH that crashes base py3.12),
#     so the safe pattern is: load the module
#     for bwa, then `unset PYTHONPATH PYTHONHOME` and run every python/Rscript via
#     the ABSOLUTE base interpreter ($PY/$RSCRIPT). The wasp step (below) does the
#     module load; setup/analyze never need bwa so they stay pure-base. -----------
# conda activate <env from environment.yml>
PY="$CONDA_PREFIX/bin/python";   [[ -x "$PY" ]]      || PY=python
RSCRIPT="$CONDA_PREFIX/bin/Rscript"; [[ -x "$RSCRIPT" ]] || RSCRIPT=Rscript
BWA_MODULE="${BWA_MODULE:-bwa/0.7.17-mil4ns7}"
BWA_IDX="${BWA_IDX:-${BWA_INDEX_ROOT:?set BWA_INDEX_ROOT (directory holding Zea/NAN_Indexes/) or BWA_IDX}/Zea/NAN_Indexes/Index_Zm_B73v5_bwa}"

LOCALINSTALL="${LOCALINSTALL:-${HOME}/LocalInstall}"
WASP="${WASP_DIR:-$LOCALINSTALL/WASP}"

NCPU="${SLURM_CPUS_PER_TASK:-16}"
CHROMS_CSV="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10"
CHROMS_WS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10"

# Per-sample label scheme (mirrors 06_44_run.sh exactly so wasp/non-wasp counts
# are comparable): B73Mo17 -> bc_key_map (AM genome_1 top-1 call); multi ->
# plate-of-origin ground truth. Multi REUSES 06_44's already-built 7-geno sites
# VCF (06_44_concordance/multi_7geno_sites.vcf.gz) — do NOT rebuild it here.
case "$SAMPLE" in
  B73Mo17_rep1|B73Mo17_rep2)
    GROUP_MODE=bc_key_map
    GENOTYPES="B73,Mo17"
    EXTRACT_VCF=$GENO/configs/Final_Mo17_relative_to_B73.vcf.gz ;;
  multiGenotypes_rep1)
    GROUP_MODE=plate
    GENOTYPES="B73,B97,Ky21,M162W,Mo18W,Oh7B,Tzi8"
    WELL_MAP=$GENO/configs/Well_to_Genotype_multiGenotypes_rep1.txt
    EXTRACT_VCF=$GENO/${SAMPLE}/diagnostics/06_44_concordance/multi_7geno_sites.vcf.gz ;;
  *) echo "ERROR: $SAMPLE not wired (supports B73Mo17_rep1/rep2, multiGenotypes_rep1)"; exit 1 ;;
esac

OUT44=$GENO/${SAMPLE}/diagnostics/06_44_concordance   # reuse bc-map + panel
OUT=$GENO/${SAMPLE}/diagnostics/06_46_wasp
SNPDIR=$OUT/wasp_snps
mkdir -p "$OUT" "$SNPDIR"

# =============================================================== STEP: setup
if [[ "$STEP" == "setup" ]]; then
  echo "[06_46] SETUP for $SAMPLE"
  # 1. WASP repo (git clone needs internet -> run on a login node)
  if [[ ! -e "$WASP/mapping/find_intersecting_snps.py" ]]; then
    echo "[06_46] cloning WASP -> $WASP"
    git clone https://github.com/bmvdgeijn/WASP.git "$WASP"
  else
    echo "[06_46] WASP present: $WASP"
  fi
  $PY -c "import pysam,numpy; print('[06_46] pysam',pysam.__version__,'numpy',numpy.__version__)"
  # 2. per-chrom SNP files for WASP: {chrom}.snps.txt.gz = pos<TAB>ref<TAB>alt
  #    biallelic SNPs only (WASP requires 1bp alleles), chr1-10.
  echo "[06_46] building WASP snp_dir (biallelic SNPs, chr1-10) -> $SNPDIR"
  keep_re='^(chr1|chr2|chr3|chr4|chr5|chr6|chr7|chr8|chr9|chr10)$'
  # single streaming pass; split by chrom; then gzip.
  rm -f "$SNPDIR"/*.snps.txt "$SNPDIR"/*.snps.txt.gz
  bcftools view -H -v snps -m2 -M2 "$EXTRACT_VCF" \
    | awk -v d="$SNPDIR" -v re="$keep_re" 'BEGIN{OFS="\t"}
        $1 ~ re { print $2,$4,$5 >> (d"/"$1".snps.txt") }'
  for c in $CHROMS_WS; do
    [[ -s "$SNPDIR/$c.snps.txt" ]] && gzip -f "$SNPDIR/$c.snps.txt" \
      && echo "[06_46]   $c: $(zcat "$SNPDIR/$c.snps.txt.gz" | wc -l) snps"
  done
  echo "[06_46] SETUP done."; exit 0
fi

# =============================================================== STEP: wasp
if [[ "$STEP" == "wasp" ]]; then
  MODE="${MODE:?wasp needs MODE=raw|clean}"
  # bwa via module (0.7.17 that mapped these BAMs), then scrub the py3.9/PYTHONPATH
  # injection so $PY/$RSCRIPT below run on clean base py3.12 (see the env note above).
  module load Bioinformatics "$BWA_MODULE" 2>/dev/null || true
  command -v samtools >/dev/null 2>&1 || module load samtools 2>/dev/null || true
  unset PYTHONPATH PYTHONHOME
  BWA="${BWA:-bwa}"
  SAMTOOLS="${SAMTOOLS:-samtools}"
  case "$MODE" in
    raw)   BAM=$MAP/ambientmapper_input/${SAMPLE}_B73v5_scifiATAC.mq10.BC.rmdup.mm.bam ;;
    clean) BAM=$AM/${SAMPLE}/clean_bams_alpha05_C0_nd/${SAMPLE}_B73v5_scifiATAC.mq10.BC.rmdup.mm.Clean.bam ;;
    *) echo "MODE must be raw|clean"; exit 1 ;;
  esac
  [[ -f "$BAM"     ]] || { echo "Missing BAM: $BAM"; exit 1; }
  [[ -f "$BAM.bai" ]] || { echo "Missing BAM index (.bai): $BAM.bai — index in a scratch dir, do not modify the data dir"; exit 1; }
  [[ -e "$WASP/mapping/find_intersecting_snps.py" ]] || { echo "WASP not set up — run STEP=setup first"; exit 1; }
  [[ -s "$SNPDIR/chr1.snps.txt.gz" ]] || { echo "snp_dir empty — run STEP=setup first"; exit 1; }
  command -v "$BWA" >/dev/null       || { echo "bwa not found ($BWA) — module load Bioinformatics $BWA_MODULE failed; set BWA=/abs/path/bwa"; exit 1; }
  command -v "$SAMTOOLS" >/dev/null  || { echo "samtools not found ($SAMTOOLS) — set SAMTOOLS=/abs/path/samtools"; exit 1; }
  [[ -f "$BWA_IDX.bwt" ]]            || { echo "bwa index missing: $BWA_IDX.bwt"; exit 1; }
  echo "[06_46] bwa=$($BWA 2>&1 | awk '/^Version/{print $2; exit}')  ($(command -v "$BWA"))"
  echo "[06_46] samtools=$($SAMTOOLS --version 2>/dev/null | awk 'NR==1{print $2}')  py=$($PY --version 2>&1)"

  WORK=$OUT/work_$MODE
  mkdir -p "$WORK"
  base=$(basename "$BAM" .bam)
  echo "[06_46] WASP sample=$SAMPLE mode=$MODE"
  echo "[06_46]   BAM=$BAM"
  echo "[06_46]   WORK=$WORK  base=$base"

  # 1. find reads overlapping SNPs; write allele-swapped remap FASTQs
  echo "[06_46] (1) find_intersecting_snps"
  $PY "$WASP/mapping/find_intersecting_snps.py" \
    --is_paired_end --is_sorted \
    --output_dir "$WORK" --snp_dir "$SNPDIR" \
    "$BAM"
  # outputs: $WORK/${base}.keep.bam  ${base}.to.remap.bam  ${base}.remap.fq1.gz  ${base}.remap.fq2.gz

  # 2. remap the allele-swapped reads with the ORIGINAL aligner/params (bwa mem -M)
  echo "[06_46] (2) bwa mem remap"
  "$BWA" mem -M -t "$NCPU" "$BWA_IDX" \
      "$WORK/${base}.remap.fq1.gz" "$WORK/${base}.remap.fq2.gz" \
    | "$SAMTOOLS" sort -@ "$NCPU" -o "$WORK/${base}.remap.sorted.bam" -
  "$SAMTOOLS" index -@ "$NCPU" "$WORK/${base}.remap.sorted.bam"

  # 3. keep only reads that remapped back to the same locus
  echo "[06_46] (3) filter_remapped_reads"
  $PY "$WASP/mapping/filter_remapped_reads.py" \
    "$WORK/${base}.to.remap.bam" \
    "$WORK/${base}.remap.sorted.bam" \
    "$WORK/${base}.remap.keep.bam"

  # 4. merge the never-remapped keepers + the passed remaps -> WASP-passed BAM
  echo "[06_46] (4) merge + sort + index"
  "$SAMTOOLS" merge -f -@ "$NCPU" "$WORK/${base}.wasp.unsorted.bam" \
    "$WORK/${base}.keep.bam" "$WORK/${base}.remap.keep.bam"
  WASP_BAM=$OUT/${SAMPLE}_${MODE}.wasp.bam
  "$SAMTOOLS" sort -@ "$NCPU" -o "$WASP_BAM" "$WORK/${base}.wasp.unsorted.bam"
  "$SAMTOOLS" index -@ "$NCPU" "$WASP_BAM"

  # sanity: the kept reads must still carry the BC tag (06_44 groups by it)
  nbc=$("$SAMTOOLS" view "$WASP_BAM" 2>/dev/null | head -100000 | grep -c 'BC:Z:' || true)
  echo "[06_46]   BC-tagged reads in first 100k of wasp.bam: $nbc"
  [[ "$nbc" -gt 0 ]] || { echo "[06_46] ERROR: BC tag lost by WASP — 06_44 grouping would fail"; exit 1; }

  # read accounting
  echo "[06_46]   read counts:"
  printf '    %-14s %s\n' input       "$("$SAMTOOLS" view -c "$BAM")"
  printf '    %-14s %s\n' keep        "$("$SAMTOOLS" view -c "$WORK/${base}.keep.bam")"
  printf '    %-14s %s\n' to.remap    "$("$SAMTOOLS" view -c "$WORK/${base}.to.remap.bam")"
  printf '    %-14s %s\n' remap.kept  "$("$SAMTOOLS" view -c "$WORK/${base}.remap.keep.bam")"
  printf '    %-14s %s\n' wasp.pass   "$("$SAMTOOLS" view -c "$WASP_BAM")"

  # 5. re-count alleles with the SAME 06_44 extractor -> counts + .bias.tsv sidecar.
  #    Group labels reuse 06_44's exact per-sample scheme (bc_key_map for B73Mo17,
  #    plate-of-origin for multi) so wasp vs non-wasp numbers stay comparable.
  echo "[06_46] (5) 06_44 extractor on wasp.bam (group_mode=$GROUP_MODE)"
  if [[ "$GROUP_MODE" == bc_key_map ]]; then
    BC_MAP=$OUT44/bc_to_genome1.tsv
    [[ -s "$BC_MAP" ]] || { echo "Missing bc-map $BC_MAP — run 06_44 prep/extract first"; exit 1; }
    GROUP_ARGS=(--group-mode bc_key_map --bc-map "$BC_MAP")
  else
    [[ -s "$WELL_MAP" ]] || { echo "Missing well-map $WELL_MAP"; exit 1; }
    GROUP_ARGS=(--group-mode plate --well-map "$WELL_MAP")
  fi
  $PY "${REPO_ROOT}/workflows/03b_variant_based_comparison/06_44_extract_allele_counts.py" \
    --bam "$WASP_BAM" --vcf "$EXTRACT_VCF" \
    "${GROUP_ARGS[@]}" \
    --bc-tag BC --min-mapq 10 --threads "$NCPU" --chroms "$CHROMS_CSV" \
    --out "$OUT/${SAMPLE}_${MODE}_wasp_allele_counts.tsv.gz"

  # free the big intermediates (keep the bams for debugging; drop the fastqs)
  rm -f "$WORK/${base}.remap.fq1.gz" "$WORK/${base}.remap.fq2.gz" "$WORK/${base}.wasp.unsorted.bam"
  echo "[06_46] WASP done mode=$MODE -> $WASP_BAM  + counts"; exit 0
fi

# =============================================================== STEP: analyze
if [[ "$STEP" == "analyze" ]]; then
  RAW=$OUT/${SAMPLE}_raw_wasp_allele_counts.tsv.gz
  CLN=$OUT/${SAMPLE}_clean_wasp_allele_counts.tsv.gz
  for f in "$RAW" "$CLN"; do [[ -s "$f" ]] || { echo "Missing wasp counts (run STEP=wasp first): $f"; exit 1; }; done

  PANEL=$OUT/genotype_panel.tsv
  if [[ ! -s "$PANEL" ]]; then
    if [[ -s "$OUT44/genotype_panel.tsv" ]]; then
      cp "$OUT44/genotype_panel.tsv" "$PANEL"       # same VCF -> identical panel
    else
      bcftools query -s "$GENOTYPES" -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$EXTRACT_VCF" > "$PANEL"
    fi
  fi

  echo "[06_46] ANALYZE (concordance on WASP-corrected counts) sample=$SAMPLE"
  $RSCRIPT "${REPO_ROOT}/workflows/03b_variant_based_comparison/06_44_concordance.R" \
    --counts-raw "$RAW" --counts-clean "$CLN" \
    --genotype-panel "$PANEL" --genotypes "$GENOTYPES" \
    --out-prefix "$OUT/${SAMPLE}_wasp" \
    --min-reads 2 --n-boot 200 --block-mb 1

  # ---- side-by-side: WASP vs non-WASP (refbias + self-concordance) -----------
  echo
  echo "=================  WASP vs non-WASP compare ($SAMPLE)  ================="
  $PY - "$SAMPLE" "$OUT" "$OUT44" <<'PY'
import sys, os, csv
sample, out, out44 = sys.argv[1], sys.argv[2], sys.argv[3]

def read_bias(path):
    d={}
    if not os.path.exists(path): return d
    with open(path) as fh:
        r=csv.DictReader(fh, delimiter='\t')
        for row in r:
            n=int(row['n_obs']); nm=int(row['sum_nm'])
            d[(row['group'],row['obs_class'])]=(n, nm/n if n else float('nan'))
    return d

print("\n-- ref-mapping bias: mean NM per allele class (want ALT->REF gap to shrink) --")
print(f"{'mode':6} {'group':6} {'class':4} {'n_obs':>12} {'mean_nm':>8}")
for tag, path in [("noWASP", f"{out44}/{sample}_raw_allele_counts.tsv.bias.tsv"),
                  ("WASP",   f"{out}/{sample}_raw_wasp_allele_counts.tsv.bias.tsv")]:
    b=read_bias(path)
    for (g,c),(n,mnm) in sorted(b.items()):
        print(f"{tag:6} {g:6} {c:4} {n:12,d} {mnm:8.3f}")

def read_selfconc(path):
    # *_self_concordance_ci.tsv from 06_44_concordance.R
    rows=[]
    if not os.path.exists(path): return rows
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            rows.append(row)
    return rows

print("\n-- self-concordance (WASP-corrected) -- see full CIs in the tsv --")
for row in read_selfconc(f"{out}/{sample}_wasp_self_concordance_ci.tsv"):
    print("  ", {k:row[k] for k in row})
print(f"\n(compare against non-WASP: {out44}/{sample}_self_concordance_ci.tsv)")
PY
  echo "[06_46] ANALYZE done -> $OUT/${SAMPLE}_wasp_*"; exit 0
fi

echo "Unknown STEP: $STEP (need setup|wasp|analyze)"; exit 1
