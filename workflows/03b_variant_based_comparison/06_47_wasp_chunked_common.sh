# 06_47_wasp_chunked_common.sh — sourced config / reuse seam for the chunked WASP family.
# =============================================================================
# NOT executable on its own. Sourced by 06_47_wasp_chunked_{prep,array,combine,analyze}.sh.
#
# WHY chunked. The monolithic 06_46_wasp_pilot.sh runs WASP genome-wide and OOM-killed
# on multiGenotypes_rep1 at stage 3 filter_remapped_reads.py (--mem=72G vs a 254 GB
# to.remap.bam — the filter holds a per-read-name dict that scales with reads-to-remap).
# Chunking WASP by chromosome is *provably behavior-preserving* (WASP flushes its
# read-pair cache at every chromosome boundary → it already processes one chrom at a
# time; remap still uses the full genome-wide bwa index), and it bounds the filter's
# peak RAM to ~the largest chrom (chr1 ≈ 16%), ~10× lower.
#
# This file provides (for the 4 step scripts that source it):
#   * env: conda base ($PY/$RSCRIPT), bwa module name, bwa index, WASP clone path
#   * per-sample case: GROUP_MODE / GENOTYPES / EXTRACT_VCF / WELL_MAP  (== 06_46)
#   * SNPDIR — REUSED from 06_46 (already-built per-chrom biallelic SNPs; no rebuild)
#   * chrom list (chr1..10) + a 1-indexed array for SLURM_ARRAY_TASK_ID
#   * path helpers (slice / per-chrom keep / per-chrom remap.keep)
#   * set_group_args() — the 06_44 extractor's per-sample group flags (combine only)
#
# Requires SAMPLE in env. MODE (raw|clean) is required by prep/array/combine and is
# left unset by analyze (which reads both modes' counts).
#
# ENV NOTE (same as 06_46): `ml Bioinformatics` breaks base py3.12
# two ways (shadows `python`, injects a py3.9 numpy via PYTHONPATH). So: only the array
# step loads the bwa module, and every step that runs python/Rscript does so via the
# ABSOLUTE base interpreter ($PY/$RSCRIPT) after `unset PYTHONPATH PYTHONHOME`. This
# file itself only conda-activates base (never loads a module).
# =============================================================================

: "${SAMPLE:?need SAMPLE (B73Mo17_rep1|B73Mo17_rep2|multiGenotypes_rep1)}"

PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"   # root of the data tree (3_Mapping/, 5_AmbientDetection/, ...)
PROJ=${PROJECT_ROOT}
GENO=$PROJ/5_Genotyping
AM=$PROJ/5_AmbientDetection
MAP=$PROJ/3_Mapping
cd "$GENO"

# --- base env: pysam/numpy + bcftools/bgzip/tabix + samtools + R (SNPRelate,...) ---
# conda activate <env from environment.yml>
PY="$CONDA_PREFIX/bin/python";        [[ -x "$PY" ]]      || PY=python
RSCRIPT="$CONDA_PREFIX/bin/Rscript";  [[ -x "$RSCRIPT" ]] || RSCRIPT=Rscript
BWA_MODULE="${BWA_MODULE:-bwa/0.7.17-mil4ns7}"   # the 0.7.17 that mapped these BAMs
BWA_IDX="${BWA_IDX:-${BWA_INDEX_ROOT:?set BWA_INDEX_ROOT (directory holding Zea/NAN_Indexes/) or BWA_IDX}/Zea/NAN_Indexes/Index_Zm_B73v5_bwa}"
LOCALINSTALL="${LOCALINSTALL:-${HOME}/LocalInstall}"
WASP="${WASP_DIR:-$LOCALINSTALL/WASP}"

NCPU="${SLURM_CPUS_PER_TASK:-16}"
CHROMS_CSV="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10"
CHROMS_WS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10"
CHROMS_ARR=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10)  # 1-indexed: ${CHROMS_ARR[$((id-1))]}

# --- per-sample label scheme (mirrors 06_46 / 06_44_run.sh exactly so wasp/non-wasp
#     counts stay comparable). B73Mo17 -> bc_key_map (AM genome_1 top-1 call);
#     multi -> plate-of-origin. Multi REUSES 06_44's already-built 7-geno sites VCF. ---
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
  *) echo "ERROR: $SAMPLE not wired (B73Mo17_rep1/rep2, multiGenotypes_rep1)"; exit 1 ;;
esac

OUT44=$GENO/${SAMPLE}/diagnostics/06_44_concordance          # reuse bc-map + panel
OUT=$GENO/${SAMPLE}/diagnostics/06_47_wasp_chunked           # NEW output tree
SNPDIR=$GENO/${SAMPLE}/diagnostics/06_46_wasp/wasp_snps      # REUSE 06_46's per-chrom SNPs
mkdir -p "$OUT"

# --- BAM resolution (raw|clean); only used when MODE is set (prep/array/combine) ---
resolve_bam() {  # $1 = raw|clean -> echoes BAM path
  case "$1" in
    raw)   echo "$MAP/ambientmapper_input/${SAMPLE}_B73v5_scifiATAC.mq10.BC.rmdup.mm.bam" ;;
    clean) echo "$AM/${SAMPLE}/clean_bams_alpha05_C0_nd/${SAMPLE}_B73v5_scifiATAC.mq10.BC.rmdup.mm.Clean.bam" ;;
    *) echo "MODE must be raw|clean" >&2; return 1 ;;
  esac
}

if [[ -n "${MODE:-}" ]]; then
  BAM="$(resolve_bam "$MODE")"
  BASE="$(basename "$BAM" .bam)"
  WORK=$OUT/work_$MODE
  SLICES=$WORK/slices
fi

# --- per-chrom path helpers. Per-chrom outputs are CANONICALLY named keep.bam /
#     remap.keep.bam inside a per-chrom dir, independent of WASP's basename scheme. ---
slice_path()  { echo "$SLICES/${BASE}.$1.bam"; }       # $1 = chrom
chrom_dir()   { echo "$WORK/$1"; }                     # $1 = chrom
chrom_keep()  { echo "$WORK/$1/keep.bam"; }            # WASP keep (renamed to canonical)
chrom_rkeep() { echo "$WORK/$1/remap.keep.bam"; }      # filter_remapped_reads output

# --- 06_44 extractor group flags (combine only). Validates inputs, sets GROUP_ARGS. ---
set_group_args() {
  if [[ "$GROUP_MODE" == bc_key_map ]]; then
    local bc_map=$OUT44/bc_to_genome1.tsv
    [[ -s "$bc_map" ]] || { echo "Missing bc-map $bc_map — run 06_44 prep/extract first"; exit 1; }
    GROUP_ARGS=(--group-mode bc_key_map --bc-map "$bc_map")
  else
    [[ -s "${WELL_MAP:-}" ]] || { echo "Missing well-map ${WELL_MAP:-<unset>}"; exit 1; }
    GROUP_ARGS=(--group-mode plate --well-map "$WELL_MAP")
  fi
}
