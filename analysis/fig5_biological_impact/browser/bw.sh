#!/bin/bash
########## BATCH Lines for Resource Request ##########
#SBATCH --time=0:30:00            # Wall clock limit
#SBATCH --nodes=1                 # Number of nodes
#SBATCH --ntasks=1                # Number of tasks
#SBATCH --cpus-per-task=5         # CPUs per task
#SBATCH --mem=5G                  # Memory per node
#SBATCH --job-name=fig5_L_bed_to_bw
#SBATCH --partition=standard
#SBATCH --output=figures/main/fig5/browser/_logs_bw/%x-%A-%a.log
#SBATCH --array=0-119
# =============================================================================
# browser/bw.sh  -  STEP 3 of the Fig 5 panel L genome-browser chain (HPC only, SLURM array):
#   per-group BED.gz -> RPM-normalised BigWig, one array task per track.
# Inputs  figures/main/fig5/browser/bed_by_type/*.bed.gz          (browser/split.sh)
#         $DATA/socrates/_data/_GenomeInfo/SM2v2_{B73v5,TAIR10}.chrs.size.txt
# Output  figures/main/fig5/browser/BWs/<track>.bw (+ _bed/, _bedgraph/ intermediates)
# Requires bedtools, wigToBigWig (UCSC), pigz in the active conda env.
# Usage   (from the repo root, on the HPC)
#   mkdir -p figures/main/fig5/browser/_logs_bw
#   sbatch analysis/fig5_biological_impact/browser/bw.sh
# =============================================================================
#
# Design notes:
#   (1) chrom sizes are resolved PER FILE from the track name (TAIR10 vs B73v5): two genomes;
#   (2) no `grep -v '^scaf'`: the split step already whitelisted chroms against the genome's
#       size file, so genomecov can never see an unknown contig;
#   (3) the array over-allocates (0-119 vs ~100 real tracks) and surplus tasks exit 0 with a
#       message instead of exit 1 -- an over-range task is expected here, not an error.
# Coordinate sort, RPM scale = 1e6 / total reads, wigToBigWig.
# SLOP: tn5 sites are extended by +/- SLOP bp. At +/-49 (a common single-bp convention) the
#   tracks looked SPARSE at browser scale, so the default here is +/-100 (~200 bp footprint per
#   insertion). Override per run with:  SLOP=49 sbatch --export=ALL,SLOP=49 browser/bw.sh
#   Changing SLOP changes every track: resubmit the WHOLE array and re-render.
# RPM is PER TRACK: each group is scaled by its own read total. Fine within a track; when
#   comparing a locus ACROSS stages remember the groups are different cell sets picked by
#   different partitions (argmax calls).
# =============================================================================

# conda activate <env from environment.yml>   # must provide bedtools, wigToBigWig, pigz

REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
cd "$REPO_ROOT"

DATA=data/processed/scifiATAC_B73_Arabidopsis
BROWSER=figures/main/fig5/browser
INPUTDIR="$BROWSER/bed_by_type"
OUTDIR="$BROWSER/BWs"
GINFO="$DATA/socrates/_data/_GenomeInfo"

export NTHREADS=${SLURM_CPUS_PER_TASK:-4}
mkdir -p "$OUTDIR" "$OUTDIR/_bed" "$OUTDIR/_bedgraph"

# build and check track list
mapfile -t SAMPLES < <(
  ls "$INPUTDIR"/*.bed.gz | xargs -n1 basename -s .bed.gz | sort
)
IDX=${SLURM_ARRAY_TASK_ID:-0}
if (( IDX >= ${#SAMPLES[@]} )); then
  echo "array task $IDX >= ${#SAMPLES[@]} tracks - nothing to do"
  exit 0
fi
SAMPLE=${SAMPLES[$IDX]}

# per-file genome -> chrom sizes
case "$SAMPLE" in
  *At_TAIR10*) CHR_SIZES="$GINFO/SM2v2_TAIR10.chrs.size.txt" ;;
  *B73_B73v5*) CHR_SIZES="$GINFO/SM2v2_B73v5.chrs.size.txt"  ;;
  *) echo "cannot resolve genome from track name: $SAMPLE"; exit 1 ;;
esac

echo "$(date): Processing $SAMPLE  (sizes: $(basename "$CHR_SIZES"))"

BED="$INPUTDIR/${SAMPLE}.bed.gz"
SORTED="$OUTDIR/_bed/${SAMPLE}.sorted.bed.gz"
BG="$OUTDIR/_bedgraph/${SAMPLE}.bedgraph"
BW="$OUTDIR/${SAMPLE}.bw"

# sanity checks
[[ -f "$BED"       ]] || { echo "Missing input BED: $BED"; exit 1; }
[[ -f "$CHR_SIZES" ]] || { echo "Missing chrom.sizes: $CHR_SIZES"; exit 1; }

# 1) extend & sort
SLOP=${SLOP:-100}
echo "$(date): Step 1: slop +/-${SLOP}bp & sort"
pigz -dc -p $NTHREADS "$BED" \
  | bedtools slop -i - -g "$CHR_SIZES" -b "$SLOP" \
  | sort --parallel=$NTHREADS -k1,1 -k2,2n \
  | pigz -p $NTHREADS > "$SORTED"

# 2) compute scale & generate RPM-normalized BedGraph
echo "$(date): Step 2: compute scale factor and bedGraph"
TOTAL_READS=$(pigz -dc -p $NTHREADS "$SORTED" | wc -l)
[[ "$TOTAL_READS" -gt 0 ]] || { echo "EMPTY track: $SAMPLE"; exit 1; }

SCALE=$(awk -v tot=$TOTAL_READS 'BEGIN{printf "%.8f", 1000000/tot}')
echo "  -> Total reads = $TOTAL_READS; scale = $SCALE"
pigz -dc -p $NTHREADS "$SORTED" \
  | bedtools genomecov -i stdin -g "$CHR_SIZES" -bg -scale $SCALE \
  > "$BG"

# 3) convert to BigWig
echo "$(date): Step 3: bigWig conversion"
wigToBigWig "$BG" "$CHR_SIZES" "$BW"
[[ -s "$BW" ]] || { echo "FAILED to write $BW"; exit 1; }

echo "$(date): Done $SAMPLE"
