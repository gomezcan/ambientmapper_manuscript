#!/bin/bash
########## BATCH Lines for Resource Request ##########
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --job-name=Step4_bgnull_SM2v2_plate
#SBATCH --partition=standard
#SBATCH --output=_logs/Step4_bgnull_SM2v2_plate_%A_%a.log
#SBATCH --array=0-5
#
# 4_3n_bgnull_plate.sh -- Step-4 (plate arm) TEST A: covariate-conditioned background-gene null
#                         for the per-cell reciprocal z produced by 4_3i.
# ---------------------------------------------------------------------------
# Answers "is a cluster's marker-type call beyond chance?" by pushing EVERY gene in the smoothed
# perkb matrix through the identical rZ pipeline, then comparing the panel against background genes
# at the SAME expression. rZ is expression-driven (via Zj) and zero-inflated, so each gene is first
# converted to its mid-rank quantile among background genes in its own expression neighbourhood.
# Two levels:
#   SET  (cluster x type)  = headline, competitive gene-set test, BH over ~n_cl x n_type
#   GENE (per marker)      = detail layer (p floor ~ 1/(2*window))
# Two conditioning modes in ONE run (they share the expensive genome-wide pass):
#   mean     = expression conditioned                    (primary)
#   mean_sd  = expression AND variability conditioned    (conservative)
#
# Validated on a synthetic fixture with known answers (6 independent replicates): planted signal
# recovered 6/6; structure-free decoys at high expression 0/6 false positives; sets lying outside
# the background expression range are marked UNTESTABLE rather than silently called.
#
#   0: SM2_At_TAIR10            Pre    3: SM2_B73_B73v5            Pre
#   1: Clean.SM2v2wd_At_TAIR10  wd     4: Clean.SM2v2wd_B73_B73v5  wd
#   2: Clean.SM2v2_At_TAIR10    nd     5: Clean.SM2v2_B73_B73v5    nd
#
# Anchor for the significance thread = At, geom, Pre + wd  ->  sbatch --array=0,1 <this>
# Full sweep (adds nd + maize, quantifies the nd over-cleaning loss) -> sbatch --array=0-5 <this>
#
# Prereq: SMOOTHED perkb matrix (4_3f_smooth_plate.sh GENE_LENGTH_NORM=perkb) + Step-2 v7 metadata
#         + the 4_3i output dir (used only for the reconciliation check).
# Env knobs: GENE_LENGTH_NORM(perkb) METRIC(geom) WINDOW(500) NDRAW(10000) FDR(0.05) MODE(mean,mean_sd)
# ---------------------------------------------------------------------------

# conda activate ambientmapper-manuscript   (environment.yml at the repo root)
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-.}"
BASE="${PROJECT_ROOT}/6_socrates"
SCRIPTS="${SCRIPTS:-${BASE}/0_scripts}"
MARK="${BASE}/_data/markers"
STEP2="${BASE}/SM2v2_plate/step2_cluster"
OUT="${BASE}/SM2v2_plate/step3_compare"
s4_3n="${SCRIPTS}/common/4_3n_background_null_rZ.R"

NORM="${GENE_LENGTH_NORM:-perkb}"      # rZ REQUIRES perkb (Zj compares genes within a cell)
METRIC="${METRIC:-geom}"               # geom = AND = primary metric for the rZ framework
WINDOW="${WINDOW:-500}"                # background genes per expression neighbourhood
NDRAW="${NDRAW:-10000}"                # set-level null draws (p floor = 1/(1+NDRAW))
FDR="${FDR:-0.05}"
MODE="${MODE:-mean,mean_sd}"
EXCLUDE="${EXCLUDE:-dividing}"         # cell-cycle: cross-tissue confound, excluded from annotation
MINMK="${MINMK:-3}"                    # min markers for a type to be set-testable
TOPK="${TOPK:-5}"                      # top-k sensitivity statistic (mean of the k highest-u markers)

case "$NORM" in raw) VAR="" ;; perkb) VAR=".perkb" ;; *) echo "ERROR: GENE_LENGTH_NORM must be raw|perkb"; exit 1 ;; esac

declare -a NAMES=( SM2_At_TAIR10 Clean.SM2v2wd_At_TAIR10 Clean.SM2v2_At_TAIR10 SM2_B73_B73v5 Clean.SM2v2wd_B73_B73v5 Clean.SM2v2_B73_B73v5 )
declare -a STAGES=( Pre wd nd Pre wd nd )
i="${SLURM_ARRAY_TASK_ID:?set via --array}"
name="${NAMES[$i]:?bad array index $i}"; stage="${STAGES[$i]}"
species="${name##*_}"

# PANEL: which marker panel the test runs on. The convention (Step-4 note) is that the UNCAPPED
# augmented panel is for INFORMATIVENESS ONLY and annotation uses the
# informative-pruned/capped panel. That matters doubly here:
#   - set test : an uncapped set (e.g. At phloem, 676 markers) is diluted by markers irrelevant to this
#                tissue, so its mean regresses to background -> loss of power. Capping at 15 keeps sets
#                comparable (3-15) and on-target.
#   - gene test: BH over 278 capped markers instead of 4587 is ~16x better multiplicity, so the
#                gene layer stops being floored out.
# Outputs are suffixed per panel so the two runs never collide.
PANEL="${PANEL:-top15}"                # top15 (annotation convention, DEFAULT) | augmented (informativeness)
case "$PANEL" in
  top15)     PSUF=".top15" ;;
  augmented) PSUF=".aug"   ;;
  *) echo "ERROR: PANEL must be top15|augmented"; exit 1 ;;
esac
case "${species}:${PANEL}" in
  B73v5:top15)      panel="${MARK}/markers.maize.informative_top15.bed" ;;
  B73v5:augmented)  panel="${MARK}/markers.maize.Marand2025_PlantscRNAdb4.bed" ;;
  TAIR10:top15)     panel="${MARK}/markers.At.informative_top15.bed" ;;
  TAIR10:augmented) panel="${MARK}/markers.At.Ecker2025_PlantscRNAdb4.bed" ;;
  *) echo "ERROR: unknown species '$species' from '$name'"; exit 1 ;;
esac

sm="${OUT}/${name}.plate${VAR}.genes.smoothed.rds"
meta="$(ls ${STEP2}/${name}.mQCv6.updated_metadata_v7.*.txt 2>/dev/null | head -1)"
cmp="${OUT}/4_3i_reciprocal_percell_${name}${VAR}/${name}.percell_rZ.cluster_mean.tsv"
adir="${OUT}/4_3n_bgnull_${name}${VAR}${PSUF}"

[[ -f "$sm" ]]                     || { echo "ERROR: missing smoothed matrix $sm (run 4_3f_smooth_plate.sh GENE_LENGTH_NORM=$NORM first)"; exit 2; }
[[ -n "${meta:-}" && -f "$meta" ]] || { echo "ERROR: no v7 metadata for $name in $STEP2"; exit 2; }
[[ -f "$panel" ]]                  || { echo "ERROR: missing panel $panel"; exit 2; }
[[ -f "$cmp" ]]                    || { echo "NOTE: no 4_3i cluster_mean for $name -- reconciliation check will be skipped"; cmp=""; }

echo "$(date): Step-4 TEST A bgnull | $name ($stage norm=$NORM metric=$METRIC mode=$MODE window=$WINDOW B=$NDRAW panel=$PANEL:$(basename "$panel"))"
Rscript "$s4_3n" "$adir" "$name" "$sm" "$meta" "$panel" "$stage" \
        LouvainClusters "$WINDOW" "$MODE" "$EXCLUDE" "$METRIC" "$FDR" "$NDRAW" "$MINMK" 1 "$cmp" "$TOPK"
echo "$(date): DONE $name -> $adir"
