#!/bin/bash
########## BATCH Lines for Resource Request ##########
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --job-name=Step4_markers_SM2v2_plate
#SBATCH --partition=standard
#SBATCH --output=_logs/Step4_markers_SM2v2_plate_%A_%a.log
#SBATCH --array=0-2
#
# 4_3_eval_plate.sh  --  Part-2 Step-4 (markers), PLATE arm.
# ---------------------------------------------------------------------------
# Plate sibling of 4_3_eval_pergenome.sh (indep). Per genome x stage: build the per-genome
# gene-body accessibility matrix (common/4_0) from the PLATE-SPLIT BED, then run the marker
# annotation (common/4_3) on the FROZEN plate Step-2 v7 clusters.
#
# WHY: the plate At Pre 5-cluster solution has cluster 1 = low depth / high
# organelle+doublet, and nd cleaning strips 90% of it (43x odds) while wd keeps it. Markers
# adjudicate whether cl1 is ambient (no coherent At identity -> nd was RIGHT / wd underclean)
# or real At (clear At cell-type call -> nd OVERCLEANED). Yardstick from indep At Pre:
# coherent = top_score>~0.5 & low entropy (e.g. cl2 At:stele 2.167); ambient = top_score~0.1,
# entropy~3. Plate At should also read dominant_species=At (vs indep's B73-contaminated object).
#
# ANNOTATION PANEL = informative_top15 (the no-cap augmented panel is for INFORMATIVENESS
# only; the mean-z call uses markers.{At,maize}.informative_top15.bed + the 4_3 min_markers=3 guard).
#
# Array (At first = the decisive cl1 question; B73 after):
#   0: SM2_At_TAIR10             Pre   TAIR10        3: SM2_B73_B73v5             Pre   B73v5
#   1: Clean.SM2v2wd_At_TAIR10   wd    TAIR10        4: Clean.SM2v2wd_B73_B73v5   wd    B73v5
#   2: Clean.SM2v2_At_TAIR10     nd    TAIR10        5: Clean.SM2v2_B73_B73v5     nd    B73v5
#   sbatch --array=0-2 ...   # At (decisive) -- default
#   sbatch --array=3-5 ...   # B73
#
# Prereq: plate Step-2 clustering done (SM2v2_plate/step2_cluster/<name>.mQCv6.updated_metadata_v7.*.txt);
#         plate-split BEDs present; informative panels + per-species GFF present.
# NOTE: check <adir>/<name>.marker_overlap_qc.tsv on first run -- ~0 overlap = GFF geneIDs != panel geneIDs.
# ---------------------------------------------------------------------------

# conda activate ambientmapper-manuscript   (environment.yml at the repo root)
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-.}"
BASE="${PROJECT_ROOT}/6_socrates"
SCRIPTS="${SCRIPTS:-${BASE}/0_scripts}"
BEDDIR="${BASE}/_data/_BED_files"
GI="${BASE}/_data/_GenomeInfo"
MARK="${BASE}/_data/markers"
STEP2="${BASE}/SM2v2_plate/step2_cluster"
OUT="${BASE}/SM2v2_plate/step3_compare"
mkdir -p "$OUT"

s4_0="${SCRIPTS}/common/4_0_build_gene_body_matrix.R"
s4_3="${SCRIPTS}/common/4_3_marker_accessibility.R"
UP=500; DOWN=500                       # gene body + 500/500 promoter window (Part-1 locked value)
NORM="${GENE_LENGTH_NORM:-raw}"        # raw (default) | perkb -> length-normalized gene activity (writes .perkb variant)
case "$NORM" in raw) VAR="" ;; perkb) VAR=".perkb" ;; *) echo "ERROR: GENE_LENGTH_NORM must be raw|perkb"; exit 1 ;; esac

declare -a NAMES=( SM2_At_TAIR10 Clean.SM2v2wd_At_TAIR10 Clean.SM2v2_At_TAIR10 SM2_B73_B73v5 Clean.SM2v2wd_B73_B73v5 Clean.SM2v2_B73_B73v5 )
declare -a STAGES=( Pre wd nd Pre wd nd )

i="${SLURM_ARRAY_TASK_ID:?set via --array}"
name="${NAMES[$i]:?bad array index $i}"
stage="${STAGES[$i]}"
species="${name##*_}"                  # trailing token -> B73v5 or TAIR10

case "$species" in
  B73v5)  ann="${GI}/B73v5.gff3";  panel="${MARK}/markers.maize.informative_top15.bed" ;;
  TAIR10) ann="${GI}/TAIR10.gff3"; panel="${MARK}/markers.At.informative_top15.bed" ;;
  *) echo "ERROR: unknown species '$species' from name '$name'"; exit 1 ;;
esac

bed="${BEDDIR}/${name}_scifiATAC.mq10.BC.rmdup.mm.tn5.bed.gz"
meta="$(ls ${STEP2}/${name}.mQCv6.updated_metadata_v7.*.txt 2>/dev/null | head -1)"

[[ -f "$bed"   ]]                  || { echo "ERROR: missing BED $bed"; exit 2; }
[[ -f "$ann"   ]]                  || { echo "ERROR: missing annotation $ann"; exit 2; }
[[ -n "${meta:-}" && -f "$meta" ]] || { echo "ERROR: no plate v7 metadata for $name in $STEP2"; exit 2; }
[[ -f "$panel" ]]                  || { echo "ERROR: missing panel $panel"; exit 2; }

mat="${OUT}/${name}.plate${VAR}.genes.sparse.rds"
adir="${OUT}/4_3_annotation_plate_${name}${VAR}"

echo "$(date): Step-4 plate markers | $name (species=$species stage=$stage norm=$NORM)"
echo "  bed=$(basename "$bed")  ann=$(basename "$ann")  meta=$(basename "$meta")  panel=$(basename "$panel")"

# (1) per-genome gene-body matrix from the plate BED (skip if present and newer than the BED)
if [[ -f "$mat" && "$mat" -nt "$bed" ]]; then
  echo " - [1/2] gene-body matrix up to date ($(basename "$mat")) -- skip"
else
  echo " - [1/2] building gene-body matrix ($NORM) -> $(basename "$mat")"
  GENE_LENGTH_NORM="$NORM" Rscript "$s4_0" "${OUT}/${name}.plate${VAR}" "$ann" "$UP" "$DOWN" "$bed"
fi

# (2) marker annotation on the frozen plate v7 clusters (informative panel + min_markers=3 default)
echo " - [2/2] marker annotation -> $(basename "$adir")"
Rscript "$s4_3" "$adir" "$name" "$mat" "$meta" "$panel" "$stage" LouvainClusters

echo "$(date): DONE $name -> $adir/${name}.cluster_annotation.tsv"
