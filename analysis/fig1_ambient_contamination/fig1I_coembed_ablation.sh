#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G
#SBATCH --job-name=fig1I_coembed_ablation
#SBATCH --partition=standard
#SBATCH --output=_logs/fig1I_coembed_ablation_%A.log
#
# fig1I_coembed_ablation.sh -- Figure 1, panel I: oracle feature ablation of the Fig 1F co-embedding.
# Runs analysis/fig1_ambient_contamination/fig1I_coembed_ablation.R on the raw per-genome Socrates
# objects with the frozen Fig 1F configuration. Takes NO arguments; every value is defaulted and
# overridable by environment variable for a documented sensitivity run.
#   mkdir -p _logs && sbatch analysis/fig1_ambient_contamination/fig1I_coembed_ablation.sh   # the control
#   MIN_BINS=20 sbatch analysis/fig1_ambient_contamination/fig1I_coembed_ablation.sh         # sensitivity run
#   GRID=1 sbatch analysis/fig1_ambient_contamination/fig1I_coembed_ablation.sh              # + the 36-config band (slower)
#   bash analysis/fig1_ambient_contamination/fig1I_coembed_ablation.sh                       # same run, foreground
# Run from the repository root. HPC only (R with Socrates, Seurat, FNN); the memory peak is the load
# of the raw B73v5 object (about 2.05M bins x 1.96M barcodes), which is subsetted and freed at once.
#
# What it does. Gives each Fig 1F barcode ONLY the bins of the genome its PLATE INDEX says it came
# from, re-embeds with the IDENTICAL frozen configuration, and re-scores the same k = 15 obs/exp
# species-mixing statistic. If the species separate, the co-projection was carried by cross-genome
# signal, not by shared cell identity. AmbientMapper is nowhere in this pipeline (raw pre-QC
# per-genome matrices in, plate index as the only mask), which is what makes the control non-circular.
#
# Inputs are the RAW step0_qc objects, not the co-embed matrix: the co-embed build zero-fills the
# genome a barcode did not pass QC on, so ablating inside that matrix would zero most
# Arabidopsis-plate barcodes for a reason that is an artifact of the build.
#
# One cell gate only, and it is reported: MIN_BINS non-zero own-genome bins (default 50). No QC is
# re-derived, and callClusters(m.clst) is used for LABELS ONLY; the mixing statistic is computed on
# the full embedded set. The ladder file reports the per-plate attrition.
#
# THE FROZEN CONFIGURATION. These values are DEFAULTED, never required, and must not be changed: the
# control is valid only while the ablated object is embedded EXACTLY as Fig 1F was, so that the
# only difference between the panels is the ablation (held fixed for cross-object comparability,
# not re-optimised). They are the values of COEMB_CFG in fig1.R:
#     pcs=20  k_near=30  min_dist=0.3  res=0.5  m.clst=50
# Requiring the caller to pass them would invite a typo that silently invalidates the control.
#
# Two different k values: KNN=30 is the UMAP neighbourhood that BUILDS the embedding; MIXK=15 is the
# mixing-statistic neighbourhood (the number quoted in the manuscript).
#
# Overridable (all optional): STAGE PCS KNN MD MIN_BINS RES MCLST MIXK OWNMIN GRID SEED REFMETA REFMINC.
# STAGE defaults to Pre (Fig 1F is the pre-clean panel).
#
# Requires the co-embed barcode set (.coembed.meta_full.tsv, built by Stage 3.1b inside
# workflows/05_qc_and_embedding/part1_indep/3_3_coembed_cluster.sh) and the step0_qc raw objects.
# This script never writes into the input directories.

set -euo pipefail
REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
cd "$REPO_ROOT"
mkdir -p _logs

# conda activate <env from environment.yml>

DATA=data/processed/scifiATAC_B73_Arabidopsis
SOC="${DATA}/socrates/SM2v2_indep"
OUT="${SOC}/coembed"
QC="${SOC}/step0_qc"
SCR=analysis/fig1_ambient_contamination
ABL=figures/main/fig1/coembed_ablation     # outputs (object, metadata, ladder, mixing summary, plots/)

STAGE="${STAGE:-Pre}"
# --- frozen Fig 1F configuration. DEFAULTED, never required, DO NOT CHANGE (COEMB_CFG in fig1.R). ---
PCS="${PCS:-20}"
KNN="${KNN:-30}"
MD="${MD:-0.3}"
RES="${RES:-0.5}"
MCLST="${MCLST:-50}"      # LABELS ONLY -- never gates the mixing statistic
# --- the one cell gate, and the readout knobs ---
MIN_BINS="${MIN_BINS:-50}"  # embeddability: non-zero own-genome bins. Reported per plate.
MIXK="${MIXK:-15}"          # mixing-statistic k
OWNMIN="${OWNMIN:-200}"     # own-genome READ floor for the second reported scope
GRID="${GRID:-0}"           # 1 -> also run the 36-configuration robustness grid
SEED="${SEED:-1}"

case "$STAGE" in
  Pre)     RAW_B="${QC}/SM2_B73v5.raw.soc.rds"
           RAW_A="${QC}/SM2_TAIR10.raw.soc.rds" ;;
  Post_wd) RAW_B="${QC}/Clean.SM2v2wd_B73v5.raw.soc.rds"
           RAW_A="${QC}/Clean.SM2v2wd_TAIR10.raw.soc.rds" ;;
  *) echo "unknown STAGE '$STAGE' (expected Pre or Post_wd)"; exit 2 ;;
esac

PREFIX="${OUT}/SM2v2_coembed_${STAGE}"
META_FULL="${PREFIX}.coembed.meta_full.tsv"

if [[ ! -s "$META_FULL" ]]; then
  echo "MISSING the Fig 1F barcode set: $META_FULL"
  echo "  It is built by Stage 3.1b, which runs inside 3_3_coembed_cluster.sh:"
  echo "  STAGE=${STAGE} PCS=${PCS} KNN=${KNN} MD=${MD} RES=${RES} sbatch workflows/05_qc_and_embedding/part1_indep/3_3_coembed_cluster.sh"
  echo "  (this script deliberately does not write into the coembed/ directory)"
  exit 2
fi
for f in "$RAW_B" "$RAW_A"; do
  [[ -s "$f" ]] || { echo "MISSING raw per-genome object: $f"; exit 2; }
done

# The unablated comparator: the Fig 1F v7 metadata at the frozen config. Its filename carries the
# clustering step's `minc_50` tag, which is that step's cell gate and is unrelated to MIN_BINS here;
# do not substitute one for the other.
REFMINC="${REFMINC:-50}"
CFG_TAG="pcs_${PCS}.k_near_${KNN}.min_dis_${MD}.minc_${REFMINC}.res_${RES}"
REFMETA="${REFMETA:-${OUT}/${STAGE}_cluster/SM2v2_coembed_${STAGE}.updated_metadata_v7.${CFG_TAG}.txt}"
if [[ ! -s "$REFMETA" ]]; then
  echo "$(date): NOTE -- unablated comparator not found, continuing without it:"
  echo "    $REFMETA"
  REFMETA="NA"
fi

mkdir -p "${ABL}/plots"

echo "$(date): [$STAGE] ORACLE FEATURE ABLATION, rebuilt from the RAW per-genome objects"
echo "  barcode set (Fig 1F)    : ${META_FULL}"
echo "  raw own-genome sources  : ${RAW_B}"
echo "                            ${RAW_A}"
echo "  frozen embedding config : pcs=${PCS} k_near=${KNN} min_dist=${MD} res=${RES} m.clst=${MCLST} (labels only)"
echo "  embeddability gate      : >= ${MIN_BINS} non-zero own-genome bins  (the ONLY cell filter, reported per plate)"
echo "  mixing statistic        : k=${MIXK} (obs/exp vs 2p(1-p)), on the EMBEDDED set; second scope own reads >= ${OWNMIN}"
echo "  robustness grid         : ${GRID}  (1 = run the 36-configuration band)"
echo "  reference (unablated)   : ${REFMETA}"
echo "  outputs                 : ${ABL}"

Rscript "${SCR}/fig1I_coembed_ablation.R" \
        "$META_FULL" "$RAW_B" "$RAW_A" \
        "$ABL" "SM2v2_coembed_${STAGE}" \
        "$PCS" "$KNN" "$MD" "$MIN_BINS" "$RES" \
        "$SEED" "$MCLST" "$MIXK" "$OWNMIN" "$GRID" "$REFMETA"

echo "$(date): DONE [$STAGE]"
echo " -> ${ABL}/SM2v2_coembed_${STAGE}.ablated.mixing_summary.*.tsv        (the headline obs/exp numbers)"
echo " -> ${ABL}/SM2v2_coembed_${STAGE}.ablated.embeddability_ladder.*.tsv  (per-plate ladder -- report it with the panel)"
echo " -> ${ABL}/SM2v2_coembed_${STAGE}.ablated.updated_metadata_v7.*.txt   (umap1/umap2 + Genome + own reads/bins + mixing)"
echo "    report obs/exp, never raw, and always beside the ladder."
