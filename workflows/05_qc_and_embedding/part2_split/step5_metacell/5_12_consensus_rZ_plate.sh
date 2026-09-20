#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --job-name=Step5_12_consensus_rZ
#SBATCH --partition=standard
#SBATCH --output=_logs/Step5_12_consensus_rZ_%A_%a.log
#SBATCH --array=0-5
#
# 5_12_consensus_rZ_plate.sh -- marker rZ on the CONSENSUS meta-cells (the plot-ready tables for
# the Fig-5 heatmap script). SLURM OPTIONAL -- light enough to run locally (sparse perkb is
# 10-70 MB):   bash 5_12_consensus_rZ_plate.sh            # all six objects
#              bash 5_12_consensus_rZ_plate.sh 0 3        # just At Pre + maize Pre
# -------------------------------------------------------------------------------------------------
# WHAT IT DOES, per object (ENGINES REUSED, never forked -- only this driver is new):
#   1. adapter: consensus/<name>.consensus.cells.tsv -> cell_to_seacell.tsv format
#      (index, SEACell="cMC-<consensus_mc>", LouvainClusters), kept in _adapters/ for provenance
#   2. 0_scripts/common/4_3p_metacell_rZ.R on that membership (raw perkb, informative_top15 panel)
#   3. 0_scripts/common/4_3z_consensus_rZ_join.R -> <name>.consensus_rZ.heatmap_input.tsv
#      (rZ long format + conf_meanF/stab/loyalty/purity attributes, n_cells reconciliation gate)
#
# WHY CONSENSUS MEMBERSHIP: a per-seed marker plot resurrects the "which seed?" indefensibility
# the consensus layer was built to close (At Pre seeds agree at ARI 0.47). The consensus partition
# with its confidence attribute is the unit downstream work can name.
#
# READING RULES: descriptive until the meta-cell GENE null runs (top_type is provisional, not
#   annotation); geom is the primary metric (Zi cap = (N-1)/sqrt(N) = 3.47 for At's k=14 -- 4_3p
#   prints it); WEIGHT by confidence, never filter; condition on n_cells cross-stage (At nd -34%).
#
# The Fig 5 consensus panel (analysis/fig5_biological_impact/fig5_I_consensus.R) consumes
# SM2v2_plate/step5_metacell/rZ_annotation/consensus/*.heatmap_input.tsv.
# -------------------------------------------------------------------------------------------------

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-.}"
BASE="${PROJECT_ROOT}/6_socrates"
S3="${BASE}/SM2v2_plate/step3_compare"
S5="${BASE}/SM2v2_plate/step5_metacell"
MK="${BASE}/_data/markers"
CONS="${S5}/consensus"
OUT="${S5}/rZ_annotation/consensus"
ADP="${OUT}/_adapters"
ENGINE_RZ="${SCRIPTS:-${BASE}/0_scripts}/common/4_3p_metacell_rZ.R"
ENGINE_JOIN="${SCRIPTS:-${BASE}/0_scripts}/common/4_3z_consensus_rZ_join.R"
TOPN="${TOPN:-6}"

declare -a NAMES=( SM2_At_TAIR10 Clean.SM2v2wd_At_TAIR10 Clean.SM2v2_At_TAIR10 \
                   SM2_B73_B73v5 Clean.SM2v2wd_B73_B73v5 Clean.SM2v2_B73_B73v5 )
declare -a STAGES=( Pre wd nd Pre wd nd )
declare -a PANELS=( At At At maize maize maize )

mkdir -p "$OUT" "$ADP"

run_one () {
  local i="$1"
  local name="${NAMES[$i]}" stage="${STAGES[$i]}" panel="${PANELS[$i]}"
  local cells="${CONS}/${name}.consensus.cells.tsv"
  local c2s="${ADP}/${name}.consensus_c2s.tsv"
  local mat="${S3}/${name}.plate.perkb.genes.sparse.rds"
  local bed="${MK}/markers.${panel}.informative_top15.bed"
  [[ -s "$cells" ]] || { echo "ERROR: missing $cells -- run 5_9 first"; return 3; }
  [[ -s "$mat"   ]] || { echo "ERROR: missing $mat"; return 3; }
  [[ -s "$bed"   ]] || { echo "ERROR: missing $bed"; return 3; }

  echo "--- [${i}] ${name} (${stage}) : adapter -> 4_3p -> 4_3z"
  Rscript - "$cells" "$c2s" <<'RS'
a <- commandArgs(trailingOnly = TRUE)
d <- read.delim(a[1], check.names = FALSE)
stopifnot(all(c("cellID", "consensus_mc", "LouvainClusters") %in% names(d)))
write.table(data.frame(index = d$cellID, SEACell = paste0("cMC-", d$consensus_mc),
                       LouvainClusters = d$LouvainClusters),
            a[2], sep = "\t", quote = FALSE, row.names = FALSE)
RS
  Rscript "$ENGINE_RZ" "$OUT" "$name" "$mat" "$c2s" "$bed" "$stage" "$TOPN"
  Rscript "$ENGINE_JOIN" "$OUT" "$CONS" "$name" "$stage" "$panel"
}

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  run_one "$SLURM_ARRAY_TASK_ID"
else
  tasks=( "$@" ); [[ ${#tasks[@]} -gt 0 ]] || tasks=( 0 1 2 3 4 5 )
  for t in "${tasks[@]}"; do run_one "$t"; done
fi
echo "$(date): 5_12 done"
