#!/usr/bin/env bash
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --job-name=QCs_1_2_1_3
#SBATCH --partition=standard
#SBATCH --output=_logs/%x-%j.log
#
# 1_2_1_3_run.sh — final stage of the chunked SM2v2 Socrates chain.
#
# Runs:
#   1_2_filter_lowQC_cells_scifiATAC_data.R  <POOL>.raw.soc.rds <POOL>
#       -> <POOL>.updated_metadata.txt
#   1_3_metaQC_scifiATAC_data.R              <POOL>.updated_metadata.txt <POOL> <DEPTH>
#       -> <POOL>.minDepth<DEPTH>.updated_metadata_v{1..6}.txt + plots/<POOL>.minDepth<DEPTH>.QC_FIGURES.pdf
#
# Both steps are skipped if their primary output is already present (sentinel
# resume — matches the reference 1_QC_scifiATAC_data.sh pattern).
#
# Submitted by submit_chunked_pipeline_SM2v2.sh with --export=ALL,POOL=...,OUTDIR=...

set -euo pipefail

# Paths are anchored on PROJECT_ROOT ($BASH_SOURCE self-location does not survive sbatch).
PROJECT_ROOT="${PROJECT_ROOT:-.}"
BASE="${PROJECT_ROOT}/6_socrates"
SCRIPTS="${SCRIPTS:-${BASE}/0_scripts}"

POOL="${POOL:?POOL not set; submit via submit_chunked_pipeline_SM2v2.sh}"
OUTDIR="${OUTDIR:?OUTDIR not set}"
DEPTH="${DEPTH:-200}"   # 200 matches the established SM2 analysis (Fig 1 objects); was 500

# conda activate ambientmapper-manuscript   (environment.yml at the repo root)
# module load macs2   (site-specific; Socrates' .preRunChecks() needs macs2 on PATH even with precomputed peaks)

step_1_2="${SCRIPTS}/common/1_2_filter_lowQC_cells_scifiATAC_data.R"
step_1_3="${SCRIPTS}/common/1_3_metaQC_scifiATAC_data.R"

mkdir -p "${BASE}/${OUTDIR}/step0_qc/plots"
cd "${BASE}/${OUTDIR}/step0_qc"

rds="${POOL}.raw.soc.rds"
meta="${POOL}.updated_metadata.txt"

[[ -s "$rds" ]] || { echo "ERROR: ${rds} missing (run 1_1c_merge_and_qc.sh first)"; exit 1; }

# ---- 1_2: filter low-QC cells + good/bad reference correlation -------------
if [ -s "$meta" ]; then
  echo " - SKIP 1_2: $meta already present ($(du -h "$meta" | cut -f1))"
else
  echo " - running 1_2 for $POOL"
  Rscript "$step_1_2" "$rds" "$POOL"
fi

# ---- 1_3: meta-QC depth cascade + 5-panel PDF ------------------------------
v6="${POOL}.minDepth${DEPTH}.updated_metadata_v6.txt"
if [ -s "$v6" ]; then
  echo " - SKIP 1_3: $v6 already present"
else
  echo " - running 1_3 for $POOL (depth_filter=${DEPTH})"
  Rscript "$step_1_3" "$meta" "$POOL" "$DEPTH"
fi

echo " - 1_2 + 1_3 done for $POOL (OUTDIR=$OUTDIR)"
ls -lh "${POOL}".*updated_metadata*.txt 2>/dev/null || true
