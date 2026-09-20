#!/usr/bin/env bash
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=180G
#SBATCH --job-name=QCs_mergeSOC
#SBATCH --partition=standard
#SBATCH --output=_logs/%x-%j.log
#
# 1_1c_merge_and_qc.sh — merge chunk soc.objs + run isCellv2.
#
# After this completes, downstream 1_2 + 1_3 can run against the merged
# <OUTDIR>/<POOL>.raw.soc.rds (drop-in for the monolithic 1_1 output).
#
# Memory: 180G is the safety floor for the in-place sparseMatrix() triplet
# rebuild inside mergeSocratesRDS — the temporary as.data.frame(summary(counts))
# expansion is the peak.

set -euo pipefail

# Paths are anchored on PROJECT_ROOT ($BASH_SOURCE self-location does not survive sbatch).
PROJECT_ROOT="${PROJECT_ROOT:-.}"
BASE="${PROJECT_ROOT}/6_socrates"
SCRIPTS="${SCRIPTS:-${BASE}/0_scripts}"
cd "${BASE}"

POOL="${POOL:?POOL not set; submit via submit_chunked_pipeline_SM2v2.sh}"
OUTDIR="${OUTDIR:?OUTDIR not set}"
N="${N:-5}"

# conda activate ambientmapper-manuscript   (environment.yml at the repo root)
# module load macs2   (site-specific; Socrates' .preRunChecks() needs macs2 on PATH even with precomputed peaks)

CHUNKDIR="${BASE}/${OUTDIR}/chunks/${POOL}"
RESULTSDIR="${BASE}/${OUTDIR}/step0_qc"

mkdir -p "$RESULTSDIR"
cd "$RESULTSDIR"

OUT="${POOL}"
RDS="${OUT}.raw.soc.rds"

if [ -s "$RDS" ]; then
  echo " - SKIP merge: $RDS already present ($(du -h "$RDS" | cut -f1))"
  exit 0
fi

# Verify all N chunk soc.objs exist.
missing=()
for i in $(seq 0 $((N - 1))); do
  f="${CHUNKDIR}/${POOL}.chunk${i}.soc.rds"
  [ -s "$f" ] || missing+=("$f")
done
if [ "${#missing[@]}" -gt 0 ]; then
  echo "ERROR: ${#missing[@]} chunk soc.obj files missing:" >&2
  printf '  %s\n' "${missing[@]}" >&2
  exit 1
fi

echo " - merging $N chunks for $POOL via mergeSocratesRDS + isCellv2"
echo "   chunkdir : $CHUNKDIR"
echo "   output   : $RESULTSDIR/$RDS"
Rscript "$SCRIPTS/common/1_1c_merge_and_qc.R" "$POOL" "$CHUNKDIR" "$OUT" "$N"

echo " - merge done."
ls -lh "$RDS"

# Post-success cleanup. Set KEEP_CHUNKS=1 to retain intermediates for debug.
if [ ! -s "$RDS" ]; then
  echo "ERROR: merged $RDS is missing or empty - leaving chunks/ in place" >&2
  exit 1
fi
if [ "${KEEP_CHUNKS:-0}" = "1" ]; then
  echo " - KEEP_CHUNKS=1 set; leaving $CHUNKDIR in place ($(du -sh "$CHUNKDIR" 2>/dev/null | cut -f1))"
else
  echo " - cleaning per-chunk intermediates in $CHUNKDIR"
  du -sh "$CHUNKDIR" 2>/dev/null | sed 's/^/   before: /'
  rm -f "$CHUNKDIR"/${POOL}.chunk*.bed.gz \
        "$CHUNKDIR"/${POOL}.chunk*.raw.before.soc.rds \
        "$CHUNKDIR"/${POOL}.chunk*.soc.rds \
        "$CHUNKDIR"/.split.done
  rmdir "$CHUNKDIR" 2>/dev/null || echo "   note: $CHUNKDIR not empty after cleanup; left in place"
  echo "   done."
fi
