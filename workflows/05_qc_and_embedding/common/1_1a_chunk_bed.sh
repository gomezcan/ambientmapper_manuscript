#!/usr/bin/env bash
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=8G
#SBATCH --job-name=QCs_chunkBED
#SBATCH --partition=standard
#SBATCH --output=_logs/%x-%j.log
#
# 1_1a_chunk_bed.sh
#
# Hash-partition a tn5 BED into N chunks by barcode (col 4). Same barcode always
# hashes to the same chunk -> chunks have disjoint barcode sets, a hard
# requirement for mergeSocratesRDS in 1_1c. Ported from
# 2_PopulationStress_maize/05_Cells_QC/0_scripts/pre_ambient/1_1a_chunk_bed.sh.
#
# Inputs (env vars):
#   POOL     sample name (e.g. SM2_At, Clean.SM2v2_B73)
#   OUTDIR   output directory under 6_socrates/ where chunks/<POOL>/ is written
#   N        number of chunks (default 5)
#
# Input BED resolved by convention:
#   _data/_BED_files/${POOL}_ZmATcombined_scifiATAC.mq10.BC.rmdup.mm.tn5.bed.gz
#
# Outputs:
#   <OUTDIR>/chunks/<POOL>/<POOL>.chunk{0..N-1}.bed.gz
#   <OUTDIR>/chunks/<POOL>/.split.done    (sentinel; re-run skips if present)

set -euo pipefail

# Paths are anchored on PROJECT_ROOT ($BASH_SOURCE self-location does not survive sbatch).
PROJECT_ROOT="${PROJECT_ROOT:-.}"
BASE="${PROJECT_ROOT}/6_socrates"
cd "${BASE}"

POOL="${POOL:?POOL not set; submit via submit_chunked_pipeline_SM2v2.sh or --export=ALL,POOL=...}"
OUTDIR="${OUTDIR:?OUTDIR not set; e.g. SM2 (PreClean) or SM2v2_clean (PostClean)}"
N="${N:-5}"

BED="${BASE}/_data/_BED_files/${POOL}_ZmATcombined_scifiATAC.mq10.BC.rmdup.mm.tn5.bed.gz"
CHUNKDIR="${BASE}/${OUTDIR}/chunks/${POOL}"
SENTINEL="${CHUNKDIR}/.split.done"

if [ ! -e "$BED" ]; then
  echo "ERROR: BED not found at $BED" >&2
  exit 1
fi

# Sentinel resume.
if [ -f "$SENTINEL" ]; then
  echo " - SKIP split: sentinel found at $SENTINEL"
  ls -lh "$CHUNKDIR"/${POOL}.chunk*.bed.gz 2>/dev/null | head
  exit 0
fi

mkdir -p "$CHUNKDIR"
echo " - splitting $BED"
echo "   N chunks   : $N"
echo "   output dir : $CHUNKDIR"
echo "   $(du -h "$BED" | cut -f1) compressed"

# Clean any partial leftovers from a prior aborted split.
rm -f "$CHUNKDIR"/${POOL}.chunk*.bed.gz

# Hash-partition stream. Polynomial hash mod 1000003 over barcode chars, then
# mod N for chunk index. Explicit close() of every gzip pipe in END so the
# subprocesses flush + exit cleanly (without this awk can hang on CPU-starved
# nodes and produce truncated gzip outputs).
gunzip -c "$BED" \
  | awk -v N="$N" -v outdir="$CHUNKDIR" -v pool="$POOL" '
    BEGIN {
      for (i = 0; i < 256; i++) ord[sprintf("%c", i)] = i
    }
    {
      bc = $4
      h = 0
      L = length(bc)
      for (i = 1; i <= L; i++) h = (h * 31 + ord[substr(bc, i, 1)]) % 1000003
      chunk = h % N
      cmd = "gzip -c > " outdir "/" pool ".chunk" chunk ".bed.gz"
      cmds[chunk] = cmd
      print | cmd
      lines[chunk]++
      total++
      if (total % 50000000 == 0) {
        printf("   ... processed %12d rows\n", total) > "/dev/stderr"
        fflush("/dev/stderr")
      }
    }
    END {
      for (i = 0; i < N; i++) {
        if (i in cmds) {
          rc = close(cmds[i])
          if (rc != 0) {
            printf("   WARN: close(chunk %02d) returned %d\n", i, rc) > "/dev/stderr"
          }
        }
      }
      total = 0
      for (i = 0; i < N; i++) {
        if (!(i in lines)) lines[i] = 0
        printf("   chunk %02d : %12d lines\n", i, lines[i]) > "/dev/stderr"
        total += lines[i]
      }
      printf("   TOTAL    : %12d lines\n", total) > "/dev/stderr"
    }
  '

echo " - chunk files:"
ls -lh "$CHUNKDIR"/${POOL}.chunk*.bed.gz

touch "$SENTINEL"
echo " - split done. sentinel: $SENTINEL"
