#!/usr/bin/env bash
# 0_04_split_reads_to_drop_SM2v2.sh
#
# Split SM2v2's per-read drop list (from AmbientMapper decontam, computed on the
# independent-genome mapping) into two target-specific drop lists, so we can
# apply them to the *combined-genome* ZmATcombined BAMs (one BAM per species
# subset).
#
# Same column convention as legacy SM2 (read_id, barcode, bc_key, allowed_set,
# winner_genome, p_as, reason):
#   - cleanAt list  = rows where allowed_set=At AND winner_genome=B73  (contam in At BAM)
#   - cleanB73 list = rows where allowed_set=B73 AND winner_genome=At  (contam in B73 BAM)
#
# Outputs land in 6_socrates/_data/ (alongside the legacy SM2 lists), so they're
# co-located with the BAM symlinks the clean-bams step reads.
#
# One-shot helper (no SBATCH needed; awk over a gzipped TSV takes minutes max).
# Run: PROJECT_ROOT=<data root> bash <repo>/workflows/04_decontamination/combined_genome/0_04_split_reads_to_drop_SM2v2.sh

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"   # root of the data tree (5_AmbientDetection/, 6_socrates/, ...)

BASE="${PROJECT_ROOT}/6_socrates"
DATA="${BASE}/_data"

IN="${PROJECT_ROOT}/5_AmbientDetection/SM2v2/decontam_with_design_alpha05_v2/SM2v2_reads_to_drop.tsv.gz"
OUT_AT="${DATA}/SM2v2_reads_to_drop.cleanAt.tsv.gz"
OUT_B73="${DATA}/SM2v2_reads_to_drop.cleanB73.tsv.gz"

[[ -f "${IN}" ]] || { echo "ERROR: missing input ${IN}"; exit 1; }

echo "[split] reads_to_drop -> per-target drop lists"
echo "  in  : ${IN}"
echo "  outAt : ${OUT_AT}"
echo "  outB73: ${OUT_B73}"

# Reads that landed in At wells but got winner=B73 -> drop from At BAM
zcat "${IN}" \
  | awk -F'\t' 'NR==1 || ($4=="At" && $5=="B73")' \
  | gzip -c > "${OUT_AT}"

# Reads that landed in B73 wells but got winner=At -> drop from B73 BAM
zcat "${IN}" \
  | awk -F'\t' 'NR==1 || ($4=="B73" && $5=="At")' \
  | gzip -c > "${OUT_B73}"

# Sanity counts (rows excluding header)
for f in "${OUT_AT}" "${OUT_B73}"; do
  printf "%s\t%d rows\n" "${f##*/}" "$(zcat "$f" | awk 'END{print NR-1}')"
done
