#!/usr/bin/env bash
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --job-name=gene_body_matrix
#SBATCH --partition=standard
#SBATCH --output=_logs/%x-%j.log
#
# Cluster-annotation phase, Task 1: build the gene-body accessibility ("gene activity") matrix the
# annotation steps consume -- gene x cell raw Tn5 counts, written as <prefix>.genes.sparse.rds
# (mirrors the maize_282 reference all_pools.raw.genes.sparse.rds).
#
# Array: 0 = PRE  (raw   SM2_*        BEDs -> SM2/step5_markers/SM2.genes.sparse.rds)
#        1 = POST (clean Clean.SM2v2_* BEDs -> SM2v2_clean/step5_markers/Clean.SM2v2.genes.sparse.rds)
#   sbatch --array=0-1 0_scripts/common/4_0_build_gene_body_matrix.sh
#
# Gene window = gene body + UP bp upstream (5', strand-aware) + DOWN bp downstream. UP/DOWN are
# tunable via env vars; first pass uses gene body + 2 kb promoter. Maize (sparse) tolerates a wide
# promoter; At (gene-dense) prefers a tighter one -- revisit per-species if At markers look smeared.
#   UP=500 DOWN=0 sbatch --array=0-1 0_scripts/common/4_0_build_gene_body_matrix.sh

set -euo pipefail
PROJECT_ROOT="${PROJECT_ROOT:-.}"
BASE="${PROJECT_ROOT}/6_socrates"
cd "$BASE"
# conda activate ambientmapper-manuscript   (environment.yml at the repo root)

GTF="_data/_GenomeInfo/ZmATcombined.gtf"
BEDDIR="_data/_BED_files"
UP="${UP:-2000}"; DOWN="${DOWN:-0}"
mkdir -p _logs

case "${SLURM_ARRAY_TASK_ID:-0}" in
  0) OUT="SM2/step5_markers/SM2"
     B_AT="${BEDDIR}/SM2_At_ZmATcombined_scifiATAC.mq10.BC.rmdup.mm.tn5.bed.gz"
     B_B73="${BEDDIR}/SM2_B73_ZmATcombined_scifiATAC.mq10.BC.rmdup.mm.tn5.bed.gz" ;;
  1) OUT="SM2v2_clean/step5_markers/Clean.SM2v2"
     B_AT="${BEDDIR}/Clean.SM2v2_At_ZmATcombined_scifiATAC.mq10.BC.rmdup.mm.tn5.bed.gz"
     B_B73="${BEDDIR}/Clean.SM2v2_B73_ZmATcombined_scifiATAC.mq10.BC.rmdup.mm.tn5.bed.gz" ;;
  *) echo "ERROR: array index must be 0 (PRE) or 1 (POST), got ${SLURM_ARRAY_TASK_ID:-unset}"; exit 1 ;;
esac

for f in "$GTF" "$B_AT" "$B_B73"; do
  [[ -s "$f" ]] || { echo "ERROR: missing input: $f"; exit 1; }
done

echo " - building gene-body matrix: OUT=${OUT}  window=+${UP}/${DOWN} bp"
Rscript "${SCRIPTS:-${BASE}/0_scripts}/common/4_0_build_gene_body_matrix.R" "$OUT" "$GTF" "$UP" "$DOWN" "$B_AT" "$B_B73"
echo " - done -> ${OUT}.genes.sparse.rds"
