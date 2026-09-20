#!/usr/bin/env bash
#SBATCH --job-name=SM2v2_clean_bams_combined
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=40G
#SBATCH --time=03:00:00
#SBATCH --array=0-1
#SBATCH --output=_logs/0_05_clean_bams_SM2v2_combined_%A_%a.log
#
# 0_05_clean_bams_SM2v2_combined.sh
#
# Apply SM2v2's reads_to_drop (computed from independent-genome decontam) to the
# *combined-genome* ZmATcombined BAMs in 3_Mapping/SM2_{At,B73}/, producing
# Clean.SM2v2_*.ZmATcombined_*.bam alongside the legacy Clean.SM2_* outputs.
#
# This keeps the joint (concatenated-reference) coordinate space of the
# combined-genome Socrates objects (Fig. 5E, F; Fig. S2, S3) while applying the
# current SM2v2 cleaning. Run from ${PROJECT_ROOT}/6_socrates (logs go to _logs/).
#
# Prerequisites:
#   - 0_04_split_reads_to_drop_SM2v2.sh has run, producing:
#       _data/SM2v2_reads_to_drop.cleanAt.tsv.gz
#       _data/SM2v2_reads_to_drop.cleanB73.tsv.gz
#
# Array task layout:
#   0: At  (input  3_Mapping/SM2_At/SM2_At_ZmATcombined_*.bam,
#           drop   _data/SM2v2_reads_to_drop.cleanAt.tsv.gz)
#   1: B73 (input  3_Mapping/SM2_B73/SM2_B73_ZmATcombined_*.bam,
#           drop   _data/SM2v2_reads_to_drop.cleanB73.tsv.gz)
#
# Output naming intentionally mirrors legacy:  Clean.<basename>.bam  — but the
# *basename* is rewritten "SM2" -> "SM2v2" so the new Clean BAMs sit next to the
# legacy Clean.SM2_* files without overwriting them.

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"   # root of the data tree (3_Mapping/, 5_AmbientDetection/, ...)

# sbatch copies the script body to a spool directory at runtime, so $BASH_SOURCE
# self-location does not work; paths are anchored on PROJECT_ROOT instead.
BASE="${PROJECT_ROOT}/6_socrates"
cd "${BASE}"

# conda activate <env from environment.yml>

DATA="${BASE}/_data"
OUTDIR="${DATA}/SM2v2_clean_bams_combined"
mkdir -p "${OUTDIR}" "${BASE}/_logs"

idx="${SLURM_ARRAY_TASK_ID:?Array index required (sbatch --array=0-1)}"

species_arr=(At B73)
species="${species_arr[$idx]}"

# Inputs (combined-genome BAMs from the original SM2 mapping run)
RAW_BAM="${PROJECT_ROOT}/3_Mapping/SM2_${species}/SM2_${species}_ZmATcombined_scifiATAC.mq10.BC.rmdup.mm.bam"
DROP="${DATA}/SM2v2_reads_to_drop.clean${species}.tsv.gz"

# Output basename: rewrite SM2 -> SM2v2 so legacy and SM2v2 outputs can coexist
in_base="$(basename "${RAW_BAM}" .bam)"            # SM2_At_ZmATcombined_...
out_base="${in_base/SM2_/SM2v2_}"                  # SM2v2_At_ZmATcombined_...
CLEAN_BAM="${OUTDIR}/Clean.${out_base}.bam"

[[ -f "${RAW_BAM}"     ]] || { echo "ERROR: missing input BAM ${RAW_BAM}"; exit 1; }
[[ -f "${RAW_BAM}.bai" ]] || { echo "ERROR: missing index ${RAW_BAM}.bai"; exit 1; }
[[ -s "${DROP}"        ]] || { echo "ERROR: missing/empty drop list ${DROP}  (run 0_04 first)"; exit 1; }

echo "============================================================"
echo "[$(date)] SM2v2 clean-bams task ${idx} (species=${species})"
echo "  raw BAM   : ${RAW_BAM} ($(du -h "${RAW_BAM}" | cut -f1))"
echo "  drop list : ${DROP}    ($(zcat "${DROP}" | awk 'END{print NR-1}') rows)"
echo "  out BAM   : ${CLEAN_BAM}"
echo "============================================================"

# ambientmapper clean-bams emits <basename>.Clean.bam in --out-dir; rename after
ambientmapper clean-bams \
  --reads-to-drop "${DROP}" \
  --bam "${RAW_BAM}" \
  --out-dir "${OUTDIR}" \
  --out-suffix .Clean.bam \
  --index

# Rename: <in_base>.Clean.bam  ->  Clean.<out_base>.bam
mv "${OUTDIR}/${in_base}.Clean.bam"     "${CLEAN_BAM}"
mv "${OUTDIR}/${in_base}.Clean.bam.bai" "${CLEAN_BAM}.bai"

[[ -s "${CLEAN_BAM}"     ]] || { echo "ERROR: clean BAM missing/empty"; exit 1; }
[[ -s "${CLEAN_BAM}.bai" ]] || { echo "ERROR: clean BAM index missing/empty"; exit 1; }

echo "[$(date)] DONE species=${species}"
ls -la "${OUTDIR}"/Clean.${out_base}.*
