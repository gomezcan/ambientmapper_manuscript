#!/bin/bash
#SBATCH --job-name=root1_sub_filter
#SBATCH --time=04:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --output=_logs/02_08_filter_%A_%a.out
#SBATCH --error=_logs/02_08_filter_%A_%a.err
#SBATCH --partition=standard
#SBATCH --array=0-51

#
# 02_08_filter_bams_subsample.sh — Filter Root1_rep1 BAMs to a balanced subsample
#
# Reads each Root1_rep1_<library>.mq10.BC.rmdup.mm.bam and writes a new BAM
# containing only reads whose BC:Z: tag is in the panel barcode list.
#
# Array layout: 26 genomes × 2 panels (A, B) = 52 tasks
#   panel  = panels[task_id / 26]
#   genome = genomes[task_id % 26]
#
# Inputs (must exist):
#   Root1_rep1/sub1k_<panel>/barcodes_by_genome.tsv
#       (produced by 02_07_subsample_root1_balanced.py)
#
# Outputs:
#   Root1_rep1/sub1k_<panel>/bams/Root1_rep1_<library>.mq10.BC.rmdup.mm.bam[.bai]

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"   # root of the data tree (3_Mapping/, 5_AmbientDetection/, ...)

cd "${PROJECT_ROOT}/5_AmbientDetection"

# conda activate <env from environment.yml>

OMP_NUM_THREADS=1
MKL_NUM_THREADS=1

SAMPLE=Root1_rep1
SRC_BAM_DIR=${PROJECT_ROOT}/3_Mapping/${SAMPLE}

# Genome → library suffix (B73 uses B73v5 in BAM filenames)
GENOMES=(B73 B97 CML103 CML228 CML247 CML277 CML322 CML333 CML52 CML69 \
         HP301 Il14H Ki11 Ki3 Ky21 M162W M37W Mo18W Ms71 NC350 NC358 \
         Oh43 Oh7B P39 Tx303 Tzi8)
LIBRARIES=(B73v5 B97 CML103 CML228 CML247 CML277 CML322 CML333 CML52 CML69 \
           HP301 Il14H Ki11 Ki3 Ky21 M162W M37W Mo18W Ms71 NC350 NC358 \
           Oh43 Oh7B P39 Tx303 Tzi8)
PANELS=(A B)
N_GENOMES=${#GENOMES[@]}

PANEL_IDX=$(( SLURM_ARRAY_TASK_ID / N_GENOMES ))
GENOME_IDX=$(( SLURM_ARRAY_TASK_ID % N_GENOMES ))
PANEL=${PANELS[$PANEL_IDX]}
GENOME=${GENOMES[$GENOME_IDX]}
LIBRARY=${LIBRARIES[$GENOME_IDX]}

PANEL_DIR=${SAMPLE}/sub1k_${PANEL}
BC_FILE=${PANEL_DIR}/barcodes_by_genome.tsv
SRC_BAM=${SRC_BAM_DIR}/${SAMPLE}_${LIBRARY}.mq10.BC.rmdup.mm.bam
OUT_BAM_DIR=${PANEL_DIR}/bams
OUT_BAM=${OUT_BAM_DIR}/${SAMPLE}_${LIBRARY}.mq10.BC.rmdup.mm.bam

mkdir -p "${OUT_BAM_DIR}"

if [[ ! -s ${BC_FILE} ]]; then
    echo "ERROR: missing ${BC_FILE}" >&2
    exit 1
fi
if [[ ! -s ${SRC_BAM} ]]; then
    echo "ERROR: missing ${SRC_BAM}" >&2
    exit 1
fi

echo "Panel  = ${PANEL}"
echo "Genome = ${GENOME} (library=${LIBRARY})"
echo "Source = ${SRC_BAM}"
echo "Output = ${OUT_BAM}"

python - <<PYEOF
import pysam, sys

bc_file  = "${BC_FILE}"
genome   = "${GENOME}"
src_bam  = "${SRC_BAM}"
out_bam  = "${OUT_BAM}"

# Load BC strings for this genome.
#
# barcodes_by_genome.tsv contains the full SAM-text form with the "BC:Z:"
# type prefix (e.g. "BC:Z:TGTTGACGATGGGCCT-Root1_rep1_B73v5"), but
# pysam's aln.get_tag("BC") returns only the VALUE of the tag — the
# "BC:Z:" prefix is stripped when pysam parses the SAM spec. So we strip
# the prefix here before building the keep set so the membership check
# matches.
BC_PREFIX = "BC:Z:"
keep = set()
with open(bc_file) as f:
    header = f.readline().rstrip("\n").split("\t")
    g_i  = header.index("genome")
    bc_i = header.index("bc_string")
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if parts[g_i] == genome:
            bc = parts[bc_i]
            if bc.startswith(BC_PREFIX):
                bc = bc[len(BC_PREFIX):]
            keep.add(bc)
print(f"  loaded {len(keep)} BC values for {genome}", flush=True)
if not keep:
    sys.exit(f"ERROR: no barcodes found for {genome} in {bc_file}")

bam_in  = pysam.AlignmentFile(src_bam, "rb")
bam_out = pysam.AlignmentFile(out_bam, "wb", template=bam_in)

n_in = n_out = n_bc_tag = 0
for aln in bam_in.fetch(until_eof=True):
    n_in += 1
    if not aln.has_tag("BC"):
        continue
    n_bc_tag += 1
    if aln.get_tag("BC") in keep:
        bam_out.write(aln)
        n_out += 1

bam_in.close()
bam_out.close()
print(f"  wrote {n_out:,} / {n_bc_tag:,} BC-tagged "
      f"/ {n_in:,} total alignments", flush=True)
if n_out == 0:
    sys.exit(f"ERROR: no alignments matched any panel barcode for {genome}. "
             f"Check BC tag format in source BAM.")
PYEOF

samtools index "${OUT_BAM}"
echo "  indexed ${OUT_BAM}.bai"
echo "Done."
