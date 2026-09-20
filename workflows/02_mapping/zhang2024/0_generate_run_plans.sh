#!/usr/bin/env bash
# 0_generate_run_plans.sh
# Generate the run_plan TSVs for scifi-demux step 2: one row = one (well, genome) pair = one SLURM array task.
# Columns match the scifi-demux genome_map: sample_base \t target_genome \t ref_path (BWA index prefix).
# Pools (well lists) come from config/Pools_scifi_*.txt, genome lists from this directory.
#
# Usage:  GENOMES_DIR=/path/to/GenomesIndex bash workflows/02_mapping/zhang2024/0_generate_run_plans.sh
# Output: workflows/02_mapping/zhang2024/run_plan.{B73Mo17_rep1,B73Mo17_rep2,multiGenotypes_rep1}.tsv
#         (the shipped copies carry the literal placeholder ${GENOMES_DIR} in ref_path)

set -euo pipefail

REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
CONFIG_DIR="${REPO_ROOT}/workflows/02_mapping/zhang2024"
POOL_DIR="${REPO_ROOT}/config"
INDEX_DIR="${GENOMES_DIR:?set GENOMES_DIR to the directory holding Zea/NAN_Indexes/Index_Zm_<genome>_bwa}/Zea/NAN_Indexes"

generate_plan() {
  local pool_file="$1"
  local genome_list="$2"
  local outfile="$3"

  echo "# sample_base	target_genome	ref_path" > "$outfile"
  while IFS= read -r well; do
    [[ -z "$well" ]] && continue
    while IFS= read -r genome; do
      [[ -z "$genome" ]] && continue
      printf '%s\t%s\t%s/Index_Zm_%s_bwa\n' \
        "$well" "$genome" "$INDEX_DIR" "$genome"
    done < "$genome_list"
  done < "$pool_file" >> "$outfile"

  local n
  n=$(grep -cv '^#' "$outfile" || true)
  echo "  $outfile  ($n tasks)"
}

echo "Generating run_plan TSVs..."

# B73Mo17 rep1: 96 wells x 2 genomes = 192 tasks
generate_plan \
  "${POOL_DIR}/Pools_scifi_B73Mo17_rep1.txt" \
  "${CONFIG_DIR}/Genome_list_scifi_B73_Mo17" \
  "${CONFIG_DIR}/run_plan.B73Mo17_rep1.tsv"

# B73Mo17 rep2: 96 wells x 2 genomes = 192 tasks
generate_plan \
  "${POOL_DIR}/Pools_scifi_B73Mo17_rep2.txt" \
  "${CONFIG_DIR}/Genome_list_scifi_B73_Mo17" \
  "${CONFIG_DIR}/run_plan.B73Mo17_rep2.tsv"

# multiGenotypes rep1: 96 wells x 7 genomes = 672 tasks
generate_plan \
  "${POOL_DIR}/Pools_scifi_multi_genotypes.txt" \
  "${CONFIG_DIR}/Genome_list_MultipleGenome" \
  "${CONFIG_DIR}/run_plan.multiGenotypes_rep1.tsv"

echo "Done."
