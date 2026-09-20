#!/usr/bin/env bash
#SBATCH --job-name=multi_dec_C0_nd
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=80G
#SBATCH --time=72:00:00
#SBATCH --output=_logs/05_03_decontam_MultiGenotype_C0_nd_%A.log
#
# 05_03_decontam_MultiGenotype_C0_nd.sh — multiGenotypes_rep1 decontam
#                                         (without_design) on C0 (4cfg_2026-05-01).
#
# 73,148 chunks, 7.3M BCs, 7 genomes. Per-chunk cost is ~2-3× the 2-genome
# datasets (more genome columns in the read-to-genome assignment table), so
# projected runtime ~50-70h vs SM2v2's 1,785 chunks/h baseline.
#
# Resource budget: 72h/80G/10cpu — double SM2v2's memory because 7-genome
# panel keeps more state in worker RAM (per the 4cfg ablation lesson: 64G
# OOM'd on multi/Cxmap_mq50 in the first 4cfg run).
#
# Output: multiGenotypes_rep1/decontam_without_design_alpha05_C0/ (reads_to_drop ->
#         05_04 clean-bams -> WASP and purity chain, Fig. 4L to N;
#         *_decontam_params.json -> Table S1).
# =============================================================================

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"   # root of the data tree (3_Mapping/, 5_AmbientDetection/, ...)

cd "${PROJECT_ROOT}/5_AmbientDetection"

# conda activate <env from environment.yml>

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export GOTO_NUM_THREADS=1

sample="multiGenotypes_rep1"

cells_calls="${sample}/genotyping_runs/4cfg_2026-05-01/C0/${sample}_cells_calls.tsv.gz"
[[ -f "$cells_calls" ]] || { echo "Missing: $cells_calls"; exit 1; }

threads="${SLURM_CPUS_PER_TASK:-10}"

assign_glob="${sample}/cell_map_ref_chunks/*_filtered.tsv.gz"
n_filtered=$(find "${sample}/cell_map_ref_chunks" -maxdepth 1 -name '*_filtered.tsv.gz' | wc -l)
[[ "$n_filtered" -gt 0 ]] || { echo "No assignment files matched: $assign_glob"; exit 1; }
echo "  found $n_filtered filtered chunk files"

alpha="0.05"
chunksize="1000000"

min_reads_post_clean="100"
min_allowed_frac_post_clean="0.90"
# Both post-clean thresholds are reporting flags in the AmbientMapper build used
# here, not enforced gates (see the stage README).

doublet_policy="top1"
indist_policy="top1"

safe_keep_delta_as="3"

ambiguous_policy_no_design="top1_rescue"

echo "============================================================"
echo "[$(date)] multiGenotypes_rep1 decontam (C0, without_design)"
echo "  cells_calls : $cells_calls"
echo "  assign_glob : $assign_glob"
echo "  walltime    : 72:00:00, mem=80G, threads=$threads"
echo "============================================================"

out_nod="${sample}/decontam_without_design_alpha05_C0"
mkdir -p "$out_nod"

ambientmapper decontam \
  --cells-calls "$cells_calls" \
  --out-dir "$out_nod" \
  --threads "$threads" \
  --assign-glob "$assign_glob" \
  --ambiguous-policy "$ambiguous_policy_no_design" \
  --doublet-policy "$doublet_policy" \
  --indist-policy "$indist_policy" \
  --decontam-alpha "$alpha" \
  --chunksize "$chunksize" \
  --min-reads-post-clean "$min_reads_post_clean" \
  --min-allowed-frac-post-clean "$min_allowed_frac_post_clean" \
  --safe-keep-delta-as "$safe_keep_delta_as"

echo ""
echo "============================================================"
echo "[$(date)] without-design DONE for ${sample}"
echo "  output: ${out_nod}"
ls -la "${out_nod}/" | head -20
echo "============================================================"
