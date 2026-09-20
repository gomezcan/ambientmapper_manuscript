#!/usr/bin/env bash
#SBATCH --job-name=sm2v2_dec_nd
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=40G
#SBATCH --time=24:00:00
#SBATCH --output=_logs/01_10b_sm2v2_decontam_nd_%A.log
#
# 01_10b_decontam_SM2v2_without_design.sh — SM2v2 step 7b: decontam (without_design).
#
# Companion to 01_10a_decontam_SM2v2_with_design.sh. The two passes write to
# different output dirs and share no state, so submit both in parallel.
#
# Walltime: same per-chunk rate as the WD pass (~1,785 chunks/h), 31,509 chunks
# -> ~17.6 h projected, 24 h gives headroom.
#
# Output: SM2v2/decontam_without_design_alpha05_v2/ (ND tables; *_decontam_params.json
#         is read by Table S1, reads_to_drop feeds the ND arm of 01_11).
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

sample="SM2v2"

cells_calls="${sample}/final/${sample}_cells_calls.tsv.gz"
[[ -f "$cells_calls" ]] || { echo "Missing: $cells_calls"; exit 1; }

threads="${SLURM_CPUS_PER_TASK:-10}"

assign_glob="${sample}/cell_map_ref_chunks/${sample}_*chunk*_filtered.tsv.gz"
ls $assign_glob >/dev/null 2>&1 || { echo "No assignment files matched: $assign_glob"; exit 1; }

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
echo "[$(date)] SM2v2 step 7b: decontam (without_design)"
echo "  cells_calls : $cells_calls"
echo "  assign_glob : $assign_glob"
echo "  walltime    : 24:00:00, mem=40G, threads=$threads"
echo "============================================================"

out_nod="${sample}/decontam_without_design_alpha05_v2"
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
