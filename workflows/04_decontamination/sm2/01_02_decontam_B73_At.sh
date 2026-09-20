#!/usr/bin/env bash
# 01_02_decontam_B73_At.sh — SM2 (combined-reference arm) decontamination, with-design (WD) and
# without-design (ND) passes on the 01_01 genotyping output.
# Input : <sample>/final/<sample>_cells_calls.tsv.gz, <sample>/cell_map_ref_chunks/*_filtered.tsv.gz,
#         configs/PlateDesign_<sample>_ATAC.txt (WD pass only)
# Output: <sample>/decontam_with_design_alpha05_v2/ (the SM2 tables read by Fig. 1B to E)
#         <sample>/decontam_without_design_alpha05_v2/ (no manuscript consumer for SM2)
# Run   : cd ${PROJECT_ROOT}/5_AmbientDetection && sbatch <repo>/workflows/04_decontamination/sm2/01_02_decontam_B73_At.sh SM2
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=40G
#SBATCH --job-name=ambientmapper_SM2_step2
#SBATCH --partition=standard
#SBATCH --output=_logs/01_02_decontam_B73_At_%A.log

if [[ $# -lt 1 ]]; then
  echo "Usage: sbatch $0 <sample_dir_name>"
  exit 2
fi
sample="$1"

PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"   # root of the data tree (3_Mapping/, 5_AmbientDetection/, ...)
cd "${PROJECT_ROOT}/5_AmbientDetection"

# conda activate <env from environment.yml>

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export GOTO_NUM_THREADS=1

cells_calls="${sample}/final/${sample}_cells_calls.tsv.gz"
[[ -f "$cells_calls" ]] || { echo "Missing: $cells_calls"; exit 1; }

threads="${SLURM_CPUS_PER_TASK:-5}"

assign_glob="${sample}/cell_map_ref_chunks/${sample}_chunk*_filtered.tsv.gz"
# fail early if glob matches nothing
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

# ----------------------------
# WITH DESIGN
# ----------------------------
design_file="configs/PlateDesign_${sample}_ATAC.txt"
[[ -f "$design_file" ]] || { echo "Missing design file: $design_file"; exit 1; }

ambiguous_policy="design_rescue"

out_design="${sample}/decontam_with_design_alpha05_v2"
mkdir -p "$out_design"

ambientmapper decontam \
  --cells-calls "$cells_calls" \
  --out-dir "$out_design" \
  --threads "$threads" \
  --assign-glob "$assign_glob" \
  --design-file "$design_file" \
  --design-bc-mode last \
  --design-bc-n 10 \
  --strict-design-drop-mismatch \
  --ambiguous-policy "$ambiguous_policy" \
  --doublet-policy "$doublet_policy" \
  --indist-policy "$indist_policy" \
  --decontam-alpha "$alpha" \
  --chunksize "$chunksize" \
  --min-reads-post-clean "$min_reads_post_clean" \
  --min-allowed-frac-post-clean "$min_allowed_frac_post_clean" \
  --safe-keep-delta-as "$safe_keep_delta_as"

echo "done with design"

# ----------------------------
# WITHOUT DESIGN
# ----------------------------
# conservative choice would be: ambiguous_policy_no_design="drop"
ambiguous_policy_no_design="top1_rescue"

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

echo "done without design"

