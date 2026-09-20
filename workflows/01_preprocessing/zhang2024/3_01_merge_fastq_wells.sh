#!/usr/bin/env bash
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --job-name=merge_wells
#SBATCH --partition=standard
#SBATCH --output=_logs/%x_%j.log

# =============================================================================
# 3_01_merge_fastq_wells.sh
#
# Two-phase merge of demuxed scifi-ATAC FASTQs:
#   Phase 1: Within-run merge (100 chunk parts -> 1 file per well per run)
#            Uses: scifi-demux step1 merge
#   Phase 2: Cross-run merge  (3 sequencing runs -> 1 file per well per rep)
#            Uses: cat (gzip concatenation)
#
# Usage (from ${PROJECT_ROOT}/2_CleanReads):
#   sbatch 3_01_merge_fastq_wells.sh
#   # or run locally:
#   bash  3_01_merge_fastq_wells.sh
#
# Input:  01_scifi_{LIB}_{RUN}/Corrected/part_*_R{1,3}.bc1.bc2_{WELL}.fastq.gz
# Output: combined/{DATASET}_rep{N}/{WELL}_R{1,3}.bc1.bc2.fastq.gz
# =============================================================================

set -euo pipefail
: "${PROJECT_ROOT:?set PROJECT_ROOT to the analysis tree that holds 2_CleanReads/}"
# conda activate <env from environment.yml>

BASEDIR="${PROJECT_ROOT}/2_CleanReads"
cd "$BASEDIR"
mkdir -p _logs

OUTROOT="combined"

# ---------------------------------------------------------------------------
# Dataset definitions:  DATASET  REP  RUN_SUFFIXES
#   Each biological replicate was deposited as 3 SRA runs (_1, _2, _3).
#   Work roots: 01_scifi_{DATASET}_rep{REP}_{RUN}/
# ---------------------------------------------------------------------------
DATASETS=(
  "B73Mo17:rep1:1,2,3"
  "B73Mo17:rep2:1,2,3"
  "multiGenotypes:rep1:1,2,3"
)

# ===========================
# Phase 1: Within-run merge
# ===========================
echo "=============================================="
echo "[$(date)] Phase 1: within-run merge (parts -> per-well)"
echo "=============================================="

for entry in "${DATASETS[@]}"; do
  IFS=':' read -r dataset rep runs <<< "$entry"
  IFS=',' read -ra run_arr <<< "$runs"

  for run in "${run_arr[@]}"; do
    workroot="01_scifi_${dataset}_${rep}_${run}"

    if [[ ! -d "$workroot" ]]; then
      echo "[WARN] Missing work_root: $workroot, skipping"
      continue
    fi

    # Skip if already merged
    if [[ -d "${workroot}/combined" ]] && [[ $(ls "${workroot}/combined/"*.fastq.gz 2>/dev/null | wc -l) -gt 0 ]]; then
      echo "[$(date)] ${workroot}/combined/ already exists, skipping phase 1"
      continue
    fi

    # Extract library name from run plan
    lib=$(awk -F'\t' 'NR==2{print $2}' "${workroot}/run_plan.step1.chunks.tsv")

    echo "[$(date)] Merging parts in ${workroot} (library=${lib})..."
    scifi-demux step1 merge --library "$lib" --work-root "$workroot"
    echo "[$(date)] Done: ${workroot}/combined/"
  done
done

# ===========================
# Phase 2: Cross-run merge
# ===========================
echo ""
echo "=============================================="
echo "[$(date)] Phase 2: cross-run merge (runs -> per-replicate)"
echo "=============================================="

for entry in "${DATASETS[@]}"; do
  IFS=':' read -r dataset rep runs <<< "$entry"
  IFS=',' read -ra run_arr <<< "$runs"

  OUTDIR="${OUTROOT}/${dataset}_${rep}"
  mkdir -p "$OUTDIR"

  # Get well list from the first run
  first_run="01_scifi_${dataset}_${rep}_${run_arr[0]}/combined"
  if [[ ! -d "$first_run" ]]; then
    echo "[ERROR] Phase 1 output missing: ${first_run}, skipping ${dataset}_${rep}"
    continue
  fi

  # Extract unique well names from combined/ filenames
  # Pattern: {WELL}_R{1,3}.bc1.bc2.fastq.gz
  wells=()
  while IFS= read -r w; do
    wells+=("$w")
  done < <(ls "$first_run"/*_R1.bc1.bc2.fastq.gz 2>/dev/null \
    | xargs -I{} basename {} \
    | sed 's/_R1\.bc1\.bc2\.fastq\.gz$//' \
    | sort -u)

  if [[ ${#wells[@]} -eq 0 ]]; then
    echo "[ERROR] No wells found in ${first_run}, skipping"
    continue
  fi

  echo "[$(date)] Merging ${#wells[@]} wells across ${#run_arr[@]} runs -> ${OUTDIR}/"

  for well in "${wells[@]}"; do
    for read in R1 R3; do
      outfile="${OUTDIR}/${well}_${read}.bc1.bc2.fastq.gz"

      # Skip if already exists
      if [[ -f "$outfile" ]]; then
        continue
      fi

      # Collect inputs from all runs
      inputs=()
      for run in "${run_arr[@]}"; do
        src="01_scifi_${dataset}_${rep}_${run}/combined/${well}_${read}.bc1.bc2.fastq.gz"
        if [[ -f "$src" ]]; then
          inputs+=("$src")
        else
          echo "[WARN] Missing: $src"
        fi
      done

      if [[ ${#inputs[@]} -eq 0 ]]; then
        echo "[WARN] No inputs for ${well} ${read}, skipping"
        continue
      fi

      # gzip files can be concatenated directly
      cat "${inputs[@]}" > "$outfile"
    done
  done

  echo "[$(date)] Done: ${OUTDIR}/ (${#wells[@]} wells x 2 reads)"
done

# ===========================
# Summary
# ===========================
echo ""
echo "=============================================="
echo "[$(date)] Merge complete. Summary:"
echo "=============================================="
for entry in "${DATASETS[@]}"; do
  IFS=':' read -r dataset rep runs <<< "$entry"
  outdir="${OUTROOT}/${dataset}_${rep}"
  if [[ -d "$outdir" ]]; then
    n_files=$(ls "$outdir"/*.fastq.gz 2>/dev/null | wc -l)
    n_wells=$(ls "$outdir"/*_R1.bc1.bc2.fastq.gz 2>/dev/null | wc -l)
    size=$(du -sh "$outdir" | cut -f1)
    echo "  ${dataset}_${rep}: ${n_wells} wells, ${n_files} files, ${size}"
  else
    echo "  ${dataset}_${rep}: MISSING"
  fi
done
