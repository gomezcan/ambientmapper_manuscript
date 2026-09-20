#!/bin/bash
########## BATCH Lines for Resource Request ##########
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --job-name=Merge_SocObjs_SM2v2
#SBATCH --partition=standard
#SBATCH --output=_logs/Step2_merge_SocObjs_SM2v2_%x-%j.log
#
# 2_0_0_Obj_integreation_SM2v2.sh
# Merge per-species Socrates objects into a single SocObj for the SM2v2 story.
# Two passes, one per stage:
#   - PreClean : <BASE>/SM2/step0_qc/{SM2_At,SM2_B73}.raw.soc.rds + .updated_metadata.txt
#   - PostClean: <BASE>/SM2v2_clean/step0_qc/{Clean.SM2v2_At,Clean.SM2v2_B73}.raw.soc.rds + .updated_metadata.txt
#
# 2_0_0_Obj_integreation.R auto-detects species (_At_ / _B73_) and stage
# (Clean. prefix → PostClean) from filenames, so no R-side changes needed.

set -euo pipefail

# conda activate ambientmapper-manuscript   (environment.yml at the repo root)

# Paths are anchored on PROJECT_ROOT ($BASH_SOURCE self-location does not survive sbatch).
PROJECT_ROOT="${PROJECT_ROOT:-.}"
BASE="${PROJECT_ROOT}/6_socrates"
SCRIPTS="${SCRIPTS:-${BASE}/0_scripts}"
cd "${BASE}"

# Uncomment whichever stage you need to merge (or both):
Rscript "${SCRIPTS}/part1_combined/2_0_0_Obj_integreation.R" SM2         SM2
Rscript "${SCRIPTS}/part1_combined/2_0_0_Obj_integreation.R" SM2v2_clean Clean.SM2v2

echo ' .. done ..'
