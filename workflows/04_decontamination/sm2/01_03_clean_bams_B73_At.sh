#!/usr/bin/env bash
# 01_03_clean_bams_B73_At.sh — SM2 (combined-reference arm) clean-bams: drop the WD reads_to_drop list from
# the combined-genome (ZmATcombined) SM2 BAMs and derive Tn5 insertion BEDs.
# Input : <sample>/decontam_with_design_alpha05_v2/<sample>_reads_to_drop.tsv.gz (from 01_02),
#         4_MappingCleaning/<sample>_{B73,At}/<sample>_*_ZmATcombined_scifiATAC.mq10.BC.rmdup.mm.bam
# Output: <sample>/clean_bams_alpha05/*.Clean.bam + *.tn5.bed.gz
# Note  : 4_MappingCleaning/ was later consolidated into 3_Mapping/; the same BAMs now sit
#         under 3_Mapping/SM2_{At,B73}/ and the Tn5 BED script under 3_Mapping/_archive/.
# Run   : cd ${PROJECT_ROOT}/5_AmbientDetection && sbatch <repo>/workflows/04_decontamination/sm2/01_03_clean_bams_B73_At.sh SM2

########## BATCH Lines for Resource Request ##########
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=40G
#SBATCH --job-name=ambientmapper_SM2_step3
#SBATCH --partition=standard
#SBATCH --output=_logs/01_03_clean_bams_B73_At_%A.log


###################################
#######   Conda / Modules   #######
###################################
# conda activate <env from environment.yml>


sample="${1:?Usage: sbatch this_script.sh <SAMPLE>}"

PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"   # root of the data tree (3_Mapping/, 5_AmbientDetection/, ...)
cd "${PROJECT_ROOT}/5_AmbientDetection"
DIR_INPUT="${PROJECT_ROOT}/4_MappingCleaning"
OUTDIR="${sample}/clean_bams_alpha05"
DROP="${sample}/decontam_with_design_alpha05_v2/${sample}_reads_to_drop.tsv.gz"
BAM_TO_BED="${PROJECT_ROOT}/4_MappingCleaning/1_6_scifi_makeTn5bed.py"

mkdir -p "$OUTDIR" "_logs"

b73_in="${DIR_INPUT}/${sample}_B73/${sample}_B73_ZmATcombined_scifiATAC.mq10.BC.rmdup.mm.bam"
at_in="${DIR_INPUT}/${sample}_At/${sample}_At_ZmATcombined_scifiATAC.mq10.BC.rmdup.mm.bam"

# If clean-bams preserves basename and appends ".<suffix>.bam"
b73_base="$(basename "$b73_in" .bam)"
at_base="$(basename "$at_in" .bam)"
b73_clean="${OUTDIR}/${b73_base}.Clean.bam"
at_clean="${OUTDIR}/${at_base}.Clean.bam"

ambientmapper clean-bams \
	--reads-to-drop "$DROP" \
	--bam "$b73_in" \
	--out-dir "$OUTDIR" \
	--out-suffix .Clean.bam

echo "[clean-bams] B73 cleaned -> $b73_clean"

ambientmapper clean-bams \
  --reads-to-drop "$DROP" \
  --bam "$at_in" \
  --out-dir "$OUTDIR" \
  --out-suffix .Clean.bam
echo "[clean-bams] At cleaned -> $at_clean"

# 2) make Tn5 BED (B73)
python "$BAM_TO_BED" "$b73_clean" \
  | sort -k1,1 -k2,2n \
  | uniq \
  > "${b73_clean%.bam}.tn5.bed"

pigz -p 5 "${b73_clean%.bam}.tn5.bed"
echo "[tn5] wrote ${b73_clean%.bam}.tn5.bed.gz"

# 3) make Tn5 BED (At)
python "$BAM_TO_BED" "$at_clean" \
  | sort -k1,1 -k2,2n \
  | uniq \
  > "${at_clean%.bam}.tn5.bed"
pigz -p 5 "${at_clean%.bam}.tn5.bed"
echo "[tn5] wrote ${at_clean%.bam}.tn5.bed.gz"


