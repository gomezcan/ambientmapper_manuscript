#!/usr/bin/env bash
#SBATCH --job-name=06_47_analyze
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=_logs/06_47_%x_%j.log
#
# Chunked WASP — STEP 4/4: concordance on the WASP-corrected counts (both modes),
# reusing 06_44_concordance.R + the 06_46 WASP-vs-noWASP compare block, unchanged.
# Thin tail — identical to 06_46's analyze, pointed at the chunked 06_47 counts.
# Run AFTER both modes' combines have landed.
#
#   sbatch --export=ALL,SAMPLE=B73Mo17_rep1 workflows/03b_variant_based_comparison/06_47_wasp_chunked_analyze.sh
# =============================================================================
set -euo pipefail

REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"   # this repository (sourced helpers, sibling scripts)
PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"   # root of the data tree (3_Mapping/, 5_AmbientDetection/, ...)
# SLURM copies the batch script to a spool dir, so $BASH_SOURCE can't find siblings;
# the shared helper is located through REPO_ROOT instead.
source "${REPO_ROOT}/workflows/03b_variant_based_comparison/06_47_wasp_chunked_common.sh"
unset PYTHONPATH PYTHONHOME   # pure base py/R; defensive

RAW=$OUT/${SAMPLE}_raw_wasp_allele_counts.tsv.gz
CLN=$OUT/${SAMPLE}_clean_wasp_allele_counts.tsv.gz
for f in "$RAW" "$CLN"; do [[ -s "$f" ]] || { echo "Missing wasp counts (run combine for both modes first): $f"; exit 1; }; done

PANEL=$OUT/genotype_panel.tsv
if [[ ! -s "$PANEL" ]]; then
  if [[ -s "$OUT44/genotype_panel.tsv" ]]; then
    cp "$OUT44/genotype_panel.tsv" "$PANEL"       # same VCF -> identical panel
  else
    bcftools query -s "$GENOTYPES" -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$EXTRACT_VCF" > "$PANEL"
  fi
fi

echo "[06_47/analyze] concordance on chunked-WASP counts sample=$SAMPLE"
$RSCRIPT "${REPO_ROOT}/workflows/03b_variant_based_comparison/06_44_concordance.R" \
  --counts-raw "$RAW" --counts-clean "$CLN" \
  --genotype-panel "$PANEL" --genotypes "$GENOTYPES" \
  --out-prefix "$OUT/${SAMPLE}_wasp" \
  --min-reads 2 --n-boot 200 --block-mb 1

# ---- side-by-side: chunked-WASP vs non-WASP (refbias + self-concordance) -----------
echo
echo "=================  chunked-WASP vs non-WASP compare ($SAMPLE)  ================="
$PY - "$SAMPLE" "$OUT" "$OUT44" <<'PY'
import sys, os, csv
sample, out, out44 = sys.argv[1], sys.argv[2], sys.argv[3]

def read_bias(path):
    d={}
    if not os.path.exists(path): return d
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            n=int(row['n_obs']); nm=int(row['sum_nm'])
            d[(row['group'],row['obs_class'])]=(n, nm/n if n else float('nan'))
    return d

print("\n-- ref-mapping bias: mean NM per allele class (want ALT->REF gap to shrink) --")
print(f"{'mode':6} {'group':6} {'class':4} {'n_obs':>12} {'mean_nm':>8}")
for tag, path in [("noWASP", f"{out44}/{sample}_raw_allele_counts.tsv.bias.tsv"),
                  ("WASP",   f"{out}/{sample}_raw_wasp_allele_counts.tsv.bias.tsv")]:
    b=read_bias(path)
    for (g,c),(n,mnm) in sorted(b.items()):
        print(f"{tag:6} {g:6} {c:4} {n:12,d} {mnm:8.3f}")

def read_selfconc(path):
    rows=[]
    if not os.path.exists(path): return rows
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'): rows.append(row)
    return rows

print("\n-- self-concordance (chunked-WASP-corrected) -- full CIs in the tsv --")
for row in read_selfconc(f"{out}/{sample}_wasp_self_concordance_ci.tsv"):
    print("  ", {k:row[k] for k in row})
print(f"\n(compare against non-WASP: {out44}/{sample}_self_concordance_ci.tsv)")
PY
echo "[06_47/analyze] done -> $OUT/${SAMPLE}_wasp_*"
