#!/bin/bash
########## BATCH Lines for Resource Request ##########
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=fig5_L_pgt
#SBATCH --partition=standard
#SBATCH --output=figures/main/fig5/browser/_logs_plots/%x-%A.log
# =============================================================================
# browser/plot.sh  -  STEP 5 of the Fig 5 panel L genome-browser chain (HPC only): pyGenomeTracks
#   figures at the audition loci, Pre | wd | nd per group. The four shipped renders of panel L are
#   Zm00001eb115210.companion_cells, zyb9.leaf_primordia, AT2G37170.bundle_sheath, SCL23.bundle_sheath.
# Inputs  figures/main/fig5/browser/plot_loci.tsv (browser/loci.R), browser/BWs/*.bw (browser/bw.sh),
#         $DATA/socrates/_data/_GenomeInfo/{B73v5,TAIR10}.gff3   (gene models)
# Output  figures/main/fig5/browser/plots/<locus>.<group>.{pdf,png} (+ plots/_ini/*.ini kept for hand-tweaking)
# Requires the pygenometracks conda env (browser/env.sh, one-time).
# Usage   (from the repo root, on the HPC)  mkdir -p figures/main/fig5/browser/_logs_plots
#         sbatch analysis/fig5_biological_impact/browser/plot.sh     # or: bash .../browser/plot.sh
# =============================================================================
#
# Per locus the figure stacks:
#     [ gene models ]             from the reference GFF3 (as a BED6 track, cached under browser/annot/)
#     [ <group> Pre / wd / nd ]   stage colours as every Fig 5 panel (#FF83FA / #43CD80 / #E69F00),
#                                 ONE SHARED y max across the three stages -> the cross-stage
#                                 comparison is real
#     [ ALLcells Pre / wd / nd ]  grey, their own shared y max (a different scale: the aggregate
#                                 dilutes any one type)
# More loci = append rows to plot_loci.tsv and resubmit; nothing else changes.
# COORDINATES COME FROM THE GFF3, not the marker BED (vnd7's marker-BED entry sits 148 kb from
#   its GFF3 gene body). The loci TSV was built from the GFF3; keep it that way.
# Track names are ARGMAX type calls (descriptive annotation) and each stage's group is a
#   different cell set from a different partition. RPM is per track.
# =============================================================================
set -euo pipefail

# conda activate pygenometracks   # the env created by browser/env.sh

REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
cd "$REPO_ROOT"

DATA=data/processed/scifiATAC_B73_Arabidopsis
BROWSER=figures/main/fig5/browser
BWDIR="$BROWSER/BWs"
LOCI="$BROWSER/plot_loci.tsv"
GINFO="$DATA/socrates/_data/_GenomeInfo"
ANNOT="$BROWSER/annot"
PLOTS="$BROWSER/plots"
mkdir -p "$ANNOT" "$PLOTS/_ini"

[ -f "$LOCI" ] || { echo "missing loci config: $LOCI (run browser/loci.R first)"; exit 1; }

# --- gene-model BED per genome, cached (bed6 from the gene rows of the GFF3;
#     a BED track sidesteps every GFF/GTF parser quirk in pyGenomeTracks) -----
gff_of() {
  case "$1" in
    B73v5)  echo "$GINFO/B73v5.gff3"  ;;
    TAIR10) echo "$GINFO/TAIR10.gff3" ;;
    *) echo "unknown genome: $1" >&2; return 1 ;;
  esac
}
for g in B73v5 TAIR10; do
  out="$ANNOT/${g}.genes.bed"
  [ -s "$out" ] && continue
  gff=$(gff_of "$g")
  [ -f "$gff" ] || { echo "missing GFF3 (place the reference gene annotation there): $gff"; exit 1; }
  awk -F'\t' -v OFS='\t' '$3 == "gene" {
    name = $9
    if (match(name, /Name=[^;]+/))    name = substr(name, RSTART + 5, RLENGTH - 5)
    else if (match($9, /ID=[^;]+/)) { name = substr($9, RSTART + 3, RLENGTH - 3)
                                      sub(/^gene:/, "", name) }
    print $1, $4 - 1, $5, name, 0, $7
  }' "$gff" | sort -k1,1 -k2,2n > "$out"
  echo "built $out ($(wc -l < "$out" | tr -d ' ') genes)"
done

# --- render every locus ------------------------------------------------------
python3 - "$LOCI" "$BWDIR" "$ANNOT" "$PLOTS" <<'PYEOF'
import csv, math, os, subprocess, sys
import pyBigWig

loci_tsv, bwdir, annot, plots = sys.argv[1:5]
ST      = [("SM2", "Pre", "#FF83FA"), ("Clean.SM2v2wd", "wd", "#43CD80"),
           ("Clean.SM2v2", "nd", "#E69F00")]          # the Fig 5 stage palette
SUFFIX  = {"B73v5": "B73_B73v5", "TAIR10": "At_TAIR10"}
GREY    = "#999999"

def region_max(files, chrom, start, end, nbins):
    """max of the BINNED means pyGenomeTracks will actually draw, +8%.
    NOT the bp-level max: that is set by single-bp spikes which bin-averaging
    flattens (measured: ccdp_1's bp max, 160 RPM, was ~2x the tallest drawn bin
    and squashed every track into the lower half of its panel). Same bin count
    as the tracks, exact (no zoom-level estimates)."""
    m = 0.0
    for f in files:
        bw = pyBigWig.open(f)
        if chrom not in bw.chroms():
            sys.exit(f"chrom {chrom!r} not in {f} - naming mismatch")
        vals = bw.stats(chrom, start, end, type="mean", nBins=nbins, exact=True)
        m = max(m, max((v or 0.0) for v in vals))
        bw.close()
    if m <= 0:
        return 1.0
    return float(f"{m * 1.08:.2g}")

with open(loci_tsv) as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))
if not rows:
    sys.exit("empty loci config")

for r in rows:
    locus, genome, group = r["locus"], r["genome"], r["group"]
    chrom = r["chrom"]; pad = int(r["pad"])
    start = max(0, int(r["start"]) - pad); end = int(r["end"]) + pad
    suffix = SUFFIX[genome]

    grp_bw = [f"{bwdir}/{p}_{suffix}.{group}.bw"   for p, _, _ in ST]
    all_bw = [f"{bwdir}/{p}_{suffix}.ALLcells.bw"  for p, _, _ in ST]
    missing = [f for f in grp_bw + all_bw if not os.path.isfile(f)]
    if missing:
        sys.exit("missing bigwig(s) - run browser/bw.sh first:\n  "
                 + "\n  ".join(missing))

    # bins tied to the window: ~8-10 bp/bin. A fixed 2000 over the ~1 kb
    # padded windows would be ~2 bp/bin and render spiky.
    nbins = max(300, min(2000, (end - start) // 8))
    ymax_grp = region_max(grp_bw, chrom, start, end, nbins)
    ymax_all = region_max(all_bw, chrom, start, end, nbins)

    ini = [f"# {locus} ({genome}) - {r.get('note','')}",
           "[x-axis]", "fontsize = 8", "", "[spacer]", "height = 0.15", ""]
    ini += [f"[genes]", "file = " + f"{annot}/{genome}.genes.bed",
            "title = genes", "height = 0.5",     # thin: one gene + neighbours
            "color = #2b2b2b",
            "border_color = #2b2b2b", "labels = true", "fontsize = 7",
            "file_type = bed", "display = stacked", "", "[spacer]",
            "height = 0.2", ""]
    for (pre, stage, col), f in zip(ST, grp_bw):
        ini += [f"[{group} {stage}]", f"file = {f}",
                f"title = {group} {stage}", "height = 1.5",
                f"color = {col}", "min_value = 0",
                f"max_value = {ymax_grp}",       # shared: honest across stages
                f"number_of_bins = {nbins}", "file_type = bigwig", ""]
    ini += ["[spacer]", "height = 0.2", ""]
    for (pre, stage, col), f in zip(ST, all_bw):
        ini += [f"[ALLcells {stage}]", f"file = {f}",
                f"title = ALLcells {stage}", "height = 1.0",
                f"color = {GREY}", "min_value = 0",
                f"max_value = {ymax_all}",       # their own shared scale
                f"number_of_bins = {nbins}", "file_type = bigwig", ""]

    ini_f = f"{plots}/_ini/{locus}.{group}.ini"
    with open(ini_f, "w") as fh:
        fh.write("\n".join(ini))

    region = f"{chrom}:{start}-{end}"
    ttl = (f"{locus} | {genome} {region} | group ymax {ymax_grp} RPM, "
           f"ALL ymax {ymax_all} RPM")
    for ext in ("pdf", "png"):
        out = f"{plots}/{locus}.{group}.{ext}"
        subprocess.run(["pyGenomeTracks", "--tracks", ini_f,
                        "--region", region, "--title", ttl,
                        "--width", "22", "--dpi", "200",
                        "--outFileName", out], check=True)
        if not (os.path.isfile(out) and os.path.getsize(out) > 1000):
            sys.exit(f"FAILED to write {out}")
    print(f"[done] {locus}.{group}  {region}  "
          f"(group ymax {ymax_grp}, ALL ymax {ymax_all})")
PYEOF

echo
echo "plots + kept .ini files -> $PLOTS/"
