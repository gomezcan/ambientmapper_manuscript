#!/usr/bin/env python3
"""
Supplementary Table S4 -- evidence for pooling weak_doublet barcodes with singlets.

Supports the allele-purity paragraph of Methods ("Allele-level validation of decontamination")
and the "Singlet" facet of Fig. 4L, which pools weak_doublet with single_clean and
dirty_singlet. The pooling is a real analysis choice and the paper asserts it in one sentence;
this table is the evidence behind that sentence.

Why the pooling is defensible
-----------------------------
AmbientMapper chooses between a one-genome and a two-genome model by BIC. A weak_doublet is a
barcode where the two-genome model won but then failed a sanity check, either a minor component
below doublet_minor_min or a BIC gap inside near_tie_margin. It is a doublet call the model
itself flags as unconvincing.

At informative sites those barcodes look like singlets:

    pairfrac  = P1 / (P1 + P2)    share of the assigned pair's reads on the dominant genome.
                                  0.5 for a genuine 50/50 doublet, ~1 for a singlet.
    top1or2   = (P1 + P2) / Ptot  fraction of informative reads the assigned pair explains.
    p_top1    = reads on the assigned genome / informative reads.

Inputs (read-only, DATA = data/processed/zhang2024)
    <dataset>/diagnostics/06_48_barcode_purity/<dataset>_weak_doublet_diag.tsv
    for B73Mo17_rep1, B73Mo17_rep2 and multiGenotypes_rep1, produced by the 06_48 barcode-purity
    step of workflows/03b_variant_based_comparison/ (06_48_purity_model.R).

The B73/Mo17 replicates are reported separately and pooled, because Fig. 4L pools them and the
0.02% figure quoted in Methods is the pooled value.

Writes TableS4 as .txt (TSV) and .xlsx into figures/tables/ by default.
Hard-fails if any published value drifts.

Usage (from the repo root)
--------------------------
    python3 analysis/tables/make_TableS4_weak_doublet_pooling.py [--outdir DIR]
"""

import argparse
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
DATA = os.path.join(REPO, "data", "processed", "zhang2024")
DEFAULT_OUTDIR = os.path.join(REPO, "figures", "tables")

# internal dataset key -> published label (presentation rule, no internal identifiers)
DATASETS = [
    ("B73Mo17_rep1", "B73/Mo17, replicate 1"),
    ("B73Mo17_rep2", "B73/Mo17, replicate 2"),
    ("multiGenotypes_rep1", "Multi-genotype"),
]
POOLED_KEYS = ("B73Mo17_rep1", "B73Mo17_rep2")
POOLED_LABEL = "B73/Mo17, replicates pooled"

HEADER = [
    "Dataset",
    "Singlet barcodes",
    "Weak-doublet barcodes",
    "Weak doublets as % of the pooled singlet class",
    "Median dominant-genome share of the assigned pair",
    "Median informative reads explained by the assigned pair",
    "Median reads on the assigned genome, singlets",
    "Median reads on the assigned genome, weak doublets",
]

# Values the manuscript states. The script exists to keep these honest.
PUBLISHED = {
    "pooled_b73mo17_pct": 0.02,      # "0.02% of the B73 and Mo17 singlet class"
    "multi_pct": 99.8,               # "99.8% of the multi-genotype singlet class"
    "multi_pairfrac_pct": 93.0,      # "a median of 93% of the reads assigned to its genome pair"
    "genuine_doublet_pairfrac": 0.5,  # "compared with the 50% expected for a genuine doublet"
}

CAPTION = (
    "Table S4. Evidence for pooling weak-doublet barcodes with singlets in the allele-purity "
    "analysis. AmbientMapper selects between a one-genome and a two-genome model by BIC. A "
    "weak doublet is a barcode for which the two-genome model won that comparison but then "
    "failed a sanity check, either a minor component below threshold or a near-tie in BIC "
    "against the one-genome model, so it is a doublet call the model itself flags as "
    "unconvincing. At informative sites, those at which the pooled genotypes do not all carry "
    "the same allele, these barcodes behave as singlets. The dominant genome takes a median of "
    "93% of the reads assigned to the barcode's genome pair in the multi-genotype library, "
    "against the 50% expected of a genuine equal doublet, and the assigned pair explains "
    "essentially all informative reads. The pooling is therefore negligible in the B73/Mo17 "
    "libraries and dominant in the multi-genotype library, so the two Singlet facets of "
    "Fig. 4L are not comparable in composition and this table records the difference."
)


def load(key):
    """Return {class: {metric: float_or_None}} for one dataset."""
    path = os.path.join(DATA, key, "diagnostics", "06_48_barcode_purity",
                        f"{key}_weak_doublet_diag.tsv")
    if not os.path.exists(path):
        sys.exit(f"missing input: {path}")
    out = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            cls = row["class"]
            out[cls] = {
                k: (float(row[k]) if row.get(k) not in (None, "") else None)
                for k in ("n", "med_p_top1", "med_top1or2", "med_pairfrac")
            }
    for required in ("singlet", "weak_doublet"):
        if required not in out:
            sys.exit(f"{path}: missing class {required!r}")
    return out


def fmt(value, digits=3):
    return "" if value is None else round(value, digits)


def build_row(label, singlet, weak):
    n_s, n_w = int(singlet["n"]), int(weak["n"])
    pct = 100.0 * n_w / (n_s + n_w)
    return [
        label,
        n_s,
        n_w,
        round(pct, 3),
        fmt(weak["med_pairfrac"]),
        fmt(weak["med_top1or2"]),
        fmt(singlet["med_p_top1"]),
        fmt(weak["med_p_top1"]),
    ]


def build_rows(loaded):
    rows = []
    for key, label in DATASETS:
        d = loaded[key]
        rows.append(build_row(label, d["singlet"], d["weak_doublet"]))

    # pooled B73/Mo17 row. Counts add; medians do not, so they are left blank rather
    # than faked by averaging two medians.
    n_s = sum(int(loaded[k]["singlet"]["n"]) for k in POOLED_KEYS)
    n_w = sum(int(loaded[k]["weak_doublet"]["n"]) for k in POOLED_KEYS)
    rows.insert(2, [POOLED_LABEL, n_s, n_w, round(100.0 * n_w / (n_s + n_w), 3), "", "", "", ""])
    return rows


def check(rows):
    """Hard-fail if any value the manuscript states has drifted."""
    by_label = {r[0]: r for r in rows}
    problems = []

    # recompute from the raw counts. Rounding a rounded value re-rounds under
    # banker's rounding and turns 0.0151 into 0.01 rather than 0.02.
    pooled_row = by_label[POOLED_LABEL]
    pooled_pct = round(100.0 * pooled_row[2] / (pooled_row[1] + pooled_row[2]), 2)
    if pooled_pct != PUBLISHED["pooled_b73mo17_pct"]:
        problems.append(
            f"pooled B73/Mo17 weak-doublet share is {pooled_pct}%, "
            f"manuscript says {PUBLISHED['pooled_b73mo17_pct']}%"
        )

    multi = by_label["Multi-genotype"]
    multi_pct = round(100.0 * multi[2] / (multi[1] + multi[2]), 1)
    if multi_pct != PUBLISHED["multi_pct"]:
        problems.append(
            f"multi-genotype weak-doublet share is {multi_pct}%, "
            f"manuscript says {PUBLISHED['multi_pct']}%"
        )

    pairfrac_pct = round(100.0 * float(multi[4]))
    if pairfrac_pct != PUBLISHED["multi_pairfrac_pct"]:
        problems.append(
            f"multi-genotype median pairfrac is {pairfrac_pct}%, "
            f"manuscript says {PUBLISHED['multi_pairfrac_pct']:.0f}%"
        )

    if float(multi[4]) <= PUBLISHED["genuine_doublet_pairfrac"]:
        problems.append(
            "multi-genotype median pairfrac is at or below 0.5, which would break the "
            "argument that weak doublets behave as singlets"
        )

    for key, label in DATASETS:
        row = by_label[label]
        if row[2] == 0:
            problems.append(f"{label}: no weak_doublet barcodes, table would be meaningless")

    if problems:
        sys.exit("SELF-CHECK FAILED\n  - " + "\n  - ".join(problems))


def write_txt(path, rows):
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(HEADER)
        writer.writerows(rows)
        fh.write("\n" + CAPTION + "\n")


def write_xlsx(path, rows):
    try:
        from openpyxl import Workbook
        from openpyxl.styles import Alignment, Font
    except ImportError:
        print("openpyxl not available, skipping .xlsx", file=sys.stderr)
        return False

    wb = Workbook()
    ws = wb.active
    ws.title = "Table S4"
    ws.append(HEADER)
    for cell in ws[1]:
        cell.font = Font(bold=True)
        cell.alignment = Alignment(wrap_text=True, vertical="top")
    for row in rows:
        ws.append(row)
    for row in ws.iter_rows(min_row=2):
        if str(row[0].value) == POOLED_LABEL:
            for cell in row:
                cell.font = Font(bold=True)
        row[3].number_format = "0.000"
        for i in (4, 5, 6, 7):
            row[i].number_format = "0.000"
    for column, width in zip("ABCDEFGH", (26, 16, 20, 22, 22, 22, 20, 22)):
        ws.column_dimensions[column].width = width
    ws.row_dimensions[1].height = 60

    caption_row = ws.max_row + 2
    ws.cell(row=caption_row, column=1, value=CAPTION).alignment = Alignment(
        wrap_text=True, vertical="top"
    )
    ws.merge_cells(start_row=caption_row, start_column=1,
                   end_row=caption_row + 8, end_column=8)
    wb.save(path)
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    loaded = {key: load(key) for key, _ in DATASETS}
    rows = build_rows(loaded)
    check(rows)
    print("self-check against the values stated in Methods: PASS")

    os.makedirs(args.outdir, exist_ok=True)
    txt = os.path.join(args.outdir, "TableS4_weak_doublet_pooling.txt")
    xlsx = os.path.join(args.outdir, "TableS4_weak_doublet_pooling.xlsx")
    write_txt(txt, rows)
    print(f"wrote {txt}")
    if write_xlsx(xlsx, rows):
        print(f"wrote {xlsx}")

    width = max(len(r[0]) for r in rows)
    print()
    for row in rows:
        pf = f"{row[4]:.3f}" if row[4] != "" else "    -"
        print(f"  {row[0]:<{width}}  singlet={row[1]:>7,}  weak={row[2]:>7,}  "
              f"{row[3]:>7.3f}%  pairfrac={pf}")


if __name__ == "__main__":
    main()
