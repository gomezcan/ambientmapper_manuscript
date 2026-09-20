#!/usr/bin/env python3
"""
Supplementary Table S1 -- library and plate design.

Cited by the main Methods, "scifi-ATAC-seq library construction": "we prepared each sample by
separating and tagmenting it across 32 independent wells per nuclei extraction, with two samples
processed per plate (Table S1)".

What the table answers
----------------------
For every scifi-ATAC-seq library in the paper: which plate indexes each sample occupies, how many
independent pools it spans, and whether the plate index identifies the sample. That last column is
the basis of the WD versus ND distinction, so it is generated from the run configuration rather
than asserted.

Plate index is NOT the same as per-well demultiplexing
------------------------------------------------------
The B73/Mo17 libraries were split per well before mapping, but that split is a parallel-processing
convenience, not an experimental design. Both genotypes are pooled before tagmentation and the
analysis treats the whole plate as one sample. The table therefore reports the whole plate for
those libraries and does not present the per-well split as a design.

Scope
-----
scifi-ATAC-seq libraries only. The maize root library is a 10x scATAC library with no combinatorial
plate, so it has no design to report and appears in Table S2 instead. The synthetic benchmark is
simulated and has no physical plate.

Design of record for the interspecies library
---------------------------------------------
`config/PlateDesign_SM2_ATAC.txt` is the two-block design (maize, Arabidopsis) that was supplied to
AmbientMapper and that matches the released data. Do NOT substitute the three-block scifi-demux
plate design (`PlateDesign_scifi_At_B73_rep1.txt`): it carries a third block that is not part of
this study, and the script hard-fails if it is passed one.

Inputs (all read-only)
----------------------
    config/PlateDesign_SM2_ATAC.txt
    config/PlateDesign_scifi_B73Mo17_rep{1,2}.txt
    config/Well_to_Genotype_multiGenotypes_rep1.txt
    data/processed/<dataset>/.../*_decontam_params.json   (one per decontamination run, see DECONTAM_PARAMS)

Writes TableS1 as .txt (TSV) and .xlsx into figures/tables/ by default.

Usage (from the repo root)
--------------------------
    python3 analysis/tables/make_TableS1_library_design.py [--outdir DIR]
"""

import argparse
import csv
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
CONFIG_DIR = os.path.join(REPO, "config")
DATA = os.path.join(REPO, "data", "processed")
DEFAULT_OUTDIR = os.path.join(REPO, "figures", "tables")

SM2_DESIGN = os.path.join(CONFIG_DIR, "PlateDesign_SM2_ATAC.txt")
B73MO17_DESIGN = os.path.join(CONFIG_DIR, "PlateDesign_scifi_B73Mo17_rep{rep}.txt")
MULTI_WELLMAP = os.path.join(CONFIG_DIR, "Well_to_Genotype_multiGenotypes_rep1.txt")

# Runs whose *_decontam_params.json decides the "design supplied" column (paths under DATA).
DECONTAM_PARAMS = {
    ("interspecies", "WD"): (
        "scifiATAC_B73_Arabidopsis/SM2v2/decontam_with_design_alpha05_v2/SM2v2_decontam_params.json"
    ),
    ("interspecies", "ND"): (
        "scifiATAC_B73_Arabidopsis/SM2v2/decontam_without_design_alpha05_v2/SM2v2_decontam_params.json"
    ),
    ("b73mo17_rep1", "ND"): (
        "zhang2024/B73Mo17_rep1/decontam_without_design_alpha05_C0/B73Mo17_rep1_decontam_params.json"
    ),
    ("b73mo17_rep2", "ND"): (
        "zhang2024/B73Mo17_rep2/decontam_without_design_alpha05_C0/B73Mo17_rep2_decontam_params.json"
    ),
    ("multigenotype", "ND"): (
        "zhang2024/multiGenotypes_rep1/decontam_without_design_alpha05_C0/"
        "multiGenotypes_rep1_decontam_params.json"
    ),
}
DECONTAM_ROOT = DATA

# Published library names. Internal identifiers are deliberately absent.
LIBRARY_LABEL = {
    "interspecies": "Maize and Arabidopsis",
    "b73mo17_rep1": "B73/Mo17, replicate 1",
    "b73mo17_rep2": "B73/Mo17, replicate 2",
    "multigenotype": "Multi-genotype",
}

SM2_SAMPLE_LABEL = {"B73": "Maize B73", "At": "Arabidopsis"}

# Verified structure. The script self-checks against this so a silent config change cannot ship.
EXPECTED = {
    "interspecies": {"samples": 2, "wells_per_sample": 32, "total": 64},
    "b73mo17_rep1": {"total": 96},
    "b73mo17_rep2": {"total": 96},
    "multigenotype": {"genotypes": 7, "pools": 8, "wells_per_pool": 12, "total": 96},
}

WELL_RE = re.compile(r"^([A-H])(\d+)(?:-(\d+))?$")
ROWS = "ABCDEFGH"


def expand_wells(spec):
    """'A1-4,B1-4' -> ['A1','A2','A3','A4','B1',...]. Raises on anything unparseable."""
    wells = []
    for token in spec.split(","):
        token = token.strip()
        if not token:
            continue
        m = WELL_RE.match(token)
        if not m:
            raise ValueError(f"unparseable well token: {token!r}")
        row, start, end = m.group(1), int(m.group(2)), m.group(3)
        end = int(end) if end else start
        if not 1 <= start <= end <= 12:
            raise ValueError(f"well token out of range: {token!r}")
        wells.extend(f"{row}{i}" for i in range(start, end + 1))
    return wells


def _join(items):
    items = list(items)
    if len(items) == 1:
        return items[0]
    return ", ".join(items[:-1]) + " and " + items[-1]


def _int_ranges(values):
    """[1,2,3,4,7] -> '1 to 4 and 7'."""
    values = sorted(values)
    groups, run = [], [values[0]]
    for v in values[1:]:
        if v == run[-1] + 1:
            run.append(v)
        else:
            groups.append(run)
            run = [v]
    groups.append(run)
    return _join(f"{g[0]} to {g[-1]}" if len(g) > 1 else str(g[0]) for g in groups)


def describe_wells(wells):
    """Compact human-readable plate position for a set of wells."""
    wells = sorted(set(wells), key=lambda w: (w[0], int(w[1:])))
    rows = sorted({w[0] for w in wells})
    cols = sorted({int(w[1:]) for w in wells})
    if set(wells) != {f"{r}{c}" for r in rows for c in cols}:
        return ", ".join(wells)
    if len(rows) == 8 and len(cols) == 12:
        return "All 96 wells, A1 to H12"
    if len(rows) == 8:
        return (
            f"Column{'s' if len(cols) > 1 else ''} {_int_ranges(cols)}, "
            f"{rows[0]}{cols[0]} to {rows[-1]}{cols[-1]}"
        )
    if len(cols) == 12:
        return f"Row{'s' if len(rows) > 1 else ''} {_join(rows)}"
    return f"Row{'s' if len(rows) > 1 else ''} {_join(rows)}, column{'s' if len(cols) > 1 else ''} {_int_ranges(cols)}"


def read_sm2_design(path):
    """Return [(sample, [wells])]. Hard-fails on the 3-block file."""
    blocks = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            sample, spec = line.split("\t")
            blocks.append((sample.strip(), expand_wells(spec)))
    if len(blocks) != 2:
        sys.exit(
            f"{path} has {len(blocks)} sample blocks, expected 2.\n"
            "  This is almost certainly the three-block scifi-demux plate design, which carries a\n"
            "  third sample that is not part of this study. Pass the two-block design supplied to\n"
            "  AmbientMapper (maize and Arabidopsis blocks only)."
        )
    return blocks


def read_perwell_design(path):
    """Per-well split file. Returns the well set. The split is parallelism, not a design."""
    wells = set()
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            _, well = line.split("\t")
            wells.add(well.strip())
    return wells


def read_multi_wellmap(path):
    """Return {pool: {'genotype':g, 'wells':[...]}} from the well-to-genotype map."""
    pools = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            pool = row["genotype_rep"].strip()
            entry = pools.setdefault(pool, {"genotype": row["genotype"].strip(), "wells": []})
            entry["wells"].append(row["well"].strip())
    return pools


def design_supplied():
    """Read design_file out of every run's params JSON. Returns {library: {mode: bool}}."""
    out = {}
    for (library, mode), rel in DECONTAM_PARAMS.items():
        path = os.path.join(DECONTAM_ROOT, rel)
        if not os.path.exists(path):
            sys.exit(f"missing run parameters: {path}")
        with open(path) as fh:
            out.setdefault(library, {})[mode] = json.load(fh).get("design_file") is not None
    return out


def check(sm2_blocks, b73mo17_wells, multi_pools, supplied):
    problems = []

    want = EXPECTED["interspecies"]
    if len(sm2_blocks) != want["samples"]:
        problems.append(f"interspecies samples: got {len(sm2_blocks)}, expected {want['samples']}")
    for sample, wells in sm2_blocks:
        if len(wells) != want["wells_per_sample"]:
            problems.append(
                f"interspecies {sample} wells: got {len(wells)}, "
                f"expected {want['wells_per_sample']}"
            )
    total = sum(len(w) for _, w in sm2_blocks)
    if total != want["total"]:
        problems.append(f"interspecies total wells: got {total}, expected {want['total']}")

    for key in ("b73mo17_rep1", "b73mo17_rep2"):
        if len(b73mo17_wells[key]) != EXPECTED[key]["total"]:
            problems.append(
                f"{key} wells: got {len(b73mo17_wells[key])}, expected {EXPECTED[key]['total']}"
            )

    want = EXPECTED["multigenotype"]
    genotypes = {i["genotype"] for i in multi_pools.values()}
    if len(multi_pools) != want["pools"]:
        problems.append(f"multi-genotype pools: got {len(multi_pools)}, expected {want['pools']}")
    if len(genotypes) != want["genotypes"]:
        problems.append(
            f"multi-genotype genotypes: got {len(genotypes)}, expected {want['genotypes']}"
        )
    for pool, info in sorted(multi_pools.items()):
        if len(info["wells"]) != want["wells_per_pool"]:
            problems.append(
                f"multi-genotype pool {pool}: got {len(info['wells'])} wells, "
                f"expected {want['wells_per_pool']}"
            )
    total = sum(len(i["wells"]) for i in multi_pools.values())
    if total != want["total"]:
        problems.append(f"multi-genotype total wells: got {total}, expected {want['total']}")

    # The WD/ND contrast is the point of the table. Assert it rather than trust directory names.
    if not supplied["interspecies"]["WD"]:
        problems.append("interspecies WD run carries no design_file, the WD/ND contrast is broken")
    if supplied["interspecies"]["ND"]:
        problems.append("interspecies ND run carries a design_file, it is not design-free")
    for library in ("b73mo17_rep1", "b73mo17_rep2", "multigenotype"):
        if supplied[library]["ND"]:
            problems.append(f"{library} run carries a design_file, expected none")

    if problems:
        sys.exit("REGRESSION against the verified library design:\n  " + "\n  ".join(problems))


HEADER = [
    "Library", "Sample or genotype", "Plate index", "Plate-index pools", "Total Tn5 wells",
    "Plate index identifies the sample", "Plate design supplied to AmbientMapper",
]

CAPTION = (
    "Table S1. Library and plate design of the scifi-ATAC-seq libraries analysed in this study. "
    "Each well of the 96-well plate carries a unique pair of Tn5 plate indexes, so the plate index "
    "recorded in a barcode identifies the well in which that nucleus was tagmented. Where a well "
    "holds material from a single sample, the plate index therefore also identifies the sample, "
    "and this map is what AmbientMapper consumes in with-design (WD) mode. In the maize and "
    "Arabidopsis library each species was tagmented across its own set of 32 wells, and the map "
    "was supplied to the WD run and withheld from the ND run, which is the comparison used "
    "throughout the paper. In the B73/Mo17 libraries both genotypes were pooled before "
    "tagmentation and the whole plate holds a single sample, so no plate index distinguishes "
    "them and only design-free assignment is possible. Those libraries were demultiplexed per "
    "well to parallelise mapping, which is a processing convenience rather than an experimental "
    "design, and every downstream analysis treats the plate as one pooled sample. In the "
    "multi-genotype library each genotype occupies its own plate rows, but that map was not "
    "supplied to any AmbientMapper run and was used only as independent ground truth to score "
    "the resulting calls. The maize root library is a 10x scATAC library with no combinatorial "
    "plate and is described in Table S2. Plate positions and well counts are generated directly "
    "from the design files used by the analysis pipeline."
)


def build_rows(sm2_blocks, b73mo17_wells, multi_pools):
    rows = []

    for sample, wells in sm2_blocks:
        rows.append([
            LIBRARY_LABEL["interspecies"], SM2_SAMPLE_LABEL.get(sample, sample),
            describe_wells(wells), 1, len(wells), "Yes",
            "Yes for the WD run, withheld for the ND run",
        ])
    rows.append([
        LIBRARY_LABEL["interspecies"], "Library total", "", len(sm2_blocks),
        sum(len(w) for _, w in sm2_blocks), "", "",
    ])

    for key in ("b73mo17_rep1", "b73mo17_rep2"):
        wells = b73mo17_wells[key]
        rows.append([
            LIBRARY_LABEL[key], "B73 and Mo17, pooled before tagmentation",
            describe_wells(wells), 1, len(wells),
            "No, one pooled sample across the whole plate", "No",
        ])

    genotypes = {}
    for pool, info in multi_pools.items():
        entry = genotypes.setdefault(info["genotype"], {"pools": 0, "wells": []})
        entry["pools"] += 1
        entry["wells"].extend(info["wells"])
    order = ["B73"] + sorted(g for g in genotypes if g != "B73")
    for genotype in order:
        info = genotypes[genotype]
        rows.append([
            LIBRARY_LABEL["multigenotype"], genotype, describe_wells(info["wells"]),
            info["pools"], len(info["wells"]), "Yes",
            "No, used only as evaluation ground truth",
        ])
    rows.append([
        LIBRARY_LABEL["multigenotype"], "Library total", "", len(multi_pools),
        sum(len(i["wells"]) for i in genotypes.values()), "", "",
    ])

    return rows


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
    ws.title = "Table S1"
    ws.append(HEADER)
    for cell in ws[1]:
        cell.font = Font(bold=True)
        cell.alignment = Alignment(wrap_text=True, vertical="bottom")
    for row in rows:
        ws.append(row)
    for row in ws.iter_rows(min_row=2):
        if str(row[1].value) == "Library total":
            for cell in row:
                cell.font = Font(bold=True)
        for cell in row[1:]:
            cell.alignment = Alignment(wrap_text=True, vertical="top")
    for column, width in zip("ABCDEFG", (24, 34, 24, 15, 15, 24, 30)):
        ws.column_dimensions[column].width = width
    ws.row_dimensions[1].height = 44

    caption_row = ws.max_row + 2
    ws.cell(row=caption_row, column=1, value=CAPTION).alignment = Alignment(
        wrap_text=True, vertical="top"
    )
    ws.merge_cells(
        start_row=caption_row, start_column=1, end_row=caption_row + 11, end_column=7
    )
    wb.save(path)
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    for path in (SM2_DESIGN, MULTI_WELLMAP):
        if not os.path.exists(path):
            sys.exit(f"missing input: {path}")

    sm2_blocks = read_sm2_design(SM2_DESIGN)
    b73mo17_wells = {}
    for rep in (1, 2):
        path = B73MO17_DESIGN.format(rep=rep)
        if not os.path.exists(path):
            sys.exit(f"missing input: {path}")
        b73mo17_wells[f"b73mo17_rep{rep}"] = read_perwell_design(path)
    multi_pools = read_multi_wellmap(MULTI_WELLMAP)
    supplied = design_supplied()

    check(sm2_blocks, b73mo17_wells, multi_pools, supplied)
    print("self-check against the verified library design: PASS")

    rows = build_rows(sm2_blocks, b73mo17_wells, multi_pools)

    os.makedirs(args.outdir, exist_ok=True)
    txt = os.path.join(args.outdir, "TableS1_library_and_plate_design.txt")
    xlsx = os.path.join(args.outdir, "TableS1_library_and_plate_design.xlsx")
    write_txt(txt, rows)
    print(f"wrote {txt}")
    if write_xlsx(xlsx, rows):
        print(f"wrote {xlsx}")

    w0 = max(len(str(r[0])) for r in rows)
    w1 = max(len(str(r[1])) for r in rows)
    w2 = max(len(str(r[2])) for r in rows)
    for row in rows:
        print(f"{row[0]:<{w0}}  {row[1]:<{w1}}  {row[2]:<{w2}}  "
              f"pools={row[3]:>2}  wells={row[4]:>3}")


if __name__ == "__main__":
    main()
