"""4_1a_extract_at_markers.py -- marker prep, step 1 (Arabidopsis).

Extracts the seedling-relevant marker sheets from the Ecker et al. 2025 Arabidopsis atlas
supplementary table (saved in the working directory as spr_MOESM3.xlsx; stored in
6_socrates/_data/markers/ as Ecker2025_NatPlants_AtAtlas_SupplTable1.xlsx) into a long TSV
(geneID, name, organ, cell_type, sheet) written as at_markers_raw.tsv, the input of 4_1b.
Run from 6_socrates/_data/markers/:  python3 4_1a_extract_at_markers.py
"""
import openpyxl, re
from collections import Counter

wb = openpyxl.load_workbook("spr_MOESM3.xlsx", read_only=True, data_only=True)
# seedling-relevant sheets (exclude Flower/Silique/MERFISH-probe sheets)
keep = ['Epidermal Markers','Guard Cell Markers','Trichome Markers','Mesophyll Markers',
        'Vasculature Markers','SAM Markers','Root Cell Type Markers','Other Cell Type Markers']
agi = re.compile(r'^AT[1-5CM]G\d{5}$', re.I)

rows = []
for ws in wb.worksheets:
    if ws.title not in keep:
        continue
    data = list(ws.iter_rows(values_only=True))
    hdr = next((i for i, r in enumerate(data[:6])
                if r and any(isinstance(c, str) and c.strip() == "Gene ID" for c in r)), None)
    if hdr is None:
        continue
    H = [str(c).strip().lower() if c is not None else "" for c in data[hdr]]
    def ci(n): return H.index(n) if n in H else None
    gi, ni, oi, ti = ci("gene id"), ci("name"), ci("organ"), ci("cell type")
    for r in data[hdr+1:]:
        if not r or gi is None or r[gi] is None:
            continue
        gid = str(r[gi]).strip().upper()
        if not agi.match(gid):
            continue
        name  = str(r[ni]).strip() if ni is not None and r[ni] else gid
        organ = str(r[oi]).strip() if oi is not None and r[oi] else "NA"
        ctype = str(r[ti]).strip() if ti is not None and r[ti] else "NA"
        rows.append((gid, name, organ, ctype, ws.title))

with open("at_markers_raw.tsv", "w") as f:
    f.write("geneID\tname\torgan\tcell_type\tsheet\n")
    for r in rows:
        f.write("\t".join(r) + "\n")

print("total marker rows:", len(rows), "| unique genes:", len(set(r[0] for r in rows)))
print("\nby sheet:")
for s, n in Counter(r[4] for r in rows).most_common():
    print(f"  {n:3d}  {s}")
print("\nby organ:")
for o, n in Counter(r[2] for r in rows).most_common():
    print(f"  {n:3d}  {o}")
print("\ncell types (top 45):")
for ct, n in Counter(r[3] for r in rows).most_common(45):
    print(f"  {n:3d}  {ct}")
