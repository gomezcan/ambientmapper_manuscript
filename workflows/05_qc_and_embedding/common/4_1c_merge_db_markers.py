#!/usr/bin/env python3
"""Merge PlantscRNAdb 4.0 high-confidence markers into the curated Ecker (At) / Marand (maize)
panels, producing the SM2v2.v2 marker set for the Step-4 per-genome annotation.

Design:
  * NO per-cell-type cap by default -- keep the FULL high-confidence candidate set and let the
    4_3 marker-accessibility z-scores rank markers empirically in the SM2v2 data ("test all";
    optional positional arg sets a cap for a trimmed variant).
  * Keep ALL canonical Marand/Ecker markers (a cap only ever ADDS db genes, never drops canonical).
  * Each gene is assigned a single primary cell type = the mapped canonical type where it attains
    its highest COSG specificity score (the 4_3 loader dedups by geneID, so one type per gene).

Sources (all in _data/markers/):
  * Gene_marker.txt              -- PlantscRNAdb 4.0 full download (33 species; hcmarker=="Yes" used)
  * markers.At.Ecker2025.bed     -- canonical At panel (TAIR AGI)
  * markers.maize.Marand2025.bed -- canonical maize panel (AGPv5 Zm00001eb)
  * B73v4_to_B73v5.tsv           -- MaizeGDB Zm00001d->Zm00001eb lift (1:1 used for the panel)
  * maize_v5_gene_coords.tsv     -- v5 gene coords for db-added maize rows (geneID chr start end)

ID compatibility: At markers are TAIR AGI (no lift). Maize PlantscRNAdb markers are v4 (Zm00001d)
and are lifted v4->v5; only unambiguous 1:1 lifts enter the panel (tandem/paralog 1:many skipped).

Usage:  python3 0_scripts/common/4_1c_merge_db_markers.py [CAP]
        CAP omitted -> no cap (full set). CAP=N -> at most N markers/type (canonical kept regardless).
"""
import sys, os
from collections import defaultdict

BASE = os.path.join(os.environ.get("PROJECT_ROOT", "."), "6_socrates")
MARK = BASE + "/_data/markers"
ECKER  = MARK + "/markers.At.Ecker2025.bed"
MARAND = MARK + "/markers.maize.Marand2025.bed"
PSRDB  = MARK + "/Gene_marker.txt"
XREF   = MARK + "/B73v4_to_B73v5.tsv"
COORDS = MARK + "/maize_v5_gene_coords.tsv"

CAP = int(sys.argv[1]) if len(sys.argv) > 1 else None
cap_eff = CAP if CAP else 10**9

# ---- PlantscRNAdb fine-type -> canonical vocab (None = DROP: not present in SM2v2 seedling) ----
AT_MAP = {
 "Bundle sheath cell":"bundle_sheath","Columella":"columella","Columella initial cell":"columella",
 "Columella root cap cell":"columella","Companion cell":"phloem_companion","Cortex":"cortex",
 "Dividing protoderm":"dividing","Early anther":None,"Early sieve element":"phloem",
 "Epidermal cell layer of the emerging lateral root":"epidermis","Epidermis":"epidermis",
 "Explant vasculature and callus founder cell":None,"Floral organ abscission zone":None,
 "Flower meristem":None,"G2/M-phase cell":"dividing","Generative nuclei":None,"Guard cell":"guard_cell",
 "Guard mother cell":"guard_cell","Hydathode":None,"Hypocotyl":None,"Initials":None,
 "Intermediate anther":None,"Late anther":None,"Lateral root cap":"lateral_root_cap",
 "Lateral root cap like cell":"lateral_root_cap","Leaf abaxial pavement cell":"epidermis",
 "Leaf adaxial pavement cell":"epidermis","Leaf epidermis":"epidermis","Leaf guard cell":"guard_cell",
 "Leaf lamina epidermis":"epidermis","Leaf pavement cell":"epidermis","Leaf preprocambium":"procambium",
 "Leaf vascular system":None,"Mesophyll cell":"mesophyll","Metaphloem":"phloem",
 "Metaphloem sieve element":"phloem","Metaxylem":"xylem","Microspore nuclei":None,"Myrosinase cell":None,
 "Non-hair root epidermal cell":"atrichoblast","Pericycle cell":"pericycle","Phloem":"phloem",
 "Phloem parenchyma cell":"phloem","Phloem pole pericycle cell":"pericycle","Pholem":"phloem",
 "Photosynthetic cell":"mesophyll","Plant embryo cotyledon":None,"Primary phloem":"phloem",
 "Primary xylem":"xylem","Procambium":"procambium","Proliferating cell":"dividing","Protophloem":"phloem",
 "Protophloem sieve element":"phloem","Protoxylem":"xylem","Quiescent center":"QC","Radicle":None,
 "Root cortex":"cortex","Root endodermis":"endodermis","Root hair cell":"trichoblast",
 "Root initial cell":None,"Root meristem":None,"Root procambium":"procambium","S-phase cell":"dividing",
 "Shoot apical meristem":"sam","Shoot meristem initial cell":"sam","Shoot parenchyma":"cortex",
 "Shoot system cortex":"cortex","Shoot system endodermis":"endodermis","Shoot system epidermis":"epidermis",
 "Shoot system vascular system":None,"Sieve element":"phloem","Sperm nuclei":None,
 "Spongy mesophyll cell":"mesophyll","Stomatal initial cell":"guard_cell","Suspensor":None,
 "Tracheary element":"xylem","Vascular initial cell":"procambium","Vegetative nuclei":None,"Xylem":"xylem",
 "Xylem parenchyma cell":"xylem","Xylem pole pericycle cell":"pericycle","Young guard cell":"guard_cell",
}
# maize: confident shoot/leaf correspondences only; root + generic-vascular + ambiguous-SAM + ear = DROP
MAIZE_MAP = {
 "Bundle sheath cell":"bundle_sheath","Mesophyll cell":"mesophyll",
 "Leaf pavement cell":"developing_pavement_cell","Subsidiary cell":"subsidiary_cell",
 "Leaf lamina epidermis":"epidermis","Meristem epidermis layer":"epidermis","Xylem":"xylem",
 "Companion cell":"companion_cells","Meristemoid":"stomatal_precursor",
 "Root stele":None,"Root epidermis":None,"Mature stele":None,"Phloem":None,"Root cortex":None,
 "Root endodermis":None,"Early stele":None,"Pericycle cell":None,"Root pith":None,"Mature endodermis":None,
 "Ear cortex":None,"Meristem boundary":None,"Leaf vascular system":None,"Meristem base":None,
 "Determinate lateral organ":None,"Adaxial meristem periphery":None,
}

def load_bed(path):
    rows = []
    with open(path) as f:
        f.readline()
        for line in f:
            line = line.rstrip("\n").rstrip("\r")
            if line:
                rows.append(line.split("\t"))
    return rows

def load_psrdb():
    at, mz = [], []
    with open(PSRDB) as f:
        f.readline()
        for line in f:
            line = line.rstrip("\n").rstrip("\r")
            if not line: continue
            p = [c.strip('"') for c in line.split("\t")]
            if len(p) < 16 or p[15] != "Yes": continue          # high-confidence only
            gene, ct, sp = p[4], p[3], p[5]
            try: cosg = float(p[9])
            except: cosg = 0.0
            if sp == "Arabidopsis thaliana": at.append((ct, gene, cosg))
            elif sp == "Zea mays":           mz.append((ct, gene, cosg))
    return at, mz

def load_lift():
    m = {}
    with open(XREF) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 2: continue
            tgt = [x for x in p[1].split(",") if x.startswith("Zm00001eb")]
            if len(tgt) == 1:                                    # unambiguous 1:1 only
                m[p[0]] = tgt[0]
    return m

def load_v5_coords():
    c = {}
    if not os.path.exists(COORDS): return c
    with open(COORDS) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 4:
                c[p[0]] = (p[1], p[2], p[3])
    return c

def build_db_by_type(records, typemap, id_ok, lift=None):
    best = {}                                                    # gene -> (cosg, canon_type)
    unmapped = set()
    for ct, gene, cosg in records:
        if lift is not None:
            if not gene.startswith("Zm00001d") or gene not in lift: continue
            gene = lift[gene]
        if not id_ok(gene): continue
        if ct not in typemap:
            unmapped.add(ct); continue
        canon = typemap[ct]
        if canon is None: continue
        if gene not in best or cosg > best[gene][0]:
            best[gene] = (cosg, canon)
    by_type = defaultdict(list)
    for gene, (cosg, canon) in best.items():
        by_type[canon].append((gene, cosg))
    for t in by_type:
        by_type[t].sort(key=lambda x: -x[1])
    return by_type, unmapped

def merge(canon_rows, by_type, cap, chrom_default, coords=None):
    canon_by_type = defaultdict(list); canon_genes = set()
    for r in canon_rows:
        canon_by_type[r[5]].append(r); canon_genes.add(r[3])
    out, stats = [], {}
    for typ in sorted(set(list(canon_by_type) + list(by_type))):
        kept = [(r, "canonical") for r in canon_by_type.get(typ, [])]   # never drop canonical
        n_can = len(kept)
        for gene, cosg in by_type.get(typ, []):
            if len(kept) >= cap: break
            if gene in canon_genes: continue
            ch, s, e = coords[gene] if (coords and gene in coords) else (chrom_default, "0", "0")
            kept.append(([ch, s, e, gene, gene, typ], "PlantscRNAdb4"))
        stats[typ] = (n_can, len(kept) - n_can, len(kept))
        out.extend(kept)
    return out

def write_bed(path, rows):
    with open(path, "w") as f:
        f.write("chr\tstart\tend\tgeneID\tname\ttype\n")
        for r, _ in rows:
            f.write("\t".join(r[:6]) + "\n")

def dump_map(path, m, target_col):
    with open(path, "w") as f:
        f.write("plantscrnadb_type\t%s\n" % target_col)
        for k in sorted(m):
            f.write("%s\t%s\n" % (k, m[k] if m[k] else "DROP"))

def summarize(tag, rows):
    by = defaultdict(lambda: [0, 0])
    for r, src in rows:
        by[r[5]][0 if src == "canonical" else 1] += 1
    lines = ["## %s   total=%d  (canonical=%d  +PlantscRNAdb=%d)  types=%d"
             % (tag, len(rows), sum(v[0] for v in by.values()),
                sum(v[1] for v in by.values()), len(by)),
             "%-32s %6s %6s %6s" % ("type", "canon", "+db", "tot")]
    for t in sorted(by, key=lambda x: -(by[x][0] + by[x][1])):
        c, d = by[t]
        lines.append("%-32s %6d %6d %6d" % (t, c, d, c + d))
    return "\n".join(lines)

# ---------------- run ----------------
ecker, marand = load_bed(ECKER), load_bed(MARAND)
at_rec, mz_rec = load_psrdb()
lift, coords = load_lift(), load_v5_coords()
is_agi = lambda g: len(g) >= 8 and g[:2] == "AT" and g[2] in "12345CM" and g[3] == "G"
is_v5  = lambda g: g.startswith("Zm00001eb")

at_by, at_un = build_db_by_type(at_rec, AT_MAP, is_agi)
mz_by, mz_un = build_db_by_type(mz_rec, MAIZE_MAP, is_v5, lift=lift)
at_rows = merge(ecker,  at_by, cap_eff, "At", coords=None)
mz_rows = merge(marand, mz_by, cap_eff, "NA", coords=coords)

write_bed(MARK + "/markers.At.Ecker2025_PlantscRNAdb4.bed", at_rows)
write_bed(MARK + "/markers.maize.Marand2025_PlantscRNAdb4.bed", mz_rows)
write_bed(MARK + "/markers.SM2v2.v2.bed", at_rows + mz_rows)
dump_map(MARK + "/celltype_map.At.PlantscRNAdb_to_Ecker.tsv", AT_MAP, "ecker_type")
dump_map(MARK + "/celltype_map.maize.PlantscRNAdb_to_Marand.tsv", MAIZE_MAP, "marand_type")

report = "\n\n".join([
    "SM2v2.v2 marker panel build report   (cap=%s)" % ("none" if CAP is None else CAP),
    summarize("Arabidopsis (TAIR10)  markers.At.Ecker2025_PlantscRNAdb4.bed", at_rows),
    summarize("Maize (B73v5)  markers.maize.Marand2025_PlantscRNAdb4.bed", mz_rows),
    "At fine-types skipped (unmapped): %s"  % (sorted(at_un) or "none"),
    "Maize fine-types skipped (unmapped): %s" % (sorted(mz_un) or "none"),
])
open(MARK + "/markers.SM2v2.v2.build_report.txt", "w").write(report + "\n")
print(report)
print("\nWROTE -> _data/markers/{markers.At.Ecker2025_PlantscRNAdb4.bed, "
      "markers.maize.Marand2025_PlantscRNAdb4.bed, markers.SM2v2.v2.bed, "
      "celltype_map.*.tsv, markers.SM2v2.v2.build_report.txt}")
