#!/usr/bin/env python3
"""Build an INFORMATIVE-PRUNED marker panel from a 4_3d informativeness table.

Rebalances the uncapped augmented panel for the cluster ANNOTATION (4_3): the mean-z call needs
comparable per-type marker counts, so keep each cell type's top-N markers by empirical in-data
specificity (max_z across clusters, from 4_3d). Freeze on the PRE object (superset) and apply the
same pruned panel to BOTH stages so Pre/Post annotation stays comparable. Informativeness itself
stays computed on the full uncapped set -- this only trims what the ANNOTATION scores over.

Usage:
  4_3e_build_informative_panel.py <informativeness.tsv> <augmented_panel.bed> <out.bed> [N=15]
    informativeness.tsv : <PreName>.marker_informativeness.tsv (from 4_3d) -- needs geneID,type,max_z
    augmented_panel.bed : the full augmented bed (source of chr/start/end/name/type + coords)
"""
import sys
from collections import defaultdict

infp, bedp, outp = sys.argv[1], sys.argv[2], sys.argv[3]
N = int(sys.argv[4]) if len(sys.argv) > 4 else 15

# top-N geneIDs per cell type by max_z
rows = []
with open(infp) as f:
    h = f.readline().rstrip("\n").split("\t"); ix = {c: i for i, c in enumerate(h)}
    for line in f:
        p = line.rstrip("\n").split("\t")
        rows.append((p[ix["type"]], p[ix["geneID"]], float(p[ix["max_z"]])))
by_type = defaultdict(list)
for t, g, z in rows:
    by_type[t].append((g, z))
keep, per_type = set(), {}
for t, lst in by_type.items():
    lst.sort(key=lambda x: -x[1])
    sel = lst[:N]
    per_type[t] = len(sel)
    for g, _ in sel:
        keep.add(g)

# subset the augmented bed (preserves header, coords, ordering)
n_out = 0
with open(bedp) as f, open(outp, "w") as o:
    o.write(f.readline())                       # header
    for line in f:
        if line.split("\t")[3] in keep:
            o.write(line); n_out += 1

print("kept %d unique genes -> %d bed rows (<=%d/type, %d types) -> %s"
      % (len(keep), n_out, N, len(per_type), outp))
