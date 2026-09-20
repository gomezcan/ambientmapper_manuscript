"""4_1b_curate_at_markers.py -- marker prep, step 2 (Arabidopsis).

Collapses the fine cell-type labels of at_markers_raw.tsv (4_1a) onto the canonical seedling
cell-type vocabulary, ranks candidates per type (canonical, then named, then sheet order) and
writes the curated panel at_markers.curated.bed (chr start end geneID name type), which becomes
markers.At.Ecker2025.bed in 6_socrates/_data/markers/ and feeds 4_1c and the 4_3* engines.
Run from 6_socrates/_data/markers/:  python3 4_1b_curate_at_markers.py
"""
import re
from collections import defaultdict, Counter

rows = []
with open("at_markers_raw.tsv") as f:
    next(f)
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) == 5:
            rows.append(tuple(p))   # geneID, name, organ, cell_type, sheet

agi = re.compile(r'^AT[1-5CM]G\d{5}$', re.I)

def collapse(ct):
    s = re.sub(r'^(Meristem|Elongation|Maturation|Differentiation)[_ ]', '', ct.strip(), flags=re.I)
    s = re.sub(r'^(Proximal|Distal)[_ ]', '', s, flags=re.I)
    low = s.lower().strip()
    rules = [   # order matters; atrichoblast BEFORE trichoblast (substring), specific before generic
        ('quiescent center', 'QC'), ('quiescent centre', 'QC'),
        ('atrichoblast', 'atrichoblast'), ('trichoblast', 'trichoblast'),
        ('columella', 'columella'),
        ('lateral root cap', 'lateral_root_cap'), ('root cap', 'lateral_root_cap'),
        ('cortex', 'cortex'), ('endodermis', 'endodermis'), ('pericycle', 'pericycle'),
        ('metaphloem & companion', 'phloem_companion'), ('companion', 'phloem_companion'),
        ('protophloem', 'phloem'), ('metaphloem', 'phloem'), ('phloem', 'phloem'),
        ('protoxylem', 'xylem'), ('metaxylem', 'xylem'), ('xylem', 'xylem'),
        ('procambr', 'procambium'),                       # catches procambium + 'procambrium' typo
        ('abaxial epiderm', 'epidermis'), ('adaxial epiderm', 'epidermis'),
        ('protoderm', 'epidermis'), ('epiderm', 'epidermis'),
        ('stomat', 'guard_cell'), ('guard', 'guard_cell'),
        ('mesophyll', 'mesophyll'), ('bundle sheath', 'bundle_sheath'),
        ('trichome', 'trichome'), ('stele', 'stele'),
        ('dividing', 'dividing'), ('g1/s', 'dividing'), ('cell cycle', 'dividing'),
        ('sam', 'sam'), ('shoot apical', 'sam'), ('meristem', 'sam'),
    ]
    for k, v in rules:
        if k in low:
            return v
    if low == 'qc':
        return 'QC'
    return low.replace(' ', '_').replace('&', 'and')

KEEP = {'epidermis','guard_cell','trichome','mesophyll','bundle_sheath','sam','dividing',
        'procambium','phloem','phloem_companion','xylem',
        'trichoblast','atrichoblast','cortex','endodermis','pericycle','QC','columella',
        'lateral_root_cap','stele'}

bytype = defaultdict(list)
for gid, name, organ, ct, sheet in rows:
    mtype = collapse(ct)
    if mtype not in KEEP:
        continue
    is_named = ((not agi.match(name)) and (name.upper() != gid.upper()) and (1 <= len(name) <= 25)
                and not re.match(r'^\d{4}-\d{2}-\d{2}', name))   # Excel turned some gene symbols into dates
    bytype[mtype].append((gid, name, is_named))

# canonical (textbook) markers per cell type -> boosted to the front, then fill with atlas named markers
CANON = {
 'QC': ['WOX5','AGL42'],
 'endodermis': ['SCR','SCARECROW','MYB36','CASP1','ESB1','SGN3','MYB68'],
 'stele': ['SHR','SHORTROOT','WOL'],
 'trichoblast': ['COBL9','EXPA7','EXP7','RHD2','RSL4','LRX1','COW1','MRH1'],
 'atrichoblast': ['GL2','GLABRA2','WER','WEREWOLF','CPC'],
 'cortex': ['CO2','CORTEX'],
 'pericycle': ['LBD16','LBD29','GATA23','SKP2B','PUCHI'],
 'columella': ['PIN3','SMB','FEZ','DRO1','LZY','ARF'],
 'lateral_root_cap': ['SMB','BRN1','BRN2','SOMBRERO'],
 'procambium': ['ATHB8','PXY','MP','MONOPTEROS','TMO5','TMO6','WOL','TDR'],
 'phloem': ['APL','NEN4','SMXL5','BRX','OPS','CALS7','SEOR1','NAC45','NAC86'],
 'phloem_companion': ['SUC2','SUCROSE','APL','AHA3'],
 'xylem': ['VND6','VND7','VND','XCP1','XCP2','IRX3','IRX5','MYB46','MYB83','CESA7'],
 'guard_cell': ['FAMA','MUTE','SPCH','KAT1','SLAC1','MYB60','HIC','FLP','GORK'],
 'mesophyll': ['CAB','CAB1','CAB3','RBCS','LHCB','LHCA','CA1','CORI3'],
 'epidermis': ['ATML1','FDH','PDF1','PDF2','DCR','ML1'],
 'trichome': ['GL1','GLABRA1','GL3','TTG1','TRY','ETC2','MYB23','CPC'],
 'sam': ['CLV3','WUS','STM','KNAT1','CUC1','CUC2','SHOOTMERISTEMLESS'],
 'bundle_sheath': ['SCL23','SULTR2','MYB76'],
 'dividing': ['CDKA','CDKB','CYCB','CYCA','PCNA','HISTONE','KNOLLE','AURORA','MCM'],
}
def _norm(s): return re.sub(r'[^A-Z0-9]', '', str(s).upper())
def canon_rank(mtype, name):
    nm = _norm(name); toks = set(re.split(r'[\s/;,_-]+', str(name).upper())) | {nm}
    for i, p in enumerate(CANON.get(mtype, [])):
        if _norm(p) == nm or p.upper() in toks:
            return i
    return 999

CAP = 10
out = []
for mtype, lst in bytype.items():
    uniq = {}                                            # gid -> (name, is_named, sheet_order)
    for order, (gid, name, isn) in enumerate(lst):
        if gid not in uniq or (isn and not uniq[gid][1]):
            uniq[gid] = (name, isn, order)
    cands = [(gid, nm, isn, od) for gid, (nm, isn, od) in uniq.items()]
    cands.sort(key=lambda x: (canon_rank(mtype, x[1]), not x[2], x[3]))  # canonical, then named, then sheet order
    for gid, nm, isn, od in cands[:CAP]:
        out.append((gid, nm if isn else gid, mtype))

out.sort(key=lambda x: (x[2], x[1]))
with open("at_markers.curated.bed", "w") as f:
    f.write("chr\tstart\tend\tgeneID\tname\ttype\n")
    for gid, name, mtype in out:
        f.write(f"At\t0\t0\t{gid}\t{name}\t{mtype}\n")

print("curated At markers:", len(out), "| unique genes:", len(set(o[0] for o in out)),
      "| cell types:", len(set(o[2] for o in out)))
print("\nby major cell type (named markers shown):")
for ct in sorted(bytype, key=lambda c: -len([o for o in out if o[2] == c])):
    genes = [o[1] for o in out if o[2] == ct]
    print(f"  {len(genes):2d}  {ct:18s} {', '.join(genes[:8])}{'...' if len(genes) > 8 else ''}")
