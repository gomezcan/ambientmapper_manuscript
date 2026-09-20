#!/usr/bin/env python3
"""
06_48_extract_barcode_block_counts.py

Barcode-resolved, 1-Mb-block-aggregated allele match counts for the
depth- & call-class-stratified purity model (06_48). Extends the 06_44
count engine from (group x site) to (barcode x block), folding founder
allele-matching into the Python pass (a per-SNP barcode table would explode;
aggregating to 1-Mb blocks keeps LD structure for the block bootstrap while
staying compact).

WHY THIS EXISTS
  06_44 emits one pooled (group x site) table -> it cannot see barcode depth
  or singlet/doublet/ambiguous status. 06_48 keeps the barcode as the unit so
  we can (a) stratify purity by AM call class x read depth, (b) test whether
  weak_doublets are genetically-clean singlets (top1 high) or real doublets
  (top1~0.5, top1or2~1.0), and (c) decompose the raw->clean gain into
  within-barcode cleaning vs read-reweighting.

============================ DEFINITIONS ======================================
Sites: KNOWN founder panel sites (reference VCF w/ per-founder GT). A site is
  "informative" iff >=1 founder is hom-ref AND >=1 is hom-alt (discriminates
  some founder pair) -- identical to 06_44. No MAF/LD/missingness filters.
Per barcode we know (from the AM cells_calls attr map): call_class, top1
  (=genome_1, the ASSIGNED genotype), top2 (=genome_2, "None" for singlets),
  raw_depth (=n_reads). Optional plate truth (multi) via --well-map.

Emitted per (barcode, 1-Mb block), for NON-ambiguous barcodes
  (single_clean / doublet / weak_doublet):
    t_top1, m_top1  : over {top1 homozygous & informative} sites,
                      t_top1 = panel-allele reads (ref+alt),
                      m_top1 = reads carrying top1's allele.
                      -> top1 self-purity = Sum m_top1 / Sum t_top1
                         (the 06_44-comparable metric; THE singlet/doublet
                          discriminator).
    t_top2, m_top2  : same, anchored to top2 (0 if no top2).
    p_tot, p_top1, p_top2 : over {top1 hom & top2 hom & top1!=top2 &
                      informative} = pair-discriminating sites,
                      p_tot   = all reads (ref+alt+other),
                      p_top1  = reads == top1 allele,
                      p_top2  = reads == top2 allele.
                      -> top1or2 (pair-consistency) = (p_top1+p_top2)/p_tot,
                         top1 pair-fraction         = p_top1/(p_top1+p_top2).
    t_truth, m_truth: plate-truth anchored (multi only; 0 otherwise).

Emitted per (barcode, block) for AMBIGUOUS barcodes, in a SEPARATE sidecar
  (they get a margin diagnostic, not a purity line): per founder F,
    t_F, m_F over {F hom & informative} sites -> concordance_F = Sum m_F/Sum t_F;
    R then reports best / second-best / margin = best - second.

CHOICES MADE:
  C1. Anchor = AM call (top1=genome_1) for ALL datasets, uniform. Plate truth
      is an ADDITIONAL column for multi, not the primary anchor. (The
      weak_doublet question is about the assigned call, so top1 is the right
      primary anchor; truth stays available for the AM-error-rate view.)
  C2. top1 self-purity denominator = panel reads at {top1-hom & informative}
      sites (ref+alt, excludes 'other'), exactly 06_44 -> directly comparable
      to the WASP concordance numbers.
  C3. pair metrics use pair-DISCRIMINATING sites only (both hom, different
      alleles) so top1or2 / pair-fraction are well-defined; p_tot includes
      'other' so top1or2 can drop below 1 for junk barcodes.
  C4. Barcode set = whatever is in --bc-attr-map (built upstream with
      n_reads>=200 & genome_1 in founders); 'empty'/'low_reads' excluded there.

OUTPUTS
  <out>.counts.tsv.gz    barcode block t_top1 m_top1 t_top2 m_top2 p_tot p_top1 p_top2 t_truth m_truth
  <out>.ambiguous.tsv.gz barcode block <F1>_t <F1>_m ... <FK>_t <FK>_m
  <out>.barcode_attrs.tsv barcode call top1 top2 truth raw_depth
Parallel by chromosome (ProcessPoolExecutor; worker re-opens the BAM),
mirroring 06_44_extract_allele_counts.py.
"""
from __future__ import annotations

import argparse
import bisect
import gzip
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict, List, Optional, Tuple

import pysam

BLOCK = 1_000_000  # 1-Mb blocks (LD-aware; matches 06_44 bootstrap block size)


# ------------------------------------------------------------- panel loading
def load_panel(path: str, genos: List[str]):
    """Return {chrom: (sorted_pos[], {pos: site})} where site =
    (ref_char, alt_char, expected_base_per_founder{list}, informative_bool).
    expected_base_per_founder[k] = ref/alt char if founder k is homozygous,
    else None (het/missing). Panel TSV: CHROM POS REF ALT GT1..GTK (no header).
    """
    op = gzip.open if path.endswith(".gz") else open
    per_chrom: Dict[str, Dict[int, tuple]] = defaultdict(dict)
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4 + len(genos):
                continue
            chrom, pos, ref, alt = f[0], int(f[1]), f[2].upper(), f[3].upper()
            if len(ref) != 1:
                continue
            alt1 = alt.split(",")[0]
            if len(alt1) != 1 or alt1 == "*":
                continue
            exp: List[Optional[str]] = []
            n_ref_hom = n_alt_hom = 0
            for k in range(len(genos)):
                gt = f[4 + k].replace("|", "/")
                if gt in ("0/0",):
                    exp.append(ref); n_ref_hom += 1
                elif gt in ("1/1",):
                    exp.append(alt1); n_alt_hom += 1
                else:                       # het / missing / multiallelic-other
                    exp.append(None)
            informative = n_ref_hom > 0 and n_alt_hom > 0
            if not informative:
                continue                    # only keep informative sites
            per_chrom[chrom][pos] = (ref, alt1, exp, informative)
    out = {}
    for c, d in per_chrom.items():
        pos_sorted = sorted(d.keys())
        out[c] = (pos_sorted, d)
    return out


# --------------------------------------------------------- attr / well maps
def _open_text(p: str):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "rt")


def load_attr_map(path: str):
    """bc_key -> (call, top1, top2, raw_depth). TSV header:
    bc_key call top1 top2 raw_depth  (top2 'None'/'' -> None)."""
    out: Dict[str, tuple] = {}
    with _open_text(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(hdr)}
        for c in ("bc_key", "call", "top1", "top2", "raw_depth"):
            if c not in ix:
                sys.exit(f"[06_48] attr-map missing column '{c}'")
        ib, ic, i1, i2, ir = (ix["bc_key"], ix["call"], ix["top1"],
                              ix["top2"], ix["raw_depth"])
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= ir:
                continue
            top2 = f[i2].strip()
            top2 = None if top2 in ("", "None", "NA") else top2
            try:
                rd = int(float(f[ir]))
            except ValueError:
                rd = -1
            out[f[ib].strip()] = (f[ic].strip(), f[i1].strip(), top2, rd)
    return out


def load_well_map(path: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        if "well" not in hdr or "genotype" not in hdr:
            sys.exit("[06_48] well-map needs 'well' and 'genotype'")
        iw, ig = hdr.index("well"), hdr.index("genotype")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > max(iw, ig) and f[iw].strip() and f[ig].strip() not in ("", "NA"):
                out[f[iw].strip()] = f[ig].strip()
    return out


# --------------------------------------------------------------- worker
_W = {}


def _init(bam, bc_tag, min_mapq, attr, well, geno_ix):
    _W.update(bam=bam, bc_tag=bc_tag, min_mapq=min_mapq,
              attr=attr, well=well, geno_ix=geno_ix)


def _truth_for(bc: str) -> Optional[str]:
    well = _W["well"]
    if well is None:
        return None
    try:
        parts = bc.split("-", 1)[1].split("_")
        if len(parts) < 4:
            return None
        return well.get(parts[-3])
    except IndexError:
        return None


def _process_chrom(item):
    chrom, pos_sorted, site_map = item
    attr, geno_ix = _W["attr"], _W["geno_ix"]
    bc_tag, min_mapq = _W["bc_tag"], _W["min_mapq"]
    K = len(geno_ix)

    # main[(bc,block)] = [t_top1,m_top1,t_top2,m_top2,p_tot,p_top1,p_top2,t_truth,m_truth]
    main: Dict[Tuple[str, int], List[int]] = defaultdict(lambda: [0] * 9)
    # amb[(bc,block)] = [t_f0,m_f0, t_f1,m_f1, ...] length 2K
    amb: Dict[Tuple[str, int], List[int]] = defaultdict(lambda: [0] * (2 * K))
    truth_of: Dict[str, str] = {}

    bam = pysam.AlignmentFile(_W["bam"], "rb")
    try:
        if chrom not in bam.references:
            return chrom, {}, {}, {}
        for read in bam.fetch(chrom):
            if (read.is_unmapped or read.is_secondary or read.is_supplementary
                    or read.is_duplicate or read.mapping_quality < min_mapq):
                continue
            rs, re = read.reference_start, read.reference_end
            if re is None or not read.has_tag(bc_tag):
                continue
            lo = bisect.bisect_left(pos_sorted, rs + 1)
            hi = bisect.bisect_right(pos_sorted, re)
            if lo == hi:
                continue
            bc_full = read.get_tag(bc_tag)
            bc_key = bc_full.split("-", 1)[0]
            a = attr.get(bc_key)
            if a is None:
                continue
            call, top1, top2, _rd = a
            seq = read.query_sequence
            if seq is None:
                continue
            ref_to_q = {rp: qp for qp, rp in read.get_aligned_pairs(matches_only=True)}

            is_amb = (call == "ambiguous")
            if not is_amb:
                i1 = geno_ix.get(top1)
                i2 = geno_ix.get(top2) if top2 else None
                itr = geno_ix.get(_truth_for(bc_full))
                if itr is not None and bc_key not in truth_of:
                    truth_of[bc_key] = _truth_for(bc_full)

            for p in pos_sorted[lo:hi]:
                q = ref_to_q.get(p - 1)
                if q is None:
                    continue
                b = seq[q].upper()
                ref, alt1, exp, _inf = site_map[p]
                if b != ref and b != alt1:
                    is_panel = False
                else:
                    is_panel = True
                blk = p // BLOCK

                if is_amb:
                    rec = amb[(bc_key, blk)]
                    for k in range(K):
                        ek = exp[k]
                        if ek is None:
                            continue
                        if is_panel:
                            rec[2 * k] += 1            # t_F (panel read at F-hom site)
                            if b == ek:
                                rec[2 * k + 1] += 1    # m_F
                    continue

                rec = main[(bc_key, blk)]
                e1 = exp[i1] if i1 is not None else None
                e2 = exp[i2] if i2 is not None else None
                # top1 self
                if e1 is not None and is_panel:
                    rec[0] += 1
                    if b == e1:
                        rec[1] += 1
                # top2 self
                if e2 is not None and is_panel:
                    rec[2] += 1
                    if b == e2:
                        rec[3] += 1
                # pair-discriminating (both hom, different alleles)
                if e1 is not None and e2 is not None and e1 != e2:
                    rec[4] += 1                          # p_tot (incl 'other')
                    if b == e1:
                        rec[5] += 1
                    elif b == e2:
                        rec[6] += 1
                # truth
                if itr is not None:
                    etr = exp[itr]
                    if etr is not None and is_panel:
                        rec[7] += 1
                        if b == etr:
                            rec[8] += 1
    finally:
        bam.close()
    return chrom, dict(main), dict(amb), truth_of


# ----------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--genotype-panel", required=True,
                    help="TSV: CHROM POS REF ALT GT1..GTK (bcftools query, no header)")
    ap.add_argument("--genotypes", required=True, help="comma list, panel GT order")
    ap.add_argument("--bc-attr-map", required=True,
                    help="TSV: bc_key call top1 top2 raw_depth (from cells_calls)")
    ap.add_argument("--well-map", help="optional plate truth: well genotype (multi)")
    ap.add_argument("--bc-tag", default="BC")
    ap.add_argument("--min-mapq", type=int, default=10)
    ap.add_argument("--threads", type=int, default=10)
    ap.add_argument("--chroms", default="", help="comma list to keep (default: all in panel)")
    ap.add_argument("--out", required=True, help="output prefix (writes .counts/.ambiguous/.barcode_attrs)")
    args = ap.parse_args()

    genos = args.genotypes.split(",")
    geno_ix = {g: i for i, g in enumerate(genos)}
    print(f"[06_48] genotypes: {genos}", file=sys.stderr)

    panel = load_panel(args.genotype_panel, genos)
    if args.chroms.strip():
        keep = set(args.chroms.split(","))
        panel = {c: v for c, v in panel.items() if c in keep}
    n_sites = sum(len(v[0]) for v in panel.values())
    print(f"[06_48] {n_sites:,} informative panel sites across {len(panel)} chrom(s)", file=sys.stderr)

    attr = load_attr_map(args.bc_attr_map)
    well = load_well_map(args.well_map) if args.well_map else None
    print(f"[06_48] {len(attr):,} barcodes in attr map; well-map: "
          f"{'yes ('+str(len(well))+' wells)' if well else 'no'}", file=sys.stderr)

    work = [(c, panel[c][0], panel[c][1]) for c in sorted(panel)]
    print(f"[06_48] {len(work)} chrom(s) x {args.threads} workers", file=sys.stderr)

    cf = gzip.open(args.out + ".counts.tsv.gz", "wt")
    af = gzip.open(args.out + ".ambiguous.tsv.gz", "wt")
    cf.write("barcode\tblock\tt_top1\tm_top1\tt_top2\tm_top2\tp_tot\tp_top1\tp_top2\tt_truth\tm_truth\n")
    af.write("barcode\tblock\t" + "\t".join(f"{g}_t\t{g}_m" for g in genos) + "\n")

    truth_all: Dict[str, str] = {}
    n_main = n_amb = 0
    with ProcessPoolExecutor(max_workers=args.threads, initializer=_init,
                             initargs=(args.bam, args.bc_tag, args.min_mapq,
                                       attr, well, geno_ix)) as ex:
        futs = [ex.submit(_process_chrom, w) for w in work]
        for fut in as_completed(futs):
            chrom, main, amb, truth_of = fut.result()
            for (bc, blk), r in main.items():
                cf.write(f"{bc}\t{chrom}_{blk}\t" + "\t".join(map(str, r)) + "\n")
            for (bc, blk), r in amb.items():
                af.write(f"{bc}\t{chrom}_{blk}\t" + "\t".join(map(str, r)) + "\n")
            truth_all.update(truth_of)
            n_main += len(main); n_amb += len(amb)
            print(f"[06_48]   {chrom}: {len(main):,} main + {len(amb):,} amb (bc,block) rows",
                  file=sys.stderr)
    cf.close(); af.close()

    with open(args.out + ".barcode_attrs.tsv", "wt") as bf:
        bf.write("barcode\tcall\ttop1\ttop2\ttruth\traw_depth\n")
        for bc, (call, top1, top2, rd) in attr.items():
            bf.write(f"{bc}\t{call}\t{top1}\t{top2 or 'None'}\t{truth_all.get(bc,'NA')}\t{rd}\n")

    print(f"[06_48] DONE — {n_main:,} main + {n_amb:,} ambiguous (bc,block) rows; "
          f"{len(attr):,} barcode attrs -> {args.out}.*", file=sys.stderr)


if __name__ == "__main__":
    main()
