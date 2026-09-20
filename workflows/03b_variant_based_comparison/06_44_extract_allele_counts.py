#!/usr/bin/env python3
"""
06_44_extract_allele_counts.py

2d (count-based) ENGINE. Scan a BAM and AGGREGATE, per (pseudo-bulk group ×
reference-VCF SNP site), the number of reads carrying REF / ALT / other. Unlike
06_35 (per-read, reservoir-sampled diagnostic), this uses ALL reads genome-wide
and emits one compact group x site count table per BAM — the input to a
reference-anchored allele concordance / correlation estimate
(06_44_concordance.R).

WHY COUNTS, NOT GENOTYPES
  A pseudo-bulk is a POOL of cells, not a diploid individual: its identity
  signal is the allele FRACTION at known founder-discriminating sites, read
  straight off pileup counts. No SNP discovery, no MAF / LD / missingness
  filters — the sites are already known (the reference VCF), so we just count.
  This is the fix for the SNPRelate GRM keyhole: that pipeline needed a dense
  diploid genotype matrix and collapsed 17M known sites to ~3.5k.

GROUP ASSIGNMENT (per read, from the BAM BC tag):
  --group-mode bc_key_map : bc_key = BC.split('-')[0]  -> group via --bc-map TSV
                            (B73Mo17: bc_key -> AM genome_1 'top-1' call)
  --group-mode plate      : well = BC suffix parts[-3]  -> group via --well-map
                            (multiGenotypes: plate-of-origin, ground truth)

The SAME group map is applied to the raw and the clean BAM, so a raw-vs-clean
delta is a fixed-barcode-set, fixed-label comparison where only the reads differ.

OUTPUT (TSV.gz): chrom  pos  ref  alt  group  n_ref  n_alt  n_other
SIDECAR (TSV)  : {out%.gz}.bias.tsv  group  obs_class  n_obs  sum_nm  sum_mapq
                 (reference-mapping-bias readout: NM / MAPQ for ref- vs alt-obs)

Parallel by chromosome (ProcessPoolExecutor; worker re-opens the BAM).
Adapted from 06_35_extract_read_alleles.py.
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


# ---------------------------------------------------------------- VCF loading

def load_snp_positions(vcf_path: str) -> Dict[str, List[Tuple[int, str, str]]]:
    """Return {chrom: [(pos_1based, ref, alt1), ...]} sorted by pos. SNVs only.

    Only the first ALT (alt1) is retained, matching the 06_35 genotype-panel
    convention so the R join (ref / alt1 / founder GT) is consistent.
    """
    out: Dict[str, List[Tuple[int, str, str]]] = defaultdict(list)
    op = gzip.open if vcf_path.endswith(".gz") else open
    with op(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            chrom, pos, ref, alt = f[0], f[1], f[3], f[4]
            if len(ref) != 1:
                continue
            alts = [a for a in alt.split(",") if len(a) == 1 and a != "*"]
            if not alts:
                continue
            out[chrom].append((int(pos), ref.upper(), alts[0].upper()))
    for c in out:
        out[c].sort(key=lambda x: x[0])
    return dict(out)


# ------------------------------------------------------------- target loading

def _open_text(p: str):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "rt")


def load_bc_map(path: str) -> Dict[str, str]:
    """Return {bc_key: group} from a TSV with header columns 'bc_key' and 'group'."""
    out: Dict[str, str] = {}
    with _open_text(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(header)}
        if "bc_key" not in ix or "group" not in ix:
            sys.exit(f"[06_44] {path} needs columns 'bc_key' and 'group'")
        ib, ig = ix["bc_key"], ix["group"]
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(ib, ig):
                continue
            g = f[ig].strip()
            if g and g != "NA":
                out[f[ib]] = g
    return out


def load_well_map(path: str, use_rep: bool = False) -> Dict[str, str]:
    """Return {well: genotype} (or genotype_rep). Same contract as 06_41."""
    col = "genotype_rep" if use_rep else "genotype"
    out: Dict[str, str] = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if "well" not in header or col not in header:
            sys.exit(f"[06_44] well-map needs columns 'well' and '{col}'")
        iw, ig = header.index("well"), header.index(col)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(iw, ig):
                continue
            well, geno = f[iw].strip(), f[ig].strip()
            if well and geno and geno != "NA":
                out[well] = geno
    return out


# ------------------------------------------------------------- worker (shared)

_W_BAM: Optional[str] = None
_W_BC_TAG: str = "BC"
_W_MODE: str = "bc_key_map"
_W_BC_MAP: Optional[Dict[str, str]] = None
_W_WELL_MAP: Optional[Dict[str, str]] = None
_W_MIN_MAPQ: int = 10


def _init_worker(bam, bc_tag, mode, bc_map, well_map, min_mapq):
    global _W_BAM, _W_BC_TAG, _W_MODE, _W_BC_MAP, _W_WELL_MAP, _W_MIN_MAPQ
    _W_BAM, _W_BC_TAG, _W_MODE = bam, bc_tag, mode
    _W_BC_MAP, _W_WELL_MAP, _W_MIN_MAPQ = bc_map, well_map, min_mapq


def _group_for(bc: str) -> Optional[str]:
    """Map a BAM BC tag value to its pseudo-bulk group, or None."""
    if _W_MODE == "bc_key_map":
        assert _W_BC_MAP is not None
        return _W_BC_MAP.get(bc.split("-", 1)[0])
    # plate: well is parts[-3] of the underscore-split BC suffix (see 06_41)
    assert _W_WELL_MAP is not None
    try:
        parts = bc.split("-", 1)[1].split("_")
        if len(parts) < 4:
            return None
        return _W_WELL_MAP.get(parts[-3])
    except IndexError:
        return None


def _process_chrom(chrom_snps: Tuple[str, List[Tuple[int, str, str]]]):
    chrom, snp_list = chrom_snps
    assert _W_BAM is not None
    snp_pos = [s[0] for s in snp_list]
    snp_info = {s[0]: (s[1], s[2]) for s in snp_list}

    # counts[(group, pos)] = [n_ref, n_alt, n_other]
    counts: Dict[Tuple[str, int], List[int]] = defaultdict(lambda: [0, 0, 0])
    # bias[(group, class)] = [n_obs, sum_nm, sum_mapq]; class in {"ref","alt"}
    bias: Dict[Tuple[str, str], List[int]] = defaultdict(lambda: [0, 0, 0])

    bam = pysam.AlignmentFile(_W_BAM, "rb")
    try:
        if chrom not in bam.references:
            return chrom, {}, {}, {}
        for read in bam.fetch(chrom):
            if (read.is_unmapped or read.is_secondary or read.is_supplementary
                    or read.is_duplicate):
                continue
            if read.mapping_quality < _W_MIN_MAPQ:
                continue
            rs, re = read.reference_start, read.reference_end
            if re is None:
                continue
            lo = bisect.bisect_left(snp_pos, rs + 1)
            hi = bisect.bisect_right(snp_pos, re)
            if lo == hi:
                continue
            if not read.has_tag(_W_BC_TAG):
                continue
            grp = _group_for(read.get_tag(_W_BC_TAG))
            if grp is None:
                continue

            seq = read.query_sequence
            if seq is None:
                continue
            ref_to_q = {rp: qp for qp, rp in read.get_aligned_pairs(matches_only=True)}
            mapq = read.mapping_quality
            nm = read.get_tag("NM") if read.has_tag("NM") else 0
            for p in snp_pos[lo:hi]:
                q = ref_to_q.get(p - 1)
                if q is None:
                    continue
                b = seq[q].upper()
                ref, alt1 = snp_info[p]
                key = (grp, p)
                if b == ref:
                    counts[key][0] += 1
                    bb = bias[(grp, "ref")]; bb[0] += 1; bb[1] += nm; bb[2] += mapq
                elif b == alt1:
                    counts[key][1] += 1
                    bb = bias[(grp, "alt")]; bb[0] += 1; bb[1] += nm; bb[2] += mapq
                else:
                    counts[key][2] += 1
    finally:
        bam.close()

    return chrom, dict(counts), snp_info, dict(bias)


# ----------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--vcf", required=True, help="Reference VCF (CHROM/POS/REF/ALT only).")
    ap.add_argument("--group-mode", choices=["bc_key_map", "plate"], required=True)
    ap.add_argument("--bc-map", help="bc_key_map mode: TSV with 'bc_key','group' columns.")
    ap.add_argument("--well-map", help="plate mode: TSV with 'well','genotype'[,'genotype_rep'].")
    ap.add_argument("--use-rep", action="store_true", help="plate: use genotype_rep column.")
    ap.add_argument("--bc-tag", default="BC")
    ap.add_argument("--min-mapq", type=int, default=10)
    ap.add_argument("--threads", type=int, default=10)
    ap.add_argument("--chroms", default="", help="Comma-separated chroms to keep (default: all in VCF).")
    ap.add_argument("--out", required=True, help="Output TSV(.gz).")
    args = ap.parse_args()

    print(f"[06_44] loading SNP positions from {args.vcf}", file=sys.stderr)
    snps = load_snp_positions(args.vcf)
    if args.chroms.strip():
        keep = set(args.chroms.split(","))
        snps = {c: v for c, v in snps.items() if c in keep}
    n_snps = sum(len(v) for v in snps.values())
    print(f"[06_44]   {n_snps:,} SNVs across {len(snps)} chromosome(s)", file=sys.stderr)

    if args.group_mode == "bc_key_map":
        if not args.bc_map:
            sys.exit("[06_44] --group-mode bc_key_map needs --bc-map")
        bc_map = load_bc_map(args.bc_map)
        well_map = None
        groups = sorted(set(bc_map.values()))
        print(f"[06_44] bc_key_map: {len(bc_map):,} barcodes -> {len(groups)} groups: "
              f"{', '.join(groups)}", file=sys.stderr)
    else:
        if not args.well_map:
            sys.exit("[06_44] --group-mode plate needs --well-map")
        well_map = load_well_map(args.well_map, use_rep=args.use_rep)
        bc_map = None
        groups = sorted(set(well_map.values()))
        print(f"[06_44] plate: {len(well_map)} wells -> {len(groups)} groups: "
              f"{', '.join(groups)}", file=sys.stderr)

    work = [(c, snps[c]) for c in sorted(snps.keys())]
    print(f"[06_44] processing {len(work)} chromosome(s) x {args.threads} workers", file=sys.stderr)

    # accumulate the small bias table in main; stream the big counts per chrom
    bias_total: Dict[Tuple[str, str], List[int]] = defaultdict(lambda: [0, 0, 0])
    grp_sites: Dict[str, int] = defaultdict(int)
    grp_reads: Dict[str, int] = defaultdict(int)

    op = gzip.open if args.out.endswith(".gz") else open
    n_rows = 0
    with op(args.out, "wt") as out:
        out.write("chrom\tpos\tref\talt\tgroup\tn_ref\tn_alt\tn_other\n")
        with ProcessPoolExecutor(
            max_workers=args.threads, initializer=_init_worker,
            initargs=(args.bam, args.bc_tag, args.group_mode, bc_map, well_map, args.min_mapq),
        ) as ex:
            futs = [ex.submit(_process_chrom, w) for w in work]
            for fut in as_completed(futs):
                chrom, counts, snp_info, bias = fut.result()
                for (grp, pos), (nr, na, no) in counts.items():
                    ref, alt1 = snp_info[pos]
                    out.write(f"{chrom}\t{pos}\t{ref}\t{alt1}\t{grp}\t{nr}\t{na}\t{no}\n")
                    grp_sites[grp] += 1
                    grp_reads[grp] += nr + na + no
                for k, v in bias.items():
                    b = bias_total[k]
                    b[0] += v[0]; b[1] += v[1]; b[2] += v[2]
                n_rows += len(counts)
                print(f"[06_44]   {chrom}: {len(counts):,} group x site rows", file=sys.stderr)

    bias_path = (args.out[:-3] if args.out.endswith(".gz") else args.out) + ".bias.tsv"
    with open(bias_path, "wt") as bf:
        bf.write("group\tobs_class\tn_obs\tsum_nm\tsum_mapq\n")
        for (grp, cls), (n, snm, smq) in sorted(bias_total.items()):
            bf.write(f"{grp}\t{cls}\t{n}\t{snm}\t{smq}\n")

    print(f"[06_44] DONE — {n_rows:,} rows -> {args.out}", file=sys.stderr)
    for g in groups:
        print(f"[06_44]   {g}: {grp_sites[g]:,} covered sites, {grp_reads[g]:,} allele obs",
              file=sys.stderr)
    print(f"[06_44] bias sidecar -> {bias_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
