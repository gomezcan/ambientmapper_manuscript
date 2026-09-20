#!/usr/bin/env python3
"""
03_03_filter_disc_peaks.py — Filter ortholog peaks by discriminative power

Selects peaks where BOTH Il14H and Ki11 have ≥N expected SNPs per read
relative to B73. These peaks produce reads that can distinguish genomes,
enabling meaningful genotyping evaluation.

Input:  synthetic/orthologs/ortholog_map.tsv + {B73,Il14H,Ki11}_peak_seqs.fa
Output: synthetic_disc/orthologs/ortholog_map.tsv + filtered FASTAs

Usage:
  python 03_03_filter_disc_peaks.py [--min-snps-per-read 1.0] [--read-length 75]
"""

import argparse
import os
import sys


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ortho-map",
                   default="synthetic/orthologs/ortholog_map.tsv",
                   help="Input ortholog map TSV (default: synthetic/orthologs/ortholog_map.tsv)")
    p.add_argument("--fasta-dir",
                   default="synthetic/orthologs",
                   help="Directory containing {B73,Il14H,Ki11}_peak_seqs.fa")
    p.add_argument("--outdir",
                   default="synthetic_disc/orthologs",
                   help="Output directory (default: synthetic_disc/orthologs)")
    p.add_argument("--min-snps-per-read", type=float, default=1.0,
                   help="Min expected SNPs per read for BOTH genomes (default: 1.0)")
    p.add_argument("--read-length", type=int, default=75,
                   help="Read length for SNP calculation (default: 75)")
    return p.parse_args()


def load_ortholog_map(path):
    """Load ortholog_map.tsv, return header + list of row dicts."""
    rows = []
    with open(path) as f:
        header = f.readline().strip()
        col_names = header.split("\t")
        for line in f:
            fields = line.strip().split("\t")
            row = dict(zip(col_names, fields))
            rows.append(row)
    return header, col_names, rows


def filter_peaks(rows, read_length, min_snps):
    """Filter to peaks with ≥min_snps expected SNPs per read in BOTH genomes."""
    kept = []
    for row in rows:
        il14h_id = float(row["Il14H_identity"])
        ki11_id = float(row["Ki11_identity"])
        il14h_snps = read_length * (1.0 - il14h_id)
        ki11_snps = read_length * (1.0 - ki11_id)
        min_val = min(il14h_snps, ki11_snps)
        if min_val >= min_snps:
            row["_il14h_snps_per_read"] = il14h_snps
            row["_ki11_snps_per_read"] = ki11_snps
            row["_min_snps_per_read"] = min_val
            kept.append(row)
    return kept


def subset_fasta(input_path, output_path, peak_ids):
    """Subset a FASTA file to only include sequences whose header matches peak_ids."""
    peak_set = set(peak_ids)
    writing = False
    n_written = 0
    with open(input_path) as fin, open(output_path, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                # Header format: >peak_id::chrom:start-end(strand)
                seq_id = line[1:].split("::")[0].strip()
                writing = seq_id in peak_set
                if writing:
                    fout.write(line)
                    n_written += 1
            elif writing:
                fout.write(line)
    return n_written


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # Load
    print(f"Loading {args.ortho_map}...")
    header, col_names, rows = load_ortholog_map(args.ortho_map)
    print(f"  {len(rows)} peaks loaded")

    # Filter
    kept = filter_peaks(rows, args.read_length, args.min_snps_per_read)
    print(f"\nFilter: ≥{args.min_snps_per_read} SNPs per {args.read_length}bp read (both genomes)")
    print(f"  Kept: {len(kept)} / {len(rows)} peaks ({len(kept)/len(rows)*100:.1f}%)")

    if len(kept) == 0:
        print("ERROR: No peaks passed filter. Lower --min-snps-per-read?")
        sys.exit(1)

    # Statistics
    il_snps = [r["_il14h_snps_per_read"] for r in kept]
    ki_snps = [r["_ki11_snps_per_read"] for r in kept]
    print(f"\n  Il14H SNPs/read: min={min(il_snps):.1f}, median={sorted(il_snps)[len(il_snps)//2]:.1f}, max={max(il_snps):.1f}")
    print(f"  Ki11  SNPs/read: min={min(ki_snps):.1f}, median={sorted(ki_snps)[len(ki_snps)//2]:.1f}, max={max(ki_snps):.1f}")

    # Estimate read pool at different coverage levels
    total_len = sum(int(r["B73_end"]) - int(r["B73_start"]) for r in kept)
    print(f"\n  Total sequence length: {total_len/1000:.1f} Kbp (avg {total_len/len(kept):.0f} bp/peak)")
    for cov in [100, 250, 500]:
        reads = total_len * cov // args.read_length
        print(f"  At {cov}x coverage: ~{reads/1000:.0f}K reads/genome")

    # Write filtered ortholog_map.tsv
    out_map = os.path.join(args.outdir, "ortholog_map.tsv")
    with open(out_map, "w") as f:
        f.write(header + "\n")
        for row in kept:
            f.write("\t".join(row[c] for c in col_names) + "\n")
    print(f"\nWrote {out_map} ({len(kept)} peaks)")

    # Subset FASTAs
    peak_ids = [r["peak_id"] for r in kept]
    genomes = ["B73", "Il14H", "Ki11"]
    for g in genomes:
        in_fa = os.path.join(args.fasta_dir, f"{g}_peak_seqs.fa")
        out_fa = os.path.join(args.outdir, f"{g}_peak_seqs.fa")
        n = subset_fasta(in_fa, out_fa, peak_ids)
        print(f"  {g}: {n} sequences written to {out_fa}")

    print(f"\nDone. Output directory: {args.outdir}")


if __name__ == "__main__":
    main()
