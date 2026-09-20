#!/usr/bin/env python3
"""
03_02_find_orthologs.py — Phase 1: Find orthologous peak regions across genomes

Given a BED of B73 peaks and minimap2 alignments of those peaks to other genomes,
parse the alignments, filter for high-quality single-hit orthologs, and extract
the orthologous sequences.

Usage:
  python 03_02_find_orthologs.py \
    --peaks synthetic/peaks/B73_peaks_filtered.bed \
    --b73-fasta /path/to/Zm-B73-REFERENCE-NAM-5.0.chrs.fa \
    --target-genomes Il14H Ki11 \
    --target-fastas /path/to/Il14H.fa /path/to/Ki11.fa \
    --sam-dir synthetic/orthologs/alignments \
    --outdir synthetic/orthologs \
    --min-coverage 0.90 \
    --min-identity 0.90
"""

import argparse
import os
import re
import subprocess
import sys
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--peaks", required=True,
                   help="Filtered B73 peaks BED (from 03_01_call_peaks.sh)")
    p.add_argument("--b73-fasta", required=True,
                   help="B73 reference FASTA (for bedtools getfasta)")
    p.add_argument("--target-genomes", nargs="+", required=True,
                   help="Target genome names (e.g., Il14H Ki11)")
    p.add_argument("--target-fastas", nargs="+", required=True,
                   help="Target genome FASTA paths (same order as --target-genomes)")
    p.add_argument("--sam-dir", required=True,
                   help="Directory containing minimap2 SAM files (one per target)")
    p.add_argument("--outdir", required=True,
                   help="Output directory for ortholog FASTAs and map")
    p.add_argument("--min-coverage", type=float, default=0.90,
                   help="Minimum alignment coverage of B73 peak (0-1, default: 0.90)")
    p.add_argument("--min-identity", type=float, default=0.90,
                   help="Minimum alignment identity (0-1, default: 0.90)")
    p.add_argument("--max-hits", type=int, default=1,
                   help="Maximum allowed hits per peak in target genome (default: 1)")
    return p.parse_args()


def load_peaks(bed_path):
    """Load BED file, return list of (chrom, start, end, peak_id)."""
    peaks = []
    with open(bed_path) as f:
        for i, line in enumerate(f):
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            chrom, start, end = cols[0], int(cols[1]), int(cols[2])
            # Use name column if present (col 4), else generate
            peak_id = cols[3] if len(cols) > 3 else f"peak_{i:05d}"
            peaks.append((chrom, start, end, peak_id))
    return peaks


def parse_sam_alignments(sam_path, peaks_by_name, min_coverage, min_identity, max_hits):
    """
    Parse minimap2 SAM output. For each query (B73 peak), find the best
    alignment in the target genome.

    Returns: dict peak_id -> {target_chrom, target_start, target_end, strand,
                               coverage, identity, mapq}
    """
    # Count hits per query to filter multi-mappers
    hit_counts = defaultdict(int)
    hits = defaultdict(list)

    with open(sam_path) as f:
        for line in f:
            if line.startswith("@"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 11:
                continue

            qname = cols[0]
            flag = int(cols[1])
            rname = cols[2]
            pos = int(cols[3])  # 1-based
            mapq = int(cols[4])
            cigar = cols[5]
            seq = cols[9]

            # Skip unmapped
            if flag & 4 or rname == "*":
                continue

            # Parse CIGAR to get alignment length on reference and query
            ref_consumed = 0
            query_consumed = 0
            query_aligned = 0
            n_matches = 0
            n_mismatches = 0

            # Parse optional tags for NM (edit distance)
            nm = 0
            for tag_col in cols[11:]:
                if tag_col.startswith("NM:i:"):
                    nm = int(tag_col.split(":")[2])
                    break

            # Parse CIGAR
            for length, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar):
                length = int(length)
                if op in ('M', '=', 'X'):
                    ref_consumed += length
                    query_consumed += length
                    query_aligned += length
                elif op == 'I':
                    query_consumed += length
                elif op == 'D':
                    ref_consumed += length
                elif op in ('S', 'H'):
                    if op == 'S':
                        query_consumed += length
                elif op == 'N':
                    ref_consumed += length

            # Get query length from peaks_by_name
            if qname not in peaks_by_name:
                continue
            query_length = peaks_by_name[qname][2] - peaks_by_name[qname][1]

            # Coverage = fraction of query aligned
            coverage = query_aligned / query_length if query_length > 0 else 0

            # Identity = (aligned bases - edit distance) / aligned bases
            identity = (query_aligned - nm) / query_aligned if query_aligned > 0 else 0

            # Strand
            strand = "-" if (flag & 16) else "+"

            # Target coordinates (0-based half-open)
            target_start = pos - 1
            target_end = target_start + ref_consumed

            hit_counts[qname] += 1
            hits[qname].append({
                "target_chrom": rname,
                "target_start": target_start,
                "target_end": target_end,
                "strand": strand,
                "coverage": coverage,
                "identity": identity,
                "mapq": mapq,
                "nm": nm,
                "ref_consumed": ref_consumed,
                "query_aligned": query_aligned,
            })

    # Filter: single best hit per peak
    results = {}
    n_multi = 0
    n_low_cov = 0
    n_low_id = 0
    n_pass = 0

    for peak_id, peak_hits in hits.items():
        if hit_counts[peak_id] > max_hits:
            n_multi += 1
            continue

        # Take best hit by identity, then coverage
        best = max(peak_hits, key=lambda h: (h["identity"], h["coverage"]))

        if best["coverage"] < min_coverage:
            n_low_cov += 1
            continue
        if best["identity"] < min_identity:
            n_low_id += 1
            continue

        results[peak_id] = best
        n_pass += 1

    return results, n_pass, n_multi, n_low_cov, n_low_id


def extract_sequences(fasta_path, regions, output_fasta):
    """
    Extract sequences from a FASTA using bedtools getfasta.
    regions: list of (chrom, start, end, name, strand)
    """
    # Write temp BED
    tmp_bed = output_fasta + ".tmp.bed"
    with open(tmp_bed, "w") as f:
        for chrom, start, end, name, strand in regions:
            f.write(f"{chrom}\t{start}\t{end}\t{name}\t0\t{strand}\n")

    # Run bedtools getfasta
    cmd = [
        "bedtools", "getfasta",
        "-fi", fasta_path,
        "-bed", tmp_bed,
        "-fo", output_fasta,
        "-name",
        "-s",  # strand-aware
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"bedtools getfasta error: {result.stderr}", file=sys.stderr)
        sys.exit(1)

    os.remove(tmp_bed)

    # Count sequences
    n_seqs = sum(1 for line in open(output_fasta) if line.startswith(">"))
    return n_seqs


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    if len(args.target_genomes) != len(args.target_fastas):
        print("Error: --target-genomes and --target-fastas must have same length",
              file=sys.stderr)
        sys.exit(1)

    # Load B73 peaks
    peaks = load_peaks(args.peaks)
    print(f"Loaded {len(peaks)} B73 peaks from {args.peaks}")

    # Build lookup by peak name (for SAM parsing)
    peaks_by_name = {}
    for chrom, start, end, peak_id in peaks:
        # bedtools getfasta names queries as chrom:start-end or the name field
        # minimap2 uses the FASTA header as query name
        # We'll use the format that bedtools getfasta -name produces
        key = f"{peak_id}::{chrom}:{start}-{end}"
        peaks_by_name[key] = (chrom, start, end)
        # Also index by just peak_id in case minimap2 truncates
        peaks_by_name[peak_id] = (chrom, start, end)
        # And by coordinate format (bedtools default without -name)
        peaks_by_name[f"{chrom}:{start}-{end}"] = (chrom, start, end)

    # Parse alignments for each target genome
    all_orthologs = {}  # genome -> {peak_id -> hit_info}
    for genome, fasta in zip(args.target_genomes, args.target_fastas):
        sam_path = os.path.join(args.sam_dir, f"B73_to_{genome}.sam")
        if not os.path.exists(sam_path):
            print(f"Error: SAM file not found: {sam_path}", file=sys.stderr)
            sys.exit(1)

        print(f"\nParsing {genome} alignments from {sam_path}...")
        orthologs, n_pass, n_multi, n_low_cov, n_low_id = parse_sam_alignments(
            sam_path, peaks_by_name, args.min_coverage, args.min_identity, args.max_hits
        )
        all_orthologs[genome] = orthologs
        print(f"  {genome}: {n_pass} pass, {n_multi} multi-hit, "
              f"{n_low_cov} low-coverage, {n_low_id} low-identity")

    # Find peaks with valid orthologs in ALL target genomes
    valid_peak_ids = set()
    for chrom, start, end, peak_id in peaks:
        # Check all possible name formats
        found_in_all = True
        for genome in args.target_genomes:
            ortho = all_orthologs[genome]
            # Try different name formats
            found = False
            for key_fmt in [peak_id,
                            f"{peak_id}::{chrom}:{start}-{end}",
                            f"{chrom}:{start}-{end}"]:
                if key_fmt in ortho:
                    found = True
                    break
            if not found:
                found_in_all = False
                break
        if found_in_all:
            valid_peak_ids.add(peak_id)

    print(f"\nPeaks with orthologs in ALL genomes: {len(valid_peak_ids)} / {len(peaks)}")

    if len(valid_peak_ids) == 0:
        print("Error: No valid ortholog triplets found. Check SAM files and name formats.",
              file=sys.stderr)
        sys.exit(1)

    # Write ortholog map
    map_path = os.path.join(args.outdir, "ortholog_map.tsv")
    with open(map_path, "w") as f:
        header = ["peak_id", "B73_chrom", "B73_start", "B73_end"]
        for genome in args.target_genomes:
            header.extend([f"{genome}_chrom", f"{genome}_start", f"{genome}_end",
                           f"{genome}_strand", f"{genome}_coverage", f"{genome}_identity"])
        f.write("\t".join(header) + "\n")

        for chrom, start, end, peak_id in peaks:
            if peak_id not in valid_peak_ids:
                continue
            row = [peak_id, chrom, str(start), str(end)]
            for genome in args.target_genomes:
                ortho = all_orthologs[genome]
                # Find the hit with matching key
                hit = None
                for key_fmt in [peak_id,
                                f"{peak_id}::{chrom}:{start}-{end}",
                                f"{chrom}:{start}-{end}"]:
                    if key_fmt in ortho:
                        hit = ortho[key_fmt]
                        break
                if hit:
                    row.extend([hit["target_chrom"], str(hit["target_start"]),
                                str(hit["target_end"]), hit["strand"],
                                f"{hit['coverage']:.4f}", f"{hit['identity']:.4f}"])
                else:
                    row.extend(["NA"] * 6)
            f.write("\t".join(row) + "\n")

    print(f"Wrote ortholog map: {map_path}")

    # Extract B73 peak sequences (only valid peaks)
    b73_regions = []
    for chrom, start, end, peak_id in peaks:
        if peak_id in valid_peak_ids:
            b73_regions.append((chrom, start, end, peak_id, "+"))

    b73_fasta = os.path.join(args.outdir, "B73_peak_seqs.fa")
    n_b73 = extract_sequences(args.b73_fasta, b73_regions, b73_fasta)
    print(f"Extracted {n_b73} B73 peak sequences -> {b73_fasta}")

    # Extract target genome sequences at orthologous positions
    for genome, fasta in zip(args.target_genomes, args.target_fastas):
        ortho = all_orthologs[genome]
        regions = []
        for chrom, start, end, peak_id in peaks:
            if peak_id not in valid_peak_ids:
                continue
            hit = None
            for key_fmt in [peak_id,
                            f"{peak_id}::{chrom}:{start}-{end}",
                            f"{chrom}:{start}-{end}"]:
                if key_fmt in ortho:
                    hit = ortho[key_fmt]
                    break
            if hit:
                regions.append((hit["target_chrom"], hit["target_start"],
                                hit["target_end"], peak_id, hit["strand"]))

        out_fasta = os.path.join(args.outdir, f"{genome}_peak_seqs.fa")
        n_seqs = extract_sequences(fasta, regions, out_fasta)
        print(f"Extracted {n_seqs} {genome} ortholog sequences -> {out_fasta}")

    # Summary stats
    print(f"\n=== Ortholog finding summary ===")
    print(f"  Input peaks:     {len(peaks)}")
    print(f"  Valid triplets:  {len(valid_peak_ids)}")
    print(f"  Yield:           {100*len(valid_peak_ids)/len(peaks):.1f}%")
    print(f"\nOutput files:")
    print(f"  {map_path}")
    print(f"  {b73_fasta}")
    for genome in args.target_genomes:
        print(f"  {os.path.join(args.outdir, f'{genome}_peak_seqs.fa')}")
    print("\nDone.")


if __name__ == "__main__":
    main()
