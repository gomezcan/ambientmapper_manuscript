#!/usr/bin/env python3
"""
03_07_extract_metrics.py — Phase 5: Extract QCMapping from synthetic BAMs

Parses BWA-mapped BAMs with pysam and writes per-genome QCMapping files in the
format ambientmapper expects (Read, BC, MAPQ, AS, NM, XAcount, frag_loc).

Bypasses `ambientmapper extract` (which requires CB tags in BAMs). Instead,
barcodes are parsed from read names: @{barcode}|{source_genome}|{art_id}

frag_loc format: chr:R1_pos:mate_pos (matching ambientmapper extract output,
needed by _friend_rescue in the assign step).

Only outputs R1 (first-in-pair) records to match one-record-per-pair convention
used by the real pipeline.

Usage:
  python 03_07_extract_metrics.py \
    --bam-dir synthetic/mapping/alpha_000 \
    --outdir synthetic/alpha_000/qc \
    --genome-map B73v5:B73 Il14H:Il14H Ki11:Ki11
"""

import argparse
import os
import sys

try:
    import pysam
except ImportError:
    print("Error: pysam required. Install with: pip install pysam", file=sys.stderr)
    sys.exit(1)


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--bam-dir", required=True,
                   help="Directory with sorted BAMs (syn_{bwa_genome}.sorted.bam)")
    p.add_argument("--outdir", required=True,
                   help="Output directory for QCMapping files")
    p.add_argument("--genome-map", nargs="+", required=True,
                   help="BWA genome name to config name mapping (e.g., B73v5:B73 Il14H:Il14H)")
    p.add_argument("--sample-suffix", default="",
                   help="Suffix to append to barcode (e.g., '-SynBench')")
    return p.parse_args()


def parse_barcode(qname):
    """Extract barcode from read name: {barcode}|{source_genome}|{art_id}"""
    parts = qname.split("|")
    if len(parts) >= 1:
        return parts[0]
    return qname


def extract_qcmapping(bam_path, output_path, sample_suffix=""):
    """
    Extract QCMapping from a BAM file.
    Outputs one row per read pair (R1 / first-in-pair only).
    """
    bam = pysam.AlignmentFile(bam_path, "rb")
    n_total = 0
    n_written = 0
    n_unmapped = 0
    n_secondary = 0
    n_r2_skip = 0

    with open(output_path, "w") as out:
        out.write("Read\tBC\tMAPQ\tAS\tNM\tXAcount\tfrag_loc\n")

        for aln in bam.fetch(until_eof=True):
            n_total += 1

            # Skip unmapped
            if aln.is_unmapped:
                n_unmapped += 1
                continue

            # Skip secondary and supplementary alignments
            if aln.is_secondary or aln.is_supplementary:
                n_secondary += 1
                continue

            # Only output R1 (first-in-pair) to match one-per-pair convention
            if aln.is_paired and aln.is_read2:
                n_r2_skip += 1
                continue

            qname = aln.query_name
            barcode = parse_barcode(qname) + sample_suffix
            mapq = aln.mapping_quality
            AS = aln.get_tag("AS") if aln.has_tag("AS") else 0
            NM = aln.get_tag("NM") if aln.has_tag("NM") else 0

            # XA tag: alternative hits; count semicolons
            if aln.has_tag("XA"):
                xa_str = aln.get_tag("XA")
                xa_count = xa_str.count(";")
            else:
                xa_count = 0

            # Fragment location: chr:R1_pos:mate_pos
            chrom = aln.reference_name or "*"
            pos = aln.reference_start
            mate_unmapped = aln.mate_is_unmapped if aln.is_paired else True
            mate_pos = aln.next_reference_start if (aln.is_paired and not mate_unmapped) else -1
            frag_loc = f"{chrom}:{pos}:{mate_pos}"

            out.write(f"{qname}\t{barcode}\t{mapq}\t{AS}\t{NM}\t{xa_count}\t{frag_loc}\n")
            n_written += 1

    bam.close()
    return n_total, n_written, n_unmapped, n_secondary, n_r2_skip


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # Parse genome mapping (BWA name -> config name)
    genome_map = {}
    for mapping in args.genome_map:
        parts = mapping.split(":")
        if len(parts) != 2:
            print(f"Error: invalid genome mapping '{mapping}', expected BWA:CONFIG",
                  file=sys.stderr)
            sys.exit(1)
        genome_map[parts[0]] = parts[1]

    print(f"Genome mapping: {genome_map}")
    print(f"BAM dir: {args.bam_dir}")
    print(f"Output:  {args.outdir}")
    print()

    for bwa_name, config_name in genome_map.items():
        bam_path = os.path.join(args.bam_dir, f"syn_{bwa_name}.sorted.bam")
        if not os.path.exists(bam_path):
            print(f"Warning: BAM not found: {bam_path}", file=sys.stderr)
            continue

        output_path = os.path.join(args.outdir, f"{config_name}_QCMapping.txt")
        print(f"[{config_name}] Extracting from {bam_path}...")

        n_total, n_written, n_unmapped, n_secondary, n_r2_skip = extract_qcmapping(
            bam_path, output_path, args.sample_suffix
        )

        print(f"  Total alignments: {n_total:,}")
        print(f"  Written (R1):     {n_written:,}")
        print(f"  Unmapped:         {n_unmapped:,}")
        print(f"  Secondary/suppl:  {n_secondary:,}")
        print(f"  R2 skipped:       {n_r2_skip:,}")
        print(f"  -> {output_path}")
        print()

    print("Done.")


if __name__ == "__main__":
    main()
