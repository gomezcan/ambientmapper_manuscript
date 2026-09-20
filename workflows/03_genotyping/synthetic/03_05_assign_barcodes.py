#!/usr/bin/env python3
"""
03_05_assign_barcodes.py — Phase 3: Assign reads to synthetic barcodes
                           (Titration design)

Produces 15 independent datasets from a shared set of template barcodes.
Cell reads are identical across all datasets; only contamination reads change.

Titration levels:
  alpha = 0 (baseline) + [0.02, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50]
  contaminants = Il14H, Ki11

Template barcodes (B73-skewed, total = 416 per dataset):
  - 300 singlets (B73, 6 depth bins × 50 per bin)
  - 96 doublets (B73+Il14H, B73+Ki11 at 4 splits × 3 depths × 4 reps)
  - 20 empty barcodes

Read name format: @{barcode}|{source_genome}|{original_art_name}/1

Input:  synthetic/reads/{B73,Il14H,Ki11}_{1,2}.fq.gz
Output: synthetic/barcoded/templates.tsv
        synthetic/barcoded/alpha_000/  (baseline)
        synthetic/barcoded/alpha_002_Il14H/
        ...
        synthetic/barcoded/alpha_050_Ki11/
        Each with: all_reads_{R1,R2}.fq.gz + truth_table.tsv
"""

import argparse
import gzip
import os
import sys
import numpy as np


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

GENOMES = ["B73", "Il14H", "Ki11"]
DEPTHS = [250, 500, 1000, 2000, 5000, 7500]
ALPHAS_NONZERO = [0.02, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50]
CONTAMINANTS = ["Il14H", "Ki11"]
DOUBLET_PAIRS = [("B73", "Il14H"), ("B73", "Ki11")]
DOUBLET_SPLITS = [0.50, 0.60, 0.70, 0.80]
DOUBLET_DEPTHS = [1000, 3000, 7500]
EMPTY_DEPTH_RANGE = (10, 50)

DEPTH_BIN_LABELS = {250: "250", 500: "500", 1000: "1K",
                    2000: "2K", 5000: "5K", 7500: "7.5K"}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--reads-dir", default="synthetic/reads",
                   help="Directory with ART FASTQs (default: synthetic/reads)")
    p.add_argument("--outdir", default="synthetic/barcoded",
                   help="Output directory (default: synthetic/barcoded)")
    p.add_argument("--seed", type=int, default=42,
                   help="Random seed (default: 42)")
    p.add_argument("--depth-jitter", type=float, default=0.10,
                   help="Relative uniform jitter on target depth (default: 0.10)")
    p.add_argument("--n-singlet-per-bin", type=int, default=50,
                   help="Singlet templates per depth bin (default: 50)")
    p.add_argument("--n-doublet-reps", type=int, default=4,
                   help="Replicates per doublet config (default: 4)")
    p.add_argument("--n-empty", type=int, default=20,
                   help="Number of empty barcodes (default: 20)")
    return p.parse_args()


# ---------------------------------------------------------------------------
# FASTQ I/O
# ---------------------------------------------------------------------------

def load_read_pool(r1_path, r2_path):
    """Load paired FASTQ files into list of (base_name, seq1, qual1, seq2, qual2)."""
    pairs = []
    with gzip.open(r1_path, "rt") as f1, gzip.open(r2_path, "rt") as f2:
        while True:
            name1 = f1.readline().strip()
            if not name1:
                break
            seq1 = f1.readline().strip()
            f1.readline()  # +
            qual1 = f1.readline().strip()

            f2.readline()  # name2 (not needed — same base name)
            seq2 = f2.readline().strip()
            f2.readline()  # +
            qual2 = f2.readline().strip()

            # Strip @ and /1 suffix to get base name
            base = name1[1:]
            if base.endswith("/1"):
                base = base[:-2]
            pairs.append((base, seq1, qual1, seq2, qual2))
    return pairs


# ---------------------------------------------------------------------------
# Template generation
# ---------------------------------------------------------------------------

def generate_templates(rng, args):
    """Generate template barcodes (shared across all datasets)."""
    templates = []
    bc_id = 0

    def add(bc_type, **kwargs):
        nonlocal bc_id
        kwargs["barcode"] = f"SYN{bc_id:06d}"
        kwargs["type"] = bc_type
        templates.append(kwargs)
        bc_id += 1

    # Singlets (B73)
    for depth in DEPTHS:
        for _ in range(args.n_singlet_per_bin):
            add("singlet", true_genome="B73",
                target_depth=depth, depth_bin=DEPTH_BIN_LABELS[depth])

    # Doublets
    for g1, g2 in DOUBLET_PAIRS:
        for rho in DOUBLET_SPLITS:
            for depth in DOUBLET_DEPTHS:
                for _ in range(args.n_doublet_reps):
                    add("doublet", true_genome=f"{g1}+{g2}",
                        genome_1=g1, genome_2=g2, rho=rho,
                        target_depth=depth, depth_bin=f"{depth}")

    # Empty
    for _ in range(args.n_empty):
        depth = int(rng.integers(EMPTY_DEPTH_RANGE[0], EMPTY_DEPTH_RANGE[1]))
        add("empty", true_genome="none",
            target_depth=depth, depth_bin="empty")

    # Apply depth jitter
    for t in templates:
        target = t["target_depth"]
        jitter = rng.uniform(1 - args.depth_jitter, 1 + args.depth_jitter)
        t["actual_depth"] = max(5, int(round(target * jitter)))

    return templates


# ---------------------------------------------------------------------------
# Pool partitioning
# ---------------------------------------------------------------------------

def partition_pools(pools, templates, rng):
    """
    Two-phase pool partitioning:
      Phase A: Allocate cell reads (fixed, shared across all datasets)
      Phase B: Reserve contamination slices (unique per dataset)

    Returns:
      cell_reads: dict barcode -> list of (name, seq1, qual1, seq2, qual2, source_genome)
      contam_slices: dict contaminant -> dict alpha -> list of read tuples
    """
    offsets = {g: 0 for g in GENOMES}

    def consume(genome, n):
        start = offsets[genome]
        end = start + n
        if end > len(pools[genome]):
            print(f"ERROR: exhausted {genome} pool at offset {start}, need {n} more",
                  file=sys.stderr)
            sys.exit(1)
        batch = pools[genome][start:end]
        offsets[genome] = end
        return batch

    # --- Phase A: Cell reads ---
    cell_reads = {}
    for t in templates:
        barcode = t["barcode"]
        d = t["actual_depth"]
        reads = []

        if t["type"] == "singlet":
            for pair in consume("B73", d):
                reads.append((*pair, "B73"))

        elif t["type"] == "doublet":
            n_g1 = int(round(t["rho"] * d))
            n_g2 = d - n_g1
            for pair in consume(t["genome_1"], n_g1):
                reads.append((*pair, t["genome_1"]))
            for pair in consume(t["genome_2"], n_g2):
                reads.append((*pair, t["genome_2"]))

        elif t["type"] == "empty":
            per_g = d // 3
            remainder = d - per_g * 3
            for i, g in enumerate(GENOMES):
                n = per_g + (1 if i < remainder else 0)
                for pair in consume(g, n):
                    reads.append((*pair, g))

        cell_reads[barcode] = reads

    print(f"  Cell reads allocated:")
    for g in GENOMES:
        print(f"    {g}: {offsets[g]:,} read pairs consumed")

    # --- Phase B: Contamination read slices ---
    # For each (alpha, contaminant), compute total contamination reads needed
    # and slice from the remaining pool
    contam_slices = {c: {} for c in CONTAMINANTS}

    for contaminant in CONTAMINANTS:
        for alpha in ALPHAS_NONZERO:
            total_contam = 0
            for t in templates:
                if t["type"] == "empty":
                    continue
                n_cell = t["actual_depth"]
                n_contam = max(0, int(round(n_cell * alpha / (1 - alpha))))
                total_contam += n_contam

            batch = consume(contaminant, total_contam)
            contam_slices[contaminant][alpha] = batch

        print(f"  {contaminant} contamination reserved: "
              f"{offsets[contaminant]:,} / {len(pools[contaminant]):,} "
              f"({offsets[contaminant]/len(pools[contaminant])*100:.1f}%)")

    return cell_reads, contam_slices


# ---------------------------------------------------------------------------
# Dataset writing
# ---------------------------------------------------------------------------

def write_dataset(templates, cell_reads, contam_reads, alpha, contaminant,
                  rng, outdir):
    """
    Write one dataset's FASTQs and truth table.

    contam_reads: list of read tuples for this (alpha, contaminant), or None for baseline.
    """
    label = "alpha_000" if alpha == 0.0 else f"alpha_{int(alpha*100):03d}_{contaminant}"
    ds_dir = os.path.join(outdir, label)
    os.makedirs(ds_dir, exist_ok=True)

    r1_path = os.path.join(ds_dir, "all_reads_R1.fq.gz")
    r2_path = os.path.join(ds_dir, "all_reads_R2.fq.gz")
    truth_path = os.path.join(ds_dir, "truth_table.tsv")

    contam_offset = 0
    total_reads = 0
    bc_meta = []  # per-barcode metadata for truth table

    with gzip.open(r1_path, "wt") as r1_out, gzip.open(r2_path, "wt") as r2_out:
        for t in templates:
            barcode = t["barcode"]

            # Start with cell reads (shared, read-only)
            reads = list(cell_reads[barcode])  # shallow copy
            n_cell = len(reads)
            n_contam = 0

            # Add contamination reads (singlets and doublets, not empties)
            if alpha > 0 and t["type"] != "empty" and contam_reads is not None:
                n_contam = max(0, int(round(t["actual_depth"] * alpha / (1 - alpha))))
                for pair in contam_reads[contam_offset:contam_offset + n_contam]:
                    reads.append((*pair, contaminant))
                contam_offset += n_contam

            # Shuffle within barcode
            idx = rng.permutation(len(reads))

            for j in idx:
                name, seq1, qual1, seq2, qual2, source = reads[j]
                tagged = f"{barcode}|{source}|{name}"
                r1_out.write(f"@{tagged}/1\n{seq1}\n+\n{qual1}\n")
                r2_out.write(f"@{tagged}/2\n{seq2}\n+\n{qual2}\n")

            total_reads += len(reads)
            bc_meta.append({
                "n_cell_reads": n_cell,
                "n_contam_reads": n_contam,
                "n_total_reads": len(reads),
            })

    # Write truth table
    write_truth_table(templates, alpha, contaminant, bc_meta, truth_path)

    return label, total_reads


def write_truth_table(templates, alpha, contaminant, bc_meta, path):
    """Write per-barcode ground truth for one dataset."""
    header = ["barcode", "type", "true_genome", "contaminant", "alpha",
              "n_cell_reads", "n_contam_reads", "n_total_reads",
              "target_depth", "actual_depth",
              "genome_1", "genome_2", "rho"]
    with open(path, "w") as f:
        f.write("\t".join(header) + "\n")
        for t, meta in zip(templates, bc_meta):
            contam_label = contaminant if alpha > 0 else "none"
            row = [
                t["barcode"],
                t["type"],
                t.get("true_genome", "none"),
                contam_label,
                f"{alpha:.2f}",
                str(meta["n_cell_reads"]),
                str(meta["n_contam_reads"]),
                str(meta["n_total_reads"]),
                str(t["target_depth"]),
                str(t["actual_depth"]),
                t.get("genome_1", "NA"),
                t.get("genome_2", "NA"),
                f"{t['rho']:.2f}" if t["type"] == "doublet" else "NA",
            ]
            f.write("\t".join(row) + "\n")


def write_templates_file(templates, path):
    """Write shared template definitions."""
    header = ["barcode", "type", "true_genome", "target_depth", "actual_depth",
              "depth_bin", "genome_1", "genome_2", "rho"]
    with open(path, "w") as f:
        f.write("\t".join(header) + "\n")
        for t in templates:
            row = [
                t["barcode"],
                t["type"],
                t.get("true_genome", "none"),
                str(t["target_depth"]),
                str(t["actual_depth"]),
                t.get("depth_bin", "NA"),
                t.get("genome_1", "NA"),
                t.get("genome_2", "NA"),
                f"{t['rho']:.2f}" if t["type"] == "doublet" else "NA",
            ]
            f.write("\t".join(row) + "\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    rng = np.random.default_rng(args.seed)
    os.makedirs(args.outdir, exist_ok=True)

    # --- Generate templates ---
    templates = generate_templates(rng, args)

    types = {}
    for t in templates:
        types[t["type"]] = types.get(t["type"], 0) + 1
    print(f"Template barcodes: {len(templates)} total")
    for tp, n in sorted(types.items()):
        print(f"  {tp}: {n}")

    total_cell = sum(t["actual_depth"] for t in templates)
    print(f"  Total cell reads: {total_cell:,}")

    # --- Write templates file ---
    templates_path = os.path.join(args.outdir, "templates.tsv")
    write_templates_file(templates, templates_path)
    print(f"  Wrote {templates_path}")

    # --- Load and shuffle read pools ---
    pools = {}
    for g in GENOMES:
        r1 = os.path.join(args.reads_dir, f"{g}_1.fq.gz")
        r2 = os.path.join(args.reads_dir, f"{g}_2.fq.gz")
        print(f"\nLoading {g} reads from {r1}...")
        pairs = load_read_pool(r1, r2)
        print(f"  Loaded {len(pairs):,} read pairs, shuffling...")
        rng.shuffle(pairs)
        pools[g] = pairs

    # --- Partition pools ---
    print(f"\nPartitioning read pools...")
    cell_reads, contam_slices = partition_pools(pools, templates, rng)

    # Free pool memory (cell_reads and contam_slices hold references to data)
    del pools

    # --- Generate datasets ---
    n_datasets = 1 + len(ALPHAS_NONZERO) * len(CONTAMINANTS)
    print(f"\nWriting {n_datasets} datasets...")

    # Dataset 1: baseline (alpha=0)
    label, n_reads = write_dataset(
        templates, cell_reads, None, 0.0, "none", rng, args.outdir)
    print(f"  {label}: {n_reads:,} read pairs")

    # Datasets 2-15: contamination titration
    for contaminant in CONTAMINANTS:
        for alpha in ALPHAS_NONZERO:
            label, n_reads = write_dataset(
                templates, cell_reads, contam_slices[contaminant][alpha],
                alpha, contaminant, rng, args.outdir)
            print(f"  {label}: {n_reads:,} read pairs")

    # --- Summary ---
    print(f"\n=== Phase 3 summary ===")
    print(f"  Templates:  {len(templates)} barcodes")
    print(f"  Datasets:   {n_datasets}")
    print(f"  Output dir: {args.outdir}")
    print(f"\n  Template types:")
    for tp in sorted(types):
        subset = [t for t in templates if t["type"] == tp]
        depths = [t["actual_depth"] for t in subset]
        print(f"    {tp}: {len(subset)} (depth {min(depths)}-{max(depths)})")
    print(f"\n  Datasets:")
    print(f"    alpha_000 (baseline)")
    for alpha in ALPHAS_NONZERO:
        pct = int(alpha * 100)
        print(f"    alpha_{pct:03d}_Il14H, alpha_{pct:03d}_Ki11")
    print(f"\nDone.")


if __name__ == "__main__":
    main()
