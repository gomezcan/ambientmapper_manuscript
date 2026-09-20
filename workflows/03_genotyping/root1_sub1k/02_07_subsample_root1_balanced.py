#!/usr/bin/env python3
"""
02_07_subsample_root1_balanced.py — Draw a depth-balanced barcode panel from
                                    Root1_rep1 for fast config validation.

Reads the per-genome `*_bc_counts.txt` files in 3_Mapping/Root1_rep1/, sums
each barcode's `nuclear` reads across all 26 NAM genomes, bins barcodes by
total depth, and draws a uniform sample within each bin.

Two panels with different seeds are produced in one run (default: A=1, B=2).

Output (per panel, under <outdir>/sub1k_<panel>/):
  barcodes.tsv             — bare 16-mer barcodes, one per line
  barcodes_by_genome.tsv   — per-genome `BC:Z:<bc>-<library>` strings
                             (one row per barcode × genome) for BAM filtering
  barcodes_with_bin.tsv    — barcode, total_nuclear_reads, depth_bin
  manifest.json            — seed, bin counts, source files

Default depth bins (across-genome nuclear read totals):
  b1: [250, 500),  b2: [500, 1K), b3: [1K, 2K),
  b4: [2K, 5K),    b5: [5K, 10K), b6: [10K, inf)

Within each bin, sample uniformly. If a bin has fewer candidates than the
target N, take all of them and log the shortfall (no over-borrowing).
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

# Root of the data tree (3_Mapping/, 5_AmbientDetection/), taken from the
# PROJECT_ROOT environment variable, falling back to the current directory.
_PROJECT_ROOT = os.environ.get("PROJECT_ROOT", os.getcwd())

BC_COUNTS_DIR_DEFAULT = f"{_PROJECT_ROOT}/3_Mapping/Root1_rep1"
OUTDIR_DEFAULT = f"{_PROJECT_ROOT}/5_AmbientDetection/Root1_rep1"
SAMPLE = "Root1_rep1"

# 26 NAM genomes — must match keys in configs/Root1_rep1.ambientmapper.json.
# The library suffix in `cellID` matches the BAM filename stem
# (Root1_rep1_<library>), which is the same string except B73 -> B73v5.
GENOMES = [
    ("B73", "B73v5"),
    ("B97", "B97"),
    ("CML103", "CML103"),
    ("CML228", "CML228"),
    ("CML247", "CML247"),
    ("CML277", "CML277"),
    ("CML322", "CML322"),
    ("CML333", "CML333"),
    ("CML52", "CML52"),
    ("CML69", "CML69"),
    ("HP301", "HP301"),
    ("Il14H", "Il14H"),
    ("Ki11", "Ki11"),
    ("Ki3", "Ki3"),
    ("Ky21", "Ky21"),
    ("M162W", "M162W"),
    ("M37W", "M37W"),
    ("Mo18W", "Mo18W"),
    ("Ms71", "Ms71"),
    ("NC350", "NC350"),
    ("NC358", "NC358"),
    ("Oh43", "Oh43"),
    ("Oh7B", "Oh7B"),
    ("P39", "P39"),
    ("Tx303", "Tx303"),
    ("Tzi8", "Tzi8"),
]

DEFAULT_BINS = [
    ("b1", 250, 500),
    ("b2", 500, 1_000),
    ("b3", 1_000, 2_000),
    ("b4", 2_000, 5_000),
    ("b5", 5_000, 10_000),
    ("b6", 10_000, None),  # None = +inf
]


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--bc-counts-dir", default=BC_COUNTS_DIR_DEFAULT,
                   help="Directory containing Root1_rep1_<library>_bc_counts.txt")
    p.add_argument("--outdir", default=OUTDIR_DEFAULT,
                   help="Parent directory; panels go in <outdir>/sub1k_<panel>/")
    p.add_argument("--n-per-bin", type=int, default=166,
                   help="Target barcodes per depth bin (default: 166 → ~1000 total)")
    p.add_argument("--panels", nargs="+", default=["A:1", "B:2"],
                   help="Panels to draw, format <name>:<seed> (default: A:1 B:2)")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def parse_cellid(cellid: str, library: str) -> str | None:
    """
    cellID format: BC:Z:<16mer>-<SAMPLE>_<library>
    e.g. BC:Z:TGTTGACGATGGGCCT-Root1_rep1_B73v5

    Strip the BC:Z: prefix and the full -<SAMPLE>_<library> suffix,
    return the bare 16-mer. Returns None on format mismatch.
    """
    if not cellid.startswith("BC:Z:"):
        return None
    body = cellid[len("BC:Z:"):]
    suffix = f"-{SAMPLE}_{library}"  # e.g. "-Root1_rep1_B73v5"
    if not body.endswith(suffix):
        return None
    return body[: -len(suffix)]


def load_per_genome_counts(bc_counts_dir: Path, library: str):
    """
    Load one *_bc_counts.txt and return dict bare_bc -> nuclear_reads.

    The files are R-exported TSVs with an unnamed row-name column at the
    front, so the header has 7 columns but each data row has 8 columns.
    pandas.read_csv handles this automatically: the first unnamed column
    becomes the DataFrame index and `cellID` / `nuclear` are addressed by
    name rather than by position.
    """
    path = bc_counts_dir / f"{SAMPLE}_{library}_bc_counts.txt"
    if not path.exists():
        sys.exit(f"ERROR: missing {path}")

    df = pd.read_csv(path, sep="\t")

    for col in ("cellID", "nuclear"):
        if col not in df.columns:
            sys.exit(f"ERROR: {path} missing '{col}' column. "
                     f"Columns: {list(df.columns)}")

    out: dict[str, int] = {}
    for cellid, nuclear in zip(df["cellID"].astype(str),
                               df["nuclear"].astype("Int64")):
        bc = parse_cellid(cellid, library)
        if bc is None:
            continue
        try:
            out[bc] = int(nuclear)
        except (ValueError, TypeError):
            continue
    return out


# ---------------------------------------------------------------------------
# Binning + sampling
# ---------------------------------------------------------------------------

def assign_bin(total_reads: int) -> str | None:
    for label, lo, hi in DEFAULT_BINS:
        if hi is None:
            if total_reads >= lo:
                return label
        elif lo <= total_reads < hi:
            return label
    return None


def sample_panel(bin_to_bcs: dict[str, list[str]],
                 n_per_bin: int,
                 seed: int):
    rng = np.random.default_rng(seed)
    chosen: list[tuple[str, str]] = []  # (barcode, bin)
    bin_counts: dict[str, dict[str, int]] = {}
    for label, _, _ in DEFAULT_BINS:
        pool = bin_to_bcs.get(label, [])
        n_target = min(n_per_bin, len(pool))
        if n_target < n_per_bin:
            print(f"  [{label}] WARNING: pool has {len(pool)} < target {n_per_bin}",
                  file=sys.stderr)
        idx = rng.choice(len(pool), size=n_target, replace=False) if pool else []
        picked = [pool[i] for i in idx]
        chosen.extend((bc, label) for bc in picked)
        bin_counts[label] = {"available": len(pool), "drawn": n_target}
    return chosen, bin_counts


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_panel(panel_dir: Path,
                chosen: list[tuple[str, str]],
                totals: dict[str, int],
                seed: int,
                bin_counts: dict[str, dict[str, int]]):
    panel_dir.mkdir(parents=True, exist_ok=True)

    # bare barcodes
    with open(panel_dir / "barcodes.tsv", "w") as f:
        for bc, _ in chosen:
            f.write(bc + "\n")

    # barcode + total + bin
    with open(panel_dir / "barcodes_with_bin.tsv", "w") as f:
        f.write("barcode\ttotal_nuclear_reads\tdepth_bin\n")
        for bc, b in chosen:
            f.write(f"{bc}\t{totals[bc]}\t{b}\n")

    # per-genome BC strings (for samtools/pysam BAM filtering)
    with open(panel_dir / "barcodes_by_genome.tsv", "w") as f:
        f.write("genome\tlibrary\tbc_string\n")
        for genome, library in GENOMES:
            for bc, _ in chosen:
                f.write(f"{genome}\t{library}\tBC:Z:{bc}-{SAMPLE}_{library}\n")

    manifest = {
        "panel_dir": str(panel_dir),
        "seed": seed,
        "n_total": len(chosen),
        "bin_counts": bin_counts,
        "bins": [
            {"label": l, "lo": lo, "hi": hi} for l, lo, hi in DEFAULT_BINS
        ],
        "sample": SAMPLE,
        "genomes": [g for g, _ in GENOMES],
    }
    with open(panel_dir / "manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    bc_dir = Path(args.bc_counts_dir)
    out_root = Path(args.outdir)

    print(f"Loading bc_counts for {len(GENOMES)} genomes from {bc_dir} ...")
    totals: dict[str, int] = defaultdict(int)
    for genome, library in GENOMES:
        per = load_per_genome_counts(bc_dir, library)
        for bc, n in per.items():
            totals[bc] += n
        print(f"  {genome:<8} ({library:<8}): {len(per):>8,} barcodes")
    print(f"\nUnion of barcodes across genomes: {len(totals):,}")

    # Bin
    bin_to_bcs: dict[str, list[str]] = defaultdict(list)
    for bc, n in totals.items():
        b = assign_bin(n)
        if b is not None:
            bin_to_bcs[b].append(bc)

    print("\nDepth bin pool sizes:")
    for label, lo, hi in DEFAULT_BINS:
        hi_s = "inf" if hi is None else f"{hi}"
        print(f"  {label} [{lo:>6}, {hi_s:>6}): {len(bin_to_bcs.get(label, [])):>7,}")

    # Draw panels
    for spec in args.panels:
        if ":" not in spec:
            sys.exit(f"ERROR: bad --panels spec '{spec}', expected name:seed")
        name, seed_str = spec.split(":", 1)
        seed = int(seed_str)
        print(f"\nPanel {name} (seed={seed}):")
        chosen, bin_counts = sample_panel(bin_to_bcs, args.n_per_bin, seed)
        panel_dir = out_root / f"sub1k_{name}"
        write_panel(panel_dir, chosen, totals, seed, bin_counts)
        print(f"  Wrote {len(chosen)} barcodes -> {panel_dir}")
        for label, info in bin_counts.items():
            print(f"    {label}: drawn={info['drawn']}/{info['available']}")

    print("\nDone.")


if __name__ == "__main__":
    main()
