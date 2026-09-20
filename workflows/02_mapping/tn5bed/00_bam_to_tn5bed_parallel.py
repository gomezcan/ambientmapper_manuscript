#!/usr/bin/env python3
"""
Per-chrom parallel BAM -> tn5.bed.gz extractor.

Sidecar for ambientmapper `clean-bams`, which does not yet expose --threads.
Mirrors the behaviour of workflows/02_mapping/legacy_sm2_root1/1_6_scifi_makeTn5bed.py + the
sort -k1,1 -k2,2n | uniq | pigz post-step, but fans out one worker per chrom.

Usage:
  python 00_bam_to_tn5bed_parallel.py \
      --bam path/to/clean.bam \
      --out path/to/clean.tn5.bed.gz \
      --threads 10 --pigz-threads 2

Output is a single sorted+uniq'd BED with columns: chrom, cut, cut+1, BC, strand.
Within-chrom order is by position; cross-chrom order follows the BAM header SQ
list, which is the same order the original linear pipeline produced.
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pysam


def emit_bed_chrom(
    bam_path: str,
    chrom: str,
    out_dir: str,
    sort_buffer: str = "1G",
) -> tuple[str, str, int, int]:
    """Emit tn5 bed for one chrom, then dedup full records, position-ordered.

    Returns (chrom, sorted_path, n_total_reads, n_emitted).
    """
    raw_path = os.path.join(out_dir, f"{chrom}.raw.bed")
    sorted_path = os.path.join(out_dir, f"{chrom}.sorted.bed")
    n_total = 0
    n_emit = 0
    bam = pysam.AlignmentFile(bam_path, "rb")
    try:
        with open(raw_path, "w") as fh:
            for read in bam.fetch(chrom):
                n_total += 1
                if read.is_unmapped:
                    continue
                try:
                    bc = read.get_tag("BC")
                except KeyError:
                    bc = "NA"
                if read.is_reverse:
                    cut = read.reference_end - 5
                    strand = "-"
                else:
                    cut = read.reference_start + 4
                    strand = "+"
                fh.write(f"{read.reference_name}\t{cut}\t{cut+1}\t{bc}\t{strand}\n")
                n_emit += 1
    finally:
        bam.close()

    # Dedup PCR-equivalent tn5 records (same chrom + cut + barcode + strand),
    # keeping output position-ordered. The uniqueness key MUST include the
    # barcode (k4) and strand (k5): `sort -u -k1,1 -k2,2n` alone dedups on
    # (chrom, cut) only, which collapses every cell sharing an insertion site
    # to a single arbitrary barcode and silently discards the rest (~75% of
    # reads in accessible chromatin). end (k3) is deterministic (cut+1).
    subprocess.run(
        [
            "sort",
            "-S", sort_buffer,
            "--parallel=1",
            "-u",
            "-k1,1",
            "-k2,2n",
            "-k4,4",
            "-k5,5",
            raw_path,
            "-o", sorted_path,
        ],
        check=True,
    )
    os.remove(raw_path)
    return (chrom, sorted_path, n_total, n_emit)


def get_chrom_order(bam_path: str) -> list[str]:
    bam = pysam.AlignmentFile(bam_path, "rb")
    try:
        return list(bam.references)
    finally:
        bam.close()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--bam", required=True, help="Input BAM (must be indexed).")
    ap.add_argument("--out", required=True, help="Output tn5.bed.gz path.")
    ap.add_argument("--threads", type=int, default=10,
                    help="Number of parallel per-chrom workers (default 10).")
    ap.add_argument("--pigz-threads", type=int, default=2,
                    help="pigz worker count (default 2).")
    ap.add_argument("--scratch-dir", default=None,
                    help="Per-chrom temp dir root. Defaults to $TMPDIR or /tmp.")
    ap.add_argument("--sort-buffer", default="1G",
                    help="Per-worker sort buffer for `sort -S` (default 1G).")
    ap.add_argument("--skip-existing", action="store_true",
                    help="If --out exists and is non-empty, skip.")
    args = ap.parse_args()

    bam_path = os.path.abspath(args.bam)
    out_path = os.path.abspath(args.out)

    if args.skip_existing and os.path.exists(out_path) and os.path.getsize(out_path) > 0:
        print(f"[tn5-parallel] OUT exists ({os.path.getsize(out_path):,} B); --skip-existing -> done.")
        return 0

    if not Path(bam_path + ".bai").exists():
        print(f"[tn5-parallel] ERROR: BAM index missing: {bam_path}.bai", file=sys.stderr)
        return 1

    chroms = get_chrom_order(bam_path)
    print(f"[tn5-parallel] BAM        : {bam_path}")
    print(f"[tn5-parallel] OUT        : {out_path}")
    print(f"[tn5-parallel] threads    : {args.threads}")
    print(f"[tn5-parallel] pigz-thr   : {args.pigz_threads}")
    print(f"[tn5-parallel] sort-buf   : {args.sort_buffer}")
    print(f"[tn5-parallel] chroms     : {len(chroms)} from BAM header")
    if len(chroms) > 10:
        print(f"[tn5-parallel]            : first10={chroms[:10]} ... last2={chroms[-2:]}")
    else:
        print(f"[tn5-parallel]            : {chroms}")

    scratch_root = args.scratch_dir or os.environ.get("TMPDIR") or "/tmp"
    Path(scratch_root).mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    with tempfile.TemporaryDirectory(prefix="tn5bed_", dir=scratch_root) as work_dir:
        print(f"[tn5-parallel] scratch    : {work_dir}")
        chrom_to_path: dict[str, str] = {}
        chrom_to_stats: dict[str, tuple[int, int]] = {}
        with ProcessPoolExecutor(max_workers=args.threads) as ex:
            futures = {
                ex.submit(emit_bed_chrom, bam_path, c, work_dir, args.sort_buffer): c
                for c in chroms
            }
            for fut in as_completed(futures):
                c = futures[fut]
                try:
                    chrom, sorted_path, n_total, n_emit = fut.result()
                except Exception as e:
                    print(f"[tn5-parallel] ERROR chrom={c}: {e!r}", file=sys.stderr)
                    raise
                chrom_to_path[chrom] = sorted_path
                chrom_to_stats[chrom] = (n_total, n_emit)
                print(f"[tn5-parallel]  [{chrom}] total={n_total:,} emit={n_emit:,} -> {sorted_path}")

        # Concat per-chrom sorted BEDs in BAM SQ order, pipe through pigz.
        tmp_out = out_path + ".tmp"
        print(f"[tn5-parallel] concat -> pigz -> {tmp_out}")
        with open(tmp_out, "wb") as fh_out:
            pigz = subprocess.Popen(
                ["pigz", "-c", "-p", str(args.pigz_threads)],
                stdin=subprocess.PIPE,
                stdout=fh_out,
            )
            assert pigz.stdin is not None
            try:
                for chrom in chroms:
                    p = chrom_to_path.get(chrom)
                    if p is None:
                        continue
                    with open(p, "rb") as fh_in:
                        shutil.copyfileobj(fh_in, pigz.stdin, length=1 << 20)
            finally:
                pigz.stdin.close()
                rc = pigz.wait()
        if rc != 0:
            print(f"[tn5-parallel] ERROR pigz exit={rc}", file=sys.stderr)
            try:
                os.remove(tmp_out)
            except FileNotFoundError:
                pass
            return rc
        os.replace(tmp_out, out_path)

    elapsed = time.time() - t0
    total_total = sum(s[0] for s in chrom_to_stats.values())
    total_emit = sum(s[1] for s in chrom_to_stats.values())
    size = os.path.getsize(out_path)
    print(f"[tn5-parallel] DONE wall={elapsed:.1f}s "
          f"reads_total={total_total:,} reads_emitted={total_emit:,} "
          f"out_bytes={size:,}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
