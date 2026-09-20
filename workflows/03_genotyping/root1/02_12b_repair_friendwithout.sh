#!/usr/bin/env bash
#SBATCH --job-name=root1_repair_fwo
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=_logs/02_12b_repair_fwo_%A.log
#
# 02_12b_repair_friendwithout.sh — repair a partial relabel from 02_12a
#
# Background:
#   The α=0.05 assign of 02_12a completed (9023/9023 chunks), but its inline
#   `rescued → ambiguous` relabel step can crash mid-loop on a truncated gzip
#   member (EOFError: Compressed file ended before the end-of-stream marker),
#   leaving cell_map_ref_chunks_alpha005_friendwithout/ incomplete (7816/9023
#   files in the run that produced the manuscript). Phase 4 friend-OFF configs
#   (S04, S10) need all 9023.
#
# What this script does:
#   - AUDIT (always): inventory cell_map_ref_chunks/ and the friendwithout
#     sibling, identify missing relabels, test gzip integrity of missing
#     source files. Writes _audit_<timestamp>.tsv into the friendwithout dir.
#   - REPAIR (only if RUN_REPAIR=1): for each non-corrupt missing source
#     file, do the rescued→ambiguous relabel and write into friendwithout.
#     Skip corrupt files with logging. Atomic write via .tmp + rename.
#
# Usage:
#   # 1. audit only (default)
#   sbatch workflows/03_genotyping/root1/02_12b_repair_friendwithout.sh
#
#   # 2. after reviewing the audit TSV, repair
#   sbatch --export=ALL,RUN_REPAIR=1 workflows/03_genotyping/root1/02_12b_repair_friendwithout.sh
#
# =============================================================================

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"   # root of the data tree (3_Mapping/, 5_AmbientDetection/, ...)

cd "${PROJECT_ROOT}/5_AmbientDetection"

# conda activate <env from environment.yml>

SAMPLE=Root1_rep1
WORKDIR=${PROJECT_ROOT}/5_AmbientDetection/Root1_rep1
SRC_DIR=${WORKDIR}/cell_map_ref_chunks
DST_DIR=${WORKDIR}/cell_map_ref_chunks_alpha005_friendwithout

RUN_REPAIR=${RUN_REPAIR:-0}
TS=$(date +%Y%m%d_%H%M%S)

echo "============================================================"
echo "[$(date)] 02_12b repair friendwithout — sample=${SAMPLE}"
echo "  src    = ${SRC_DIR}"
echo "  dst    = ${DST_DIR}"
echo "  mode   = $([[ "${RUN_REPAIR}" == "1" ]] && echo REPAIR || echo AUDIT_ONLY)"
echo "============================================================"

[[ -d "${SRC_DIR}" ]] || { echo "ERROR: src dir not found: ${SRC_DIR}" >&2; exit 1; }
[[ -d "${DST_DIR}" ]] || { echo "ERROR: dst dir not found: ${DST_DIR}" >&2; exit 1; }

python - <<PYEOF
import gzip
import json
import os
import sys
import time
from glob import glob

SRC = "${SRC_DIR}"
DST = "${DST_DIR}"
RUN_REPAIR = int("${RUN_REPAIR}")
AUDIT_TSV = os.path.join(DST, "_audit_${TS}.tsv")

t0 = time.time()
src_files = sorted(glob(os.path.join(SRC, "*_filtered.tsv.gz")))
dst_basenames = {os.path.basename(p) for p in glob(os.path.join(DST, "*_filtered.tsv.gz"))}
missing = [p for p in src_files if os.path.basename(p) not in dst_basenames]

print(f"  src filtered files     : {len(src_files):>6,}")
print(f"  dst filtered files     : {len(dst_basenames):>6,}")
print(f"  missing in dst         : {len(missing):>6,}")
sys.stdout.flush()

# AUDIT — test gzip integrity by streaming each missing source file.
print("  audit: testing gzip integrity of missing source files...")
sys.stdout.flush()

ok_sources, corrupt_sources = [], []
for i, src_f in enumerate(missing, start=1):
    if i % 500 == 0:
        print(f"    audited {i:,}/{len(missing):,} ({time.time()-t0:.0f}s)")
        sys.stdout.flush()
    try:
        with gzip.open(src_f, "rt") as fh:
            for _ in fh:
                pass
        ok_sources.append(src_f)
    except (EOFError, OSError) as e:
        corrupt_sources.append((src_f, type(e).__name__, str(e)))

print(f"  audit complete         : ok={len(ok_sources):,}  corrupt={len(corrupt_sources):,}  ({time.time()-t0:.0f}s)")
sys.stdout.flush()

with open(AUDIT_TSV, "w") as fh:
    fh.write("status\tpath\terror_type\terror\n")
    for p in ok_sources:
        fh.write(f"ok_missing\t{p}\t\t\n")
    for p, et, e in corrupt_sources:
        fh.write(f"corrupt\t{p}\t{et}\t{e}\n")
print(f"  wrote audit            : {AUDIT_TSV}")
sys.stdout.flush()

if corrupt_sources:
    print("\n  CORRUPT SOURCE FILES:")
    for p, et, e in corrupt_sources[:20]:
        print(f"    {os.path.basename(p)}  [{et}: {e}]")
    if len(corrupt_sources) > 20:
        print(f"    ... and {len(corrupt_sources)-20} more (see audit TSV)")
    sys.stdout.flush()

if not RUN_REPAIR:
    print("\n  AUDIT_ONLY mode — no relabel performed.")
    print("  To run the repair, resubmit with: sbatch --export=ALL,RUN_REPAIR=1 ...")
    sys.exit(0)

# REPAIR — relabel ok_sources only. Skip corrupt with logging.
print(f"\n  REPAIR: relabeling {len(ok_sources):,} ok_missing files (skipping {len(corrupt_sources):,} corrupt)")
sys.stdout.flush()

n_done, n_relabeled, n_reads = 0, 0, 0
for i, src_f in enumerate(ok_sources, start=1):
    base = os.path.basename(src_f)
    out_f = os.path.join(DST, base)
    tmp_f = out_f + ".tmp"
    try:
        with gzip.open(src_f, "rt") as fin, gzip.open(tmp_f, "wt") as fout:
            header = fin.readline()
            fout.write(header)
            cols = header.rstrip("\n").split("\t")
            i_cls = cols.index("assigned_class")
            for line in fin:
                parts = line.rstrip("\n").split("\t")
                n_reads += 1
                if parts[i_cls] == "rescued":
                    parts[i_cls] = "ambiguous"
                    n_relabeled += 1
                fout.write("\t".join(parts) + "\n")
        os.rename(tmp_f, out_f)
        n_done += 1
    except Exception as e:
        if os.path.exists(tmp_f):
            os.unlink(tmp_f)
        print(f"    FAIL during repair: {base}  ({type(e).__name__}: {e})")
    if i % 200 == 0:
        print(f"    repaired {i:,}/{len(ok_sources):,} ({time.time()-t0:.0f}s)")
        sys.stdout.flush()

# Final accounting
final_dst = len({os.path.basename(p) for p in glob(os.path.join(DST, "*_filtered.tsv.gz"))})
print(f"\n  REPAIR done            : relabeled={n_done:,}  reads={n_reads:,}  rescued→ambiguous={n_relabeled:,}")
print(f"  final dst count        : {final_dst:,}/{len(src_files):,} (expected {len(src_files)-len(corrupt_sources)})")

summary = {
    "timestamp"          : "${TS}",
    "src_files"          : len(src_files),
    "dst_before"         : len(dst_basenames),
    "missing"            : len(missing),
    "ok_sources"         : len(ok_sources),
    "corrupt_sources"    : len(corrupt_sources),
    "repaired"           : n_done,
    "reads_seen"         : n_reads,
    "rescued_relabeled"  : n_relabeled,
    "dst_after"          : final_dst,
    "audit_tsv"          : AUDIT_TSV,
}
with open(os.path.join(DST, "_repair_summary_${TS}.json"), "w") as fh:
    json.dump(summary, fh, indent=2)
print(f"  wrote summary          : _repair_summary_${TS}.json")
PYEOF

echo ""
echo "============================================================"
echo "[$(date)] 02_12b done"
echo "============================================================"
