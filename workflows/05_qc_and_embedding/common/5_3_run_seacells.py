#!/usr/bin/env python
"""
5_3_run_seacells.py -- Step 5.3: build meta-cells with SEACells from the Step-5.1b bundle.

Consumes the portable bundle written by 5_1b_export_for_seacells.R (MatrixMarket + TSVs), so the R
side never needs anndata and this side never needs R:
    <p>.mtx.gz    cells x ACRs, raw Tn5 insertion counts
    <p>.svd.tsv   the frozen Step-2 embedding      -> obsm['X_svd']
    <p>.obs.tsv   LouvainClusters + QC             -> obs
    <p>.var.tsv   ACR ids                          -> var

WHAT SEACELLS DOES AND DOES NOT DO FOR US
  It solves meta-cell CONSTRUCTION (archetypal analysis on diffusion components, built for scATAC
  sparsity). It takes ONE AnnData and therefore cannot know that Pre/wd/nd are meant to be COMPARED.
  Measured: At nd cells carry +37% median counts vs Pre, so at fixed n_SEACells the nd
  meta-cells come out deeper for free -- a stage-dependent bandwidth pointing the same way as the
  cleaning effect under study (the `nla` trap). RARE_TO is therefore applied HERE, after
  aggregation: every meta-cell is subsampled down to a common fragment budget.

PURITY IS SCORED ON CLUSTER ID, NOT CELL TYPE -- deliberately. The annotation is the open question
this layer exists to test, so scoring against it would be circular.

n_SEACells: SEACells' published default is one meta-cell per ~75 cells. Our own sweep found the
choice barely matters for maize (adjusted purity 0.876 -> 0.797 from k=10 to k=100) and that At
cannot support many meta-cells at all. Pass --cells-per-seacell to override.

NULL MODE (for the SVD-shuffle null; its engine and driver are not part of the shipped workflow).
  --shuffle-svd K   permute EACH COLUMN of X_svd independently (RNG seed K) before the kernel is
                    built. NOTE: THE KERNEL IS BUILT ON X_svd, NOT ON THE ACR MATRIX. X_svd is the
                    frozen Step-2 *tile* SVD from 5_1b; the ACR matrix only ever feeds aggregation.
                    Shuffling the ACR counts would therefore return an IDENTICAL partition and
                    report perfect reproducibility as an artifact. Shuffle the embedding.
  --membership-only skip aggregation/rarefaction/evaluate/h5ad; emit the partition only.
  Both default OFF and the deliverable path is unchanged when they are.

Usage:
  python 5_3_run_seacells.py --bundle <prefix> --out <outdir> [--cells-per-seacell 75]
                             [--rarefy-to N] [--seed 42] [--n-waypoint-eigs 10] [--n-neighbors 15]
                             [--shuffle-svd K] [--membership-only]
"""
import argparse, os, sys, gzip, warnings
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp
import anndata as ad

warnings.filterwarnings("ignore")


def log(msg):
    print(f"[5_3] {msg}", flush=True)


def load_bundle(prefix):
    mtx = prefix + ".mtx.gz"
    with gzip.open(mtx, "rb") as fh:
        X = sio.mmread(fh).tocsr()                      # cells x ACRs
    obs = pd.read_csv(prefix + ".obs.tsv", sep="\t")
    var = pd.read_csv(prefix + ".var.tsv", sep="\t")
    svd = pd.read_csv(prefix + ".svd.tsv", sep="\t")

    obs = obs.set_index("cellID")
    var = var.set_index("acrID")
    svd = svd.set_index("cellID")
    assert list(svd.index) == list(obs.index), "svd.tsv and obs.tsv row order differ"
    assert X.shape == (obs.shape[0], var.shape[0]), f"matrix {X.shape} vs obs/var {obs.shape[0]}x{var.shape[0]}"

    A = ad.AnnData(X=X, obs=obs, var=var)
    A.obsm["X_svd"] = svd.values.astype(np.float64)
    A.raw = A                                            # summarize_by_SEACell reads .raw by default
    A.obs["LouvainClusters"] = A.obs["LouvainClusters"].astype(str).astype("category")
    return A


def rarefy_rows(M, target, rng):
    """Subsample each row's counts down to `target` total. Rows already <= target are untouched.
    Multivariate-hypergeometric per row (sampling without replacement from the row's counts)."""
    M = sp.csr_matrix(M)
    out = M.copy().astype(np.int64)
    tot = np.asarray(M.sum(axis=1)).ravel()
    n_down = 0
    for i in range(M.shape[0]):
        t = tot[i]
        if t <= target:
            continue
        s, e = M.indptr[i], M.indptr[i + 1]
        counts = M.data[s:e].astype(np.int64)
        # draw `target` reads without replacement from this meta-cell's pooled reads
        picked = rng.multivariate_hypergeometric(counts, int(target))
        out.data[s:e] = picked
        n_down += 1
    out.eliminate_zeros()
    return out, n_down


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bundle", required=True, help="prefix written by 5_1b (no extension)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--cells-per-seacell", type=float, default=75.0,
                    help="SEACells default heuristic is 1 meta-cell per ~75 cells")
    ap.add_argument("--n-seacells", type=int, default=0, help="explicit count; overrides --cells-per-seacell")
    ap.add_argument("--rarefy-to", type=int, default=0,
                    help="common fragment budget per meta-cell; 0 = off. Pass the SAME value to every "
                         "stage of a comparison -- this is what SEACells cannot do for us.")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--n-waypoint-eigs", type=int, default=10)
    ap.add_argument("--n-neighbors", type=int, default=15)
    ap.add_argument("--max-iter", type=int, default=100)
    ap.add_argument("--min-iter", type=int, default=10)
    ap.add_argument("--shuffle-svd", type=int, default=0,
                    help="NULL MODE. >0 = permute EACH COLUMN of X_svd independently with this RNG "
                         "seed before building the kernel, destroying cell-cell structure while "
                         "preserving every PC's marginal distribution. The grouping that comes back "
                         "is what SEACells produces from structureless data = the floor. 0 = off. "
                         "Kept independent of --seed so the SEACells init can be held fixed while "
                         "only the shuffle varies.")
    ap.add_argument("--membership-only", action="store_true",
                    help="NULL MODE. Skip aggregation, rarefaction, evaluate, soft assignments and "
                         "the h5ad; write only cell_to_seacell.tsv + seacell_run.tsv. Null draws "
                         "need the partition and nothing else -- scoring happens in 4_3t.")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    name = os.path.basename(args.bundle)
    rng = np.random.default_rng(args.seed)
    np.random.seed(args.seed)

    A = load_bundle(args.bundle)

    # ---- NULL MODE: column-wise shuffle of the embedding ------------------------
    # THE KERNEL IS BUILT ON X_svd (see build_kernel_on below), NOT on the ACR matrix -- X_svd is
    # the FROZEN STEP-2 TILE SVD written by 5_1b. So randomizing the ACR counts would leave the
    # grouping BYTE-IDENTICAL and report perfect reproducibility as an artifact. The embedding is
    # the object that has to move.
    #   Column-wise (per-PC) permutation, not a row shuffle: a row shuffle relabels cells and
    # leaves the point cloud -- hence the partition -- unchanged up to naming.
    #   KNOWN CONSEQUENCE, measure it downstream: PC1 usually tracks depth in scATAC, so shuffling
    # it independently destroys the null's depth structure. 4_3t reports null-vs-real
    # sd(log10 depth) and conditions the rZ null on group SIZE for exactly this reason.
    if args.shuffle_svd > 0:
        srng = np.random.default_rng(args.shuffle_svd)
        Z = A.obsm["X_svd"].copy()
        for j in range(Z.shape[1]):
            Z[:, j] = Z[srng.permutation(Z.shape[0]), j]
        A.obsm["X_svd"] = Z
        log(f"NULL MODE: X_svd column-wise shuffled ({Z.shape[0]} cells x {Z.shape[1]} PCs), "
            f"shuffle_seed={args.shuffle_svd}. Grouping below is the STRUCTURELESS FLOOR.")
        if args.rarefy_to > 0:
            log("WARNING: --rarefy-to with --shuffle-svd is meaningless (rarefaction acts after "
                "assignment and never touches the rZ path). Ignoring is not automatic -- it will "
                "still run; prefer --membership-only for null draws.")
    n_sc = args.n_seacells if args.n_seacells > 0 else max(5, int(round(A.n_obs / args.cells_per_seacell)))
    log(f"{name}: {A.n_obs} cells x {A.n_vars} ACRs | n_SEACells={n_sc} "
        f"(~{A.n_obs/n_sc:.0f} cells each) | seed={args.seed}")
    if n_sc < 10:
        log(f"WARNING: only {n_sc} meta-cells. Our sweep showed <=20 meta-cells gives estimates whose "
            f"between-seed spread swamps any stage effect. Treat as a documented limit case.")

    from SEACells.core import SEACells as SC
    import SEACells.core as core
    import SEACells.evaluate as ev

    model = SC(A, build_kernel_on="X_svd", n_SEACells=n_sc, use_gpu=False, verbose=True,
               n_waypoint_eigs=args.n_waypoint_eigs, n_neighbors=args.n_neighbors)
    log("constructing kernel ...")
    model.construct_kernel_matrix()
    log("initializing archetypes ...")
    model.initialize_archetypes()
    # SEACells RAISES RuntimeWarning (it does not warn) when it reaches max_iter without converging,
    # which would kill the SLURM job rather than return the partial fit. The model state at that
    # point is still usable -- assignments exist, they just have not converged. So: retry once with
    # 3x the iterations, and if it still will not converge, proceed with a LOUD flag recorded in the
    # run TSV rather than losing the run. Never let this fail silently either way.
    log(f"fitting (max_iter={args.max_iter}) ...")
    converged = True
    try:
        model.fit(max_iter=args.max_iter, min_iter=args.min_iter)
    except RuntimeWarning:
        retry = args.max_iter * 3
        log(f"NOT CONVERGED at {args.max_iter} iters -- retrying with {retry}")
        try:
            model.fit(max_iter=retry, min_iter=args.min_iter)
            log(f"converged on retry ({retry} iters)")
        except RuntimeWarning:
            converged = False
            log(f"*** STILL NOT CONVERGED after {retry} iters -- proceeding with the partial fit. "
                f"Treat these meta-cells as provisional and re-run with --max-iter higher. ***")

    hard = model.get_hard_assignments()
    A.obs["SEACell"] = hard["SEACell"].astype(str).values
    n_real = A.obs["SEACell"].nunique()
    log(f"fit done: {n_real} non-empty meta-cells")

    # ---- membership-only short circuit (null draws) ------------------------------
    # Everything below aggregates, rarefies, evaluates and writes an h5ad. A null draw needs the
    # partition and nothing else: 4_3t re-aggregates the REAL perkb gene matrix over this
    # membership itself. Skipping saves the h5ad (per-draw storage) and the evaluate pass.
    if args.membership_only:
        A.obs[["SEACell", "LouvainClusters"]].to_csv(
            os.path.join(args.out, f"{name}.cell_to_seacell.tsv"), sep="\t")
        pd.DataFrame({
            "bundle": [args.bundle], "n_cells": [A.n_obs], "n_acrs": [A.n_vars],
            "n_seacells_requested": [n_sc], "n_seacells_realized": [n_real],
            "cells_per_seacell": [A.n_obs / max(n_real, 1)],
            "median_frags": [np.nan], "rarefy_to": [args.rarefy_to], "seed": [args.seed],
            "shuffle_svd": [args.shuffle_svd], "membership_only": [True],
            "converged": [converged],
        }).to_csv(os.path.join(args.out, f"{name}.seacell_run.tsv"), sep="\t", index=False)
        log(f"membership-only: wrote {name}.cell_to_seacell.tsv ({n_real} groups) -- done")
        return

    # ---- aggregate raw counts per meta-cell -------------------------------------
    agg = core.summarize_by_SEACell(A, SEACells_label="SEACell",
                                    celltype_label="LouvainClusters", summarize_layer="raw")
    tot = np.asarray(sp.csr_matrix(agg.X).sum(axis=1)).ravel()
    log(f"pooled fragments/meta-cell: median {np.median(tot):.0f}  "
        f"min {tot.min():.0f}  max {tot.max():.0f}")

    # ---- rarefaction: the stage-matching step SEACells cannot do -----------------
    if args.rarefy_to > 0:
        below = int((tot < args.rarefy_to).sum())
        Xr, n_down = rarefy_rows(agg.X, args.rarefy_to, rng)
        agg.layers["raw_counts"] = sp.csr_matrix(agg.X)
        agg.X = Xr
        log(f"rarefied to {args.rarefy_to}: {n_down} meta-cells downsampled, "
            f"{below} already below budget (left as-is -- report them)")
        agg.obs["below_budget"] = (tot < args.rarefy_to)
    agg.obs["n_frags_pooled"] = tot

    # ---- SEACells' own quality metrics ------------------------------------------
    # purity is scored on CLUSTER ID (see module docstring): the annotation is the open question.
    try:
        comp = ev.compactness(A, low_dim_embedding="X_svd", SEACells_label="SEACell")
        sep = ev.separation(A, low_dim_embedding="X_svd", nth_nbr=1, SEACells_label="SEACell")
        pur = ev.compute_celltype_purity(A, "LouvainClusters")
        qc = pd.concat([comp, sep, pur], axis=1)
        qc.to_csv(os.path.join(args.out, f"{name}.seacell_qc.tsv"), sep="\t")
        log(f"compactness median {float(np.nanmedian(comp.values)):.4f} (lower=tighter) | "
            f"separation median {float(np.nanmedian(sep.values)):.4f} (higher=better)")
        pcol = [c for c in pur.columns if "purity" in c.lower()]
        if pcol:
            log(f"cluster purity median {float(np.nanmedian(pur[pcol[0]].values)):.4f}")
    except Exception as e:
        log(f"WARNING: SEACells evaluate step failed ({type(e).__name__}: {e}) -- continuing")

    # ---- outputs -----------------------------------------------------------------
    A.obs[["SEACell", "LouvainClusters"]].to_csv(
        os.path.join(args.out, f"{name}.cell_to_seacell.tsv"), sep="\t")
    soft = model.get_soft_assignments()
    np.save(os.path.join(args.out, f"{name}.soft_assignments.npy"), np.asarray(soft[0]))
    agg.write_h5ad(os.path.join(args.out, f"{name}.metacells.h5ad"))
    pd.DataFrame({
        "bundle": [args.bundle], "n_cells": [A.n_obs], "n_acrs": [A.n_vars],
        "n_seacells_requested": [n_sc], "n_seacells_realized": [n_real],
        "cells_per_seacell": [A.n_obs / max(n_real, 1)],
        "median_frags": [float(np.median(tot))], "rarefy_to": [args.rarefy_to], "seed": [args.seed],
        "shuffle_svd": [args.shuffle_svd], "membership_only": [False],
        "converged": [converged],
    }).to_csv(os.path.join(args.out, f"{name}.seacell_run.tsv"), sep="\t", index=False)
    if not converged:
        log("REMINDER: converged=False in the run TSV -- do not report these meta-cells as final.")
    log(f"wrote -> {args.out}/{name}.{{metacells.h5ad, cell_to_seacell.tsv, seacell_qc.tsv, seacell_run.tsv}}")


if __name__ == "__main__":
    sys.exit(main())
