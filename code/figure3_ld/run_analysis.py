#!/usr/bin/env python3
"""
run_analysis.py — FULL re-run of the LD-block stratified null analysis behind
Figure 3 A/B (PONE-D-26-36293 revision), from the raw inputs to the final tables,
followed by a discrepancy report against the frozen tables.

This driver is a faithful port of the originally executed analysis (executed
notebook of session sess_bdbf02ac8d9f, cells 0-2, 13, 32-33, 72-74, 87, 97-101,
107, 112-115, 122-131; drivers regen_v5.py / partb_v5.py). The statistical core
(sampler.py / sampler_v5.py) is copied VERBATIM from the frozen analysis package
(MD5-locked, see config.json). The gene-wise historical bootstrap reuses the
shipped repository script code/03_fig3_bootstrap_fdr5.py (byte-identical to the
script that produced the submitted results) via import — no re-implementation.

Chain:
  1. MD5 gates on every input (fail fast, no processing on mismatch).
  2. Cohort/pool/universe via the shipped script functions (byte-parity path):
     7,707-gene universe, 1,022 -> 1,016 risk set (6 shortfall-bin drops),
     3,605-gene background pool.
  3. LD-block attachment from data/figure3_ld/inputs/gene_block_map.tsv
     (dedup keep-max-Z; midpoint assignment was performed at package build time).
  4. Eight-axis L2G>=0.5 membership (5 primary + 3 exploratory), shipped-script
     logic extended to 8 axes; primary-5 sets asserted identical to shipped.
  5. Four null families:
       a. gene-wise gPS-matched bootstrap, B=1,000 (HISTORICAL / submitted spec;
          rng = default_rng(20260527), shipped code path, extended to 8 axes);
       b. PRIMARY: LD-block standard-gPS null, B=10,000, 8 axes
          (seed SeedSequence(20260527).spawn(6)[2]);
       c. LD-block size-only null (no gPS conditioning), B=10,000
          (seed spawn(16)[6]);
       d. exact gPS^(-d) nulls for the 5 primary axes, B=10,000 each
          (seeds spawn(16)[7..11]); per-gene counted disease sets from the L2G
          ledger at its native score floor, gPS^(-d) = gPS - |D_i ∩ axis_ids|
          clipped at 0, bins/targets rebuilt per axis.
  6. Final tables (written to results/figure3_ld/recomputed/):
       Table_S_block_vs_gene_bootstrap.tsv            (8 axes, gene-wise vs block)
       Table_final_five_axis_sensitivity.tsv          (5 axes x 4 specifications)
       Table_final_five_axis_exact_gps_minus_d.tsv    (5 axes, exact exclusion)
       Figure3B_main_source_data.tsv                  (5 axes, primary spec)
       null_axis_percentages_block_gps_matched.tsv    (10,000 x 8 null matrix)
  7. Comparison vs data/figure3_ld/frozen/ (SHA256-gated): gene-set identity,
     observed %, null mean/SD/q2.5/q97.5, b_le/b_ge, p, q per axis/spec.
     Tolerances follow the storage precision of each frozen table (4 dp for
     Table_S, 6 significant figures elsewhere; the null matrix is compared
     bit-for-bit). ALL discrepancies are written to
     results/figure3_ld/recomputed/comparison_report.md and the process exits
     non-zero. No parameter is ever adjusted to force a match.

Usage:
  python code/figure3_ld/run_analysis.py --full     # from the repository root
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import importlib.util
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as _st
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
CFG = json.loads((HERE / "config.json").read_text())

sys.path.insert(0, str(HERE))
from sampler import build_sampler, run_block_null            # noqa: E402
from sampler_v5 import run_block_null_rec                     # noqa: E402

# ----------------------------------------------------------------------------- constants
BASE_SEED = CFG["seed"]["base_seed"]
N_TARGET = CFG["cohort"]["n_target"]
SIZE_WINDOW = 25
B_GW = CFG["replicates"]["genewise"]
B_BLK = CFG["replicates"]["block"]

AXES8_IDS = {k: v for k, v in CFG["axes"]["disease_ids"].items()}
DISEASES8 = list(AXES8_IDS.keys())          # 5 primary (v4 order) + AF, Gout, OA
DISEASES5 = DISEASES8[:5]
FIVE = DISEASES5

DISPLAY = {"T2D": "Type 2 diabetes", "CAD": "Coronary artery disease",
           "Hypertension": "Hypertension", "Asthma": "Asthma",
           "Prostate Ca": "Prostate cancer", "Atrial Fib": "Atrial fibrillation",
           "Gout": "Gout", "OA": "Osteoarthritis"}

OUTDIR = REPO / "results" / "figure3_ld" / "recomputed"
FROZEN = REPO / "data" / "figure3_ld" / "frozen"


# ----------------------------------------------------------------------------- helpers
def file_md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def md5_gates() -> None:
    """Fail-fast MD5 verification of every input and of the verbatim code copies."""
    bad = []
    for group in ("repo_inputs", "ld_inputs", "code_md5"):
        for rel, want in CFG["inputs"][group].items():
            p = REPO / rel
            if not p.exists():
                bad.append(f"MISSING {rel}")
                continue
            got = file_md5(p)
            if got != want:
                bad.append(f"MD5 MISMATCH {rel}: expected {want}, got {got}")
            else:
                print(f"  MD5 OK  {rel}")
    if bad:
        for b in bad:
            print("ERROR:", b, file=sys.stderr)
        sys.exit("MD5 gates failed — halting before any computation.")


def import_shipped():
    """Import the shipped historical script (code/03_fig3_bootstrap_fdr5.py)."""
    spec = importlib.util.spec_from_file_location(
        "shipped", REPO / "code" / "03_fig3_bootstrap_fdr5.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def parse_ids(x):
    if pd.isna(x):
        return []
    s = str(x).strip()
    if s.startswith("["):
        try:
            return ast.literal_eval(s)
        except Exception:
            return []
    return [s]


def cp_count(k, n, alpha=0.05):
    """Clopper-Pearson interval on k/n (executed notebook, cells 112/115/122)."""
    lo = _st.beta.ppf(alpha / 2, k, n - k + 1) if k > 0 else 0.0
    hi = _st.beta.ppf(1 - alpha / 2, k + 1, n - k) if k < n else 1.0
    return lo, hi


def bh_q(pvec):
    """Benjamini-Hochberg q-values (executed notebook, cell 115)."""
    pvec = np.asarray(pvec, float)
    n = len(pvec)
    order = np.argsort(pvec)
    ranks = np.empty(n, int)
    ranks[order] = np.arange(1, n + 1)
    q = pvec * n / ranks
    qo = np.minimum.accumulate(q[order][::-1])[::-1]
    out = np.empty(n)
    out[order] = qo
    return np.minimum(out, 1.0)


def stats_from_pcts(null, obs, B):
    """Uniform R2#7 convention (executed notebook, cell 129)."""
    nm = null.mean()
    nsd = null.std(ddof=1)
    q_lo, q_hi = np.quantile(null, [0.025, 0.975])
    b_le = int((null <= obs).sum())
    b_ge = int((null >= obs).sum())
    p_two = min(1.0, 2 * min((b_le + 1) / (B + 1), (b_ge + 1) / (B + 1)))
    b_min = min(b_le, b_ge)
    pi_lo, pi_hi = cp_count(b_min, B)
    return dict(null_mean=nm, null_sd=nsd, emp_lo=q_lo, emp_hi=q_hi,
                resid=obs - nm, resid_lo=obs - q_hi, resid_hi=obs - q_lo,
                b_le=b_le, b_ge=b_ge, p_two=p_two,
                p_mc_lo=min(1.0, 2 * pi_lo), p_mc_hi=min(1.0, 2 * pi_hi))


def tol6(x):
    """Half of the 6th significant digit — storage tolerance of the frozen TSVs
    (executed notebook, cell 124)."""
    if x == 0:
        return 5e-7
    return 0.5 * 10 ** (math.floor(math.log10(abs(x))) - 5)


def build_realised_gps_table(gene_ids, sizes_b, pool_z1, scz_target):
    """Realised gene-level gPS matching diagnostic (39 risk-occupied bins).

    Transferred VERBATIM from the executed notebook (sess_bdbf02ac8d9f, cells 75-76).
    `gene_ids` is the (B, max_size) -1-padded array of realised pool-gene positional
    indices recorded by run_block_null_rec; `sizes_b` is the recorder's per-set sizes.
    No parameter is adjusted; rel_dev uses the RAW null mean (code wins over the
    data-dictionary note that wrongly described it as mean_minus_obs / null_mean).
    """
    pool = pool_z1.reset_index(drop=True)
    pool_gps_int = pool['gPS'].astype(int).to_numpy()

    risk = scz_target
    risk_gps_int = risk['gPS'].astype(int).to_numpy()
    tgt = pd.Series(risk_gps_int).value_counts().sort_index()
    bins = tgt.index.to_numpy()                      # 39 risk-occupied bins

    # FULL integer grid over pool support (pool may occupy bins the risk set does not)
    all_bins = np.arange(0, max(pool_gps_int.max(), bins.max()) + 1)

    B = gene_ids.shape[0]
    gps_flat = np.where(gene_ids >= 0, pool_gps_int[np.clip(gene_ids, 0, None)], -1)
    Hfull = np.zeros((B, len(all_bins)), dtype=np.int32)
    for b in range(B):
        row = gps_flat[b][gene_ids[b] >= 0]
        Hfull[b] = np.bincount(row, minlength=len(all_bins))[:len(all_bins)]

    assert np.all(Hfull.sum(1) == sizes_b), 'histogram totals != recorded sizes'

    risk_mask = np.isin(all_bins, bins)
    Hr = Hfull[:, risk_mask]                       # null counts on the 39 risk-occupied bins
    obs = tgt.to_numpy()

    rows = []
    for j, g in enumerate(bins):
        col = Hr[:, j]
        m = col.mean()
        rows.append(dict(gps_bin=int(g), observed=int(obs[j]),
                         null_mean=round(m, 3), null_sd=round(col.std(ddof=1), 3),
                         null_p2_5=np.percentile(col, 2.5), null_p97_5=np.percentile(col, 97.5),
                         null_min=col.min(), null_max=col.max(),
                         mean_minus_obs=round(m - obs[j], 3),
                         rel_dev=round((m - obs[j]) / obs[j], 4) if obs[j] > 0 else np.nan,
                         obs_within_95=bool(np.percentile(col, 2.5) <= obs[j] <= np.percentile(col, 97.5))))
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------- pipeline
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--full", action="store_true",
                    help="run the complete re-analysis (only mode)")
    args = ap.parse_args()
    if not args.full:
        ap.error("only --full is supported (fast figure-only mode lives in "
                 "plot_figure3.py)")

    t_start = time.time()
    OUTDIR.mkdir(parents=True, exist_ok=True)

    print("== 1. MD5 gates ==")
    md5_gates()

    shipped = import_shipped()

    # ------------------------------------------------------------------ setup
    print("\n== 2. Cohort / pool / universe (shipped code path) ==")
    gene_dedup, scz_primary, pool_z1 = shipped.build_cohort_and_pool()
    scz_target = shipped.derive_bootstrap_target(scz_primary, pool_z1)
    l2g_by = shipped.build_l2g_membership(shipped.INPUTS["l2g_diseases_full"]["path"])

    assert len(gene_dedup) == CFG["cohort"]["n_universe"], len(gene_dedup)
    assert len(pool_z1) == CFG["cohort"]["n_pool"], len(pool_z1)
    assert len(scz_primary) == CFG["cohort"]["n_risk_primary"], len(scz_primary)
    assert len(scz_target) == N_TARGET, len(scz_target)
    dropped = sorted(set(scz_primary.ensembl_id) - set(scz_target.ensembl_id))
    assert dropped == shipped.EXPECTED_DROPPED_IDS
    assert dropped == sorted(CFG["cohort"]["shortfall_dropped_ids"])
    print(f"  universe {len(gene_dedup)} | pool {len(pool_z1)} | "
          f"risk {len(scz_primary)} -> {len(scz_target)} (6 shortfall-bin drops match)")

    # ------------------------------------------------------------------ blocks
    print("\n== 3. LD-block attachment (gene_block_map.tsv, dedup keep-max-Z) ==")
    gbm = pd.read_csv(REPO / "data/figure3_ld/inputs/gene_block_map.tsv", sep="\t")
    gbm_d = gbm.sort_values("magma_z", ascending=False).drop_duplicates("ensembl_id",
                                                                        keep="first")
    flag_uni = set(gbm.loc[gbm["in_analysis_universe"] == True, "ensembl_id"])  # noqa: E712
    uni_ids = set(gene_dedup.ensembl_id)
    assert flag_uni == uni_ids, f"universe flag mismatch: {len(flag_uni ^ uni_ids)}"
    blk = gbm_d.set_index("ensembl_id")["block_id"]
    gene_dedup["block_id"] = gene_dedup["ensembl_id"].map(blk)
    scz_target["block_id"] = scz_target["ensembl_id"].map(blk)
    pool_z1["block_id"] = pool_z1["ensembl_id"].map(blk)
    assert gene_dedup["block_id"].notna().all()
    print(f"  cohort blocks {scz_target['block_id'].nunique()} | "
          f"pool blocks {pool_z1['block_id'].nunique()}")

    # ------------------------------------------------------- 8-axis membership
    print("\n== 4. Axis membership (L2G >= 0.5, 8 axes) ==")
    l2g = pd.read_csv(shipped.INPUTS["l2g_diseases_full"]["path"])
    l2g05 = l2g[l2g["score"] >= CFG["axes"]["l2g_threshold"]].copy()
    l2g05["pid"] = l2g05["diseaseIds"].apply(parse_ids)
    ex05 = l2g05[["geneId", "pid"]].explode("pid").rename(columns={"pid": "diseaseId"})
    axis_members8 = {d: set(ex05[ex05["diseaseId"].isin(ids)]["geneId"]) & uni_ids
                     for d, ids in AXES8_IDS.items()}
    for d in DISEASES5:  # primary-5 sets must equal the shipped-script sets
        assert axis_members8[d] == (l2g_by[d] & uni_ids), d
    axis_n8 = {d: len(s) for d, s in axis_members8.items()}
    for d, n in CFG["axes"]["axis_n_in_universe"].items():
        assert axis_n8[d] == n, (d, axis_n8[d], n)
    print("  axis sizes (in universe):", axis_n8)

    # observed percentages (8 axes)
    target_ids = set(scz_target["ensembl_id"])
    scz_obs8 = {d: 100.0 * len(target_ids & axis_members8[d]) / N_TARGET
                for d in DISEASES8}
    obs8 = [scz_obs8[d] for d in DISEASES8]
    assert abs(obs8[0] - 13.3858) < 1e-3
    print("  observed %:", {d: round(scz_obs8[d], 4) for d in DISEASES8})

    # ------------------------------------------------- gene-wise null (1k, 8x)
    print("\n== 5a. Gene-wise gPS-matched null, B=1,000 (historical spec) ==")
    # shipped 5-axis run (identical code path as the submitted results) ...
    scz_obs5, boot_pct5 = shipped.run_bootstrap(scz_target, pool_z1, l2g_by)
    # NOTE (provenance): the repository's committed historical
    # results/figure3_bootstrap_stats.tsv is the SUBMISSION-TIME gene-wise table.
    # It differs slightly from the final re-executed null (e.g. T2D CI lo
    # 0.0984 vs 0.1969, Prostate Ca bg_mean 6.0348 vs 6.0386, centered-p
    # convention vs R2#7). The authoritative gene-wise reference for the
    # revision is the frozen sensitivity table's submitted_genewise_gps_matched
    # rows; the null regenerated here matches those rows exactly (verified in
    # step 7). We therefore do NOT assert equality with the historical TSV.
    stats_check = shipped.compute_stats(scz_obs5, boot_pct5)
    print("  shipped 5-axis gene-wise null regenerated (authoritative comparison "
          "against the frozen sensitivity table happens in step 7)")

    # ... extended to 8 axes with the identical RNG stream (regen_v5 HALT-B path)
    target_gps_freq = scz_target["gPS"].astype(int).value_counts().to_dict()
    pool_by_gps = {int(g): pool_z1.loc[pool_z1["gPS"].astype(int) == int(g),
                                       "ensembl_id"].to_numpy()
                   for g in target_gps_freq}
    rng = np.random.default_rng(BASE_SEED)
    bins_sorted = sorted(target_gps_freq.keys())
    boot8 = np.zeros((B_GW, len(DISEASES8)))
    boot_ids = np.empty((B_GW, N_TARGET), dtype=object)
    for it in range(B_GW):
        sampled = []
        for g in bins_sorted:
            sampled.extend(rng.choice(pool_by_gps[g], size=target_gps_freq[g],
                                      replace=False))
        sampled_set = set(sampled)
        boot_ids[it] = np.array(sampled, dtype=object)
        for j, d in enumerate(DISEASES8):
            boot8[it, j] = 100.0 * len(sampled_set & axis_members8[d]) / N_TARGET
    assert np.array_equal(boot8[:, :5], boot_pct5), "8-axis extension drifted from shipped stream"
    print(f"  gene-wise null: {B_GW} sets x 8 axes; first 5 axes bit-identical to shipped")

    # --------------------------------------------- PRIMARY block null (10k, 8x)
    print("\n== 5b. PRIMARY LD-block standard-gPS null, B=10,000, 8 axes ==")
    S8 = build_sampler(pool_z1, scz_target, axis_members8, DISEASES8)
    risk_loggps = np.log2(scz_target["gPS"].to_numpy() + 1)
    ss = np.random.SeedSequence(BASE_SEED)
    step_names = ["step2_clumping", "step3_genewise_ks", "step3_block_10k",
                  "step3_block_1k", "step4_jackknife", "step6_diag"]
    SEEDS = dict(zip(step_names, ss.spawn(len(step_names))))
    SP = np.random.SeedSequence(BASE_SEED).spawn(16)
    assert SP[2].spawn_key == (2,) and SEEDS["step3_block_10k"].spawn_key == (2,)

    t0 = time.time()
    sizes_b, means_b, kss_b, pcts_b8, redraws_b, borrowed_b = run_block_null(
        SEEDS["step3_block_10k"], B_BLK, S8["target_freq"], N_TARGET, S8, risk_loggps,
        label="block10k")
    print(f"  10k null done in {time.time() - t0:.0f}s: sizes {sizes_b.mean():.1f}"
          f"+/-{sizes_b.std():.1f} [{sizes_b.min()},{sizes_b.max()}], "
          f"redraws {redraws_b}, borrowed/set {borrowed_b.mean():.1f}")

    # -------------------------------------- recorder pass (Amendment 1, gated) ----
    # run_block_null_rec is run as an ADDITIONAL pass with the SAME seed; the primary
    # path above (run_block_null) is left untouched so the verified 80/80 byte-identity
    # is never put at risk. The recorder returns identical null sets PLUS per-set
    # gene_ids. Before those gene_ids are used for anything we GATE on bit-identity to
    # the primary null (named check "recorder_bit_identity"); on any drift we halt and
    # report the magnitude — no parameter is adjusted.
    print("\n== 5b-rec. Recorder pass for the realised-gene-id diagnostic (gated) ==")
    t0 = time.time()
    rec_primary = run_block_null_rec(
        SEEDS["step3_block_10k"], B_BLK, S8["target_freq"], N_TARGET, S8, risk_loggps,
        compute_ks=True, label="block10k-rec", progress_every=5000)
    d_sizes = int(np.abs(rec_primary["sizes"] - sizes_b).max())
    d_means = float(np.abs(rec_primary["means"] - means_b).max())
    d_pcts = float(np.abs(rec_primary["pcts"] - pcts_b8).max())
    recorder_ok = (d_sizes == 0) and (d_means == 0.0) and (d_pcts == 0.0)
    RECORDER_CHECK = (
        "- recorder_bit_identity: additional `run_block_null_rec` pass vs primary "
        "`run_block_null` (same seed, spawn_key (2,)) — max|Δsizes| "
        f"{d_sizes}, max|Δmeans| {d_means:.3e}, max|Δpcts| {d_pcts:.3e} -> "
        + ("PASS (bit-identical; recorded gene_ids are the primary null's)"
           if recorder_ok else "**FAIL** (recorder drifted)"))
    print(f"  recorder pass done in {time.time() - t0:.0f}s | {RECORDER_CHECK}")
    if not recorder_ok:
        (OUTDIR / "comparison_report.md").write_text(
            "# Comparison: recomputed vs frozen LD-block tables\n\n"
            "## Recorder bit-identity gate (Amendment 1)\n" + RECORDER_CHECK + "\n\n"
            "The recorder pass drifted from the primary null, so its gene_ids cannot "
            "be trusted as the primary null's realised sets. Halting before the "
            "realised-gene-id table is built. No parameter was adjusted.\n")
        print("ERROR: recorder drifted from the primary null — see comparison_report.md. "
              "Halting; no parameter adjusted.", file=sys.stderr)
        return 1

    # ------------------------------------------------------ size-only (10k, 8x)
    print("\n== 5c. LD-block size-only null (no gPS conditioning), B=10,000 ==")
    S_size = dict(S8)
    S_size["bin_blocks"] = [np.arange(S8["B_pool"])]
    S_size["bins_sorted"] = [0]
    S_size["bin_logvals"] = np.array([0.0])
    S_size["borrow_order"] = [np.array([], dtype=int)]
    S_size["target_freq"] = {0: N_TARGET}
    S_size["empty_bins"] = []
    t0 = time.time()
    rec_size = run_block_null_rec(SP[6], B_BLK, {0: N_TARGET}, N_TARGET, S_size,
                                  risk_loggps, compute_ks=False, label="sizeonly",
                                  progress_every=5000)
    print(f"  size-only done in {time.time() - t0:.0f}s: sizes "
          f"{rec_size['sizes'].mean():.1f}+/-{rec_size['sizes'].std():.1f}, "
          f"redraws {rec_size['redraws']}, borrowed {int(rec_size['borrowed'].sum())} (expect 0)")

    # ------------------------------------------- exact gPS^(-d) nulls (5 axes)
    print("\n== 5d. Exact gPS^(-d) nulls, B=10,000 per primary axis ==")
    # per-gene counted disease sets from the L2G ledger at its native floor
    l2g_full = l2g.copy()
    l2g_full["disease_list"] = l2g_full["diseaseIds"].apply(parse_ids)
    expl = l2g_full[["geneId", "disease_list"]].explode("disease_list")
    gene_counted = expl.groupby("geneId")["disease_list"].apply(lambda s: set(s)).to_dict()

    # provenance gate: implemented gPS == |counted disease set| on pool+risk
    impl = pd.concat([pool_z1[["ensembl_id", "gPS"]], scz_target[["ensembl_id", "gPS"]]])
    mm = [g for g, gps in impl.itertuples(index=False) if len(gene_counted[g]) != gps]
    assert len(mm) == 0, f"implemented gPS != |D_gene| for {len(mm)} genes"
    print(f"  provenance gate: implemented gPS == |counted ledger set| for all "
          f"{len(impl)} pool+risk genes")

    exact_c = {ax: {g: len(ds & set(AXES8_IDS[ax])) for g, ds in gene_counted.items()}
               for ax in DISEASES8}

    # exact vs L2G-approx subtraction bookkeeping (notebook cells 97/99)
    oa_ids = set(AXES8_IDS["OA"])
    oa_count05 = (ex05[ex05["diseaseId"].isin(oa_ids)]
                  .groupby("geneId")["diseaseId"].nunique().to_dict())
    genes_pr = impl["ensembl_id"].to_numpy()
    is_risk = np.array([g in target_ids for g in genes_pr])
    ev_rows = []
    for ax in DISEASES8:
        if ax == "OA":
            approx = np.array([oa_count05.get(g, 0) for g in genes_pr])
        else:
            mem = axis_members8[ax]
            approx = np.array([1 if g in mem else 0 for g in genes_pr])
        exact = np.array([exact_c[ax][g] for g in genes_pr])
        diff = exact - approx
        ev_rows.append(dict(axis=ax, n_poolrisk_diff=int((diff != 0).sum()),
                            n_exact_gt=int((diff > 0).sum()),
                            n_exact_lt=int((diff < 0).sum()),
                            extra_subtracted_total=int(diff.sum()),
                            risk_n_diff=int((diff[is_risk] != 0).sum()),
                            risk_extra_subtracted=int(diff[is_risk].sum())))
    ev = pd.DataFrame(ev_rows)
    assert (ev["n_exact_lt"] == 0).all(), "exact must subtract >= approx"
    ev.to_csv(OUTDIR / "phase6_exact_vs_approx_contrib.tsv", sep="\t", index=False)

    def build_loo_inputs(c_pool, c_risk):
        pool_d = pool_z1.copy()
        pool_d["gPS"] = np.clip(pool_z1["gPS"].to_numpy() - c_pool, 0, None)
        risk_d = scz_target.copy()
        risk_d["gPS"] = np.clip(scz_target["gPS"].to_numpy() - c_risk, 0, None)
        S_d = build_sampler(pool_d, risk_d, axis_members8, DISEASES8)
        rld = np.log2(risk_d["gPS"].to_numpy() + 1)
        return pool_d, risk_d, S_d, rld

    def occ_from_ids(gene_ids, blk_pool):
        return np.array([np.unique(blk_pool[row[row >= 0]]).size for row in gene_ids])

    exact_runs = {}
    for i, ax in enumerate(FIVE):
        t0 = time.time()
        c_pool = pool_z1["ensembl_id"].map(exact_c[ax]).to_numpy()
        c_risk = scz_target["ensembl_id"].map(exact_c[ax]).to_numpy()
        pool_d, risk_d, S_d, rld = build_loo_inputs(c_pool, c_risk)
        assert sum(S_d["target_freq"].values()) == N_TARGET
        R = run_block_null_rec(SP[7 + i], B_BLK, S_d["target_freq"], N_TARGET, S_d, rld,
                               compute_ks=False, label=f"exact_{ax}",
                               progress_every=5000)
        occ = occ_from_ids(R["gene_ids"], pool_d["block_id"].to_numpy())
        exact_runs[ax] = dict(pcts=R["pcts"], sizes=R["sizes"], means=R["means"],
                              redraws=R["redraws"],
                              borrowed_mean=float(R["borrowed"].mean()), occ=occ,
                              n_bins=len(S_d["bins_sorted"]),
                              target_mean=float(rld.mean()),
                              gps_pool=pool_d["gPS"].to_numpy(),
                              gps_risk=risk_d["gPS"].to_numpy())
        print(f"  EXACT {ax}: sizes {R['sizes'].mean():.1f}+/-{R['sizes'].std():.1f}, "
              f"redraws {R['redraws']}, borrowed/set {exact_runs[ax]['borrowed_mean']:.1f}, "
              f"bins {exact_runs[ax]['n_bins']}, occ {occ.mean():.1f} "
              f"[{np.percentile(occ, 2.5):.0f},{np.percentile(occ, 97.5):.0f}], "
              f"{time.time() - t0:.0f}s", flush=True)
        del R

    # ================================================================ tables
    print("\n== 6. Assembling final tables ==")

    # ---- Table_S_block_vs_gene_bootstrap.tsv (8 axes; notebook cells 32-33) ----
    def stats_r27(null, obs):
        B = len(null)
        resid = obs - null.mean()
        lo = obs - np.percentile(null, 97.5)
        hi = obs - np.percentile(null, 2.5)
        p_lo = (int((null <= obs).sum()) + 1) / (B + 1)
        p_up = (int((null >= obs).sum()) + 1) / (B + 1)
        return resid, lo, hi, min(1.0, 2 * min(p_lo, p_up))

    def p_centered(null, obs):
        res = obs - null
        m = res.mean()
        c = res - m
        return (1.0 + np.sum(np.abs(c) >= abs(m))) / (len(null) + 1)

    rows = []
    for j, d in enumerate(DISEASES8):
        obs = scz_obs8[d]
        gr, glo, ghi, gp = stats_r27(boot8[:, j], obs)
        br, blo, bhi, bp = stats_r27(pcts_b8[:, j], obs)
        rows.append(dict(axis=d, primary=d in DISEASES5, n_L2G=axis_n8[d],
                         obs_pct=round(obs, 4),
                         gw_resid=gr, gw_lo=glo, gw_hi=ghi, gw_p=gp,
                         gw_p_ctr=p_centered(boot8[:, j], obs),
                         bw_resid=br, bw_lo=blo, bw_hi=bhi, bw_p=bp,
                         bw_p_ctr=p_centered(pcts_b8[:, j], obs),
                         gw_nullsd=boot8[:, j].std(), bw_nullsd=pcts_b8[:, j].std()))
    df8 = pd.DataFrame(rows)
    df8["sd_inflation"] = df8["bw_nullsd"] / df8["gw_nullsd"]
    for col, tag in [("gw_p", "gw"), ("gw_p_ctr", "gw_ctr"),
                     ("bw_p", "bw"), ("bw_p_ctr", "bw_ctr")]:
        q5 = np.full(8, np.nan)
        q5[:5] = multipletests(df8[col].values[:5], method="fdr_bh")[1]
        df8[f"{tag}_q_BH5"] = q5
        df8[f"{tag}_q_BH8"] = multipletests(df8[col].values, method="fdr_bh")[1]
    table_s = df8[["axis", "primary", "n_L2G", "obs_pct",
                   "gw_resid", "gw_lo", "gw_hi", "gw_p", "gw_q_BH5", "gw_q_BH8",
                   "gw_p_ctr", "gw_ctr_q_BH5", "gw_ctr_q_BH8",
                   "bw_resid", "bw_lo", "bw_hi", "bw_p", "bw_q_BH5", "bw_q_BH8",
                   "gw_nullsd", "bw_nullsd", "sd_inflation"]].copy()
    for c in table_s.columns:
        if table_s[c].dtype == float:
            table_s[c] = table_s[c].round(4)
    table_s.to_csv(OUTDIR / "Table_S_block_vs_gene_bootstrap.tsv", sep="\t", index=False)
    print("  Table_S_block_vs_gene_bootstrap.tsv (8 rows)")

    # ---- per-specification stats (notebook cells 129-130) ----
    gw_rows = {ax: stats_from_pcts(boot8[:, j], obs8[j], B_GW)
               for j, ax in enumerate(DISEASES5)}
    so_rows = {ax: stats_from_pcts(rec_size["pcts"][:, j], obs8[j], B_BLK)
               for j, ax in enumerate(DISEASES5)}
    std_rows = {ax: stats_from_pcts(pcts_b8[:, j], obs8[j], B_BLK)
                for j, ax in enumerate(DISEASES5)}

    # exact-run stats (notebook cells 122/124: ddof=1)
    verify_rows = []
    for i, ax in enumerate(FIVE):
        null = exact_runs[ax]["pcts"][:, i]   # axis i's own null under exact gPS^(-d)
        st = stats_from_pcts(null, obs8[i], B_BLK)
        occ = exact_runs[ax]["occ"]
        st.update(occ_mean=occ.mean(), occ_lo=np.quantile(occ, 0.025),
                  occ_hi=np.quantile(occ, 0.975), n_bins=exact_runs[ax]["n_bins"],
                  redraws=exact_runs[ax]["redraws"],
                  borrowed_mean=exact_runs[ax]["borrowed_mean"])
        verify_rows.append(st)

    # ---- Table_final_five_axis_sensitivity.tsv (5 x 4; cell 130) ----
    spec_map = [("submitted_genewise_gps_matched", gw_rows, B_GW,
                 "submitted specification (gene-wise resampling, B=1,000); superseded by the LD-block primary"),
                ("block_size_only", so_rows, B_BLK,
                 "alternative null specification: LD-block resampling with block-size matching only (no gPS conditioning)"),
                ("block_standard_gps", std_rows, B_BLK,
                 "PRIMARY: LD-block resampling with standard-gPS matching"),
                ("block_exact_gps_minus_d", {ax: verify_rows[i] for i, ax in enumerate(FIVE)},
                 B_BLK, "primary sensitivity: LD-block resampling with exact axis-excluded gPS^(-d) from the source ledger")]
    recs = []
    for spec, srows, Bspec, note in spec_map:
        ps = np.array([srows[ax]["p_two"] for ax in FIVE])
        qs = bh_q(ps)
        for j, ax in enumerate(FIVE):
            r = srows[ax]
            b_min = int(min(r["b_le"], r["b_ge"]))
            pi_lo, pi_hi = cp_count(b_min, Bspec)
            v_lo, v_hi = ps.copy(), ps.copy()
            v_lo[j] = min(1.0, 2 * pi_lo)
            v_hi[j] = min(1.0, 2 * pi_hi)
            recs.append(dict(axis=ax, specification=spec, B=Bspec,
                             observed_pct=obs8[DISEASES8.index(ax)],
                             null_mean=r["null_mean"], null_sd=r["null_sd"],
                             resid_pp=r["resid"], ref95_lo=r["resid_lo"],
                             ref95_hi=r["resid_hi"], b_le=r["b_le"], b_ge=r["b_ge"],
                             p_two=r["p_two"], p_mc_lo=r["p_mc_lo"], p_mc_hi=r["p_mc_hi"],
                             q_BH5=qs[j], q_BH5_mc_lo=bh_q(v_lo)[j],
                             q_BH5_mc_hi=bh_q(v_hi)[j],
                             floor_flag=("b_min=0: p at replicate floor" if b_min == 0 else ""),
                             note=note))
    sens = pd.DataFrame(recs)
    sens.to_csv(OUTDIR / "Table_final_five_axis_sensitivity.tsv", sep="\t", index=False,
                float_format="%.6g")
    print("  Table_final_five_axis_sensitivity.tsv (20 rows)")

    # ---- Table_final_five_axis_exact_gps_minus_d.tsv (5 rows; cells 123-125) ----
    gps_univ = gene_dedup.set_index("ensembl_id")["gPS"]
    risk_ids = target_ids
    pool_ids = set(pool_z1["ensembl_id"])
    sub_docs = []
    for i, ax in enumerate(FIVE):
        ids = AXES8_IDS[ax]
        contrib = exact_c[ax]
        changed_univ = [g for g in uni_ids if contrib.get(g, 0) > 0]
        changed_risk = [g for g in risk_ids if contrib.get(g, 0) > 0]
        changed_pool = [g for g in pool_ids if contrib.get(g, 0) > 0]
        amounts = pd.Series([contrib[g] for g in changed_univ])
        members = set(axis_members8[ax])
        member_not_counted = sum(1 for g in members if contrib.get(g, 0) == 0)
        before = gps_univ.copy()
        after = gps_univ - pd.Series({g: contrib.get(g, 0) for g in uni_ids}
                                     ).reindex(gps_univ.index).fillna(0)
        sub_docs.append(dict(
            axis=ax, disease_ids_removed="|".join(ids),
            n_genes_changed_universe=len(changed_univ),
            n_genes_changed_risk=len(changed_risk),
            n_genes_changed_pool=len(changed_pool),
            subtraction_amounts=sorted(amounts.unique().tolist()),
            member_not_counted=member_not_counted,
            gps_univ_mean_before=before.mean(), gps_univ_mean_after=after.mean(),
            gps_univ_sd_before=before.std(), gps_univ_sd_after=after.std(),
            gps_univ_median_before=before.median(), gps_univ_median_after=after.median(),
            gps_univ_max_before=int(before.max()), gps_univ_max_after=int(after.max())))
    subdoc5 = pd.DataFrame(sub_docs)
    assert all(subdoc5["member_not_counted"] == 0), "superset property violated"

    exact5 = pd.DataFrame(verify_rows)
    p5 = exact5["p_two"].values
    q5 = bh_q(p5)
    q5_mc_lo, q5_mc_hi = [], []
    for i in range(5):
        b_min = int(min(exact5.loc[i, "b_le"], exact5.loc[i, "b_ge"]))
        pi_lo, pi_hi = cp_count(b_min, B_BLK)
        p_lo, p_hi = min(1.0, 2 * pi_lo), min(1.0, 2 * pi_hi)
        v_lo, v_hi = p5.copy(), p5.copy()
        v_lo[i], v_hi[i] = p_lo, p_hi
        q5_mc_lo.append(bh_q(v_lo)[i])
        q5_mc_hi.append(bh_q(v_hi)[i])
    exact5["q_BH5"] = q5
    exact5["q_BH5_mc_lo"] = q5_mc_lo
    exact5["q_BH5_mc_hi"] = q5_mc_hi
    exact5["floor_flag"] = ["b_min=0: p at 10,000-replicate floor"
                            if min(r["b_le"], r["b_ge"]) == 0 else ""
                            for _, r in exact5.iterrows()]
    exact5["axis"] = FIVE

    m = exact5.merge(subdoc5, on="axis")
    match_rows = []
    for i, ax in enumerate(FIVE):
        dm = exact_runs[ax]["means"] - exact_runs[ax]["target_mean"]
        match_rows.append(dict(axis=ax,
                               match_maxabs_dmean=float(np.abs(dm).max()),
                               match_pct_within_005=float(100 * (np.abs(dm) < 0.05).mean())))
    m = m.merge(pd.DataFrame(match_rows), on="axis")
    m["subtraction_per_gene"] = m["subtraction_amounts"].apply(
        lambda a: "1 (all changed genes)" if a == [1] else str(a))
    ev5 = ev[ev["axis"].isin(FIVE)][["axis", "n_poolrisk_diff", "risk_n_diff"]].rename(
        columns={"n_poolrisk_diff": "n_exact_beyond_approx_poolrisk",
                 "risk_n_diff": "n_exact_beyond_approx_risk"})
    m = m.merge(ev5, on="axis")
    m = m.rename(columns={"observed": "observed_pct", "resid": "resid_pp",
                          "resid_lo": "ref95_lo", "resid_hi": "ref95_hi",
                          "n_genes_changed_universe": "n_genes_changed_universe_7707",
                          "n_genes_changed_risk": "n_genes_changed_risk_1016",
                          "n_genes_changed_pool": "n_genes_changed_pool_3605"})
    m["observed_pct"] = [obs8[DISEASES8.index(ax)] for ax in m["axis"]]
    final_cols = ["axis", "disease_ids_removed",
                  "n_genes_changed_universe_7707", "n_genes_changed_risk_1016",
                  "n_genes_changed_pool_3605", "subtraction_per_gene",
                  "n_exact_beyond_approx_poolrisk", "n_exact_beyond_approx_risk",
                  "gps_univ_mean_before", "gps_univ_mean_after",
                  "gps_univ_sd_before", "gps_univ_sd_after",
                  "gps_univ_median_before", "gps_univ_median_after",
                  "gps_univ_max_before", "gps_univ_max_after",
                  "observed_pct", "null_mean", "null_sd", "resid_pp",
                  "ref95_lo", "ref95_hi", "b_le", "b_ge", "p_two",
                  "p_mc_lo", "p_mc_hi", "q_BH5", "q_BH5_mc_lo", "q_BH5_mc_hi",
                  "floor_flag", "occ_mean", "occ_lo", "occ_hi",
                  "match_maxabs_dmean", "match_pct_within_005",
                  "n_bins", "redraws", "borrowed_mean"]
    m = m[final_cols]
    m.to_csv(OUTDIR / "Table_final_five_axis_exact_gps_minus_d.tsv", sep="\t",
             index=False, float_format="%.6g")
    print("  Table_final_five_axis_exact_gps_minus_d.tsv (5 rows)")

    # ---- Figure3B_main_source_data.tsv (5 rows, primary spec) ----
    # Provenance fidelity: the frozen Figure3B table was derived from the
    # *written* sensitivity table (%.6g storage), not from full-precision
    # intermediates. Re-read the file just written so the derived band columns
    # reproduce the frozen values exactly (e.g. CAD 2.85433 - (-2.03325) =
    # 4.88758, not the full-precision 4.88759).
    sens_stored = pd.read_csv(OUTDIR / "Table_final_five_axis_sensitivity.tsv",
                              sep="\t")
    prim = sens_stored[sens_stored.specification == "block_standard_gps"].set_index("axis")
    f3b_rows = []
    for ax in FIVE:
        f = prim.loc[ax]
        f3b_rows.append(dict(
            axis=ax, axis_display=DISPLAY[ax], observed_pct=f.observed_pct,
            null_mean_pct=f.null_mean, null_sd=f.null_sd,
            null_band_lo_q0_025=f.observed_pct - f.ref95_hi,
            null_band_hi_q0_975=f.observed_pct - f.ref95_lo,
            resid_pp=f.resid_pp, p_two=f.p_two, q_BH5=f.q_BH5,
            floor_flag=f.floor_flag))
    f3b = pd.DataFrame(f3b_rows)
    f3b_out = f3b.rename(columns={"null_band_lo_q0_025": "null_band_lo_q0.025",
                                  "null_band_hi_q0_975": "null_band_hi_q0.975"})
    # frozen file written with pandas default formatting (full double repr,
    # e.g. 7.586360000000001) — no float_format, unlike the %.6g tables
    f3b_out.to_csv(OUTDIR / "Figure3B_main_source_data.tsv", sep="\t", index=False)
    print("  Figure3B_main_source_data.tsv (5 rows)")

    # ---- null matrix (10,000 x 8) ----
    null_mat = pd.DataFrame(pcts_b8, columns=DISEASES8)
    null_mat.insert(0, "set_id", np.arange(len(null_mat)))
    null_mat.to_csv(OUTDIR / "null_axis_percentages_block_gps_matched.tsv", sep="\t",
                    index=False, float_format="%.6f")
    print("  null_axis_percentages_block_gps_matched.tsv (10,000 x 9)")

    # ---- Table_realised_genelevel_gps.tsv (39 occupied bins; cells 75-76) ----
    # Built from the recorder pass's gene_ids (bit-identity to the primary null gated
    # in step 5b-rec). Written with pandas default formatting (no float_format), exactly
    # as the frozen file was.
    perbin = build_realised_gps_table(rec_primary["gene_ids"], rec_primary["sizes"],
                                      pool_z1, scz_target)
    perbin.to_csv(OUTDIR / "Table_realised_genelevel_gps.tsv", sep="\t", index=False)
    n_out95 = int((~perbin["obs_within_95"]).sum())
    print(f"  Table_realised_genelevel_gps.tsv ({len(perbin)} rows; "
          f"{n_out95} bins outside the null 95% envelope)")

    # ============================================================ comparison
    print("\n== 7. Comparison vs frozen tables ==")
    report = compare_to_frozen(boot8, pcts_b8, rec_size, exact_runs, scz_obs8,
                               axis_members8, uni_ids, recorder_check=RECORDER_CHECK)

    # Historical-TVS provenance note: committed submission-time gene-wise table
    # vs the final re-executed gene-wise null (informational, not a gate).
    hist = pd.read_csv(REPO / "results" / "figure3_bootstrap_stats.tsv", sep="\t")
    hist_note = ["", "## Historical gene-wise TSV (informational, NOT a gate)",
                 "The repository's committed `results/figure3_bootstrap_stats.tsv` is the "
                 "submission-time table. Differences vs the final re-executed gene-wise "
                 "null (this run, matching the frozen sensitivity table exactly):"]
    for j, ax in enumerate(DISEASES5):
        null = boot8[:, j]
        obs = scz_obs8[ax]
        h = hist[hist["axis"] == ax].iloc[0]
        ci_lo_new = obs - np.quantile(null, 0.975)
        hist_note.append(
            f"- {ax}: bg_mean {h['bg_mean_pct']} -> {null.mean():.4f}; "
            f"CI lo {h['ci_lo']} -> {ci_lo_new:.4f}; "
            f"p_two_sided {h['p_two_sided']} (centered convention) vs "
            f"{min(1.0, 2 * min(int((null <= obs).sum()) + 1, int((null >= obs).sum()) + 1) / (B_GW + 1)):.6f} (R2#7)")
    hist_note.append("The frozen revision tables are authoritative; the historical TSV is "
                     "kept unchanged as the record of the submitted analysis.")
    report += "\n".join(hist_note)
    report_path = OUTDIR / "comparison_report.md"
    report_path.write_text(report)
    print(report)
    print(f"\nTotal runtime {time.time() - t_start:.0f}s. Report: {report_path}")

    n_fail = report.count("**FAIL**")
    if n_fail:
        print(f"\n{n_fail} comparison checks FAILED — see comparison_report.md. "
              f"No parameter was adjusted; investigate before use.", file=sys.stderr)
        return 1
    print("\nAll comparison checks passed.")
    return 0


# ----------------------------------------------------------------------------- comparison
def compare_to_frozen(boot8, pcts_b8, rec_size, exact_runs, scz_obs8, axis_members8,
                      uni_ids, recorder_check="") -> str:
    """Compare recomputed outputs against the SHA256-gated frozen tables.

    Tolerances: Table_S stored at 4 dp; other tables at 6 significant figures
    (tol6); the null matrix is compared bit-for-bit before formatting. The realised
    gene-level gPS table is compared by file sha256 (byte-identity) with a per-cell
    pass/fail tally.
    """
    lines = ["# Comparison: recomputed vs frozen LD-block tables", ""]
    if recorder_check:
        lines.append("## Recorder bit-identity gate (Amendment 1)")
        lines.append(recorder_check)
        lines.append("")
    lines.append("Frozen inputs verified against `data/figure3_ld/frozen/SHA256SUMS.txt` "
                 "before comparison.")
    ok = True

    # 0. frozen integrity
    import hashlib as hl
    sums = {}
    for ln in (FROZEN / "SHA256SUMS.txt").read_text().splitlines():
        h, fn = ln.split()
        sums[fn] = h
    for fn, h in sums.items():
        got = hl.sha256((FROZEN / fn).read_bytes()).hexdigest()
        if got != h:
            lines.append(f"- **FAIL** frozen file integrity: {fn}")
            ok = False
    lines.append("- frozen SHA256 integrity: PASS" if ok else "")

    def cmp_table(name, recomputed_path, float_cols, int_cols, tol_fun, extra_note=""):
        nonlocal ok
        fz = pd.read_csv(FROZEN / name, sep="\t")
        rc = pd.read_csv(recomputed_path, sep="\t")
        lines.append(f"\n## {name}")
        if list(fz.columns) != list(rc.columns):
            lines.append(f"- **FAIL** column mismatch: frozen {list(fz.columns)} vs "
                         f"recomputed {list(rc.columns)}")
            ok = False
            return
        if len(fz) != len(rc):
            lines.append(f"- **FAIL** row count {len(fz)} vs {len(rc)}")
            ok = False
            return
        worst = []
        for c in float_cols:
            d = np.abs(fz[c].to_numpy(float) - rc[c].to_numpy(float))
            tol = np.array([tol_fun(v) for v in fz[c].to_numpy(float)])
            if (d > tol).any():
                i = int(np.argmax(d - tol))
                worst.append(f"**FAIL** {c}: max|diff| {d.max():.3e} at row {i} "
                             f"(frozen {fz[c].iloc[i]}, recomputed {rc[c].iloc[i]})")
                ok = False
            else:
                worst.append(f"PASS {c}: max|diff| {d.max():.3e} within storage tolerance")
        for c in int_cols:
            neq = (fz[c].to_numpy() != rc[c].to_numpy())
            if neq.any():
                worst.append(f"**FAIL** {c}: {int(neq.sum())} unequal entries")
                ok = False
            else:
                worst.append(f"PASS {c}: identical")
        for c in fz.columns:
            if c in float_cols or c in int_cols:
                continue
            fz_s = fz[c].fillna("").astype(str).to_numpy()
            rc_s = rc[c].fillna("").astype(str).to_numpy()
            if (fz_s != rc_s).any():
                worst.append(f"**FAIL** text column {c}: differs")
                ok = False
        lines.extend("  - " + w for w in worst)
        if extra_note:
            lines.append(f"  - {extra_note}")

    # Table_S (4 dp storage)
    cmp_table("Table_S_block_vs_gene_bootstrap.tsv",
              OUTDIR / "Table_S_block_vs_gene_bootstrap.tsv",
              float_cols=["obs_pct", "gw_resid", "gw_lo", "gw_hi", "gw_p",
                          "gw_q_BH5", "gw_q_BH8", "gw_p_ctr", "gw_ctr_q_BH5",
                          "gw_ctr_q_BH8", "bw_resid", "bw_lo", "bw_hi", "bw_p",
                          "bw_q_BH5", "bw_q_BH8", "gw_nullsd", "bw_nullsd",
                          "sd_inflation"],
              int_cols=["n_L2G"], tol_fun=lambda v: 5e-5)

    # sensitivity (6 sig figs)
    cmp_table("Table_final_five_axis_sensitivity.tsv",
              OUTDIR / "Table_final_five_axis_sensitivity.tsv",
              float_cols=["observed_pct", "null_mean", "null_sd", "resid_pp",
                          "ref95_lo", "ref95_hi", "p_two", "p_mc_lo", "p_mc_hi",
                          "q_BH5", "q_BH5_mc_lo", "q_BH5_mc_hi"],
              int_cols=["B", "b_le", "b_ge"], tol_fun=tol6)

    # exact (6 sig figs)
    cmp_table("Table_final_five_axis_exact_gps_minus_d.tsv",
              OUTDIR / "Table_final_five_axis_exact_gps_minus_d.tsv",
              float_cols=["gps_univ_mean_before", "gps_univ_mean_after",
                          "gps_univ_sd_before", "gps_univ_sd_after",
                          "gps_univ_median_before", "gps_univ_median_after",
                          "observed_pct", "null_mean", "null_sd", "resid_pp",
                          "ref95_lo", "ref95_hi", "p_two", "p_mc_lo", "p_mc_hi",
                          "q_BH5", "q_BH5_mc_lo", "q_BH5_mc_hi",
                          "occ_mean", "occ_lo", "occ_hi",
                          "match_maxabs_dmean", "match_pct_within_005",
                          "borrowed_mean"],
              int_cols=["n_genes_changed_universe_7707", "n_genes_changed_risk_1016",
                        "n_genes_changed_pool_3605", "n_exact_beyond_approx_poolrisk",
                        "n_exact_beyond_approx_risk", "gps_univ_max_before",
                        "gps_univ_max_after", "b_le", "b_ge", "n_bins", "redraws"],
              tol_fun=tol6)

    # Figure3B source (6 sig figs)
    cmp_table("Figure3B_main_source_data.tsv",
              OUTDIR / "Figure3B_main_source_data.tsv",
              float_cols=["observed_pct", "null_mean_pct", "null_sd",
                          "null_band_lo_q0.025", "null_band_hi_q0.975", "resid_pp",
                          "p_two", "q_BH5"],
              int_cols=[], tol_fun=tol6)

    # null matrix: bit-for-bit before formatting
    lines.append("\n## null_axis_percentages_block_gps_matched.tsv")
    fz_null = pd.read_csv(FROZEN / "null_axis_percentages_block_gps_matched.tsv",
                          sep="\t")
    fz_vals = fz_null[DISEASES8].to_numpy()
    if fz_vals.shape != pcts_b8.shape:
        lines.append(f"- **FAIL** shape {fz_vals.shape} vs {pcts_b8.shape}")
        ok = False
    else:
        # frozen stored at %.6f -> round recomputed to the same storage
        # precision and require exact equality (bit-identical at stored
        # precision). A raw full-precision diff can reach exactly 5e-7 (half
        # the last stored digit) even when the underlying values are identical.
        diff = np.abs(fz_vals - np.round(pcts_b8, 6))
        bit = diff.max() == 0.0
        lines.append(f"- max|diff| after %.6f rounding: {diff.max():.3e} -> "
                     + ("PASS (bit-identical at stored precision)" if bit
                        else "**FAIL**"))
        ok = ok and bit

    # realised gene-level gPS table: file byte-identity (sha256) + per-cell tally.
    # The frozen file was written with pandas default formatting (no float_format),
    # so the standard is exact equality, not a storage tolerance.
    lines.append("\n## Table_realised_genelevel_gps.tsv")
    fz_real = FROZEN / "Table_realised_genelevel_gps.tsv"
    rc_real = OUTDIR / "Table_realised_genelevel_gps.tsv"
    fz_hash = hl.sha256(fz_real.read_bytes()).hexdigest()
    rc_hash = hl.sha256(rc_real.read_bytes()).hexdigest()
    byte_ok = fz_hash == rc_hash
    lines.append(f"- byte-identity sha256: recomputed `{rc_hash}` vs frozen "
                 f"`{fz_hash}` -> " + ("PASS" if byte_ok else "**FAIL**"))
    fz_r = pd.read_csv(fz_real, sep="\t", dtype=str)
    rc_r = pd.read_csv(rc_real, sep="\t", dtype=str)
    if list(fz_r.columns) == list(rc_r.columns) and len(fz_r) == len(rc_r):
        fz_v = fz_r.fillna("").to_numpy()
        rc_v = rc_r.fillna("").to_numpy()
        match = fz_v == rc_v
        n_pass = int(match.sum())
        n_tot = int(match.size)
        lines.append(f"- per-cell tally: {n_pass}/{n_tot} cells identical "
                     f"({len(fz_r)} rows x {len(fz_r.columns)} columns)")
        if not match.all():
            ok = False
            bad = np.argwhere(~match)
            for r, c in bad[:25]:
                lines.append(f"  - **FAIL** cell (row {r}, column "
                             f"'{fz_r.columns[c]}'): frozen {fz_v[r, c]!r} vs "
                             f"recomputed {rc_v[r, c]!r}")
            if len(bad) > 25:
                lines.append(f"  - ... and {len(bad) - 25} further differing cells")
    else:
        lines.append(f"- **FAIL** shape/columns differ: frozen {fz_r.shape} "
                     f"{list(fz_r.columns)} vs recomputed {rc_r.shape} "
                     f"{list(rc_r.columns)}")
        ok = False
    ok = ok and byte_ok

    # gene-set identity summary
    lines.append("\n## Gene sets")
    lines.append(f"- universe n={len(uni_ids)}; risk n={N_TARGET}; pool n=3,605 "
                 "(asserted during setup; HALT on mismatch)")
    lines.append("- axis membership sets asserted identical to the shipped-script "
                 "sets for the 5 primary axes (HALT on mismatch)")
    lines.append("- observed percentages: "
                 + ", ".join(f"{d} {scz_obs8[d]:.4f}" for d in DISEASES8))

    lines.append("\n## Conclusion")
    lines.append("All recomputed values match the frozen tables within storage "
                 "precision." if ok else
                 "DISCREPANCIES PRESENT — see FAIL lines above. No parameter was "
                 "adjusted to force agreement.")
    return "\n".join(lines)


if __name__ == "__main__":
    raise SystemExit(main())
