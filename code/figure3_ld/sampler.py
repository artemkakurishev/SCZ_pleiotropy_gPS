"""Stratified block sampler for the LD-aware gPS-matched null (PONE-D-26-36293, R1#3).

Design (user-adjudicated Decision 2):
- Unit of resampling = LDetect EUR LD block; matching target = gPS (integer bins).
- Blocks binned by mean log2(gPS+1) of their pool genes, assigned to the nearest
  risk-set integer-gPS bin in log2(gPS+1) space.
- Per replicate: within each bin, shuffle blocks and take toward the bin's gene-count
  target (risk set's per-bin count) using UNBIASED RANDOMIZED ROUNDING at the crossing
  block: include it with probability (t - cum) / block_size, so E[count] = target.
- NEIGHBOUR BORROWING (documented deviation forced by block-level supply): mid/high
  gPS bins have fewer pool blocks than the risk target requires, because high-gPS
  low-Z genes sit in mixed blocks whose means bin low. When a bin's own blocks are
  exhausted, the unmet remainder is drawn from the nearest bin(s) in log2(gPS+1)
  space with remaining supply. Borrowed gene counts are recorded per set.
- Size window: accept if |total - target_n| <= 25, else redraw (size is uncorrelated
  with the axis statistic). NO KS or match-quality rejection (amendment b).
"""
import time

import numpy as np
import pandas as pd
from scipy.stats import ks_2samp


def build_sampler(pool_z1: pd.DataFrame, scz_target: pd.DataFrame,
                  axis_members: dict, diseases: list) -> dict:
    pool = pool_z1.reset_index(drop=True)
    pool_block_codes, pool_block_uniq = pd.factorize(pool['block_id'])
    B_pool = len(pool_block_uniq)
    blk_sizes = np.bincount(pool_block_codes)
    loggps_pool = np.log2(pool['gPS'].to_numpy() + 1)
    blk_mean_loggps = np.bincount(pool_block_codes, weights=loggps_pool) / blk_sizes
    genes_flat = np.argsort(pool_block_codes, kind='stable')
    offsets = np.zeros(B_pool + 1, dtype=int)
    offsets[1:] = np.cumsum(blk_sizes)

    target_freq = scz_target['gPS'].astype(int).value_counts().to_dict()
    bins_sorted = sorted(target_freq.keys())
    bin_logvals = np.log2(np.array(bins_sorted) + 1)
    assign = np.argmin(np.abs(blk_mean_loggps[:, None] - bin_logvals[None, :]), axis=1)
    bin_blocks = [np.where(assign == bi)[0] for bi in range(len(bins_sorted))]
    empty_bins = [g for bi, g in enumerate(bins_sorted) if len(bin_blocks[bi]) == 0]

    # static borrow order: for each bin, other bins sorted by log-space distance
    nb = len(bins_sorted)
    borrow_order = []
    for bi in range(nb):
        dist = np.abs(bin_logvals - bin_logvals[bi])
        dist[bi] = np.inf
        borrow_order.append(np.argsort(dist, kind='stable'))

    pool_eid = pool['ensembl_id'].to_numpy()
    axis_mat = np.column_stack([np.isin(pool_eid, list(axis_members[d])) for d in diseases])

    return dict(pool=pool, B_pool=B_pool, blk_sizes=blk_sizes, loggps_pool=loggps_pool,
                blk_mean_loggps=blk_mean_loggps, genes_flat=genes_flat, offsets=offsets,
                target_freq=target_freq, bins_sorted=bins_sorted, bin_logvals=bin_logvals,
                bin_blocks=bin_blocks, empty_bins=empty_bins,
                empty_shortfall=int(sum(target_freq[g] for g in empty_bins)),
                borrow_order=borrow_order,
                axis_mat=axis_mat, pool_block_uniq=pool_block_uniq)


def _draw_run(rng, S, order, ptr, bi, t):
    """Take blocks from bin bi's shuffled order toward gene target t (randomized rounding
    at the crossing block: include it w.p. (t - cum)/sz, so the run is UNBIASED around t —
    symmetric under/overshoot). Returns (chosen block codes, genes achieved, exhausted).
    `exhausted` is True only when the bin's block supply ran out with cum < t; rounding
    shortfalls at a rejected crossing block are NOT exhaustion and must not be topped up
    (topping those up while keeping accepted-crossing overshoots biases set size upward)."""
    chosen = []
    cum = 0
    blk_sizes = S['blk_sizes']
    obi = order[bi]
    while cum < t and ptr[bi] < len(obi):
        c = obi[ptr[bi]]
        sz = int(blk_sizes[c])
        if cum + sz < t:
            chosen.append(c); ptr[bi] += 1; cum += sz
        else:
            if rng.random() < (t - cum) / sz:
                chosen.append(c); ptr[bi] += 1; cum += sz
            break
    exhausted = bool(cum < t and ptr[bi] >= len(obi))
    return chosen, cum, exhausted


def draw_stratified(rng, S, tgt_freq):
    """One stratified block sample. Returns (gene index array, n genes borrowed across bins).
    Borrowing triggers only on supply exhaustion (structural deficit), from the nearest
    bin(s) in log2(gPS+1) space with remaining blocks."""
    nb = len(S['bins_sorted'])
    order = [rng.permutation(bl) for bl in S['bin_blocks']]
    ptr = [0] * nb
    chosen = []
    n_borrowed = 0
    for bi, g in enumerate(S['bins_sorted']):
        t = tgt_freq[g]
        if t <= 0:
            continue
        ch, cum, exhausted = _draw_run(rng, S, order, ptr, bi, t)
        chosen.extend(ch)
        deficit = t - cum
        if exhausted and deficit > 0:
            # Borrow from nearest bins with remaining supply. SINGLE crossing decision
            # across the whole borrow sequence (NO re-roll): take donor blocks while
            # sz < deficit, then ONE randomized-rounding crossing (include w.p.
            # deficit/sz), then STOP whether accepted or not. Re-rolling the crossing
            # across donor bins until acceptance would bias E[borrowed] above deficit.
            stop = False
            for bj in S['borrow_order'][bi]:
                while ptr[bj] < len(order[bj]):
                    c = order[bj][ptr[bj]]
                    sz = int(S['blk_sizes'][c])
                    if sz < deficit:
                        chosen.append(c); ptr[bj] += 1; n_borrowed += sz; deficit -= sz
                    else:
                        if rng.random() < deficit / sz:
                            chosen.append(c); ptr[bj] += 1; n_borrowed += sz
                        stop = True
                        break
                if stop:
                    break
    if not chosen:
        return np.array([], dtype=int), 0
    cb = np.array(chosen)
    off, gf = S['offsets'], S['genes_flat']
    return np.concatenate([gf[off[c]:off[c + 1]] for c in cb]), n_borrowed


def run_block_null(seed, n_sets: int, tgt_freq: dict, target_n: int, S: dict,
                   risk_loggps: np.ndarray, compute_ks: bool = True,
                   label: str = '', progress_every: int = 2000):
    rng = np.random.default_rng(seed)
    n_axes = S['axis_mat'].shape[1]
    sizes = np.zeros(n_sets, int)
    means = np.zeros(n_sets)
    kss = np.zeros(n_sets)
    borrowed = np.zeros(n_sets, int)
    pcts = np.zeros((n_sets, n_axes))
    redraws = 0
    i = 0
    t0 = time.time()
    while i < n_sets:
        idx, nb_genes = draw_stratified(rng, S, tgt_freq)
        if abs(len(idx) - target_n) > 25:
            redraws += 1
            if redraws > 50 * n_sets:
                raise RuntimeError(f'redraw storm: {redraws} redraws at accepted set {i}')
            continue
        lg = S['loggps_pool'][idx]
        sizes[i] = len(idx)
        means[i] = lg.mean()
        borrowed[i] = nb_genes
        if compute_ks:
            kss[i] = ks_2samp(lg, risk_loggps).statistic
        pcts[i] = 100.0 * S['axis_mat'][idx].sum(axis=0) / len(idx)
        i += 1
        if progress_every and i % progress_every == 0:
            print(f'  [{label}] {i}/{n_sets} sets | {time.time() - t0:.0f}s | redraws {redraws}',
                  flush=True)
    return sizes, means, kss, pcts, redraws, borrowed
