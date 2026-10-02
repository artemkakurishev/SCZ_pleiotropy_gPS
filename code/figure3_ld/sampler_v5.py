"""Instrumented copy of v4 sampler.py — v5 regeneration with RNG-free recording.

CRITICAL CONSTRAINT (amendment 1): the RNG call sequence must be IDENTICAL to v4's
sampler.py — same rng calls, same order, same arguments. All recording is pure
bookkeeping of already-drawn values (list appends, array writes). No recording
statement consumes or alters the RNG stream. Bit-identity to v4 pcts_b is verified
as a HALT condition by the driver (regen_v5.py).

Additions vs v4 (all RNG-free):
- draw_stratified_rec: returns per-bin chosen block codes and a borrow-event log
  (recipient_bin, donor_bin, n_genes) in addition to (gene idx, n_borrowed).
- run_block_null_rec: records per-set gene indices (padded int32), per-set draw
  attempts (for the redraw-cap report, item 2.5), and concatenates borrow logs.
"""
import time

import numpy as np

from sampler import build_sampler, _draw_run  # unchanged v4 builder + helper


def draw_stratified_rec(rng, S, tgt_freq):
    """Exact copy of v4 draw_stratified with recording appended.

    Returns (gene_idx, n_borrowed, chosen_per_bin, borrow_events) where
    chosen_per_bin[bi] = list of block codes chosen during bin bi's processing
    (own blocks first, then borrowed), and borrow_events = list of
    (recipient_bin_index, donor_bin_index, n_genes) tuples.
    """
    nb = len(S['bins_sorted'])
    order = [rng.permutation(bl) for bl in S['bin_blocks']]
    ptr = [0] * nb
    chosen = []
    chosen_per_bin = [[] for _ in range(nb)]
    borrow_events = []
    n_borrowed = 0
    for bi, g in enumerate(S['bins_sorted']):
        t = tgt_freq[g]
        if t <= 0:
            continue
        ch, cum, exhausted = _draw_run(rng, S, order, ptr, bi, t)
        chosen.extend(ch)
        chosen_per_bin[bi].extend(ch)
        deficit = t - cum
        if exhausted and deficit > 0:
            stop = False
            for bj in S['borrow_order'][bi]:
                while ptr[bj] < len(order[bj]):
                    c = order[bj][ptr[bj]]
                    sz = int(S['blk_sizes'][c])
                    if sz < deficit:
                        chosen.append(c); ptr[bj] += 1; n_borrowed += sz; deficit -= sz
                        chosen_per_bin[bi].append(c)
                        borrow_events.append((bi, bj, sz))
                    else:
                        if rng.random() < deficit / sz:
                            chosen.append(c); ptr[bj] += 1; n_borrowed += sz
                            chosen_per_bin[bi].append(c)
                            borrow_events.append((bi, bj, sz))
                        stop = True
                        break
                if stop:
                    break
    if not chosen:
        return np.array([], dtype=int), 0, chosen_per_bin, borrow_events
    cb = np.array(chosen)
    off, gf = S['offsets'], S['genes_flat']
    return np.concatenate([gf[off[c]:off[c + 1]] for c in cb]), n_borrowed, chosen_per_bin, borrow_events


def run_block_null_rec(seed, n_sets: int, tgt_freq: dict, target_n: int, S: dict,
                       risk_loggps: np.ndarray, compute_ks: bool = True,
                       label: str = '', progress_every: int = 2000):
    """Mirror of v4 run_block_null (same accept/reject logic, same rng stream via
    draw_stratified_rec) + RNG-free recording of gene IDs, attempts, borrow flow."""
    from scipy.stats import ks_2samp
    rng = np.random.default_rng(seed)
    n_axes = S['axis_mat'].shape[1]
    sizes = np.zeros(n_sets, int)
    means = np.zeros(n_sets)
    kss = np.zeros(n_sets)
    borrowed = np.zeros(n_sets, int)
    attempts = np.zeros(n_sets, int)
    pcts = np.zeros((n_sets, n_axes))
    max_size = target_n + 25
    gene_ids = np.full((n_sets, max_size), -1, dtype=np.int32)
    all_borrow_events = []   # list per set
    redraws = 0
    i = 0
    t0 = time.time()
    while i < n_sets:
        idx, nb_genes, cpb, bev = draw_stratified_rec(rng, S, tgt_freq)
        attempts[i] += 1
        if abs(len(idx) - target_n) > 25:
            redraws += 1
            if redraws > 50 * n_sets:
                raise RuntimeError(f'redraw storm: {redraws} redraws at accepted set {i}')
            continue
        lg = S['loggps_pool'][idx]
        sizes[i] = len(idx)
        means[i] = lg.mean()
        borrowed[i] = nb_genes
        gene_ids[i, :len(idx)] = idx
        all_borrow_events.append(bev)
        if compute_ks:
            kss[i] = ks_2samp(lg, risk_loggps).statistic
        pcts[i] = 100.0 * S['axis_mat'][idx].sum(axis=0) / len(idx)
        i += 1
        if progress_every and i % progress_every == 0:
            print(f'  [{label}] {i}/{n_sets} sets | {time.time() - t0:.0f}s | redraws {redraws}',
                  flush=True)
    return dict(sizes=sizes, means=means, kss=kss, borrowed=borrowed, attempts=attempts,
                pcts=pcts, gene_ids=gene_ids, borrow_events=all_borrow_events,
                redraws=redraws)
