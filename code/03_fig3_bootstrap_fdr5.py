#!/usr/bin/env python3
"""
03_fig3_bootstrap_fdr5.py
=========================

End-to-end script that reproduces ``bootstrap_stats_5axis.tsv`` byte-identically
from the three documented raw inputs in ``inputs/``.

Method (one paragraph)
----------------------
1. Cohort: genes with MAGMA gene-level Z >= ``Z_THRESHOLD`` (3.09), after
   region-exclusion (chr X; MHC chr6:25-35 Mb; 17q21 chr17:43.5-44.9 Mb) and
   dedup-on-Ensembl keep-max-Z. 6 outlier-bin SCZ genes are dropped to make
   per-bin without-replacement sampling feasible. Emergent N = 1016.
2. Background pool: cohort universe with Z <= ``POOL_THRESHOLD`` (1.0).
3. L2G axis membership: a gene "belongs" to an axis if at least one credible-set
   row has ``score >= L2G_THRESHOLD`` (0.5) AND its parsed ``diseaseIds`` list
   intersects the axis EFO/MONDO ID set.
4. Bootstrap: B = ``N_ITER`` (1000) iterations. Each iteration draws, per gPS
   bin and without replacement, the same number of genes as the cohort has in
   that bin; per axis, computes bg%_b = #sampled members / N_TARGET * 100.
5. Stats per axis: residual_b = obs% - bg%_b; reports mean, 2.5/97.5 percentile
   CI, p_gt_zero, and the centered-permutation two-sided p:
       p_two = (1 + #{|residual_b - mean(residual)| >= |mean(residual)|}) / (B+1)
6. FDR: Benjamini-Hochberg computed across the 5 targeted clinical axes.

Reproducibility
---------------
* ``SEED = 20260527`` (NumPy ``default_rng``); deterministic sorted-bin iteration.
* MD5s of all three raw inputs are asserted before any data processing.
"""
from __future__ import annotations

import ast
import hashlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

# ─────────────────────────────────────────────────────────────────────────────
# Constants (declared a priori; do NOT tune to match the shipped TSV)
# ─────────────────────────────────────────────────────────────────────────────
SEED            = 20260527
N_ITER          = 1000
Z_THRESHOLD     = 3.09     # SCZ cohort MAGMA Z >= this
L2G_THRESHOLD   = 0.5      # L2G credible-set link >= this
POOL_THRESHOLD  = 1.0      # background pool MAGMA Z <= this

# 5 disease axes used in the supplementary figure
AXES = {
    'T2D':          {'MONDO_0005148'},
    'CAD':          {'EFO_0001645'},
    'Hypertension': {'EFO_0000537'},
    'Asthma':       {'MONDO_0004979'},
    'Prostate Ca':  {'EFO_0001663'},
}
DISEASES = list(AXES.keys())

# Region exclusions for MAGMA cohort
EXCL_MHC_CHR    = '6'
EXCL_MHC_START  = 25_000_000
EXCL_MHC_STOP   = 35_000_000
EXCL_17Q_CHR    = '17'
EXCL_17Q_START  = 43_500_000
EXCL_17Q_STOP   = 44_900_000

# Expected emergent cohort size (post bin-shortfall exclusion)
EXPECTED_N_TARGET = 1016

# Expected 6 dropped Ensembl IDs
EXPECTED_DROPPED_IDS = sorted([
    'ENSG00000183527',
    'ENSG00000132394',
    'ENSG00000196628',
    'ENSG00000112182',
    'ENSG00000187323',
    'ENSG00000138821',
])

# ─────────────────────────────────────────────────────────────────────────────
# Input paths and recorded MD5s
# ─────────────────────────────────────────────────────────────────────────────
HERE = Path(__file__).resolve().parent
INPUTS = {
    'gene_annotation_table': {
        'path': HERE.parent / 'inputs' / 'gene_annotation_table.tsv',
        'md5':  '6b308a7d300f9de33f6f458fdab9c79a',
    },
    'gps_source': {
        'path': HERE.parent / 'inputs' / 'disease_ta_index_pandas.csv',
        'md5':  '03b7c1e7fc211cd69c6d00141bff20b0',
    },
    'l2g_diseases_full': {
        'path': HERE.parent / 'inputs' / 'l2g_diseases_full.csv',
        'md5':  'f426e77b1bb71e4136d2577dc4dd6215',
    },
}

OUT_TSV = HERE.parent / 'results' / 'figure3_bootstrap_stats.tsv'

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────
def file_md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()

def assert_inputs() -> None:
    for name, meta in INPUTS.items():
        if not meta['path'].exists():
            sys.exit(f"ERROR: missing input file: {meta['path']}")
        got = file_md5(meta['path'])
        if got != meta['md5']:
            sys.exit(
                f"ERROR: MD5 mismatch for {meta['path'].name}:\n"
                f"  expected: {meta['md5']}\n  observed: {got}\n"
            )

def excl_regions(df: pd.DataFrame) -> pd.DataFrame:
    df = df[df['chr'].astype(str) != 'X'].copy()
    df = df[~((df['chr'].astype(str) == EXCL_MHC_CHR) &
              (df['start'] >= EXCL_MHC_START) &
              (df['stop']  <= EXCL_MHC_STOP))].copy()
    df = df[~((df['chr'].astype(str) == EXCL_17Q_CHR) &
              (df['start'] >= EXCL_17Q_START) &
              (df['stop']  <= EXCL_17Q_STOP))].copy()
    return df

def parse_disease_ids(x) -> list[str]:
    if pd.isna(x):
        return []
    s = str(x).strip()
    if s.startswith('['):
        try:
            return ast.literal_eval(s)
        except Exception:
            return []
    return [s]

def build_l2g_membership(l2g_path: Path) -> dict[str, set[str]]:
    l2g = pd.read_csv(l2g_path)
    l2g = l2g[l2g['score'] >= L2G_THRESHOLD].copy()
    l2g['diseaseIds_parsed'] = l2g['diseaseIds'].apply(parse_disease_ids)
    exploded = (
        l2g[['geneId', 'diseaseIds_parsed']]
        .explode('diseaseIds_parsed')
        .rename(columns={'diseaseIds_parsed': 'diseaseId'})
    )
    return {
        axis: set(exploded[exploded['diseaseId'].isin(ids)]['geneId'])
        for axis, ids in AXES.items()
    }

# ─────────────────────────────────────────────────────────────────────────────
# Stage 1 — cohort + pool
# ─────────────────────────────────────────────────────────────────────────────
def build_cohort_and_pool() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    gene_raw = pd.read_csv(INPUTS['gene_annotation_table']['path'], sep='\t')
    gps_src  = pd.read_csv(INPUTS['gps_source']['path'], usecols=['geneId', 'uniqueDiseases'])

    gene = (
        gene_raw.drop(columns=['gPS'], errors='ignore')
                .merge(gps_src.rename(columns={'geneId': 'ensembl_id', 'uniqueDiseases': 'gPS'}),
                       on='ensembl_id', how='left')
    )
    gene['gPS'] = gene['gPS'].fillna(0).astype(int)
    gene_clean = excl_regions(gene)
    
    gene_dedup = (gene_clean.sort_values('magma_z', ascending=False)
                            .drop_duplicates('ensembl_id', keep='first'))

    scz_primary = gene_dedup[gene_dedup['magma_z'] >= Z_THRESHOLD].copy()
    pool_z1     = gene_dedup[gene_dedup['magma_z'] <= POOL_THRESHOLD].copy()
    return gene_dedup, scz_primary, pool_z1

# ─────────────────────────────────────────────────────────────────────────────
# Stage 1b — bin-coverage pre-flight
# ─────────────────────────────────────────────────────────────────────────────
def derive_bootstrap_target(scz_primary: pd.DataFrame, pool_z1: pd.DataFrame) -> pd.DataFrame:
    target_need = scz_primary['gPS'].astype(int).value_counts().sort_index()
    pool_have   = pool_z1['gPS'].astype(int).value_counts()

    shortfall = []
    for gps_bin, need in target_need.items():
        have = int(pool_have.get(gps_bin, 0))
        if have < need:
            shortfall.append((int(gps_bin), int(need), have, need - have))

    to_drop = set()
    for gps_bin, _need, _have, short in shortfall:
        rows = (scz_primary[scz_primary['gPS'].astype(int) == gps_bin]
                .sort_values('magma_z', ascending=True))
        to_drop.update(rows.head(short)['ensembl_id'].tolist())

    scz_target = scz_primary[~scz_primary['ensembl_id'].isin(to_drop)].copy()
    assert len(scz_target) == EXPECTED_N_TARGET, f"Emergent N != {EXPECTED_N_TARGET}"
    return scz_target

# ─────────────────────────────────────────────────────────────────────────────
# Stage 2 & 3 — L2G membership & bootstrap
# ─────────────────────────────────────────────────────────────────────────────
def run_bootstrap(scz_target: pd.DataFrame, pool_z1: pd.DataFrame, 
                  l2g_by_disease: dict[str, set[str]]) -> tuple[dict[str, float], np.ndarray]:
    target_ids  = set(scz_target['ensembl_id'])
    n_target    = len(target_ids)

    scz_obs_pct = {d: 100.0 * len(target_ids & l2g_by_disease[d]) / n_target for d in DISEASES}

    target_gps_freq = scz_target['gPS'].astype(int).value_counts().to_dict()
    pool_by_gps = {
        int(g): pool_z1.loc[pool_z1['gPS'].astype(int) == int(g), 'ensembl_id'].to_numpy()
        for g in target_gps_freq
    }

    rng = np.random.default_rng(SEED)
    bins_sorted = sorted(target_gps_freq.keys())
    boot_pct = np.zeros((N_ITER, len(DISEASES)))

    for it in range(N_ITER):
        sampled: list[str] = []
        for g in bins_sorted:
            sampled.extend(rng.choice(pool_by_gps[g], size=target_gps_freq[g], replace=False))
        sampled_set = set(sampled)
        for j, d in enumerate(DISEASES):
            boot_pct[it, j] = 100.0 * len(sampled_set & l2g_by_disease[d]) / n_target

    return scz_obs_pct, boot_pct

# ─────────────────────────────────────────────────────────────────────────────
# Stage 4 — Stats, FDR, Verdict
# ─────────────────────────────────────────────────────────────────────────────
def compute_stats(scz_obs_pct: dict[str, float], boot_pct: np.ndarray) -> pd.DataFrame:
    scz_arr   = np.array([scz_obs_pct[d] for d in DISEASES])
    residuals = scz_arr[None, :] - boot_pct

    rows = []
    for j, d in enumerate(DISEASES):
        res      = residuals[:, j]
        mean_res = res.mean()
        ci_lo, ci_hi = np.percentile(res, [2.5, 97.5])
        p_gt_zero = (res > 0).mean()
        
        centered = res - res.mean()
        p_two    = (1.0 + np.sum(np.abs(centered) >= np.abs(mean_res))) / (N_ITER + 1)
        
        rows.append({
            'axis':           d,
            'scz_obs_pct':    round(scz_obs_pct[d], 4),
            'bg_mean_pct':    round(boot_pct[:, j].mean(), 4),
            'residual_mean':  round(mean_res, 4),
            'ci_lo':          round(ci_lo, 4),
            'ci_hi':          round(ci_hi, 4),
            'p_gt_zero':      round(p_gt_zero, 4),
            'p_two_sided':    round(p_two, 4),
        })
    stats_df = pd.DataFrame(rows)

    # Calculate BH-FDR directly across the 5 axes
    _, q_vals, _, _ = multipletests(stats_df['p_two_sided'].values, alpha=0.05, method='fdr_bh')
    stats_df['q_bh'] = np.round(q_vals, 4)

    def verdict(row) -> str:
        q = row['q_bh']
        if q < 0.05 and row['residual_mean'] < 0:
            return 'FDR depletion'
        if q < 0.05 and row['residual_mean'] > 0:
            return 'FDR excess'
        return 'borderline / null'

    stats_df['verdict'] = stats_df.apply(verdict, axis=1)
    return stats_df

# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────
def main() -> int:
    assert_inputs()
    _, scz_primary, pool_z1 = build_cohort_and_pool()
    scz_target = derive_bootstrap_target(scz_primary, pool_z1)
    l2g_by     = build_l2g_membership(INPUTS['l2g_diseases_full']['path'])
    
    scz_obs, boot_pct = run_bootstrap(scz_target, pool_z1, l2g_by)
    stats_df = compute_stats(scz_obs, boot_pct)

    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    stats_df.to_csv(OUT_TSV, sep='\t', index=False, float_format='%.4f')
    return 0

if __name__ == '__main__':
    raise SystemExit(main())