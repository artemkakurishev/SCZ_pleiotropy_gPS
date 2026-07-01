#!/usr/bin/env python3
"""
03_fig3_make_figure.py
======================

Render the polished supplementary FDR-5 forest figure from
``bootstrap_stats_5axis.tsv``. Simplified version without suptitles
and footnotes.
"""
from __future__ import annotations

import ast
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['font.family']  = ['Liberation Sans', 'Arimo', 'DejaVu Sans']
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['svg.hashsalt'] = 'supp_fdr5_2026'

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
TSV  = HERE.parent / 'results' / 'figure3_bootstrap_stats.tsv'
OUT_PNG = HERE.parent / 'results' / 'figure3_bootstrap_fdr5.png'
OUT_PDF = HERE.parent / 'results' / 'figure3_bootstrap_fdr5.pdf'
OUT_SVG = HERE.parent / 'results' / 'figure3_bootstrap_fdr5.svg'

CLINICAL_ORDER = ['T2D', 'CAD', 'Hypertension', 'Asthma', 'Prostate Ca']

AXIS_5_IDS = {
    'T2D':          {'MONDO_0005148'},
    'CAD':          {'EFO_0001645'},
    'Hypertension': {'EFO_0000537'},
    'Asthma':       {'MONDO_0004979'},
    'Prostate Ca':  {'EFO_0001663'},
}
L2G_THRESHOLD = 0.5

COLOR_BY_VERDICT = {
    'FDR excess':         '#1f78b4',
    'FDR depletion':      '#d73027',
    'borderline / null':  '#7f7f7f',
}

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

def compute_all_gene_background_pct() -> dict[str, float]:
    gene_raw = pd.read_csv(HERE.parent / 'inputs' / 'gene_annotation_table.tsv', sep='\t')
    gps_src  = pd.read_csv(HERE.parent / 'inputs' / 'disease_ta_index_pandas.csv',
                           usecols=['geneId', 'uniqueDiseases'])
    g = (gene_raw.drop(columns=['gPS'], errors='ignore')
                 .merge(gps_src.rename(columns={'geneId': 'ensembl_id',
                                                'uniqueDiseases': 'gPS'}),
                        on='ensembl_id', how='left'))
    g['gPS'] = g['gPS'].fillna(0).astype(int)
    mask = (
        (g['chr'].astype(str) != 'X') &
        ~((g['chr'].astype(str) == '6')  & (g['start'] >= 25_000_000) & (g['stop'] <= 35_000_000)) &
        ~((g['chr'].astype(str) == '17') & (g['start'] >= 43_500_000) & (g['stop'] <= 44_900_000))
    )
    g_dedup = (g[mask].sort_values('magma_z', ascending=False)
                       .drop_duplicates('ensembl_id', keep='first'))
    all_ids = set(g_dedup['ensembl_id'])
    n_all   = len(all_ids)

    l2g = pd.read_csv(HERE.parent / 'inputs' / 'l2g_diseases_full.csv')
    l2g = l2g[l2g['score'] >= L2G_THRESHOLD].copy()
    l2g['_dis'] = l2g['diseaseIds'].apply(parse_disease_ids)
    exploded = l2g[['geneId', '_dis']].explode('_dis').rename(columns={'_dis': 'diseaseId'})
    out: dict[str, float] = {}
    for axis, ids in AXIS_5_IDS.items():
        members = set(exploded[exploded['diseaseId'].isin(ids)]['geneId'])
        out[axis] = 100.0 * len(all_ids & members) / n_all
    return out

def render() -> None:
    stats_df = pd.read_csv(TSV, sep='\t')
    df = stats_df.set_index('axis').loc[CLINICAL_ORDER].reset_index()
    bg_all = compute_all_gene_background_pct()
    df['raw_delta'] = df.apply(lambda r: r['scz_obs_pct'] - bg_all[r['axis']], axis=1)

    XLIM = (-5.0, 6.5)

    fig, (axA, axB) = plt.subplots(
        1, 2, figsize=(12, 4.5), sharey=True, sharex=True,
        gridspec_kw={'width_ratios': [1, 1.15], 'wspace': 0.18},
    )
    y = np.arange(len(df))[::-1]

    for i, (_, r) in enumerate(df.iterrows()):
        col = COLOR_BY_VERDICT[r['verdict']]
        axA.barh(y[i], r['raw_delta'], color=col, edgecolor='none', height=0.65, alpha=0.95)
        
        bg_pct  = r['scz_obs_pct'] - r['raw_delta']
        scz_pct = r['scz_obs_pct']
        side_label = f"SCZ {scz_pct:.1f}%  BG {bg_pct:.1f}%"
        gap = 0.25
        if r['raw_delta'] >= 0:
            axA.text(r['raw_delta'] + gap, y[i], side_label, va='center', ha='left', fontsize=8.5, color='#222222')
        else:
            axA.text(r['raw_delta'] - gap, y[i], side_label, va='center', ha='right', fontsize=8.5, color='#222222')

    axA.axvline(0, color='black', linewidth=0.6, alpha=0.6)
    axA.set_yticks(y)
    axA.set_yticklabels(df['axis'].tolist(), fontweight='bold')
    axA.set_xlabel('SCZ% − all-gene background % (pp)', fontsize=10)
    axA.set_title('A   Raw enrichment (SCZ vs all-gene background)', fontsize=11, loc='left', pad=8, fontweight='bold')
    axA.spines['top'].set_visible(False)
    axA.spines['right'].set_visible(False)
    axA.tick_params(axis='both', labelsize=10)

    for i, (_, r) in enumerate(df.iterrows()):
        yi  = y[i]
        col = COLOR_BY_VERDICT[r['verdict']]
        err = [[r['residual_mean'] - r['ci_lo']], [r['ci_hi'] - r['residual_mean']]]
        axB.errorbar(r['residual_mean'], yi, xerr=err, fmt='o', color=col, ecolor=col,
                     capsize=4, elinewidth=1.6, markersize=8, markeredgecolor='black', markeredgewidth=0.6)
        
        # Обновлено: используем q_bh (соответствует обновленному bootstrap_fdr5.py)
        q = r['q_bh']
        qstr = f"q={q:.3f}" if q >= 0.001 else f"q={q:.4f}"
        label = f"{r['residual_mean']:+.2f} [{r['ci_lo']:+.2f}, {r['ci_hi']:+.2f}] {qstr}"
        
        axB.text(r['ci_hi'] + 0.20, yi, label, va='center', ha='left', fontsize=8.0, color=col)

    axB.axvline(0, color='black', linestyle='--', linewidth=0.8, alpha=0.6)
    axB.set_xlabel('gPS-matched residual (pp) ± 95% CI', fontsize=10)
    axB.set_title('B   Pleiotropy-adjusted residual + FDR verdict', fontsize=11, loc='left', pad=8, fontweight='bold')
    axB.spines['top'].set_visible(False)
    axB.spines['right'].set_visible(False)
    axB.tick_params(axis='both', labelsize=10)

    axA.set_xlim(*XLIM)
    axB.set_xlim(*XLIM)

    fig.subplots_adjust(left=0.10, right=0.97, top=0.90, bottom=0.18)

    handles = [
        mpatches.Patch(facecolor=COLOR_BY_VERDICT['FDR excess'], edgecolor='black', label='FDR excess'),
        mpatches.Patch(facecolor=COLOR_BY_VERDICT['FDR depletion'], edgecolor='black', label='FDR depletion'),
        mpatches.Patch(facecolor=COLOR_BY_VERDICT['borderline / null'], edgecolor='black', label='borderline / null'),
    ]
    fig.legend(handles=handles, loc='center', ncol=3, frameon=False, fontsize=9.5, bbox_to_anchor=(0.5, 0.06))

    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=300)
    fig.savefig(OUT_PDF)
    fig.savefig(OUT_SVG)
    plt.close(fig)

if __name__ == '__main__':
    render()