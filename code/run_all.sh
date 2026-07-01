#!/usr/bin/env bash
#
# run_all.sh — One-script-per-figure pipeline runner for the PeerJ deposit.
#
# Each python script:
#   * fails fast on missing or MD5-mismatched inputs,
#   * writes its figure(s) to ../results/figureN_*.pdf|svg|png and its
#     underlying numeric data to ../results/figureN_*_data.tsv (or, for
#     Figure 3, ../results/figure3_bootstrap_stats.tsv).
#
# Usage:
#   cd code
#   bash run_all.sh
#
# Requirements: python 3.11, the pinned package set in metadata/requirements.txt.
#
# Note: figure 2 (LOEUF × gPS hexbin) is currently MD5-gated and will exit 1
#       with an instruction message until the manuscript LOEUF table has been
#       staged at inputs/loeuf_gnomad_v41.tsv and its MD5 registered in both
#       code/02_fig2_hexbin_gradient.py and metadata/inputs_manifest.tsv.
#
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

PYTHON="${PYTHON:-python3}"

echo "[$(date +%H:%M:%S)] Figure 1 — Gene-property forest (9 traits, Model A)"
$PYTHON 01_fig1_gene_property.py

echo "[$(date +%H:%M:%S)] Figure 2 — LOEUF × gPS hexbin (SCZ − BD ΔZ)"
if ! $PYTHON 02_fig2_hexbin_gradient.py ; then
    echo "  WARNING: figure 2 skipped (LOEUF not yet provisioned). Continuing." >&2
fi

echo "[$(date +%H:%M:%S)] Figure 3 — gPS-matched bootstrap statistics"
$PYTHON 03_fig3_bootstrap_fdr5.py

echo "[$(date +%H:%M:%S)] Figure 3 — render bootstrap figure"
$PYTHON 03_fig3_make_figure.py

echo "[$(date +%H:%M:%S)] Figure 4 — SCHEMA targets across the gPS × Z landscape"
$PYTHON 04_fig4_schema_targets.py

echo "[$(date +%H:%M:%S)] Figure 5 — LOO gPS sensitivity (SCZ vs BD)"
$PYTHON 05_fig5_loo_sensitivity.py

echo "[$(date +%H:%M:%S)] Pipeline complete. Outputs in ../results/."
