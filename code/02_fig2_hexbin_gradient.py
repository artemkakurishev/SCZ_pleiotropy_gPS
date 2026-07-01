#!/usr/bin/env python3
"""
02_fig2_hexbin_gradient.py
==========================

Render Figure 2 — the constraint × pleiotropy landscape: a hexbin map
showing mean ΔZ (Z_SCZ − Z_BIP) across LOEUF (gnomAD v4.1) on the x-axis
and log₂(gPS + 1) on the y-axis.

The plot is restricted to protein-coding genes with gPS ≥ 1 AND non-null
LOEUF, expected N ≈ 7,514 per the published manuscript.

Hexbin specification (FROZEN — matches published Figure 2)
----------------------------------------------------------
* gridsize = 22
* mincnt   = 5
* cmap     = 'RdBu_r' (red SCZ-enriched, blue BD-enriched)
* C        = Z_SCZ - Z_BIP (mean per bin via reduce_C_function=np.mean)
* symmetric vmin/vmax around 0 → vmax = max(|mean ΔZ|) clipped at 0.6
* Expected ~214 bins after mincnt filtering.

Inputs
------
* inputs/gene_annotation_table.tsv      — MD5 6b308a7d300f9de33f6f458fdab9c79a
* inputs/bip_gene_annotation_table.tsv  — MD5 c944322ea7c6a489a7229411162761c0
* inputs/loeuf_gnomad_v41.tsv           — MD5 58c78cec40fb791268aeda6df346fbf4

Outputs
-------
* results/figure2_hexbin_gradient.{pdf,svg,png}
* results/figure2_hexbin_gradient_data.tsv — the joined per-gene table
  actually rendered (gene_symbol, ensembl_id, gPS, LOEUF, Z_SCZ, Z_BIP, ΔZ).

Reproducibility notes
---------------------
* MD5 fail-fast on every input. The LOEUF table is the Ensembl-anchored
  MANE-select (with canonical fallback) view of gnomAD v4.1
  constraint_metrics, built from the public Broad GCS bucket; its MD5 is
  registered below.
* SCZ and BIP gene-level Z-scores are merged on `ensembl_id`.
* LOEUF is merged on ensembl_id (the primary key shared across all three
  inputs). The script does not currently use any symbol-based fallback;
  any LOEUF row whose ensembl_id is absent from the MAGMA tables is
  dropped by the inner merge.
"""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = ["Liberation Sans", "Arimo", "DejaVu Sans"]
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["svg.hashsalt"] = "fig2_hexbin_2026"

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ─────────────────────────────────────────────────────────────────────────────
# Paths
# ─────────────────────────────────────────────────────────────────────────────
HERE = Path(__file__).resolve().parent
INPUTS_DIR = HERE.parent / "inputs"
RESULTS_DIR = HERE.parent / "results"

LOEUF_EXPECTED_MD5 = "58c78cec40fb791268aeda6df346fbf4"

INPUTS = {
    "gene_annotation_scz": {
        "path": INPUTS_DIR / "gene_annotation_table.tsv",
        "md5": "6b308a7d300f9de33f6f458fdab9c79a",
    },
    "gene_annotation_bip": {
        "path": INPUTS_DIR / "bip_gene_annotation_table.tsv",
        "md5": "c944322ea7c6a489a7229411162761c0",
    },
    "loeuf_gnomad_v41": {
        "path": INPUTS_DIR / "loeuf_gnomad_v41.tsv",
        "md5": LOEUF_EXPECTED_MD5,
    },
}

OUT_PDF = RESULTS_DIR / "figure2_hexbin_gradient.pdf"
OUT_SVG = RESULTS_DIR / "figure2_hexbin_gradient.svg"
OUT_PNG = RESULTS_DIR / "figure2_hexbin_gradient.png"
OUT_TSV = RESULTS_DIR / "figure2_hexbin_gradient_data.tsv"

# ─────────────────────────────────────────────────────────────────────────────
# FROZEN published Figure 2 spec
# ─────────────────────────────────────────────────────────────────────────────
GPS_MIN = 1                     # Model A restriction: gPS >= 1
GRIDSIZE = 22
MINCNT = 5
CMAP = "RdBu_r"
VMAX_CLIP = 0.6                 # symmetric colorbar limit
EXPECTED_N_GENES = 7514         # per manuscript
EXPECTED_N_BINS_APPROX = 214    # per manuscript


# ─────────────────────────────────────────────────────────────────────────────
# MD5 gate
# ─────────────────────────────────────────────────────────────────────────────
def _file_md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def assert_inputs() -> None:
    for name, meta in INPUTS.items():
        p: Path = meta["path"]
        if not p.exists():
            sys.exit(f"ERROR: missing input file: {p}")
        got = _file_md5(p)
        if got != meta["md5"]:
            sys.exit(
                f"ERROR: MD5 mismatch for {p.name}:\n"
                f"  expected: {meta['md5']}\n  observed: {got}\n"
            )


# ─────────────────────────────────────────────────────────────────────────────
# Data assembly
# ─────────────────────────────────────────────────────────────────────────────
def _load_loeuf(path: Path) -> pd.DataFrame:
    """Load a LOEUF table. Accepts gnomAD v4.1 standard column layouts.

    Required columns (any of these names will be recognised):
      - gene ID:     'ensembl_id' | 'ensembl_gene_id' | 'gene_id' | 'transcript'
      - gene symbol: 'gene_symbol' | 'gene' | 'hgnc_symbol'  (optional, used as fallback)
      - LOEUF:       'loeuf' | 'lof.oe_ci.upper' | 'oe_lof_upper'

    The output DataFrame has columns: ensembl_id, gene_symbol (or NaN), LOEUF.
    """
    df = pd.read_csv(path, sep="\t" if path.suffix in (".tsv", ".txt") else ",")
    cols_lower = {c.lower(): c for c in df.columns}

    def _pick(candidates: list[str]) -> str | None:
        for c in candidates:
            if c.lower() in cols_lower:
                return cols_lower[c.lower()]
        return None

    ens_col = _pick([
        "ensembl_id", "ensembl_gene_id", "gene_id", "transcript", "transcript_id"
    ])
    sym_col = _pick(["gene_symbol", "gene", "hgnc_symbol", "symbol"])
    loeuf_col = _pick(["loeuf", "lof.oe_ci.upper", "oe_lof_upper"])

    if loeuf_col is None:
        sys.exit(
            f"ERROR: LOEUF column not found in {path.name}. "
            f"Looked for: loeuf, lof.oe_ci.upper, oe_lof_upper. "
            f"Columns present: {list(df.columns)}"
        )

    keep = {}
    if ens_col is not None:
        keep["ensembl_id"] = df[ens_col].astype(str).str.split(".").str[0]
    else:
        keep["ensembl_id"] = pd.NA
    if sym_col is not None:
        keep["gene_symbol"] = df[sym_col].astype(str)
    else:
        keep["gene_symbol"] = pd.NA
    keep["LOEUF"] = pd.to_numeric(df[loeuf_col], errors="coerce")

    out = pd.DataFrame(keep).dropna(subset=["LOEUF"])
    # If multiple rows per gene, keep the smallest LOEUF (most constrained transcript)
    if ens_col is not None and out["ensembl_id"].duplicated().any():
        out = out.sort_values("LOEUF").drop_duplicates("ensembl_id", keep="first")
    return out.reset_index(drop=True)


def _load_magma(path: Path, ensembl_col: str = "ensembl_id",
                z_col: str = "magma_z", gps_col: str = "gPS") -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    keep_cols = [ensembl_col, z_col, gps_col]
    if "gene_symbol" in df.columns:
        keep_cols.append("gene_symbol")
    elif "hgnc_symbol" in df.columns:
        keep_cols.append("hgnc_symbol")
    sub = df[[c for c in keep_cols if c in df.columns]].copy()
    sub = sub.rename(columns={ensembl_col: "ensembl_id"})
    sub["ensembl_id"] = sub["ensembl_id"].astype(str).str.split(".").str[0]
    return sub


def assemble_plot_frame() -> pd.DataFrame:
    df_scz = _load_magma(INPUTS["gene_annotation_scz"]["path"])
    df_bip = _load_magma(INPUTS["gene_annotation_bip"]["path"])
    df_loeuf = _load_loeuf(INPUTS["loeuf_gnomad_v41"]["path"])

    df_scz = df_scz.rename(columns={"magma_z": "Z_SCZ"})
    df_bip = df_bip.rename(columns={"magma_z": "Z_BIP"})

    # Merge SCZ + BIP on ensembl_id
    merged = pd.merge(
        df_scz[["ensembl_id", "Z_SCZ", "gPS"]],
        df_bip[["ensembl_id", "Z_BIP"]],
        on="ensembl_id", how="inner",
    )

    # Merge LOEUF (ensembl_id is the primary key)
    merged = pd.merge(
        merged,
        df_loeuf[["ensembl_id", "LOEUF", "gene_symbol"]],
        on="ensembl_id", how="inner",
    )

    # Restrict to Model A
    merged = merged[merged["gPS"] >= GPS_MIN].copy()
    merged = merged.dropna(subset=["Z_SCZ", "Z_BIP", "LOEUF", "gPS"])

    # Derived columns
    merged["log2_gPS1"] = np.log2(merged["gPS"] + 1.0)
    merged["deltaZ"] = merged["Z_SCZ"] - merged["Z_BIP"]

    return merged.reset_index(drop=True)


# ─────────────────────────────────────────────────────────────────────────────
# Render
# ─────────────────────────────────────────────────────────────────────────────
def render(merged: pd.DataFrame) -> None:
    n_genes = len(merged)
    # Soft check (the published N is 7,514; we don't fail-fast if the LOEUF
    # source differs by a few hundred genes, just warn).
    if abs(n_genes - EXPECTED_N_GENES) > 200:
        sys.stderr.write(
            f"WARNING: rendered gene count {n_genes} differs from expected "
            f"{EXPECTED_N_GENES} by more than 200. Check LOEUF source.\n"
        )

    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    out_cols = ["ensembl_id", "gene_symbol", "gPS", "log2_gPS1",
                "LOEUF", "Z_SCZ", "Z_BIP", "deltaZ"]
    merged[out_cols].to_csv(OUT_TSV, sep="\t", index=False, float_format="%.6g")

    # ── Figure ───────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(9.0, 6.5))
    fig.subplots_adjust(left=0.10, right=0.95, top=0.92, bottom=0.18)

    # Compute symmetric color limits from mean ΔZ per cell. We perform a dry
    # hexbin to find the per-bin mean magnitudes, then plot with VMAX_CLIP as
    # the ceiling so the colorbar is bounded as in the published figure.
    hb_dry = ax.hexbin(
        merged["LOEUF"], merged["log2_gPS1"],
        C=merged["deltaZ"], reduce_C_function=np.mean,
        gridsize=GRIDSIZE, mincnt=MINCNT, cmap=CMAP,
    )
    cell_means = hb_dry.get_array()
    vmax_data = float(np.nanmax(np.abs(cell_means)))
    vmax = min(vmax_data, VMAX_CLIP)
    hb_dry.remove()

    hb = ax.hexbin(
        merged["LOEUF"], merged["log2_gPS1"],
        C=merged["deltaZ"], reduce_C_function=np.mean,
        gridsize=GRIDSIZE, mincnt=MINCNT, cmap=CMAP,
        vmin=-vmax, vmax=vmax,
    )

    n_bins = (hb.get_array() != 0).sum()
    if abs(n_bins - EXPECTED_N_BINS_APPROX) > 30:
        sys.stderr.write(
            f"WARNING: rendered bin count {n_bins} differs from expected "
            f"{EXPECTED_N_BINS_APPROX} by more than 30. Check gridsize/mincnt.\n"
        )

    # Axes
    ax.set_xlabel("LOEUF (gnomAD v4.1)", fontsize=11.5, fontweight="bold")
    ax.set_ylabel(r"$\log_{2}(gPS + 1)$", fontsize=11.5, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", labelsize=10)

    # Titles, matching the published figure
    fig.text(
        0.525, 0.96,
        "SCZ–BD contrast is broadly uniform across constraint and pleiotropy",
        ha="center", va="center", fontsize=12.5, fontweight="bold",
    )
    fig.text(
        0.525, 0.92,
        "SCZ signal exceeds BD throughout; the constrained × pleiotropic corner is not a privileged regime",
        ha="center", va="center", fontsize=10, style="italic", color="#444444",
    )

    # Colorbar (horizontal, bottom)
    cbar_ax = fig.add_axes([0.18, 0.07, 0.66, 0.025])
    cbar = fig.colorbar(hb, cax=cbar_ax, orientation="horizontal")
    cbar.set_label("Mean ΔZ per hex (Z_SCZ − Z_BIP)", fontsize=10)
    cbar.ax.tick_params(labelsize=9)

    # Below-cbar caption: BD-enriched | SCZ-enriched anchors
    fig.text(0.18, 0.035, "BD-enriched", ha="left", va="center",
             fontsize=9, color="#0027C2", fontweight="bold")
    fig.text(0.84, 0.035, "SCZ-enriched", ha="right", va="center",
             fontsize=9, color="#9B0F0F", fontweight="bold")

    fig.savefig(OUT_PNG, dpi=300)
    fig.savefig(OUT_PDF)
    fig.savefig(OUT_SVG)
    plt.close(fig)


def main() -> int:
    assert_inputs()
    merged = assemble_plot_frame()
    render(merged)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
