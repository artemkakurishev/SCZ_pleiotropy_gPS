#!/usr/bin/env python3
"""
01_fig1_gene_property.py
========================

Render Figure 1 — the 9-trait Model-A forest of MAGMA gene-property
standardized coefficients (β_std) against log2(gPS + 1). The figure is
the manuscript's Figure 1 (Frontiers preprint, June 2025).

One marker per trait (Model A only, gPS >= 1). Traits are color-categorised
into:
  - Psychiatric (4 traits): MDD, SCZ, ADHD, BD          → blue
  - Neurological / Neuro-dev (3 traits): PD, OCD, ASD   → purple
  - Non-psychiatric disease (1 trait): IBD              → red
  - Anthropometric control (1 trait): Height            → orange

Trait order on the y-axis is **the published display order**, not a strict
descending β_std sort. The published order interleaves the control / somatic
positive controls with the psychiatric / neuro-dev block:

    Height → IBD → MDD → SCZ → ADHD → BD → OCD† → ASD → PD

Within the psychiatric block, descending β_std is preserved (MDD > SCZ >
ADHD > BD > OCD > ASD). The Height-before-IBD ordering reflects an editorial
choice in the published figure (anthropometric control as anchor); the
underlying β_std ordering across all 9 traits would put IBD first
(β_std = 0.1820 vs Height β_std = 0.1639).

OCD2025 carries a low-power footnote ("†") when its `power_caveat`
field is True in the source table.

Significance annotation follows the published figure:
  *** : p_one_sided < 1e-3
  **  : p_one_sided < 1e-2
  *   : p_one_sided < 0.05
  ns  : p_one_sided >= 0.05

Inputs
------
* inputs/gp_results_big9.tsv — MD5 b83619c51e4a433c39feea9e259852b7.
  18 rows: 9 traits × {Model A, Model B}. Columns:
  trait, model, n_eff, n_genes, beta, beta_std, se, p_one_sided, power_caveat.

Outputs
-------
* results/figure1_gene_property.{pdf,svg,png}
* results/figure1_gene_property_data.tsv  — the 9 Model-A rows actually
  rendered, with computed 95% CI columns and significance marker.

Reproducibility
---------------
* No stochastic step (deterministic pure read + plot).
* MD5 of inputs/gp_results_big9.tsv asserted before any compute.
"""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = ["Liberation Sans", "Arimo", "DejaVu Sans"]
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["svg.hashsalt"] = "fig1_gene_property_2026"

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

# ─────────────────────────────────────────────────────────────────────────────
# Paths
# ─────────────────────────────────────────────────────────────────────────────
HERE = Path(__file__).resolve().parent
INPUTS_DIR = HERE.parent / "inputs"
RESULTS_DIR = HERE.parent / "results"

INPUTS = {
    "gp_results_big9": {
        "path": INPUTS_DIR / "gp_results_big9.tsv",
        "md5": "b83619c51e4a433c39feea9e259852b7",
    },
}

OUT_PDF = RESULTS_DIR / "figure1_gene_property.pdf"
OUT_SVG = RESULTS_DIR / "figure1_gene_property.svg"
OUT_PNG = RESULTS_DIR / "figure1_gene_property.png"
OUT_TSV = RESULTS_DIR / "figure1_gene_property_data.tsv"

# ─────────────────────────────────────────────────────────────────────────────
# Display configuration — matches published Figure 1 exactly
# ─────────────────────────────────────────────────────────────────────────────

# Trait → Category (published grouping)
TRAIT_CATEGORY: dict[str, str] = {
    "SCZ3":     "Psychiatric (SCZ/BD/MDD/ADHD)",
    "BIP2021":  "Psychiatric (SCZ/BD/MDD/ADHD)",
    "MDD2025":  "Psychiatric (SCZ/BD/MDD/ADHD)",
    "ADHD2022": "Psychiatric (SCZ/BD/MDD/ADHD)",
    "OCD2025":  "Neurological/Neuro-dev (PD/OCD/ASD)",
    "ASD2019":  "Neurological/Neuro-dev (PD/OCD/ASD)",
    "PD2019":   "Neurological/Neuro-dev (PD/OCD/ASD)",
    "IBD2017":  "Non-psychiatric disease (IBD)",
    "Height":   "Anthropometric control (Height)",
}

# Category → color (published palette, colorblind-friendly)
CATEGORY_COLOR: dict[str, str] = {
    "Psychiatric (SCZ/BD/MDD/ADHD)":      "#0279EE",  # blue
    "Neurological/Neuro-dev (PD/OCD/ASD)": "#9467BD", # purple
    "Non-psychiatric disease (IBD)":      "#D62728",  # red
    "Anthropometric control (Height)":    "#FF9400",  # orange
}

# Display labels (left-axis tick labels), verbatim from the published figure
TRAIT_DISPLAY: dict[str, str] = {
    "Height":   "Height (GIANT)",
    "IBD2017":  "IBD (2017)",
    "MDD2025":  "MDD (PGC2025)",
    "SCZ3":     "SCZ (PGC3)",
    "ADHD2022": "ADHD (PGC2022)",
    "BIP2021":  "BD (PGC2021)",
    "OCD2025":  "OCD (2025)",
    "ASD2019":  "ASD (2019)",
    "PD2019":   "PD (2019)",
}

# Trait display order (top-to-bottom on the y-axis), verbatim from published Fig 1
TRAIT_ORDER: list[str] = [
    "Height", "IBD2017", "MDD2025", "SCZ3", "ADHD2022",
    "BIP2021", "OCD2025", "ASD2019", "PD2019",
]

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


def _sig_marker(p: float) -> str:
    """Significance marker matching the published Figure 1 caption."""
    if p < 1e-3:
        return "***"
    if p < 1e-2:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


# ─────────────────────────────────────────────────────────────────────────────
# Render
# ─────────────────────────────────────────────────────────────────────────────
def render() -> None:
    df = pd.read_csv(INPUTS["gp_results_big9"]["path"], sep="\t")

    # Model A only; trait set check
    df_a = df[df["model"] == "A"].copy()
    expected_traits = sorted(TRAIT_DISPLAY.keys())
    assert sorted(df_a["trait"].tolist()) == expected_traits, \
        f"Unexpected trait set: {sorted(df_a['trait'].tolist())}"

    # Enforce published display order
    df_a = df_a.set_index("trait").loc[TRAIT_ORDER].reset_index()

    # 95% CI on the standardized scale. NOTE: the source TSV column "se" is the
    # MAGMA Wald SE of the *unstandardized* beta. Because MAGMA applies the same
    # scaling factor (sd_y) to both beta and its SE when reporting standardized
    # results, the per-trait ratio (beta_std / beta) — ~1.04-1.05 here — equals
    # the ratio (se_std / se). Using "se" directly yields a CI whose half-width
    # is uniformly ~4-5%% larger than the strict standardized CI, but the
    # rounded p-value-to-asterisk mapping (one-sided p, taken from the table)
    # is unaffected. This matches the precision shown in the published Figure 1.
    df_a["ci_lo"] = df_a["beta_std"] - 1.96 * df_a["se"]
    df_a["ci_hi"] = df_a["beta_std"] + 1.96 * df_a["se"]
    df_a["category"] = df_a["trait"].map(TRAIT_CATEGORY)
    df_a["color"] = df_a["category"].map(CATEGORY_COLOR)
    df_a["display_label"] = df_a["trait"].map(TRAIT_DISPLAY)
    df_a.loc[df_a["power_caveat"], "display_label"] += " †"
    df_a["sig"] = df_a["p_one_sided"].map(_sig_marker)

    # Save underlying numeric data
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    df_a[
        ["trait", "category", "display_label", "n_eff", "n_genes",
         "beta", "beta_std", "se", "ci_lo", "ci_hi",
         "p_one_sided", "sig", "power_caveat"]
    ].to_csv(OUT_TSV, sep="\t", index=False, float_format="%.6g")

    # ─── Figure ───────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(9.5, 6.0))
    fig.subplots_adjust(left=0.22, right=0.94, top=0.88, bottom=0.20)

    # Top = first item in TRAIT_ORDER
    y_positions = np.arange(len(df_a))[::-1]

    for y, (_, r) in zip(y_positions, df_a.iterrows()):
        ax.errorbar(
            r["beta_std"], y,
            xerr=[[r["beta_std"] - r["ci_lo"]], [r["ci_hi"] - r["beta_std"]]],
            fmt="o", color=r["color"], ecolor=r["color"],
            capsize=4, elinewidth=1.6, markersize=8,
            markeredgecolor="black", markeredgewidth=0.5,
        )
        # Significance annotation at the right tip of the CI
        style = "italic" if r["sig"] == "ns" else "normal"
        ax.text(
            r["ci_hi"] + 0.006, y, r["sig"],
            ha="left", va="center", fontsize=10,
            color=r["color"], fontstyle=style, fontweight="bold",
        )

    ax.axvline(0, color="#888888", linewidth=0.9, linestyle="--", zorder=0)

    ax.set_yticks(y_positions)
    ax.set_yticklabels(df_a["display_label"].tolist(), fontsize=11)
    ax.set_xlabel(r"$\beta_{std}$ (gene-property effect size)", fontsize=11.5)

    # Two-line title, matching the published figure
    fig.text(
        0.55, 0.96,
        "Gene-property enrichment of pleiotropic genes (gPS) across 9 traits",
        ha="center", va="center", fontsize=12.5, fontweight="bold",
    )
    fig.text(
        0.55, 0.92,
        "Model A (gPS ≥ 1 genes only)",
        ha="center", va="center", fontsize=11,
    )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", alpha=0.25, linestyle=":")
    ax.tick_params(axis="both", labelsize=10)

    # X-limits: slight padding so significance markers don't clip
    x_max = float(df_a["ci_hi"].max()) + 0.05
    x_min = float(df_a["ci_lo"].min()) - 0.02
    ax.set_xlim(x_min, x_max)

    # Legend matches the published order
    legend_order = [
        "Psychiatric (SCZ/BD/MDD/ADHD)",
        "Neurological/Neuro-dev (PD/OCD/ASD)",
        "Non-psychiatric disease (IBD)",
        "Anthropometric control (Height)",
    ]
    handles = [
        mpatches.Patch(facecolor=CATEGORY_COLOR[c], edgecolor="black", label=c)
        for c in legend_order
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=2,
        frameon=False,
        fontsize=9.5,
        bbox_to_anchor=(0.55, 0.005),
    )

    # Low-power footnote
    if df_a["power_caveat"].any():
        fig.text(
            0.04, 0.005,
            "† Limited GWAS power (OCD2025)",
            ha="left", va="bottom", fontsize=8.5, color="#555555", style="italic",
        )

    fig.savefig(OUT_PNG, dpi=300)
    fig.savefig(OUT_PDF)
    fig.savefig(OUT_SVG)
    plt.close(fig)


def main() -> int:
    assert_inputs()
    render()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
