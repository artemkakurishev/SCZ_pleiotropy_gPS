#!/usr/bin/env python3
"""
05_fig5_loo_sensitivity.py
==========================

Render Figure 5 — the gPS leave-one-out (LOO) sensitivity line plot for
SCZ3 vs BD2021 across four progressively broader disease-removal tiers.

Plot layout
-----------
Two lines (SCZ3 blue, BD 2021 orange), each carrying four points at
x = Baseline / Core psych (−19) / All psych (−61) / Psych + CNS (−216).
Y = β_std with 95% CI errorbars. Markers are filled when p_onesided < 0.05
and open (white-filled, colored edge) when p_onesided ≥ 0.05.

Bottom annotations show gPS variance retained vs baseline (100%, 94%,
91%, 79%), computed as 100 - gps_variance_reduction_pct from the source
table.

Inputs
------
* inputs/loo_gps_sensitivity_frontiers.tsv — MD5
  16b15719d497627f8cc8dee92f4ddedf. Columns: trait, scenario, n_genes,
  beta, beta_std, SE, p_onesided, attenuation_pct, gps_variance_reduction_pct.

Outputs
-------
* results/figure5_loo_sensitivity.{pdf,svg,png}
* results/figure5_loo_sensitivity_data.tsv — the rendered 8 rows
  (2 traits × 4 scenarios) with computed display columns.
"""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = ["Liberation Sans", "Arimo", "DejaVu Sans"]
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["svg.hashsalt"] = "fig5_loo_2026"

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import pandas as pd

# ─────────────────────────────────────────────────────────────────────────────
# Paths
# ─────────────────────────────────────────────────────────────────────────────
HERE = Path(__file__).resolve().parent
INPUTS_DIR = HERE.parent / "inputs"
RESULTS_DIR = HERE.parent / "results"

INPUTS = {
    "loo_gps_sensitivity": {
        "path": INPUTS_DIR / "loo_gps_sensitivity_frontiers.tsv",
        "md5": "16b15719d497627f8cc8dee92f4ddedf",
    },
}

OUT_PDF = RESULTS_DIR / "figure5_loo_sensitivity.pdf"
OUT_SVG = RESULTS_DIR / "figure5_loo_sensitivity.svg"
OUT_PNG = RESULTS_DIR / "figure5_loo_sensitivity.png"
OUT_TSV = RESULTS_DIR / "figure5_loo_sensitivity_data.tsv"

# ─────────────────────────────────────────────────────────────────────────────
# Display configuration — matches published Figure 5 verbatim
# ─────────────────────────────────────────────────────────────────────────────

# Scenario ordering left-to-right on the x-axis
SCENARIO_ORDER: list[str] = ["original", "noPsych_core", "noPsych_OT", "noPsychNoNerv"]

# Bottom-row tick labels (two-line, with the disease-count delta in parentheses)
SCENARIO_TICK_LABELS: dict[str, str] = {
    "original":      "Baseline\n(full gPS)",
    "noPsych_core":  "Core psych\n(−19)",
    "noPsych_OT":    "All psych\n(−61)",
    "noPsychNoNerv": "Psych + CNS\n(−216)",
}

# Trait → display label / color
TRAIT_LABEL: dict[str, str] = {
    "SCZ3":    "SCZ3",
    "BIP2021": "BD 2021",
}
TRAIT_COLOR: dict[str, str] = {
    "SCZ3":    "#0279EE",    # blue
    "BIP2021": "#FF9400",    # orange
}

# Significance threshold for filled vs open markers
P_THRESHOLD = 0.05

# Horizontal jitter for the two traits at each x so their CIs don't overplot
JITTER = 0.05

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
# Render
# ─────────────────────────────────────────────────────────────────────────────
def render() -> None:
    df = pd.read_csv(INPUTS["loo_gps_sensitivity"]["path"], sep="\t")

    # Schema sanity
    expected_cols = {
        "trait", "scenario", "n_genes", "beta", "beta_std", "SE",
        "p_onesided", "attenuation_pct", "gps_variance_reduction_pct",
    }
    missing = expected_cols - set(df.columns)
    if missing:
        sys.exit(f"ERROR: missing columns in LOO table: {missing}")

    expected_traits = {"SCZ3", "BIP2021"}
    if set(df["trait"].unique()) != expected_traits:
        sys.exit(f"ERROR: unexpected trait set: {sorted(df['trait'].unique())}")
    if set(df["scenario"].unique()) != set(SCENARIO_ORDER):
        sys.exit(f"ERROR: unexpected scenario set: {sorted(df['scenario'].unique())}")

    # Derived columns
    df["ci_lo"] = df["beta_std"] - 1.96 * df["SE"]
    df["ci_hi"] = df["beta_std"] + 1.96 * df["SE"]
    df["filled"] = df["p_onesided"] < P_THRESHOLD
    df["variance_retained_pct"] = 100.0 - df["gps_variance_reduction_pct"]

    # Save numeric outputs (sorted for stability)
    df_out = df.sort_values(["trait", "scenario"]).copy()
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    df_out[
        ["trait", "scenario", "n_genes", "beta", "beta_std", "SE",
         "ci_lo", "ci_hi", "p_onesided", "filled",
         "attenuation_pct", "gps_variance_reduction_pct",
         "variance_retained_pct"]
    ].to_csv(OUT_TSV, sep="\t", index=False, float_format="%.6g")

    # ── Figure ───────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(9.0, 5.5))
    fig.subplots_adjust(left=0.10, right=0.95, top=0.91, bottom=0.20)

    x_pos = np.arange(len(SCENARIO_ORDER), dtype=float)
    x_for_trait: dict[str, np.ndarray] = {
        "SCZ3":    x_pos - JITTER,
        "BIP2021": x_pos + JITTER,
    }

    for trait in ["SCZ3", "BIP2021"]:
        color = TRAIT_COLOR[trait]
        # Order rows by SCENARIO_ORDER
        sub = (df.set_index(["trait", "scenario"])
                 .loc[trait]
                 .reindex(SCENARIO_ORDER)
                 .reset_index())
        xs = x_for_trait[trait]
        ys = sub["beta_std"].to_numpy()
        yerr = np.array([sub["beta_std"] - sub["ci_lo"],
                         sub["ci_hi"] - sub["beta_std"]])

        # Line connecting the points
        ax.plot(xs, ys, color=color, linewidth=1.8, zorder=2)

        # Errorbar caps + markers (filled vs open by significance)
        for x, y, lo, hi, filled in zip(
            xs,
            ys,
            yerr[0],
            yerr[1],
            sub["filled"].to_numpy(),
        ):
            ax.errorbar(
                x, y, yerr=[[lo], [hi]],
                ecolor=color, elinewidth=1.4, capsize=4, capthick=1.2,
                fmt="none", zorder=2,
            )
            if filled:
                ax.scatter(x, y, s=70, c=color, edgecolor=color,
                           linewidth=1.4, zorder=4)
            else:
                ax.scatter(x, y, s=70, facecolors="white", edgecolor=color,
                           linewidth=1.6, zorder=4)

    # Zero reference line
    ax.axhline(0, color="#888888", linewidth=0.8, linestyle=":", zorder=1)

    # X-axis ticks
    ax.set_xticks(x_pos)
    ax.set_xticklabels(
        [SCENARIO_TICK_LABELS[s] for s in SCENARIO_ORDER],
        fontsize=10,
    )
    ax.set_xlim(-0.5, len(SCENARIO_ORDER) - 0.5)

    # Y-axis
    ax.set_ylabel(r"Gene-property coefficient $\beta_{std}$ (95% CI)",
                  fontsize=11)

    # Top title
    ax.set_title(
        "gPS leave-one-out sensitivity across removal tiers",
        fontsize=12.5, fontweight="bold", pad=10,
    )

    # Variance-retained annotations along the bottom (above the x-tick labels)
    # Use Axes coordinates for x, data coordinates for y.
    var_ret = (
        df_out.set_index(["trait", "scenario"])
              .xs("SCZ3")
              .reindex(SCENARIO_ORDER)["variance_retained_pct"]
              .to_numpy()
    )
    # SCZ and BD share the same gPS variance reduction (the gPS table is
    # recomputed once per scenario, then both traits regressed on it). Sanity
    # check.
    var_ret_bip = (
        df_out.set_index(["trait", "scenario"])
              .xs("BIP2021")
              .reindex(SCENARIO_ORDER)["variance_retained_pct"]
              .to_numpy()
    )
    if not np.allclose(var_ret, var_ret_bip, atol=1e-6):
        sys.exit("ERROR: gps_variance_reduction_pct differs between SCZ3 and BIP2021; "
                 "the column is expected to be a property of the gPS recomputation, not trait.")

    # Place pct labels just above the tick-label band
    y_anchor = ax.get_ylim()[0] + 0.005
    for x, pct in zip(x_pos, var_ret):
        ax.text(
            x, y_anchor, f"{pct:.0f}%",
            ha="center", va="bottom", fontsize=10.5, fontweight="bold",
            color="#444444",
        )

    # Italic caption for the variance-retained row (placed under the x-ticks)
    fig.text(
        0.525, 0.10,
        "gPS variance retained (% of baseline)",
        ha="center", va="center", fontsize=9, style="italic", color="#555555",
    )

    # Sub-x-axis label (below the tick labels)
    fig.text(
        0.525, 0.025,
        "LOO removal tier (gPS recomputed after disease removal)",
        ha="center", va="center", fontsize=10,
    )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", labelsize=10)

    # ── Legend ───────────────────────────────────────────────────────────
    legend_handles = [
        mlines.Line2D([], [], color=TRAIT_COLOR["SCZ3"], marker="o",
                      markersize=8, linewidth=2, label=TRAIT_LABEL["SCZ3"]),
        mlines.Line2D([], [], color=TRAIT_COLOR["BIP2021"], marker="o",
                      markersize=8, linewidth=2, label=TRAIT_LABEL["BIP2021"]),
        mlines.Line2D([], [], color="black", marker="o", linestyle="None",
                      markersize=8, label="p < 0.05"),
        mlines.Line2D([], [], color="black", marker="o", linestyle="None",
                      markerfacecolor="white", markeredgewidth=1.5,
                      markersize=8, label="n.s. (p ≥ 0.05)"),
    ]
    ax.legend(
        handles=legend_handles,
        loc="upper right",
        fontsize=9.5, frameon=False,
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
