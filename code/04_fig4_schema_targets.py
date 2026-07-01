#!/usr/bin/env python3
"""
04_fig4_schema_targets.py
=========================

Render Figure 4 — SCHEMA rare-variant target stratification across the
gPS × MAGMA-Z(SCZ) landscape, with bubble area proportional to SCHEMA PTV
odds ratio.

Plotting set
------------
The source table `schema_pav_gps_table.tsv` carries the 34 SCHEMA
FDR < 0.1 PTV genes. Per the user-confirmed brief, we exclude genes
absent from MAGMA (i.e. `scz_source == 'absent'`); this drops 5 genes
(GRIA3, SLF2, H1-4, MAGEC1, EIF2S3 — documented in
`inputs/SCHEMA_absent_from_MAGMA.md`) and leaves 29 plottable genes.

Aesthetic conventions
---------------------
* X-axis: Genetic Pleiotropy Score (gPS). Genes with `tsepilov_zone ==
  'not_in_gPS'` are anchored in a deterministic, low-discrepancy horizontal
  spread between 0.30 and 0.95 (left of the gPS=1 tick), matching the
  published display.
* Y-axis: MAGMA Z-score (SCZ).
* Bubble area ∝ SCHEMA OR (PTV); s = OR_PTV * 30 (legend shows OR=2/10/20).
  OR_PTV == inf gets a 99th-percentile fallback size so it stays on-canvas.
* Pleiotropy zone colors (matching the published legend):
    Optimal (gPS 2-5)        : #75A025 (green)
    High pleiotropy (gPS ≥ 6) : #FF9400 (orange)
    Low pleiotropy (gPS = 1)  : #0279EE (blue)
    Not in gPS table          : #BBBBBB (grey)
* SETD1A: bolded label with a leading ★, marking it as the single PAV-
  fraction-high coding gene (PAV-fraction = 0.83) in the published set.
* The "Optimal pleiotropy zone" band (x = 1.5..5.5) is shaded light green.

Inputs
------
* inputs/schema_pav_gps_table.tsv — MD5 b9c5f67a83510d0d87f94981914f1a40.

Outputs
-------
* results/figure4_schema_targets.{pdf,svg,png}
* results/figure4_schema_targets_data.tsv  — the 29 plottable genes
  with the computed display columns (including jittered x and zone).
"""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = ["Liberation Sans", "Arimo", "DejaVu Sans"]
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["svg.hashsalt"] = "fig4_schema_2026"

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import numpy as np
import pandas as pd

try:
    from adjustText import adjust_text
    HAS_ADJUST = True
except ImportError:
    HAS_ADJUST = False

# ─────────────────────────────────────────────────────────────────────────────
# Paths
# ─────────────────────────────────────────────────────────────────────────────
HERE = Path(__file__).resolve().parent
INPUTS_DIR = HERE.parent / "inputs"
RESULTS_DIR = HERE.parent / "results"

INPUTS = {
    "schema_pav_gps_table": {
        "path": INPUTS_DIR / "schema_pav_gps_table.tsv",
        "md5": "b9c5f67a83510d0d87f94981914f1a40",
    },
}

OUT_PDF = RESULTS_DIR / "figure4_schema_targets.pdf"
OUT_SVG = RESULTS_DIR / "figure4_schema_targets.svg"
OUT_PNG = RESULTS_DIR / "figure4_schema_targets.png"
OUT_TSV = RESULTS_DIR / "figure4_schema_targets_data.tsv"

# ─────────────────────────────────────────────────────────────────────────────
# Aesthetic configuration
# ─────────────────────────────────────────────────────────────────────────────
ZONE_COLOR: dict[str, str] = {
    "optimal_2to5":    "#75A025",   # green
    "high_gPS6plus":   "#FF9400",   # orange
    "low_gPS1":        "#0279EE",   # blue
    "not_in_gPS":      "#BBBBBB",   # grey
}
ZONE_LABEL: dict[str, str] = {
    "optimal_2to5":    "Optimal (gPS 2–5)",
    "high_gPS6plus":   "High pleiotropy (gPS ≥ 6)",
    "low_gPS1":        "Low pleiotropy (gPS = 1)",
    "not_in_gPS":      "Not in gPS table",
}
ZONE_ORDER: list[str] = ["optimal_2to5", "high_gPS6plus", "low_gPS1", "not_in_gPS"]

OR_REF_SIZES: list[int] = [2, 10, 20]   # legend reference odds-ratios
BUBBLE_SCALE = 30.0                      # s = OR * BUBBLE_SCALE
INF_OR_FALLBACK = 25.0                   # used in lieu of inf for bubble sizing

# Deterministic horizontal spread for `not_in_gPS` genes (avoids vertical stack)
NOT_IN_GPS_X_LO = 0.30
NOT_IN_GPS_X_HI = 0.95

# Optimal-zone band (light green, behind everything)
OPTIMAL_BAND_X_LO = 1.5
OPTIMAL_BAND_X_HI = 5.5

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


def _deterministic_spread(n: int, lo: float, hi: float) -> np.ndarray:
    """Deterministic, low-discrepancy spread of n points in [lo, hi]."""
    return np.linspace(lo, hi, n)


# ─────────────────────────────────────────────────────────────────────────────
# Render
# ─────────────────────────────────────────────────────────────────────────────
def render() -> None:
    df = pd.read_csv(INPUTS["schema_pav_gps_table"]["path"], sep="\t")

    # ── Filter: SCHEMA FDR<0.1 PTV genes, present in MAGMA ───────────────
    n_all = len(df)
    df = df[df["schema_category"] == "PTV"].copy()
    df = df[df["schema_fdr"] < 0.1].copy()
    df = df[df["scz_source"] != "absent"].copy()
    n_plot = len(df)
    if n_all != 34:
        sys.exit(f"ERROR: expected 34 SCHEMA rows, got {n_all}")
    if n_plot != 29:
        sys.exit(f"ERROR: expected 29 plottable rows, got {n_plot}")

    # Sort stably by gene_symbol so the not_in_gPS spread is deterministic
    df = df.sort_values("gene_symbol").reset_index(drop=True)

    # ── X coordinate: gPS for in_gPS, jittered spread for not_in_gPS ────
    not_in = df["tsepilov_zone"] == "not_in_gPS"
    xs = df["gPS"].astype(float).to_numpy(copy=True)
    spread = _deterministic_spread(int(not_in.sum()), NOT_IN_GPS_X_LO, NOT_IN_GPS_X_HI)
    xs[not_in.values] = spread
    df["x"] = xs
    df["y"] = df["magma_z_scz"].astype(float)

    # Bubble area
    or_for_size = df["OR_ptv"].replace([np.inf, -np.inf], np.nan).fillna(INF_OR_FALLBACK)
    df["bubble_s"] = or_for_size.clip(upper=INF_OR_FALLBACK) * BUBBLE_SCALE

    # Zone color
    df["color"] = df["tsepilov_zone"].map(ZONE_COLOR)
    assert df["color"].notna().all(), \
        f"Unmapped pleiotropy zone(s): {df.loc[df['color'].isna(), 'tsepilov_zone'].unique()}"

    # PAV-fraction lead marker.
    # The upstream source table does not carry a numeric PAV-fraction; per the
    # published manuscript §3.4 the only PAV-fraction > 0.5 lead in the set is
    # SETD1A (PAV-fraction = 0.83). All other genes are PAV-fraction 0-0.5.
    df["pav_lead"] = df["gene_symbol"].eq("SETD1A")

    # Save underlying numeric data
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    out_cols = [
        "gene_symbol", "ensembl_id", "gPS", "tsepilov_zone",
        "magma_z_scz", "magma_z_bip", "OR_ptv", "schema_p", "schema_fdr",
        "has_PAV", "pav_lead", "x", "y", "bubble_s", "color",
    ]
    df[out_cols].to_csv(OUT_TSV, sep="\t", index=False, float_format="%.6g")

    # ── Figure ───────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(11.5, 7.0))
    fig.subplots_adjust(left=0.07, right=0.78, top=0.95, bottom=0.10)

    # Fixed axis limits BEFORE plotting so the optimal-band label has a stable
    # vertical anchor.
    x_min, x_max = -0.5, 25.0
    y_min, y_max = -2.5, 7.0
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)

    # Optimal-zone shaded band (drawn first; sits behind everything else)
    ax.axvspan(
        OPTIMAL_BAND_X_LO, OPTIMAL_BAND_X_HI,
        color="#75A025", alpha=0.12, zorder=0,
    )

    # Optimal-zone label at the top of the band, inside axes
    ax.text(
        (OPTIMAL_BAND_X_LO + OPTIMAL_BAND_X_HI) / 2.0,
        y_max - 0.35,
        "Optimal\npleiotropy zone",
        ha="center", va="top", fontsize=9, color="#557018",
        style="italic", zorder=1,
    )

    # Zero reference line
    ax.axhline(0, color="#888888", linewidth=0.8, linestyle="--", zorder=1)

    # Bubbles, one zone at a time so the legend works cleanly
    for zone in ZONE_ORDER:
        sub = df[df["tsepilov_zone"] == zone]
        if sub.empty:
            continue
        ax.scatter(
            sub["x"], sub["y"],
            s=sub["bubble_s"], c=ZONE_COLOR[zone],
            alpha=0.75, edgecolor="black", linewidth=0.7,
            zorder=2,
        )

    # SETD1A PAV-lead star marker overlay
    pav_lead_rows = df[df["pav_lead"]]
    for _, r in pav_lead_rows.iterrows():
        ax.scatter(
            r["x"], r["y"], marker="*",
            s=160, c="black", zorder=4,
        )

    # Gene labels (slight upward offset, adjustText handles collisions)
    texts = []
    for _, r in df.iterrows():
        label = r["gene_symbol"]
        if r["pav_lead"]:
            label = f"{label} ★"
        weight = "bold" if r["pav_lead"] else "normal"
        texts.append(
            ax.text(
                r["x"], r["y"], label,
                fontsize=8.5, ha="center", va="center",
                fontweight=weight, zorder=5,
            )
        )

    # NOTE: adjust_text iterates until label-collision convergence and is
    # not bit-deterministic across runs (slight pixel/layout drift between
    # invocations). The underlying numeric data (figure4_schema_targets_data.tsv)
    # IS bit-stable; the rendered figures are visually-faithful but not
    # MD5-stable. This matches the brief's "numeric/vector parity, NOT
    # pixel parity" requirement.
    if HAS_ADJUST:
        adjust_text(
            texts,
            ax=ax,
            expand=(1.15, 1.25),
            force_text=(0.4, 0.7),
            force_static=(0.05, 0.1),
            only_move={"texts": "xy"},
            arrowprops=dict(arrowstyle="-", color="grey", lw=0.45, alpha=0.5),
        )

    # Axes
    ax.set_xlabel("Genetic Pleiotropy Score (gPS)", fontsize=11)
    ax.set_ylabel("MAGMA Z-score (SCZ)", fontsize=11)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", labelsize=10)

    # ── Legend: pleiotropy zones (top right) ─────────────────────────────
    zone_handles = [
        mpatches.Patch(facecolor=ZONE_COLOR[z], edgecolor="black", label=ZONE_LABEL[z])
        for z in ZONE_ORDER
    ]
    leg1 = ax.legend(
        handles=zone_handles,
        title="Pleiotropy zone",
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        fontsize=9, title_fontsize=9.5, frameon=False,
    )
    ax.add_artist(leg1)

    # ── Legend: SCHEMA OR bubble sizes (middle right) ────────────────────
    size_handles = []
    for OR in OR_REF_SIZES:
        size_handles.append(
            mlines.Line2D(
                [], [], color="white", marker="o", linestyle="None",
                markerfacecolor="#DDDDDD", markeredgecolor="black",
                markersize=np.sqrt(OR * BUBBLE_SCALE),
                label=f"OR = {OR}",
            )
        )
    leg2 = ax.legend(
        handles=size_handles,
        title="SCHEMA OR (PTV)",
        loc="center left",
        bbox_to_anchor=(1.02, 0.55),
        fontsize=9, title_fontsize=9.5, frameon=False,
        labelspacing=1.6,
    )
    ax.add_artist(leg2)

    # ── Legend: PAV-fraction marker (bottom right) ───────────────────────
    pav_handles = [
        mlines.Line2D(
            [], [], color="black", marker="*", linestyle="None",
            markersize=12, label="PAV-fraction ≥ 0.5\n(fine-mapped PAV lead)",
        ),
        mlines.Line2D(
            [], [], color="white", marker="o", linestyle="None",
            markerfacecolor="#DDDDDD", markeredgecolor="black",
            markersize=8, label="PAV-fraction 0–0.5",
        ),
    ]
    ax.legend(
        handles=pav_handles,
        title="PAV-fraction (VEP, PP≥0.5)",
        loc="lower left",
        bbox_to_anchor=(1.02, 0.05),
        fontsize=9, title_fontsize=9.5, frameon=False,
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
