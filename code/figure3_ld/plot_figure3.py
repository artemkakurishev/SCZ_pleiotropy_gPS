#!/usr/bin/env python
"""
plot_figure3.py — PONE-D-26-36293 revision: combined two-panel Figure 3.

  Panel A: the five-disease inferential plot under the PRIMARY LD-block,
           standard-gPS reference model (delta scale in percentage points;
           solid zero reference line; q_BH5 column).
  Panel B: T2D across three LD-block reference models on an ABSOLUTE percentage
           scale (block_size_only / block_standard_gps / block_exact_gps_minus_d).
           Grey bands = central 95% of each model's empirical null distribution
           of T2D annotation percentages; internal black tick = null mean;
           dashed vertical line = observed T2D percentage (fixed across models).

PLOTTING ONLY. No resampling, no matching, no p/q computation. Every plotted value
is read from hash-verified tables; the only derived quantities are deterministic
arithmetic re-expressions of table cells:

  Panel A:
    null_mean_pct       = observed_pct - bw_resid
    [Q0.025, Q0.975](N) = [observed_pct - bw_hi, observed_pct - bw_lo]
    delta_pp            = observed_pct - null_mean_pct
    centred band        = [Q0.025(N) - null_mean, Q0.975(N) - null_mean]
  Panel B (empirically verified below against the null matrix for
  block_standard_gps):
    ref95_lo = observed_pct - null_quantile_97.5
    ref95_hi = observed_pct - null_quantile_2.5
    => absolute band   = [observed_pct - ref95_hi, observed_pct - ref95_lo]

Bands are empirical-null reference ranges, NOT confidence intervals; they are not
forced symmetric and are not centred on the observed values. q-values are the
existing Benjamini-Hochberg adjustments across the five primary disease axes,
computed separately within each specification (NOT recomputed across the three
T2D rows). The exact-exclusion row uses the disease-ledger exclusion result
(MONDO_0005148 removed from the gPS ledger), not the L2G>=0.5 approximation.

Two modes:
  FAST (default):     python code/figure3_ld/plot_figure3.py
      Reads the frozen tables in data/figure3_ld/frozen/ (SHA256-verified) and
      builds Fig3.tif in seconds.
  FROM RECOMPUTED:    python code/figure3_ld/plot_figure3.py --from-recomputed
      Reads results/figure3_ld/recomputed/ produced by
      `python code/figure3_ld/run_analysis.py --full` (which already compared
      those tables against the frozen ones and halted on discrepancy).

Font: Arial with Arimo/DejaVu Sans fallbacks (PLOS-permitted). Arial lacks
Unicode superscript glyphs, so scientific notation uses mathtext with an Arial
fontset. Where genuine Arial is not installed, Arimo (metric-compatible) is used
and the figure geometry is unchanged.

Export: .pdf/.svg + 600-dpi opaque PNG; Pillow -> RGB LZW TIFF (dpi tag 600),
staged on local disk then copied (S3 mounts cannot seek). TIFF, typography
(8-12 pt) and stroke (>= 0.75 pt) audits halt on violation.

Adapted from the verified revision figure package make_figure3_AB.py; the
assertion battery is preserved in full.
"""

from pathlib import Path
import argparse
import hashlib
import json
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.transforms import blended_transform_factory
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
FROZEN = REPO / "data" / "figure3_ld" / "frozen"
RECOMPUTED = REPO / "results" / "figure3_ld" / "recomputed"
OUTDIR = REPO / "results" / "figure3_ld"

# ----------------------------------------------------------------------------- constants
matplotlib.rcParams.update({
    # Arial (PLOS-permitted); Arimo/DejaVu Sans as declared fallbacks.
    "font.family": ["Arial", "Arimo", "DejaVu Sans"],
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold",
    "mathtext.cal": "Arial:italic",
    "mathtext.default": "regular",
    "font.size": 9,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 9,
    "legend.fontsize": 8,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "lines.linewidth": 1.0,
    "patch.linewidth": 0.75,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "savefig.transparent": False,
})

DPI = 600
UNIVERSE = 7707
RISK_N = 1016
B_REPL = 10000

AXES5 = ["T2D", "CAD", "Hypertension", "Asthma", "Prostate Ca"]
AXES8 = AXES5 + ["Atrial Fib", "Gout", "OA"]
DISPLAY = {
    "T2D": "Type 2 diabetes",
    "CAD": "Coronary artery disease",
    "Hypertension": "Hypertension",
    "Asthma": "Asthma",
    "Prostate Ca": "Prostate cancer",
    "Atrial Fib": "Atrial fibrillation",
    "Gout": "Gout",
    "OA": "Osteoarthritis",
}

C_OBS = "#000000"   # observed: prominent black
C_BAND = "#7A7A7A"  # empirical-null reference band (light grey fill at alpha 0.35)
C_ZERO = "#333333"  # Panel A zero reference line

# Panel B: three LD-block reference models, display order top -> bottom
MODELS = ["block_size_only", "block_standard_gps", "block_exact_gps_minus_d"]
MODEL_LABELS = {
    "block_size_only": "Without gPS",
    "block_standard_gps": "With gPS",
    "block_exact_gps_minus_d": "With gPS excluding T2D",
}
# validation-only expected absolute central-95% null ranges (approximate)
EXPECTED_BANDS = {"block_size_only": (8.54, 12.22),
                  "block_standard_gps": (10.87, 13.96),
                  "block_exact_gps_minus_d": (11.32, 14.44)}
EXPECTED_Q_DISP = {"block_size_only": "0.004", "block_standard_gps": "0.272",
                   "block_exact_gps_minus_d": "0.510"}

# Frozen spot values for Panel A (second-layer assertions)
FROZEN_Q5 = {"T2D": 0.272223, "CAD": 0.00049995, "Hypertension": 0.29917,
             "Asthma": 0.0023331, "Prostate Ca": 0.00049995}
FROZEN_P = {"T2D": 0.217778, "CAD": 0.00019998, "Hypertension": 0.29917,
            "Asthma": 0.00139986, "Prostate Ca": 0.00019998}


# ----------------------------------------------------------------------------- source gate
def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def gate_frozen() -> None:
    """Fast mode: verify the frozen tables against their SHA256 manifest."""
    sums = {}
    for ln in (FROZEN / "SHA256SUMS.txt").read_text().splitlines():
        h, fn = ln.split()
        sums[fn] = h
    bad = [fn for fn, h in sums.items() if sha256(FROZEN / fn) != h]
    if bad:
        sys.exit(f"ERROR: frozen-table SHA256 mismatch: {bad}")
    print(f"Frozen-table SHA256 gate: PASS ({len(sums)} files)")


def gate_recomputed() -> None:
    """Recomputed mode: the re-run must have completed with a clean comparison."""
    rep = RECOMPUTED / "comparison_report.md"
    if not rep.exists():
        sys.exit("ERROR: results/figure3_ld/recomputed/ not found — run "
                 "`python code/figure3_ld/run_analysis.py --full` first.")
    if "**FAIL**" in rep.read_text():
        sys.exit("ERROR: recomputed tables failed the frozen comparison — see "
                 "results/figure3_ld/recomputed/comparison_report.md; refusing to "
                 "plot unverified values.")
    print("Recomputed-table comparison gate: PASS (comparison_report.md clean)")


# ----------------------------------------------------------------------------- formatting
def fmt_sci(v: float) -> str:
    """2-significant-figure scientific notation with a mathtext exponent (Arial
    lacks Unicode superscript glyphs)."""
    mant, exp = f"{v:.1e}".split("e")
    return f"{mant}×$10^{{{int(exp)}}}$"


def fmt_q(q: float) -> str:
    """Panel A q formatting: 3 dp when >= 0.01, else 2-sig-fig scientific.
    Never prints a small nonzero value as zero."""
    return f"{q:.3f}" if q >= 0.01 else fmt_sci(q)


# ----------------------------------------------------------------------------- loading
def load_panel_a(src: Path) -> pd.DataFrame:
    """Panel A values + the established assertion battery."""
    t8 = pd.read_csv(src / "Table_S_block_vs_gene_bootstrap.tsv", sep="\t")
    f3b = pd.read_csv(src / "Figure3B_main_source_data.tsv", sep="\t")
    null_mat = pd.read_csv(src / "null_axis_percentages_block_gps_matched.tsv", sep="\t")

    assert len(t8) == 8 and set(t8["axis"]) == set(AXES8), "8-axis table shape/content"
    assert len(f3b) == 5 and set(f3b["axis"]) == set(AXES5), "Figure3B source shape/content"
    assert null_mat.shape == (B_REPL, 9), f"null matrix shape {null_mat.shape}"

    rows = []
    for _, r in t8.iterrows():
        a = r["axis"]
        obs = float(r["obs_pct"])
        mu = obs - float(r["bw_resid"])
        q025 = obs - float(r["bw_hi"])
        q975 = obs - float(r["bw_lo"])
        rows.append(dict(
            axis=a, axis_display=DISPLAY[a], primary=bool(r["primary"]),
            observed_pct=obs, null_mean_pct=mu, null_q025=q025, null_q975=q975,
            band_lo_c=q025 - mu, band_hi_c=q975 - mu, delta_pp=obs - mu,
            p_two=float(r["bw_p"]),
            q_BH5=float(r["bw_q_BH5"]) if pd.notna(r["bw_q_BH5"]) else np.nan,
        ))
    df = pd.DataFrame(rows).set_index("axis").loc[AXES8].reset_index()
    d5 = df[df["primary"]].reset_index(drop=True)

    # 1. cross-consistency with the Figure3B source (5 axes).
    #    Tolerances are rounding-aware (8-axis table stores 4 dp).
    f3b_i = f3b.set_index("axis")
    for _, g in d5.iterrows():
        a, f = g["axis"], f3b_i.loc[g["axis"]]
        assert abs(g["observed_pct"] - f["observed_pct"]) < 6e-5, f"{a} obs mismatch"
        assert abs(g["null_mean_pct"] - f["null_mean_pct"]) < 2e-4, f"{a} null mean"
        assert abs(g["null_q025"] - f["null_band_lo_q0.025"]) < 2e-4, f"{a} q025"
        assert abs(g["null_q975"] - f["null_band_hi_q0.975"]) < 2e-4, f"{a} q975"
        assert abs(g["delta_pp"] - f["resid_pp"]) < 2e-4, f"{a} resid"
        assert abs(g["p_two"] - f["p_two"]) < 6e-5, f"{a} p"
        assert abs(g["q_BH5"] - f["q_BH5"]) < 6e-5, f"{a} q_BH5"

    # 2. direct verification against the 10,000-replicate null matrix
    for _, g in d5.iterrows():
        col = null_mat[g["axis"]].to_numpy()
        assert abs(col.mean() - g["null_mean_pct"]) < 1e-3, f'{g["axis"]} matrix mean'
        assert abs(np.quantile(col, 0.025) - g["null_q025"]) < 1e-3, f'{g["axis"]} matrix q025'
        assert abs(np.quantile(col, 0.975) - g["null_q975"]) < 1e-3, f'{g["axis"]} matrix q975'
        assert g["null_q025"] < g["null_mean_pct"] < g["null_q975"], f'{g["axis"]} band order'

    # 3. frozen spot values (second layer)
    for _, g in d5.iterrows():
        a = g["axis"]
        assert abs(g["p_two"] - FROZEN_P[a]) < 6e-5, f"{a} p frozen"
        assert abs(g["q_BH5"] - FROZEN_Q5[a]) < 6e-5, f"{a} q_BH5 frozen"
        assert abs(g["delta_pp"] - (g["observed_pct"] - g["null_mean_pct"])) < 1e-12
        assert abs(g["band_lo_c"] - (g["null_q025"] - g["null_mean_pct"])) < 1e-12
        assert abs(g["band_hi_c"] - (g["null_q975"] - g["null_mean_pct"])) < 1e-12
    return d5


def load_panel_b(src: Path, null_mat) -> pd.DataFrame:
    """T2D rows for the three LD-block reference models + validation battery."""
    sens = pd.read_csv(src / "Table_final_five_axis_sensitivity.tsv", sep="\t")
    exact = pd.read_csv(src / "Table_final_five_axis_exact_gps_minus_d.tsv", sep="\t")

    pb = (sens[(sens.axis == "T2D") & (sens.specification.isin(MODELS))]
          .set_index("specification").loc[MODELS].reset_index())
    assert len(pb) == 3, f"expected exactly 3 T2D model rows, got {len(pb)}"
    assert pb.observed_pct.nunique() == 1, "observed % must be identical across models"
    assert (pb.B == B_REPL).all(), "each model must use 10,000 replicates"

    # exact disease-ledger exclusion, not the L2G>=0.5 approximation:
    # residual must be +0.523 pp (approximate-exclusion result was ~+0.59 pp)
    r_ex = pb.loc[pb.specification == "block_exact_gps_minus_d"].iloc[0]
    assert abs(r_ex.resid_pp - 0.522595) < 1e-4, f"exact-exclusion resid {r_ex.resid_pp}"

    # cross-check against the dedicated exact-exclusion table
    ex = exact[exact.axis == "T2D"].iloc[0]
    assert ex.disease_ids_removed == "MONDO_0005148", "T2D ledger removal (MONDO_0005148)"
    for col in ("observed_pct", "null_mean", "resid_pp", "ref95_lo", "ref95_hi", "q_BH5"):
        assert abs(ex[col] - r_ex[col]) < 1e-6, f"exact-table cross-check failed: {col}"

    # documented absolute-band transformation:
    # ref95_lo = obs - q97.5 ; ref95_hi = obs - q2.5  =>  band = [obs-ref95_hi, obs-ref95_lo]
    pb["null_lo"] = pb.observed_pct - pb.ref95_hi
    pb["null_hi"] = pb.observed_pct - pb.ref95_lo
    for m, (lo, hi) in EXPECTED_BANDS.items():
        r = pb.loc[pb.specification == m].iloc[0]
        assert abs(r.null_lo - lo) < 0.01 and abs(r.null_hi - hi) < 0.01, \
            f"{m} band [{r.null_lo:.4f}, {r.null_hi:.4f}] vs expected ~[{lo}, {hi}]"
        assert r.null_lo < r.null_mean < r.null_hi, f"{m} null mean outside band"

    # empirical verification of the transformation for block_standard_gps against
    # the 10,000-replicate null matrix (absolute T2D percentages)
    r_std = pb.loc[pb.specification == "block_standard_gps"].iloc[0]
    col = null_mat["T2D"].to_numpy()
    assert abs(np.quantile(col, 0.025) - r_std.null_lo) < 2e-3, "matrix q025 vs table"
    assert abs(np.quantile(col, 0.975) - r_std.null_hi) < 2e-3, "matrix q975 vs table"
    assert abs(col.mean() - r_std.null_mean) < 2e-3, "matrix mean vs table"

    # displayed q-values: existing BH across the five primary axes within each model
    for m, qs in EXPECTED_Q_DISP.items():
        q = float(pb.loc[pb.specification == m, "q_BH5"].iloc[0])
        assert f"{q:.3f}" == qs, f"{m} q_BH5 {q} displays as {q:.3f}, expected {qs}"

    # geometric expectations for the observed line
    obs = float(pb.observed_pct.iloc[0])
    r_size = pb.loc[pb.specification == "block_size_only"].iloc[0]
    assert obs > r_size.null_hi, "observed must fall right of the size-only band"
    for m in ("block_standard_gps", "block_exact_gps_minus_d"):
        r = pb.loc[pb.specification == m].iloc[0]
        assert r.null_lo < obs < r.null_hi, f"observed must fall inside the {m} band"

    pb["display_label"] = pb.specification.map(MODEL_LABELS)
    return pb


# ----------------------------------------------------------------------------- drawing
def _rows(n):
    return np.arange(n)[::-1]  # first row on top


def _style_ax(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", width=0.8)


def draw_panel_a(ax, axq, d5):
    """Panel A grammar: black point = observed - matched-null mean; grey band =
    central 95% of the centred empirical null (not forced symmetric, not a CI);
    thin zero line; single right column 'FDR q' (q_BH5)."""
    n = len(d5)
    y = _rows(n)
    for yi, (_, r) in zip(y, d5.iterrows()):
        ax.barh(yi, r["band_hi_c"] - r["band_lo_c"], left=r["band_lo_c"], height=0.52,
                facecolor=C_BAND, alpha=0.35, edgecolor="#4D4D4D", linewidth=0.75,
                zorder=2)
    ax.axvline(0, color=C_ZERO, linewidth=0.9, zorder=3)  # above bands, below points
    for yi, (_, r) in zip(y, d5.iterrows()):
        ax.plot(r["delta_pp"], yi, marker="o", markersize=7.0,
                markerfacecolor=C_OBS, markeredgecolor=C_OBS, markeredgewidth=0.9,
                zorder=4, linestyle="none")
    ax.set_yticks(y)
    ax.set_yticklabels(d5["axis_display"], fontsize=11)
    ax.set_xlim(-4.35, 2.35)
    ax.set_xticks([-4, -2, 0, 2])
    ax.tick_params(axis="x", labelsize=10)
    ax.set_xlabel("Overlap difference (percentage points)", fontsize=11)
    ax.set_ylim(-0.62, n - 1 + 1.00)
    _style_ax(ax)
    ytop = n - 1 + 0.70
    ax.text(-4.20, ytop, "Less than expected", ha="left", va="center",
            fontsize=10, color="#4D4D4D")
    ax.text(2.20, ytop, "More than expected", ha="right", va="center",
            fontsize=10, color="#4D4D4D")

    axq.set_xlim(0, 1)
    axq.set_ylim(-0.62, n - 1 + 1.00)
    axq.axis("off")
    axq.text(0.5, ytop, "FDR q", ha="center", va="center", fontsize=11,
             fontweight="bold")
    for yi, (_, r) in zip(y, d5.iterrows()):
        axq.text(0.5, yi, fmt_q(r["q_BH5"]), ha="center", va="center", fontsize=10.5)
    return ax, axq


def draw_panel_b(ax, axq, pb, title=True):
    """Panel B grammar: T2D under three LD-block reference models, ABSOLUTE scale.

    Grey band = central 95% of the model's empirical null distribution of T2D
    annotation percentages; short black tick = null mean; dashed vertical line =
    observed T2D percentage (fixed), labelled once above the plotting region.
    """
    n = len(pb)
    y = _rows(n)
    obs = float(pb.observed_pct.iloc[0])
    for yi, (_, r) in zip(y, pb.iterrows()):
        ax.barh(yi, r["null_hi"] - r["null_lo"], left=r["null_lo"], height=0.52,
                facecolor=C_BAND, alpha=0.35, edgecolor="#4D4D4D", linewidth=0.75,
                zorder=2)
        # null mean: short solid vertical black tick spanning the band height
        ax.plot([r["null_mean"], r["null_mean"]], [yi - 0.26, yi + 0.26],
                color=C_OBS, linewidth=1.4, zorder=3, solid_capstyle="butt")
    # observed: one dashed vertical black line across all rows
    ax.axvline(obs, color=C_OBS, linewidth=1.0, linestyle=(0, (4, 2)), zorder=4)
    ax.set_yticks(y)
    ax.set_yticklabels(pb["display_label"], fontsize=11)
    ax.set_xlim(7.9, 15.0)
    ax.set_xticks([8, 10, 12, 14])
    ax.tick_params(axis="x", labelsize=10)
    ax.set_xlabel("T2D-annotated genes (%)", fontsize=11)
    ax.set_ylim(-0.62, n - 1 + 1.00)
    _style_ax(ax)

    # panel title (left) and single observed-line label (at the line's x),
    # both above the plotting region on the same baseline -> no collision
    trans = blended_transform_factory(ax.transData, ax.transAxes)
    if title:
        ax.text(0.0, 1.045, "T2D across reference models", transform=ax.transAxes,
                ha="left", va="bottom", fontsize=11, fontweight="bold")
    ax.text(obs, 1.045, f"Observed: {obs:.2f}%", transform=trans,
            ha="center", va="bottom", fontsize=10)

    # single aligned numerical column: FDR q (q_BH5 within each specification)
    ytop = n - 1 + 0.70
    axq.set_xlim(0, 1)
    axq.set_ylim(-0.62, n - 1 + 1.00)
    axq.axis("off")
    axq.text(0.5, ytop, "FDR q", ha="center", va="center", fontsize=11,
             fontweight="bold")
    for yi, (_, r) in zip(y, pb.iterrows()):
        axq.text(0.5, yi, f"{r['q_BH5']:.3f}", ha="center", va="center", fontsize=10.5)
    return ax, axq


# ----------------------------------------------------------------------------- audits
def audit_figure(fig, name):
    """All rendered text 8-12 pt; all drawn strokes >= 0.75 pt. Halt on violation."""
    sizes, widths = [], []
    for t in fig.findobj(matplotlib.text.Text):
        if t.get_text().strip():
            sizes.append(float(t.get_size()))
    for ln in fig.findobj(Line2D):
        lw = float(ln.get_linewidth())
        if lw > 0 and (ln.get_linestyle() not in ("None", "none", "")):
            widths.append(lw)
        mew = float(ln.get_markeredgewidth())
        if mew > 0 and ln.get_marker() not in ("None", "none", "", None):
            widths.append(mew)
    for p in fig.findobj(matplotlib.patches.Patch):
        lw = float(p.get_linewidth())
        if lw > 0:
            widths.append(lw)
    for ax in fig.axes:
        for sp in ax.spines.values():
            if sp.get_visible():
                widths.append(float(sp.get_linewidth()))
    fmin, fmax = min(sizes), max(sizes)
    wmin = min(widths)
    assert 8.0 <= fmin <= 12.0, f"{name}: font min {fmin} pt outside [8, 12]"
    assert 8.0 <= fmax <= 12.0, f"{name}: font max {fmax} pt outside [8, 12]"
    assert wmin >= 0.75, f"{name}: stroke min {wmin} pt < 0.75"
    return {"figure": name, "font_min_pt": fmin, "font_max_pt": fmax,
            "stroke_min_pt": round(wmin, 3)}


def audit_fonts(fig, name):
    """Every rendered text artist must resolve to Arial (or declared fallbacks)."""
    from matplotlib import font_manager as fm
    allowed = {"Arial", "Arimo", "DejaVu Sans"}
    used = {}
    for t in fig.findobj(matplotlib.text.Text):
        if not t.get_text().strip():
            continue
        path = fm.findfont(t.get_fontproperties(), fallback_to_default=True)
        fam = fm.FontProperties(fname=path).get_name()
        used[fam] = used.get(fam, 0) + 1
    assert set(used) <= allowed, f"{name}: unexpected fonts {set(used) - allowed}"
    return {f"fonts_{name}": used}


# ----------------------------------------------------------------------------- export
def save_figure(fig, stem, png_name=None, make_tif=True):
    """matplotlib -> pdf/svg + 600-dpi opaque PNG; Pillow -> RGB LZW TIFF; verify."""
    import shutil
    from PIL import Image

    for ext in ("pdf", "svg"):
        fig.savefig(OUTDIR / f"{stem}.{ext}", bbox_inches="tight", pad_inches=0.03)

    png = OUTDIR / (png_name or f"{stem}_{DPI}dpi.png")
    fig.savefig(png, dpi=DPI, transparent=False, facecolor="white",
                bbox_inches="tight", pad_inches=0.03)

    rec = {"figure": stem}
    if make_tif:
        # Pillow seeks back to update the TIFF header -> stage on local disk, then copy
        tif_tmp = Path("/tmp") / f"{stem}.tif"
        with Image.open(png) as im:
            im.convert("RGB").save(tif_tmp, compression="tiff_lzw", dpi=(DPI, DPI))
        tif = OUTDIR / f"{stem}.tif"
        shutil.copy(tif_tmp, tif)
        mb = tif.stat().st_size / 1e6
        with Image.open(tif) as im:
            w_px, h_px = im.size
            dpi_x, dpi_y = im.info.get("dpi", (None, None))
            dpi_x, dpi_y = int(round(float(dpi_x))), int(round(float(dpi_y)))
            comp = im.info.get("compression", "")
            n_pages = getattr(im, "n_frames", 1)
            mode, bands = im.mode, im.getbands()
        rec.update({
            "file": tif.name, "width_in": round(w_px / dpi_x, 3),
            "height_in": round(h_px / dpi_y, 3), "dpi_x": dpi_x, "dpi_y": dpi_y,
            "mode": mode, "n_channels": len(bands), "has_alpha": ("A" in bands),
            "compression": comp, "n_pages": n_pages, "file_size_MB": round(mb, 3),
        })
        assert 2.63 <= rec["width_in"] <= 7.5, f'{stem}: width {rec["width_in"]} in'
        assert rec["height_in"] <= 8.75, f'{stem}: height {rec["height_in"]} in'
        assert 300 <= dpi_x <= 600 and dpi_x == dpi_y, f"{stem}: dpi {dpi_x}/{dpi_y}"
        assert mode == "RGB" and len(bands) == 3 and not rec["has_alpha"]
        assert comp == "tiff_lzw" and n_pages == 1 and mb < 10
    return rec


# ----------------------------------------------------------------------------- renderers
def _panel_legend_handles(kind):
    if kind == "A":
        return [
            Line2D([0], [0], marker="o", linestyle="none", markersize=7.0,
                   markerfacecolor=C_OBS, markeredgecolor=C_OBS, markeredgewidth=0.9,
                   label="Observed difference"),
            Patch(facecolor=C_BAND, alpha=0.35, edgecolor="#4D4D4D", linewidth=0.75,
                  label="95% null range"),
        ]
    return [
        Patch(facecolor=C_BAND, alpha=0.35, edgecolor="#4D4D4D", linewidth=0.75,
              label="95% null range"),
        Line2D([0], [0], marker="|", linestyle="none", markersize=9,
               markeredgewidth=1.4, markeredgecolor=C_OBS, label="Null mean"),
    ]


def render_combined(d5, pb):
    """Two vertically stacked panels sharing identical horizontal geometry so the
    plotting-region edges and right-hand FDR-q columns align across panels."""
    fig = plt.figure(figsize=(7.1, 6.0))
    subA, subB = fig.subfigures(2, 1, height_ratios=[2.90, 2.55], hspace=0.10)
    for sub in (subA, subB):
        sub.set_facecolor("white")

    # Panel A (identical inner geometry to Panel B -> aligned q columns)
    gsA = subA.add_gridspec(1, 2, width_ratios=[3.55, 0.95], left=0.225, right=0.99,
                            top=0.925, bottom=0.32, wspace=0.03)
    axA, axAq = subA.add_subplot(gsA[0]), subA.add_subplot(gsA[1])
    draw_panel_a(axA, axAq, d5)
    subA.legend(handles=_panel_legend_handles("A"), loc="lower left", frameon=False,
                ncol=2, fontsize=10, bbox_to_anchor=(0.225, 0.05),
                handletextpad=0.5, columnspacing=1.2)
    subA.text(0.022, 0.965, "A", ha="left", va="top", fontsize=12, fontweight="bold")

    # Panel B
    gsB = subB.add_gridspec(1, 2, width_ratios=[3.55, 0.95], left=0.225, right=0.99,
                            top=0.855, bottom=0.345, wspace=0.03)
    axB, axBq = subB.add_subplot(gsB[0]), subB.add_subplot(gsB[1])
    draw_panel_b(axB, axBq, pb)
    subB.legend(handles=_panel_legend_handles("B"), loc="lower left", frameon=False,
                ncol=2, fontsize=10, bbox_to_anchor=(0.225, 0.03),
                handletextpad=0.5, columnspacing=1.2)
    subB.text(0.022, 0.955, "B", ha="left", va="top", fontsize=12, fontweight="bold")
    return fig


def render_panel_b_standalone(pb):
    fig = plt.figure(figsize=(7.1, 2.35))
    gs = fig.add_gridspec(1, 2, width_ratios=[3.55, 0.95], left=0.225, right=0.99,
                          top=0.80, bottom=0.335, wspace=0.03)
    ax, axq = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    draw_panel_b(ax, axq, pb)
    fig.legend(handles=_panel_legend_handles("B"), loc="lower left", frameon=False,
               ncol=2, fontsize=10, bbox_to_anchor=(0.225, 0.03),
               handletextpad=0.5, columnspacing=1.2)
    return fig


# ----------------------------------------------------------------------------- table
def write_plotting_table(pb, mode_label):
    prov = ("data/figure3_ld/frozen/Table_final_five_axis_sensitivity.tsv "
            "(SHA-256 03a82ba0…); exact-exclusion row cross-checked against "
            "Table_final_five_axis_exact_gps_minus_d.tsv (SHA-256 1c5375d1…); "
            "block_standard_gps band verified against the null matrix "
            "null_axis_percentages_block_gps_matched.tsv (SHA-256 35731bcb…); "
            f"plotted from {mode_label}")
    out = pd.DataFrame({
        "specification": pb.specification,
        "display_label": pb.display_label,
        "observed_pct": pb.observed_pct,
        "null_mean_pct": pb.null_mean,
        "null_q025_abs": pb.null_lo.round(5),
        "null_q975_abs": pb.null_hi.round(5),
        "q_BH5": pb.q_BH5,
        "n_replicates": pb.B.astype(int),
        "band_definition": ("central 95% of the empirical null distribution of T2D "
                            "annotation percentages (absolute scale); not a confidence interval"),
        "q_definition": ("Benjamini-Hochberg across the five primary disease axes, "
                         "computed separately within each specification"),
        "source_file": prov,
    })
    out.to_csv(OUTDIR / "Fig3B_T2D_reference_models_source_data.tsv", sep="\t",
               index=False)
    return out


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--from-recomputed", action="store_true",
                    help="plot from results/figure3_ld/recomputed/ instead of the "
                         "frozen tables")
    args = ap.parse_args()

    OUTDIR.mkdir(parents=True, exist_ok=True)
    if args.from_recomputed:
        gate_recomputed()
        src = RECOMPUTED
        mode_label = "recomputed tables (run_analysis.py --full)"
    else:
        gate_frozen()
        src = FROZEN
        mode_label = "frozen tables (data/figure3_ld/frozen/)"

    d5 = load_panel_a(src)
    print("Panel A assertion battery: PASS (cross-checks, null-matrix quantiles, "
          "spot values)")
    null_mat = pd.read_csv(src / "null_axis_percentages_block_gps_matched.tsv", sep="\t")
    pb = load_panel_b(src, null_mat)
    print("Panel B validation battery: PASS (3 rows, identical observed %, B=10,000, "
          "exact-ledger exclusion, band transformation, q display, geometry)")

    verification = []
    fig = render_combined(d5, pb)
    verification.append({**audit_figure(fig, "Fig3"), **audit_fonts(fig, "Fig3"),
                         **save_figure(fig, "Fig3", png_name="Fig3_preview.png")})
    plt.close(fig)
    print(f"  Fig3: rendered and verified {verification[-1]}")

    figB = render_panel_b_standalone(pb)
    verification.append({**audit_figure(figB, "Fig3B_T2D_reference_models"),
                         **audit_fonts(figB, "Fig3B_T2D_reference_models"),
                         **save_figure(figB, "Fig3B_T2D_reference_models",
                                       png_name="Fig3B_T2D_reference_models.png",
                                       make_tif=False)})
    plt.close(figB)
    print(f"  Fig3B standalone: rendered and verified {verification[-1]}")

    tab = write_plotting_table(pb, mode_label)
    print(f"  plotting table: {len(tab)} rows -> Fig3B_T2D_reference_models_source_data.tsv")

    (OUTDIR / "verification_tiff.json").write_text(json.dumps(verification, indent=2))
    print("\nVerification summary:")
    for rec in verification:
        print(" ", json.dumps(rec))


if __name__ == "__main__":
    main()
