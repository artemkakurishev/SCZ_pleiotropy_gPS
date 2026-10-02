# Reproducibility Pipeline — Pleiotropy-Matched Framework for SCZ/BD Risk Genes

Refactored "one script = one figure" Python 3.11 pipeline accompanying the
manuscript:

> Kurishev AO, Marilovtseva EV, Semina EV, Chaika YA, Golimbet VE.
> *A Pleiotropy-Matched Framework Extends Cross-Disease Interpretation of
> Schizophrenia and Bipolar Disorder Risk Genes.* (Frontiers preprint,
> 2025-06-26; PeerJ submission in preparation.)

The original analysis lived in a monolithic Jupyter notebook
(`execution_trace_magma-worker.ipynb`, 573 cells across multiple machine
restarts) and a pair of frozen bootstrap scripts. This deposit replaces
that with six small, MD5-gated Python scripts — one per figure — that fail
fast on missing or altered inputs and emit both vector (PDF/SVG) and
raster (PNG, 300 dpi) renderings together with the underlying numeric
data tables.

## Structure

```
pleiotropy_deposit/
├── code/
│   ├── 01_fig1_gene_property.py        Figure 1 — 9-trait β_std forest
│   ├── 02_fig2_hexbin_gradient.py      Figure 2 — LOEUF × gPS hexbin (SCZ−BD ΔZ)
│   ├── 03_fig3_bootstrap_fdr5.py       Figure 3 — gPS-matched bootstrap stats (HISTORICAL)
│   ├── 03_fig3_make_figure.py          Figure 3 — render bootstrap figure (HISTORICAL)
│   ├── figure3_ld/                     Figure 3 A/B — LD-block stratified null (CURRENT)
│   │   ├── run_analysis.py                 full re-run + frozen-table comparison
│   │   ├── plot_figure3.py                 figure rendering (fast / from-recomputed)
│   │   ├── config.json                     all analysis constants, transferred verbatim
│   │   ├── sampler.py                      LD-block sampler (original, unchanged)
│   │   └── sampler_v5.py                   v5 sampler extensions (original, unchanged)
│   ├── 04_fig4_schema_targets.py       Figure 4 — SCHEMA × gPS × MAGMA Z bubble plot
│   ├── 05_fig5_loo_sensitivity.py      Figure 5 — LOO line plot (SCZ vs BD)
│   ├── magma_pipeline.md               MAGMA commands as executed for the deposit
│   └── run_all.sh                      Wrapper: runs every script in order
├── inputs/
│   ├── gene_annotation_table.tsv               (SCZ3 MAGMA gene-level + gPS)
│   ├── bip_gene_annotation_table.tsv           (BIP2021 MAGMA gene-level + gPS)
│   ├── schema_pav_gps_table.tsv                (SCHEMA FDR<0.1 PTV, n=34)
│   ├── gp_results_big9.tsv                     (MAGMA gene-property, 9 traits × 2 models)
│   ├── loo_gps_sensitivity_frontiers.tsv       (Frontiers Suppl Table S3)
│   ├── removed_disease_ledger.tsv              (Frontiers Suppl Table S4)
│   ├── disease_ta_index_pandas.csv             (Open Targets / Gentropy gPS table)
│   ├── l2g_diseases_full.csv                   (Open Targets L2G credible-set predictions)
│   ├── loeuf_gnomad_v41.tsv                    gnomAD v4.1 LOEUF (Ensembl-anchored, for Figure 2)
│   ├── EXTERNAL_INPUTS_README.md               (per-file provenance walk-through)
│   └── SCHEMA_absent_from_MAGMA.md             (the 5 SCHEMA genes absent from MAGMA)
├── data/figure3_ld/
│   ├── inputs/                         (LD-block inputs: gene_block_map, ldetect_EUR_blocks,
│   │                                    NCBI37.3.gene.loc + provenance.md)
│   └── frozen/                         (frozen revision tables + SHA256SUMS.txt)
├── docs/
│   └── figure3_ld_methods.md           (Figure 3 A/B method as executed + commands)
├── metadata/
│   ├── inputs_manifest.tsv             (per-input provenance + MD5 ground truth)
│   └── requirements.txt                (full transitive freeze; see ../requirements.txt for top-level pins)
├── results/                            (auto-generated; one .pdf, .svg, .png + data.tsv per figure)
│   └── figure3_ld/                     (Figure 3 A/B: Fig3.tif/pdf/svg + recomputed/)
├── MANIFEST.tsv                        (per-file MD5 + role manifest for the deposit)
├── README.md                           (this file)
├── requirements.txt                    (Python 3.11 pin set; see metadata/ for full freeze)
└── requirements_frozen.txt             (verbatim deposit-time pip freeze, preserved)
```

## Quickstart

```bash
cd pleiotropy_deposit
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
cd code
bash run_all.sh
```

Outputs land in `results/`. Figures are emitted as PDF (vector), SVG
(editable text), and PNG (300 dpi raster). Underlying numeric data is
written as a tab-separated table next to each figure, e.g.
`figure1_gene_property_data.tsv`.

## Figure 3: which calculation was used before, and which creates it now

Figure 3 changed between submission and revision. Both pipelines are kept.

- **Historical (submission):** `code/03_fig3_bootstrap_fdr5.py` +
  `code/03_fig3_make_figure.py` — a **gene-wise** gPS-matched bootstrap
  (1,000 replicates), writing `results/figure3_bootstrap_*`. These files are
  preserved unchanged as the record of the submitted analysis. The committed
  `results/figure3_bootstrap_stats.tsv` is the submission-time table and
  differs slightly from the final re-executed gene-wise null (see
  `docs/figure3_ld_methods.md`).
- **Current (revision, creates Figure 3 A/B):** `code/figure3_ld/` — an
  **LD-block stratified** null (LDetect EUR blocks, 10,000 replicates) with
  three reference models (size-only / standard gPS / exact gPS excluding the
  axis disease). Full method: `docs/figure3_ld_methods.md`.

Two modes (from the repository root):

```bash
# (a) Fast — render Figure 3 from the frozen, SHA256-verified tables (seconds)
python code/figure3_ld/plot_figure3.py

# (b) Full — re-run the whole analysis, compare to the frozen tables, then render
python code/figure3_ld/run_analysis.py --full
python code/figure3_ld/plot_figure3.py --from-recomputed
```

The full re-run writes `results/figure3_ld/recomputed/` and halts on any
discrepancy with the frozen tables beyond storage precision; the
`--from-recomputed` figure mode refuses to plot from an unverified
recomputation. Executed end-to-end for this deposit, all six recomputed
tables match the frozen ones — five byte-identical at storage precision, and
the realised gene-level gPS diagnostic table
(`Table_realised_genelevel_gps.tsv`, 39 occupied bins) byte-identical
(429/429 cells). Figure 3 outputs land in
`results/figure3_ld/` (`Fig3.tif/.pdf/.svg`, a standalone Panel B, and the
Panel-B source table).

## MD5 fail-fast policy

Every script declares its inputs with an `INPUTS` dict that pins each
file to an expected MD5. Before any compute happens, `assert_inputs()`
verifies file existence and hash. On mismatch, the script exits 1 with a
diff message naming the file and the expected vs observed hash. This
prevents silent runs against an accidentally-substituted upstream file.

To accept a deliberate upstream change:
1. Compute the new MD5 with `md5sum path/to/file`.
2. Update the corresponding row in `metadata/inputs_manifest.tsv`.
3. Update the `INPUTS[...]['md5']` literal in the affected script.
4. Re-run.

## Reproducibility & determinism

| Output | Bit-stable across runs? |
|---|---|
| `figure3_bootstrap_stats.tsv` | **Yes** (NumPy `default_rng(seed=20260527)`, fixed N_ITER=1000) |
| `figure3_ld/recomputed/*.tsv` | **Yes** — byte-identical to `data/figure3_ld/frozen/` (SeedSequence(20260527), fixed replicate counts) |
| `figure3_ld/Fig3.{tif,pdf,svg}` | **Values yes, pixels no.** Plotted numbers are assertion-locked to the tables; raster pixels depend on the installed Arial-metric font (Arimo here), and PDF/SVG bytes embed a build timestamp. |
| `figure3_ld/Fig3B_..._source_data.tsv` | **Numeric content yes.** The `source_file` provenance string records the rendering mode (frozen vs recomputed); all numeric cells are identical across modes. |
| All other `*_data.tsv`        | **Yes** (deterministic from inputs)                          |
| `figure3_bootstrap_fdr5.{pdf,png,svg}` | **Yes** (no label repulsion, fixed layout) |
| `figure1_*`, `figure5_*` figures | **Yes** (deterministic layout)                           |
| `figure4_schema_targets.{pdf,png,svg}` | **No, by design.** `adjustText` converges to slightly different label layouts across runs. Visual content and numeric data are identical; only label positions drift by a few pixels. |
| `figure2_*` | **Yes** (no random ops; deterministic hexbin layout)        |

The brief required "numeric/vector parity, NOT pixel parity"; Figure 4's
label-repulsion drift falls within that allowance.

Files that intentionally drift are flagged in `MANIFEST.tsv` with
`strand=visual_only` (currently just the three Fig 4 image files). The MANIFEST snapshot still records the current MD5s of these files, but a downstream automated MD5 audit should treat `visual_only` rows as informational, not contractual. All TSVs, the bootstrap stats, scripts, and other rendered figures remain MD5-stable.

## Figure-specific reproduction notes

**Figure 1** — Trait order is the published display order
(Height → IBD → MDD → SCZ → ADHD → BD → OCD† → ASD → PD), not a strict
descending-β_std sort. Rationale: `code/01_fig1_gene_property.py` docstring.

**Figure 4** — 29 of the 34 SCHEMA FDR<0.1 PTV genes are plotted; the 5
without a MAGMA Z-score cannot be positioned and are removed by
`schema_category=='PTV' & schema_fdr<0.1 & scz_source!='absent'`.
Excluded genes: `inputs/SCHEMA_absent_from_MAGMA.md`.

**Figure 5** — `SE` in `loo_gps_sensitivity_frontiers.tsv` (Suppl.
Table S3) is already the SE of `beta_std`, so CIs are `beta_std ± 1.96 × SE`
with no rescaling. This differs from Figure 1, where `se` is on the
*unstandardized* β scale (see the `01_fig1_gene_property.py` docstring).

**Figure 2 (LOEUF input)** — After the SCZ ∩ BIP ∩ LOEUF merge, 7,653 genes
are plotted vs the manuscript's ≈7,514; the +139 difference is the price of
including non-MANE canonical-Ensembl rows and falls inside the script's
±200-gene tolerance. Full derivation of `loeuf_gnomad_v41.tsv`:
`inputs/EXTERNAL_INPUTS_README.md`.

## Contact

Artemiy O. Kurishev — Russian Mental Health Research Center / Engelhardt
Institute of Molecular Biology, Moscow.
