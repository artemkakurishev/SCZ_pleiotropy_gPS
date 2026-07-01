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
│   ├── 03_fig3_bootstrap_fdr5.py       Figure 3 — gPS-matched bootstrap stats
│   ├── 03_fig3_make_figure.py          Figure 3 — render bootstrap figure
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
├── metadata/
│   ├── inputs_manifest.tsv             (per-input provenance + MD5 ground truth)
│   └── requirements.txt                (full transitive freeze; see ../requirements.txt for top-level pins)
├── results/                            (auto-generated; one .pdf, .svg, .png + data.tsv per figure)
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
| All other `*_data.tsv`        | **Yes** (deterministic from inputs)                          |
| `figure3_bootstrap_fdr5.{pdf,png,svg}` | **Yes** (no label repulsion, fixed layout) |
| `figure1_*`, `figure5_*` figures | **Yes** (deterministic layout)                           |
| `figure4_schema_targets.{pdf,png,svg}` | **No, by design.** `adjustText` converges to slightly different label layouts across runs. Visual content and numeric data are identical; only label positions drift by a few pixels. |
| `figure2_*` | **Yes** (no random ops; deterministic hexbin layout)        |

The brief required "numeric/vector parity, NOT pixel parity"; Figure 4's
label-repulsion drift falls within that allowance.

Files that intentionally drift are flagged in `MANIFEST.tsv` with
`strand=visual_only` (currently just the three Fig 4 image files). The MANIFEST snapshot still records the current MD5s of these files, but a downstream automated MD5 audit should treat `visual_only` rows as informational, not contractual. All TSVs, the bootstrap stats, scripts, and other rendered figures remain MD5-stable.

## Notes on Figure 1 trait order

The 9-trait forest in Figure 1 displays traits in the **published manuscript
order** (Height → IBD → MDD → SCZ → ADHD → BD → OCD† → ASD → PD), which is
NOT strictly descending β_std (a strict descending sort would put IBD before
Height, since β_std(IBD)=0.182 > β_std(Height)=0.164). The published order
groups the polygenic positive control (Height) first, then the non-psychiatric
disease (IBD), then descending β_std within the psychiatric block, and finally
the neurological/neuro-developmental traits. This matches the published
Figure 1 verbatim; see `code/01_fig1_gene_property.py` docstring for the
explicit ordering decision.

## Notes on Figure 4 filtering

The published Figure 4 panel renders 29 of the 34 SCHEMA FDR<0.1 PTV genes;
the 5 SCHEMA genes absent from MAGMA (`GRIA3`, `SLF2`, `H1-4`, `MAGEC1`,
`EIF2S3`) are excluded — they have no gene-level Z-score to position
against. The filter is `schema_category=='PTV' & schema_fdr<0.1 &
scz_source!='absent'` and is documented in
`inputs/SCHEMA_absent_from_MAGMA.md`.

## Notes on Figure 5 SE column

`loo_gps_sensitivity_frontiers.tsv` (Frontiers Supplementary Table S3)
reports `SE` as the standard error of `beta_std` directly (i.e. it is
already on the standardized scale), so 95% CIs are computed as
`beta_std ± 1.96 × SE` without any rescaling. This differs from Figure 1's
`gp_results_big9.tsv`, where `se` is the SE of the unstandardized
β (see the `01_fig1_gene_property.py` docstring for the per-trait scaling
note).

## Notes on the LOEUF input (Figure 2)

`inputs/loeuf_gnomad_v41.tsv` is a tidy gene-level slice of the public
`gnomad.v4.1.constraint_metrics.tsv` table
(`gs://gcp-public-data--gnomad/release/4.1/constraint/`, ETag
`14df4b2acb581fcbbb2a82a3a555fd35`, dated 2024-04-18). The slice keeps
the Ensembl-anchored branch of the gnomAD file (gene_id starts with
`ENSG`), filters to MANE-select with a canonical-Ensembl fallback for the
1,142 non-MANE genes, drops the 957 genes that carry any `constraint_flags`
(equivalent to the v4.1.1 browser-side filter), and renames columns to
`ensembl_id`, `gene_symbol`, `loeuf`, `loeuf_lower`, `loeuf_rank`,
`loeuf_decile`, `lof_oe`, `pli`, `lof_obs`, `lof_exp`. The result is
17,666 unique Ensembl gene rows covering 96.7% of the SCZ MAGMA gene
list. After the SCZ ∩ BIP ∩ LOEUF inner merge in Figure 2, 7,653 genes
are plotted (vs the manuscript's reported ≈7,514; the +139 difference
is the price of including the 185 non-MANE canonical Ensembl rows and
falls inside the script's ±200-gene soft-warning tolerance).

## Contact

Artemiy O. Kurishev — Russian Mental Health Research Center / Engelhardt
Institute of Molecular Biology, Moscow.
