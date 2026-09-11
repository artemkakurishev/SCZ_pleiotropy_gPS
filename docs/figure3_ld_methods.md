# Figure 3 A/B — LD-block stratified null analysis (methods as executed)

This document records the method **as actually executed** for the revision of
PONE-D-26-36293, the exact commands that reproduce it, and the provenance of
every input. It accompanies `code/figure3_ld/` and supersedes the gene-wise
bootstrap for Figure 3 (see "Which calculation creates Figure 3" below).

All constants below were transferred from the originally executed analysis code
(`code/figure3_ld/sampler.py`, MD5 `90a0ad89aba97e0b735c074fa9968270`;
`code/figure3_ld/sampler_v5.py`, MD5 `302eceb92675ef883356c0e8b74721dd`;
shipped `code/03_fig3_bootstrap_fdr5.py`, MD5 `d75c7dd6b2a2fa7d8003860eb4616b11`)
and are pinned in `code/figure3_ld/config.json`. **No parameter was tuned to
match the frozen outputs.**

## Which calculation creates Figure 3

| | Historical (kept, not deleted) | Current (creates Figure 3 A/B) |
|---|---|---|
| Script | `code/03_fig3_bootstrap_fdr5.py` + `code/03_fig3_make_figure.py` | `code/figure3_ld/run_analysis.py` + `code/figure3_ld/plot_figure3.py` |
| Null unit | individual genes | LDetect EUR LD blocks |
| Replicates | 1,000 | 10,000 |
| Output | `results/figure3_bootstrap_*` | `results/figure3_ld/` |
| Status | submission-time analysis, preserved unchanged | revision analysis |

The old gene-wise bootstrap is retained as the historical record of the
submitted analysis. The committed `results/figure3_bootstrap_stats.tsv` is the
**submission-time** table; it differs slightly from the final re-executed
gene-wise null (the `submitted_genewise_gps_matched` row of the frozen
sensitivity table), e.g. T2D null mean 11.7336 → 11.7338, reference-interval
lower bound 0.0984 → 0.1969, p 0.032 (centred convention) → 0.037962 (R2#7
convention), q 0.040 → 0.0475; Prostate Ca null mean 6.0348 → 6.0386. The
frozen revision tables are authoritative; the historical TSV is kept unchanged.
These differences are itemised in
`results/figure3_ld/recomputed/comparison_report.md` (informational section,
not a gate).

## Commands

Two modes, both from the repository root:

**(a) Fast — Figure 3 from the frozen tables (seconds).**
Verifies the frozen tables against `data/figure3_ld/frozen/SHA256SUMS.txt`,
then renders. No resampling is performed.

```bash
python code/figure3_ld/plot_figure3.py
```

**(b) Full — re-run the entire analysis, compare to frozen, then render
(~1–2 min).** Recomputes every null from the raw inputs, writes
`results/figure3_ld/recomputed/`, compares each table against the frozen
version (halting on any discrepancy beyond storage precision), and only then
is the figure allowed to read the recomputed tables.

```bash
python code/figure3_ld/run_analysis.py --full
python code/figure3_ld/plot_figure3.py --from-recomputed
```

`plot_figure3.py --from-recomputed` refuses to plot if
`results/figure3_ld/recomputed/comparison_report.md` is absent or contains a
`**FAIL**` line, so an unverified recomputation can never silently become a
figure.

## Method as executed

**Genome build.** GRCh37/hg19 throughout (autosomes identical between the two).
LDetect EUR blocks and NCBI37.3 gene locations are both GRCh37. Open Targets
(GRCh38) coordinates were never used for block assignment.

**Universe and cohorts.** Autosomal protein-coding genes, excluding chrX, the
MHC (chr6:25–35 Mb) and 17q21.31 (chr17:43.5–44.9 Mb); duplicate Ensembl IDs
resolved by keeping the maximum MAGMA Z. Universe n = 7,707. Risk cohort:
MAGMA Z ≥ 3.09 → 1,022 genes, of which 6 fall in shortfall bins and are
dropped (ENSG00000183527, ENSG00000132394, ENSG00000196628, ENSG00000112182,
ENSG00000187323, ENSG00000138821), leaving the target n = 1,016. Background
pool: Z ≤ 1.0 → **n = 3,605** (the manuscript text printed 3,614; the executed
value is 3,605 and is what all frozen and recomputed numbers use).

**Disease-axis membership (L2G ≥ 0.5).** A gene belongs to an axis if it has an
Open Targets L2G score ≥ 0.5 whose parsed `diseaseIds` intersect the axis ID
set. This is the **disease-membership threshold** and is distinct from the gPS
counting rule below. Axis IDs: T2D = {MONDO_0005148}; CAD = {EFO_0001645};
Hypertension = {EFO_0000537}; Asthma = {MONDO_0004979}; Prostate Ca =
{EFO_0001663}; Atrial Fib = {EFO_0000275}; Gout = {EFO_0004274}; OA =
{EFO_0004616, EFO_1000786, MONDO_0005178}. Axis sizes in the universe: T2D 869,
CAD 272, Hypertension 356, Asthma 291, Prostate Ca 382 (primary five);
Atrial Fib 323, Gout 279, OA 411 (exploratory three).

**gPS counting rule.** gPS is the `uniqueDiseases` column of
`disease_ta_index_pandas.csv` (Open Targets / Gentropy), i.e. the count of
distinct diseases per gene at that table's native ≈0.1 score floor — **not**
the count at L2G ≥ 0.5. A provenance gate in the original run verified that the
implemented gPS equals the per-gene counted-disease-set size for every pool and
risk gene.

**LD blocks and gene-to-block assignment.** Berisa & Pickrell (2016) LDetect
EUR partition (`fourier_ls-all.bed`, bitbucket.org/nygcresearch/ldetect-data;
raw BED MD5 `10167be5c416895f27aedf7485db05ae`), 1,703 autosomal blocks. Each
gene is assigned to the block containing the **midpoint of its MAGMA window**
(the `gene_annotation_table.tsv` start/stop columns are already the 35 kb
upstream / 10 kb downstream window bounds — verified 7,793/7,793 against
`NCBI37.3.gene.loc`), with a nearest-block fallback within 1 Mb across gaps.
100% of the 7,707 universe genes are assigned. Observed risk genes occupy 433
blocks; the pool occupies 1,145 blocks.

**Block sampler.** Blocks are binned by the mean log2(gPS+1) of their pool
genes and assigned to the nearest risk-set integer-gPS bin in log2(gPS+1)
space. A crossing block is included with probability (target − achieved) /
block size (unbiased randomized rounding). Neighbour borrowing occurs only on
supply exhaustion, from the nearest bin(s), with a single crossing decision
across the borrow sequence. A draw is accepted if |total genes − 1016| ≤ 25,
else redrawn (cap 50 × n_sets). There is no KS gate (KS is computed as a
diagnostic only).

**The three Panel B reference models** (all LD-block, B = 10,000):
1. `block_size_only` — size-matched only, no gPS conditioning.
2. `block_standard_gps` — gPS-matched (PRIMARY).
3. `block_exact_gps_minus_d` — exact disease-ledger exclusion: per-gene counted
   disease sets are taken from `l2g_diseases_full.csv` exploded with no
   additional score filter, the axis disease IDs are removed
   (gPS^(−d)_i = |D_i \ axis_ids|, clipped at 0), bins and per-bin gene-count
   targets are rebuilt, and 10,000 replicates are drawn per axis. For T2D the
   removed ID is MONDO_0005148; this is the exact exclusion, not the L2G ≥ 0.5
   approximation (which gave a residual of ≈ +0.59 pp vs the exact +0.523 pp).

**Seeds.** Base `SeedSequence(20260527)`. Because spawn-key identity makes
`spawn(6)[i] == spawn(16)[i]` for i < 6: gene-wise null =
`default_rng(20260527)`; block standard-gPS 10k = `spawn(6)[2]` (=
`spawn(16)[2]`); size-only 10k = `spawn(16)[6]`; exact gPS^(−d) 10k =
`spawn(16)[7..11]` for T2D, CAD, Hypertension, Asthma, Prostate Ca.

**p-values and FDR.** Two-sided p = min(1, 2·min(b_le + 1, b_ge + 1)/(B + 1)),
where b_le / b_ge count null replicates ≤ / ≥ the observed statistic; the floor
at B = 10,000 is 2×10⁻⁴ (never 0). FDR is Benjamini–Hochberg across the **5
primary axes, computed separately within each specification** (q_BH5); the
exploratory axes are excluded from the primary family, and q is **not**
recomputed across the three T2D rows of Panel B. Monte-Carlo intervals are
Clopper–Pearson 95% intervals on min(b_le, b_ge)/B propagated through the
p formula and the BH step (p_mc_lo/hi, q_BH5_mc_lo/hi).

## Panel A / Panel B interval definitions

The grey bands are **central 95% empirical-null reference intervals** — the
2.5th to 97.5th percentile of the matched null distribution. They are **not
confidence intervals**, are not forced to be symmetric, and are not centred on
the observed values.

- **Panel A** (delta scale, percentage points): point = observed % − matched
  null mean; band = the centred null quantiles [Q0.025(N) − null_mean,
  Q0.975(N) − null_mean]; solid zero reference line; right column = q_BH5.
- **Panel B** (absolute scale): band = [observed % − ref95_hi, observed % −
  ref95_lo] = the central 95% of the model's null distribution of T2D
  annotation percentages; internal black tick = null mean; dashed vertical line
  = observed T2D percentage (13.39%, fixed across the three models); right
  column = q_BH5 within each specification.

## Why the T2D conclusion changes: mean shift, not null widening

The move from the gene-wise to the LD-block null changes the T2D result, but
**not** because the null distribution widens. The null standard deviation is
nearly constant across the gPS-conditioned specifications; what moves is the
null **mean**:

| Specification | Null mean (%) | Null SD (%) | Residual (pp) | q_BH5 |
|---|---|---|---|---|
| gene-wise + gPS (historical, B=1,000) | 11.7338 | 0.7845 | +1.652 | 0.0475 |
| block, size-only (no gPS) | 10.3558 | 0.9456 | +3.030 | 0.0040 |
| block + standard gPS (PRIMARY) | 12.4100 | 0.7863 | +0.976 | 0.2722 |
| block + exact gPS^(−T2D) | 12.8632 | 0.7931 | +0.523 | 0.5097 |

Reading the decomposition: conditioning on LD-block structure *without* gPS
(size-only) lowers the null mean to 10.36 % and slightly widens the SD to
0.9456, so the observed 13.39 % looks *more* enriched (+3.03 pp, q = 0.004).
Adding standard gPS matching raises the null mean to 12.41 % at essentially
unchanged SD (0.7863), cutting the residual to +0.98 pp (q = 0.272). Excluding
T2D from the gPS ledger raises the null mean further to 12.86 % (SD 0.7931),
leaving +0.52 pp (q = 0.510). The attenuation is therefore driven by the
upward shift of the null mean as block structure and then T2D-ledger exclusion
are added — the variance barely changes. It would be incorrect to describe the
conclusion change as "the LD correction widened the null"; it did not.

## Matching diagnostics, occupancy, and filter audit

These diagnostics mirror the authoritative methods prose
(`methods_final_ld_gps.md`, sha256
`7acbd3ad64dd10695824890ecb94130e3d4c6ddb60464103240f6cdd29e02378`); the
realised gene-level values are reproduced cell-for-cell in
`Table_realised_genelevel_gps.tsv` (see "Reproduction result").

**Matching adequacy.** Across the 10,000 primary-null sets, the realised
gene-level gPS mean matched the risk set to within ± 0.05 log2(gPS+1) units in
**100%** of sets (requirement ≥ 95%). At the distribution level the match is
approximate: the median total-variation distance between the realised-null and
risk-set gPS distributions was **0.165**, with **17 of 39** occupied bins lying
outside the 95% bin-level envelope and a net upward borrowing of **16.1%** from
neighbouring bins. The null is therefore described as matched on **mean** gPS,
not on the full gPS distribution.

**Block occupancy.** The observed risk genes occupy **433** LD blocks, whereas
primary-null sets occupy **398.0** blocks on average (95% range **377–419**).
An indicative, linearity-based estimate places the effect of this occupancy
difference on the T2D residual at approximately **+0.45 percentage points**
toward the null; no qualitative five-axis conclusion is affected.

**Acceptance-filter audit.** Removing the ± 25-gene acceptance filter changes
the null standard deviations by ≤ 1.18%; the filter has no material effect on
inference.

## Reproduction result (this deposit)

`python code/figure3_ld/run_analysis.py --full` was executed end-to-end. All
**six** recomputed tables match the frozen tables in `data/figure3_ld/frozen/`
— five byte-identical at storage precision and the sixth,
`Table_realised_genelevel_gps.tsv`, **byte-identical** (file sha256
`b968d9f4…`; 429/429 cells identical across the 39 occupied gPS bins). Every
comparison check PASSes (0 FAIL; see
`results/figure3_ld/recomputed/comparison_report.md`).

The realised gene-level table is built from the per-set gene identities
recorded by an **additional** `run_block_null_rec` pass run alongside the
untouched primary `run_block_null` (same seed, spawn-key `(2,)`). A named gate
(`recorder_bit_identity`) asserts the recorder's null matrix is bit-identical
to the primary's (max|Δsizes| 0, max|Δmeans| 0, max|Δpcts| 0) before the
recorded gene identities are used; on any drift the run halts and reports the
magnitude. The primary path is never modified to obtain the diagnostic.

The figure rendered from the recomputed tables is identical in plotted values
to the figure rendered from the frozen tables (the assertion battery locks
every plotted number); the TIFF is not byte-identical to the submission figure
only because genuine Arial is not installed in this environment and the
metric-compatible Arimo is used — geometry shifts by ≤ 4 px, all plotted
values unchanged.

## Inputs and provenance

Derived LD inputs (this deposit, `data/figure3_ld/inputs/`):
- `gene_block_map.tsv` (MD5 `1b95b0819be4ed2f2d4e2853d57aef90`) — per-gene
  block assignment.
- `ldetect_EUR_blocks.tsv` (MD5 `929ffdc4b95beb6ba54df88f904aa0d1`) — 1,703
  blocks with assigned block_id EUR_0001..EUR_1703.
- `NCBI37.3.gene.loc` (MD5 `7f6ffd6abb0b02c94a629e74c036d5f6`) — MAGMA gene
  locations, GRCh37.

See `data/figure3_ld/inputs/provenance.md` for the full source/build record and
`data/figure3_ld/frozen/SHA256SUMS.txt` for the frozen-table hashes. The three
shared manuscript inputs (`gene_annotation_table.tsv`,
`disease_ta_index_pandas.csv`, `l2g_diseases_full.csv`) are documented in
`metadata/inputs_manifest.tsv`. No individual-level genetic or clinical data
are included.
