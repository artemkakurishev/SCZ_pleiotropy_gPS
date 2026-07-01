# External inputs required by the pipeline

The nine files listed in `metadata/inputs_manifest.tsv` are the **only**
inputs consumed by the five `code/0X_*.py` scripts. Every script asserts
its required MD5s at the top of `main()` and exits non-zero before any
computation if a file is missing or has the wrong hash (fail-fast policy).

## Bit-for-bit preserved core inputs (4)

- `gene_annotation_table.tsv` — SCZ MAGMA gene-level Z + gPS + tier.
  md5 = `6b308a7d300f9de33f6f458fdab9c79a`. Used by Fig 2, Fig 3 (frozen).
- `bip_gene_annotation_table.tsv` — BIP MAGMA gene-level Z + gPS + tier.
  md5 = `c944322ea7c6a489a7229411162761c0`. Used by Fig 2.
- `schema_pav_gps_table.tsv` — SCHEMA × gPS × PAV-fraction table.
  md5 = `b9c5f67a83510d0d87f94981914f1a40`. Used by Fig 4.
- `loeuf_gnomad_v41.tsv` — Per-Ensembl gnomAD v4.1 LOEUF table.
  md5 = `58c78cec40fb791268aeda6df346fbf4`. Used by Fig 2.
  Extracted from gnomAD v4.1 `constraint_metrics.tsv`: MANE-select
  transcript preferred, canonical as fallback; flag-filtered rows removed;
  17,666 gene rows. Columns: `ensembl_id`, `gene_symbol`, `transcript_id`,
  `loeuf`, `loeuf_lower`, `loeuf_rank`, `loeuf_decile`, `lof_oe`, `pli`,
  `lof_obs`, `lof_exp`. Fig 2 consumes only `ensembl_id` and `loeuf`.

## Supplementary inputs extracted from the deposit / supplementary archive (3)

- `gp_results_big9.tsv` — 9-trait Model-A/B MAGMA gene-property table with
  `beta_std` and `p_one_sided`. md5 = `b83619c51e4a433c39feea9e259852b7`.
  Used by Fig 1.
- `loo_gps_sensitivity_frontiers.tsv` — Leave-one-out sensitivity for SCZ3
  and BIP2021 across four scenarios with `gps_variance_reduction_pct`.
  md5 = `16b15719d497627f8cc8dee92f4ddedf`. Used by Fig 5.
- `removed_disease_ledger.tsv` — One row per disease term removed in each
  LOO scenario. md5 = `8f50071e5ad3e832b66abbb403a583aa`. Audit-only
  (cited by Fig 5 in the script docstring; not numerically used by render).

## External, MD5-gated user-supplied inputs (2)

- `disease_ta_index_pandas.csv` — Gentropy gPS source (commit `fec6427`).
  md5 = `03b7c1e7fc211cd69c6d00141bff20b0`. Used by Fig 3 (frozen).
- `l2g_diseases_full.csv` — Open Targets 26.03 L2G credible-set predictions.
  md5 = `f426e77b1bb71e4136d2577dc4dd6215`. Used by Fig 3 (frozen).
