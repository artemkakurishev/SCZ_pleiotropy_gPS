# SCHEMA Genes Absent from the MAGMA Gene-Property Set

Five of the 34 SCHEMA PTV-significant genes (FDR < 0.1, Singh et al. 2022,
PMID 35396580) have no L2G credible-set prioritization in Open Targets
release 26.03 and are therefore absent from the MAGMA gene-property
analysis. They are listed below for full transparency.

Figure 4 (`code/04_fig4_schema_targets.py`) deliberately excludes these five
genes from the bubble plot because the y-axis (`magma_z_scz`) is undefined
for them. The script is gated by the MD5 of `schema_pav_gps_table.tsv`
(`b9c5f67a83510d0d87f94981914f1a40`); the absent rows in that source file
are identified by `scz_source == "absent"` and are excluded **before**
plotting, producing the 29-bubble figure that matches the published Figure 4.

| gene_symbol | ensembl_id      | schema_p   | schema_fdr | schema_OR | reason        |
|-------------|-----------------|------------|------------|-----------|---------------|
| GRIA3       | ENSG00000125675 | 5.98e-07   | 0.00155    | inf       | not_in_MAGMA  |
| SLF2        | ENSG00000119906 | 2.5e-05    | 0.0245     | 3.05      | not_in_MAGMA  |
| H1-4        | ENSG00000168298 | 5.84e-05   | 0.0447     | 7.03      | not_in_MAGMA  |
| MAGEC1      | ENSG00000155495 | 6.18e-05   | 0.0447     | 2.25      | not_in_MAGMA  |
| EIF2S3      | ENSG00000130741 | 8.23e-05   | 0.0467     | inf       | not_in_MAGMA  |

Source of truth: `inputs/schema_pav_gps_table.tsv`, rows where
`scz_source == "absent"`. A standalone supplementary
table with these five entries is also retained at
`data/intermediate/supp_table_S4_schema_absent_from_magma.tsv` in the
deposit for citation in the manuscript.
