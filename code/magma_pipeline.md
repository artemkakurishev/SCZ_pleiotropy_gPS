# MAGMA commands as executed
Source: /mnt/results/ledger_v1/ledger_3_magma.tsv, column `command_full` (27 runs). MAGMA commands as executed; transcribed verbatim, unedited.

```bash
magma --bfile ref/g1000_eur --gene-annot scz3_annot.genes.annot --pval scz3_snp_pval.txt ncol=N --gene-model snp-wise=mean --out scz3_gene
```

```bash
magma --gene-results scz3_gene.genes.raw --gene-covar gene_property_log_gps_nonzero.txt --model direction=pos --out scz3_gene_property_modelA
```

```bash
magma --gene-results scz3_gene.genes.raw --gene-covar gene_property_log_gps_zeros.txt --model direction=pos --out scz3_gene_property_modelB
```

```bash
magma --bfile ref/g1000_eur --gene-annot bip2021_annot.genes.annot --pval bip2021_snp_pval.txt ncol=N --gene-model snp-wise=mean --out bip2021_gene
```

```bash
magma --gene-results bip2021_gene.genes.raw --gene-covar gene_property_log_gps_nonzero.txt --model direction=pos --out bip2021_gene_property_modelA
```

```bash
magma --gene-results bip2021_gene.genes.raw --gene-covar gene_property_log_gps_zeros.txt --model direction=pos --out bip2021_gene_property_modelB
```

```bash
magma --bfile ref/g1000_eur --gene-annot mdd2025_annot.genes.annot --pval mdd2025_snp_pval.txt ncol=N --gene-model snp-wise=mean --out mdd2025_gene
```

```bash
magma --gene-results mdd2025_gene.genes.raw --gene-covar mdd_fixedN_covar_nonzero.txt --model direction=pos --out mdd2025_gene_property_modelA
```

```bash
magma --gene-results mdd2025_gene.genes.raw --gene-covar mdd_fixedN_covar_zeros.txt --model direction=pos --out mdd2025_gene_property_modelB
```

```bash
magma --bfile ref/g1000_eur --gene-annot height_annot.genes.annot --pval height_snp_pval.txt ncol=N --gene-model snp-wise=mean --out height_gene
```

```bash
magma --gene-results height_gene.genes.raw --gene-covar height_cov_A.txt (gene_property_log_gps_nonzero.txt) --model direction=pos --out height_gene_property_modelA
```

```bash
magma --gene-results height_gene.genes.raw --gene-covar height_cov_B.txt (gene_property_log_gps_zeros.txt) --model direction=pos --out height_gene_property_modelB
```

```bash
magma --bfile ref/g1000_eur --gene-annot ibd2017_annot.genes.annot --pval ibd2017_snp_pval.txt ncol=N --gene-model snp-wise=mean --out ibd2017_gene
```

```bash
magma --gene-results ibd2017_gene.genes.raw --gene-covar gene_property_log_gps_nonzero.txt --model direction=pos --out ibd2017_gene_property_modelA
```

```bash
magma --gene-results ibd2017_gene.genes.raw --gene-covar gene_property_log_gps_zeros.txt --model direction=pos --out ibd2017_gene_property_modelB
```

```bash
magma --bfile ref/g1000_eur --gene-annot adhd2022_annot.genes.annot --pval adhd2022_snp_pval.txt ncol=N --gene-model snp-wise=mean --out adhd2022_gene
```

```bash
magma --gene-results adhd2022_gene.genes.raw --gene-covar gene_property_log_gps_nonzero.txt --model direction=pos --out adhd2022_gene_property_modelA
```

```bash
magma --gene-results adhd2022_gene.genes.raw --gene-covar gene_property_log_gps_zeros.txt --model direction=pos --out adhd2022_gene_property_modelB
```

```bash
magma --bfile ref/g1000_eur --gene-annot pd2019_annot.genes.annot --pval pd2019_snp_pval.txt ncol=N --gene-model snp-wise=mean --out pd2019_gene
```

```bash
magma --gene-results pd2019_gene.genes.raw --gene-covar gene_property_log_gps_nonzero.txt --model direction=pos --out pd2019_gene_property_modelA
```

```bash
magma --gene-results pd2019_gene.genes.raw --gene-covar gene_property_log_gps_zeros.txt --model direction=pos --out pd2019_gene_property_modelB
```

```bash
magma --bfile ref/g1000_eur --gene-annot ocd2025_annot.genes.annot --pval ocd2025_snp_pval.txt ncol=N --gene-model snp-wise=mean --out ocd2025_gene
```

```bash
magma --gene-results ocd2025_gene.genes.raw --gene-covar gene_property_log_gps_nonzero.txt --model direction=pos --out ocd2025_gene_property_modelA
```

```bash
magma --gene-results ocd2025_gene.genes.raw --gene-covar gene_property_log_gps_zeros.txt --model direction=pos --out ocd2025_gene_property_modelB
```

```bash
magma --bfile ref/g1000_eur --gene-annot asd2019_annot.genes.annot --pval asd2019_snp_pval.txt ncol=N --gene-model snp-wise=mean --out asd2019_gene
```

```bash
magma --gene-results asd2019_gene.genes.raw --gene-covar gene_property_log_gps_nonzero.txt --model direction=pos --out asd2019_gene_property_modelA
```

```bash
magma --gene-results asd2019_gene.genes.raw --gene-covar gene_property_log_gps_zeros.txt --model direction=pos --out asd2019_gene_property_modelB
```
