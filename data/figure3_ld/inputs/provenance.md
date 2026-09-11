# LD block + gene-location provenance (data/figure3_ld/inputs/)

## LDetect EUR blocks (`ldetect_EUR_blocks.tsv`)
- Source: Berisa & Pickrell (2016, Bioinformatics 32:283-5), LDetect EUR partition
  `EUR/fourier_ls-all.bed`, downloaded from bitbucket.org/nygcresearch/ldetect-data.
- Raw BED MD5: `10167be5c416895f27aedf7485db05ae`; derived TSV (block_id assigned
  EUR_0001..EUR_1703 in chr/start order) MD5: `929ffdc4b95beb6ba54df88f904aa0d1`.
- 1,703 blocks, autosomal, hg19/GRCh37 coordinates, median 1.52 Mb.
- Assembly note: hg19 and GRCh37 are identical on autosomes at block scale; the EUR
  partition matches the g1000_eur MAGMA reference panel used for the SCZ gene-level
  analysis. Open Targets (GRCh38) coordinates were never used for block assignment.

## MAGMA gene locations (`NCBI37.3.gene.loc`)
- Source: canonical CTG/SURFsara MAGMA distribution,
  https://vu.data.surfsara.nl/index.php/s/Pj2orwuF2JYyKxq/download
  (link published on the official MAGMA page, cncr.nl/research/magma/; the legacy
  ctg.cncr.nl URLs now 301-redirect). NCBI 37.3 gene locations = GRCh37 build.
- MD5: `7f6ffd6abb0b02c94a629e74c036d5f6`; 19,427 genes; columns
  entrez_id, chr, start, stop, strand, symbol.

## Gene-to-block map (`gene_block_map.tsv`)
- MD5: `1b95b0819be4ed2f2d4e2853d57aef90`. One row per analysis-universe gene with
  its assigned LDetect EUR block_id.
- Assignment rule: block containing the midpoint of the gene's MAGMA window; the
  `gene_annotation_table.tsv` start/stop columns are already the 35 kb upstream /
  10 kb downstream window bounds (see cross-check below). Nearest-block fallback
  within 1 Mb across gaps. 100% of the 7,707 universe genes assigned.

## Build cross-check
- The analysis table `gene_annotation_table.tsv` (MD5 `6b308a7d300f9de33f6f458fdab9c79a`)
  start/stop coordinates equal `NCBI37.3.gene.loc` body coordinates extended by the
  MAGMA window convention (strand-aware 35 kb upstream / 10 kb downstream) for
  **7,793/7,793 genes exactly** (100% start+stop agreement), and for all 7,707/7,707
  analysis-universe genes. This confirms (i) both files are GRCh37, (ii) the table's
  start/stop columns already ARE the MAGMA 35/10 window bounds.
- Consequence: the MAGMA window for any gene is `[table_start, table_stop]`; the gene
  body is `[table_start + 35,000, table_stop - 10,000]` on the + strand
  (`[table_start + 10,000, table_stop - 35,000]` on the - strand); TSS = body start (+)
  / body stop (-).
- Strand agreement between the NCBI .loc file and the Ensembl GRCh37 lookup used in the
  v4 step-6 diagnostic: 99.91% (7 of 7,707 genes discordant; .loc strand is
  authoritative here because it is the MAGMA annotation source).

## Correction to the v4 step-6 diagnostic (recorded for transparency)
- v4's assignment-rule diagnostic treated the table start/stop as the gene *body* and
  added a further 35/10 kb, i.e. it double-extended the window, and placed the "TSS"
  35 kb outside the true body. Recomputed with the correct convention:
  TSS-vs-midpoint block reassignment 193/7,707 (v4 reported 329); MAGMA-window spans
  >1 block in 638/7,707 genes (v4 reported 835; corrected distribution
  {1: 7069, 2: 624, 3: 12, 4: 2}); risk-cohort subsets 41 and 102 (v4: 52 and 124).
- Direction of bias: v4 OVERESTIMATED assignment ambiguity; the correction strengthens
  the diagnostic's conclusion (midpoint assignment retained). No primary result changes:
  block_id was always assigned by window midpoint, which is unaffected.
