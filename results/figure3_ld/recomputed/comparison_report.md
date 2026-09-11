# Comparison: recomputed vs frozen LD-block tables

## Recorder bit-identity gate (Amendment 1)
- recorder_bit_identity: additional `run_block_null_rec` pass vs primary `run_block_null` (same seed, spawn_key (2,)) — max|Δsizes| 0, max|Δmeans| 0.000e+00, max|Δpcts| 0.000e+00 -> PASS (bit-identical; recorded gene_ids are the primary null's)

Frozen inputs verified against `data/figure3_ld/frozen/SHA256SUMS.txt` before comparison.
- frozen SHA256 integrity: PASS

## Table_S_block_vs_gene_bootstrap.tsv
  - PASS obs_pct: max|diff| 0.000e+00 within storage tolerance
  - PASS gw_resid: max|diff| 0.000e+00 within storage tolerance
  - PASS gw_lo: max|diff| 0.000e+00 within storage tolerance
  - PASS gw_hi: max|diff| 0.000e+00 within storage tolerance
  - PASS gw_p: max|diff| 0.000e+00 within storage tolerance
  - PASS gw_q_BH5: max|diff| nan within storage tolerance
  - PASS gw_q_BH8: max|diff| 0.000e+00 within storage tolerance
  - PASS gw_p_ctr: max|diff| 0.000e+00 within storage tolerance
  - PASS gw_ctr_q_BH5: max|diff| nan within storage tolerance
  - PASS gw_ctr_q_BH8: max|diff| 0.000e+00 within storage tolerance
  - PASS bw_resid: max|diff| 0.000e+00 within storage tolerance
  - PASS bw_lo: max|diff| 0.000e+00 within storage tolerance
  - PASS bw_hi: max|diff| 0.000e+00 within storage tolerance
  - PASS bw_p: max|diff| 0.000e+00 within storage tolerance
  - PASS bw_q_BH5: max|diff| nan within storage tolerance
  - PASS bw_q_BH8: max|diff| 0.000e+00 within storage tolerance
  - PASS gw_nullsd: max|diff| 0.000e+00 within storage tolerance
  - PASS bw_nullsd: max|diff| 0.000e+00 within storage tolerance
  - PASS sd_inflation: max|diff| 0.000e+00 within storage tolerance
  - PASS n_L2G: identical

## Table_final_five_axis_sensitivity.tsv
  - PASS observed_pct: max|diff| 0.000e+00 within storage tolerance
  - PASS null_mean: max|diff| 0.000e+00 within storage tolerance
  - PASS null_sd: max|diff| 0.000e+00 within storage tolerance
  - PASS resid_pp: max|diff| 0.000e+00 within storage tolerance
  - PASS ref95_lo: max|diff| 0.000e+00 within storage tolerance
  - PASS ref95_hi: max|diff| 0.000e+00 within storage tolerance
  - PASS p_two: max|diff| 0.000e+00 within storage tolerance
  - PASS p_mc_lo: max|diff| 0.000e+00 within storage tolerance
  - PASS p_mc_hi: max|diff| 0.000e+00 within storage tolerance
  - PASS q_BH5: max|diff| 0.000e+00 within storage tolerance
  - PASS q_BH5_mc_lo: max|diff| 0.000e+00 within storage tolerance
  - PASS q_BH5_mc_hi: max|diff| 0.000e+00 within storage tolerance
  - PASS B: identical
  - PASS b_le: identical
  - PASS b_ge: identical

## Table_final_five_axis_exact_gps_minus_d.tsv
  - PASS gps_univ_mean_before: max|diff| 0.000e+00 within storage tolerance
  - PASS gps_univ_mean_after: max|diff| 0.000e+00 within storage tolerance
  - PASS gps_univ_sd_before: max|diff| 0.000e+00 within storage tolerance
  - PASS gps_univ_sd_after: max|diff| 0.000e+00 within storage tolerance
  - PASS gps_univ_median_before: max|diff| 0.000e+00 within storage tolerance
  - PASS gps_univ_median_after: max|diff| 0.000e+00 within storage tolerance
  - PASS observed_pct: max|diff| 0.000e+00 within storage tolerance
  - PASS null_mean: max|diff| 0.000e+00 within storage tolerance
  - PASS null_sd: max|diff| 0.000e+00 within storage tolerance
  - PASS resid_pp: max|diff| 0.000e+00 within storage tolerance
  - PASS ref95_lo: max|diff| 0.000e+00 within storage tolerance
  - PASS ref95_hi: max|diff| 0.000e+00 within storage tolerance
  - PASS p_two: max|diff| 0.000e+00 within storage tolerance
  - PASS p_mc_lo: max|diff| 0.000e+00 within storage tolerance
  - PASS p_mc_hi: max|diff| 0.000e+00 within storage tolerance
  - PASS q_BH5: max|diff| 0.000e+00 within storage tolerance
  - PASS q_BH5_mc_lo: max|diff| 0.000e+00 within storage tolerance
  - PASS q_BH5_mc_hi: max|diff| 0.000e+00 within storage tolerance
  - PASS occ_mean: max|diff| 0.000e+00 within storage tolerance
  - PASS occ_lo: max|diff| 0.000e+00 within storage tolerance
  - PASS occ_hi: max|diff| 0.000e+00 within storage tolerance
  - PASS match_maxabs_dmean: max|diff| 0.000e+00 within storage tolerance
  - PASS match_pct_within_005: max|diff| 0.000e+00 within storage tolerance
  - PASS borrowed_mean: max|diff| 0.000e+00 within storage tolerance
  - PASS n_genes_changed_universe_7707: identical
  - PASS n_genes_changed_risk_1016: identical
  - PASS n_genes_changed_pool_3605: identical
  - PASS n_exact_beyond_approx_poolrisk: identical
  - PASS n_exact_beyond_approx_risk: identical
  - PASS gps_univ_max_before: identical
  - PASS gps_univ_max_after: identical
  - PASS b_le: identical
  - PASS b_ge: identical
  - PASS n_bins: identical
  - PASS redraws: identical

## Figure3B_main_source_data.tsv
  - PASS observed_pct: max|diff| 0.000e+00 within storage tolerance
  - PASS null_mean_pct: max|diff| 0.000e+00 within storage tolerance
  - PASS null_sd: max|diff| 0.000e+00 within storage tolerance
  - PASS null_band_lo_q0.025: max|diff| 0.000e+00 within storage tolerance
  - PASS null_band_hi_q0.975: max|diff| 0.000e+00 within storage tolerance
  - PASS resid_pp: max|diff| 0.000e+00 within storage tolerance
  - PASS p_two: max|diff| 0.000e+00 within storage tolerance
  - PASS q_BH5: max|diff| 0.000e+00 within storage tolerance

## null_axis_percentages_block_gps_matched.tsv
- max|diff| after %.6f rounding: 0.000e+00 -> PASS (bit-identical at stored precision)

## Table_realised_genelevel_gps.tsv
- byte-identity sha256: recomputed `b968d9f449b9b1f71afa5e4cfee898743b1b92137f2035a3a8848bcc0b7d361c` vs frozen `b968d9f449b9b1f71afa5e4cfee898743b1b92137f2035a3a8848bcc0b7d361c` -> PASS
- per-cell tally: 429/429 cells identical (39 rows x 11 columns)

## Gene sets
- universe n=7707; risk n=1016; pool n=3,605 (asserted during setup; HALT on mismatch)
- axis membership sets asserted identical to the shipped-script sets for the 5 primary axes (HALT on mismatch)
- observed percentages: T2D 13.3858, CAD 2.8543, Hypertension 5.4134, Asthma 3.7402, Prostate Ca 3.0512, Atrial Fib 4.7244, Gout 3.4449, OA 6.8898

## Conclusion
All recomputed values match the frozen tables within storage precision.
## Historical gene-wise TSV (informational, NOT a gate)
The repository's committed `results/figure3_bootstrap_stats.tsv` is the submission-time table. Differences vs the final re-executed gene-wise null (this run, matching the frozen sensitivity table exactly):
- T2D: bg_mean 11.7336 -> 11.7338; CI lo 0.0984 -> 0.1969; p_two_sided 0.032 (centered convention) vs 0.037962 (R2#7)
- CAD: bg_mean 4.9469 -> 4.9469; CI lo -3.0536 -> -3.0536; p_two_sided 0.001 (centered convention) vs 0.001998 (R2#7)
- Hypertension: bg_mean 4.8043 -> 4.8043; CI lo -0.3937 -> -0.3937; p_two_sided 0.2478 (centered convention) vs 0.281718 (R2#7)
- Asthma: bg_mean 5.2933 -> 5.2933; CI lo -2.5591 -> -2.5591; p_two_sided 0.001 (centered convention) vs 0.001998 (R2#7)
- Prostate Ca: bg_mean 6.0348 -> 6.0386; CI lo -4.2323 -> -4.2323; p_two_sided 0.001 (centered convention) vs 0.001998 (R2#7)
The frozen revision tables are authoritative; the historical TSV is kept unchanged as the record of the submitted analysis.