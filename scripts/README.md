# CAMP IBS Pipeline Scripts

This folder is organized around one active modular workflow plus a few focused plotting/reporting helpers.

## Active Modular Workflow

Run through `wheat_ibs_modular_pipeline.sh`.

- `step00_vcf_to_ibs.sh`: optional VCF to PLINK IBS generation.
- `step01_prepare_matrix.R`: read `.mibs/.mibs.id`, keep exact map-ID matches, and write a symmetric IBS matrix.
- `step02_build_reference_cluster.R`: build strict DNA duplicate clusters from the 2group reference, separately for the anchor group and B25.
- `step02b_build_reference_context.R`: build the 2group Z23/B25 bidirectional match context used by downstream RNA rescue.
- `step03_pairwise_matching.R`: match every query group against the anchor group in map order.
- `step04_diagnose_low_ibs.R`: draw low-IBS cross heatmaps for samples below threshold.
- `step05_cluster_rescue_and_summary.R`: combine pairwise matching, 2group cluster reference, and 2group match context into final sample decisions.

## Focused 2group DNA QC

- `two_group_internal_heatmaps.R`: publication/checking heatmaps for Z23, B25, expected Z23-B25 pairs, and full B25 x Z23 cross matrices.
- `two_group_ibs_density_thresholds.R`: expected-pair versus random/non-pair IBS density, ECDF, boxplot, and q99 background checks.
- `run_2cluster.R`: batch runner for 2group heatmaps and density plots.
- `dna_bidirectional_fine_classification.R`: standalone detailed 2group bidirectional classification from heatmap outputs.
- `dna_bidirectional_match_scatter.R`: scatter/barplot summary from curated Excel diagnosis sheets.

## Reporting Helpers

- `dna_diagnosis_summary.R`: summary plots/tables from diagnosis spreadsheets.
- `dna_ibs_density_from_xlsx.R`: spreadsheet-based IBS density summaries.
- `plot_group_heatmaps.R`: standalone heatmap export from existing IBS files.
- `samples.R`: older wrapper for DNA/RNA joint diagnosis; kept for compatibility while the modular steps remain the primary workflow.
- `audit_2group_bidirectional_rules.R`: rule-audit helper for checking historical classifications.

## Decision Principle

The 2group reference is now represented by two separate evidence layers:

- `cluster_table`: high-similarity DNA neighborhoods defined by IBS >= 0.99.
- `reference_context`: row-wise Z23/B25 bidirectional match status and confidence.

The 5group/RNA workflow uses Z23 as the primary identity anchor, then uses B25 and clusters only as rescue/explanation evidence. B25 rescue is flagged for manual review if the corresponding 2group reference context is `True_mismatch` or `No_data`.
