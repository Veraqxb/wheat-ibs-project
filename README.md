# Wheat IBS Project Record

This folder keeps a local record of the wheat IBS identification project, including:

- workflow notes
- decision logic for sample filtering and retention
- code snapshots used during the project
- config template examples

## Structure

- `docs/project_log.md`
  Project summary, workflow, threshold logic, and sample decision rules.
- `docs/server_notes.md`
  Recommended server-side directory layout and execution notes.
- `scripts/`
  Code snapshots for the pipeline and report scripts.
- `configs/`
  Config template snapshot.

## New Modular Pipeline

The repository now includes a modular IBS matching workflow with two entry modes:

- `MODE=vcf`
  Read chromosome-level VCF files, filter, compute PLINK IBS, and export `.mibs/.mibs.id`
- `MODE=ibs`
  Reuse existing `.mibs/.mibs.id` and run downstream matching directly

Main entry:

- `bin/camp-ibs`
- `scripts/wheat_ibs_modular_pipeline.sh`
- `scripts/plot_group_heatmaps.R`

Steps:

- `scripts/two_group_internal_heatmaps.R`
- `scripts/run_2cluster.R`
- `scripts/two_group_ibs_density_thresholds.R`
- `scripts/samples.R`
- `scripts/dna_diagnosis_summary.R`
- `scripts/step00_vcf_to_ibs.sh`
- `scripts/step01_prepare_matrix.R`
- `scripts/step02_build_reference_cluster.R`
- `scripts/step02b_build_reference_context.R`
- `scripts/step03_pairwise_matching.R`
- `scripts/step04_diagnose_low_ibs.R`
- `scripts/step05_cluster_rescue_and_summary.R`

Project-level script roles:

- `two_group_internal_heatmaps.R` performs focused `2group` DNA reference QC, including Z23/B25 duplicate clusters, assignment tables, expected-pair diagonal heatmap, full cross-group heatmap, and structure heatmaps.
- `run_2cluster.R` batch-runs `two_group_internal_heatmaps.R` and the 2group density-threshold analysis from either a TSV config or automatic `ibs_file/` discovery.
- `two_group_ibs_density_thresholds.R` compares Z23-within, B25-within, full between-group, and expected-pair IBS distributions, then writes density/ECDF/boxplot outputs plus a threshold recommendation table.
- `samples.R` is the DNA-to-RNA joint diagnosis wrapper around the modular pipeline.
- `dna_diagnosis_summary.R` summarizes TSV diagnosis outputs.

Core design:

- The first column of the sample map is used as the anchor reference group
- Anchor-group clusters are defined by `IBS >= 0.99`
- The same-ploidy `2group` dataset can be used as the reference cluster library for the downstream `5group` dataset
- In `2group` reference mode, both `Z23` and `B25` clusters are built, their row-wise overlap is exported, and the combined cluster table is used for downstream rescue
- The `2group` reference also exports a row-wise bidirectional Z23/B25 match context table, so downstream RNA rescue can know whether the corresponding B25 reference is trusted, review-only, or excluded
- Every downstream group is matched against the anchor in map order
- Sample names are matched by exact string equality only; no fuzzy ID rescue is used when reading the map against `.mibs.id`
- RNA groups use Z23-first matching and B25 as rescue
- B25 rescue never silently overrides Z23. If the 2group reference context is `True_mismatch` or `No_data`, B25 rescue is written as `REVIEW_B25_CONTEXT_RISK` instead of an automatic rescued sample.
- Low-IBS / swapped samples are checked against anchor high-similarity clusters
- `DETECTION_MODE=simple` can be used to disable B25 direct rescue and scan later groups only against the first-column cluster reference
- RNA scan summaries include the within-group top 1% IBS background value and its difference from the 0.90 threshold, written to `*_threshold_check.log`

Formal same-ploidy configs:

- `configs/C2_modular.config.sh`
- `configs/C4_modular.config.sh`
- `configs/C6_modular.config.sh`

These configs contain both:

- `2group` reference inputs for duplicate-cluster library building
- `5group` query inputs for downstream sample identification

So in practice you only need to run three config files, one per ploidy.

Recommended usage:

Packaged command-line wrapper from existing `.mibs/.mibs.id` files:

```bash
wheat_ibs_project/bin/camp-ibs \
  --mode ibs \
  --work-root results/C2_packaged \
  --prefix C2_packaged \
  --map maps/c2_id_map.txt \
  --reference-mibs ibs_file/C2_2groups_renamed.mibs \
  --reference-id ibs_file/C2_2groups_renamed.mibs.id \
  --query-mibs ibs_file/C2_5groups.final_qc.mibs \
  --query-id ibs_file/C2_5groups.final_qc.mibs.id \
  --anchor-col Z23 \
  --secondary-col B25 \
  --rna-groups "TC SC FC"
```

Use `--dry-run` to write the generated config without executing the pipeline.

Config-file mode:

```bash
bash scripts/wheat_ibs_modular_pipeline.sh configs/C2_modular.config.sh
bash scripts/wheat_ibs_modular_pipeline.sh configs/C4_modular.config.sh
bash scripts/wheat_ibs_modular_pipeline.sh configs/C6_modular.config.sh
```

Packaging design:

- `docs/package_design.md`
- `examples/example_ibs_paired.sh`
- `examples/example_vcf_paired.sh`

Focused 2group reference QC:

```bash
Rscript scripts/run_2cluster.R \
  --ibs-dir ibs_file \
  --map-dir maps \
  --info-file camp_info.txt \
  --outdir all_2group_internal_heatmaps
```

This single command generates Z23/B25 internal cluster heatmaps, expected Z23-B25 heatmaps, duplicate-cluster tables, and 2group IBS density-threshold plots.

To write the three same-ploidy config files directly on the server so they stay identical to the repository templates:

```bash
bash scripts/write_server_modular_configs.sh
```

Standalone heatmap export from existing IBS files:

```bash
Rscript scripts/plot_group_heatmaps.R \
  --reference-mibs /path/to/2group.mibs \
  --reference-id /path/to/2group.mibs.id \
  --query-mibs /path/to/5group.mibs \
  --query-id /path/to/5group.mibs.id \
  --map /path/to/id_map.txt \
  --outdir /path/to/heatmaps \
  --prefix C2_heatmaps \
  --anchor-col Z23 \
  --secondary-col B25 \
  --rna-groups "TC FC SC"
```

This script exports internal group heatmaps and ordered cross-group heatmaps. In cross-group plots, the anchor/reference group is placed on the x-axis, the query group is placed on the y-axis, and expected one-to-one pairs are outlined on the diagonal.

## Current Project Scope

The project focuses on:

- VCF filtering from chromosome-level VCF files
- IBS calculation with PLINK
- within-group repeatability checks
- between-group identity checks
- Z23-based decision logic for transcriptome samples

## Main Datasets

- Random-site resequencing VCFs for `Z23` and `B25`
- Exon / genotype / transcriptome VCFs for `Z23`, `B25`, `TC`, `SC`, `FC`

## Practical Analysis Principle

- Use random-site DNA VCF as the primary identity evidence
- Use exon / RNA-based VCF as supporting evidence
- Use `best match + second best + margin` rather than a single IBS value
- Use `Z23` as the first identity baseline for transcriptome sample decisions
