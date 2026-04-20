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

- `scripts/wheat_ibs_modular_pipeline.sh`

Steps:

- `scripts/step00_vcf_to_ibs.sh`
- `scripts/step01_prepare_matrix.R`
- `scripts/step02_build_reference_cluster.R`
- `scripts/step03_pairwise_matching.R`
- `scripts/step04_diagnose_low_ibs.R`
- `scripts/step05_cluster_rescue_and_summary.R`

Core design:

- The first column of the sample map is used as the anchor reference group
- Anchor-group clusters are defined by `IBS > 0.99`
- Every downstream group is matched against the anchor in map order
- Sample names are matched by exact string equality only; no fuzzy ID rescue is used when reading the map against `.mibs.id`
- RNA groups use Z23-first matching and B25 as rescue
- Low-IBS / swapped samples are checked against anchor high-similarity clusters

Example configs:

- `configs/wheat_ibs_modular_ibs_example.config.sh`
  Directly read existing `.mibs/.mibs.id` and run the full matching workflow
- `configs/wheat_ibs_modular_vcf_example.config.sh`
  Start from chromosome-level VCF files, compute IBS, then run the full matching workflow

Recommended usage:

```bash
bash scripts/wheat_ibs_modular_pipeline.sh configs/wheat_ibs_modular_ibs_example.config.sh
```

```bash
bash scripts/wheat_ibs_modular_pipeline.sh configs/wheat_ibs_modular_vcf_example.config.sh
```

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
