# CAMP DNA/RNA Sample Identity QC Workflow

This project implements an IBS-based sample identity QC workflow for CAMP population data.

## Entry Modes

The workflow supports two reproducible entry modes.

- `VCF mode`: starts from chromosome-level VCF files, filters variants, computes PLINK IBS, and exports `.mibs/.mibs.id`.
- `IBS mode`: starts from existing `.mibs`, `.mibs.id`, and sample map files.

Both modes converge on the same downstream matching logic.

## Biological Design

- DNA reference groups are `Z23` and `B25`.
- RNA groups are `TC`, `FC`, and `SC`.
- The same-ploidy `2group` data are used first to establish reliable DNA duplicate clusters.
- The same-ploidy `5group` or `6group` data are then diagnosed against the DNA reference.

## Decision Rules

- DNA duplicate clusters are defined by threshold clusters using `IBS >= 0.99`.
- `NA` values are preserved in displayed heatmaps.
- `NA` values do not participate in threshold-based cluster decisions.
- Hierarchical clustering is only a structure visualization aid.
- Identity decisions use explicit IBS thresholds, map-defined expected pairs, and threshold-cluster rescue logic.
- RNA matching uses a looser threshold than DNA, but DNA duplicate clusters remain strict.

## Matrix Handling

The pipeline keeps two matrix concepts separate whenever `NA` values matter.

- Diagnosis matrix: preserves `NA` and is used for all threshold decisions.
- Plotting/order matrix: may use temporary imputation only to compute display order for structure plots.

## Main Script Roles

- `two_group_internal_heatmaps.R`: focused `2group` DNA reference QC.
- `run_2cluster.R`: batch runner for multiple `2group` inputs.
- `samples.R`: DNA-to-RNA joint diagnosis entry point using the modular workflow.
- `dna_diagnosis_summary.R`: compact summary of final diagnosis outputs.
- `wheat_ibs_modular_pipeline.sh`: canonical full pipeline entry.

## Outputs

The workflow prefers reproducible TSV and PDF outputs:

- duplicate pair table
- cluster table
- assignment table
- map-order reference table
- expected-pair IBS table
- full cross-group IBS table
- pairwise matching tables
- low-IBS sample table
- final matching table
- group-level matching summary
- threshold-cluster heatmaps
- hierarchical structure heatmaps
- ordered cross-group heatmaps
