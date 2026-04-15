# Project Log

## Goal

Build a reusable workflow for wheat germplasm identity checking from VCF filtering to IBS calculation and final sample decisions.

## Analysis Layers

### 1. Random-site resequencing VCF

Use as the main dataset for identity determination because it is closer to genome-wide random fingerprinting and is less affected by expression bias.

Main tasks:

- evaluate within-group reproducibility across years
- evaluate whether mapped `Z23 <-> B25` pairs are consistent
- identify mislabeling, swaps, or duplicate mismatches

### 2. Exon / genotype / transcriptome VCF

Use as supporting evidence rather than the sole decision basis.

Main tasks:

- evaluate whether `TC`, `SC`, `FC` samples are consistent with expected DNA samples
- verify whether non-matching transcriptome samples align better with `B25`
- support decisions about follow-up resequencing

## Filtering Logic

### Core VCF processing

1. Collect only chromosome-level files: `chr*.vcf.gz`
2. Concatenate per-chromosome VCFs using `bcftools concat -f vcf.list`
3. Keep bi-allelic SNPs only
4. Keep variants with `MAC >= 2`
5. Convert to PLINK format
6. Test candidate `--geno` thresholds before fixing the final value
7. Remove high-heterozygosity SNPs carefully
8. Calculate IBS matrix with `plink --distance ibs square`

## Why `--geno` Must Be Tested

`B25` contains shallow sequencing samples. A strict `--geno` threshold can remove too many SNPs and make IBS unstable.

Recommended test values:

- `0.2`
- `0.3`
- `0.4`
- `0.5`

## Choosing the Final `--geno`

Do not choose the threshold that keeps the most SNPs by default.

Choose the most stringent threshold that still:

- retains enough SNPs for stable IBS
- preserves expected pair matching
- keeps `margin = best - second_best` reasonably large

Suggested ranking:

1. Exclude thresholds with too few remaining SNPs
2. Prefer thresholds with the highest number of correct expected matches
3. If ties remain, prefer the larger average margin
4. If still tied, prefer the stricter threshold

## IBS Decision Variables

For each sample, calculate:

- `ibs_expected`
- `ibs_best`
- `second_best`
- `margin`
- `status`

Definitions:

- `ibs_expected`: IBS to the expected paired sample in the map
- `ibs_best`: highest IBS to any candidate in the comparison group
- `second_best`: second highest IBS
- `margin`: `ibs_best - second_best`
- `status`: whether expected pair equals best pair

## Group-Level Logic

### Within-group checks

Use self-comparison matrices for:

- `Z23`
- `B25`
- `TC`
- `SC`
- `FC`

Purpose:

- confirm that repeated or related samples cluster correctly
- identify outlier samples within the same group

### Between-group checks

Use mapped pair comparisons such as:

- `B25 vs Z23`
- `TC vs Z23`
- `SC vs Z23`
- `FC vs Z23`

Purpose:

- verify expected pair matching
- identify better-than-expected matches to alternate groups

## Transcriptome Decision Logic

Use `Z23` as the primary reference for transcriptome samples.

For each `TC` / `SC` / `FC` sample:

1. Check whether it matches the expected `Z23` sample
2. If not, check whether it matches the corresponding `B25` sample
3. If it matches `B25`, keep it for later resequencing
4. If it matches neither `Z23` nor `B25`, mark it as a removal candidate

## Recommended Final Labels

- `KEEP_Z23`
- `KEEP_B25`
- `REVIEW`
- `REMOVE`

## Recommended Output Files

- `*_pair_summary.tsv`
- `*_summary.tsv`
- `*_IBS_matrix.tsv`
- `*_mismatch_samples.tsv`
- `*_keep_for_reseq.tsv`
- `*_remove_candidates.tsv`
- overview heatmap PDFs
- issue-only heatmap PDFs

## Visualization Strategy

Large full heatmaps should not be the only result.

Recommended figures:

1. overview heatmap
2. issue-only heatmap
3. expected-vs-best comparison plot
4. margin ranking plot
5. clustering or MDS plot based on `1 - IBS`

For heatmaps:

- use `RdYlBu`
- mark `IBS = 1` with strongest highlight
- optionally mark `IBS >= 0.99` with a weaker highlight

