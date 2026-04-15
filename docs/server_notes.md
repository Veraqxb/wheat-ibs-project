# Server Notes

## Recommended Server Root

Example:

`/data1/xuebing/Z25_bam`

## Recommended Pipeline Folder

`/data1/xuebing/Z25_bam/05_IBS_pipeline`

Suggested layout:

- `scripts/`
- `configs/`
- `maps/`
- `results/`

## Practical Execution Pattern

For the current server project, chromosome VCFs already exist in:

- `*_all/VCF`

Therefore the preferred run mode is:

- start from chromosome-level VCFs
- concatenate VCFs
- filter variants
- calculate IBS
- run report generation

## Important Implementation Notes

- collect only `chr*.vcf.gz`
- avoid including old merged VCF files in the concat step
- use `bcftools concat -f vcf.list`
- treat `.csi` and `.tbi` as optional index types
- use `Rscript` rather than `bash` for `.R` files

## Map File Placement

Recommended location:

`/data1/xuebing/Z25_bam/05_IBS_pipeline/maps/`

Map file requirements:

- tab-delimited
- first row is header
- group names must exactly match config values
- sample IDs must exactly match PLINK `.mibs.id` sample names

## Current Known Practical Issues

- `GROUP_X` / `GROUP_Y` must match actual map headers exactly
- `plink --distance ibs square` outputs a full square matrix, not lower triangle
- heatmaps with hundreds of samples should be treated as overview figures, not sole evidence

