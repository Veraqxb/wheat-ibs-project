# TODO

## Workflow

- Add a server-ready `VCF-only` pipeline snapshot matching the latest debugging results
- Add `--geno` scan support for `0.2 / 0.3 / 0.4 / 0.5`
- Add automatic recommendation of the final `--geno` threshold
- Add resume logic so report generation can be rerun without repeating concat/filter/plink

## Reporting

- Add mismatch sample export
- Add keep/remove candidate export
- Add issue-only heatmap output
- Add expected-vs-best IBS comparison plot
- Add margin ranking plot
- Add clustering or MDS output based on `1 - IBS`

## Decision Support

- Separate default rules for DNA random-site data and exon/RNA data
- Add `KEEP_Z23 / KEEP_B25 / REVIEW / REMOVE` classification into report outputs
- Add summary table for transcriptome samples that fail Z23 matching but match B25

## Documentation

- Add server deployment steps as a copy-paste shell block
- Add sample map format examples for `Z23/B25` and `TC/SC/FC`
- Add notes on heatmap interpretation for large sample sets

