#!/usr/bin/env bash
set -euo pipefail

# Run C2/C4/C6 from a single manifest table.
# Execute from the project root that contains wheat_ibs_project/.
# First edit wheat_ibs_project/examples/manifest_ibs.tsv so paths point to
# your maps/ and ibs_file/ directories.

wheat_ibs_project/bin/camp-ibs run-manifest \
  --manifest wheat_ibs_project/examples/manifest_ibs.tsv \
  --work-root results/manifest_run \
  --anchor-col Z23 \
  --secondary-col B25 \
  --rna-groups "TC SC FC" \
  --dna-threshold 0.99 \
  --rna-threshold 0.90
