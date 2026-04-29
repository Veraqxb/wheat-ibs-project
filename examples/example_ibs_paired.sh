#!/usr/bin/env bash
set -euo pipefail

# Run the packaged CAMP IBS workflow from existing .mibs/.mibs.id files.
# Execute from the project root that contains ibs_file/ and maps/.

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
  --rna-groups "TC SC FC" \
  --dna-threshold 0.99 \
  --rna-threshold 0.90
