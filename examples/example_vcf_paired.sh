#!/usr/bin/env bash
set -euo pipefail

# Run the packaged CAMP IBS workflow from paired 2group/5group VCF directories.
# Update tool paths and VCF directories before running.

wheat_ibs_project/bin/camp-ibs \
  --mode vcf \
  --work-root results/C2_from_vcf \
  --prefix C2_from_vcf \
  --map maps/c2_id_map.txt \
  --project-root /data1/xuebing/Z25_bam \
  --reference-vcf-dir /data1/xuebing/Z25_bam/C2_all/VCF \
  --query-vcf-dir /data1/xuebing/Z25_bam/C2_5group/VCF \
  --bcftools /data1/xuebing/Z25_bam/bcftools \
  --bgzip /data1/xuebing/Z25_bam/bgzip \
  --tabix /data1/xuebing/Z25_bam/tabix \
  --plink /data1/xuebing/Z25_bam/plink \
  --parallel /data1/xuebing/Z25_bam/parallel \
  --threads 8 \
  --reference-max-geno 0.3 \
  --query-max-geno 0.4
