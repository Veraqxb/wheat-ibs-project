#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash scripts/write_server_modular_configs.sh [target_config_dir]

Default target_config_dir:
  /data1/xuebing/Z25_bam/05_IBS/configs

This script writes the three same-ploidy modular config files:
  C2_modular.config.sh
  C4_modular.config.sh
  C6_modular.config.sh
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

TARGET_DIR="${1:-/data1/xuebing/Z25_bam/05_IBS/configs}"
mkdir -p "$TARGET_DIR"

cat > "${TARGET_DIR}/C2_modular.config.sh" <<'EOF'
MODE="ibs"

WORK_ROOT="/data1/xuebing/Z25_bam/05_IBS/results/C2_modular"
PREFIX="C2_modular"
MAP_FILE="/data1/xuebing/Z25_bam/05_IBS/maps/c2_id_map.txt"

ANCHOR_COL="Z23"
SECONDARY_COL="B25"
RNA_GROUPS="TC SC FC"
DNA_THRESHOLD=0.99
RNA_THRESHOLD=0.90
IBS_ZMIN=0.7
IBS_ZMAX=1.0

REFERENCE_PREFIX="C2_2group_reference"
QUERY_PREFIX="C2_5group_query"
REFERENCE_ALL_GROUPS="Z23 B25"
QUERY_ALL_GROUPS="Z23 B25 TC SC FC"
REFERENCE_MAX_GENO=0.3
QUERY_MAX_GENO=0.4

# MODE=ibs inputs
REFERENCE_IBS_FILE="/data1/xuebing/Z25_bam/05_IBS/ibs_file/C2_2groups_renamed.mibs"
REFERENCE_IBS_ID_FILE="/data1/xuebing/Z25_bam/05_IBS/ibs_file/C2_2groups_renamed.mibs.id"
QUERY_IBS_FILE="/data1/xuebing/Z25_bam/05_IBS/ibs_file/C2_5groups.final_qc.mibs"
QUERY_IBS_ID_FILE="/data1/xuebing/Z25_bam/05_IBS/ibs_file/C2_5groups.final_qc.mibs.id"

# MODE=vcf inputs
PROJECT_ROOT="/data1/xuebing/Z25_bam"
REFERENCE_VCF_SOURCE_DIR="/data1/xuebing/Z25_bam/C2_all/VCF"
QUERY_VCF_SOURCE_DIR="/data1/xuebing/Z25_bam/C2_5group/VCF"
BCFTOOLS_BIN="${PROJECT_ROOT}/bcftools"
BGZIP_BIN="${PROJECT_ROOT}/bgzip"
TABIX_BIN="${PROJECT_ROOT}/tabix"
PLINK_BIN="${PROJECT_ROOT}/plink"
PARALLEL_BIN="${PROJECT_ROOT}/parallel"
RSCRIPT_BIN="${PROJECT_ROOT}/Rscript"
THREADS_PARALLEL=8
CHR_SET=42
CHR_LIST="005 006 011 012 017 018 023 024 029 030 035 036 041 042"
MIN_MAC=2
HIGH_HET_THRESHOLD=0.05
COPY_MODE="link"
EOF

cat > "${TARGET_DIR}/C4_modular.config.sh" <<'EOF'
MODE="ibs"

WORK_ROOT="/data1/xuebing/Z25_bam/05_IBS/results/C4_modular"
PREFIX="C4_modular"
MAP_FILE="/data1/xuebing/Z25_bam/05_IBS/maps/c4_id_map.txt"

ANCHOR_COL="Z23"
SECONDARY_COL="B25"
RNA_GROUPS="TC SC FC"
DNA_THRESHOLD=0.99
RNA_THRESHOLD=0.90
IBS_ZMIN=0.7
IBS_ZMAX=1.0

REFERENCE_PREFIX="C4_2group_reference"
QUERY_PREFIX="C4_5group_query"
REFERENCE_ALL_GROUPS="Z23 B25"
QUERY_ALL_GROUPS="Z23 B25 TC SC FC"
REFERENCE_MAX_GENO=0.3
QUERY_MAX_GENO=0.4

# MODE=ibs inputs
REFERENCE_IBS_FILE="/data1/xuebing/Z25_bam/05_IBS/ibs_file/C4_2groups_renamed.mibs"
REFERENCE_IBS_ID_FILE="/data1/xuebing/Z25_bam/05_IBS/ibs_file/C4_2groups_renamed.mibs.id"
QUERY_IBS_FILE="/data1/xuebing/Z25_bam/05_IBS/ibs_file/C4_5group_final_qc.mibs"
QUERY_IBS_ID_FILE="/data1/xuebing/Z25_bam/05_IBS/ibs_file/C4_5group_final_qc.mibs.id"

# MODE=vcf inputs
PROJECT_ROOT="/data1/xuebing/Z25_bam"
REFERENCE_VCF_SOURCE_DIR="/data1/xuebing/Z25_bam/C4_all/VCF"
QUERY_VCF_SOURCE_DIR="/data1/xuebing/Z25_bam/C4_5group/VCF"
BCFTOOLS_BIN="${PROJECT_ROOT}/bcftools"
BGZIP_BIN="${PROJECT_ROOT}/bgzip"
TABIX_BIN="${PROJECT_ROOT}/tabix"
PLINK_BIN="${PROJECT_ROOT}/plink"
PARALLEL_BIN="${PROJECT_ROOT}/parallel"
RSCRIPT_BIN="${PROJECT_ROOT}/Rscript"
THREADS_PARALLEL=8
CHR_SET=42
CHR_LIST="001 002 003 004 007 008 009 010 013 014 015 016 019 020 021 022 025 026 027 028 031 032 033 034 037 038 039 040"
MIN_MAC=2
HIGH_HET_THRESHOLD=0.05
COPY_MODE="link"
EOF

cat > "${TARGET_DIR}/C6_modular.config.sh" <<'EOF'
MODE="ibs"

WORK_ROOT="/data1/xuebing/Z25_bam/05_IBS/results/C6_modular"
PREFIX="C6_modular"
MAP_FILE="/data1/xuebing/Z25_bam/05_IBS/maps/c6_id_map.txt"

ANCHOR_COL="Z23"
SECONDARY_COL="B25"
RNA_GROUPS="TC SC FC"
DNA_THRESHOLD=0.99
RNA_THRESHOLD=0.90
IBS_ZMIN=0.7
IBS_ZMAX=1.0

REFERENCE_PREFIX="C6_2group_reference"
QUERY_PREFIX="C6_5group_query"
REFERENCE_ALL_GROUPS="Z23 B25"
QUERY_ALL_GROUPS="Z23 B25 TC SC FC"
REFERENCE_MAX_GENO=0.3
QUERY_MAX_GENO=0.4

# MODE=ibs inputs
REFERENCE_IBS_FILE="/data1/xuebing/Z25_bam/05_IBS/ibs_file/C6_2groups_final_qc.mibs"
REFERENCE_IBS_ID_FILE="/data1/xuebing/Z25_bam/05_IBS/ibs_file/C6_2groups_final_qc.mibs.id"
QUERY_IBS_FILE="/data1/xuebing/Z25_bam/05_IBS/ibs_file/C6_5group_final.mibs"
QUERY_IBS_ID_FILE="/data1/xuebing/Z25_bam/05_IBS/ibs_file/C6_5group_final.mibs.id"

# MODE=vcf inputs
PROJECT_ROOT="/data1/xuebing/Z25_bam"
REFERENCE_VCF_SOURCE_DIR="/data1/xuebing/Z25_bam/C6_all/VCF"
QUERY_VCF_SOURCE_DIR="/data1/xuebing/Z25_bam/C6_5group/VCF"
BCFTOOLS_BIN="${PROJECT_ROOT}/bcftools"
BGZIP_BIN="${PROJECT_ROOT}/bgzip"
TABIX_BIN="${PROJECT_ROOT}/tabix"
PLINK_BIN="${PROJECT_ROOT}/plink"
PARALLEL_BIN="${PROJECT_ROOT}/parallel"
RSCRIPT_BIN="${PROJECT_ROOT}/Rscript"
THREADS_PARALLEL=8
CHR_SET=42
CHR_LIST="001 002 003 004 005 006 007 008 009 010 011 012 013 014 015 016 017 018 019 020 021 022 023 024 025 026 027 028 029 030 031 032 033 034 035 036 037 038 039 040 041 042"
MIN_MAC=2
HIGH_HET_THRESHOLD=0.05
COPY_MODE="link"
EOF

echo "Wrote:"
printf '  %s\n' \
  "${TARGET_DIR}/C2_modular.config.sh" \
  "${TARGET_DIR}/C4_modular.config.sh" \
  "${TARGET_DIR}/C6_modular.config.sh"
