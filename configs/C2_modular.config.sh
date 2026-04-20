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
