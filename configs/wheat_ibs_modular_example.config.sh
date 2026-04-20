MODE="ibs"

WORK_ROOT="/data1/xuebing/Z25_bam/05_IBS_pipeline/results/C2_modular"
PREFIX="C2_modular"
MAP_FILE="/data1/xuebing/Z25_bam/05_IBS_pipeline/maps/c2_id_map.txt"

# MODE=ibs
IBS_FILE="/data1/xuebing/Z25_bam/05_IBS_pipeline/results/C2_5group/04_plink/C2_5groups.final_qc.mibs"
IBS_ID_FILE="/data1/xuebing/Z25_bam/05_IBS_pipeline/results/C2_5group/04_plink/C2_5groups.final_qc.mibs.id"

# MODE=vcf
PROJECT_ROOT="/data1/xuebing/Z25_bam"
VCF_SOURCE_DIR="/data1/xuebing/Z25_bam/C2_5group/VCF"
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
MAX_GENO=0.4
HIGH_HET_THRESHOLD=0.05
COPY_MODE="link"

ANCHOR_COL="Z23"
SECONDARY_COL="B25"
RNA_GROUPS="TC SC FC"
ALL_GROUPS="Z23 B25 TC SC FC"
DNA_THRESHOLD=0.99
RNA_THRESHOLD=0.90
IBS_ZMIN=0.7
IBS_ZMAX=1.0
