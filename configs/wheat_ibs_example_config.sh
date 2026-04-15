FASTCALL_JAR="/path/to/TIGER_F3_20251121.jar"
REFERENCE_DIR="/data/public_data/wheat/00_genome/reference/v1.0/byChr"
BAM_INFO="/path/to/C4_bam_info_206.txt"
LIB_DIR="/path/to/v4_random_50k_lib"
CHR_LIST="001 002 003 004 007 008 009 010 013 014 015 016 019 020 021 022 025 026 027 028 031 032 033 034 037 038 039 040"

THREADS_FASTCALL=36
THREADS_PARALLEL=8
JAVA_XMX="100g"

SAMTOOLS_BIN="/path/to/samtools"
BCFTOOLS_BIN="/path/to/bcftools"
BGZIP_BIN="/path/to/bgzip"
TABIX_BIN="/path/to/tabix"
PLINK_BIN="/path/to/plink"
PARALLEL_BIN="/path/to/parallel"
RSCRIPT_BIN="/usr/bin/Rscript"

WORK_ROOT="/path/to/project_run"
PREFIX="C4_2groups"
CHR_SET=42

MAP_FILE="/path/to/C4_map_2group.txt"
GROUP_Y="B25"
GROUP_X="23"

FASTCALL_F=30
FASTCALL_G=20
FASTCALL_H=0.05
FASTCALL_E=0
MIN_MAC=2
MAX_GENO=0.2
HIGH_HET_THRESHOLD=0.05
IBS_ZMIN=0.7
IBS_ZMAX=1.0
