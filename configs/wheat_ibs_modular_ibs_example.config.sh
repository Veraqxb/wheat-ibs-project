MODE="ibs"

# Output root for the full modular matching workflow
WORK_ROOT="/data1/xuebing/Z25_bam/05_IBS_pipeline/results/C2_modular_ibs"
PREFIX="C2_modular_ibs"
MAP_FILE="/data1/xuebing/Z25_bam/05_IBS_pipeline/maps/c2_id_map.txt"

# Existing IBS inputs
IBS_FILE="/data1/xuebing/Z25_bam/05_IBS_pipeline/results/C2_5group/04_plink/C2_5groups.final_qc.mibs"
IBS_ID_FILE="/data1/xuebing/Z25_bam/05_IBS_pipeline/results/C2_5group/04_plink/C2_5groups.final_qc.mibs.id"

# Downstream matching behavior
RSCRIPT_BIN="/data1/xuebing/Z25_bam/Rscript"
ANCHOR_COL="Z23"
SECONDARY_COL="B25"
RNA_GROUPS="TC SC FC"
DNA_THRESHOLD=0.99
RNA_THRESHOLD=0.90
IBS_ZMIN=0.7
IBS_ZMAX=1.0
