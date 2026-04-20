#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash wheat_ibs_modular_pipeline.sh <config.sh>

Required config:
  MODE=vcf|ibs
  WORK_ROOT
  PREFIX
  MAP_FILE

For MODE=vcf also require:
  Either:
    VCF_SOURCE_DIR
  Or paired reference/query inputs:
    REFERENCE_VCF_SOURCE_DIR
    QUERY_VCF_SOURCE_DIR

For MODE=ibs also require:
  Either:
    IBS_FILE
    IBS_ID_FILE
  Or paired reference/query inputs:
    REFERENCE_IBS_FILE
    REFERENCE_IBS_ID_FILE
    QUERY_IBS_FILE
    QUERY_IBS_ID_FILE

Optional:
  RSCRIPT_BIN        default Rscript
  ANCHOR_COL         default first map column
  SECONDARY_COL      default B25
  RNA_GROUPS         default "TC SC FC"
  DNA_THRESHOLD      default 0.99
  RNA_THRESHOLD      default 0.90
  IBS_ZMIN           default 0.7
  IBS_ZMAX           default 1.0
EOF
}

require_var() {
  local var_name="$1"
  if [[ -z "${!var_name:-}" ]]; then
    echo "Missing required config variable: ${var_name}" >&2
    exit 1
  fi
}

[[ $# -eq 1 ]] || { usage; exit 1; }
CONFIG_FILE="$1"
[[ -f "$CONFIG_FILE" ]] || { echo "Config file not found: $CONFIG_FILE" >&2; exit 1; }
source "$CONFIG_FILE"

require_var MODE
require_var WORK_ROOT
require_var PREFIX
require_var MAP_FILE

RSCRIPT_BIN="${RSCRIPT_BIN:-Rscript}"
ANCHOR_COL="${ANCHOR_COL:-}"
SECONDARY_COL="${SECONDARY_COL:-B25}"
RNA_GROUPS="${RNA_GROUPS:-TC SC FC}"
DNA_THRESHOLD="${DNA_THRESHOLD:-0.99}"
RNA_THRESHOLD="${RNA_THRESHOLD:-0.90}"
IBS_ZMIN="${IBS_ZMIN:-0.7}"
IBS_ZMAX="${IBS_ZMAX:-1.0}"
REFERENCE_PREFIX="${REFERENCE_PREFIX:-${PREFIX}_2group}"
QUERY_PREFIX="${QUERY_PREFIX:-${PREFIX}_5group}"
REFERENCE_ALL_GROUPS="${REFERENCE_ALL_GROUPS:-Z23 B25}"
QUERY_ALL_GROUPS="${QUERY_ALL_GROUPS:-Z23 B25 ${RNA_GROUPS}}"
REFERENCE_MAX_GENO="${REFERENCE_MAX_GENO:-0.3}"
QUERY_MAX_GENO="${QUERY_MAX_GENO:-0.4}"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REF_ROOT="${WORK_ROOT}/reference_2group"
QUERY_ROOT="${WORK_ROOT}/query_5group"
STEP01_REF_DIR="${REF_ROOT}/step01_prepare_matrix"
STEP02_DIR="${REF_ROOT}/step02_build_reference_cluster"
STEP01_QUERY_DIR="${QUERY_ROOT}/step01_prepare_matrix"
STEP03_DIR="${QUERY_ROOT}/step03_pairwise_matching"
STEP04_DIR="${QUERY_ROOT}/step04_diagnose_low_ibs"
STEP05_DIR="${WORK_ROOT}/step05_cluster_rescue_and_summary"
TMP_DIR="${WORK_ROOT}/tmp"

mkdir -p "$STEP01_REF_DIR" "$STEP02_DIR" "$STEP01_QUERY_DIR" "$STEP03_DIR" "$STEP04_DIR" "$STEP05_DIR" "$TMP_DIR"

run_step00_with_temp_config() {
  local temp_config="$1"
  shift
  cat > "$temp_config" <<EOF
PROJECT_ROOT="${PROJECT_ROOT}"
VCF_SOURCE_DIR="$1"
WORK_ROOT="$2"
PREFIX="$3"
MAP_FILE="${MAP_FILE}"
BCFTOOLS_BIN="${BCFTOOLS_BIN}"
BGZIP_BIN="${BGZIP_BIN}"
TABIX_BIN="${TABIX_BIN}"
PLINK_BIN="${PLINK_BIN}"
PARALLEL_BIN="${PARALLEL_BIN}"
THREADS_PARALLEL="${THREADS_PARALLEL:-8}"
CHR_SET="${CHR_SET:-42}"
CHR_LIST="${CHR_LIST:-}"
MIN_MAC="${MIN_MAC:-2}"
MAX_GENO="$4"
HIGH_HET_THRESHOLD="${HIGH_HET_THRESHOLD:-0.05}"
COPY_MODE="${COPY_MODE:-link}"
ALL_GROUPS="$5"
EOF
  bash "${SCRIPT_DIR}/step00_vcf_to_ibs.sh" "$temp_config"
}

PAIRED_MODE=false
if [[ -n "${REFERENCE_IBS_FILE:-}" || -n "${REFERENCE_VCF_SOURCE_DIR:-}" ]]; then
  PAIRED_MODE=true
fi

if [[ "$MODE" == "vcf" ]]; then
  require_var PROJECT_ROOT
  require_var BCFTOOLS_BIN
  require_var BGZIP_BIN
  require_var TABIX_BIN
  require_var PLINK_BIN
  require_var PARALLEL_BIN
  if [[ "$PAIRED_MODE" == true ]]; then
    require_var REFERENCE_VCF_SOURCE_DIR
    require_var QUERY_VCF_SOURCE_DIR
    run_step00_with_temp_config "${TMP_DIR}/reference_step00.config.sh" "$REFERENCE_VCF_SOURCE_DIR" "$REF_ROOT" "$REFERENCE_PREFIX" "$REFERENCE_MAX_GENO" "$REFERENCE_ALL_GROUPS"
    run_step00_with_temp_config "${TMP_DIR}/query_step00.config.sh" "$QUERY_VCF_SOURCE_DIR" "$QUERY_ROOT" "$QUERY_PREFIX" "$QUERY_MAX_GENO" "$QUERY_ALL_GROUPS"
    REFERENCE_IBS_FILE="${REF_ROOT}/step00_vcf_to_ibs/04_plink/${REFERENCE_PREFIX}.final_qc.mibs"
    REFERENCE_IBS_ID_FILE="${REF_ROOT}/step00_vcf_to_ibs/04_plink/${REFERENCE_PREFIX}.final_qc.mibs.id"
    QUERY_IBS_FILE="${QUERY_ROOT}/step00_vcf_to_ibs/04_plink/${QUERY_PREFIX}.final_qc.mibs"
    QUERY_IBS_ID_FILE="${QUERY_ROOT}/step00_vcf_to_ibs/04_plink/${QUERY_PREFIX}.final_qc.mibs.id"
  else
    require_var VCF_SOURCE_DIR
    bash "${SCRIPT_DIR}/step00_vcf_to_ibs.sh" "$CONFIG_FILE"
    IBS_FILE="${WORK_ROOT}/step00_vcf_to_ibs/04_plink/${PREFIX}.final_qc.mibs"
    IBS_ID_FILE="${WORK_ROOT}/step00_vcf_to_ibs/04_plink/${PREFIX}.final_qc.mibs.id"
  fi
elif [[ "$MODE" == "ibs" ]]; then
  if [[ "$PAIRED_MODE" == true ]]; then
    require_var REFERENCE_IBS_FILE
    require_var REFERENCE_IBS_ID_FILE
    require_var QUERY_IBS_FILE
    require_var QUERY_IBS_ID_FILE
  else
    require_var IBS_FILE
    require_var IBS_ID_FILE
  fi
else
  echo "Unsupported MODE: ${MODE}" >&2
  exit 1
fi

if [[ "$PAIRED_MODE" == true ]]; then
  "$RSCRIPT_BIN" "${SCRIPT_DIR}/step01_prepare_matrix.R" \
    --mibs "$REFERENCE_IBS_FILE" \
    --id "$REFERENCE_IBS_ID_FILE" \
    --map "$MAP_FILE" \
    --outdir "$STEP01_REF_DIR" \
    --prefix "$REFERENCE_PREFIX"

  "$RSCRIPT_BIN" "${SCRIPT_DIR}/step02_build_reference_cluster.R" \
    --matrix "${STEP01_REF_DIR}/${REFERENCE_PREFIX}_ibs_matrix.tsv" \
    --map "$MAP_FILE" \
    --outdir "$STEP02_DIR" \
    --prefix "$REFERENCE_PREFIX" \
    --anchor-col "$ANCHOR_COL" \
    --zmin "$IBS_ZMIN" \
    --zmax "$IBS_ZMAX"

  "$RSCRIPT_BIN" "${SCRIPT_DIR}/step01_prepare_matrix.R" \
    --mibs "$QUERY_IBS_FILE" \
    --id "$QUERY_IBS_ID_FILE" \
    --map "$MAP_FILE" \
    --outdir "$STEP01_QUERY_DIR" \
    --prefix "$QUERY_PREFIX"

  "$RSCRIPT_BIN" "${SCRIPT_DIR}/step03_pairwise_matching.R" \
    --matrix "${STEP01_QUERY_DIR}/${QUERY_PREFIX}_ibs_matrix.tsv" \
    --map "$MAP_FILE" \
    --outdir "$STEP03_DIR" \
    --prefix "$QUERY_PREFIX" \
    --anchor-col "$ANCHOR_COL" \
    --rna-groups "$RNA_GROUPS" \
    --dna-threshold "$DNA_THRESHOLD" \
    --rna-threshold "$RNA_THRESHOLD" \
    --zmin "$IBS_ZMIN" \
    --zmax "$IBS_ZMAX"

  "$RSCRIPT_BIN" "${SCRIPT_DIR}/step04_diagnose_low_ibs.R" \
    --matrix "${STEP01_QUERY_DIR}/${QUERY_PREFIX}_ibs_matrix.tsv" \
    --map "$MAP_FILE" \
    --pairwise-dir "$STEP03_DIR" \
    --outdir "$STEP04_DIR" \
    --prefix "$QUERY_PREFIX" \
    --anchor-col "$ANCHOR_COL" \
    --zmin "$IBS_ZMIN" \
    --zmax "$IBS_ZMAX"

  "$RSCRIPT_BIN" "${SCRIPT_DIR}/step05_cluster_rescue_and_summary.R" \
    --matrix "${STEP01_QUERY_DIR}/${QUERY_PREFIX}_ibs_matrix.tsv" \
    --map "$MAP_FILE" \
    --cluster-table "${STEP02_DIR}/${REFERENCE_PREFIX}_cluster_table.tsv" \
    --pairwise-dir "$STEP03_DIR" \
    --outdir "$STEP05_DIR" \
    --prefix "$PREFIX" \
    --anchor-col "$ANCHOR_COL" \
    --secondary-col "$SECONDARY_COL" \
    --rna-groups "$RNA_GROUPS" \
    --dna-threshold "$DNA_THRESHOLD" \
    --rna-threshold "$RNA_THRESHOLD"
else
  "$RSCRIPT_BIN" "${SCRIPT_DIR}/step01_prepare_matrix.R" \
    --mibs "$IBS_FILE" \
    --id "$IBS_ID_FILE" \
    --map "$MAP_FILE" \
    --outdir "${WORK_ROOT}/step01_prepare_matrix" \
    --prefix "$PREFIX"

  "$RSCRIPT_BIN" "${SCRIPT_DIR}/step02_build_reference_cluster.R" \
    --matrix "${WORK_ROOT}/step01_prepare_matrix/${PREFIX}_ibs_matrix.tsv" \
    --map "$MAP_FILE" \
    --outdir "${WORK_ROOT}/step02_build_reference_cluster" \
    --prefix "$PREFIX" \
    --anchor-col "$ANCHOR_COL" \
    --zmin "$IBS_ZMIN" \
    --zmax "$IBS_ZMAX"

  "$RSCRIPT_BIN" "${SCRIPT_DIR}/step03_pairwise_matching.R" \
    --matrix "${WORK_ROOT}/step01_prepare_matrix/${PREFIX}_ibs_matrix.tsv" \
    --map "$MAP_FILE" \
    --outdir "${WORK_ROOT}/step03_pairwise_matching" \
    --prefix "$PREFIX" \
    --anchor-col "$ANCHOR_COL" \
    --rna-groups "$RNA_GROUPS" \
    --dna-threshold "$DNA_THRESHOLD" \
    --rna-threshold "$RNA_THRESHOLD" \
    --zmin "$IBS_ZMIN" \
    --zmax "$IBS_ZMAX"

  "$RSCRIPT_BIN" "${SCRIPT_DIR}/step04_diagnose_low_ibs.R" \
    --matrix "${WORK_ROOT}/step01_prepare_matrix/${PREFIX}_ibs_matrix.tsv" \
    --map "$MAP_FILE" \
    --pairwise-dir "${WORK_ROOT}/step03_pairwise_matching" \
    --outdir "${WORK_ROOT}/step04_diagnose_low_ibs" \
    --prefix "$PREFIX" \
    --anchor-col "$ANCHOR_COL" \
    --zmin "$IBS_ZMIN" \
    --zmax "$IBS_ZMAX"

  "$RSCRIPT_BIN" "${SCRIPT_DIR}/step05_cluster_rescue_and_summary.R" \
    --matrix "${WORK_ROOT}/step01_prepare_matrix/${PREFIX}_ibs_matrix.tsv" \
    --map "$MAP_FILE" \
    --cluster-table "${WORK_ROOT}/step02_build_reference_cluster/${PREFIX}_cluster_table.tsv" \
    --pairwise-dir "${WORK_ROOT}/step03_pairwise_matching" \
    --outdir "${WORK_ROOT}/step05_cluster_rescue_and_summary" \
    --prefix "$PREFIX" \
    --anchor-col "$ANCHOR_COL" \
    --secondary-col "$SECONDARY_COL" \
    --rna-groups "$RNA_GROUPS" \
    --dna-threshold "$DNA_THRESHOLD" \
    --rna-threshold "$RNA_THRESHOLD"
fi

echo "Modular pipeline completed for ${PREFIX}"
