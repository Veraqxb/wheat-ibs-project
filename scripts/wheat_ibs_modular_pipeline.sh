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
  VCF_SOURCE_DIR and tool paths needed by step00_vcf_to_ibs.sh

For MODE=ibs also require:
  IBS_FILE
  IBS_ID_FILE

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

[[ $# -eq 1 ]] || { usage; exit 1; }
CONFIG_FILE="$1"
[[ -f "$CONFIG_FILE" ]] || { echo "Config file not found: $CONFIG_FILE" >&2; exit 1; }
source "$CONFIG_FILE"

: "${MODE:?missing MODE}"
: "${WORK_ROOT:?missing WORK_ROOT}"
: "${PREFIX:?missing PREFIX}"
: "${MAP_FILE:?missing MAP_FILE}"

RSCRIPT_BIN="${RSCRIPT_BIN:-Rscript}"
ANCHOR_COL="${ANCHOR_COL:-}"
SECONDARY_COL="${SECONDARY_COL:-B25}"
RNA_GROUPS="${RNA_GROUPS:-TC SC FC}"
DNA_THRESHOLD="${DNA_THRESHOLD:-0.99}"
RNA_THRESHOLD="${RNA_THRESHOLD:-0.90}"
IBS_ZMIN="${IBS_ZMIN:-0.7}"
IBS_ZMAX="${IBS_ZMAX:-1.0}"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
STEP01_DIR="${WORK_ROOT}/step01_prepare_matrix"
STEP02_DIR="${WORK_ROOT}/step02_build_reference_cluster"
STEP03_DIR="${WORK_ROOT}/step03_pairwise_matching"
STEP04_DIR="${WORK_ROOT}/step04_diagnose_low_ibs"
STEP05_DIR="${WORK_ROOT}/step05_cluster_rescue_and_summary"

mkdir -p "$STEP01_DIR" "$STEP02_DIR" "$STEP03_DIR" "$STEP04_DIR" "$STEP05_DIR"

if [[ "$MODE" == "vcf" ]]; then
  bash "${SCRIPT_DIR}/step00_vcf_to_ibs.sh" "$CONFIG_FILE"
  IBS_FILE="${WORK_ROOT}/step00_vcf_to_ibs/04_plink/${PREFIX}.final_qc.mibs"
  IBS_ID_FILE="${WORK_ROOT}/step00_vcf_to_ibs/04_plink/${PREFIX}.final_qc.mibs.id"
elif [[ "$MODE" == "ibs" ]]; then
  : "${IBS_FILE:?missing IBS_FILE for MODE=ibs}"
  : "${IBS_ID_FILE:?missing IBS_ID_FILE for MODE=ibs}"
else
  echo "Unsupported MODE: ${MODE}" >&2
  exit 1
fi

"$RSCRIPT_BIN" "${SCRIPT_DIR}/step01_prepare_matrix.R" \
  --mibs "$IBS_FILE" \
  --id "$IBS_ID_FILE" \
  --map "$MAP_FILE" \
  --outdir "$STEP01_DIR" \
  --prefix "$PREFIX"

"$RSCRIPT_BIN" "${SCRIPT_DIR}/step02_build_reference_cluster.R" \
  --matrix "${STEP01_DIR}/${PREFIX}_ibs_matrix.tsv" \
  --map "$MAP_FILE" \
  --outdir "$STEP02_DIR" \
  --prefix "$PREFIX" \
  --anchor-col "$ANCHOR_COL" \
  --zmin "$IBS_ZMIN" \
  --zmax "$IBS_ZMAX"

"$RSCRIPT_BIN" "${SCRIPT_DIR}/step03_pairwise_matching.R" \
  --matrix "${STEP01_DIR}/${PREFIX}_ibs_matrix.tsv" \
  --map "$MAP_FILE" \
  --outdir "$STEP03_DIR" \
  --prefix "$PREFIX" \
  --anchor-col "$ANCHOR_COL" \
  --rna-groups "$RNA_GROUPS" \
  --dna-threshold "$DNA_THRESHOLD" \
  --rna-threshold "$RNA_THRESHOLD" \
  --zmin "$IBS_ZMIN" \
  --zmax "$IBS_ZMAX"

"$RSCRIPT_BIN" "${SCRIPT_DIR}/step04_diagnose_low_ibs.R" \
  --matrix "${STEP01_DIR}/${PREFIX}_ibs_matrix.tsv" \
  --map "$MAP_FILE" \
  --pairwise-dir "$STEP03_DIR" \
  --outdir "$STEP04_DIR" \
  --prefix "$PREFIX" \
  --anchor-col "$ANCHOR_COL" \
  --zmin "$IBS_ZMIN" \
  --zmax "$IBS_ZMAX"

"$RSCRIPT_BIN" "${SCRIPT_DIR}/step05_cluster_rescue_and_summary.R" \
  --matrix "${STEP01_DIR}/${PREFIX}_ibs_matrix.tsv" \
  --map "$MAP_FILE" \
  --cluster-table "${STEP02_DIR}/${PREFIX}_cluster_table.tsv" \
  --pairwise-dir "$STEP03_DIR" \
  --outdir "$STEP05_DIR" \
  --prefix "$PREFIX" \
  --anchor-col "$ANCHOR_COL" \
  --secondary-col "$SECONDARY_COL" \
  --rna-groups "$RNA_GROUPS" \
  --dna-threshold "$DNA_THRESHOLD" \
  --rna-threshold "$RNA_THRESHOLD"

echo "Modular pipeline completed for ${PREFIX}"
