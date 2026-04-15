#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash wheat_ibs_pipeline_server.sh <config.sh>

Required variables in config:
  PROJECT_ROOT
  VCF_SOURCE_DIR
  WORK_ROOT
  PREFIX
  MAP_FILE
  GROUP_Y
  GROUP_X
  BCFTOOLS_BIN
  BGZIP_BIN
  TABIX_BIN
  PLINK_BIN
  PARALLEL_BIN
  RSCRIPT_BIN

Optional variables:
  THREADS_PARALLEL    default 8
  CHR_SET             default 42
  CHR_LIST            optional chromosome whitelist, e.g. "001 002 003"
  GROUP_MODE          optional label such as 2group or 5group
  PLOIDY_TAG          optional label such as C2/C4/C6
  ALL_GROUPS          optional group list for record keeping
  MIN_MAC             default 2
  MAX_GENO            default 0.2
  HIGH_HET_THRESHOLD  default 0.05
  IBS_ZMIN            default 0.7
  IBS_ZMAX            default 1.0
  COPY_MODE           default copy; can be copy or link
EOF
}

[[ $# -eq 1 ]] || { usage; exit 1; }

CONFIG_FILE="$1"
[[ -f "$CONFIG_FILE" ]] || { echo "Config file not found: $CONFIG_FILE" >&2; exit 1; }
source "$CONFIG_FILE"

: "${PROJECT_ROOT:?missing PROJECT_ROOT}"
: "${VCF_SOURCE_DIR:?missing VCF_SOURCE_DIR}"
: "${WORK_ROOT:?missing WORK_ROOT}"
: "${PREFIX:?missing PREFIX}"
: "${MAP_FILE:?missing MAP_FILE}"
: "${GROUP_Y:?missing GROUP_Y}"
: "${GROUP_X:?missing GROUP_X}"
: "${BCFTOOLS_BIN:?missing BCFTOOLS_BIN}"
: "${BGZIP_BIN:?missing BGZIP_BIN}"
: "${TABIX_BIN:?missing TABIX_BIN}"
: "${PLINK_BIN:?missing PLINK_BIN}"
: "${PARALLEL_BIN:?missing PARALLEL_BIN}"
: "${RSCRIPT_BIN:?missing RSCRIPT_BIN}"

THREADS_PARALLEL="${THREADS_PARALLEL:-8}"
CHR_SET="${CHR_SET:-42}"
MIN_MAC="${MIN_MAC:-2}"
MAX_GENO="${MAX_GENO:-0.2}"
HIGH_HET_THRESHOLD="${HIGH_HET_THRESHOLD:-0.05}"
IBS_ZMIN="${IBS_ZMIN:-0.7}"
IBS_ZMAX="${IBS_ZMAX:-1.0}"
COPY_MODE="${COPY_MODE:-copy}"
CHR_LIST="${CHR_LIST:-}"
GROUP_MODE="${GROUP_MODE:-}"
PLOIDY_TAG="${PLOIDY_TAG:-}"
ALL_GROUPS="${ALL_GROUPS:-}"

PIPELINE_DIR="$WORK_ROOT"
INPUT_DIR="${PIPELINE_DIR}/01_input_vcf"
MERGE_DIR="${PIPELINE_DIR}/02_merged_vcf"
FILTER_DIR="${PIPELINE_DIR}/03_filter"
PLINK_DIR="${PIPELINE_DIR}/04_plink"
REPORT_DIR="${PIPELINE_DIR}/05_report"
LOG_DIR="${PIPELINE_DIR}/logs"

PIPELINE_LOG="${LOG_DIR}/${PREFIX}.pipeline.log"
CONCAT_ERR_LOG="${LOG_DIR}/${PREFIX}.concat.err"
MERGED_VCF="${MERGE_DIR}/${PREFIX}.merged.vcf.gz"
BI_VCF="${FILTER_DIR}/${PREFIX}.bi_mac${MIN_MAC}.vcf.gz"
PLINK_RAW="${PLINK_DIR}/${PREFIX}.bi_mac${MIN_MAC}"
PLINK_QC1="${PLINK_DIR}/${PREFIX}.qc_step1"
PLINK_FINAL="${PLINK_DIR}/${PREFIX}.final_qc"

mkdir -p "$INPUT_DIR" "$MERGE_DIR" "$FILTER_DIR" "$PLINK_DIR" "$REPORT_DIR" "$LOG_DIR"

log() {
  echo "[$(date '+%F %T')] $*" | tee -a "$PIPELINE_LOG"
}

require_file() {
  [[ -f "$1" ]] || { echo "Required file missing: $1" >&2; exit 1; }
}

require_exec() {
  [[ -x "$1" ]] || { echo "Required executable missing: $1" >&2; exit 1; }
}

run_cmd() {
  log "$*"
  "$@"
}

for exe in "$BCFTOOLS_BIN" "$BGZIP_BIN" "$TABIX_BIN" "$PLINK_BIN" "$PARALLEL_BIN" "$RSCRIPT_BIN"; do
  require_exec "$exe"
done

[[ -d "$VCF_SOURCE_DIR" ]] || { echo "VCF source dir not found: $VCF_SOURCE_DIR" >&2; exit 1; }
require_file "$MAP_FILE"

log "Pipeline started for ${PREFIX}"
log "Source VCF directory: ${VCF_SOURCE_DIR}"
log "Working directory: ${WORK_ROOT}"
[[ -n "$PLOIDY_TAG" ]] && log "Ploidy tag: ${PLOIDY_TAG}"
[[ -n "$GROUP_MODE" ]] && log "Group mode: ${GROUP_MODE}"
[[ -n "$ALL_GROUPS" ]] && log "All groups: ${ALL_GROUPS}"

collect_expected_vcfs() {
  local search_dir="$1"
  local pattern="$2"

  if [[ -n "$CHR_LIST" ]]; then
    local chr
    for chr in $CHR_LIST; do
      find "$search_dir" -maxdepth 1 -type f -name "chr${chr}.${pattern}" | sort
    done
  else
    find "$search_dir" -maxdepth 1 -type f -name "chr*.${pattern}" | sort
  fi
}

prepare_vcfs() {
  log "Preparing input VCF files"

  mapfile -t plain_vcfs < <(collect_expected_vcfs "$VCF_SOURCE_DIR" "vcf")
  mapfile -t gz_vcfs < <(collect_expected_vcfs "$VCF_SOURCE_DIR" "vcf.gz")

  if [[ ${#plain_vcfs[@]} -eq 0 && ${#gz_vcfs[@]} -eq 0 ]]; then
    echo "No expected chromosome VCF files found in $VCF_SOURCE_DIR" >&2
    exit 1
  fi

  if [[ ${#plain_vcfs[@]} -gt 0 ]]; then
    log "Compressing plain VCF files in source directory"
    printf '%s\n' "${plain_vcfs[@]}" | "$PARALLEL_BIN" -j "$THREADS_PARALLEL" "$BGZIP_BIN" -f {}
  fi

  mapfile -t gz_vcfs < <(collect_expected_vcfs "$VCF_SOURCE_DIR" "vcf.gz")
  [[ ${#gz_vcfs[@]} -gt 0 ]] || { echo "No expected chromosome .vcf.gz files found after compression" >&2; exit 1; }

  log "Indexing source VCF files"
  printf '%s\n' "${gz_vcfs[@]}" | "$PARALLEL_BIN" -j "$THREADS_PARALLEL" "$BCFTOOLS_BIN" index -f {}

  log "Collecting VCF files into analysis directory"
  rm -f "${INPUT_DIR}"/*.vcf.gz "${INPUT_DIR}"/*.vcf.gz.csi "${INPUT_DIR}"/*.vcf.gz.tbi "${INPUT_DIR}/vcf.list" 2>/dev/null || true

  for vcf in "${gz_vcfs[@]}"; do
    base=$(basename "$vcf")

    if [[ "$COPY_MODE" == "link" ]]; then
      ln -sf "$vcf" "${INPUT_DIR}/${base}"
      if [[ -f "${vcf}.csi" ]]; then
        ln -sf "${vcf}.csi" "${INPUT_DIR}/${base}.csi"
      fi
      if [[ -f "${vcf}.tbi" ]]; then
        ln -sf "${vcf}.tbi" "${INPUT_DIR}/${base}.tbi"
      fi
    else
      cp "$vcf" "${INPUT_DIR}/${base}"
      if [[ -f "${vcf}.csi" ]]; then
        cp "${vcf}.csi" "${INPUT_DIR}/${base}.csi"
      fi
      if [[ -f "${vcf}.tbi" ]]; then
        cp "${vcf}.tbi" "${INPUT_DIR}/${base}.tbi"
      fi
    fi
  done
}

concat_vcfs() {
  log "Concatenating chromosome VCF files"

  (
    cd "$INPUT_DIR"
    : > vcf.list
    if [[ -n "$CHR_LIST" ]]; then
      local chr
      for chr in $CHR_LIST; do
        [[ -f "chr${chr}.vcf.gz" ]] && printf 'chr%s.vcf.gz\n' "$chr" >> vcf.list
      done
    else
      find . -maxdepth 1 -type f -name "chr*.vcf.gz" | sed 's|^./||' | sort > vcf.list
    fi
    [[ -s vcf.list ]] || { echo "No chromosome VCF files found in $INPUT_DIR" >&2; exit 1; }

    "$BCFTOOLS_BIN" concat -f vcf.list -Oz -o "$MERGED_VCF" 2> "$CONCAT_ERR_LOG"
  )

  run_cmd "$BCFTOOLS_BIN" index -f "$MERGED_VCF"
}

filter_vcf() {
  log "Filtering SNPs: bi-allelic and MAC >= ${MIN_MAC}"
  run_cmd "$BCFTOOLS_BIN" view -m2 -M2 -v snps -i "MAC>=${MIN_MAC}" "$MERGED_VCF" -Oz -o "$BI_VCF"
  run_cmd "$TABIX_BIN" -f -p vcf "$BI_VCF"
}

run_plink_qc() {
  log "Converting VCF to PLINK bed"
  run_cmd "$PLINK_BIN" --vcf "$BI_VCF" --make-bed --chr-set "$CHR_SET" --out "$PLINK_RAW" --double-id

  log "Filtering variants by missing rate: --geno ${MAX_GENO}"
  run_cmd "$PLINK_BIN" --bfile "$PLINK_RAW" --chr-set "$CHR_SET" --geno "$MAX_GENO" --make-bed --out "$PLINK_QC1"

  log "Calculating Hardy-Weinberg statistics"
  run_cmd "$PLINK_BIN" --bfile "$PLINK_QC1" --chr-set "$CHR_SET" --hardy --out "${PLINK_DIR}/${PREFIX}.hardy"

  awk -v thr="$HIGH_HET_THRESHOLD" 'NR>1 && $7 > thr {print $2}' "${PLINK_DIR}/${PREFIX}.hardy.hwe" > "${PLINK_DIR}/${PREFIX}.high_het_snps.txt"

  log "Removing high-heterozygosity SNPs"
  run_cmd "$PLINK_BIN" --bfile "$PLINK_QC1" --chr-set "$CHR_SET" --exclude "${PLINK_DIR}/${PREFIX}.high_het_snps.txt" --make-bed --out "$PLINK_FINAL"

  log "Calculating IBS matrix"
  run_cmd "$PLINK_BIN" --bfile "$PLINK_FINAL" --distance ibs square --chr-set "$CHR_SET" --out "$PLINK_FINAL"
}

run_report() {
  log "Generating IBS report and heatmaps"
  run_cmd "$RSCRIPT_BIN" "$(dirname "$0")/wheat_ibs_report.R" \
    --mibs "${PLINK_FINAL}.mibs" \
    --id "${PLINK_FINAL}.mibs.id" \
    --map "$MAP_FILE" \
    --group-y "$GROUP_Y" \
    --group-x "$GROUP_X" \
    --outdir "$REPORT_DIR" \
    --prefix "$PREFIX" \
    --zmin "$IBS_ZMIN" \
    --zmax "$IBS_ZMAX"
}

prepare_vcfs
concat_vcfs
filter_vcf
run_plink_qc
run_report

log "Pipeline completed successfully"
