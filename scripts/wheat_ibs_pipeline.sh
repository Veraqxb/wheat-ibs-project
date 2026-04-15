#!/usr/bin/env bash

set -Eeuo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash wheat_ibs_pipeline.sh <config.sh>

Config variables expected in <config.sh>:
  FASTCALL_JAR           TIGER_F3 jar path
  REFERENCE_DIR          directory containing chrXX.fa.gz
  BAM_INFO               bam info table for FastCall3
  LIB_DIR                directory containing chrXX-50k.Lib.gz
  CHR_LIST               quoted list, e.g. "001 002 003"
  THREADS_FASTCALL       FastCall3 thread count
  THREADS_PARALLEL       bgzip/index parallel jobs
  JAVA_XMX               e.g. 100g
  SAMTOOLS_BIN           samtools path used by FastCall3
  BCFTOOLS_BIN           bcftools path
  BGZIP_BIN              bgzip path
  TABIX_BIN              tabix path
  PLINK_BIN              plink path
  PARALLEL_BIN           GNU parallel path
  RSCRIPT_BIN            Rscript path
  WORK_ROOT              output root dir
  PREFIX                 output prefix
  CHR_SET                e.g. 42 for wheat
  MAP_FILE               tab-delimited sample map
  GROUP_Y                y-axis group column name
  GROUP_X                x-axis group column name
Optional:
  FASTCALL_F             default 30
  FASTCALL_G             default 20
  FASTCALL_H             default 0.05
  FASTCALL_E             default 0
  FASTCALL_K_SUBDIR      default all_calls
  MIN_MAC                default 2
  MAX_GENO               default 0.2
  HIGH_HET_THRESHOLD     default 0.05
  IBS_ZMIN               default 0.7
  IBS_ZMAX               default 1.0
EOF
}

[[ $# -eq 1 ]] || { usage; exit 1; }

CONFIG_FILE="$1"
[[ -f "$CONFIG_FILE" ]] || { echo "Config not found: $CONFIG_FILE" >&2; exit 1; }
source "$CONFIG_FILE"

: "${FASTCALL_JAR:?missing FASTCALL_JAR}"
: "${REFERENCE_DIR:?missing REFERENCE_DIR}"
: "${BAM_INFO:?missing BAM_INFO}"
: "${LIB_DIR:?missing LIB_DIR}"
: "${CHR_LIST:?missing CHR_LIST}"
: "${THREADS_FASTCALL:?missing THREADS_FASTCALL}"
: "${THREADS_PARALLEL:?missing THREADS_PARALLEL}"
: "${JAVA_XMX:?missing JAVA_XMX}"
: "${SAMTOOLS_BIN:?missing SAMTOOLS_BIN}"
: "${BCFTOOLS_BIN:?missing BCFTOOLS_BIN}"
: "${BGZIP_BIN:?missing BGZIP_BIN}"
: "${TABIX_BIN:?missing TABIX_BIN}"
: "${PLINK_BIN:?missing PLINK_BIN}"
: "${PARALLEL_BIN:?missing PARALLEL_BIN}"
: "${RSCRIPT_BIN:?missing RSCRIPT_BIN}"
: "${WORK_ROOT:?missing WORK_ROOT}"
: "${PREFIX:?missing PREFIX}"
: "${CHR_SET:?missing CHR_SET}"
: "${MAP_FILE:?missing MAP_FILE}"
: "${GROUP_Y:?missing GROUP_Y}"
: "${GROUP_X:?missing GROUP_X}"

FASTCALL_F="${FASTCALL_F:-30}"
FASTCALL_G="${FASTCALL_G:-20}"
FASTCALL_H="${FASTCALL_H:-0.05}"
FASTCALL_E="${FASTCALL_E:-0}"
FASTCALL_K_SUBDIR="${FASTCALL_K_SUBDIR:-all_calls}"
MIN_MAC="${MIN_MAC:-2}"
MAX_GENO="${MAX_GENO:-0.2}"
HIGH_HET_THRESHOLD="${HIGH_HET_THRESHOLD:-0.05}"
IBS_ZMIN="${IBS_ZMIN:-0.7}"
IBS_ZMAX="${IBS_ZMAX:-1.0}"

FASTCALL_OUT_DIR="${WORK_ROOT}/${PREFIX}_FastCall3"
FASTCALL_LOG_DIR="${WORK_ROOT}/${PREFIX}_logs"
VCF_DIR="${WORK_ROOT}/${PREFIX}_VCF"
FILTER_DIR="${WORK_ROOT}/${PREFIX}_vcf_filter"
PLINK_DIR="${WORK_ROOT}/${PREFIX}_plink"
REPORT_DIR="${WORK_ROOT}/${PREFIX}_report"
PIPELINE_LOG="${WORK_ROOT}/${PREFIX}.pipeline.log"
MERGED_VCF="${VCF_DIR}/${PREFIX}.merged.vcf.gz"
BI_VCF="${FILTER_DIR}/${PREFIX}.bi_mac${MIN_MAC}.vcf.gz"
PLINK_RAW="${PLINK_DIR}/${PREFIX}.bi_mac${MIN_MAC}"
PLINK_QC1="${PLINK_DIR}/${PREFIX}.qc_step1"
PLINK_FINAL="${PLINK_DIR}/${PREFIX}.final_qc"

mkdir -p "$WORK_ROOT" "$FASTCALL_OUT_DIR" "$FASTCALL_LOG_DIR" "$VCF_DIR" "$FILTER_DIR" "$PLINK_DIR" "$REPORT_DIR"

log() {
  echo "[$(date '+%F %T')] $*" | tee -a "$PIPELINE_LOG"
}

require_file() {
  [[ -f "$1" ]] || { echo "Required file missing: $1" >&2; exit 1; }
}

require_cmd() {
  [[ -x "$1" ]] || { echo "Required executable missing or not executable: $1" >&2; exit 1; }
}

run_cmd() {
  log "$*"
  "$@"
}

log "Starting pipeline for ${PREFIX}"

for exe in "$SAMTOOLS_BIN" "$BCFTOOLS_BIN" "$BGZIP_BIN" "$TABIX_BIN" "$PLINK_BIN" "$PARALLEL_BIN" "$RSCRIPT_BIN"; do
  require_cmd "$exe"
done
require_file "$FASTCALL_JAR"
require_file "$BAM_INFO"
require_file "$MAP_FILE"

run_fastcall_chr() {
  local chr="$1"
  local chr_num
  chr_num=$((10#$chr))
  local ref="${REFERENCE_DIR}/chr${chr}.fa.gz"
  local lib="${LIB_DIR}/chr${chr}-50k.Lib.gz"
  local chr_log="${FASTCALL_LOG_DIR}/chr${chr}.log"

  require_file "$ref"
  require_file "$lib"

  log "FastCall3 scan: chr${chr}"
  java -Xmx"${JAVA_XMX}" -jar "$FASTCALL_JAR" \
    -app FastCall3 \
    -mod scan \
    -a "$ref" \
    -b "$BAM_INFO" \
    -c "$lib" \
    -d "$chr_num" \
    -e "$FASTCALL_E" \
    -f "$FASTCALL_F" \
    -g "$FASTCALL_G" \
    -h "$FASTCALL_H" \
    -i "$SAMTOOLS_BIN" \
    -j "$THREADS_FASTCALL" \
    -k "$FASTCALL_OUT_DIR" \
    >"$chr_log" 2>&1
}

collect_vcfs() {
  find "$FASTCALL_OUT_DIR" -type f \( -name "*.vcf" -o -name "*.vcf.gz" \) | sort
}

compress_and_index_vcfs() {
  log "Compressing plain VCF files"
  mapfile -t plain_vcfs < <(find "$FASTCALL_OUT_DIR" -type f -name "*.vcf" | sort)
  if [[ ${#plain_vcfs[@]} -gt 0 ]]; then
    printf '%s\n' "${plain_vcfs[@]}" | "$PARALLEL_BIN" -j "$THREADS_PARALLEL" "$BGZIP_BIN" -f {}
  fi

  log "Indexing VCF.GZ files"
  mapfile -t gz_vcfs < <(find "$FASTCALL_OUT_DIR" -type f -name "*.vcf.gz" | sort)
  [[ ${#gz_vcfs[@]} -gt 0 ]] || { echo "No VCF.gz files found in $FASTCALL_OUT_DIR" >&2; exit 1; }
  printf '%s\n' "${gz_vcfs[@]}" | "$PARALLEL_BIN" -j "$THREADS_PARALLEL" "$BCFTOOLS_BIN" index -f {}
}

copy_vcfs_to_workdir() {
  log "Collecting per-chromosome VCFs"
  mapfile -t gz_vcfs < <(find "$FASTCALL_OUT_DIR" -type f -name "*.vcf.gz" | sort)
  [[ ${#gz_vcfs[@]} -gt 0 ]] || { echo "No VCF.gz files to collect" >&2; exit 1; }
  for vcf in "${gz_vcfs[@]}"; do
    cp "$vcf" "$VCF_DIR/"
    [[ -f "${vcf}.csi" ]] && cp "${vcf}.csi" "$VCF_DIR/"
    [[ -f "${vcf}.tbi" ]] && cp "${vcf}.tbi" "$VCF_DIR/"
  done
}

concat_vcfs() {
  log "Concatenating chromosome VCFs"
  (
    cd "$VCF_DIR"
    mapfile -t local_vcfs < <(find . -maxdepth 1 -type f -name "*.vcf.gz" | sort)
    [[ ${#local_vcfs[@]} -gt 0 ]] || { echo "No copied VCF.gz files in $VCF_DIR" >&2; exit 1; }
    "$BCFTOOLS_BIN" concat "${local_vcfs[@]}" -Oz -o "$MERGED_VCF"
  )
  run_cmd "$BCFTOOLS_BIN" index -f "$MERGED_VCF"
}

filter_vcf() {
  log "Filtering bi-allelic SNPs with MAC >= ${MIN_MAC}"
  run_cmd "$BCFTOOLS_BIN" view -m2 -M2 -v snps -i "MAC>=${MIN_MAC}" "$MERGED_VCF" -Oz -o "$BI_VCF"
  run_cmd "$TABIX_BIN" -f -p vcf "$BI_VCF"
}

run_plink_qc() {
  log "Converting to PLINK"
  run_cmd "$PLINK_BIN" --vcf "$BI_VCF" --make-bed --chr-set "$CHR_SET" --out "$PLINK_RAW" --double-id

  log "Filtering SNP missingness with --geno ${MAX_GENO}"
  run_cmd "$PLINK_BIN" --bfile "$PLINK_RAW" --chr-set "$CHR_SET" --geno "$MAX_GENO" --make-bed --out "$PLINK_QC1"

  log "Computing HWE summary for high-heterozygosity SNP removal"
  run_cmd "$PLINK_BIN" --bfile "$PLINK_QC1" --chr-set "$CHR_SET" --hardy --out "${PLINK_DIR}/${PREFIX}.hardy"

  awk -v thr="$HIGH_HET_THRESHOLD" 'NR>1 && $7 > thr {print $2}' "${PLINK_DIR}/${PREFIX}.hardy.hwe" > "${PLINK_DIR}/${PREFIX}.high_het_snps.txt"

  log "Removing high-heterozygosity SNPs"
  run_cmd "$PLINK_BIN" --bfile "$PLINK_QC1" --chr-set "$CHR_SET" --exclude "${PLINK_DIR}/${PREFIX}.high_het_snps.txt" --make-bed --out "$PLINK_FINAL"

  log "Calculating IBS"
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
    --zmax "$IBS_ZMAX" \
    --het-threshold "$HIGH_HET_THRESHOLD"
}

for chr in $CHR_LIST; do
  run_fastcall_chr "$chr"
done

compress_and_index_vcfs
copy_vcfs_to_workdir
concat_vcfs
filter_vcf
run_plink_qc
run_report

log "Pipeline completed successfully"
