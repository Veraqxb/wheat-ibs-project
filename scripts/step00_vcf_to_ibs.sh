#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash step00_vcf_to_ibs.sh <config.sh>

Config variables:
  PROJECT_ROOT
  VCF_SOURCE_DIR
  WORK_ROOT
  PREFIX
  MAP_FILE
  BCFTOOLS_BIN
  BGZIP_BIN
  TABIX_BIN
  PLINK_BIN
  PARALLEL_BIN

Optional:
  THREADS_PARALLEL   default 8
  CHR_SET            default 42
  CHR_LIST           optional chromosome whitelist
  MIN_MAC            default 2
  MAX_GENO           default 0.3
  HIGH_HET_THRESHOLD default 0.05
  COPY_MODE          default link
  ALL_GROUPS         optional map groups for missingness profiles
  GROUP_X/GROUP_Y/SECONDARY_GROUP and *_MATCH_* settings for keep-file generation
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
: "${BCFTOOLS_BIN:?missing BCFTOOLS_BIN}"
: "${BGZIP_BIN:?missing BGZIP_BIN}"
: "${TABIX_BIN:?missing TABIX_BIN}"
: "${PLINK_BIN:?missing PLINK_BIN}"
: "${PARALLEL_BIN:?missing PARALLEL_BIN}"

THREADS_PARALLEL="${THREADS_PARALLEL:-8}"
CHR_SET="${CHR_SET:-42}"
CHR_LIST="${CHR_LIST:-}"
MIN_MAC="${MIN_MAC:-2}"
MAX_GENO="${MAX_GENO:-0.3}"
HIGH_HET_THRESHOLD="${HIGH_HET_THRESHOLD:-0.05}"
COPY_MODE="${COPY_MODE:-link}"
ALL_GROUPS="${ALL_GROUPS:-}"
GROUP_X="${GROUP_X:-}"
GROUP_Y="${GROUP_Y:-}"
SECONDARY_GROUP="${SECONDARY_GROUP:-}"
GROUP_X_MATCH_COL="${GROUP_X_MATCH_COL:-}"
GROUP_Y_MATCH_COL="${GROUP_Y_MATCH_COL:-}"
SECONDARY_GROUP_MATCH_COL="${SECONDARY_GROUP_MATCH_COL:-}"
GROUP_X_MATCH_MODE="${GROUP_X_MATCH_MODE:-direct}"
GROUP_Y_MATCH_MODE="${GROUP_Y_MATCH_MODE:-direct}"
SECONDARY_GROUP_MATCH_MODE="${SECONDARY_GROUP_MATCH_MODE:-direct}"

STEP_DIR="${WORK_ROOT}/step00_vcf_to_ibs"
INPUT_DIR="${STEP_DIR}/01_input_vcf"
MERGE_DIR="${STEP_DIR}/02_merged_vcf"
FILTER_DIR="${STEP_DIR}/03_filter"
PLINK_DIR="${STEP_DIR}/04_plink"
LOG_DIR="${STEP_DIR}/logs"
mkdir -p "$INPUT_DIR" "$MERGE_DIR" "$FILTER_DIR" "$PLINK_DIR" "$LOG_DIR"

PIPELINE_LOG="${LOG_DIR}/${PREFIX}.step00.log"
MERGED_VCF="${MERGE_DIR}/${PREFIX}.merged.vcf.gz"
FILTERED_VCF="${FILTER_DIR}/${PREFIX}.bi_mac${MIN_MAC}.vcf.gz"
PLINK_RAW="${PLINK_DIR}/${PREFIX}.bi_mac${MIN_MAC}"
PLINK_QC1="${PLINK_DIR}/${PREFIX}.qc_step1"
PLINK_FINAL="${PLINK_DIR}/${PREFIX}.final_qc"
MISSING_PREFIX="${PLINK_DIR}/${PREFIX}.missing"

log() {
  echo "[$(date '+%F %T')] $*" | tee -a "$PIPELINE_LOG"
}

run_cmd() {
  log "$*"
  "$@"
}

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

contains_word() {
  local needle="$1"
  shift || true
  local item
  for item in "$@"; do
    [[ "$item" == "$needle" ]] && return 0
  done
  return 1
}

get_group_match_settings() {
  local group_name="$1"
  local match_col=""
  local match_mode="direct"
  if [[ "$group_name" == "$GROUP_X" ]]; then
    match_col="$GROUP_X_MATCH_COL"
    match_mode="$GROUP_X_MATCH_MODE"
  elif [[ "$group_name" == "$GROUP_Y" ]]; then
    match_col="$GROUP_Y_MATCH_COL"
    match_mode="$GROUP_Y_MATCH_MODE"
  elif [[ -n "$SECONDARY_GROUP" && "$group_name" == "$SECONDARY_GROUP" ]]; then
    match_col="$SECONDARY_GROUP_MATCH_COL"
    match_mode="$SECONDARY_GROUP_MATCH_MODE"
  fi
  printf '%s\t%s\n' "$match_col" "$match_mode"
}

build_group_keep_file() {
  local group_name="$1"
  local keep_file="$2"
  local settings
  settings="$(get_group_match_settings "$group_name")"
  local match_col="${settings%%$'\t'*}"
  local match_mode="${settings#*$'\t'}"

  awk -v group_col="$group_name" \
      -v match_col="$match_col" \
      -v match_mode="$match_mode" \
      -v fam_file="${PLINK_RAW}.fam" '
    BEGIN {
      FS = OFS = "\t"
      while ((getline < fam_file) > 0) {
        present[$2] = 1
      }
      close(fam_file)
    }
    NR == 1 {
      for (i = 1; i <= NF; i++) idx[$i] = i
      if (!(group_col in idx)) exit 2
      if (match_col != "" && !(match_col in idx)) exit 3
      next
    }
    {
      if (match_col != "") sample_id = $(idx[match_col]); else sample_id = $(idx[group_col])
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", sample_id)
      if (sample_id == "" || sample_id == "NA" || sample_id == "-") next
      if (match_mode == "b25_to_2") sub(/^B25C2_/, "2_", sample_id)
      if (sample_id in present) print sample_id, sample_id
    }
  ' "$MAP_FILE" | awk '!seen[$0]++' > "$keep_file"
}

run_group_missingness() {
  local group_name="$1"
  local keep_file="${PLINK_DIR}/${PREFIX}.${group_name}.keep"
  local miss_out="${MISSING_PREFIX}.${group_name}"
  build_group_keep_file "$group_name" "$keep_file"
  [[ -s "$keep_file" ]] || { log "Skipping missingness for ${group_name}: no matching samples"; return; }
  [[ -f "${miss_out}.lmiss" ]] && { log "Missingness for ${group_name} already exists"; return; }
  run_cmd "$PLINK_BIN" --bfile "$PLINK_RAW" --chr-set "$CHR_SET" --keep "$keep_file" --missing --out "$miss_out"
}

log "step00_vcf_to_ibs started for ${PREFIX}"

mapfile -t plain_vcfs < <(collect_expected_vcfs "$VCF_SOURCE_DIR" "vcf")
mapfile -t gz_vcfs < <(collect_expected_vcfs "$VCF_SOURCE_DIR" "vcf.gz")
if [[ ${#plain_vcfs[@]} -gt 0 ]]; then
  printf '%s\n' "${plain_vcfs[@]}" | "$PARALLEL_BIN" -j "$THREADS_PARALLEL" "$BGZIP_BIN" -f {}
fi
mapfile -t gz_vcfs < <(collect_expected_vcfs "$VCF_SOURCE_DIR" "vcf.gz")
[[ ${#gz_vcfs[@]} -gt 0 ]] || { echo "No chromosome VCF files found in $VCF_SOURCE_DIR" >&2; exit 1; }

missing_index_vcfs=()
for vcf in "${gz_vcfs[@]}"; do
  if [[ ! -f "${vcf}.csi" && ! -f "${vcf}.tbi" ]]; then
    missing_index_vcfs+=("$vcf")
  fi
done
if [[ ${#missing_index_vcfs[@]} -gt 0 ]]; then
  printf '%s\n' "${missing_index_vcfs[@]}" | "$PARALLEL_BIN" -j "$THREADS_PARALLEL" "$BCFTOOLS_BIN" index -f {}
fi

rm -f "${INPUT_DIR}"/*.vcf.gz "${INPUT_DIR}"/*.vcf.gz.csi "${INPUT_DIR}"/*.vcf.gz.tbi "${INPUT_DIR}/vcf.list" 2>/dev/null || true
for vcf in "${gz_vcfs[@]}"; do
  base=$(basename "$vcf")
  if [[ "$COPY_MODE" == "link" ]]; then
    ln -sf "$vcf" "${INPUT_DIR}/${base}"
    [[ -f "${vcf}.csi" ]] && ln -sf "${vcf}.csi" "${INPUT_DIR}/${base}.csi"
    [[ -f "${vcf}.tbi" ]] && ln -sf "${vcf}.tbi" "${INPUT_DIR}/${base}.tbi"
  else
    cp "$vcf" "${INPUT_DIR}/${base}"
    [[ -f "${vcf}.csi" ]] && cp "${vcf}.csi" "${INPUT_DIR}/${base}.csi"
    [[ -f "${vcf}.tbi" ]] && cp "${vcf}.tbi" "${INPUT_DIR}/${base}.tbi"
  fi
done

if [[ ! -f "$MERGED_VCF" ]]; then
  (
    cd "$INPUT_DIR"
    if [[ -n "$CHR_LIST" ]]; then
      : > vcf.list
      for chr in $CHR_LIST; do
        [[ -f "chr${chr}.vcf.gz" ]] && printf 'chr%s.vcf.gz\n' "$chr" >> vcf.list
      done
    else
      find . -maxdepth 1 -type f -name "chr*.vcf.gz" | sed 's|^./||' | sort > vcf.list
    fi
    [[ -s vcf.list ]] || { echo "No VCF files collected into ${INPUT_DIR}" >&2; exit 1; }
    run_cmd "$BCFTOOLS_BIN" concat --threads "$THREADS_PARALLEL" -f vcf.list -Oz -o "$MERGED_VCF"
  )
  run_cmd "$BCFTOOLS_BIN" index -f "$MERGED_VCF"
fi

if [[ ! -f "$FILTERED_VCF" ]]; then
  run_cmd "$BCFTOOLS_BIN" view --threads "$THREADS_PARALLEL" -m2 -M2 -v snps -i "MAC>=${MIN_MAC}" "$MERGED_VCF" -Oz -o "$FILTERED_VCF"
  run_cmd "$TABIX_BIN" -f -p vcf "$FILTERED_VCF"
fi

if [[ ! -f "${PLINK_RAW}.bed" ]]; then
  run_cmd "$PLINK_BIN" --vcf "$FILTERED_VCF" --make-bed --chr-set "$CHR_SET" --out "$PLINK_RAW" --double-id
fi

groups=()
for group_name in $ALL_GROUPS; do
  contains_word "$group_name" "${groups[@]}" || groups+=("$group_name")
done
[[ -n "$GROUP_X" ]] && contains_word "$GROUP_X" "${groups[@]}" || groups+=("$GROUP_X")
[[ -n "$GROUP_Y" ]] && contains_word "$GROUP_Y" "${groups[@]}" || groups+=("$GROUP_Y")
if [[ -n "$SECONDARY_GROUP" ]]; then
  contains_word "$SECONDARY_GROUP" "${groups[@]}" || groups+=("$SECONDARY_GROUP")
fi
for group_name in "${groups[@]}"; do
  [[ -n "$group_name" ]] && run_group_missingness "$group_name"
done

if [[ ! -f "${PLINK_QC1}.bed" ]]; then
  run_cmd "$PLINK_BIN" --bfile "$PLINK_RAW" --chr-set "$CHR_SET" --geno "$MAX_GENO" --make-bed --out "$PLINK_QC1"
fi

if [[ ! -f "${PLINK_DIR}/${PREFIX}.hardy.hwe" ]]; then
  run_cmd "$PLINK_BIN" --bfile "$PLINK_QC1" --chr-set "$CHR_SET" --hardy --out "${PLINK_DIR}/${PREFIX}.hardy"
fi

awk -v thr="$HIGH_HET_THRESHOLD" 'NR>1 && $7 > thr {print $2}' "${PLINK_DIR}/${PREFIX}.hardy.hwe" > "${PLINK_DIR}/${PREFIX}.high_het_snps.txt"

if [[ ! -f "${PLINK_FINAL}.bed" ]]; then
  run_cmd "$PLINK_BIN" --bfile "$PLINK_QC1" --chr-set "$CHR_SET" --exclude "${PLINK_DIR}/${PREFIX}.high_het_snps.txt" --make-bed --out "$PLINK_FINAL"
fi

if [[ ! -f "${PLINK_FINAL}.mibs" || ! -f "${PLINK_FINAL}.mibs.id" ]]; then
  run_cmd "$PLINK_BIN" --bfile "$PLINK_FINAL" --distance ibs square --chr-set "$CHR_SET" --out "$PLINK_FINAL"
fi

log "step00_vcf_to_ibs completed successfully"
