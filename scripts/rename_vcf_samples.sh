#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash rename_vcf_samples.sh <vcf_dir> <rename_map.tsv> [bcftools_bin]

Arguments:
  vcf_dir          directory containing chr*.vcf.gz
  rename_map.tsv   two-column tab-delimited file: old_sample<TAB>new_sample
  bcftools_bin     optional, default: bcftools

This script:
  1. checks sample names in each chr*.vcf.gz
  2. applies bcftools reheader using the rename table
  3. writes renamed VCFs into <vcf_dir>/renamed/
  4. indexes the renamed VCFs

Notes:
  - keep original VCF files unchanged
  - use the renamed directory as VCF_SOURCE_DIR in downstream IBS analysis
EOF
}

[[ $# -ge 2 ]] || { usage; exit 1; }

VCF_DIR="$1"
RENAME_MAP="$2"
BCFTOOLS_BIN="${3:-bcftools}"
OUT_DIR="${VCF_DIR%/}/renamed"

[[ -d "$VCF_DIR" ]] || { echo "VCF directory not found: $VCF_DIR" >&2; exit 1; }
[[ -f "$RENAME_MAP" ]] || { echo "Rename map not found: $RENAME_MAP" >&2; exit 1; }
[[ -x "$BCFTOOLS_BIN" || "$BCFTOOLS_BIN" == "bcftools" ]] || { echo "bcftools not executable: $BCFTOOLS_BIN" >&2; exit 1; }

mkdir -p "$OUT_DIR"

echo "[1/4] Checking rename map format..."
awk 'NF < 2 {bad=1} END {exit bad}' "$RENAME_MAP" || {
  echo "Rename map must be two-column tab-delimited: old_sample new_sample" >&2
  exit 1
}

echo "[2/4] Validating sample names against VCFs..."
first_vcf="$(find "$VCF_DIR" -maxdepth 1 -type f -name 'chr*.vcf.gz' | sort | head -n 1)"
[[ -n "$first_vcf" ]] || { echo "No chr*.vcf.gz found in $VCF_DIR" >&2; exit 1; }

"$BCFTOOLS_BIN" query -l "$first_vcf" > "${OUT_DIR}/current_samples.txt"
cut -f1 "$RENAME_MAP" > "${OUT_DIR}/rename_old.txt"

missing_count="$(grep -Fvx -f "${OUT_DIR}/current_samples.txt" "${OUT_DIR}/rename_old.txt" | wc -l | tr -d ' ')"
if [[ "$missing_count" != "0" ]]; then
  echo "These old sample IDs are not present in ${first_vcf}:" >&2
  grep -Fvx -f "${OUT_DIR}/current_samples.txt" "${OUT_DIR}/rename_old.txt" >&2
  exit 1
fi

echo "[3/4] Reheadering chromosome VCF files..."
find "$VCF_DIR" -maxdepth 1 -type f -name 'chr*.vcf.gz' | sort | while read -r vcf; do
  base="$(basename "$vcf")"
  out_vcf="${OUT_DIR}/${base}"
  echo "  - ${base}"
  "$BCFTOOLS_BIN" reheader -s "$RENAME_MAP" "$vcf" -o "$out_vcf"
  "$BCFTOOLS_BIN" index -f "$out_vcf"
done

echo "[4/4] Done."
echo "Renamed VCF directory: ${OUT_DIR}"
echo "Use this directory as VCF_SOURCE_DIR in the IBS config."
