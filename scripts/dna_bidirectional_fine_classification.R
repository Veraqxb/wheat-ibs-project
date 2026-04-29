#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(key, default = NULL) {
  idx <- match(key, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[[idx + 1]]
}

root_dir <- get_arg("--root", getwd())
input_dir <- get_arg("--input-dir", file.path(root_dir, "two_group_heatmap_check", "current"))
outdir <- get_arg("--outdir", file.path(root_dir, "bidirectional_fine_classification_260427"))
threshold <- as.numeric(get_arg("--threshold", "0.99"))

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%F %T")), sprintf(...), "\n", sep = "")
}

standardize_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "N/A", "NaN", "NULL", "null", "-")] <- NA_character_
  x
}

display_id_from_anchor <- function(anchor_id, ploidy) {
  x <- standardize_id(anchor_id)
  out <- sub("^Z23", "", x)
  out <- ifelse(startsWith(out, ploidy), out, paste0(ploidy, out))
  out
}

read_cross_matrix <- function(path) {
  dt <- fread(path, data.table = FALSE, check.names = FALSE)
  sample_col <- names(dt)[1]
  rn <- standardize_id(dt[[sample_col]])
  mat <- as.matrix(dt[, -1, drop = FALSE])
  suppressWarnings(storage.mode(mat) <- "numeric")
  rownames(mat) <- rn
  colnames(mat) <- standardize_id(colnames(mat))
  mat
}

find_group_dir <- function(ploidy) {
  candidates <- c(
    file.path(input_dir, ploidy),
    file.path(input_dir, paste0(ploidy, "_2group")),
    file.path(input_dir, paste0(ploidy, "_2groups"))
  )
  pair_files <- file.path(candidates, paste0(ploidy, "_2group_pair_summary.tsv"))
  hit <- candidates[file.exists(pair_files)]
  if (length(hit) > 0) return(hit[1])

  recursive_hit <- list.files(
    input_dir,
    pattern = paste0("^", ploidy, "_2group_pair_summary\\.tsv$"),
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(recursive_hit) == 0) stop("No pair_summary found for ", ploidy, " under ", input_dir)
  dirname(recursive_hit[1])
}

best_name <- function(values, names_vec) {
  if (all(is.na(values))) return(NA_character_)
  names_vec[which.max(replace(values, is.na(values), -Inf))]
}

best_value <- function(values) {
  if (all(is.na(values))) return(NA_real_)
  max(values, na.rm = TRUE)
}

classify_one <- function(ploidy) {
  group_dir <- find_group_dir(ploidy)
  prefix <- paste0(ploidy, "_2group")
  pair_path <- file.path(group_dir, paste0(prefix, "_pair_summary.tsv"))
  cross_path <- file.path(group_dir, paste0(prefix, "_B25_x_Z23_cross_matrix.tsv"))

  if (!file.exists(pair_path)) stop("Missing pair summary: ", pair_path)
  if (!file.exists(cross_path)) stop("Missing cross matrix: ", cross_path)

  pair <- fread(pair_path)
  cross <- read_cross_matrix(cross_path)

  required_cols <- c("row_index", "anchor_id", "secondary_id", "expected_pair_ibs")
  missing_cols <- setdiff(required_cols, names(pair))
  if (length(missing_cols) > 0) {
    stop("Missing columns in ", pair_path, ": ", paste(missing_cols, collapse = ", "))
  }

  pair[, anchor_id := standardize_id(anchor_id)]
  pair[, secondary_id := standardize_id(secondary_id)]

  # B25 -> Z23: each B25 row scans all Z23 columns.
  best_z23_by_b25 <- apply(cross, 1, best_name, names_vec = colnames(cross))
  best_z23_ibs_by_b25 <- apply(cross, 1, best_value)

  # Z23 -> B25: each Z23 column scans all B25 rows.
  best_b25_by_z23 <- apply(cross, 2, best_name, names_vec = rownames(cross))
  best_b25_ibs_by_z23 <- apply(cross, 2, best_value)

  pair[, best_Z23 := best_z23_by_b25[secondary_id]]
  pair[, best_Z23_ibs := best_z23_ibs_by_b25[secondary_id]]
  pair[, best_B25 := best_b25_by_z23[anchor_id]]
  pair[, best_B25_ibs := best_b25_ibs_by_z23[anchor_id]]

  pair[, ID := display_id_from_anchor(anchor_id, ploidy)]
  pair[, expected_high := !is.na(expected_pair_ibs) & expected_pair_ibs >= threshold]
  pair[, b25_to_z23_exact := expected_high &
         !is.na(best_Z23) &
         best_Z23 == anchor_id &
         !is.na(best_Z23_ibs) &
         best_Z23_ibs >= threshold]
  pair[, z23_to_b25_exact := expected_high &
         !is.na(best_B25) &
         best_B25 == secondary_id &
         !is.na(best_B25_ibs) &
         best_B25_ibs >= threshold]
  pair[, best_shifted_high := (!is.na(best_Z23_ibs) & best_Z23_ibs >= threshold) |
         (!is.na(best_B25_ibs) & best_B25_ibs >= threshold)]

  pair[, dna_diagnosis_fine := fifelse(
    is.na(expected_pair_ibs),
    "No_data",
    fifelse(
      b25_to_z23_exact & z23_to_b25_exact,
      "Exact_match",
      fifelse(
        b25_to_z23_exact & !z23_to_b25_exact,
        "B25_to_Z23_only_exact",
        fifelse(
          !b25_to_z23_exact & z23_to_b25_exact,
          "Z23_to_B25_only_exact",
          fifelse(
            expected_high & best_shifted_high,
            "Expected_high_but_best_shifted",
            fifelse(
              !expected_high & best_shifted_high,
              "True_mismatch",
              "low_ibs_match"
            )
          )
        )
      )
    )
  )]

  # Short aliases are convenient for Excel reports, while the long class keeps
  # the scan direction explicit.
  pair[, dna_diagnosis_short := fifelse(
    dna_diagnosis_fine == "B25_to_Z23_only_exact",
    "B25_only_exact",
    fifelse(
      dna_diagnosis_fine == "Z23_to_B25_only_exact",
      "Z23_only_exact",
      dna_diagnosis_fine
    )
  )]

  pair[, ploidy := ploidy]
  pair[, source_dir := group_dir]

  out_cols <- c(
    "ID", "ploidy", "row_index",
    "anchor_id", "secondary_id",
    "best_Z23", "best_Z23_ibs",
    "best_B25", "best_B25_ibs",
    "expected_pair_ibs",
    "b25_to_z23_exact", "z23_to_b25_exact",
    "dna_diagnosis_fine", "dna_diagnosis_short",
    "source_dir"
  )
  out <- pair[, ..out_cols]

  setorder(out, row_index)
  fwrite(
    out,
    file.path(outdir, paste0(ploidy, "_bidirectional_fine_classification.tsv")),
    sep = "\t",
    quote = FALSE
  )
  out
}

all_dt <- rbindlist(lapply(c("C2", "C4", "C6"), classify_one), fill = TRUE)
fwrite(
  all_dt,
  file.path(outdir, "all_bidirectional_fine_classification.tsv"),
  sep = "\t",
  quote = FALSE
)

summary_dt <- all_dt[, .N, by = .(ploidy, dna_diagnosis_fine, dna_diagnosis_short)]
setorder(summary_dt, ploidy, dna_diagnosis_fine)
fwrite(
  summary_dt,
  file.path(outdir, "all_bidirectional_fine_classification_summary.tsv"),
  sep = "\t",
  quote = FALSE
)

attention_dt <- all_dt[!dna_diagnosis_fine %in% c("Exact_match", "low_ibs_match")]
fwrite(
  attention_dt,
  file.path(outdir, "attention_samples_bidirectional_fine.tsv"),
  sep = "\t",
  quote = FALSE
)

msg("Wrote all table: %s", file.path(outdir, "all_bidirectional_fine_classification.tsv"))
msg("Wrote summary: %s", file.path(outdir, "all_bidirectional_fine_classification_summary.tsv"))
msg("Wrote attention table: %s", file.path(outdir, "attention_samples_bidirectional_fine.tsv"))
