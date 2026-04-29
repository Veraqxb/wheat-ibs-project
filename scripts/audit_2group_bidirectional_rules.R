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

root_dir <- get_arg("--root", "/Users/veraqiu/Desktop/CAMP_sample_identify")
outdir <- get_arg("--outdir", file.path(root_dir, "two_group_rule_audit"))
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

find_latest_group_dir <- function(root_dir, ploidy) {
  base <- file.path(root_dir, "two_group_heatmap_check")
  candidates <- list.dirs(base, recursive = FALSE, full.names = TRUE)
  candidates <- candidates[grepl(paste0("^", ploidy, "($|_)"), basename(candidates))]
  pair_files <- file.path(candidates, paste0(ploidy, "_2group_pair_summary.tsv"))
  candidates <- candidates[file.exists(pair_files)]
  if (length(candidates) == 0) stop("No pair_summary directory found for ", ploidy)
  info <- file.info(pair_files[file.exists(pair_files)])
  candidates[which.max(info$mtime)]
}

audit_one <- function(ploidy) {
  group_dir <- find_latest_group_dir(root_dir, ploidy)
  prefix <- paste0(ploidy, "_2group")
  pair_path <- file.path(group_dir, paste0(prefix, "_pair_summary.tsv"))
  cross_path <- file.path(group_dir, paste0(prefix, "_B25_x_Z23_cross_matrix.tsv"))

  if (!file.exists(pair_path)) stop("Missing pair summary: ", pair_path)
  if (!file.exists(cross_path)) stop("Missing cross matrix: ", cross_path)

  pair <- fread(pair_path)
  cross <- read_cross_matrix(cross_path)

  pair[, anchor_id := standardize_id(anchor_id)]
  pair[, secondary_id := standardize_id(secondary_id)]
  pair[, best_anchor_raw := standardize_id(best_anchor_raw)]

  # B25 -> Z23: each row scans all Z23 anchors.
  best_anchor_recalc <- apply(cross, 1, function(x) {
    if (all(is.na(x))) return(NA_character_)
    colnames(cross)[which.max(replace(x, is.na(x), -Inf))]
  })
  best_anchor_ibs_recalc <- apply(cross, 1, function(x) {
    if (all(is.na(x))) return(NA_real_)
    max(x, na.rm = TRUE)
  })

  # Z23 -> B25: each column scans all B25 samples.
  best_secondary_recalc <- apply(cross, 2, function(x) {
    if (all(is.na(x))) return(NA_character_)
    rownames(cross)[which.max(replace(x, is.na(x), -Inf))]
  })
  best_secondary_ibs_recalc <- apply(cross, 2, function(x) {
    if (all(is.na(x))) return(NA_real_)
    max(x, na.rm = TRUE)
  })

  pair[, best_anchor_recalc := best_anchor_recalc[secondary_id]]
  pair[, best_anchor_ibs_recalc := best_anchor_ibs_recalc[secondary_id]]
  pair[, best_secondary_recalc := best_secondary_recalc[anchor_id]]
  pair[, best_secondary_ibs_recalc := best_secondary_ibs_recalc[anchor_id]]

  pair[, expected_high := !is.na(expected_pair_ibs) & expected_pair_ibs >= threshold]
  pair[, b25_to_z23_exact := expected_high &
         !is.na(best_anchor_recalc) &
         best_anchor_recalc == anchor_id &
         !is.na(best_anchor_ibs_recalc) &
         best_anchor_ibs_recalc >= threshold]
  pair[, z23_to_b25_exact := expected_high &
         !is.na(best_secondary_recalc) &
         best_secondary_recalc == secondary_id &
         !is.na(best_secondary_ibs_recalc) &
         best_secondary_ibs_recalc >= threshold]
  pair[, bidirectional_status := fifelse(
    is.na(expected_pair_ibs), "No_data",
    fifelse(b25_to_z23_exact & z23_to_b25_exact, "Bidirectional_exact",
      fifelse(b25_to_z23_exact & !z23_to_b25_exact, "Only_B25_to_Z23_exact",
        fifelse(!b25_to_z23_exact & z23_to_b25_exact, "Only_Z23_to_B25_exact",
          fifelse(expected_high, "Expected_high_but_best_shifted", "Low_expected_pair")
        )
      )
    )
  )]

  pair[, exact_rule_overlap := fifelse(
    diagnosis_label == "Exact_match" & bidirectional_status != "Bidirectional_exact",
    "Current_Exact_but_not_bidirectional",
    "OK"
  )]
  pair[, ploidy := ploidy]
  pair[, source_dir := group_dir]

  out_cols <- c(
    "ploidy", "row_index", "anchor_id", "secondary_id",
    "expected_pair_ibs",
    "best_anchor_raw", "best_anchor_ibs",
    "best_anchor_recalc", "best_anchor_ibs_recalc",
    "best_secondary_recalc", "best_secondary_ibs_recalc",
    "diagnosis_label", "expected_high",
    "b25_to_z23_exact", "z23_to_b25_exact",
    "bidirectional_status", "exact_rule_overlap", "source_dir"
  )
  pair[, ..out_cols]
}

all_audit <- rbindlist(lapply(c("C2", "C4", "C6"), audit_one), fill = TRUE)

fwrite(
  all_audit,
  file.path(outdir, "all_2group_bidirectional_rule_audit.tsv"),
  sep = "\t",
  quote = FALSE
)

summary_dt <- all_audit[, .N, by = .(ploidy, diagnosis_label, bidirectional_status, exact_rule_overlap)]
setorder(summary_dt, ploidy, diagnosis_label, bidirectional_status)
fwrite(
  summary_dt,
  file.path(outdir, "all_2group_bidirectional_rule_summary.tsv"),
  sep = "\t",
  quote = FALSE
)

problem_dt <- all_audit[exact_rule_overlap == "Current_Exact_but_not_bidirectional"]
fwrite(
  problem_dt,
  file.path(outdir, "current_exact_but_not_bidirectional.tsv"),
  sep = "\t",
  quote = FALSE
)

msg("Wrote audit table: %s", file.path(outdir, "all_2group_bidirectional_rule_audit.tsv"))
msg("Wrote summary: %s", file.path(outdir, "all_2group_bidirectional_rule_summary.tsv"))
msg("Current Exact_match but not bidirectional: %d", nrow(problem_dt))
