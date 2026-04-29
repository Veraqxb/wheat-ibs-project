#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- if (length(script_arg) == 0 || is.na(script_arg)) getwd() else sub("^--file=", "", script_arg)
script_path <- gsub("~\\+~", " ", script_path, fixed = FALSE)
source(file.path(dirname(normalizePath(script_path)), "ibs_common.R"))

opt <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("matrix", "map", "outdir", "prefix")
missing <- required[!required %in% names(opt)]
if (length(missing) > 0) stop("Missing arguments: ", paste(missing, collapse = ", "))

opt[["anchor-col"]] <- opt[["anchor-col"]] %||% ""
opt[["secondary-col"]] <- opt[["secondary-col"]] %||% "B25"
opt[["threshold"]] <- opt[["threshold"]] %||% "0.99"

dir.create(opt[["outdir"]], recursive = TRUE, showWarnings = FALSE)

map_df <- read_sample_map(opt[["map"]])
anchor_col <- if (nzchar(opt[["anchor-col"]])) opt[["anchor-col"]] else names(map_df)[1]
secondary_col <- opt[["secondary-col"]]
threshold <- as.numeric(opt[["threshold"]])

if (!(anchor_col %in% names(map_df))) stop("Anchor column not found in map: ", anchor_col)
if (!(secondary_col %in% names(map_df))) stop("Secondary column not found in map: ", secondary_col)

mat_df <- read.table(opt[["matrix"]], header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
ibs_mat <- as.matrix(mat_df[, -1, drop = FALSE])
mode(ibs_mat) <- "numeric"
rownames(ibs_mat) <- mat_df[[1]]
colnames(ibs_mat) <- colnames(mat_df)[-1]

anchor_all <- standardize_id(map_df[[anchor_col]])
secondary_all <- standardize_id(map_df[[secondary_col]])

anchor_ref <- unique(anchor_all[!is.na(anchor_all) & anchor_all %in% colnames(ibs_mat)])
secondary_ref <- unique(secondary_all[!is.na(secondary_all) & secondary_all %in% rownames(ibs_mat)])
if (length(anchor_ref) == 0) stop("No valid anchor IDs found in IBS matrix for ", anchor_col)
if (length(secondary_ref) == 0) stop("No valid secondary IDs found in IBS matrix for ", secondary_col)

best_match <- function(query_id, candidate_ids) {
  out <- list(id = NA_character_, ibs = NA_real_, second_id = NA_character_, second_ibs = NA_real_)
  if (is.na(query_id) || !(query_id %in% rownames(ibs_mat))) return(out)
  candidate_ids <- candidate_ids[candidate_ids %in% colnames(ibs_mat)]
  if (length(candidate_ids) == 0) return(out)
  vals <- as.numeric(ibs_mat[query_id, candidate_ids, drop = TRUE])
  names(vals) <- candidate_ids
  if (all(is.na(vals))) return(out)
  ord <- order(vals, decreasing = TRUE, na.last = TRUE)
  out$id <- names(vals)[ord[1]]
  out$ibs <- vals[ord[1]]
  if (length(ord) >= 2) {
    out$second_id <- names(vals)[ord[2]]
    out$second_ibs <- vals[ord[2]]
  }
  out
}

lookup_ibs <- function(row_id, col_id) {
  if (is.na(row_id) || is.na(col_id)) return(NA_real_)
  if (!(row_id %in% rownames(ibs_mat)) || !(col_id %in% colnames(ibs_mat))) return(NA_real_)
  as.numeric(ibs_mat[row_id, col_id])
}

classify_reference <- function(expected_ibs, anchor_id, secondary_id, best_z23, best_z23_ibs, best_b25, best_b25_ibs) {
  if (is.na(anchor_id) || is.na(secondary_id) || is.na(expected_ibs)) return("No_data")

  expected_high <- expected_ibs >= threshold
  b25_to_z23_exact <- expected_high && !is.na(best_z23) && identical(best_z23, anchor_id) && !is.na(best_z23_ibs) && best_z23_ibs >= threshold
  z23_to_b25_exact <- expected_high && !is.na(best_b25) && identical(best_b25, secondary_id) && !is.na(best_b25_ibs) && best_b25_ibs >= threshold
  any_high_shift <- (!is.na(best_z23_ibs) && best_z23_ibs >= threshold) || (!is.na(best_b25_ibs) && best_b25_ibs >= threshold)

  if (b25_to_z23_exact && z23_to_b25_exact) return("Exact_match")
  if (b25_to_z23_exact) return("B25_to_Z23_only_exact")
  if (z23_to_b25_exact) return("Z23_to_B25_only_exact")
  if (expected_high) return("Expected_high_but_best_shifted")
  if (any_high_shift) return("True_mismatch")
  "low_ibs_match"
}

short_class <- function(x) {
  ifelse(
    x == "B25_to_Z23_only_exact", "B25_only_exact",
    ifelse(x == "Z23_to_B25_only_exact", "Z23_only_exact", x)
  )
}

reference_confidence <- function(x) {
  ifelse(
    x == "Exact_match", "trusted",
    ifelse(x %in% c("True_mismatch", "No_data"), "exclude", "review")
  )
}

rows <- lapply(seq_len(nrow(map_df)), function(i) {
  anchor_id <- anchor_all[i]
  secondary_id <- secondary_all[i]
  expected_ibs <- lookup_ibs(secondary_id, anchor_id)
  z23_res <- best_match(secondary_id, anchor_ref)
  b25_res <- best_match(anchor_id, secondary_ref)
  fine <- classify_reference(
    expected_ibs = expected_ibs,
    anchor_id = anchor_id,
    secondary_id = secondary_id,
    best_z23 = z23_res$id,
    best_z23_ibs = z23_res$ibs,
    best_b25 = b25_res$id,
    best_b25_ibs = b25_res$ibs
  )

  data.frame(
    row_index = i,
    anchor_group = anchor_col,
    secondary_group = secondary_col,
    anchor_id = anchor_id,
    secondary_id = secondary_id,
    expected_pair_ibs = expected_ibs,
    best_Z23 = z23_res$id,
    best_Z23_ibs = z23_res$ibs,
    second_Z23 = z23_res$second_id,
    second_Z23_ibs = z23_res$second_ibs,
    best_B25 = b25_res$id,
    best_B25_ibs = b25_res$ibs,
    second_B25 = b25_res$second_id,
    second_B25_ibs = b25_res$second_ibs,
    dna_diagnosis_fine = fine,
    dna_diagnosis_short = short_class(fine),
    reference_confidence = reference_confidence(fine),
    stringsAsFactors = FALSE
  )
})

context_df <- do.call(rbind, rows)

write.table(
  context_df,
  file = file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_reference_context.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

summary_df <- as.data.frame(
  table(
    dna_diagnosis_short = context_df$dna_diagnosis_short,
    reference_confidence = context_df$reference_confidence
  ),
  stringsAsFactors = FALSE
)
names(summary_df)[3] <- "count"
summary_df <- summary_df[summary_df$count > 0, , drop = FALSE]
write.table(
  summary_df,
  file = file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_reference_context_summary.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

attention_df <- context_df[context_df$reference_confidence != "trusted", , drop = FALSE]
write.table(
  attention_df,
  file = file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_reference_context_attention.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Built 2group reference context for", nrow(context_df), "map rows\n")
