#!/usr/bin/env Rscript

source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "ibs_common.R"))

opt <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("matrix", "map", "cluster-table", "pairwise-dir", "outdir", "prefix")
missing <- required[!required %in% names(opt)]
if (length(missing) > 0) stop("Missing arguments: ", paste(missing, collapse = ", "))

opt[["anchor-col"]] <- opt[["anchor-col"]] %||% ""
opt[["secondary-col"]] <- opt[["secondary-col"]] %||% "B25"
opt[["rna-groups"]] <- opt[["rna-groups"]] %||% "TC SC FC"
opt[["dna-threshold"]] <- opt[["dna-threshold"]] %||% "0.99"
opt[["rna-threshold"]] <- opt[["rna-threshold"]] %||% "0.90"
opt[["detection-mode"]] <- opt[["detection-mode"]] %||% "full"

dir.create(opt[["outdir"]], recursive = TRUE, showWarnings = FALSE)
map_df <- read_sample_map(opt[["map"]])
anchor_col <- if (nzchar(opt[["anchor-col"]])) opt[["anchor-col"]] else names(map_df)[1]
secondary_col <- opt[["secondary-col"]]
rna_groups <- unique(strsplit(opt[["rna-groups"]], "[[:space:]]+")[[1]])
rna_groups <- rna_groups[nzchar(rna_groups)]

mat_df <- read.table(opt[["matrix"]], header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
ibs_mat <- as.matrix(mat_df[, -1, drop = FALSE]); mode(ibs_mat) <- "numeric"
rownames(ibs_mat) <- mat_df[[1]]
colnames(ibs_mat) <- colnames(mat_df)[-1]

cluster_df <- read.table(opt[["cluster-table"]], header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
if (!("sample_id" %in% names(cluster_df))) stop("Cluster table missing sample_id column")
if (!("reference_group" %in% names(cluster_df))) cluster_df$reference_group <- anchor_col

pair_files <- list.files(opt[["pairwise-dir"]], pattern = "_pairwise\\.tsv$", full.names = TRUE)
pair_tables <- lapply(pair_files, function(x) read.table(x, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE))
names(pair_tables) <- sub("^.*_([^_/]+)_vs_.*_pairwise\\.tsv$", "\\1", basename(pair_files))

if (!(secondary_col %in% names(pair_tables))) {
  warning("Secondary reference group pairwise table not found for ", secondary_col, "; RNA rescue will use direct IBS lookup only.")
}

cluster_lookup <- split(cluster_df, paste(cluster_df$reference_group, cluster_df$sample_id, sep = "||"))
cluster_members_for <- function(reference_group, sample_id) {
  hit <- cluster_lookup[[paste(reference_group, sample_id, sep = "||")]]
  if (is.null(hit) || nrow(hit) == 0) return(character(0))
  members <- unique(unlist(strsplit(hit$cluster_members[1], ";", fixed = TRUE)))
  members[nzchar(members)]
}

compute_secondary_match <- function(query_id, ref_ids, threshold) {
  out <- list(expected_ibs = NA_real_, best_match = NA_character_, best_ibs = NA_real_, match_type = "NO_DATA", match_info = "No valid IBS data")
  if (is.na(query_id) || !(query_id %in% rownames(ibs_mat))) {
    out$match_info <- "Query sample not found in IBS matrix"
    return(out)
  }
  ref_ids <- ref_ids[!is.na(ref_ids) & ref_ids %in% colnames(ibs_mat)]
  if (length(ref_ids) == 0) {
    out$match_info <- "No valid secondary-reference samples in IBS matrix"
    return(out)
  }
  vals <- as.numeric(ibs_mat[query_id, ref_ids, drop = TRUE])
  names(vals) <- ref_ids
  if (all(is.na(vals))) return(out)
  ord <- order(vals, decreasing = TRUE, na.last = TRUE)
  out$best_match <- names(vals)[ord[1]]
  out$best_ibs <- vals[ord[1]]
  if (!is.na(out$best_ibs) && out$best_ibs > threshold) {
    out$match_type <- "MATCH"
    out$match_info <- paste("Secondary best match above threshold:", out$best_match)
  } else {
    out$match_type <- "LOW_IBS"
    out$match_info <- "Secondary best match below threshold"
  }
  out
}

simple_rows <- list()
summary_rows <- list()

for (group_name in names(map_df)[-1]) {
  pair_df <- pair_tables[[group_name]]
  if (is.null(pair_df) || nrow(pair_df) == 0) next

  group_threshold <- pick_threshold(group_name, rna_groups = rna_groups, dna_threshold = as.numeric(opt[["dna-threshold"]]), rna_threshold = as.numeric(opt[["rna-threshold"]]))
  group_rows <- list()

  for (i in seq_len(nrow(pair_df))) {
    row <- pair_df[i, , drop = FALSE]
    sample_id <- if ("sample_id" %in% names(row)) row$sample_id else row$query_id
    match_type <- row$match_type
    match_info <- row$match_info
    final_ibs <- row$pair_ibs
    matched_id <- row$anchor_id

    if (group_name %in% rna_groups) {
      secondary_expected <- if (secondary_col %in% names(map_df)) standardize_id(map_df[[secondary_col]][row$row_index]) else NA_character_
      secondary_res <- compute_secondary_match(row$query_id, secondary_expected, group_threshold)

      if (identical(match_type, "MATCH")) {
        match_type <- "MATCH_Z23"
        match_info <- paste("Matched expected", anchor_col)
        matched_id <- row$anchor_id
        final_ibs <- row$pair_ibs
      } else if (opt[["detection-mode"]] != "simple" && identical(secondary_res$match_type, "MATCH")) {
        match_type <- "RESCUED_B25"
        match_info <- paste("Unmatched to", anchor_col, "but rescued by", secondary_col, secondary_res$best_match)
        matched_id <- secondary_res$best_match
        final_ibs <- secondary_res$best_ibs
      } else {
        cluster_members <- cluster_members_for(anchor_col, row$anchor_id)
        secondary_cluster_members <- if (opt[["detection-mode"]] != "simple" && !is.na(secondary_expected)) cluster_members_for(secondary_col, secondary_expected) else character(0)
        neighbor_candidates <- c(row$best_anchor, if (opt[["detection-mode"]] != "simple") secondary_res$best_match else NA_character_)
        rescue_hit <- intersect(cluster_members, neighbor_candidates)
        secondary_rescue_hit <- intersect(secondary_cluster_members, neighbor_candidates)
        if (length(rescue_hit) > 0) {
          match_type <- "RESCUED_CLUSTER"
          match_info <- paste("Rescued by anchor-cluster neighbor", rescue_hit[1])
          matched_id <- rescue_hit[1]
          if (!is.na(row$best_anchor_ibs) && row$best_anchor == rescue_hit[1]) {
            final_ibs <- row$best_anchor_ibs
          } else if (!is.na(secondary_res$best_ibs) && secondary_res$best_match == rescue_hit[1]) {
            final_ibs <- secondary_res$best_ibs
          }
        } else if (length(secondary_rescue_hit) > 0) {
          match_type <- "RESCUED_CLUSTER"
          match_info <- paste("Rescued by", secondary_col, "cluster neighbor", secondary_rescue_hit[1])
          matched_id <- secondary_rescue_hit[1]
          if (!is.na(secondary_res$best_ibs) && secondary_res$best_match == secondary_rescue_hit[1]) {
            final_ibs <- secondary_res$best_ibs
          } else if (!is.na(row$best_anchor_ibs) && row$best_anchor == secondary_rescue_hit[1]) {
            final_ibs <- row$best_anchor_ibs
          }
        } else if (identical(match_type, "SWAPPED")) {
          match_type <- "SWAPPED"
          match_info <- paste("Best anchor is", row$best_anchor, "instead of expected", row$anchor_id)
          matched_id <- row$best_anchor
          final_ibs <- row$best_anchor_ibs
        } else if (identical(match_type, "NO_DATA")) {
          match_type <- "NO_DATA"
          match_info <- "No usable IBS values against Z23/B25"
          matched_id <- NA_character_
          final_ibs <- NA_real_
        } else {
          match_type <- "DROP"
          match_info <- paste("Unmatched to", anchor_col, "and", secondary_col)
          matched_id <- row$best_anchor
          final_ibs <- row$best_anchor_ibs
        }
      }
    } else {
      if (identical(match_type, "MATCH")) {
        match_type <- "MATCH_Z23"
        match_info <- paste("Expected anchor matched", anchor_col)
      } else if (identical(match_type, "SWAPPED")) {
        match_type <- "SWAPPED"
        match_info <- paste("Best anchor is", row$best_anchor)
        matched_id <- row$best_anchor
        final_ibs <- row$best_anchor_ibs
      } else if (identical(match_type, "NO_DATA")) {
        match_type <- "NO_DATA"
      } else {
        match_type <- "DROP"
        match_info <- paste("Expected and best anchor below", sprintf("%.2f", group_threshold))
      }
    }

    group_rows[[i]] <- data.frame(
      sample_id = sample_id,
      group_name = group_name,
      matched_id = matched_id,
      ibs = final_ibs,
      match_type = match_type,
      match_info = match_info,
      stringsAsFactors = FALSE
    )
  }

  group_df <- do.call(rbind, group_rows)
  write.table(group_df, file = file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_", group_name, "_matching_simple.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  simple_rows[[group_name]] <- group_df
  summary_rows[[group_name]] <- as.data.frame(
    table(
      group_name = rep(group_name, nrow(group_df)),
      match_type = group_df$match_type
    ),
    stringsAsFactors = FALSE
  )
}

full_simple <- do.call(rbind, simple_rows)
full_simple$sample_id <- ifelse(is.na(full_simple$sample_id), "", full_simple$sample_id)
write.table(full_simple, file = file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_full_matching_table.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

summary_df <- do.call(rbind, summary_rows)
names(summary_df)[3] <- "count"
write.table(summary_df, file = file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_matching_summary.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

low_df <- full_simple[full_simple$match_type %in% c("SWAPPED", "DROP", "NO_DATA"), , drop = FALSE]
write.table(low_df, file = file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_low_ibs_sample_table.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

rna_mismatch_df <- full_simple[full_simple$group_name %in% rna_groups & full_simple$match_type == "DROP", , drop = FALSE]
write.table(rna_mismatch_df, file = file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_rna_complete_mismatch.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

draw_category_barplot(summary_df, file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_matching_summary_categories.pdf")), paste(opt[["prefix"]], "final matching categories"))

cat("Built simplified final matching outputs for", nrow(full_simple), "rows\n")
