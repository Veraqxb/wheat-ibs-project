#!/usr/bin/env Rscript

source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "ibs_common.R"))

opt <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("matrix", "map", "outdir", "prefix")
missing <- required[!required %in% names(opt)]
if (length(missing) > 0) stop("Missing arguments: ", paste(missing, collapse = ", "))

opt[["anchor-col"]] <- opt[["anchor-col"]] %||% ""
opt[["rna-groups"]] <- opt[["rna-groups"]] %||% "TC SC FC"
opt[["dna-threshold"]] <- opt[["dna-threshold"]] %||% "0.99"
opt[["rna-threshold"]] <- opt[["rna-threshold"]] %||% "0.90"
opt[["zmin"]] <- opt[["zmin"]] %||% "0.7"
opt[["zmax"]] <- opt[["zmax"]] %||% "1.0"

dir.create(opt[["outdir"]], recursive = TRUE, showWarnings = FALSE)
map_df <- read_sample_map(opt[["map"]])
anchor_col <- if (nzchar(opt[["anchor-col"]])) opt[["anchor-col"]] else names(map_df)[1]
rna_groups <- unique(strsplit(opt[["rna-groups"]], "[[:space:]]+")[[1]])
rna_groups <- rna_groups[nzchar(rna_groups)]

mat_df <- read.table(opt[["matrix"]], header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
ibs_mat <- as.matrix(mat_df[, -1, drop = FALSE]); mode(ibs_mat) <- "numeric"
rownames(ibs_mat) <- mat_df[[1]]
colnames(ibs_mat) <- colnames(mat_df)[-1]

anchor_ids_all <- standardize_id(map_df[[anchor_col]])
anchor_ids <- anchor_ids_all[!is.na(anchor_ids_all) & anchor_ids_all %in% rownames(ibs_mat)]
if (length(anchor_ids) == 0) stop("No anchor-group samples matched IBS matrix for column: ", anchor_col)

summary_list <- list()
all_low <- list()

for (group_name in names(map_df)[-1]) {
  group_ids_all <- standardize_id(map_df[[group_name]])
  sample_idx <- seq_len(nrow(map_df))
  pair_df <- data.frame(
    row_index = sample_idx,
    sample_id = group_ids_all,
    group_name = group_name,
    anchor_id = anchor_ids_all,
    query_id = group_ids_all,
    stringsAsFactors = FALSE
  )

  group_threshold <- pick_threshold(group_name, rna_groups = rna_groups, dna_threshold = as.numeric(opt[["dna-threshold"]]), rna_threshold = as.numeric(opt[["rna-threshold"]]))
  pair_df$threshold <- group_threshold
  pair_df$pair_ibs <- NA_real_
  pair_df$best_anchor <- NA_character_
  pair_df$best_anchor_ibs <- NA_real_
  pair_df$second_anchor <- NA_character_
  pair_df$second_anchor_ibs <- NA_real_
  pair_df$third_anchor <- NA_character_
  pair_df$third_anchor_ibs <- NA_real_
  pair_df$match_type <- "NO_DATA"
  pair_df$match_info <- "No valid IBS data"

  for (i in seq_len(nrow(pair_df))) {
    qid <- pair_df$query_id[i]
    aid <- pair_df$anchor_id[i]
    if (is.na(qid) || !(qid %in% rownames(ibs_mat))) {
      pair_df$match_info[i] <- "Query sample not found in IBS matrix"
      next
    }

    vals <- as.numeric(ibs_mat[qid, anchor_ids, drop = TRUE])
    names(vals) <- anchor_ids
    pair_df$pair_ibs[i] <- if (!is.na(aid) && aid %in% anchor_ids) ibs_mat[qid, aid] else NA_real_

    if (all(is.na(vals))) {
      pair_df$match_info[i] <- "All anchor IBS values are NA"
      next
    }

    ord <- order(vals, decreasing = TRUE, na.last = TRUE)
    if (length(ord) >= 1) {
      pair_df$best_anchor[i] <- names(vals)[ord[1]]
      pair_df$best_anchor_ibs[i] <- vals[ord[1]]
    }
    if (length(ord) >= 2) {
      pair_df$second_anchor[i] <- names(vals)[ord[2]]
      pair_df$second_anchor_ibs[i] <- vals[ord[2]]
    }
    if (length(ord) >= 3) {
      pair_df$third_anchor[i] <- names(vals)[ord[3]]
      pair_df$third_anchor_ibs[i] <- vals[ord[3]]
    }

    if (!is.na(pair_df$pair_ibs[i]) && pair_df$pair_ibs[i] > group_threshold) {
      pair_df$match_type[i] <- "MATCH"
      pair_df$match_info[i] <- paste("Expected anchor passes threshold:", aid)
    } else if (!is.na(pair_df$best_anchor_ibs[i]) && pair_df$best_anchor_ibs[i] > group_threshold) {
      pair_df$match_type[i] <- "SWAPPED"
      pair_df$match_info[i] <- paste("Best anchor above threshold:", pair_df$best_anchor[i])
    } else if (!is.na(pair_df$pair_ibs[i]) || !is.na(pair_df$best_anchor_ibs[i])) {
      pair_df$match_type[i] <- "LOW_IBS"
      pair_df$match_info[i] <- "Expected and best anchor below threshold"
    }
  }

  pair_out <- file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_", group_name, "_vs_", anchor_col, "_pairwise.tsv"))
  write.table(pair_df, file = pair_out, sep = "\t", quote = FALSE, row.names = FALSE)

  valid_query_ids <- unique(pair_df$query_id[!is.na(pair_df$query_id) & pair_df$query_id %in% rownames(ibs_mat)])
  if (length(valid_query_ids) >= 2) {
    within_mat <- ibs_mat[valid_query_ids, valid_query_ids, drop = FALSE]
    draw_heatmap(within_mat, file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_", group_name, "_within_heatmap.pdf")), paste(opt[["prefix"]], group_name, "within-group heatmap"), zlim = c(as.numeric(opt[["zmin"]]), as.numeric(opt[["zmax"]])))
    dup_df <- build_similarity_clusters(valid_query_ids, ibs_mat, threshold = 0.99)
    dup_ids <- unique(dup_df$sample_id[dup_df$cluster_size > 1])
    if (length(dup_ids) >= 2) {
      dup_mat <- ibs_mat[dup_ids, dup_ids, drop = FALSE]
      draw_heatmap(dup_mat, file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_", group_name, "_duplicate_heatmap.pdf")), paste(opt[["prefix"]], group_name, "duplicate heatmap"), zlim = c(as.numeric(opt[["zmin"]]), as.numeric(opt[["zmax"]])))
    }
    internal_vals <- get_upper_triangle_values(within_mat)
  } else {
    internal_vals <- numeric(0)
  }

  draw_density_plot(
    internal_vals = internal_vals,
    pair_vals = pair_df$pair_ibs,
    out_file = file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_", group_name, "_density.pdf")),
    title = paste(opt[["prefix"]], group_name, "internal vs expected-pair IBS"),
    internal_label = paste0(group_name, "_internal"),
    pair_label = "paired_1to1",
    threshold = group_threshold
  )

  group_summary <- as.data.frame(table(group_name = group_name, match_type = factor(pair_df$match_type, levels = c("MATCH", "SWAPPED", "LOW_IBS", "NO_DATA"))), stringsAsFactors = FALSE)
  names(group_summary)[3] <- "count"
  group_summary$threshold <- group_threshold
  summary_list[[group_name]] <- group_summary
  all_low[[group_name]] <- pair_df[pair_df$match_type %in% c("SWAPPED", "LOW_IBS", "NO_DATA"), , drop = FALSE]
}

summary_df <- do.call(rbind, summary_list)
write.table(summary_df, file = file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_group_matching_summary.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

low_df <- do.call(rbind, all_low)
if (!is.null(low_df) && nrow(low_df) > 0) {
  write.table(low_df, file = file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_low_ibs_samples.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
}

draw_category_barplot(summary_df, file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_group_matching_summary.pdf")), paste(opt[["prefix"]], "matching categories by group"))

cat("Completed pairwise matching for", length(names(map_df)) - 1, "groups against", anchor_col, "\n")
