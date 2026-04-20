#!/usr/bin/env Rscript

source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "ibs_common.R"))

opt <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("matrix", "map", "outdir", "prefix")
missing <- required[!required %in% names(opt)]
if (length(missing) > 0) stop("Missing arguments: ", paste(missing, collapse = ", "))

opt[["anchor-col"]] <- opt[["anchor-col"]] %||% ""
opt[["cluster-threshold"]] <- opt[["cluster-threshold"]] %||% "0.99"
opt[["zmin"]] <- opt[["zmin"]] %||% "0.7"
opt[["zmax"]] <- opt[["zmax"]] %||% "1.0"

dir.create(opt[["outdir"]], recursive = TRUE, showWarnings = FALSE)
map_df <- read_sample_map(opt[["map"]])
anchor_col <- if (nzchar(opt[["anchor-col"]])) opt[["anchor-col"]] else names(map_df)[1]
if (!(anchor_col %in% names(map_df))) stop("Anchor column not found in map: ", anchor_col)

mat_df <- read.table(opt[["matrix"]], header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
ibs_mat <- as.matrix(mat_df[, -1, drop = FALSE]); mode(ibs_mat) <- "numeric"
rownames(ibs_mat) <- mat_df[[1]]
colnames(ibs_mat) <- colnames(mat_df)[-1]

anchor_ids <- unique(standardize_id(map_df[[anchor_col]]))
anchor_ids <- anchor_ids[!is.na(anchor_ids) & anchor_ids %in% rownames(ibs_mat)]
if (length(anchor_ids) == 0) stop("No anchor-group samples matched IBS matrix for column: ", anchor_col)

cluster_df <- build_similarity_clusters(anchor_ids, ibs_mat, threshold = as.numeric(opt[["cluster-threshold"]]))
cluster_df$anchor_group <- anchor_col
write.table(cluster_df, file = file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_cluster_table.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

ordered_anchor <- unique(c(
  cluster_df$sample_id[order(cluster_df$cluster_id, cluster_df$sample_id)],
  setdiff(anchor_ids, cluster_df$sample_id)
))
ordered_anchor <- ordered_anchor[ordered_anchor %in% anchor_ids]
full_mat <- ibs_mat[ordered_anchor, ordered_anchor, drop = FALSE]
draw_heatmap(full_mat, file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_anchor_cluster_heatmap.pdf")), paste(opt[["prefix"]], anchor_col, "cluster heatmap"), zlim = c(as.numeric(opt[["zmin"]]), as.numeric(opt[["zmax"]])))
draw_heatmap(full_mat, file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_anchor_full_heatmap.pdf")), paste(opt[["prefix"]], anchor_col, "full within-group heatmap"), zlim = c(as.numeric(opt[["zmin"]]), as.numeric(opt[["zmax"]])))

dup_ids <- unique(cluster_df$sample_id[cluster_df$cluster_size > 1])
if (length(dup_ids) >= 2) {
  dup_mat <- ibs_mat[dup_ids, dup_ids, drop = FALSE]
  draw_heatmap(dup_mat, file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_anchor_duplicate_heatmap.pdf")), paste(opt[["prefix"]], anchor_col, "duplicate-only heatmap"), zlim = c(as.numeric(opt[["zmin"]]), as.numeric(opt[["zmax"]])))
}

cat("Built anchor clusters for", anchor_col, "with", length(anchor_ids), "samples\n")
