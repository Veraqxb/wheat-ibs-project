#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
if (is.na(script_arg)) {
  stop("Cannot determine script path for plot_group_heatmaps.R")
}
script_dirname <- dirname(normalizePath(sub("^--file=", "", script_arg)))
source(file.path(script_dirname, "ibs_common.R"))

args <- parse_args(commandArgs(trailingOnly = TRUE))

required <- c("reference-mibs", "reference-id", "query-mibs", "query-id", "map", "outdir", "prefix")
missing_args <- required[!required %in% names(args)]
if (length(missing_args) > 0) {
  stop("Missing arguments: ", paste(missing_args, collapse = ", "))
}

anchor_col <- args[["anchor-col"]] %||% "Z23"
secondary_col <- args[["secondary-col"]] %||% "B25"
rna_groups <- strsplit(args[["rna-groups"]] %||% "TC FC SC", "[[:space:],]+")[[1]]
rna_groups <- rna_groups[nzchar(rna_groups)]
zmin <- as.numeric(args[["zmin"]] %||% "0.7")
zmax <- as.numeric(args[["zmax"]] %||% "1.0")
threshold <- as.numeric(args[["threshold"]] %||% "0.9")

dir.create(args[["outdir"]], recursive = TRUE, showWarnings = FALSE)

map_df <- read_sample_map(args[["map"]])
reference_mat <- read_mibs_matrix(args[["reference-mibs"]], args[["reference-id"]])
query_mat <- read_mibs_matrix(args[["query-mibs"]], args[["query-id"]])

threshold_palette <- function(zlim, threshold = 0.9) {
  threshold <- min(max(threshold, zlim[1]), zlim[2])
  lower_breaks <- seq(zlim[1], threshold, length.out = 31)
  upper_breaks <- seq(threshold, zlim[2], length.out = 91)
  breaks <- c(lower_breaks, upper_breaks[-1])
  lower_cols <- colorRampPalette(c("#F7F7F7", "#FFFFE5", "#FFF7BC"))(length(lower_breaks) - 1)
  upper_cols <- colorRampPalette(c("#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026"))(length(upper_breaks) - 1)
  list(colors = c(lower_cols, upper_cols), breaks = breaks)
}

draw_threshold_heatmap <- function(sub_mat, file, title, xlab = "", ylab = "", zlim = c(0.7, 1.0), threshold = 0.9) {
  if (is.null(sub_mat) || nrow(sub_mat) == 0 || ncol(sub_mat) == 0) return(invisible(NULL))

  pal <- threshold_palette(zlim, threshold)
  nx <- ncol(sub_mat)
  ny <- nrow(sub_mat)
  pdf(file, width = max(8, nx * 0.5 + 3), height = max(8, ny * 0.5 + 3))
  par(mar = c(12, 12, 4, 2))
  image_mat <- t(sub_mat)[, ny:1, drop = FALSE]
  image(
    x = seq_len(nx),
    y = seq_len(ny),
    z = image_mat,
    col = pal$colors,
    breaks = pal$breaks,
    axes = FALSE,
    main = title,
    xlab = xlab,
    ylab = ylab
  )
  axis(1, at = seq_len(nx), labels = colnames(sub_mat), las = 2, cex.axis = 0.7)
  axis(2, at = seq_len(ny), labels = rev(rownames(sub_mat)), las = 1, cex.axis = 0.7)
  abline(h = seq(0.5, ny + 0.5, by = 1), col = "grey82")
  abline(v = seq(0.5, nx + 0.5, by = 1), col = "grey82")
  box(col = "grey50")
  for (i in seq_len(nx)) {
    for (j in seq_len(ny)) {
      val <- sub_mat[j, i]
      if (is.na(val)) next
      plot_y <- ny - j + 1
      text_col <- if (val >= 0.97) "white" else if (val >= threshold) "black" else "#4D4D4D"
      text(i, plot_y, sprintf("%.3f", val), cex = 0.55, col = text_col, font = ifelse(val >= 0.99, 2, 1))
    }
  }
  dev.off()
}

get_group_ids <- function(map_df, col_name, ibs_mat) {
  if (!(col_name %in% names(map_df))) return(character(0))
  ids <- standardize_id(map_df[[col_name]])
  ids <- ids[!is.na(ids) & ids %in% rownames(ibs_mat)]
  unique(ids)
}

build_expected_pair_matrix <- function(map_df, x_col, y_col, ibs_mat) {
  if (!(x_col %in% names(map_df)) || !(y_col %in% names(map_df))) {
    return(NULL)
  }
  x_ids <- standardize_id(map_df[[x_col]])
  y_ids <- standardize_id(map_df[[y_col]])
  keep <- !is.na(x_ids) & !is.na(y_ids) & x_ids %in% colnames(ibs_mat) & y_ids %in% rownames(ibs_mat)
  if (!any(keep)) return(NULL)

  x_ids <- x_ids[keep]
  y_ids <- y_ids[keep]
  pair_mat <- matrix(NA_real_, nrow = length(y_ids), ncol = length(x_ids))
  for (i in seq_along(y_ids)) {
    pair_mat[i, i] <- ibs_mat[y_ids[i], x_ids[i]]
  }
  rownames(pair_mat) <- y_ids
  colnames(pair_mat) <- x_ids
  pair_mat
}

draw_internal_if_available <- function(mat, map_df, group_col, prefix_root, dataset_label) {
  ids <- get_group_ids(map_df, group_col, mat)
  if (length(ids) < 2) return(invisible(NULL))
  sub_mat <- mat[ids, ids, drop = FALSE]
  draw_threshold_heatmap(
    sub_mat,
    file.path(args[["outdir"]], paste0(prefix_root, "_", dataset_label, "_", group_col, "_internal_heatmap.pdf")),
    paste(args[["prefix"]], dataset_label, group_col, "internal IBS"),
    zlim = c(zmin, zmax),
    threshold = threshold
  )
}

draw_pair_if_available <- function(mat, map_df, x_col, y_col, prefix_root, dataset_label) {
  pair_mat <- build_expected_pair_matrix(map_df, x_col, y_col, mat)
  if (is.null(pair_mat)) return(invisible(NULL))
  draw_threshold_heatmap(
    pair_mat,
    file.path(args[["outdir"]], paste0(prefix_root, "_", dataset_label, "_", y_col, "_vs_", x_col, "_pair_heatmap.pdf")),
    paste(args[["prefix"]], dataset_label, y_col, "vs", x_col, "expected 1-to-1 IBS"),
    xlab = x_col,
    ylab = y_col,
    zlim = c(zmin, zmax),
    threshold = threshold
  )
}

draw_internal_if_available(reference_mat, map_df, anchor_col, args[["prefix"]], "reference2group")
draw_internal_if_available(reference_mat, map_df, secondary_col, args[["prefix"]], "reference2group")
draw_pair_if_available(reference_mat, map_df, anchor_col, secondary_col, args[["prefix"]], "reference2group")

draw_internal_if_available(query_mat, map_df, anchor_col, args[["prefix"]], "query5group")
draw_internal_if_available(query_mat, map_df, secondary_col, args[["prefix"]], "query5group")

for (group_name in rna_groups) {
  draw_pair_if_available(query_mat, map_df, anchor_col, group_name, args[["prefix"]], "query5group")
  draw_pair_if_available(query_mat, map_df, secondary_col, group_name, args[["prefix"]], "query5group")
}

cat("Completed standalone heatmap export for", args[["prefix"]], "\n")
