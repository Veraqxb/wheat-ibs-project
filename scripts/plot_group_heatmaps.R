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
zmin <- as.numeric(args[["zmin"]] %||% "0.0")
zmax <- as.numeric(args[["zmax"]] %||% "1.0")
threshold <- as.numeric(args[["threshold"]] %||% "0.9")

dir.create(args[["outdir"]], recursive = TRUE, showWarnings = FALSE)

map_df <- read_sample_map(args[["map"]])
reference_mat <- read_mibs_matrix(args[["reference-mibs"]], args[["reference-id"]])
query_mat <- read_mibs_matrix(args[["query-mibs"]], args[["query-id"]])

threshold_palette <- function(zlim, threshold = 0.9) {
  threshold <- min(max(threshold, zlim[1]), zlim[2])
  lower_breaks <- seq(zlim[1], threshold, length.out = 31)
  upper_breaks <- seq(threshold, zlim[2], length.out = 121)
  breaks <- c(lower_breaks, upper_breaks[-1])
  lower_cols <- colorRampPalette(c("#F7F7F7", "#FFFFE5", "#FFF7BC"))(length(lower_breaks) - 1)
  upper_cols <- colorRampPalette(c("#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026", "#7F0000"))(length(upper_breaks) - 1)
  list(colors = c(lower_cols, upper_cols), breaks = breaks)
}

draw_threshold_heatmap <- function(sub_mat, file, title, xlab = "", ylab = "", zlim = c(0.0, 1.0), threshold = 0.9, mark_diagonal = FALSE, na_col = "#000000") {
  if (is.null(sub_mat) || nrow(sub_mat) == 0 || ncol(sub_mat) == 0) return(invisible(NULL))

  pal <- threshold_palette(zlim, threshold)
  nx <- ncol(sub_mat)
  ny <- nrow(sub_mat)
  pdf(file, width = max(8, nx * 0.5 + 3), height = max(8, ny * 0.5 + 3))
  par(mar = c(12, 12, 4, 2))
  image_mat <- t(sub_mat)[, ny:1, drop = FALSE]
  image_fill <- image_mat
  image_fill[is.na(image_fill)] <- zlim[1]
  image(
    x = seq_len(nx),
    y = seq_len(ny),
    z = image_fill,
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
  na_idx <- which(is.na(image_mat), arr.ind = TRUE)
  if (nrow(na_idx) > 0) {
    rect(
      xleft = na_idx[, 1] - 0.5,
      ybottom = na_idx[, 2] - 0.5,
      xright = na_idx[, 1] + 0.5,
      ytop = na_idx[, 2] + 0.5,
      col = na_col,
      border = NA
    )
  }
  if (mark_diagonal) {
    diag_n <- min(nx, ny)
    rect(
      xleft = seq_len(diag_n) - 0.5,
      ybottom = ny - seq_len(diag_n) + 0.5,
      xright = seq_len(diag_n) + 0.5,
      ytop = ny - seq_len(diag_n) + 1.5,
      border = "#252525",
      lwd = 1.4
    )
  }
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

build_display_slots <- function(raw_ids, valid_ids) {
  raw_ids <- standardize_id(raw_ids)
  if (length(raw_ids) == 0) {
    return(data.frame(actual_id = character(0), display_id = character(0), stringsAsFactors = FALSE))
  }

  display <- character(length(raw_ids))
  seen <- list()
  na_count <- 0L
  for (i in seq_along(raw_ids)) {
    cur <- raw_ids[[i]]
    if (is.na(cur) || !(cur %in% valid_ids)) {
      na_count <- na_count + 1L
      display[[i]] <- sprintf("NA_%03d", na_count)
    } else {
      seen[[cur]] <- (seen[[cur]] %||% 0L) + 1L
      display[[i]] <- if (seen[[cur]] > 1L) sprintf("%s#%d", cur, seen[[cur]]) else cur
    }
  }
  data.frame(actual_id = raw_ids, display_id = display, stringsAsFactors = FALSE)
}

build_internal_map_matrix <- function(map_df, group_col, ibs_mat) {
  if (!(group_col %in% names(map_df))) {
    return(NULL)
  }

  slots <- build_display_slots(map_df[[group_col]], rownames(ibs_mat))
  if (nrow(slots) == 0) return(NULL)

  out <- matrix(NA_real_, nrow = nrow(slots), ncol = nrow(slots))
  rownames(out) <- slots$display_id
  colnames(out) <- slots$display_id

  for (i in seq_len(nrow(slots))) {
    for (j in seq_len(nrow(slots))) {
      a <- slots$actual_id[[i]]
      b <- slots$actual_id[[j]]
      if (!is.na(a) && !is.na(b) && a %in% rownames(ibs_mat) && b %in% colnames(ibs_mat)) {
        out[i, j] <- ibs_mat[a, b]
      }
    }
  }
  out
}

build_ordered_cross_matrix <- function(map_df, x_col, y_col, ibs_mat) {
  if (!(x_col %in% names(map_df)) || !(y_col %in% names(map_df))) {
    return(NULL)
  }
  x_slots <- build_display_slots(map_df[[x_col]], colnames(ibs_mat))
  y_slots <- build_display_slots(map_df[[y_col]], rownames(ibs_mat))
  if (nrow(x_slots) == 0 || nrow(y_slots) == 0) return(NULL)

  out <- matrix(NA_real_, nrow = nrow(y_slots), ncol = nrow(x_slots))
  rownames(out) <- y_slots$display_id
  colnames(out) <- x_slots$display_id

  for (i in seq_len(nrow(y_slots))) {
    for (j in seq_len(nrow(x_slots))) {
      y_id <- y_slots$actual_id[[i]]
      x_id <- x_slots$actual_id[[j]]
      if (!is.na(y_id) && !is.na(x_id) && y_id %in% rownames(ibs_mat) && x_id %in% colnames(ibs_mat)) {
        out[i, j] <- ibs_mat[y_id, x_id]
      }
    }
  }
  out
}

draw_internal_if_available <- function(mat, map_df, group_col, prefix_root, dataset_label) {
  sub_mat <- build_internal_map_matrix(map_df, group_col, mat)
  if (is.null(sub_mat) || nrow(sub_mat) == 0) return(invisible(NULL))
  draw_threshold_heatmap(
    sub_mat,
    file.path(args[["outdir"]], paste0(prefix_root, "_", dataset_label, "_", group_col, "_internal_heatmap.pdf")),
    paste(args[["prefix"]], dataset_label, group_col, "internal IBS"),
    zlim = c(zmin, zmax),
    threshold = threshold
  )
}

draw_pair_if_available <- function(mat, map_df, x_col, y_col, prefix_root, dataset_label) {
  pair_mat <- build_ordered_cross_matrix(map_df, x_col, y_col, mat)
  if (is.null(pair_mat)) return(invisible(NULL))
  draw_threshold_heatmap(
    pair_mat,
    file.path(args[["outdir"]], paste0(prefix_root, "_", dataset_label, "_", y_col, "_vs_", x_col, "_pair_heatmap.pdf")),
    paste(args[["prefix"]], dataset_label, y_col, "vs", x_col, "expected 1-to-1 IBS"),
    xlab = x_col,
    ylab = y_col,
    zlim = c(zmin, zmax),
    threshold = threshold,
    mark_diagonal = TRUE
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
