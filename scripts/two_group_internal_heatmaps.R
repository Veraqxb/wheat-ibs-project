#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(pheatmap)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
if (is.na(script_arg)) stop("Cannot determine script path for two_group_internal_heatmaps.R")
script_path <- sub("^--file=", "", script_arg)
# Rscript can encode spaces in --file paths as "~+~" on some systems.
script_path <- gsub("~\\+~", " ", script_path, fixed = FALSE)
script_dirname <- dirname(normalizePath(script_path))
source(file.path(script_dirname, "ibs_common.R"))

opt <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("mibs", "id", "map", "outdir", "prefix")
missing <- required[!required %in% names(opt)]
if (length(missing) > 0) stop("Missing arguments: ", paste(missing, collapse = ", "))

anchor_col <- opt[["anchor-col"]] %||% "Z23"
secondary_col <- opt[["secondary-col"]] %||% "B25"
cluster_threshold <- as.numeric(opt[["cluster-threshold"]] %||% "0.99")
zmin <- as.numeric(opt[["zmin"]] %||% "0.7")
zmax <- as.numeric(opt[["zmax"]] %||% "1.0")
info_file <- opt[["info"]] %||% NA_character_

dir.create(opt[["outdir"]], recursive = TRUE, showWarnings = FALSE)

map_df <- read_sample_map(opt[["map"]])
ibs_mat <- read_mibs_matrix(opt[["mibs"]], opt[["id"]])
info_df <- NULL
if (!is.na(info_file) && file.exists(info_file)) {
  info_df <- read.table(info_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, fill = TRUE)
  for (nm in names(info_df)) info_df[[nm]] <- standardize_id(info_df[[nm]])
}

if (!(anchor_col %in% names(map_df))) stop("Missing anchor column in map: ", anchor_col)
if (!(secondary_col %in% names(map_df))) stop("Missing secondary column in map: ", secondary_col)

build_display_slots <- function(raw_ids, valid_ids = NULL, missing_prefix = "NA") {
  raw_ids <- standardize_id(raw_ids)
  labels <- character(length(raw_ids))
  resolved_ids <- rep(NA_character_, length(raw_ids))
  seen <- list()
  na_count <- 0L

  for (i in seq_along(raw_ids)) {
    id <- raw_ids[i]
    valid <- !is.na(id) && id != "" && (is.null(valid_ids) || id %in% valid_ids)
    if (valid) {
      seen[[id]] <- (seen[[id]] %||% 0L) + 1L
      suffix <- if (seen[[id]] > 1L) paste0("#", seen[[id]]) else ""
      labels[i] <- paste0(id, suffix)
      resolved_ids[i] <- id
    } else {
      na_count <- na_count + 1L
      labels[i] <- sprintf("%s_%03d", missing_prefix, na_count)
      resolved_ids[i] <- NA_character_
    }
  }

  data.frame(
    slot_index = seq_along(raw_ids),
    raw_id = raw_ids,
    resolved_id = resolved_ids,
    display_id = labels,
    stringsAsFactors = FALSE
  )
}

build_internal_map_matrix <- function(map_df, col_name, ibs_mat) {
  slots <- build_display_slots(map_df[[col_name]], valid_ids = rownames(ibs_mat), missing_prefix = paste0(col_name, "_NA"))
  n <- nrow(slots)
  out <- matrix(NA_real_, nrow = n, ncol = n)
  rownames(out) <- slots$display_id
  colnames(out) <- slots$display_id

  valid_idx <- which(!is.na(slots$resolved_id))
  for (i in valid_idx) {
    for (j in valid_idx) {
      out[i, j] <- ibs_mat[slots$resolved_id[i], slots$resolved_id[j]]
    }
  }
  # Internal self-comparisons for real samples should be displayed as IBS=1.
  # Missing map slots remain NA, so they are still shown with the NA color.
  if (length(valid_idx) > 0) {
    out[cbind(valid_idx, valid_idx)] <- 1
  }
  list(matrix = out, slots = slots)
}

build_cross_display_matrix <- function(map_df, anchor_col, secondary_col, ibs_mat) {
  anchor_slots <- build_display_slots(map_df[[anchor_col]], valid_ids = colnames(ibs_mat), missing_prefix = paste0(anchor_col, "_NA"))
  secondary_slots <- build_display_slots(map_df[[secondary_col]], valid_ids = rownames(ibs_mat), missing_prefix = paste0(secondary_col, "_NA"))

  out <- matrix(NA_real_, nrow = nrow(secondary_slots), ncol = nrow(anchor_slots))
  rownames(out) <- secondary_slots$display_id
  colnames(out) <- anchor_slots$display_id

  valid_rows <- which(!is.na(secondary_slots$resolved_id))
  valid_cols <- which(!is.na(anchor_slots$resolved_id))
  for (i in valid_rows) {
    for (j in valid_cols) {
      out[i, j] <- ibs_mat[secondary_slots$resolved_id[i], anchor_slots$resolved_id[j]]
    }
  }

  list(matrix = out, anchor_slots = anchor_slots, secondary_slots = secondary_slots)
}

build_display_number_matrix <- function(mat) {
  nm <- matrix("", nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  nm[!is.na(mat)] <- sprintf("%.3f", mat[!is.na(mat)])
  nm
}

heatmap_palette <- function(zmin, zmax) {
  colors <- c(
    colorRampPalette(c("#2166AC", "#67A9CF"))(40),  # <0.90 blue
    colorRampPalette(c("#D9F0A3", "#A6D96A"))(20),  # 0.90-0.95 light green
    colorRampPalette(c("#FFF7BC", "#FEE391"))(20),  # 0.95-0.99 light yellow
    colorRampPalette(c("#FB6A4A", "#CB181D"))(20)   # >=0.99 red
  )
  # pheatmap requires length(breaks) == length(colors) + 1.
  # Drop the first break of each later segment to avoid over-indexing
  # high IBS values such as 1.0 into the NA color.
  breaks <- c(
    seq(zmin, 0.90, length.out = 41),
    seq(0.900001, 0.95, length.out = 21)[-1],
    seq(0.950001, 0.99, length.out = 21)[-1],
    seq(0.990001, zmax + 1e-06, length.out = 21)[-1]
  )
  list(colors = colors, breaks = breaks, na_col = "#4D4D4D")
}

build_group_annotation <- function(sample_ids, info_dt, group_name) {
  if (is.null(info_dt) || nrow(info_dt) == 0) return(NULL)
  if (!("GROUP" %in% names(info_dt))) return(NULL)

  key_col <- if (group_name == "Z23") {
    if ("CAMP编号" %in% names(info_dt)) "CAMP编号" else if (group_name %in% names(info_dt)) group_name else NULL
  } else if (group_name == "B25") {
    if ("B25" %in% names(info_dt)) "B25" else if (group_name %in% names(info_dt)) group_name else NULL
  } else if (group_name %in% names(info_dt)) {
    group_name
  } else {
    NULL
  }

  if (is.null(key_col)) return(NULL)
  tmp <- unique(info_dt[, c(key_col, "GROUP"), drop = FALSE])
  colnames(tmp) <- c("sample_id", "GROUP")
  tmp$sample_id <- standardize_id(tmp$sample_id)
  tmp$GROUP <- standardize_id(tmp$GROUP)
  tmp <- tmp[!is.na(tmp$sample_id) & tmp$sample_id != "", , drop = FALSE]
  # camp_info may contain repeated sample rows; keep the first annotation per ID
  # so pheatmap receives unique row names.
  tmp <- tmp[!duplicated(tmp$sample_id), , drop = FALSE]
  ann <- data.frame(sample_id = unique(sample_ids), stringsAsFactors = FALSE)
  ann <- merge(ann, tmp, by = "sample_id", all.x = TRUE, sort = FALSE)
  ann$GROUP[is.na(ann$GROUP) | ann$GROUP == ""] <- "Unknown"
  rownames(ann) <- ann$sample_id
  data.frame(GROUP = ann$GROUP, row.names = ann$sample_id, stringsAsFactors = FALSE)
}

make_group_colors <- function(groups) {
  groups <- unique(groups)
  palette_pool <- c(
    "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
    "#56B4E9", "#F0E442", "#999999", "#332288", "#88CCEE",
    "#44AA99", "#117733", "#DDCC77", "#CC6677", "#882255", "#AA4499"
  )
  if (length(groups) > length(palette_pool)) {
    palette_pool <- grDevices::colorRampPalette(palette_pool)(length(groups))
  }
  cols <- palette_pool[seq_along(groups)]
  names(cols) <- groups
  list(GROUP = cols)
}

build_slot_annotation <- function(slots, info_df, group_name) {
  ann_valid <- build_group_annotation(slots$resolved_id[!is.na(slots$resolved_id)], info_df, group_name)
  ann_full <- data.frame(
    GROUP = rep("Missing", nrow(slots)),
    row.names = slots$display_id,
    stringsAsFactors = FALSE
  )
  if (!is.null(ann_valid)) {
    valid_ann_idx <- match(rownames(ann_valid), slots$resolved_id)
    valid_ann_idx <- valid_ann_idx[!is.na(valid_ann_idx)]
    ann_full$GROUP[valid_ann_idx] <- ann_valid$GROUP
  }
  ann_full
}

map_order <- data.frame(
  row_index = seq_len(nrow(map_df)),
  anchor_id = standardize_id(map_df[[anchor_col]]),
  secondary_id = standardize_id(map_df[[secondary_col]]),
  stringsAsFactors = FALSE
)
write.table(map_order, file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_map_order.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

group_ids <- function(col_name) {
  ids <- standardize_id(map_df[[col_name]])
  ids[!is.na(ids) & ids %in% rownames(ibs_mat)]
}

duplicate_pairs <- function(ids, group_name) {
  ids <- unique(ids)
  if (length(ids) < 2) {
    return(data.frame(group = character(0), sample_a = character(0), sample_b = character(0), ibs = numeric(0), stringsAsFactors = FALSE))
  }
  sub_mat <- ibs_mat[ids, ids, drop = FALSE]
  idx <- which(upper.tri(sub_mat) & !is.na(sub_mat) & sub_mat >= cluster_threshold, arr.ind = TRUE)
  if (nrow(idx) == 0) {
    return(data.frame(group = character(0), sample_a = character(0), sample_b = character(0), ibs = numeric(0), stringsAsFactors = FALSE))
  }
  data.frame(
    group = group_name,
    sample_a = rownames(sub_mat)[idx[, 1]],
    sample_b = colnames(sub_mat)[idx[, 2]],
    ibs = sub_mat[idx],
    stringsAsFactors = FALSE
  )
}

ordered_for_structure <- function(sub_mat) {
  # Temporary imputation is used only for display ordering. The original
  # diagnosis matrix is returned unchanged for plotting.
  if (nrow(sub_mat) < 3) return(rownames(sub_mat))
  order_mat <- sub_mat
  diag(order_mat) <- 1
  finite_vals <- order_mat[is.finite(order_mat)]
  fill_value <- if (length(finite_vals) > 0) median(finite_vals, na.rm = TRUE) else 0
  order_mat[is.na(order_mat)] <- fill_value
  dist_mat <- as.dist(1 - order_mat)
  hc <- hclust(dist_mat, method = "average")
  rownames(order_mat)[hc$order]
}

heat_cfg <- heatmap_palette(zmin, zmax)

pheatmap_with_consistent_style <- function(mat, file, main, annotation_row = NULL, annotation_col = NULL,
                                           annotation_colors = NULL, cluster_rows = FALSE, cluster_cols = FALSE,
                                           width_scale = 0.35, height_scale = 0.35,
                                           force_numbers = FALSE) {
  nmax <- max(nrow(mat), ncol(mat))
  if (force_numbers) {
    plot_width <- max(10, min(180, 4 + ncol(mat) * max(width_scale, 0.31)))
    plot_height <- max(10, min(180, 4 + nrow(mat) * max(height_scale, 0.31)))
  } else {
    plot_width <- max(8, min(72, 4 + ncol(mat) * width_scale))
    plot_height <- max(8, min(72, 4 + nrow(mat) * height_scale))
  }
  label_font <- if (nmax <= 80) {
    5
  } else if (nmax <= 220) {
    3
  } else if (nmax <= 650) {
    1.7
  } else {
    1.2
  }
  number_font <- if (nmax <= 30) {
    8
  } else if (nmax <= 80) {
    5
  } else if (nmax <= 220) {
    2.2
  } else if (nmax <= 650) {
    if (force_numbers) 2.6 else 1.2
  } else {
    if (force_numbers) 1.2 else 0.8
  }
  # The explicit *_with_IBS_values.pdf outputs must always contain IBS labels.
  # Non-value versions stay clean for structure-level viewing.
  show_numbers <- force_numbers || nmax <= 80
  disp_mat <- if (show_numbers) build_display_number_matrix(mat) else FALSE
  pdf(file, width = plot_width, height = plot_height)
  pheatmap(
    mat,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    annotation_row = annotation_row,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    display_numbers = disp_mat,
    fontsize = label_font,
    fontsize_row = label_font,
    fontsize_col = label_font,
    fontsize_number = number_font,
    angle_col = 90,
    color = heat_cfg$colors,
    breaks = heat_cfg$breaks,
    na_col = heat_cfg$na_col,
    border_color = NA,
    main = main
  )
  dev.off()
}

draw_group_outputs <- function(col_name) {
  ids <- unique(group_ids(col_name))
  if (length(ids) == 0) return(NULL)

  map_mat_info <- build_internal_map_matrix(map_df, col_name, ibs_mat)
  write.table(
    map_mat_info$slots,
    file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_", col_name, "_map_slots.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  ann_map <- build_group_annotation(map_mat_info$slots$resolved_id[!is.na(map_mat_info$slots$resolved_id)], info_df, col_name)
  ann_map_full <- NULL
  ann_colors <- NULL
  if (!is.null(ann_map)) {
    ann_map_full <- data.frame(GROUP = rep("Missing", nrow(map_mat_info$slots)), row.names = map_mat_info$slots$display_id, stringsAsFactors = FALSE)
    valid_ann_idx <- match(rownames(ann_map), map_mat_info$slots$resolved_id)
    valid_ann_idx <- valid_ann_idx[!is.na(valid_ann_idx)]
    ann_map_full$GROUP[valid_ann_idx] <- ann_map$GROUP
    ann_colors <- make_group_colors(c(ann_map_full$GROUP, "Missing"))
    if (!("Missing" %in% names(ann_colors$GROUP))) ann_colors$GROUP <- c(ann_colors$GROUP, Missing = "#000000")
  }

  if (nrow(map_mat_info$matrix) > 0) {
    pheatmap_with_consistent_style(
      map_mat_info$matrix,
      file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_", col_name, "_internal_map_order_heatmap.pdf")),
      paste(opt[["prefix"]], col_name, "internal heatmap (map order)"),
      annotation_row = ann_map_full,
      annotation_col = ann_map_full,
      annotation_colors = ann_colors
    )
  }

  cluster_df <- build_similarity_clusters(ids, ibs_mat, threshold = cluster_threshold)
  cluster_df$reference_group <- col_name
  write.table(cluster_df, file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_", col_name, "_cluster_table.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

  assignment_df <- cluster_df[, c("reference_group", "sample_id", "cluster_id", "cluster_size", "cluster_members", "neighbor_ids", "max_neighbor_ibs"), drop = FALSE]
  write.table(assignment_df, file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_", col_name, "_assignment_table.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

  sub_mat <- ibs_mat[ids, ids, drop = FALSE]
  cluster_order <- cluster_df$sample_id[order(ifelse(is.na(cluster_df$cluster_id), "singleton", cluster_df$cluster_id), cluster_df$sample_id)]
  cluster_order <- unique(cluster_order[cluster_order %in% rownames(sub_mat)])
  ann_mat <- build_group_annotation(cluster_order, info_df, col_name)
  ann_colors_cluster <- if (!is.null(ann_mat)) make_group_colors(ann_mat$GROUP) else NULL
  if (length(cluster_order) > 1) {
    draw_heatmap(sub_mat[cluster_order, cluster_order, drop = FALSE], file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_", col_name, "_threshold_cluster_heatmap.pdf")), paste(opt[["prefix"]], col_name, "IBS threshold-cluster heatmap"), zlim = c(zmin, zmax))
    if (!is.null(ann_mat)) {
      pheatmap_with_consistent_style(
        sub_mat[cluster_order, cluster_order, drop = FALSE],
        file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_", col_name, "_threshold_cluster_heatmap_with_GROUP.pdf")),
        paste(opt[["prefix"]], col_name, "IBS threshold-cluster heatmap + GROUP"),
        annotation_row = ann_mat,
        annotation_col = ann_mat,
        annotation_colors = ann_colors_cluster
      )
    }
  }

  h_order <- ordered_for_structure(sub_mat)
  if (length(h_order) > 1) {
    draw_heatmap(sub_mat[h_order, h_order, drop = FALSE], file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_", col_name, "_hierarchical_structure_heatmap.pdf")), paste(opt[["prefix"]], col_name, "hierarchical structure heatmap"), zlim = c(zmin, zmax))
    ann_h <- build_group_annotation(h_order, info_df, col_name)
    if (!is.null(ann_h)) {
      pheatmap_with_consistent_style(
        sub_mat[h_order, h_order, drop = FALSE],
        file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_", col_name, "_hierarchical_structure_heatmap_with_GROUP.pdf")),
        paste(opt[["prefix"]], col_name, "hierarchical structure heatmap + GROUP"),
        annotation_row = ann_h,
        annotation_col = ann_h,
        annotation_colors = make_group_colors(ann_h$GROUP)
      )
    }
    data.frame(sample_id = h_order, group = col_name, taxa = col_name, stringsAsFactors = FALSE) |>
      write.table(file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_", col_name, "_hierarchical_order.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  }

  cluster_df
}

anchor_clusters <- draw_group_outputs(anchor_col)
secondary_clusters <- draw_group_outputs(secondary_col)
all_clusters <- do.call(rbind, Filter(Negate(is.null), list(anchor_clusters, secondary_clusters)))
write.table(all_clusters, file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_cluster_table.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

dup_df <- rbind(
  duplicate_pairs(group_ids(anchor_col), anchor_col),
  duplicate_pairs(group_ids(secondary_col), secondary_col)
)
write.table(dup_df, file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_duplicate_pair_table.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

anchor_ids <- standardize_id(map_df[[anchor_col]])
secondary_ids <- standardize_id(map_df[[secondary_col]])

cross_info <- build_cross_display_matrix(map_df, anchor_col, secondary_col, ibs_mat)
cross_mat <- cross_info$matrix
ann_cross_row <- build_slot_annotation(cross_info$secondary_slots, info_df, secondary_col)
ann_cross_col <- build_slot_annotation(cross_info$anchor_slots, info_df, anchor_col)
ann_cross_colors <- make_group_colors(c(ann_cross_row$GROUP, ann_cross_col$GROUP))
if (!("Missing" %in% names(ann_cross_colors$GROUP))) {
  ann_cross_colors$GROUP <- c(ann_cross_colors$GROUP, Missing = "#000000")
}

pair_table <- data.frame(
  row_index = seq_len(nrow(map_df)),
  anchor_id = anchor_ids,
  secondary_id = secondary_ids,
  anchor_display_id = cross_info$anchor_slots$display_id,
  secondary_display_id = cross_info$secondary_slots$display_id,
  expected_pair_ibs = mapply(function(a, b) {
    if (is.na(a) || is.na(b) || !(a %in% colnames(ibs_mat)) || !(b %in% rownames(ibs_mat))) return(NA_real_)
    ibs_mat[b, a]
  }, anchor_ids, secondary_ids),
  stringsAsFactors = FALSE
)

pair_table$best_anchor <- apply(cross_mat, 1, function(x) {
  if (all(is.na(x))) return(NA_character_)
  names(x)[order(x, decreasing = TRUE, na.last = TRUE)[1]]
})
pair_table$best_anchor_raw <- cross_info$anchor_slots$raw_id[match(pair_table$best_anchor, cross_info$anchor_slots$display_id)]
pair_table$best_anchor_ibs <- mapply(function(sec, best) {
  if (is.na(sec) || is.na(best) || !(sec %in% rownames(ibs_mat)) || !(best %in% colnames(ibs_mat))) return(NA_real_)
  ibs_mat[sec, best]
}, pair_table$secondary_id, pair_table$best_anchor_raw)
pair_table$diagnosis_label <- ifelse(
  is.na(pair_table$expected_pair_ibs),
  "No_data",
  ifelse(pair_table$expected_pair_ibs >= cluster_threshold & pair_table$best_anchor_raw == pair_table$anchor_id,
         "Exact_match",
         ifelse(pair_table$expected_pair_ibs < cluster_threshold & !is.na(pair_table$best_anchor_ibs) & pair_table$best_anchor_ibs >= cluster_threshold,
                "True_mismatch",
                "low_ibs_match"))
)
write.table(pair_table, file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_expected_pair_ibs.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cbind(sample = rownames(cross_mat), cross_mat), file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_", secondary_col, "_x_", anchor_col, "_cross_matrix.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

pair_summary <- pair_table[, c("row_index", "anchor_id", "secondary_id", "expected_pair_ibs", "best_anchor_raw", "best_anchor_ibs", "diagnosis_label")]
write.table(pair_summary, file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_pair_summary.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

diag_mat <- matrix(NA_real_, nrow = nrow(cross_mat), ncol = ncol(cross_mat), dimnames = dimnames(cross_mat))
for (i in seq_len(min(nrow(diag_mat), ncol(diag_mat)))) diag_mat[i, i] <- cross_mat[i, i]
pheatmap_with_consistent_style(
  diag_mat,
  file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_", secondary_col, "_vs_", anchor_col, "_expected_diagonal_heatmap.pdf")),
  paste(opt[["prefix"]], secondary_col, "vs", anchor_col, "expected pair diagonal"),
  annotation_row = ann_cross_row,
  annotation_col = ann_cross_col,
  annotation_colors = ann_cross_colors,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  width_scale = 0.30,
  height_scale = 0.30,
  force_numbers = FALSE
)
pheatmap_with_consistent_style(
  diag_mat,
  file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_", secondary_col, "_vs_", anchor_col, "_expected_diagonal_heatmap_with_IBS_values.pdf")),
  paste(opt[["prefix"]], secondary_col, "vs", anchor_col, "expected pair diagonal + IBS values"),
  annotation_row = ann_cross_row,
  annotation_col = ann_cross_col,
  annotation_colors = ann_cross_colors,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  width_scale = 0.30,
  height_scale = 0.30,
  force_numbers = TRUE
)
pheatmap_with_consistent_style(
  cross_mat,
  file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_", secondary_col, "_x_", anchor_col, "_full_cross_heatmap.pdf")),
  paste(opt[["prefix"]], secondary_col, "x", anchor_col, "full cross-group IBS"),
  annotation_row = ann_cross_row,
  annotation_col = ann_cross_col,
  annotation_colors = ann_cross_colors,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  width_scale = 0.28,
  height_scale = 0.28,
  force_numbers = FALSE
)
pheatmap_with_consistent_style(
  cross_mat,
  file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_", secondary_col, "_x_", anchor_col, "_full_cross_heatmap_with_IBS_values.pdf")),
  paste(opt[["prefix"]], secondary_col, "x", anchor_col, "full cross-group IBS + IBS values"),
  annotation_row = ann_cross_row,
  annotation_col = ann_cross_col,
  annotation_colors = ann_cross_colors,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  width_scale = 0.28,
  height_scale = 0.28,
  force_numbers = TRUE
)

cat("Completed 2group DNA internal and cross-group QC for", opt[["prefix"]], "\n")
