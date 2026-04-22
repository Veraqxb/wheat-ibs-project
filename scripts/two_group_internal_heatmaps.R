#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
if (is.na(script_arg)) stop("Cannot determine script path for two_group_internal_heatmaps.R")
script_dirname <- dirname(normalizePath(sub("^--file=", "", script_arg)))
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

dir.create(opt[["outdir"]], recursive = TRUE, showWarnings = FALSE)

map_df <- read_sample_map(opt[["map"]])
ibs_mat <- read_mibs_matrix(opt[["mibs"]], opt[["id"]])

if (!(anchor_col %in% names(map_df))) stop("Missing anchor column in map: ", anchor_col)
if (!(secondary_col %in% names(map_df))) stop("Missing secondary column in map: ", secondary_col)

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

draw_group_outputs <- function(col_name) {
  ids <- unique(group_ids(col_name))
  if (length(ids) == 0) return(NULL)

  cluster_df <- build_similarity_clusters(ids, ibs_mat, threshold = cluster_threshold)
  cluster_df$reference_group <- col_name
  write.table(cluster_df, file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_", col_name, "_cluster_table.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

  assignment_df <- cluster_df[, c("reference_group", "sample_id", "cluster_id", "cluster_size", "cluster_members", "neighbor_ids", "max_neighbor_ibs"), drop = FALSE]
  write.table(assignment_df, file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_", col_name, "_assignment_table.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

  sub_mat <- ibs_mat[ids, ids, drop = FALSE]
  cluster_order <- cluster_df$sample_id[order(ifelse(is.na(cluster_df$cluster_id), "singleton", cluster_df$cluster_id), cluster_df$sample_id)]
  cluster_order <- unique(cluster_order[cluster_order %in% rownames(sub_mat)])
  if (length(cluster_order) > 1) {
    draw_heatmap(sub_mat[cluster_order, cluster_order, drop = FALSE], file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_", col_name, "_threshold_cluster_heatmap.pdf")), paste(opt[["prefix"]], col_name, "IBS threshold-cluster heatmap"), zlim = c(zmin, zmax))
  }

  h_order <- ordered_for_structure(sub_mat)
  if (length(h_order) > 1) {
    draw_heatmap(sub_mat[h_order, h_order, drop = FALSE], file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_", col_name, "_hierarchical_structure_heatmap.pdf")), paste(opt[["prefix"]], col_name, "hierarchical structure heatmap"), zlim = c(zmin, zmax))
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
keep_pair <- !is.na(anchor_ids) & !is.na(secondary_ids) & anchor_ids %in% colnames(ibs_mat) & secondary_ids %in% rownames(ibs_mat)
pair_table <- data.frame(
  row_index = which(keep_pair),
  anchor_id = anchor_ids[keep_pair],
  secondary_id = secondary_ids[keep_pair],
  expected_pair_ibs = mapply(function(a, b) ibs_mat[b, a], anchor_ids[keep_pair], secondary_ids[keep_pair]),
  stringsAsFactors = FALSE
)

cross_mat <- ibs_mat[secondary_ids[keep_pair], anchor_ids[keep_pair], drop = FALSE]
rownames(cross_mat) <- secondary_ids[keep_pair]
colnames(cross_mat) <- anchor_ids[keep_pair]

pair_table$best_anchor <- apply(cross_mat, 1, function(x) {
  if (all(is.na(x))) return(NA_character_)
  names(x)[order(x, decreasing = TRUE, na.last = TRUE)[1]]
})
pair_table$best_anchor_ibs <- mapply(function(sec, best) if (is.na(best)) NA_real_ else ibs_mat[sec, best], pair_table$secondary_id, pair_table$best_anchor)
pair_table$diagnosis_label <- ifelse(
  is.na(pair_table$expected_pair_ibs),
  "No_data",
  ifelse(pair_table$expected_pair_ibs >= cluster_threshold & pair_table$best_anchor == pair_table$anchor_id,
         "Exact_match",
         ifelse(pair_table$expected_pair_ibs < cluster_threshold & !is.na(pair_table$best_anchor_ibs) & pair_table$best_anchor_ibs >= cluster_threshold,
                "True_mismatch",
                "low_ibs_match"))
)
write.table(pair_table, file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_expected_pair_ibs.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cbind(sample = rownames(cross_mat), cross_mat), file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_", secondary_col, "_x_", anchor_col, "_cross_matrix.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

diag_mat <- matrix(NA_real_, nrow = nrow(cross_mat), ncol = ncol(cross_mat), dimnames = dimnames(cross_mat))
diag(diag_mat) <- diag(cross_mat)
draw_heatmap(diag_mat, file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_", secondary_col, "_vs_", anchor_col, "_expected_diagonal_heatmap.pdf")), paste(opt[["prefix"]], secondary_col, "vs", anchor_col, "expected pair diagonal"), zlim = c(zmin, zmax))
draw_heatmap(cross_mat, file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_", secondary_col, "_x_", anchor_col, "_full_cross_heatmap.pdf")), paste(opt[["prefix"]], secondary_col, "x", anchor_col, "full cross-group IBS"), zlim = c(zmin, zmax))

cat("Completed 2group DNA internal and cross-group QC for", opt[["prefix"]], "\n")
