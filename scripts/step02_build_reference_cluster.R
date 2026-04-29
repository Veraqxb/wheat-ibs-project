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
opt[["secondary-col"]] <- opt[["secondary-col"]] %||% ""
opt[["cluster-threshold"]] <- opt[["cluster-threshold"]] %||% "0.99"
opt[["zmin"]] <- opt[["zmin"]] %||% "0.7"
opt[["zmax"]] <- opt[["zmax"]] %||% "1.0"

dir.create(opt[["outdir"]], recursive = TRUE, showWarnings = FALSE)
map_df <- read_sample_map(opt[["map"]])
anchor_col <- if (nzchar(opt[["anchor-col"]])) opt[["anchor-col"]] else names(map_df)[1]
secondary_col <- opt[["secondary-col"]]
if (!(anchor_col %in% names(map_df))) stop("Anchor column not found in map: ", anchor_col)
if (nzchar(secondary_col) && !(secondary_col %in% names(map_df))) stop("Secondary column not found in map: ", secondary_col)

mat_df <- read.table(opt[["matrix"]], header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
ibs_mat <- as.matrix(mat_df[, -1, drop = FALSE]); mode(ibs_mat) <- "numeric"
rownames(ibs_mat) <- mat_df[[1]]
colnames(ibs_mat) <- colnames(mat_df)[-1]

anchor_ids <- unique(standardize_id(map_df[[anchor_col]]))
anchor_ids <- anchor_ids[!is.na(anchor_ids) & anchor_ids %in% rownames(ibs_mat)]
if (length(anchor_ids) == 0) stop("No anchor-group samples matched IBS matrix for column: ", anchor_col)

cluster_df <- build_similarity_clusters(anchor_ids, ibs_mat, threshold = as.numeric(opt[["cluster-threshold"]]))
cluster_df$reference_group <- anchor_col
all_cluster_df <- cluster_df
write.table(cluster_df, file = file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_", anchor_col, "_cluster_table.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

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

intersection_df <- data.frame()
if (nzchar(secondary_col)) {
  secondary_ids <- unique(standardize_id(map_df[[secondary_col]]))
  secondary_ids <- secondary_ids[!is.na(secondary_ids) & secondary_ids %in% rownames(ibs_mat)]
  if (length(secondary_ids) > 0) {
    secondary_cluster_df <- build_similarity_clusters(secondary_ids, ibs_mat, threshold = as.numeric(opt[["cluster-threshold"]]))
    secondary_cluster_df$reference_group <- secondary_col
    all_cluster_df <- rbind(all_cluster_df, secondary_cluster_df)
    write.table(secondary_cluster_df, file = file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_", secondary_col, "_cluster_table.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

    ordered_secondary <- unique(c(
      secondary_cluster_df$sample_id[order(secondary_cluster_df$cluster_id, secondary_cluster_df$sample_id)],
      setdiff(secondary_ids, secondary_cluster_df$sample_id)
    ))
    ordered_secondary <- ordered_secondary[ordered_secondary %in% secondary_ids]
    secondary_mat <- ibs_mat[ordered_secondary, ordered_secondary, drop = FALSE]
    draw_heatmap(secondary_mat, file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_", secondary_col, "_cluster_heatmap.pdf")), paste(opt[["prefix"]], secondary_col, "cluster heatmap"), zlim = c(as.numeric(opt[["zmin"]]), as.numeric(opt[["zmax"]])))

    secondary_dup_ids <- unique(secondary_cluster_df$sample_id[secondary_cluster_df$cluster_size > 1])
    if (length(secondary_dup_ids) >= 2) {
      secondary_dup_mat <- ibs_mat[secondary_dup_ids, secondary_dup_ids, drop = FALSE]
      draw_heatmap(secondary_dup_mat, file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_", secondary_col, "_duplicate_heatmap.pdf")), paste(opt[["prefix"]], secondary_col, "duplicate-only heatmap"), zlim = c(as.numeric(opt[["zmin"]]), as.numeric(opt[["zmax"]])))
    }

    intersection_rows <- lapply(seq_len(nrow(map_df)), function(i) {
      anchor_id_i <- standardize_id(map_df[[anchor_col]][i])
      secondary_id_i <- standardize_id(map_df[[secondary_col]][i])
      if (is.na(anchor_id_i) || is.na(secondary_id_i)) return(NULL)
      if (!(anchor_id_i %in% cluster_df$sample_id) || !(secondary_id_i %in% secondary_cluster_df$sample_id)) return(NULL)
      anchor_hit <- cluster_df[cluster_df$sample_id == anchor_id_i, , drop = FALSE]
      secondary_hit <- secondary_cluster_df[secondary_cluster_df$sample_id == secondary_id_i, , drop = FALSE]
      anchor_members <- unique(unlist(strsplit(anchor_hit$cluster_members[1], ";", fixed = TRUE)))
      secondary_members <- unique(unlist(strsplit(secondary_hit$cluster_members[1], ";", fixed = TRUE)))
      shared_members <- intersect(anchor_members, secondary_members)
      data.frame(
        row_index = i,
        anchor_id = anchor_id_i,
        secondary_id = secondary_id_i,
        anchor_cluster_id = anchor_hit$cluster_id[1],
        secondary_cluster_id = secondary_hit$cluster_id[1],
        shared_member_count = length(shared_members),
        shared_members = if (length(shared_members) > 0) paste(shared_members, collapse = ";") else "",
        stringsAsFactors = FALSE
      )
    })
    intersection_df <- do.call(rbind, Filter(Negate(is.null), intersection_rows))
    if (!is.null(intersection_df) && nrow(intersection_df) > 0) {
      write.table(intersection_df, file = file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_cluster_intersection.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
    }
  }
}

write.table(all_cluster_df, file = file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_cluster_table.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

cat("Built reference clusters for", anchor_col, "with", length(anchor_ids), "samples\n")
