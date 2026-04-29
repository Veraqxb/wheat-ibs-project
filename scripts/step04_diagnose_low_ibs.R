#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- if (length(script_arg) == 0 || is.na(script_arg)) getwd() else sub("^--file=", "", script_arg)
script_path <- gsub("~\\+~", " ", script_path, fixed = FALSE)
source(file.path(dirname(normalizePath(script_path)), "ibs_common.R"))

opt <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("matrix", "map", "pairwise-dir", "outdir", "prefix")
missing <- required[!required %in% names(opt)]
if (length(missing) > 0) stop("Missing arguments: ", paste(missing, collapse = ", "))

opt[["anchor-col"]] <- opt[["anchor-col"]] %||% ""
opt[["zmin"]] <- opt[["zmin"]] %||% "0.7"
opt[["zmax"]] <- opt[["zmax"]] %||% "1.0"

dir.create(opt[["outdir"]], recursive = TRUE, showWarnings = FALSE)
map_df <- read_sample_map(opt[["map"]])
anchor_col <- if (nzchar(opt[["anchor-col"]])) opt[["anchor-col"]] else names(map_df)[1]

mat_df <- read.table(opt[["matrix"]], header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
ibs_mat <- as.matrix(mat_df[, -1, drop = FALSE]); mode(ibs_mat) <- "numeric"
rownames(ibs_mat) <- mat_df[[1]]
colnames(ibs_mat) <- colnames(mat_df)[-1]

pair_files <- list.files(opt[["pairwise-dir"]], pattern = "_pairwise\\.tsv$", full.names = TRUE)
if (length(pair_files) == 0) stop("No pairwise tables found in ", opt[["pairwise-dir"]])

all_low <- list()
for (pair_file in pair_files) {
  pair_df <- read.table(pair_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  low_df <- pair_df[pair_df$match_type %in% c("SWAPPED", "LOW_IBS", "NO_DATA"), , drop = FALSE]
  if (nrow(low_df) == 0) next
  group_name <- unique(low_df$group_name)
  group_name <- group_name[!is.na(group_name)][1]
  if (is.na(group_name) || !nzchar(group_name)) {
    group_name <- sub(".*_([^_/]+)_vs_.*_pairwise\\.tsv$", "\\1", basename(pair_file))
  }

  query_ids <- unique(low_df$query_id[!is.na(low_df$query_id) & low_df$query_id %in% rownames(ibs_mat)])
  anchor_ids <- unique(c(low_df$anchor_id, low_df$best_anchor, low_df$second_anchor, low_df$third_anchor))
  anchor_ids <- anchor_ids[!is.na(anchor_ids) & anchor_ids %in% colnames(ibs_mat)]
  if (length(query_ids) > 0 && length(anchor_ids) > 0) {
    sub_mat <- ibs_mat[query_ids, anchor_ids, drop = FALSE]
    draw_heatmap(
      sub_mat,
      file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_", group_name, "_low_ibs_heatmap.pdf")),
      paste(opt[["prefix"]], group_name, "low-IBS cross heatmap"),
      zlim = c(as.numeric(opt[["zmin"]]), as.numeric(opt[["zmax"]]))
    )
  }
  all_low[[group_name]] <- low_df
}

if (length(all_low) > 0) {
  low_all <- do.call(rbind, all_low)
  write.table(low_all, file = file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_low_ibs_all.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
}

cat("Diagnosed low-IBS samples from", length(pair_files), "pairwise tables\n")
