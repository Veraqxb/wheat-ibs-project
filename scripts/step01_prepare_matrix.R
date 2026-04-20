#!/usr/bin/env Rscript

source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "ibs_common.R"))

opt <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("mibs", "id", "map", "outdir", "prefix")
missing <- required[!required %in% names(opt)]
if (length(missing) > 0) stop("Missing arguments: ", paste(missing, collapse = ", "))

dir.create(opt[["outdir"]], recursive = TRUE, showWarnings = FALSE)
map_df <- read_sample_map(opt[["map"]])
ibs_mat <- read_mibs_matrix(opt[["mibs"]], opt[["id"]])

map_out <- map_df
for (col in names(map_out)) {
  map_out[[paste0(col, "_in_matrix")]] <- ifelse(!is.na(map_out[[col]]) & map_out[[col]] %in% rownames(ibs_mat), "YES", "NO")
}

write_matrix_tsv(ibs_mat, file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_ibs_matrix.tsv")))
write.table(map_out, file = file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_map_standardized.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

meta_df <- data.frame(
  prefix = opt[["prefix"]],
  sample_count = nrow(ibs_mat),
  map_rows = nrow(map_df),
  map_columns = ncol(map_df),
  anchor_group = names(map_df)[1],
  stringsAsFactors = FALSE
)
write.table(meta_df, file = file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_matrix_metadata.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

cat("Prepared IBS matrix for", opt[["prefix"]], "\n")
