#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
if (is.na(script_arg)) stop("Cannot determine script path for run_2cluster.R")
script_dirname <- dirname(normalizePath(sub("^--file=", "", script_arg)))
source(file.path(script_dirname, "ibs_common.R"))

opt <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("config")
missing <- required[!required %in% names(opt)]
if (length(missing) > 0) stop("Missing arguments: ", paste(missing, collapse = ", "))

cfg <- read.table(opt[["config"]], header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
required_cols <- c("mibs", "id", "map", "outdir", "prefix")
missing_cols <- setdiff(required_cols, names(cfg))
if (length(missing_cols) > 0) stop("Batch config missing columns: ", paste(missing_cols, collapse = ", "))

runner <- file.path(script_dirname, "two_group_internal_heatmaps.R")
for (i in seq_len(nrow(cfg))) {
  row <- cfg[i, , drop = FALSE]
  args <- c(
    runner,
    "--mibs", row$mibs,
    "--id", row$id,
    "--map", row$map,
    "--outdir", row$outdir,
    "--prefix", row$prefix,
    "--anchor-col", if ("anchor_col" %in% names(row)) row$anchor_col else "Z23",
    "--secondary-col", if ("secondary_col" %in% names(row)) row$secondary_col else "B25",
    "--cluster-threshold", if ("cluster_threshold" %in% names(row)) as.character(row$cluster_threshold) else "0.99"
  )
  message("Running 2group QC: ", row$prefix)
  status <- system2("Rscript", args)
  if (!identical(status, 0L)) stop("2group QC failed for prefix: ", row$prefix)
}

cat("Completed batch 2group cluster QC for", nrow(cfg), "datasets\n")
