#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
if (is.na(script_arg)) stop("Cannot determine script path for run_2cluster.R")
script_dirname <- dirname(normalizePath(sub("^--file=", "", script_arg)))
runner <- file.path(script_dirname, "two_group_internal_heatmaps.R")
density_runner <- file.path(script_dirname, "two_group_ibs_density_thresholds.R")
source(file.path(script_dirname, "ibs_common.R"))

parse_bool <- function(x) {
  if (is.logical(x)) return(x)
  x <- tolower(trimws(as.character(x)))
  x %in% c("true", "t", "1", "yes", "y")
}

msg <- function(fmt, ...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%F %T"), sprintf(fmt, ...)))

get_arg <- function(opt, nm, default = NA_character_) {
  if (!nm %in% names(opt) || is.null(opt[[nm]]) || identical(opt[[nm]], "")) default else opt[[nm]]
}

get_map_file <- function(mibs_path, map_dir) {
  f <- basename(mibs_path)
  if (grepl("^C2", f, ignore.case = TRUE)) {
    file.path(map_dir, "c2_id_map.txt")
  } else if (grepl("^C4", f, ignore.case = TRUE)) {
    file.path(map_dir, "c4_id_map.txt")
  } else if (grepl("^C6", f, ignore.case = TRUE)) {
    file.path(map_dir, "c6_id_map.txt")
  } else {
    NA_character_
  }
}

build_cmd <- function(row, defaults = list()) {
  args <- c(
    runner,
    "--mibs", row[["mibs"]],
    "--id", row[["id"]],
    "--map", row[["map"]],
    "--outdir", row[["outdir"]],
    "--prefix", row[["prefix"]],
    "--anchor-col", row[["anchor_col"]],
    "--secondary-col", row[["secondary_col"]],
    "--cluster-threshold", row[["cluster_threshold"]],
    "--zmin", row[["zmin"]],
    "--zmax", row[["zmax"]]
  )

  if (!is.na(row[["info"]]) && nzchar(row[["info"]])) {
    args <- c(args, "--info", row[["info"]])
  } else if (!is.na(defaults$info_file) && nzchar(defaults$info_file)) {
    args <- c(args, "--info", defaults$info_file)
  }

  args
}

build_density_cmd <- function(row) {
  c(
    density_runner,
    "--mibs", row[["mibs"]],
    "--id", row[["id"]],
    "--map", row[["map"]],
    "--outdir", file.path(row[["outdir"]], "density_thresholds"),
    "--prefix", row[["prefix"]],
    "--anchor-col", row[["anchor_col"]],
    "--secondary-col", row[["secondary_col"]],
    "--xmin", row[["zmin"]],
    "--xmax", row[["zmax"]]
  )
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))

config_file <- get_arg(opt, "config", NA_character_)
ibs_dir <- get_arg(opt, "ibs-dir", "ibs_file")
map_dir <- get_arg(opt, "map-dir", "maps")
info_file <- get_arg(opt, "info-file", "camp_info.txt")
main_outdir <- get_arg(opt, "outdir", "all_2group_internal_heatmaps")
default_anchor <- get_arg(opt, "anchor-col", "Z23")
default_secondary <- get_arg(opt, "secondary-col", "B25")
default_threshold <- get_arg(opt, "cluster-threshold", "0.99")
default_zmin <- get_arg(opt, "zmin", "0.7")
default_zmax <- get_arg(opt, "zmax", "1.0")
default_rscript <- get_arg(opt, "rscript", "Rscript")
run_density <- parse_bool(get_arg(opt, "density", "TRUE"))

dir.create(main_outdir, recursive = TRUE, showWarnings = FALSE)

run_dataset <- function(row, defaults) {
  row <- as.list(row)
  row$anchor_col <- row$anchor_col %||% defaults$anchor_col
  row$secondary_col <- row$secondary_col %||% defaults$secondary_col
  row$cluster_threshold <- row$cluster_threshold %||% defaults$cluster_threshold
  row$zmin <- row$zmin %||% defaults$zmin
  row$zmax <- row$zmax %||% defaults$zmax
  row$info <- row$info %||% defaults$info_file

  if (is.na(row$mibs) || !file.exists(row$mibs)) stop("Missing mibs file: ", row$mibs)
  if (is.na(row$id) || !file.exists(row$id)) stop("Missing id file: ", row$id)
  if (is.na(row$map) || !file.exists(row$map)) stop("Missing map file: ", row$map)
  if (!is.na(row$info) && nzchar(row$info) && !file.exists(row$info)) stop("Missing info file: ", row$info)

  dir.create(row$outdir, recursive = TRUE, showWarnings = FALSE)

  cmd_args <- build_cmd(row, defaults)
  msg("Running: %s", row$prefix)
  status <- system2(default_rscript, cmd_args)
  if (!identical(status, 0L)) stop("2group QC failed for prefix: ", row$prefix)

  if (run_density) {
    msg("Running density threshold analysis: %s", row$prefix)
    density_status <- system2(default_rscript, build_density_cmd(row))
    if (!identical(density_status, 0L)) stop("2group density threshold analysis failed for prefix: ", row$prefix)
  }
}

defaults <- list(
  anchor_col = default_anchor,
  secondary_col = default_secondary,
  cluster_threshold = default_threshold,
  zmin = default_zmin,
  zmax = default_zmax,
  info_file = info_file
)

if (!is.na(config_file) && file.exists(config_file)) {
  cfg <- fread(config_file, header = TRUE, sep = "\t", fill = TRUE, data.table = FALSE, check.names = FALSE)
  required_cols <- c("mibs", "id", "map", "outdir", "prefix")
  missing_cols <- setdiff(required_cols, names(cfg))
  if (length(missing_cols) > 0) stop("Batch config missing columns: ", paste(missing_cols, collapse = ", "))

  if (!("anchor_col" %in% names(cfg))) cfg$anchor_col <- defaults$anchor_col
  if (!("secondary_col" %in% names(cfg))) cfg$secondary_col <- defaults$secondary_col
  if (!("cluster_threshold" %in% names(cfg))) cfg$cluster_threshold <- defaults$cluster_threshold
  if (!("zmin" %in% names(cfg))) cfg$zmin <- defaults$zmin
  if (!("zmax" %in% names(cfg))) cfg$zmax <- defaults$zmax
  if (!("info" %in% names(cfg))) cfg$info <- defaults$info_file

  for (i in seq_len(nrow(cfg))) {
    run_dataset(cfg[i, , drop = FALSE], defaults)
  }
  msg("Completed batch 2group QC for %d datasets", nrow(cfg))
  quit(status = 0)
}

mibs_files <- list.files(ibs_dir, pattern = "\\.mibs$", full.names = TRUE)
mibs_files <- mibs_files[grepl("2group|2groups", basename(mibs_files), ignore.case = TRUE)]
if (length(mibs_files) == 0) stop("No 2group .mibs files found in: ", ibs_dir)

msg("Found %d candidate 2group mibs files.", length(mibs_files))

for (mibs_file in mibs_files) {
  id_file <- paste0(mibs_file, ".id")
  map_file <- get_map_file(mibs_file, map_dir)
  prefix <- sub("\\.mibs$", "", basename(mibs_file))

  row <- data.frame(
    mibs = mibs_file,
    id = id_file,
    map = map_file,
    outdir = file.path(main_outdir, prefix),
    prefix = prefix,
    anchor_col = defaults$anchor_col,
    secondary_col = defaults$secondary_col,
    cluster_threshold = defaults$cluster_threshold,
    zmin = defaults$zmin,
    zmax = defaults$zmax,
    info = defaults$info_file,
    stringsAsFactors = FALSE
  )

  if (!file.exists(id_file)) {
    msg("Skip %s: missing id file", prefix)
    next
  }
  if (is.na(map_file) || !file.exists(map_file)) {
    msg("Skip %s: missing map file", prefix)
    next
  }
  if (!is.na(info_file) && nzchar(info_file) && !file.exists(info_file)) {
    msg("Skip %s: missing info file", prefix)
    next
  }

  run_dataset(row, defaults)
}

msg("All done.")
