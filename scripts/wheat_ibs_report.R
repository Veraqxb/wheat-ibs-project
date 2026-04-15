#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(x) {
  res <- list()
  i <- 1
  while (i <= length(x)) {
    key <- x[[i]]
    if (!startsWith(key, "--")) stop("Invalid argument: ", key)
    if (i == length(x)) stop("Missing value for ", key)
    res[[sub("^--", "", key)]] <- x[[i + 1]]
    i <- i + 2
  }
  res
}

opt <- parse_args(args)
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

required <- c("mibs", "id", "map", "group-y", "group-x", "outdir", "prefix", "zmin", "zmax")
missing <- required[!required %in% names(opt)]
if (length(missing) > 0) {
  stop("Missing arguments: ", paste(missing, collapse = ", "))
}

opt[["group-y-match-col"]] <- opt[["group-y-match-col"]] %||% ""
opt[["group-x-match-col"]] <- opt[["group-x-match-col"]] %||% ""
opt[["group-y-match-mode"]] <- opt[["group-y-match-mode"]] %||% "direct"
opt[["group-x-match-mode"]] <- opt[["group-x-match-mode"]] %||% "direct"

dir.create(opt[["outdir"]], recursive = TRUE, showWarnings = FALSE)

cat("[1/6] Reading IBS matrix and sample map...\n")

ids_raw <- read.table(opt[["id"]], colClasses = "character", stringsAsFactors = FALSE)
if (ncol(ids_raw) < 2) {
  stop("IBS id file format is invalid: expected at least 2 columns in ", opt[["id"]])
}

ids <- ids_raw$V2
n <- length(ids)
if (n == 0) {
  stop("No sample IDs found in ", opt[["id"]])
}

map_df <- read.table(
  opt[["map"]],
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fill = TRUE,
  na.strings = c("", "NA", " ", "-")
)

if (ncol(map_df) == 0) {
  stop("Map file has no columns: ", opt[["map"]])
}

available_cols <- colnames(map_df)
if (!(opt[["group-y"]] %in% available_cols)) {
  stop(
    "Missing map column: ", opt[["group-y"]],
    "\nAvailable columns: ", paste(available_cols, collapse = ", ")
  )
}
if (!(opt[["group-x"]] %in% available_cols)) {
  stop(
    "Missing map column: ", opt[["group-x"]],
    "\nAvailable columns: ", paste(available_cols, collapse = ", ")
  )
}

if (nzchar(opt[["group-y-match-col"]]) && !(opt[["group-y-match-col"]] %in% available_cols)) {
  stop(
    "Missing map column for group-y-match-col: ", opt[["group-y-match-col"]],
    "\nAvailable columns: ", paste(available_cols, collapse = ", ")
  )
}
if (nzchar(opt[["group-x-match-col"]]) && !(opt[["group-x-match-col"]] %in% available_cols)) {
  stop(
    "Missing map column for group-x-match-col: ", opt[["group-x-match-col"]],
    "\nAvailable columns: ", paste(available_cols, collapse = ", ")
  )
}

mibs_vec <- scan(opt[["mibs"]], quiet = TRUE)
expected_tri <- n * (n + 1) / 2
expected_square <- n * n

if (length(mibs_vec) == expected_square) {
  cat("Detected square IBS matrix format.\n")
  mat <- matrix(mibs_vec, nrow = n, ncol = n, byrow = TRUE)
} else if (length(mibs_vec) == expected_tri) {
  cat("Detected lower-triangle IBS matrix format.\n")
  mat <- matrix(0, n, n)
  idx <- 1
  for (i in seq_len(n)) {
    for (j in seq_len(i)) {
      mat[i, j] <- mibs_vec[idx]
      mat[j, i] <- mibs_vec[idx]
      idx <- idx + 1
    }
  }
} else {
  stop(
    "IBS matrix size mismatch: expected either ",
    expected_square, " (square) or ",
    expected_tri, " (lower triangle) values for ",
    n, " samples, but got ", length(mibs_vec),
    ". File: ", opt[["mibs"]]
  )
}

rownames(mat) <- ids
colnames(mat) <- ids

rdylbu_base <- c(
  "#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090",
  "#FFFFBF",
  "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4", "#313695"
)
rdylbu_cols <- colorRampPalette(rev(rdylbu_base))(100)

zmin <- as.numeric(opt[["zmin"]])
zmax <- as.numeric(opt[["zmax"]])

apply_match_mode <- function(x, mode) {
  if (mode == "direct") return(x)
  if (mode == "b25_to_2") return(sub("^B25C2_", "2_", x))
  stop("Unsupported match mode: ", mode)
}

resolve_group_ids <- function(map_df, logical_col, match_col, match_mode) {
  logical_ids <- map_df[[logical_col]]
  raw_match_ids <- if (nzchar(match_col)) map_df[[match_col]] else logical_ids
  match_ids <- apply_match_mode(raw_match_ids, match_mode)
  list(label = logical_ids, match = match_ids)
}

draw_heatmap <- function(sub_mat, file, title, xlab = "", ylab = "", show_values = TRUE) {
  nx <- ncol(sub_mat)
  ny <- nrow(sub_mat)

  if (nx == 0 || ny == 0) {
    warning("Skip heatmap because matrix is empty: ", file)
    return(invisible(NULL))
  }

  pdf(file, width = max(8, nx * 0.5 + 3), height = max(8, ny * 0.5 + 3))
  par(mar = c(12, 12, 4, 2))

  image_mat <- t(sub_mat)[, ny:1, drop = FALSE]
  image(
    1:nx, 1:ny, image_mat,
    col = rdylbu_cols,
    axes = FALSE,
    xlab = xlab,
    ylab = ylab,
    main = title,
    zlim = c(zmin, zmax)
  )

  axis(1, at = 1:nx, labels = colnames(sub_mat), las = 2, cex.axis = 0.7)
  axis(2, at = 1:ny, labels = rev(rownames(sub_mat)), las = 1, cex.axis = 0.7)
  abline(h = seq(0.5, ny + 0.5, by = 1), col = "grey75")
  abline(v = seq(0.5, nx + 0.5, by = 1), col = "grey75")
  box(col = "grey50")

  for (i in seq_len(nx)) {
    for (j in seq_len(ny)) {
      val <- sub_mat[j, i]
      plot_y <- ny - j + 1

      if (is.na(val)) {
        rect(i - 0.5, plot_y - 0.5, i + 0.5, plot_y + 0.5, col = "grey85", border = "grey60", lwd = 0.8)
        if (show_values) {
          text(i, plot_y, "NA", cex = 0.55, col = "black")
        }
        next
      }

      if (isTRUE(all.equal(val, 1, tolerance = 1e-8))) {
        rect(i - 0.5, plot_y - 0.5, i + 0.5, plot_y + 0.5, border = "black", lwd = 2.8)
        points(i, plot_y, pch = 8, cex = 1.1, col = "black")
      } else if (val >= 0.99) {
        rect(i - 0.5, plot_y - 0.5, i + 0.5, plot_y + 0.5, border = "black", lwd = 1.6)
      }

      if (show_values) {
        text_col <- ifelse(val > 0.94 | val < 0.78, "white", "black")
        text_cex <- 0.6
        font_face <- 1
        if (isTRUE(all.equal(val, 1, tolerance = 1e-8))) {
          text_col <- "black"
          text_cex <- 0.78
          font_face <- 2
        } else if (val >= 0.99) {
          text_col <- "black"
          text_cex <- 0.68
          font_face <- 2
        }
        text(i, plot_y, sprintf("%.2f", val), cex = text_cex, col = text_col, font = font_face)
      }
    }
  }

  dev.off()
}

cat("[2/6] Building cross-group matrix...\n")

group_y <- resolve_group_ids(map_df, opt[["group-y"]], opt[["group-y-match-col"]], opt[["group-y-match-mode"]])
group_x <- resolve_group_ids(map_df, opt[["group-x"]], opt[["group-x-match-col"]], opt[["group-x-match-mode"]])

keep_idx <- !is.na(group_y$label) &
            !is.na(group_x$label) &
            !is.na(group_y$match) &
            !is.na(group_x$match) &
            (group_y$match %in% ids) &
            (group_x$match %in% ids)

valid_map <- map_df[keep_idx, , drop = FALSE]
if (nrow(valid_map) == 0) {
  stop(
    "No valid matched pairs found between map and IBS IDs.\n",
    "Group Y: ", opt[["group-y"]], "\n",
    "Group X: ", opt[["group-x"]], "\n",
    "Tip: check whether sample names in map file exactly match those in ", opt[["id"]]
  )
}

valid_group_y <- resolve_group_ids(valid_map, opt[["group-y"]], opt[["group-y-match-col"]], opt[["group-y-match-mode"]])
valid_group_x <- resolve_group_ids(valid_map, opt[["group-x"]], opt[["group-x-match-col"]], opt[["group-x-match-mode"]])

rows_y_label <- valid_group_y$label
rows_y_match <- valid_group_y$match
cols_x_label <- valid_group_x$label
cols_x_match <- valid_group_x$match

sub_mat <- mat[rows_y_match, cols_x_match, drop = FALSE]
rownames(sub_mat) <- rows_y_label
colnames(sub_mat) <- cols_x_label

if (nrow(sub_mat) == 0 || ncol(sub_mat) == 0) {
  stop("Cross-group IBS sub-matrix is empty after filtering.")
}

cross_pdf <- file.path(
  opt[["outdir"]],
  paste0(opt[["prefix"]], "_", opt[["group-y"]], "_vs_", opt[["group-x"]], "_IBS.pdf")
)
draw_heatmap(
  sub_mat,
  cross_pdf,
  paste(opt[["prefix"]], "IBS:", opt[["group-y"]], "vs", opt[["group-x"]]),
  paste(opt[["group-x"]], "Samples (X)"),
  paste(opt[["group-y"]], "Samples (Y)"),
  show_values = nrow(sub_mat) <= 80 && ncol(sub_mat) <= 80
)

cat("[3/6] Summarising best-match pairs...\n")

pair_summary <- data.frame(
  sample_y = rows_y_label,
  sample_y_match = rows_y_match,
  expected_x = cols_x_label,
  expected_x_match = cols_x_match,
  ibs_expected = NA_real_,
  best_x = NA_character_,
  ibs_best = NA_real_,
  second_best = NA_real_,
  margin = NA_real_,
  status = NA_character_,
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(pair_summary))) {
  sample_y <- pair_summary$sample_y[i]
  sample_y_match <- pair_summary$sample_y_match[i]
  expected_x <- pair_summary$expected_x[i]
  expected_x_match <- pair_summary$expected_x_match[i]

  vals <- as.numeric(mat[sample_y_match, cols_x_match])
  names(vals) <- cols_x_label
  if (all(is.na(vals))) {
    pair_summary$ibs_expected[i] <- NA_real_
    pair_summary$best_x[i] <- NA_character_
    pair_summary$ibs_best[i] <- NA_real_
    pair_summary$second_best[i] <- NA_real_
    pair_summary$margin[i] <- NA_real_
    pair_summary$status[i] <- "NO_DATA"
    next
  }
  ord <- order(vals, decreasing = TRUE, na.last = TRUE)
  best_idx <- ord[1]
  second_idx <- if (length(ord) >= 2) ord[2] else ord[1]

  pair_summary$ibs_expected[i] <- unname(mat[sample_y_match, expected_x_match])
  pair_summary$best_x[i] <- names(vals)[best_idx]
  pair_summary$ibs_best[i] <- vals[best_idx]
  pair_summary$second_best[i] <- vals[second_idx]
  pair_summary$margin[i] <- ifelse(is.na(vals[best_idx]) || is.na(vals[second_idx]), NA_real_, vals[best_idx] - vals[second_idx])
  pair_summary$status[i] <- ifelse(is.na(pair_summary$best_x[i]), "NO_DATA", ifelse(expected_x == pair_summary$best_x[i], "MATCH", "MISMATCH"))
}

write.table(
  pair_summary,
  file = file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_pair_summary.tsv")),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

summary_df <- data.frame(
  prefix = opt[["prefix"]],
  group_y = opt[["group-y"]],
  group_x = opt[["group-x"]],
  total_pairs_in_map = nrow(map_df),
  valid_pairs_used = nrow(valid_map),
  matched_pair_count = sum(pair_summary$status == "MATCH"),
  mismatched_pair_count = sum(pair_summary$status == "MISMATCH"),
  no_data_pair_count = sum(pair_summary$status == "NO_DATA"),
  stringsAsFactors = FALSE
)

write.table(
  summary_df,
  file = file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_summary.tsv")),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

cat("[4/6] Drawing within-group heatmaps...\n")

for (g in unique(c(opt[["group-y"]], opt[["group-x"]]))) {
  if (g == opt[["group-y"]]) {
    g_resolved <- resolve_group_ids(map_df, g, opt[["group-y-match-col"]], opt[["group-y-match-mode"]])
  } else if (g == opt[["group-x"]]) {
    g_resolved <- resolve_group_ids(map_df, g, opt[["group-x-match-col"]], opt[["group-x-match-mode"]])
  } else {
    g_resolved <- resolve_group_ids(map_df, g, "", "direct")
  }

  keep_g <- !is.na(g_resolved$label) & !is.na(g_resolved$match) & (g_resolved$match %in% ids)
  g_df <- data.frame(label = g_resolved$label[keep_g], match = g_resolved$match[keep_g], stringsAsFactors = FALSE)
  g_df <- g_df[!duplicated(g_df$match), , drop = FALSE]

  if (nrow(g_df) < 2) {
    cat("Skip self-heatmap for group", g, ": fewer than 2 valid samples.\n")
    next
  }

  self_mat <- mat[g_df$match, g_df$match, drop = FALSE]
  rownames(self_mat) <- g_df$label
  colnames(self_mat) <- g_df$label
  pdf_file <- file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_Self_", g, ".pdf"))
  draw_heatmap(
    self_mat,
    pdf_file,
    paste(opt[["prefix"]], "Internal IBS:", g),
    show_values = nrow(g_df) <= 80
  )
}

cat("[5/6] Writing full IBS matrix...\n")

write.table(
  cbind(sample = rownames(mat), mat),
  file = file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_IBS_matrix.tsv")),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

cat("[6/6] Report completed successfully.\n")
cat("Available map columns:", paste(available_cols, collapse = ", "), "\n")
cat("Valid pairs used:", nrow(valid_map), "\n")
cat("MATCH:", sum(pair_summary$status == "MATCH"), "\n")
cat("MISMATCH:", sum(pair_summary$status == "MISMATCH"), "\n")
cat("NO_DATA:", sum(pair_summary$status == "NO_DATA"), "\n")
