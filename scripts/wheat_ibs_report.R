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
required <- c("mibs", "id", "map", "group-y", "group-x", "outdir", "prefix", "zmin", "zmax")
missing <- required[!required %in% names(opt)]
if (length(missing) > 0) {
  stop("Missing arguments: ", paste(missing, collapse = ", "))
}

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

group_y_vals <- map_df[[opt[["group-y"]]]]
group_x_vals <- map_df[[opt[["group-x"]]]]

keep_idx <- !is.na(group_y_vals) &
            !is.na(group_x_vals) &
            (group_y_vals %in% ids) &
            (group_x_vals %in% ids)

valid_map <- map_df[keep_idx, , drop = FALSE]
if (nrow(valid_map) == 0) {
  stop(
    "No valid matched pairs found between map and IBS IDs.\n",
    "Group Y: ", opt[["group-y"]], "\n",
    "Group X: ", opt[["group-x"]], "\n",
    "Tip: check whether sample names in map file exactly match those in ", opt[["id"]]
  )
}

rows_y <- valid_map[[opt[["group-y"]]]]
cols_x <- valid_map[[opt[["group-x"]]]]
sub_mat <- mat[rows_y, cols_x, drop = FALSE]

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
  sample_y = rows_y,
  expected_x = cols_x,
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
  expected_x <- pair_summary$expected_x[i]

  vals <- as.numeric(mat[sample_y, cols_x])
  names(vals) <- cols_x
  ord <- order(vals, decreasing = TRUE, na.last = TRUE)
  best_idx <- ord[1]
  second_idx <- if (length(ord) >= 2) ord[2] else ord[1]

  pair_summary$ibs_expected[i] <- unname(mat[sample_y, expected_x])
  pair_summary$best_x[i] <- names(vals)[best_idx]
  pair_summary$ibs_best[i] <- vals[best_idx]
  pair_summary$second_best[i] <- vals[second_idx]
  pair_summary$margin[i] <- vals[best_idx] - vals[second_idx]
  pair_summary$status[i] <- ifelse(expected_x == pair_summary$best_x[i], "MATCH", "MISMATCH")
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
  g_samples <- map_df[[g]]
  g_samples <- g_samples[!is.na(g_samples) & g_samples %in% ids]
  g_samples <- unique(g_samples)

  if (length(g_samples) < 2) {
    cat("Skip self-heatmap for group", g, ": fewer than 2 valid samples.\n")
    next
  }

  self_mat <- mat[g_samples, g_samples, drop = FALSE]
  pdf_file <- file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_Self_", g, ".pdf"))
  draw_heatmap(
    self_mat,
    pdf_file,
    paste(opt[["prefix"]], "Internal IBS:", g),
    show_values = length(g_samples) <= 80
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
