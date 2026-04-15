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
required <- c("mibs", "id", "map", "group-y", "group-x", "outdir", "prefix", "zmin", "zmax", "het-threshold")
missing <- required[!required %in% names(opt)]
if (length(missing) > 0) {
  stop("Missing arguments: ", paste(missing, collapse = ", "))
}

dir.create(opt[["outdir"]], recursive = TRUE, showWarnings = FALSE)

cat("[1/5] Reading IBS matrix and sample map...\n")
ids_raw <- read.table(opt[["id"]], colClasses = "character", stringsAsFactors = FALSE)
ids <- ids_raw$V2
n <- length(ids)

map_df <- read.table(
  opt[["map"]],
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fill = TRUE,
  na.strings = c("", "NA", " ", "-")
)

if (!(opt[["group-y"]] %in% colnames(map_df))) stop("Missing map column: ", opt[["group-y"]])
if (!(opt[["group-x"]] %in% colnames(map_df))) stop("Missing map column: ", opt[["group-x"]])

mibs_vec <- scan(opt[["mibs"]], quiet = TRUE)
mat <- matrix(0, n, n)
idx <- 1
for (i in seq_len(n)) {
  for (j in seq_len(i)) {
    mat[i, j] <- mibs_vec[idx]
    mat[j, i] <- mibs_vec[idx]
    idx <- idx + 1
  }
}
rownames(mat) <- ids
colnames(mat) <- ids

sci_colors <- colorRampPalette(rev(c(
  "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf",
  "#e0f3f8", "#abd9e9", "#74add1", "#4575b4"
)))(100)

zmin <- as.numeric(opt[["zmin"]])
zmax <- as.numeric(opt[["zmax"]])

draw_heatmap <- function(sub_mat, file, title, xlab = "", ylab = "") {
  nx <- ncol(sub_mat)
  ny <- nrow(sub_mat)
  pdf(file, width = max(8, nx * 0.5 + 3), height = max(8, ny * 0.5 + 3))
  par(mar = c(12, 12, 4, 2))
  image_mat <- t(sub_mat)[, ny:1, drop = FALSE]
  image(
    1:nx, 1:ny, image_mat,
    col = sci_colors, axes = FALSE,
    xlab = xlab, ylab = ylab, main = title,
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
      text_col <- ifelse(val > 0.94 | val < 0.78, "white", "black")
      text(i, ny - j + 1, sprintf("%.2f", val), cex = 0.6, col = text_col)
    }
  }
  dev.off()
}

cat("[2/5] Building cross-group matrix...\n")
keep_idx <- !is.na(map_df[[opt[["group-y"]]]]) &
            !is.na(map_df[[opt[["group-x"]]]]) &
            (map_df[[opt[["group-y"]]]] %in% ids) &
            (map_df[[opt[["group-x"]]]] %in% ids)

valid_map <- map_df[keep_idx, , drop = FALSE]
rows_y <- valid_map[[opt[["group-y"]]]]
cols_x <- valid_map[[opt[["group-x"]]]]
sub_mat <- mat[rows_y, cols_x, drop = FALSE]

cross_pdf <- file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_", opt[["group-y"]], "_vs_", opt[["group-x"]], "_IBS.pdf"))
draw_heatmap(
  sub_mat,
  cross_pdf,
  paste(opt[["prefix"]], "IBS:", opt[["group-y"]], "vs", opt[["group-x"]]),
  paste(opt[["group-x"]], "Samples (X)"),
  paste(opt[["group-y"]], "Samples (Y)")
)

cat("[3/5] Summarising best-match pairs...\n")
pair_summary <- data.frame(
  sample_y = rownames(sub_mat),
  expected_x = colnames(sub_mat),
  ibs_expected = diag(sub_mat),
  best_x = NA_character_,
  ibs_best = NA_real_,
  second_best = NA_real_,
  margin = NA_real_,
  status = NA_character_,
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(sub_mat))) {
  vals <- as.numeric(sub_mat[i, ])
  ord <- order(vals, decreasing = TRUE)
  best_idx <- ord[1]
  second_idx <- if (length(ord) >= 2) ord[2] else ord[1]
  pair_summary$best_x[i] <- colnames(sub_mat)[best_idx]
  pair_summary$ibs_best[i] <- vals[best_idx]
  pair_summary$second_best[i] <- vals[second_idx]
  pair_summary$margin[i] <- vals[best_idx] - vals[second_idx]
  pair_summary$status[i] <- ifelse(pair_summary$expected_x[i] == pair_summary$best_x[i], "MATCH", "MISMATCH")
}

write.table(
  pair_summary,
  file = file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_pair_summary.tsv")),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

cat("[4/5] Drawing within-group heatmaps...\n")
for (g in unique(c(opt[["group-y"]], opt[["group-x"]]))) {
  g_samples <- map_df[[g]]
  g_samples <- g_samples[!is.na(g_samples) & g_samples %in% ids]
  g_samples <- unique(g_samples)
  if (length(g_samples) < 2) next
  self_mat <- mat[g_samples, g_samples, drop = FALSE]
  pdf_file <- file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_Self_", g, ".pdf"))
  draw_heatmap(self_mat, pdf_file, paste(opt[["prefix"]], "Internal IBS:", g))
}

cat("[5/5] Writing full IBS matrix...\n")
write.table(
  cbind(sample = rownames(mat), mat),
  file = file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_IBS_matrix.tsv")),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

cat("Report completed.\n")
