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
opt[["plot-mode"]] <- opt[["plot-mode"]] %||% "dna"
opt[["alert-threshold"]] <- opt[["alert-threshold"]] %||% "0.85"
opt[["match-threshold"]] <- opt[["match-threshold"]] %||% if (opt[["plot-mode"]] == "rna") "0.90" else "0.99"

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
match_threshold <- as.numeric(opt[["match-threshold"]])
alert_threshold <- as.numeric(opt[["alert-threshold"]])

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

write_static_html_report <- function(pair_summary, summary_df, out_file) {
  row_color <- function(status) {
    if (status == "MATCH") return("#d9ead3")
    if (status == "SWAPPED") return("#fce5cd")
    if (status == "MISMATCH") return("#f4cccc")
    if (status == "NO_DATA") return("#d9d9d9")
    "#ffffff"
  }

  con <- file(out_file, "w")
  on.exit(close(con), add = TRUE)
  writeLines("<html><head><meta charset='utf-8'><style>body{font-family:Arial,sans-serif;margin:18px} table{border-collapse:collapse;font-size:12px;margin-bottom:18px} th,td{border:1px solid #999;padding:4px 6px} th{background:#f0f0f0} h2{margin-top:22px}</style></head><body>", con)
  writeLines("<h1>IBS Report Summary</h1>", con)
  writeLines("<h2>Overall Summary</h2>", con)
  writeLines("<table>", con)
  writeLines("<tr>" %+% paste(sprintf("<th>%s</th>", names(summary_df)), collapse = "") %+% "</tr>", con)
  writeLines("<tr>" %+% paste(sprintf("<td>%s</td>", summary_df[1, ]), collapse = "") %+% "</tr>", con)
  writeLines("</table>", con)
  writeLines("<h2>Per-sample Pair Summary</h2>", con)
  writeLines("<table>", con)
  writeLines("<tr>" %+% paste(sprintf("<th>%s</th>", names(pair_summary)), collapse = "") %+% "</tr>", con)
  for (i in seq_len(nrow(pair_summary))) {
    color <- row_color(pair_summary$status[i])
    vals <- ifelse(is.na(pair_summary[i, ]), "", as.character(pair_summary[i, ]))
    writeLines("<tr style='background:" %+% color %+% "'>" %+% paste(sprintf("<td>%s</td>", vals), collapse = "") %+% "</tr>", con)
  }
  writeLines("</table></body></html>", con)
}

`%+%` <- function(a, b) paste0(a, b)

write_interactive_html_report <- function(pair_summary, summary_df, out_file) {
  if (!requireNamespace("DT", quietly = TRUE) ||
      !requireNamespace("htmltools", quietly = TRUE) ||
      !requireNamespace("htmlwidgets", quietly = TRUE)) {
    write_static_html_report(pair_summary, summary_df, out_file)
    return(invisible(FALSE))
  }

  decision_colors <- c(MATCH = "#d9ead3", SWAPPED = "#fce5cd", MISMATCH = "#f4cccc", NO_DATA = "#d9d9d9")
  row_callback <- paste0(
    "function(row, data) {",
    "var status = data[", which(names(pair_summary) == "status") - 1, "];",
    "if (status === 'MATCH') { $(row).css({'background-color':'", decision_colors["MATCH"], "'}); }",
    "else if (status === 'SWAPPED') { $(row).css({'background-color':'", decision_colors["SWAPPED"], "'}); }",
    "else if (status === 'MISMATCH') { $(row).css({'background-color':'", decision_colors["MISMATCH"], "'}); }",
    "else if (status === 'NO_DATA') { $(row).css({'background-color':'", decision_colors["NO_DATA"], "'}); }",
    "}"
  )

  summary_widget <- DT::datatable(
    summary_df,
    rownames = FALSE,
    filter = "none",
    options = list(dom = "t", paging = FALSE, ordering = FALSE, autoWidth = TRUE)
  )

  pair_widget <- DT::datatable(
    pair_summary,
    rownames = FALSE,
    filter = "top",
    extensions = c("Buttons"),
    options = list(
      dom = "Bfrtip",
      buttons = c("copy", "csv", "excel"),
      pageLength = 25,
      autoWidth = TRUE,
      scrollX = TRUE,
      rowCallback = DT::JS(row_callback)
    )
  )

  page <- htmltools::tagList(
    htmltools::tags$head(
      htmltools::tags$meta(charset = "utf-8"),
      htmltools::tags$style(htmltools::HTML("
        body { font-family: Arial, sans-serif; margin: 18px; }
        h1, h2 { margin-bottom: 10px; }
        .section { margin-bottom: 24px; }
      "))
    ),
    htmltools::tags$h1("Interactive IBS Report"),
    htmltools::tags$div(class = "section",
      htmltools::tags$h2("Overall Summary"),
      summary_widget
    ),
    htmltools::tags$div(class = "section",
      htmltools::tags$h2("Per-sample Pair Summary"),
      htmltools::tags$p("Use the search box, column filters, sorting, and export buttons to inspect sample results."),
      pair_widget
    )
  )

  htmlwidgets::saveWidget(page, file = out_file, selfcontained = TRUE)
  invisible(TRUE)
}

get_upper_triangle_values <- function(m) {
  if (is.null(m) || nrow(m) < 2 || ncol(m) < 2) return(numeric(0))
  m[lower.tri(m, diag = FALSE)]
}

build_distribution_df <- function(mat, valid_group_y, valid_group_x, pair_summary) {
  y_ids <- unique(valid_group_y$match)
  x_ids <- unique(valid_group_x$match)

  y_ids <- y_ids[!is.na(y_ids) & y_ids %in% rownames(mat)]
  x_ids <- x_ids[!is.na(x_ids) & x_ids %in% rownames(mat)]

  y_internal <- numeric(0)
  x_internal <- numeric(0)

  if (length(y_ids) >= 2) {
    y_mat <- mat[y_ids, y_ids, drop = FALSE]
    y_internal <- get_upper_triangle_values(y_mat)
    y_internal <- y_internal[!is.na(y_internal)]
  }

  if (length(x_ids) >= 2) {
    x_mat <- mat[x_ids, x_ids, drop = FALSE]
    x_internal <- get_upper_triangle_values(x_mat)
    x_internal <- x_internal[!is.na(x_internal)]
  }

  paired_vals <- pair_summary$ibs_expected
  paired_vals <- paired_vals[!is.na(paired_vals)]

  cross_nonpair <- numeric(0)
  if (length(valid_group_y$match) > 0 && length(valid_group_x$match) > 0) {
    cross_mat <- mat[valid_group_y$match, valid_group_x$match, drop = FALSE]
    diag(cross_mat) <- NA
    cross_nonpair <- as.numeric(cross_mat)
    cross_nonpair <- cross_nonpair[!is.na(cross_nonpair)]
  }

  dist_df <- data.frame(
    IBS = c(y_internal, x_internal, paired_vals, cross_nonpair),
    Category = c(
      rep(paste0(opt[["group-y"]], "_internal"), length(y_internal)),
      rep(paste0(opt[["group-x"]], "_internal"), length(x_internal)),
      rep("paired_1to1", length(paired_vals)),
      rep("cross_nonpaired", length(cross_nonpair))
    ),
    stringsAsFactors = FALSE
  )

  dist_df$Category <- factor(
    dist_df$Category,
    levels = c(
      paste0(opt[["group-y"]], "_internal"),
      paste0(opt[["group-x"]], "_internal"),
      "paired_1to1",
      "cross_nonpaired"
    )
  )
  dist_df
}

build_raincloud_df <- function(pair_summary, plot_mode, group_x_label, group_y_label) {
  if (plot_mode == "dna") {
    rain_df <- data.frame(
      Sample = rep(pair_summary$sample_y, 3),
      Type = factor(
        rep(c(paste0(group_x_label, "_internal"), "paired_match", paste0(group_y_label, "_internal")), each = nrow(pair_summary)),
        levels = c(paste0(group_x_label, "_internal"), "paired_match", paste0(group_y_label, "_internal"))
      ),
      IBS = c(pair_summary$ibs_x_background, pair_summary$ibs_expected, pair_summary$ibs_y_background),
      Match_Type = rep(pair_summary$match_type, 3),
      stringsAsFactors = FALSE
    )
  } else {
    rain_df <- data.frame(
      Sample = rep(pair_summary$sample_y, 2),
      Type = factor(
        rep(c("internal_background", "paired_match"), each = nrow(pair_summary)),
        levels = c("internal_background", "paired_match")
      ),
      IBS = c(pair_summary$ibs_y_background, pair_summary$ibs_expected),
      Match_Type = rep(pair_summary$match_type, 2),
      stringsAsFactors = FALSE
    )
  }
  rain_df[!is.na(rain_df$IBS), , drop = FALSE]
}

draw_match_type_heatmap <- function(pair_summary, out_file, title) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || nrow(pair_summary) == 0) {
    return(invisible(NULL))
  }

  type_levels <- c(
    "1_Unique_Match",
    "2_Clonal_Match",
    "3_Swapped_Mismatch",
    "4_True_Mismatch",
    "5_No_Data"
  )
  type_labels <- c(
    "Unique",
    "Clonal",
    "Swapped",
    "TrueMismatch",
    "NoData"
  )
  type_colors <- c(
    "1_Unique_Match" = "#0072B2",
    "2_Clonal_Match" = "#009E73",
    "3_Swapped_Mismatch" = "#D55E00",
    "4_True_Mismatch" = "#BDBDBD",
    "5_No_Data" = "#000000"
  )

  plot_df <- do.call(
    rbind,
    lapply(seq_len(nrow(pair_summary)), function(i) {
      data.frame(
        sample_y = pair_summary$sample_y[i],
        category = factor(type_levels, levels = type_levels),
        fill_group = ifelse(type_levels == pair_summary$match_type[i], type_levels, "inactive"),
        label = ifelse(
          type_levels == pair_summary$match_type[i],
          paste0(pair_summary$best_x[i] %||% "", ifelse(is.na(pair_summary$ibs_best[i]), "", sprintf("\n%.3f", pair_summary$ibs_best[i]))),
          ""
        ),
        stringsAsFactors = FALSE
      )
    })
  )

  plot_df$sample_y <- factor(plot_df$sample_y, levels = rev(pair_summary$sample_y))
  fill_values <- c(type_colors, inactive = "#FFFFFF")

  pdf(out_file, width = 9, height = max(6, nrow(pair_summary) * 0.26 + 2.5))
  print(
    ggplot2::ggplot(plot_df, ggplot2::aes(x = category, y = sample_y, fill = fill_group)) +
      ggplot2::geom_tile(color = "grey80", linewidth = 0.3) +
      ggplot2::geom_text(
        data = plot_df[plot_df$fill_group != "inactive", , drop = FALSE],
        ggplot2::aes(label = label),
        size = 2.4,
        lineheight = 0.9,
        color = "white",
        fontface = "bold"
      ) +
      ggplot2::scale_fill_manual(values = fill_values, guide = "none") +
      ggplot2::scale_x_discrete(labels = type_labels) +
      ggplot2::labs(
        title = title,
        subtitle = "Each row marks the active DNA pairing class for one sample",
        x = "",
        y = ""
      ) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
        plot.subtitle = ggplot2::element_text(hjust = 0.5),
        axis.text.x = ggplot2::element_text(face = "bold"),
        axis.text.y = ggplot2::element_text(size = 8),
        panel.grid = ggplot2::element_blank()
      )
  )
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
  second_x = NA_character_,
  second_best = NA_real_,
  third_x = NA_character_,
  third_best = NA_real_,
  margin = NA_real_,
  high_match_count = NA_integer_,
  high_match_ids = NA_character_,
  duplicate_x_ids = NA_character_,
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
    pair_summary$second_x[i] <- NA_character_
    pair_summary$second_best[i] <- NA_real_
    pair_summary$third_x[i] <- NA_character_
    pair_summary$third_best[i] <- NA_real_
    pair_summary$margin[i] <- NA_real_
    pair_summary$high_match_count[i] <- 0L
    pair_summary$high_match_ids[i] <- ""
    pair_summary$duplicate_x_ids[i] <- ""
    pair_summary$status[i] <- "NO_DATA"
    next
  }
  ord <- order(vals, decreasing = TRUE, na.last = TRUE)
  best_idx <- ord[1]
  second_idx <- if (length(ord) >= 2) ord[2] else ord[1]
  third_idx <- if (length(ord) >= 3) ord[3] else second_idx
  high_idx <- ord[!is.na(vals[ord]) & vals[ord] >= match_threshold]
  high_labels <- names(vals)[high_idx]
  high_pairs <- if (length(high_idx) > 0) paste(sprintf("%s(%.4f)", names(vals)[high_idx], vals[high_idx]), collapse = ";") else ""

  pair_summary$ibs_expected[i] <- unname(mat[sample_y_match, expected_x_match])
  pair_summary$best_x[i] <- names(vals)[best_idx]
  pair_summary$ibs_best[i] <- vals[best_idx]
  pair_summary$second_x[i] <- names(vals)[second_idx]
  pair_summary$second_best[i] <- vals[second_idx]
  pair_summary$third_x[i] <- names(vals)[third_idx]
  pair_summary$third_best[i] <- vals[third_idx]
  pair_summary$margin[i] <- ifelse(is.na(vals[best_idx]) || is.na(vals[second_idx]), NA_real_, vals[best_idx] - vals[second_idx])
  pair_summary$high_match_count[i] <- length(high_idx)
  pair_summary$high_match_ids[i] <- high_pairs
  pair_summary$duplicate_x_ids[i] <- if (length(high_labels) > 1) paste(high_labels, collapse = ";") else ""
  pair_summary$status[i] <- ifelse(is.na(pair_summary$best_x[i]), "NO_DATA", ifelse(expected_x == pair_summary$best_x[i], "MATCH", "MISMATCH"))
}

pair_summary$ibs_y_background <- NA_real_
pair_summary$ibs_x_background <- NA_real_
for (i in seq_len(nrow(pair_summary))) {
  sy <- pair_summary$sample_y_match[i]
  sx <- pair_summary$expected_x_match[i]
  other_y <- pair_summary$sample_y_match[-i]
  other_y <- other_y[!is.na(other_y) & other_y %in% rownames(mat)]
  if (length(other_y) > 0 && sy %in% rownames(mat)) {
    pair_summary$ibs_y_background[i] <- mean(mat[sy, other_y], na.rm = TRUE)
  }
  other_x <- pair_summary$expected_x_match[-i]
  other_x <- other_x[!is.na(other_x) & other_x %in% rownames(mat)]
  if (length(other_x) > 0 && sx %in% rownames(mat)) {
    pair_summary$ibs_x_background[i] <- mean(mat[sx, other_x], na.rm = TRUE)
  }
}

pair_summary$match_type <- NA_character_

for (i in seq_len(nrow(pair_summary))) {
  if (pair_summary$status[i] == "NO_DATA") {
    pair_summary$match_type[i] <- "5_No_Data"
  } else if (!is.na(pair_summary$ibs_expected[i]) &&
             pair_summary$ibs_expected[i] >= match_threshold &&
             pair_summary$best_x[i] == pair_summary$expected_x[i]) {
    if (!is.na(pair_summary$high_match_count[i]) && pair_summary$high_match_count[i] > 1) {
      pair_summary$match_type[i] <- "2_Clonal_Match"
    } else {
      pair_summary$match_type[i] <- "1_Unique_Match"
    }
  } else if (!is.na(pair_summary$ibs_best[i]) && pair_summary$ibs_best[i] >= match_threshold) {
    pair_summary$match_type[i] <- "3_Swapped_Mismatch"
  } else {
    pair_summary$match_type[i] <- "4_True_Mismatch"
  }
}

pair_summary$status <- ifelse(
  pair_summary$match_type %in% c("1_Unique_Match", "2_Clonal_Match"), "MATCH",
  ifelse(pair_summary$match_type == "3_Swapped_Mismatch", "SWAPPED",
    ifelse(pair_summary$match_type == "5_No_Data", "NO_DATA", "MISMATCH")
  )
)

summary_df <- data.frame(
  prefix = opt[["prefix"]],
  group_y = opt[["group-y"]],
  group_x = opt[["group-x"]],
  total_pairs_in_map = nrow(map_df),
  valid_pairs_used = nrow(valid_map),
  matched_pair_count = sum(pair_summary$status == "MATCH", na.rm = TRUE),
  swapped_pair_count = sum(pair_summary$status == "SWAPPED", na.rm = TRUE),
  mismatched_pair_count = sum(pair_summary$status == "MISMATCH", na.rm = TRUE),
  no_data_pair_count = sum(pair_summary$status == "NO_DATA", na.rm = TRUE),
  unique_match_count = sum(pair_summary$match_type == "1_Unique_Match", na.rm = TRUE),
  clonal_match_count = sum(pair_summary$match_type == "2_Clonal_Match", na.rm = TRUE),
  swapped_mismatch_count = sum(pair_summary$match_type == "3_Swapped_Mismatch", na.rm = TRUE),
  true_mismatch_count = sum(pair_summary$match_type == "4_True_Mismatch", na.rm = TRUE),
  stringsAsFactors = FALSE
)

write.table(
  pair_summary,
  file = file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_pair_summary.tsv")),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)
write.table(
  summary_df,
  file = file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_summary.tsv")),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)
write.table(
  pair_summary[pair_summary$match_type == "3_Swapped_Mismatch", , drop = FALSE],
  file = file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_swapped_mismatch_samples.tsv")),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)
write.table(
  pair_summary[pair_summary$match_type == "4_True_Mismatch", , drop = FALSE],
  file = file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_true_mismatch_samples.tsv")),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)
write.table(
  pair_summary[pair_summary$match_type == "5_No_Data", , drop = FALSE],
  file = file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_no_data_samples.tsv")),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)
draw_match_type_heatmap(
  pair_summary,
  file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_DNA_match_type_heatmap.pdf")),
  paste(opt[["prefix"]], "DNA Pairing Classes")
)
interactive_ok <- write_interactive_html_report(
  pair_summary,
  summary_df,
  file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_report_interactive.html"))
)
write_static_html_report(pair_summary, summary_df, file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_report.html")))

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

  dup_mat <- self_mat
  diag(dup_mat) <- NA
  dup_idx <- which(dup_mat >= match_threshold, arr.ind = TRUE)
  if (nrow(dup_idx) > 0) {
    dup_ids <- unique(rownames(dup_mat)[dup_idx[, 1]])
    dup_sub <- self_mat[dup_ids, dup_ids, drop = FALSE]
    draw_heatmap(
      dup_sub,
      file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_Internal_Duplicates_", g, ".pdf")),
      paste(opt[["prefix"]], "Potential Internal Duplicates:", g),
      show_values = TRUE
    )
  }
}

cat("[5/6] Writing full IBS matrix...\n")

write.table(
  cbind(sample = rownames(mat), mat),
  file = file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_IBS_matrix.tsv")),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

mismatch_df <- pair_summary[pair_summary$status == "MISMATCH", , drop = FALSE]
if (nrow(mismatch_df) > 0) {
  mismatch_rows <- mismatch_df$sample_y_match
  mismatch_cols <- unique(c(mismatch_df$expected_x_match, cols_x_match[match(mismatch_df$best_x, cols_x_label)]))
  mismatch_rows <- mismatch_rows[!is.na(mismatch_rows) & mismatch_rows %in% rownames(mat)]
  mismatch_cols <- mismatch_cols[!is.na(mismatch_cols) & mismatch_cols %in% colnames(mat)]
  if (length(mismatch_rows) > 0 && length(mismatch_cols) > 0) {
    mismatch_mat <- mat[mismatch_rows, mismatch_cols, drop = FALSE]
    rownames(mismatch_mat) <- mismatch_df$sample_y[match(mismatch_rows, mismatch_df$sample_y_match)]
    mapped_cols <- cols_x_label[match(colnames(mismatch_mat), cols_x_match)]
    colnames(mismatch_mat) <- ifelse(is.na(mapped_cols), colnames(mismatch_mat), mapped_cols)
    draw_heatmap(
      mismatch_mat,
      file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_Mismatch_Audit.pdf")),
      paste(opt[["prefix"]], "Mismatch Audit"),
      show_values = TRUE
    )
  }
}

swapped_df <- pair_summary[pair_summary$status == "SWAPPED", , drop = FALSE]
if (nrow(swapped_df) > 0) {
  swapped_rows <- swapped_df$sample_y_match
  swapped_cols <- unique(c(swapped_df$expected_x_match, cols_x_match[match(swapped_df$best_x, cols_x_label)]))
  swapped_rows <- swapped_rows[!is.na(swapped_rows) & swapped_rows %in% rownames(mat)]
  swapped_cols <- swapped_cols[!is.na(swapped_cols) & swapped_cols %in% colnames(mat)]
  if (length(swapped_rows) > 0 && length(swapped_cols) > 0) {
    swapped_mat <- mat[swapped_rows, swapped_cols, drop = FALSE]
    rownames(swapped_mat) <- swapped_df$sample_y[match(swapped_rows, swapped_df$sample_y_match)]
    mapped_cols <- cols_x_label[match(colnames(swapped_mat), cols_x_match)]
    colnames(swapped_mat) <- ifelse(is.na(mapped_cols), colnames(swapped_mat), mapped_cols)
    draw_heatmap(
      swapped_mat,
      file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_Swapped_Audit.pdf")),
      paste(opt[["prefix"]], "Swapped Audit"),
      show_values = TRUE
    )
  }
}

if (requireNamespace("ggplot2", quietly = TRUE)) {
  cat("[6/8] Drawing raincloud and density distribution plots...\n")
  library(ggplot2)
  has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)

  rain_df <- build_raincloud_df(pair_summary, opt[["plot-mode"]], opt[["group-x"]], opt[["group-y"]])

  if (nrow(rain_df) > 0) {
    c_unique <- sum(pair_summary$match_type == "1_Unique_Match", na.rm = TRUE)
    c_clone_in <- sum(pair_summary$match_type == "2_Clonal_Match", na.rm = TRUE)
    c_clone_out <- sum(pair_summary$match_type == "3_Swapped_Mismatch", na.rm = TRUE)
    c_mismatch <- sum(pair_summary$match_type == "4_True_Mismatch", na.rm = TRUE)
    c_nodata <- sum(pair_summary$match_type == "5_No_Data", na.rm = TRUE)

    stat_text <- sprintf(
      "Sample Pairing Audit\n(IBS Threshold = %.2f)\n------------------------------\n[1] Unique Match (Blue): %d\n[2] Clonal Match (Green): %d\n[3] Swapped Mismatch (Orange): %d\n[4] True Mismatch (Grey): %d\n[5] No Data (Black): %d",
      match_threshold, c_unique, c_clone_in, c_clone_out, c_mismatch, c_nodata
    )

    df_true_mismatch <- rain_df[rain_df$Match_Type == "4_True_Mismatch", , drop = FALSE]
    df_nodata <- rain_df[rain_df$Match_Type == "5_No_Data", , drop = FALSE]
    df_colored <- rain_df[rain_df$Match_Type %in% c("1_Unique_Match", "2_Clonal_Match", "3_Swapped_Mismatch"), , drop = FALSE]
    df_alert <- rain_df[rain_df$IBS < alert_threshold & rain_df$Type == "paired_match", , drop = FALSE]

    min_y <- min(0.8, min(rain_df$IBS, na.rm = TRUE))
    max_y <- max(rain_df$IBS, na.rm = TRUE)
    y_range <- max_y - min_y

    p_rain <- ggplot() +
      geom_violin(data = rain_df, aes(x = Type, y = IBS, fill = Type), trim = FALSE, alpha = 0.2, color = NA, width = 0.6) +
      geom_boxplot(data = rain_df, aes(x = Type, y = IBS, fill = Type), width = 0.15, outlier.shape = NA, alpha = 0.4, color = "black") +
      geom_line(data = df_true_mismatch, aes(x = Type, y = IBS, group = Sample), color = "grey80", alpha = 0.5, linewidth = 0.4) +
      geom_line(data = df_nodata, aes(x = Type, y = IBS, group = Sample), color = "black", alpha = 0.45, linewidth = 0.45, linetype = "dashed") +
      geom_line(data = df_colored, aes(x = Type, y = IBS, group = Sample, color = Match_Type), linewidth = 0.9, alpha = 0.8) +
      geom_point(data = df_true_mismatch, aes(x = Type, y = IBS), fill = "grey75", color = "grey85", size = 2, alpha = 0.5, shape = 21) +
      geom_point(data = df_nodata, aes(x = Type, y = IBS), fill = "black", color = "black", size = 2.1, alpha = 0.55, shape = 21) +
      geom_point(data = df_colored, aes(x = Type, y = IBS, fill = Type), size = 2.5, alpha = 0.9, shape = 21, color = "white") +
      geom_hline(yintercept = alert_threshold, color = "#d73027", linetype = "dashed", linewidth = 0.8, alpha = 0.7) +
      annotate("text", x = median(seq_along(levels(rain_df$Type))), y = alert_threshold, label = sprintf("IBS = %.2f Alert Baseline", alert_threshold), color = "#d73027", vjust = -0.6, fontface = "bold", size = 4)

    if (nrow(df_alert) > 0) {
      if (has_ggrepel) {
        p_rain <- p_rain + ggrepel::geom_text_repel(
          data = df_alert, aes(x = Type, y = IBS, label = Sample),
          color = "#d73027", fontface = "bold", size = 3.5,
          nudge_x = -0.15, direction = "y", segment.color = "grey50"
        )
      } else {
        p_rain <- p_rain + geom_text(
          data = df_alert, aes(x = Type, y = IBS, label = Sample),
          color = "#d73027", fontface = "bold", size = 3.5, vjust = 1.5, hjust = 1.1
        )
      }
    }

    p_rain <- p_rain +
      annotate("label", x = if (opt[["plot-mode"]] == "dna") 3.55 else 2.8, y = max_y + (y_range * 0.04), label = stat_text,
               hjust = 0, vjust = 1, fill = "#f8f9fa", color = "black",
               label.size = 0.8, fontface = "bold", size = 4.3, alpha = 0.9) +
      theme_bw(base_size = 15) +
      labs(
        title = paste("Raincloud IBS Audit:", opt[["prefix"]]),
        subtitle = if (opt[["plot-mode"]] == "dna") "DNA mode: Z23 internal vs paired match vs B25 internal" else "RNA mode: internal background vs paired match",
        x = "", y = "Identity By State (IBS) Score"
      ) +
      scale_fill_manual(values = if (opt[["plot-mode"]] == "dna") {
        vals <- c("#4575b4", "#fdae61", "#74add1")
        names(vals) <- levels(rain_df$Type)
        vals
      } else {
        vals <- c("#4575b4", "#fdae61")
        names(vals) <- levels(rain_df$Type)
        vals
      }) +
      scale_color_manual(values = c(
        "1_Unique_Match" = "#0072B2",
        "2_Clonal_Match" = "#009E73",
        "3_Swapped_Mismatch" = "#D55E00"
      )) +
      coord_cartesian(ylim = c(min(0.7, min_y - (y_range * 0.05)), max(1.0, max_y + (y_range * 0.05))), clip = "off") +
      theme(
        panel.border = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.8),
        legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "grey30", size = 11),
        axis.text.x = element_text(face = "bold", color = "black", size = 13),
        panel.grid.major.x = element_blank(),
        plot.margin = margin(t = 30, r = 260, b = 15, l = 15)
      )

    rain_out <- file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_density_distribution.pdf"))
    pdf(rain_out, width = 11.5, height = 7.5)
    print(p_rain)
    dev.off()
  }

  dist_df <- build_distribution_df(mat, valid_group_y, valid_group_x, pair_summary)
  if (nrow(dist_df) > 0) {
    if (opt[["plot-mode"]] == "rna") {
      dist_df <- dist_df[dist_df$Category %in% c(paste0(opt[["group-y"]], "_internal"), "paired_1to1", "cross_nonpaired"), , drop = FALSE]
      dist_df$Category <- factor(dist_df$Category, levels = c(paste0(opt[["group-y"]], "_internal"), "paired_1to1", "cross_nonpaired"))
      density_cols <- c("#4575b4", "#d73027", "#fdae61")
      names(density_cols) <- c(paste0(opt[["group-y"]], "_internal"), "paired_1to1", "cross_nonpaired")
    } else {
      density_cols <- c("#4575b4", "#74add1", "#d73027", "#fdae61")
      names(density_cols) <- c(
        paste0(opt[["group-y"]], "_internal"),
        paste0(opt[["group-x"]], "_internal"),
        "paired_1to1",
        "cross_nonpaired"
      )
    }

    p_density <- ggplot(dist_df, aes(x = IBS, color = Category, fill = Category)) +
      geom_density(alpha = 0.18, linewidth = 1.2, adjust = 1.1) +
      geom_vline(xintercept = match_threshold, linetype = "dashed", color = "grey35", linewidth = 0.7) +
      geom_vline(xintercept = alert_threshold, linetype = "dotted", color = "#d73027", linewidth = 0.7) +
      scale_color_manual(values = density_cols) +
      scale_fill_manual(values = density_cols) +
      labs(
        title = paste("IBS Density Comparison:", opt[["prefix"]]),
        subtitle = if (opt[["plot-mode"]] == "dna") "DNA mode: within-group, paired 1-to-1, and cross-group non-paired IBS distributions" else "RNA mode: internal, paired 1-to-1, and cross-group non-paired IBS distributions",
        x = "IBS",
        y = "Density",
        color = "Category",
        fill = "Category"
      ) +
      theme_bw(base_size = 14) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "right",
        panel.grid.minor = element_blank()
      ) +
      coord_cartesian(xlim = c(min(0.7, min(dist_df$IBS, na.rm = TRUE)), 1.0))

    density_out <- file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_ibs_density_comparison.pdf"))
    pdf(density_out, width = 9, height = 6.5)
    print(p_density)
    dev.off()
  }
}

cat("[6/6] Report completed successfully.\n")
cat("Available map columns:", paste(available_cols, collapse = ", "), "\n")
cat("Valid pairs used:", nrow(valid_map), "\n")
cat("MATCH:", sum(pair_summary$status == "MATCH"), "\n")
cat("SWAPPED:", sum(pair_summary$status == "SWAPPED"), "\n")
cat("MISMATCH:", sum(pair_summary$status == "MISMATCH"), "\n")
cat("NO_DATA:", sum(pair_summary$status == "NO_DATA"), "\n")
