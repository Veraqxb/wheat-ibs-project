#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(data.table)
  library(ggplot2)
})

parse_args <- function(args) {
  res <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[i]
    if (startsWith(key, "--")) {
      val <- if (i < length(args) && !startsWith(args[i + 1], "--")) args[i + 1] else TRUE
      res[[sub("^--", "", key)]] <- val
      if (!identical(val, TRUE)) i <- i + 1
    }
    i <- i + 1
  }
  res
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
`%||%` <- function(x, y) if (is.null(x)) y else x
default_xlsx <- "/Users/veraqiu/Desktop/CAMP_sample_identify/CAMP_DNA_check_withCluster_260426.xlsx"
if (!file.exists(default_xlsx)) {
  default_xlsx <- "/Users/veraqiu/Documents/New project/CAMP_DNA_check_withCluster_260426.xlsx"
}
xlsx_file <- opt[["xlsx"]] %||% default_xlsx
sheet_name <- opt[["sheet"]] %||% "C2"
outdir <- opt[["outdir"]] %||% file.path(dirname(xlsx_file), "bidirectional_scatter_260427")
supp_file <- opt[["supplement"]] %||% NA_character_
ibs_dir <- opt[["ibs-dir"]] %||% file.path(dirname(xlsx_file), "ibs_file")
map_dir <- opt[["map-dir"]] %||% file.path(dirname(xlsx_file), "maps")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

clean_char <- function(x) {
  x <- as.character(x)
  x <- gsub("\\u00A0", " ", x)
  x <- gsub("\\u3000", " ", x)
  x <- gsub("^\\ufeff", "", x)
  trimws(x)
}

clean_numeric <- function(x) {
  x <- clean_char(x)
  x[x %in% c("", "NA", "N/A", "NaN", "NULL", "null", "-")] <- NA
  suppressWarnings(as.numeric(x))
}

read_mibs_matrix <- function(mibs_file, id_file) {
  ids <- fread(id_file, header = FALSE, colClasses = "character")
  if (ncol(ids) < 2) stop("ID file must contain at least two columns: ", id_file)
  sample_ids <- clean_char(ids[[2]])
  n <- length(sample_ids)

  mat_df <- read.table(mibs_file, header = FALSE, fill = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  mat_df <- mat_df[, colSums(is.na(mat_df)) < nrow(mat_df), drop = FALSE]
  mat <- as.matrix(mat_df)
  if (ncol(mat) == n + 1) mat <- mat[, -1, drop = FALSE]
  if (nrow(mat) != n || ncol(mat) != n) {
    stop(sprintf("Dimension mismatch for %s: matrix %dx%d, ids=%d", basename(mibs_file), nrow(mat), ncol(mat), n))
  }
  suppressWarnings(storage.mode(mat) <- "numeric")
  if (any(is.na(mat[upper.tri(mat)]))) {
    mat_t <- t(mat)
    mat[upper.tri(mat)] <- mat_t[upper.tri(mat)]
  }
  rownames(mat) <- sample_ids
  colnames(mat) <- sample_ids
  mat
}

get_ploidy_files <- function(sheet_name, ibs_dir, map_dir) {
  tag <- toupper(sheet_name)
  if (tag == "C2") {
    mibs <- file.path(ibs_dir, "C2_2groups_renamed.mibs")
    map <- file.path(map_dir, "c2_id_map.txt")
  } else if (tag == "C4") {
    mibs <- file.path(ibs_dir, "C4_2groups_renamed.mibs")
    map <- file.path(map_dir, "c4_id_map.txt")
  } else if (tag == "C6") {
    mibs <- file.path(ibs_dir, "C6_2groups_final_qc.mibs")
    map <- file.path(map_dir, "c6_id_map.txt")
  } else {
    stop("Unsupported sheet/ploidy: ", sheet_name)
  }
  list(mibs = mibs, id = paste0(mibs, ".id"), map = map)
}

read_sample_map <- function(map_file) {
  dt <- fread(map_file, header = TRUE, fill = TRUE, colClasses = "character")
  setnames(dt, clean_char(names(dt)))
  for (nm in names(dt)) dt[[nm]] <- clean_char(dt[[nm]])
  dt
}

fill_from_ibs_if_available <- function(dt, sheet_name, ibs_dir, map_dir) {
  files <- get_ploidy_files(sheet_name, ibs_dir, map_dir)
  if (!file.exists(files$mibs) || !file.exists(files$id) || !file.exists(files$map)) {
    msg("Skip IBS backfill: missing mibs/id/map for %s", sheet_name)
    return(dt)
  }

  mat <- read_mibs_matrix(files$mibs, files$id)
  map_dt <- read_sample_map(files$map)
  z23_col <- if ("Z23" %in% names(map_dt)) "Z23" else if ("CAMP编号" %in% names(map_dt)) "CAMP编号" else NA_character_
  if (is.na(z23_col) || !("B25" %in% names(map_dt))) {
    msg("Skip IBS backfill: map file lacks Z23/CAMP编号 or B25 columns.")
    return(dt)
  }

  pair_map <- unique(map_dt[, .(expected_Z23 = get(z23_col), expected_B25 = B25)])
  pair_map <- pair_map[!is.na(expected_Z23) & expected_Z23 != ""]
  dt <- merge(dt, pair_map, by = "expected_Z23", all.x = TRUE, sort = FALSE, suffixes = c("", "_map"))
  if ("expected_B25_map" %in% names(dt)) {
    dt[is.na(expected_B25) | expected_B25 == "", expected_B25 := expected_B25_map]
    dt[, expected_B25_map := NULL]
  }

  dt[, best_Z23_ibs_from_matrix := mapply(function(b25, z23) {
    if (is.na(b25) || is.na(z23) || !(b25 %in% rownames(mat)) || !(z23 %in% colnames(mat))) return(NA_real_)
    as.numeric(mat[b25, z23])
  }, expected_B25, best_Z23)]

  dt[, best_B25_ibs_from_matrix := mapply(function(z23, b25) {
    if (is.na(z23) || is.na(b25) || !(z23 %in% rownames(mat)) || !(b25 %in% colnames(mat))) return(NA_real_)
    as.numeric(mat[z23, b25])
  }, expected_Z23, best_B25)]

  dt[, expected_pair_ibs_from_matrix := mapply(function(z23, b25) {
    if (is.na(z23) || is.na(b25) || !(z23 %in% rownames(mat)) || !(b25 %in% colnames(mat))) return(NA_real_)
    as.numeric(mat[z23, b25])
  }, expected_Z23, expected_B25)]

  dt[is.na(best_Z23_ibs), best_Z23_ibs := best_Z23_ibs_from_matrix]
  dt[is.na(best_B25_ibs), best_B25_ibs := best_B25_ibs_from_matrix]
  dt[is.na(expected_pair_ibs), expected_pair_ibs := expected_pair_ibs_from_matrix]
  dt
}

msg <- function(fmt, ...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%F %T"), sprintf(fmt, ...)))
}

theme_pub <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey30"),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      axis.line = element_line(linewidth = 0.5),
      axis.ticks = element_line(linewidth = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    )
}

extract_suffix <- function(x) {
  x <- clean_char(x)
  sub("^[A-Za-z]+[0-9]*", "", x)
}

extract_suffix_for_expected <- function(id, z23_prefix, b25_prefix) {
  id <- clean_char(id)
  out <- id
  if (!is.na(z23_prefix) && nzchar(z23_prefix)) {
    out <- ifelse(startsWith(out, z23_prefix), sub(paste0("^", z23_prefix), "", out), out)
  }
  if (!is.na(b25_prefix) && nzchar(b25_prefix)) {
    out <- ifelse(startsWith(out, b25_prefix), sub(paste0("^", b25_prefix), "", out), out)
  }
  out
}

infer_prefix <- function(x) {
  x <- clean_char(x)
  pref <- sub("^([A-Za-z]+[0-9]*).*", "\\1", x)
  pref[pref == x | pref == ""] <- NA_character_
  pref
}

infer_expected_ids <- function(dt) {
  dt <- copy(dt)
  z23_pref <- na.omit(unique(infer_prefix(dt$best_Z23)))
  b25_pref <- na.omit(unique(infer_prefix(dt$best_B25)))

  z23_prefix <- if (length(z23_pref) > 0) z23_pref[1] else ifelse(grepl("^Z23", dt$ID[1]), "Z23", "")
  b25_prefix <- if (length(b25_pref) > 0) b25_pref[1] else ifelse(grepl("^B25", dt$ID[1]), "B25", "")

  dt[, id_suffix := extract_suffix_for_expected(ID, z23_prefix, b25_prefix)]
  dt[, expected_Z23 := ifelse(startsWith(ID, z23_prefix), ID, paste0(z23_prefix, id_suffix))]
  dt[, expected_B25 := ifelse(startsWith(ID, b25_prefix), ID, paste0(b25_prefix, id_suffix))]
  dt
}

sheet_df <- as.data.table(read_excel(xlsx_file, sheet = sheet_name, col_names = TRUE))
setnames(sheet_df, clean_char(names(sheet_df)))
for (nm in names(sheet_df)) sheet_df[[nm]] <- clean_char(sheet_df[[nm]])

note_cols <- grep("^note(\\.\\.\\.[0-9]+)?$", names(sheet_df), value = TRUE)
if ("...8" %in% names(sheet_df)) setnames(sheet_df, "...8", "cluster_note")
if (!("cluster_note" %in% names(sheet_df)) && length(note_cols) > 0) {
  setnames(sheet_df, note_cols[1], "cluster_note")
}
if (!("cluster_note" %in% names(sheet_df))) sheet_df[, cluster_note := NA_character_]

required_cols <- c("ID", "Taxa", "best_Z23", "best_B25", "best_B25_ibs", "expected_pair_ibs", "dna_diagnosis", "match_type")
missing_cols <- setdiff(required_cols, names(sheet_df))
if (length(missing_cols) > 0) stop("Missing columns in xlsx: ", paste(missing_cols, collapse = ", "))

sheet_df <- sheet_df[!is.na(ID) & ID != "" & ID != "ID"]
sheet_df[, best_B25_ibs := clean_numeric(best_B25_ibs)]
sheet_df[, expected_pair_ibs := clean_numeric(expected_pair_ibs)]

if ("best_Z23_ibs" %in% names(sheet_df)) {
  sheet_df[, best_Z23_ibs := clean_numeric(best_Z23_ibs)]
} else {
  sheet_df[, best_Z23_ibs := NA_real_]
}

if (!is.na(supp_file) && file.exists(supp_file)) {
  supp_df <- fread(supp_file, sep = "\t", header = TRUE, fill = TRUE, colClasses = "character")
  setnames(supp_df, clean_char(names(supp_df)))
  for (nm in names(supp_df)) supp_df[[nm]] <- clean_char(supp_df[[nm]])
  if ("best_Z23_ibs" %in% names(supp_df)) {
    if (!("ID" %in% names(supp_df))) stop("Supplement file must contain ID column when best_Z23_ibs is provided.")
    supp_df[, best_Z23_ibs := clean_numeric(best_Z23_ibs)]
    sheet_df <- merge(
      sheet_df,
      supp_df[, .(ID, best_Z23_ibs_supp = best_Z23_ibs)],
      by = "ID",
      all.x = TRUE,
      sort = FALSE
    )
    sheet_df[is.na(best_Z23_ibs) & !is.na(best_Z23_ibs_supp), best_Z23_ibs := best_Z23_ibs_supp]
    sheet_df[, best_Z23_ibs_supp := NULL]
  }
}

sheet_df <- infer_expected_ids(sheet_df)
sheet_df <- fill_from_ibs_if_available(sheet_df, sheet_name, ibs_dir, map_dir)
sheet_df[, row_index := seq_len(.N)]
sheet_df[, sample_label := factor(ID, levels = ID)]
sheet_df[, match_type := clean_char(match_type)]
sheet_df[is.na(match_type) | match_type == "", match_type := "NA"]

sheet_df[, z23_match_ok := !is.na(best_Z23) & best_Z23 == expected_Z23]
sheet_df[, b25_match_ok := !is.na(best_B25) & best_B25 == expected_B25]
sheet_df[, mismatch_flag := !(z23_match_ok & b25_match_ok)]
sheet_df[, delta_Z23 := best_Z23_ibs - expected_pair_ibs]
sheet_df[, delta_B25 := best_B25_ibs - expected_pair_ibs]
sheet_df[, delta_abs := abs(fifelse(is.na(delta_Z23), 0, delta_Z23)) + abs(fifelse(is.na(delta_B25), 0, delta_B25))]
sheet_df[, mismatch_type := fifelse(
  z23_match_ok & b25_match_ok, "Both_match",
  fifelse(!z23_match_ok & !b25_match_ok, "Both_mismatch", "One_side_mismatch")
)]
sheet_df[, expected_class := fifelse(
  dna_diagnosis %in% c("Exact_match", "low_ibs_match"), "match",
  fifelse(grepl("mismatch", dna_diagnosis, ignore.case = TRUE), "mismatch", "adjust")
)]
sheet_df[, match_type_plot := factor(match_type, levels = unique(c("match", "mismatch", "adjust", "NA", sort(setdiff(unique(match_type), c("match", "mismatch", "adjust", "NA"))))))]

write.table(
  sheet_df,
  file.path(outdir, paste0(sheet_name, "_bidirectional_scatter_input.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

msg("Rows loaded: %d", nrow(sheet_df))
msg("Rows with best_Z23_ibs available: %d", sum(!is.na(sheet_df$best_Z23_ibs)))
msg("Rows with best_B25_ibs available: %d", sum(!is.na(sheet_df$best_B25_ibs)))

bright_cols <- c(
  "match" = "#1F78B4",
  "mismatch" = "#E31A1C",
  "adjust" = "#33A02C"
)

match_type_cols_base <- c(
  "match" = "#1F78B4",
  "mismatch" = "#E31A1C",
  "adjust" = "#33A02C",
  "NA" = "#9E9E9E"
)

make_match_type_colors <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "NA"
  levels_use <- unique(c(names(match_type_cols_base), setdiff(unique(x), names(match_type_cols_base))))
  cols_use <- match_type_cols_base
  extra_levels <- setdiff(levels_use, names(cols_use))
  if (length(extra_levels) > 0) {
    extra_cols <- grDevices::colorRampPalette(c("#6B6B6B", "#CFCFCF"))(length(extra_levels))
    names(extra_cols) <- extra_levels
    cols_use <- c(cols_use, extra_cols)
  }
  list(levels = levels_use, cols = cols_use[levels_use])
}

match_type_info <- make_match_type_colors(sheet_df$match_type)
sheet_df[, match_type_plot := factor(match_type, levels = match_type_info$levels)]

plot_match_type_barplot <- function(df, title, outfile_base) {
  count_dt <- copy(df)[, .N, by = .(match_type = as.character(match_type_plot))]
  count_dt <- merge(
    data.table(match_type = levels(df$match_type_plot)),
    count_dt,
    by = "match_type",
    all.x = TRUE,
    sort = FALSE
  )
  count_dt[is.na(N), N := 0L]
  count_dt[, match_type_plot := factor(match_type, levels = levels(df$match_type_plot))]
  fwrite(count_dt[, .(match_type, N)], file.path(outdir, paste0(outfile_base, ".tsv")), sep = "\t")

  y_max <- max(count_dt$N, na.rm = TRUE)
  if (!is.finite(y_max) || y_max == 0) y_max <- 1

  p <- ggplot(count_dt, aes(x = match_type_plot, y = N, fill = match_type_plot)) +
    geom_col(width = 0.68, color = "grey20", linewidth = 0.25) +
    geom_text(aes(label = N), vjust = -0.35, size = 4.2, fontface = "bold") +
    scale_fill_manual(values = match_type_info$cols, drop = FALSE) +
    expand_limits(y = y_max * 1.16) +
    labs(
      x = "match_type",
      y = "Sample count",
      fill = "match_type",
      title = title
    ) +
    theme_pub(base_size = 13) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)
    )

  ggsave(file.path(outdir, paste0(outfile_base, ".pdf")), p, width = 6.4, height = 4.8)
  ggsave(file.path(outdir, paste0(outfile_base, ".png")), p, width = 6.4, height = 4.8, dpi = 300)
}

plot_taxa_match_type_barplot <- function(df, title, outfile_base) {
  plot_dt <- copy(df)
  plot_dt[, Taxa := clean_char(Taxa)]
  plot_dt[is.na(Taxa) | Taxa == "", Taxa := "Unknown"]
  plot_dt[, match_type := as.character(match_type_plot)]

  count_dt <- plot_dt[, .N, by = .(Taxa, match_type)]
  taxa_order <- plot_dt[, .N, by = Taxa][order(-N, Taxa)]$Taxa
  full_dt <- CJ(
    Taxa = taxa_order,
    match_type = levels(plot_dt$match_type_plot),
    unique = TRUE
  )
  count_dt <- merge(full_dt, count_dt, by = c("Taxa", "match_type"), all.x = TRUE, sort = FALSE)
  count_dt[is.na(N), N := 0L]
  count_dt[, Taxa := factor(Taxa, levels = taxa_order)]
  count_dt[, match_type_plot := factor(match_type, levels = levels(plot_dt$match_type_plot))]
  count_dt[, label := ifelse(N > 0, as.character(N), "")]

  count_wide <- dcast(count_dt[, .(Taxa = as.character(Taxa), match_type, N)], Taxa ~ match_type, value.var = "N", fill = 0)
  count_total <- count_dt[, .(total = sum(N)), by = .(Taxa = as.character(Taxa))]
  count_wide <- merge(count_total, count_wide, by = "Taxa", all.x = TRUE, sort = FALSE)

  fwrite(count_dt[, .(Taxa = as.character(Taxa), match_type, N)], file.path(outdir, paste0(outfile_base, "_long.tsv")), sep = "\t")
  fwrite(count_wide, file.path(outdir, paste0(outfile_base, "_wide.tsv")), sep = "\t")

  p <- ggplot(count_dt, aes(x = Taxa, y = N, fill = match_type_plot)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.25) +
    geom_text(
      aes(label = label),
      position = position_stack(vjust = 0.5),
      size = 3.4,
      fontface = "bold",
      color = "black"
    ) +
    scale_fill_manual(values = match_type_info$cols, drop = FALSE) +
    labs(
      x = "Taxa",
      y = "Sample count",
      fill = "match_type",
      title = title
    ) +
    theme_pub(base_size = 13) +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1)
    )

  plot_width <- max(7.2, min(14, 1.2 * length(taxa_order) + 3.5))
  ggsave(file.path(outdir, paste0(outfile_base, ".pdf")), p, width = plot_width, height = 5.4)
  ggsave(file.path(outdir, paste0(outfile_base, ".png")), p, width = plot_width, height = 5.4, dpi = 300)
}

plot_y_eq_x_scatter <- function(df, y_col, label_flag_col, best_col, title, y_label, outfile_base,
                                label_mode = "id_to_best", label_source_col = "ID") {
  plot_df <- copy(df)[!is.na(expected_pair_ibs) & !is.na(get(y_col))]
  if (nrow(plot_df) == 0) return(invisible(NULL))

  plot_df[, label_text := ifelse(
    get(label_flag_col),
    if (label_mode == "id") ID else paste0(get(label_source_col), "->", get(best_col)),
    ""
  )]
  label_df <- plot_df[get(label_flag_col) == TRUE]
  lim_low <- max(0.65, min(c(plot_df$expected_pair_ibs, plot_df[[y_col]]), na.rm = TRUE) - 0.01)
  lim_high <- min(1.05, max(1.02, max(c(plot_df$expected_pair_ibs, plot_df[[y_col]]), na.rm = TRUE) + 0.03))

  label_layer <- if (nrow(label_df) == 0) {
    NULL
  } else if (requireNamespace("ggrepel", quietly = TRUE)) {
    ggrepel::geom_text_repel(
      data = label_df,
      aes(label = label_text),
      color = "#E31A1C",
      size = 2.5,
      box.padding = 0.22,
      point.padding = 0.12,
      min.segment.length = 0,
      max.overlaps = 40,
      seed = 1,
      show.legend = FALSE
    )
  } else {
    geom_text(
      data = label_df,
      aes(label = label_text),
      color = "#E31A1C",
      size = 2.5,
      nudge_y = 0.006,
      check_overlap = TRUE,
      show.legend = FALSE
    )
  }

  p <- ggplot(plot_df, aes(x = expected_pair_ibs, y = get(y_col))) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.55, color = "grey35") +
    geom_vline(xintercept = 0.90, linetype = "dotted", linewidth = 0.4, color = "grey55") +
    geom_hline(yintercept = 0.90, linetype = "dotted", linewidth = 0.4, color = "grey55") +
    geom_vline(xintercept = 0.99, linetype = "dotted", linewidth = 0.4, color = "#B2182B") +
    geom_hline(yintercept = 0.99, linetype = "dotted", linewidth = 0.4, color = "#B2182B") +
    geom_point(aes(color = match_type_plot), size = 3, alpha = 0.86) +
    label_layer +
    scale_color_manual(values = match_type_info$cols, drop = FALSE) +
    coord_cartesian(xlim = c(lim_low, 1.02), ylim = c(lim_low, lim_high), clip = "off") +
    labs(
      x = "expected_pair_ibs",
      y = y_label,
      color = "match_type",
      title = title
    ) +
    theme_pub(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      plot.subtitle = element_blank(),
      plot.margin = margin(t = 12, r = 20, b = 16, l = 12)
    )

  ggsave(file.path(outdir, paste0(outfile_base, ".pdf")), p, width = 7.8, height = 6.2)
  ggsave(file.path(outdir, paste0(outfile_base, ".png")), p, width = 7.8, height = 6.2, dpi = 300)
}

plot_expected_pair_only <- function(df, title, outfile_base) {
  plot_df <- copy(df)[!is.na(expected_pair_ibs)]
  if (nrow(plot_df) == 0) return(invisible(NULL))

  plot_df[, plot_index := seq_len(.N)]
  plot_df[, label_text := ifelse(match_type == "mismatch", ID, "")]
  label_df <- plot_df[match_type == "mismatch"]
  y_low <- max(0.65, min(plot_df$expected_pair_ibs, na.rm = TRUE) - 0.02)

  label_layer <- if (nrow(label_df) == 0) {
    NULL
  } else if (requireNamespace("ggrepel", quietly = TRUE)) {
    ggrepel::geom_text_repel(
      data = label_df,
      aes(label = label_text),
      color = "#E31A1C",
      size = 2.7,
      box.padding = 0.25,
      point.padding = 0.15,
      min.segment.length = 0,
      max.overlaps = 80,
      seed = 2,
      show.legend = FALSE
    )
  } else {
    geom_text(
      data = label_df,
      aes(label = label_text),
      color = "#E31A1C",
      size = 2.7,
      nudge_y = 0.006,
      check_overlap = TRUE,
      show.legend = FALSE
    )
  }

  p <- ggplot(plot_df, aes(x = plot_index, y = expected_pair_ibs)) +
    geom_hline(yintercept = 0.90, linetype = "dotted", linewidth = 0.45, color = "grey55") +
    geom_hline(yintercept = 0.99, linetype = "dotted", linewidth = 0.55, color = "#B2182B") +
    geom_point(aes(color = match_type_plot), size = 3, alpha = 0.86) +
    label_layer +
    scale_color_manual(values = match_type_info$cols, drop = FALSE) +
    coord_cartesian(ylim = c(y_low, 1.02), clip = "off") +
    labs(
      x = "Sample order in table",
      y = "expected_pair_ibs",
      color = "match_type",
      title = title
    ) +
    theme_pub(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      plot.subtitle = element_blank(),
      plot.margin = margin(t = 12, r = 24, b = 16, l = 12)
    )

  ggsave(file.path(outdir, paste0(outfile_base, ".pdf")), p, width = 8.2, height = 5.8)
  ggsave(file.path(outdir, paste0(outfile_base, ".png")), p, width = 8.2, height = 5.8, dpi = 300)
}

plot_expected_pair_vs_best_b25 <- function(df, title, outfile_base) {
  plot_df <- copy(df)[!is.na(expected_pair_ibs) & !is.na(best_B25_ibs)]
  if (nrow(plot_df) == 0) return(invisible(NULL))

  plot_df[, match_type_plot := factor(match_type, levels = match_type_info$levels)]
  plot_df[, label_text := ifelse(match_type == "mismatch", ID, "")]
  label_df <- plot_df[label_text != ""]
  lim_low <- max(0.65, min(c(plot_df$expected_pair_ibs, plot_df$best_B25_ibs), na.rm = TRUE) - 0.01)
  lim_high <- min(1.05, max(1.02, max(c(plot_df$expected_pair_ibs, plot_df$best_B25_ibs), na.rm = TRUE) + 0.03))

  label_layer <- if (nrow(label_df) == 0) {
    NULL
  } else if (requireNamespace("ggrepel", quietly = TRUE)) {
    ggrepel::geom_text_repel(
      data = label_df,
      aes(label = label_text),
      color = "#E31A1C",
      size = 2.7,
      box.padding = 0.25,
      point.padding = 0.15,
      min.segment.length = 0,
      max.overlaps = 80,
      seed = 3,
      show.legend = FALSE
    )
  } else {
    geom_text(
      data = label_df,
      aes(label = label_text),
      color = "#E31A1C",
      size = 2.7,
      nudge_y = 0.006,
      check_overlap = TRUE,
      show.legend = FALSE
    )
  }

  p <- ggplot(plot_df, aes(x = expected_pair_ibs, y = best_B25_ibs)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.55, color = "grey35") +
    geom_vline(xintercept = 0.90, linetype = "dotted", linewidth = 0.4, color = "grey55") +
    geom_hline(yintercept = 0.90, linetype = "dotted", linewidth = 0.4, color = "grey55") +
    geom_vline(xintercept = 0.99, linetype = "dotted", linewidth = 0.4, color = "#B2182B") +
    geom_hline(yintercept = 0.99, linetype = "dotted", linewidth = 0.4, color = "#B2182B") +
    geom_point(aes(color = match_type_plot), size = 3, alpha = 0.86) +
    label_layer +
    scale_color_manual(values = match_type_info$cols, drop = FALSE) +
    coord_cartesian(xlim = c(lim_low, 1.02), ylim = c(lim_low, lim_high), clip = "off") +
    labs(
      x = "Expected pair IBS",
      y = "Best B25 IBS",
      color = "match_type",
      title = title,
      subtitle = "Dashed line: y = x; dotted lines: IBS = 0.90 and 0.99"
    ) +
    theme_pub(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey25"),
      plot.margin = margin(t = 12, r = 24, b = 16, l = 12)
    )

  ggsave(file.path(outdir, paste0(outfile_base, ".pdf")), p, width = 7.8, height = 6.2)
  ggsave(file.path(outdir, paste0(outfile_base, ".png")), p, width = 7.8, height = 6.2, dpi = 300)
}

sheet_df[, b25_best_mismatch := !b25_match_ok]
sheet_df[, z23_best_mismatch := !z23_match_ok]
sheet_df[, best_top1_ibs := pmax(best_Z23_ibs, best_B25_ibs, na.rm = TRUE)]
sheet_df[is.infinite(best_top1_ibs), best_top1_ibs := NA_real_]
sheet_df[, any_best_mismatch := z23_best_mismatch | b25_best_mismatch]
sheet_df[, best_top1_label := fifelse(
  z23_best_mismatch & b25_best_mismatch, paste0(best_Z23, "/", best_B25),
  fifelse(z23_best_mismatch, best_Z23, fifelse(b25_best_mismatch, best_B25, ""))
)]

plot_match_type_barplot(
  sheet_df,
  title = paste0(sheet_name, ": match_type summary"),
  outfile_base = paste0(sheet_name, "_00_match_type_barplot")
)

plot_taxa_match_type_barplot(
  sheet_df,
  title = paste0(sheet_name, ": match_type summary by Taxa"),
  outfile_base = paste0(sheet_name, "_00b_taxa_match_type_barplot")
)

plot_y_eq_x_scatter(
  sheet_df,
  y_col = "best_B25_ibs",
  label_flag_col = "b25_best_mismatch",
  best_col = "best_B25",
  title = paste0(sheet_name, ": Z23 to B25 best IBS"),
  y_label = "best_B25_ibs",
  outfile_base = paste0(sheet_name, "_01_Z23_to_B25_best_ibs"),
  label_source_col = "ID"
)

plot_y_eq_x_scatter(
  sheet_df,
  y_col = "best_Z23_ibs",
  label_flag_col = "z23_best_mismatch",
  best_col = "best_Z23",
  title = paste0(sheet_name, ": B25 to Z23 best IBS"),
  y_label = "best_Z23_ibs",
  outfile_base = paste0(sheet_name, "_02_B25_to_Z23_best_ibs"),
  label_source_col = "expected_B25"
)

sheet_df[, expected_pair_mismatch := match_type == "mismatch"]

plot_y_eq_x_scatter(
  sheet_df,
  y_col = "best_top1_ibs",
  label_flag_col = "expected_pair_mismatch",
  best_col = "best_top1_label",
  title = paste0(sheet_name, ": Expected pair vs best IBS"),
  y_label = "best top1 IBS",
  outfile_base = paste0(sheet_name, "_03_expected_pair_vs_best_ibs"),
  label_mode = "id"
)

plot_expected_pair_only(
  sheet_df,
  title = paste0(sheet_name, ": Expected pair IBS only"),
  outfile_base = paste0(sheet_name, "_04_expected_pair_ibs_only")
)

plot_expected_pair_vs_best_b25(
  sheet_df,
  title = paste0(sheet_name, ": Expected vs best-match IBS"),
  outfile_base = paste0(sheet_name, "_05_expected_pair_vs_best_B25_ibs")
)

# Bidirectional check table
check_dt <- copy(sheet_df)[, .(
  ID,
  Taxa,
  expected_Z23,
  best_Z23,
  best_Z23_ibs,
  z23_match_ok,
  expected_B25,
  best_B25,
  best_B25_ibs,
  b25_match_ok,
  delta_Z23,
  delta_B25,
  delta_abs,
  expected_pair_ibs,
  expected_class,
  match_type,
  dna_diagnosis,
  cluster_note,
  mismatch_type
)]

check_dt[, bidirectional_status := fifelse(
  z23_match_ok & b25_match_ok, "Both_match",
  fifelse(!z23_match_ok & !b25_match_ok, "Both_mismatch", "One_side_mismatch")
)]

write.table(
  check_dt,
  file.path(outdir, paste0(sheet_name, "_bidirectional_check_table.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

summary_dt <- check_dt[, .N, by = .(bidirectional_status, match_type, dna_diagnosis)]
write.table(
  summary_dt,
  file.path(outdir, paste0(sheet_name, "_bidirectional_check_summary.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

abnormal_dt <- check_dt[
  bidirectional_status != "Both_match" | (is.finite(delta_abs) & delta_abs > 0),
  .(
    ID,
    Taxa,
    best_Z23,
    best_Z23_ibs,
    delta_Z23,
    best_B25,
    best_B25_ibs,
    delta_B25,
    delta_abs,
    expected_pair_ibs,
    match_type,
    dna_diagnosis,
    bidirectional_status,
    cluster_note
  )
]
if (nrow(abnormal_dt) > 0) {
  setorder(abnormal_dt, -delta_abs, ID)
  write.table(
    abnormal_dt,
    file.path(outdir, paste0(sheet_name, "_abnormal_samples_simple.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

msg("Output dir: %s", outdir)
