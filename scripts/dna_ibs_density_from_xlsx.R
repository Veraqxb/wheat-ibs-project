#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(data.table)
  library(ggplot2)
  library(scales)
})

# =========================
# Input
# =========================
xlsx_file <- "/Users/veraqiu/Desktop/CAMP_samples_2026/CAMP_DNA_check_260423.xlsx"
sheets_use <- c("C2", "C4", "C6")
outdir <- "/Users/veraqiu/Desktop/CAMP_samples_2026/dna_diagnosis_summary/ibs_density"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# =========================
# Diagnosis order and colors
# =========================
diag_levels <- c(
  "Exact_match",
  "low_ibs_match",
  "NA",
  "B25_only_exact",
  "Z23_only_exact",
  "Z23_cluster",
  "True_mismatch"
)

diag_cols <- c(
  "Exact_match" = "#0072B2",
  "low_ibs_match" = "#CC79A7",
  "NA" = "#BDBDBD",
  "B25_only_exact" = "#009E73",
  "Z23_only_exact" = "#56B4E9",
  "Z23_cluster" = "#E69F00",
  "True_mismatch" = "#D55E00"
)

metric_cols <- c(
  "expected_pair_ibs" = "#4E79A7",
  "best_B25_ibs" = "#D55E00"
)

theme_pub <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey30"),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      axis.line = element_line(linewidth = 0.5),
      axis.ticks = element_line(linewidth = 0.5),
      legend.key.size = unit(0.55, "cm")
    )
}

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

clean_diagnosis <- function(x) {
  x <- clean_char(x)
  x[is.na(x) | x %in% c("", "N/A", "NaN", "NULL", "null", "-")] <- "NA"
  x
}

prepare_diag_levels_cols <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "NA"
  extra_levels <- setdiff(unique(x), diag_levels)
  diag_levels_all <- unique(c(diag_levels, extra_levels))

  diag_cols_use <- diag_cols
  if (length(extra_levels) > 0) {
    extra_cols <- rep("#8C8C8C", length(extra_levels))
    names(extra_cols) <- extra_levels
    diag_cols_use <- c(diag_cols_use, extra_cols)
  }

  list(levels = diag_levels_all, cols = diag_cols_use)
}

safe_quantile <- function(x, probs) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(rep(NA_real_, length(probs)))
  as.numeric(quantile(x, probs = probs, na.rm = TRUE, names = FALSE))
}

read_one_sheet <- function(file, sheet_name) {
  dt <- as.data.table(suppressMessages(read_excel(file, sheet = sheet_name, col_names = TRUE)))
  setnames(dt, clean_char(names(dt)))

  keep_cols <- c("ID", "Taxa", "best_Z23", "best_B25", "best_B25_ibs", "expected_pair_ibs", "dna_diagnosis")
  keep_cols <- intersect(keep_cols, names(dt))
  dt <- dt[, ..keep_cols]

  for (j in names(dt)) {
    dt[[j]] <- clean_char(dt[[j]])
  }

  if ("ID" %in% names(dt)) {
    dt <- dt[!is.na(ID) & ID != "" & ID != "ID"]
  }
  if ("best_B25_ibs" %in% names(dt)) {
    dt <- dt[is.na(best_B25_ibs) | best_B25_ibs != "best_B25_ibs"]
    dt[, best_B25_ibs := clean_numeric(best_B25_ibs)]
  } else {
    dt[, best_B25_ibs := NA_real_]
  }
  if ("expected_pair_ibs" %in% names(dt)) {
    dt <- dt[is.na(expected_pair_ibs) | expected_pair_ibs != "expected_pair_ibs"]
    dt[, expected_pair_ibs := clean_numeric(expected_pair_ibs)]
  } else {
    dt[, expected_pair_ibs := NA_real_]
  }
  if ("dna_diagnosis" %in% names(dt)) {
    dt <- dt[is.na(dna_diagnosis) | dna_diagnosis != "dna_diagnosis"]
    dt[, dna_diagnosis := clean_diagnosis(dna_diagnosis)]
  } else {
    dt[, dna_diagnosis := "NA"]
  }

  if (!("Taxa" %in% names(dt))) {
    dt[, Taxa := "Unknown"]
  }
  dt[is.na(Taxa) | Taxa == "", Taxa := "Unknown"]
  dt[, ploidy := sheet_name]
  dt[, delta := best_B25_ibs - expected_pair_ibs]
  dt[, is_complete := !is.na(best_B25_ibs) & !is.na(expected_pair_ibs)]
  dt
}

build_long_ibs <- function(dt) {
  rbindlist(list(
    dt[!is.na(expected_pair_ibs), .(
      ID, Taxa, ploidy, dna_diagnosis, metric = "expected_pair_ibs", ibs = expected_pair_ibs
    )],
    dt[!is.na(best_B25_ibs), .(
      ID, Taxa, ploidy, dna_diagnosis, metric = "best_B25_ibs", ibs = best_B25_ibs
    )]
  ), fill = TRUE)
}

make_ibs_summary <- function(long_dt) {
  probs <- c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1)
  long_dt[, {
    qv <- safe_quantile(ibs, probs)
    .(
      n = .N,
      min = qv[1],
      q01 = qv[2],
      q05 = qv[3],
      q25 = qv[4],
      median = qv[5],
      q75 = qv[6],
      q95 = qv[7],
      q99 = qv[8],
      max = qv[9],
      mean = mean(ibs, na.rm = TRUE),
      sd = sd(ibs, na.rm = TRUE)
    )
  }, by = .(ploidy, metric, dna_diagnosis)]
}

plot_density_overlay <- function(long_dt, group_name, outdir_group) {
  plot_dt <- copy(long_dt)
  if (nrow(plot_dt) == 0) return(invisible(NULL))

  p <- ggplot(plot_dt, aes(x = ibs, color = metric, fill = metric)) +
    geom_density(alpha = 0.12, linewidth = 1) +
    geom_vline(xintercept = 0.90, linetype = "dashed", linewidth = 0.55, color = "black") +
    geom_vline(xintercept = 0.99, linetype = "dashed", linewidth = 0.55, color = "#B2182B") +
    scale_color_manual(values = metric_cols) +
    scale_fill_manual(values = metric_cols) +
    coord_cartesian(xlim = c(max(0.70, min(plot_dt$ibs, na.rm = TRUE) - 0.01), 1.00)) +
    labs(
      x = "IBS",
      y = "Density",
      color = "Metric",
      fill = "Metric",
      title = paste0(group_name, ": IBS density distribution"),
      subtitle = "Blue = expected pair IBS; orange = best B25 IBS"
    ) +
    theme_pub()

  ggsave(file.path(outdir_group, paste0(group_name, "_ibs_density_overlay.pdf")), p, width = 8.5, height = 5.6)
  ggsave(file.path(outdir_group, paste0(group_name, "_ibs_density_overlay.png")), p, width = 8.5, height = 5.6, dpi = 300)
}

plot_ecdf_overlay <- function(long_dt, group_name, outdir_group) {
  plot_dt <- copy(long_dt)
  if (nrow(plot_dt) == 0) return(invisible(NULL))

  p <- ggplot(plot_dt, aes(x = ibs, color = metric)) +
    stat_ecdf(linewidth = 1) +
    geom_vline(xintercept = 0.90, linetype = "dashed", linewidth = 0.55, color = "black") +
    geom_vline(xintercept = 0.99, linetype = "dashed", linewidth = 0.55, color = "#B2182B") +
    scale_color_manual(values = metric_cols) +
    coord_cartesian(xlim = c(max(0.70, min(plot_dt$ibs, na.rm = TRUE) - 0.01), 1.00)) +
    labs(
      x = "IBS",
      y = "Cumulative proportion",
      color = "Metric",
      title = paste0(group_name, ": IBS cumulative distribution"),
      subtitle = "Threshold guides at 0.90 and 0.99"
    ) +
    theme_pub()

  ggsave(file.path(outdir_group, paste0(group_name, "_ibs_ecdf_overlay.pdf")), p, width = 8.5, height = 5.6)
  ggsave(file.path(outdir_group, paste0(group_name, "_ibs_ecdf_overlay.png")), p, width = 8.5, height = 5.6, dpi = 300)
}

plot_density_by_diagnosis <- function(long_dt, group_name, outdir_group, diag_levels_all, diag_cols_use) {
  plot_dt <- long_dt[metric == "expected_pair_ibs"]
  if (nrow(plot_dt) == 0) return(invisible(NULL))
  plot_dt[, dna_diagnosis := factor(as.character(dna_diagnosis), levels = diag_levels_all)]

  p <- ggplot(plot_dt, aes(x = ibs, color = dna_diagnosis, fill = dna_diagnosis)) +
    geom_density(alpha = 0.10, linewidth = 0.9) +
    geom_vline(xintercept = 0.90, linetype = "dashed", linewidth = 0.55, color = "black") +
    geom_vline(xintercept = 0.99, linetype = "dashed", linewidth = 0.55, color = "#B2182B") +
    scale_color_manual(values = diag_cols_use, drop = FALSE) +
    scale_fill_manual(values = diag_cols_use, drop = FALSE) +
    coord_cartesian(xlim = c(max(0.70, min(plot_dt$ibs, na.rm = TRUE) - 0.01), 1.00)) +
    labs(
      x = "Expected pair IBS",
      y = "Density",
      color = "DNA diagnosis",
      fill = "DNA diagnosis",
      title = paste0(group_name, ": Expected-pair IBS by diagnosis")
    ) +
    theme_pub()

  ggsave(file.path(outdir_group, paste0(group_name, "_expected_pair_density_by_diagnosis.pdf")), p, width = 9.2, height = 5.8)
  ggsave(file.path(outdir_group, paste0(group_name, "_expected_pair_density_by_diagnosis.png")), p, width = 9.2, height = 5.8, dpi = 300)
}

plot_boxplot_by_metric <- function(long_dt, group_name, outdir_group) {
  plot_dt <- copy(long_dt)
  if (nrow(plot_dt) == 0) return(invisible(NULL))
  plot_dt[, metric := factor(metric, levels = c("expected_pair_ibs", "best_B25_ibs"))]

  p <- ggplot(plot_dt, aes(x = metric, y = ibs, fill = metric)) +
    geom_boxplot(outlier.size = 0.5, width = 0.68) +
    geom_hline(yintercept = 0.90, linetype = "dashed", linewidth = 0.55, color = "black") +
    geom_hline(yintercept = 0.99, linetype = "dashed", linewidth = 0.55, color = "#B2182B") +
    scale_fill_manual(values = metric_cols) +
    coord_cartesian(ylim = c(max(0.70, min(plot_dt$ibs, na.rm = TRUE) - 0.01), 1.00)) +
    labs(
      x = "",
      y = "IBS",
      fill = "Metric",
      title = paste0(group_name, ": IBS boxplot comparison")
    ) +
    theme_pub()

  ggsave(file.path(outdir_group, paste0(group_name, "_ibs_boxplot.pdf")), p, width = 6.4, height = 5.5)
  ggsave(file.path(outdir_group, paste0(group_name, "_ibs_boxplot.png")), p, width = 6.4, height = 5.5, dpi = 300)
}

plot_scatter_expected_vs_best <- function(dt, group_name, outdir_group, diag_levels_all, diag_cols_use) {
  plot_dt <- copy(dt)[is_complete == TRUE]
  if (nrow(plot_dt) == 0) return(invisible(NULL))
  plot_dt[, dna_diagnosis := factor(as.character(dna_diagnosis), levels = diag_levels_all)]

  p <- ggplot(plot_dt, aes(x = expected_pair_ibs, y = best_B25_ibs, color = dna_diagnosis)) +
    geom_point(size = 2.8, alpha = 0.85) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.5, color = "grey45") +
    geom_vline(xintercept = 0.90, linetype = "dotted", linewidth = 0.5, color = "grey35") +
    geom_hline(yintercept = 0.90, linetype = "dotted", linewidth = 0.5, color = "grey35") +
    geom_vline(xintercept = 0.99, linetype = "dotted", linewidth = 0.5, color = "#B2182B") +
    geom_hline(yintercept = 0.99, linetype = "dotted", linewidth = 0.5, color = "#B2182B") +
    scale_color_manual(values = diag_cols_use, drop = FALSE) +
    coord_cartesian(
      xlim = c(max(0.70, min(plot_dt$expected_pair_ibs, na.rm = TRUE) - 0.01), 1.00),
      ylim = c(max(0.70, min(plot_dt$best_B25_ibs, na.rm = TRUE) - 0.01), 1.00)
    ) +
    labs(
      x = "Expected pair IBS",
      y = "Best B25 IBS",
      color = "DNA diagnosis",
      title = paste0(group_name, ": Expected pair vs best B25")
    ) +
    theme_pub()

  ggsave(file.path(outdir_group, paste0(group_name, "_expected_vs_best_scatter.pdf")), p, width = 7.6, height = 6.0)
  ggsave(file.path(outdir_group, paste0(group_name, "_expected_vs_best_scatter.png")), p, width = 7.6, height = 6.0, dpi = 300)
}

plot_one_group <- function(dt, group_name) {
  outdir_group <- file.path(outdir, group_name)
  dir.create(outdir_group, recursive = TRUE, showWarnings = FALSE)

  level_info <- prepare_diag_levels_cols(dt$dna_diagnosis)
  diag_levels_all <- level_info$levels
  diag_cols_use <- level_info$cols

  long_dt <- build_long_ibs(dt)
  fwrite(dt, file.path(outdir_group, paste0(group_name, "_cleaned.tsv")), sep = "\t", quote = FALSE, na = "NA")
  fwrite(long_dt, file.path(outdir_group, paste0(group_name, "_ibs_long.tsv")), sep = "\t", quote = FALSE, na = "NA")
  fwrite(make_ibs_summary(long_dt), file.path(outdir_group, paste0(group_name, "_ibs_density_summary.tsv")), sep = "\t", quote = FALSE, na = "NA")

  plot_density_overlay(long_dt, group_name, outdir_group)
  plot_ecdf_overlay(long_dt, group_name, outdir_group)
  plot_density_by_diagnosis(long_dt, group_name, outdir_group, diag_levels_all, diag_cols_use)
  plot_boxplot_by_metric(long_dt, group_name, outdir_group)
  plot_scatter_expected_vs_best(dt, group_name, outdir_group, diag_levels_all, diag_cols_use)
}

all_list <- list()

for (sheet_name in sheets_use) {
  message("Processing sheet: ", sheet_name)
  dt <- read_one_sheet(xlsx_file, sheet_name)
  all_list[[sheet_name]] <- dt
  plot_one_group(dt, sheet_name)
}

all_dt <- rbindlist(all_list, fill = TRUE)
if (nrow(all_dt) > 0) {
  all_long <- build_long_ibs(all_dt)
  level_info <- prepare_diag_levels_cols(all_dt$dna_diagnosis)
  diag_levels_all <- level_info$levels
  diag_cols_use <- level_info$cols

  fwrite(all_dt, file.path(outdir, "all_ploidy_cleaned.tsv"), sep = "\t", quote = FALSE, na = "NA")
  fwrite(all_long, file.path(outdir, "all_ploidy_ibs_long.tsv"), sep = "\t", quote = FALSE, na = "NA")
  fwrite(make_ibs_summary(all_long), file.path(outdir, "all_ploidy_ibs_density_summary.tsv"), sep = "\t", quote = FALSE, na = "NA")

  p_all_density <- ggplot(all_long, aes(x = ibs, color = metric, fill = metric)) +
    geom_density(alpha = 0.12, linewidth = 1) +
    geom_vline(xintercept = 0.90, linetype = "dashed", linewidth = 0.55, color = "black") +
    geom_vline(xintercept = 0.99, linetype = "dashed", linewidth = 0.55, color = "#B2182B") +
    scale_color_manual(values = metric_cols) +
    scale_fill_manual(values = metric_cols) +
    facet_wrap(~ploidy, nrow = 1, scales = "fixed") +
    coord_cartesian(xlim = c(max(0.70, min(all_long$ibs, na.rm = TRUE) - 0.01), 1.00)) +
    labs(
      x = "IBS",
      y = "Density",
      color = "Metric",
      fill = "Metric",
      title = "CAMP IBS density distributions across C2, C4 and C6"
    ) +
    theme_pub()

  ggsave(file.path(outdir, "all_ploidy_ibs_density_faceted.pdf"), p_all_density, width = 12.2, height = 5.6)
  ggsave(file.path(outdir, "all_ploidy_ibs_density_faceted.png"), p_all_density, width = 12.2, height = 5.6, dpi = 300)

  p_all_ecdf <- ggplot(all_long, aes(x = ibs, color = metric)) +
    stat_ecdf(linewidth = 1) +
    geom_vline(xintercept = 0.90, linetype = "dashed", linewidth = 0.55, color = "black") +
    geom_vline(xintercept = 0.99, linetype = "dashed", linewidth = 0.55, color = "#B2182B") +
    scale_color_manual(values = metric_cols) +
    facet_wrap(~ploidy, nrow = 1, scales = "fixed") +
    coord_cartesian(xlim = c(max(0.70, min(all_long$ibs, na.rm = TRUE) - 0.01), 1.00)) +
    labs(
      x = "IBS",
      y = "Cumulative proportion",
      color = "Metric",
      title = "CAMP IBS cumulative distributions across C2, C4 and C6"
    ) +
    theme_pub()

  ggsave(file.path(outdir, "all_ploidy_ibs_ecdf_faceted.pdf"), p_all_ecdf, width = 12.2, height = 5.6)
  ggsave(file.path(outdir, "all_ploidy_ibs_ecdf_faceted.png"), p_all_ecdf, width = 12.2, height = 5.6, dpi = 300)

  p_diag <- ggplot(
    all_long[metric == "expected_pair_ibs"],
    aes(x = ibs, color = factor(dna_diagnosis, levels = diag_levels_all), fill = factor(dna_diagnosis, levels = diag_levels_all))
  ) +
    geom_density(alpha = 0.10, linewidth = 0.9) +
    geom_vline(xintercept = 0.90, linetype = "dashed", linewidth = 0.55, color = "black") +
    geom_vline(xintercept = 0.99, linetype = "dashed", linewidth = 0.55, color = "#B2182B") +
    scale_color_manual(values = diag_cols_use, drop = FALSE) +
    scale_fill_manual(values = diag_cols_use, drop = FALSE) +
    facet_wrap(~ploidy, nrow = 1, scales = "fixed") +
    coord_cartesian(xlim = c(max(0.70, min(all_long$ibs, na.rm = TRUE) - 0.01), 1.00)) +
    labs(
      x = "Expected pair IBS",
      y = "Density",
      color = "DNA diagnosis",
      fill = "DNA diagnosis",
      title = "Expected-pair IBS by diagnosis across ploidy groups"
    ) +
    theme_pub()

  ggsave(file.path(outdir, "all_ploidy_expected_pair_density_by_diagnosis.pdf"), p_diag, width = 12.5, height = 5.8)
  ggsave(file.path(outdir, "all_ploidy_expected_pair_density_by_diagnosis.png"), p_diag, width = 12.5, height = 5.8, dpi = 300)
}

cat("Done.\n")
cat("Output dir:", outdir, "\n")
