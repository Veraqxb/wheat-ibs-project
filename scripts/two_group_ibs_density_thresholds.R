#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
if (is.na(script_arg)) stop("Cannot determine script path for two_group_ibs_density_thresholds.R")
script_dirname <- dirname(normalizePath(sub("^--file=", "", script_arg)))
source(file.path(script_dirname, "ibs_common.R"))

opt <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("mibs", "id", "map", "outdir")
missing <- required[!required %in% names(opt)]
if (length(missing) > 0) stop("Missing arguments: ", paste(missing, collapse = ", "))

prefix <- opt[["prefix"]] %||% sub("\\.mibs$", "", basename(opt[["mibs"]]))
anchor_col <- opt[["anchor-col"]] %||% "Z23"
secondary_col <- opt[["secondary-col"]] %||% "B25"
xmin <- as.numeric(opt[["xmin"]] %||% "0.70")
xmax <- as.numeric(opt[["xmax"]] %||% "1.00")

dir.create(opt[["outdir"]], recursive = TRUE, showWarnings = FALSE)

msg <- function(fmt, ...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%F %T"), sprintf(fmt, ...)))
}

get_anchor_col <- function(map_df, preferred_col) {
  if (preferred_col %in% names(map_df)) return(preferred_col)
  if ("Z23" %in% names(map_df)) return("Z23")
  if ("CAMP编号" %in% names(map_df)) return("CAMP编号")
  stop("No anchor/Z23-like column found in map file")
}

extract_within_values <- function(mat, ids, label) {
  ids <- unique(standardize_id(ids))
  ids <- ids[!is.na(ids) & ids %in% rownames(mat)]
  if (length(ids) < 2) {
    return(data.table(type = character(), ibs = numeric()))
  }

  sub_m <- mat[ids, ids, drop = FALSE]
  vals <- sub_m[upper.tri(sub_m)]
  vals <- vals[!is.na(vals)]
  data.table(type = label, ibs = vals)
}

extract_between_values <- function(mat, row_ids, col_ids, label) {
  row_ids <- unique(standardize_id(row_ids))
  col_ids <- unique(standardize_id(col_ids))
  row_ids <- row_ids[!is.na(row_ids) & row_ids %in% rownames(mat)]
  col_ids <- col_ids[!is.na(col_ids) & col_ids %in% colnames(mat)]

  if (length(row_ids) == 0 || length(col_ids) == 0) {
    return(data.table(type = character(), ibs = numeric()))
  }

  vals <- as.vector(mat[row_ids, col_ids, drop = FALSE])
  vals <- vals[!is.na(vals)]
  data.table(type = label, ibs = vals)
}

extract_expected_pairs <- function(mat, map_df, anchor_col, secondary_col) {
  if (!(anchor_col %in% names(map_df)) || !(secondary_col %in% names(map_df))) {
    return(data.table())
  }

  anchor <- standardize_id(map_df[[anchor_col]])
  secondary <- standardize_id(map_df[[secondary_col]])
  keep <- !is.na(anchor) & !is.na(secondary)
  anchor <- anchor[keep]
  secondary <- secondary[keep]

  ibs_vals <- mapply(function(a, b) {
    if (a %in% colnames(mat) && b %in% rownames(mat)) {
      as.numeric(mat[b, a])
    } else {
      NA_real_
    }
  }, anchor, secondary)

  anchor_self_vals <- vapply(anchor, function(a) {
    if (a %in% rownames(mat) && a %in% colnames(mat)) as.numeric(mat[a, a]) else NA_real_
  }, numeric(1))

  secondary_self_vals <- vapply(secondary, function(b) {
    if (b %in% rownames(mat) && b %in% colnames(mat)) as.numeric(mat[b, b]) else NA_real_
  }, numeric(1))

  data.table(
    pair_index = seq_along(anchor),
    expected_anchor = anchor,
    expected_secondary = secondary,
    ibs = ibs_vals,
    anchor_self_ibs = anchor_self_vals,
    secondary_self_ibs = secondary_self_vals
  )
}

safe_quantile <- function(x, probs) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(rep(NA_real_, length(probs)))
  as.numeric(quantile(x, probs = probs, na.rm = TRUE, names = FALSE))
}

make_summary <- function(dt) {
  probs <- c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1)
  dt[, {
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
      max = qv[9]
    )
  }, by = type]
}

make_threshold_recommendation <- function(expected_vals, between_vals, anchor_vals, secondary_vals) {
  exp_q01 <- safe_quantile(expected_vals, 0.01)
  exp_q05 <- safe_quantile(expected_vals, 0.05)
  bg_q95 <- safe_quantile(between_vals, 0.95)
  bg_q99 <- safe_quantile(between_vals, 0.99)
  anchor_q99 <- safe_quantile(anchor_vals, 0.99)
  secondary_q99 <- safe_quantile(secondary_vals, 0.99)

  data.table(
    metric = c(
      "expected_q01",
      "expected_q05",
      "between_bg_q95",
      "between_bg_q99",
      paste0(tolower(anchor_col), "_within_q99"),
      paste0(tolower(secondary_col), "_within_q99"),
      "suggested_loose_lower",
      "suggested_strict_lower",
      "separation_gap_q05_vs_bgq99"
    ),
    value = c(
      exp_q01,
      exp_q05,
      bg_q95,
      bg_q99,
      anchor_q99,
      secondary_q99,
      max(0.90, bg_q99, na.rm = TRUE),
      0.99,
      exp_q05 - bg_q99
    )
  )
}

type_cols <- c(
  "Anchor_within" = "#4E79A7",
  "Secondary_within" = "#59A14F",
  "Between_all" = "#E15759",
  "Expected_pair" = "#F28E2B"
)

line_cols <- c(
  "Fixed_0.90" = "black",
  "Background_q99" = "#7B3294",
  "Strict_0.99" = "#B2182B"
)

line_types <- c(
  "Fixed_0.90" = "dashed",
  "Background_q99" = "dotdash",
  "Strict_0.99" = "dashed"
)

msg("Reading files...")
mat <- read_mibs_matrix(opt[["mibs"]], opt[["id"]])
map_df <- read_sample_map(opt[["map"]])

anchor_col <- get_anchor_col(map_df, anchor_col)
if (!(secondary_col %in% names(map_df))) stop("Missing secondary column in map: ", secondary_col)

anchor_ids <- standardize_id(map_df[[anchor_col]])
secondary_ids <- standardize_id(map_df[[secondary_col]])

msg("Extracting IBS distributions...")
anchor_within <- extract_within_values(mat, anchor_ids, "Anchor_within")
secondary_within <- extract_within_values(mat, secondary_ids, "Secondary_within")
between_all <- extract_between_values(mat, secondary_ids, anchor_ids, "Between_all")

expected_dt <- extract_expected_pairs(mat, map_df, anchor_col, secondary_col)
expected_dist <- expected_dt[!is.na(ibs), .(type = "Expected_pair", ibs)]

plot_dt <- rbindlist(list(anchor_within, secondary_within, between_all, expected_dist), fill = TRUE)
plot_dt[, type := factor(type, levels = names(type_cols))]

fwrite(plot_dt, file.path(opt[["outdir"]], paste0(prefix, "_ibs_density_input_values.tsv")), sep = "\t")
fwrite(expected_dt, file.path(opt[["outdir"]], paste0(prefix, "_expected_pair_ibs.tsv")), sep = "\t")

summary_dt <- make_summary(plot_dt)
fwrite(summary_dt, file.path(opt[["outdir"]], paste0(prefix, "_ibs_density_summary.tsv")), sep = "\t")

threshold_dt <- make_threshold_recommendation(
  expected_vals = expected_dist$ibs,
  between_vals = between_all$ibs,
  anchor_vals = anchor_within$ibs,
  secondary_vals = secondary_within$ibs
)
fwrite(threshold_dt, file.path(opt[["outdir"]], paste0(prefix, "_ibs_threshold_recommendation.tsv")), sep = "\t")

bg_q99 <- threshold_dt[metric == "between_bg_q99", value]
loose_thr <- threshold_dt[metric == "suggested_loose_lower", value]
strict_thr <- threshold_dt[metric == "suggested_strict_lower", value]

line_dt <- data.table(
  x = c(0.90, bg_q99, 0.99),
  line_type = factor(c("Fixed_0.90", "Background_q99", "Strict_0.99"), levels = names(line_cols)),
  label = c("0.90", paste0("BG q99 = ", sprintf("%.4f", bg_q99)), "0.99")
)
line_dt <- line_dt[!is.na(x)]

p1 <- ggplot(plot_dt, aes(x = ibs, color = type, fill = type)) +
  geom_density(alpha = 0.12, linewidth = 1, na.rm = TRUE) +
  geom_vline(
    data = line_dt,
    aes(xintercept = x, linetype = line_type, color = line_type),
    linewidth = 0.7,
    show.legend = TRUE
  ) +
  scale_color_manual(values = c(type_cols, line_cols), breaks = c(names(type_cols), names(line_cols))) +
  scale_fill_manual(values = type_cols) +
  scale_linetype_manual(values = line_types) +
  coord_cartesian(xlim = c(xmin, xmax)) +
  theme_classic(base_size = 13) +
  labs(
    x = "IBS",
    y = "Density",
    color = "Distribution / Threshold",
    fill = "Distribution",
    linetype = "Threshold",
    title = paste0(prefix, " 2group IBS density distributions"),
    subtitle = paste0("Background-aware threshold (between-group q99) = ", sprintf("%.4f", bg_q99), "; strict threshold = 0.99")
  )

ggsave(file.path(opt[["outdir"]], paste0(prefix, "_ibs_density_distributions.pdf")), p1, width = 9.5, height = 5.8)
ggsave(file.path(opt[["outdir"]], paste0(prefix, "_ibs_density_distributions.png")), p1, width = 9.5, height = 5.8, dpi = 300)

p2 <- ggplot(plot_dt, aes(x = ibs, color = type)) +
  stat_ecdf(linewidth = 1, na.rm = TRUE) +
  geom_vline(
    data = line_dt,
    aes(xintercept = x, linetype = line_type, color = line_type),
    linewidth = 0.7,
    show.legend = TRUE
  ) +
  scale_color_manual(values = c(type_cols, line_cols), breaks = c(names(type_cols), names(line_cols))) +
  scale_linetype_manual(values = line_types) +
  coord_cartesian(xlim = c(xmin, xmax)) +
  theme_classic(base_size = 13) +
  labs(
    x = "IBS",
    y = "Cumulative proportion",
    color = "Distribution / Threshold",
    linetype = "Threshold",
    title = paste0(prefix, " 2group IBS cumulative distributions"),
    subtitle = paste0("Background-aware threshold (between-group q99) = ", sprintf("%.4f", bg_q99))
  )

ggsave(file.path(opt[["outdir"]], paste0(prefix, "_ibs_ecdf_distributions.pdf")), p2, width = 9.5, height = 5.8)
ggsave(file.path(opt[["outdir"]], paste0(prefix, "_ibs_ecdf_distributions.png")), p2, width = 9.5, height = 5.8, dpi = 300)

p3 <- ggplot(plot_dt, aes(x = type, y = ibs, fill = type)) +
  geom_boxplot(outlier.size = 0.5, width = 0.65, na.rm = TRUE) +
  geom_hline(yintercept = 0.90, linetype = "dashed", linewidth = 0.6, color = "black") +
  geom_hline(yintercept = bg_q99, linetype = "dotdash", linewidth = 0.7, color = "#7B3294") +
  geom_hline(yintercept = 0.99, linetype = "dashed", linewidth = 0.6, color = "#B2182B") +
  scale_fill_manual(values = type_cols) +
  coord_cartesian(ylim = c(xmin, xmax)) +
  theme_classic(base_size = 13) +
  labs(
    x = "",
    y = "IBS",
    fill = "Distribution",
    title = paste0(prefix, " 2group IBS boxplot comparison"),
    subtitle = paste0("Background-aware threshold = ", sprintf("%.4f", bg_q99))
  ) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

ggsave(file.path(opt[["outdir"]], paste0(prefix, "_ibs_boxplot_comparison.pdf")), p3, width = 8.8, height = 5.8)
ggsave(file.path(opt[["outdir"]], paste0(prefix, "_ibs_boxplot_comparison.png")), p3, width = 8.8, height = 5.8, dpi = 300)

self_scatter_dt <- expected_dt[
  !is.na(anchor_self_ibs) & !is.na(secondary_self_ibs),
  .(pair_index, expected_anchor, expected_secondary, anchor_self_ibs, secondary_self_ibs)
]
fwrite(self_scatter_dt, file.path(opt[["outdir"]], paste0(prefix, "_self_identity_scatter_input.tsv")), sep = "\t")

if (nrow(self_scatter_dt) > 0) {
  p4 <- ggplot(self_scatter_dt, aes(x = anchor_self_ibs, y = secondary_self_ibs)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.8, color = "#7F7F7F") +
    geom_point(size = 2.6, alpha = 0.85, color = "#4E79A7") +
    coord_cartesian(xlim = c(xmin, xmax), ylim = c(xmin, xmax)) +
    theme_classic(base_size = 13) +
    labs(
      x = paste0(anchor_col, " self IBS"),
      y = paste0(secondary_col, " self IBS"),
      title = paste0(prefix, " self-identity y=x scatter"),
      subtitle = "Each point is one expected anchor-secondary pair; self-vs-self values are shown against the y=x reference"
    )

  ggsave(file.path(opt[["outdir"]], paste0(prefix, "_self_identity_y_eq_x_scatter.pdf")), p4, width = 7.2, height = 6.2)
  ggsave(file.path(opt[["outdir"]], paste0(prefix, "_self_identity_y_eq_x_scatter.png")), p4, width = 7.2, height = 6.2, dpi = 300)
}

pair_scatter_dt <- copy(expected_dt)
pair_scatter_dt[, match_status := fifelse(
  is.na(ibs), "No_data",
  fifelse(ibs >= strict_thr, ">=0.99", fifelse(ibs >= loose_thr, paste0(">=", sprintf("%.3f", loose_thr)), "<loose_threshold"))
)]
pair_scatter_dt[, pair_label := paste0(expected_anchor, " vs ", expected_secondary)]
fwrite(pair_scatter_dt, file.path(opt[["outdir"]], paste0(prefix, "_one_to_one_pair_scatter_input.tsv")), sep = "\t")

if (nrow(pair_scatter_dt) > 0) {
  pair_cols <- c("#B2182B", "#F28E2B", "#4E79A7", "#7F7F7F")
  names(pair_cols) <- c(">=0.99", paste0(">=", sprintf("%.3f", loose_thr)), "<loose_threshold", "No_data")

  p5 <- ggplot(pair_scatter_dt, aes(x = pair_index, y = ibs, color = match_status)) +
    geom_point(size = 2.8, alpha = 0.9, na.rm = TRUE) +
    geom_hline(yintercept = 0.90, linetype = "dashed", linewidth = 0.6, color = "black") +
    geom_hline(yintercept = bg_q99, linetype = "dotdash", linewidth = 0.7, color = "#7B3294") +
    geom_hline(yintercept = 0.99, linetype = "dashed", linewidth = 0.6, color = "#B2182B") +
    scale_color_manual(values = pair_cols, drop = FALSE) +
    coord_cartesian(ylim = c(xmin, xmax)) +
    theme_classic(base_size = 13) +
    labs(
      x = "Expected pair order",
      y = "Expected pair IBS",
      color = "Pair class",
      title = paste0(prefix, " one-to-one expected-pair scatter"),
      subtitle = paste0("Map-defined ", anchor_col, "-", secondary_col, " one-to-one pairs")
    )

  ggsave(file.path(opt[["outdir"]], paste0(prefix, "_one_to_one_pair_scatter.pdf")), p5, width = 9.5, height = 5.8)
  ggsave(file.path(opt[["outdir"]], paste0(prefix, "_one_to_one_pair_scatter.png")), p5, width = 9.5, height = 5.8, dpi = 300)
}

msg("Done.")
msg("Background-aware threshold (between-group q99): %.6f", bg_q99)
msg("Loose threshold: %.6f", loose_thr)
msg("Strict threshold: %.6f", strict_thr)
msg("Output dir: %s", opt[["outdir"]])
