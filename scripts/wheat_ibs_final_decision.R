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

`%||%` <- function(x, y) {
  if (is.null(x) || identical(x, "")) y else x
}

opt <- parse_args(args)
required <- c("map", "dna-pair-summary", "outdir", "prefix", "ploidy")
missing <- required[!required %in% names(opt)]
if (length(missing) > 0) {
  stop("Missing arguments: ", paste(missing, collapse = ", "))
}

opt[["z23-col"]] <- opt[["z23-col"]] %||% "Z23"
opt[["b25-col"]] <- opt[["b25-col"]] %||% "B25"
opt[["tc-col"]] <- opt[["tc-col"]] %||% "TC"
opt[["sc-col"]] <- opt[["sc-col"]] %||% "SC"
opt[["fc-col"]] <- opt[["fc-col"]] %||% "FC"
opt[["dna-ibs-matrix"]] <- opt[["dna-ibs-matrix"]] %||% ""

dir.create(opt[["outdir"]], recursive = TRUE, showWarnings = FALSE)

read_pair_summary <- function(path, sample_col, expected_col) {
  if (!nzchar(path) || !file.exists(path)) return(NULL)
  df <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  required_cols <- c("sample_y", "expected_x", "ibs_expected", "best_x", "ibs_best", "second_best", "margin", "status")
  miss <- required_cols[!required_cols %in% colnames(df)]
  if (length(miss) > 0) {
    stop("Pair summary missing columns in ", path, ": ", paste(miss, collapse = ", "))
  }
  names(df)[names(df) == "sample_y"] <- sample_col
  names(df)[names(df) == "expected_x"] <- expected_col
  df
}

extract_sample_id <- function(...) {
  vals <- list(...)
  for (v in vals) {
    if (length(v) == 0 || is.na(v) || !nzchar(v)) next
    m <- sub(".*_(\\d+)$", "\\1", v)
    if (!identical(m, v)) return(m)
  }
  NA_character_
}

coalesce_chr <- function(...) {
  vals <- list(...)
  out <- vals[[1]]
  for (i in 2:length(vals)) {
    fill <- is.na(out) | out == ""
    out[fill] <- vals[[i]][fill]
  }
  out
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

needed_map_cols <- c(opt[["z23-col"]], opt[["b25-col"]])
missing_map <- needed_map_cols[!needed_map_cols %in% colnames(map_df)]
if (length(missing_map) > 0) {
  stop("Map file missing required columns: ", paste(missing_map, collapse = ", "))
}

for (col in c(opt[["tc-col"]], opt[["sc-col"]], opt[["fc-col"]])) {
  if (!(col %in% colnames(map_df))) {
    map_df[[col]] <- NA_character_
  }
}

base_df <- data.frame(
  sample_id = vapply(
    seq_len(nrow(map_df)),
    function(i) extract_sample_id(
      map_df[[opt[["z23-col"]]]][i],
      map_df[[opt[["b25-col"]]]][i],
      map_df[[opt[["tc-col"]]]][i],
      map_df[[opt[["sc-col"]]]][i],
      map_df[[opt[["fc-col"]]]][i]
    ),
    character(1)
  ),
  ploidy = opt[["ploidy"]],
  dataset = opt[["prefix"]],
  expected_z23 = map_df[[opt[["z23-col"]]]],
  expected_b25 = map_df[[opt[["b25-col"]]]],
  expected_tc = map_df[[opt[["tc-col"]]]],
  expected_sc = map_df[[opt[["sc-col"]]]],
  expected_fc = map_df[[opt[["fc-col"]]]],
  stringsAsFactors = FALSE
)

dna_df <- read_pair_summary(opt[["dna-pair-summary"]], "expected_b25", "expected_z23")
base_df <- merge(base_df, dna_df, by = c("expected_b25", "expected_z23"), all.x = TRUE, sort = FALSE)
names(base_df)[names(base_df) == "best_x"] <- "dna_best_match"
names(base_df)[names(base_df) == "ibs_expected"] <- "dna_ibs_expected"
names(base_df)[names(base_df) == "ibs_best"] <- "dna_ibs_best"
names(base_df)[names(base_df) == "second_best"] <- "dna_second_best"
names(base_df)[names(base_df) == "margin"] <- "dna_margin"
names(base_df)[names(base_df) == "status"] <- "dna_status"

join_rna_summary <- function(df, opt_name, source_name, expected_col) {
  path <- opt[[opt_name]] %||% ""
  if (!nzchar(path)) return(df)
  rna_df <- read_pair_summary(path, expected_col, "rna_expected_target")
  keep_cols <- c(expected_col, "ibs_expected", "best_x", "ibs_best", "second_best", "margin", "status")
  rna_df <- rna_df[, keep_cols, drop = FALSE]
  suffix <- tolower(source_name)
  names(rna_df) <- c(
    expected_col,
    paste0(suffix, "_ibs_expected"),
    paste0(suffix, "_best_match"),
    paste0(suffix, "_ibs_best"),
    paste0(suffix, "_second_best"),
    paste0(suffix, "_margin"),
    paste0(suffix, "_status")
  )
  merge(df, rna_df, by = expected_col, all.x = TRUE, sort = FALSE)
}

for (spec in list(
  list(opt = "tc-z23-summary", name = "TC_Z23", col = "expected_tc"),
  list(opt = "sc-z23-summary", name = "SC_Z23", col = "expected_sc"),
  list(opt = "fc-z23-summary", name = "FC_Z23", col = "expected_fc"),
  list(opt = "tc-b25-summary", name = "TC_B25", col = "expected_tc"),
  list(opt = "sc-b25-summary", name = "SC_B25", col = "expected_sc"),
  list(opt = "fc-b25-summary", name = "FC_B25", col = "expected_fc")
)) {
  base_df <- join_rna_summary(base_df, spec$opt, spec$name, spec$col)
}

resolve_rna_source <- function(df, src) {
  z_col <- tolower(paste0(src, "_z23_status"))
  b_col <- tolower(paste0(src, "_b25_status"))
  expected_col <- paste0("expected_", tolower(src))
  out <- data.frame(
    source = src,
    expected_sample = df[[expected_col]],
    z23_status = if (z_col %in% names(df)) df[[z_col]] else NA_character_,
    b25_status = if (b_col %in% names(df)) df[[b_col]] else NA_character_,
    z23_best = if (tolower(paste0(src, "_z23_best_match")) %in% names(df)) df[[tolower(paste0(src, "_z23_best_match"))]] else NA_character_,
    b25_best = if (tolower(paste0(src, "_b25_best_match")) %in% names(df)) df[[tolower(paste0(src, "_b25_best_match"))]] else NA_character_,
    z23_ibs = if (tolower(paste0(src, "_z23_ibs_expected")) %in% names(df)) df[[tolower(paste0(src, "_z23_ibs_expected"))]] else NA_real_,
    b25_ibs = if (tolower(paste0(src, "_b25_ibs_expected")) %in% names(df)) df[[tolower(paste0(src, "_b25_ibs_expected"))]] else NA_real_,
    z23_margin = if (tolower(paste0(src, "_z23_margin")) %in% names(df)) df[[tolower(paste0(src, "_z23_margin"))]] else NA_real_,
    b25_margin = if (tolower(paste0(src, "_b25_margin")) %in% names(df)) df[[tolower(paste0(src, "_b25_margin"))]] else NA_real_,
    stringsAsFactors = FALSE
  )

  out$rna_status <- ifelse(
    out$z23_status == "MATCH", "MATCH_Z23",
    ifelse(out$b25_status == "MATCH", "MATCH_B25",
      ifelse((is.na(out$z23_status) | out$z23_status == "NO_DATA") & (is.na(out$b25_status) | out$b25_status == "NO_DATA"), "NO_DATA", "MISMATCH")
    )
  )
  out$rna_best_match <- ifelse(out$rna_status == "MATCH_Z23", out$z23_best, ifelse(out$rna_status == "MATCH_B25", out$b25_best, coalesce_chr(out$z23_best, out$b25_best)))
  out$rna_ibs_expected <- ifelse(out$rna_status == "MATCH_B25", out$b25_ibs, out$z23_ibs)
  out$rna_margin <- ifelse(out$rna_status == "MATCH_B25", out$b25_margin, out$z23_margin)
  out
}

rna_sources <- c("tc", "sc", "fc")
rna_details <- lapply(rna_sources, function(src) resolve_rna_source(base_df, src))
names(rna_details) <- toupper(rna_sources)

base_df$rna_source <- apply(
  do.call(cbind, lapply(rna_details, function(x) x$rna_status)),
  1,
  function(x) {
    good <- names(rna_details)[x %in% c("MATCH_Z23", "MATCH_B25")]
    if (length(good) == 0) NA_character_ else paste(good, collapse = ",")
  }
)

base_df$rna_status <- apply(
  do.call(cbind, lapply(rna_details, function(x) x$rna_status)),
  1,
  function(x) {
    if (any(x == "MATCH_Z23", na.rm = TRUE)) {
      "MATCH_Z23"
    } else if (any(x == "MATCH_B25", na.rm = TRUE)) {
      "MATCH_B25"
    } else if (all(is.na(x) | x == "NO_DATA")) {
      "NO_DATA"
    } else {
      "MISMATCH"
    }
  }
)

base_df$rna_best_match <- apply(
  do.call(cbind, lapply(rna_details, function(x) x$rna_best_match)),
  1,
  function(x) {
    x <- x[!is.na(x) & x != ""]
    if (length(x) == 0) NA_character_ else x[1]
  }
)

base_df$rna_ibs_expected <- apply(
  do.call(cbind, lapply(rna_details, function(x) x$rna_ibs_expected)),
  1,
  function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[!is.na(x)]
    if (length(x) == 0) NA_real_ else max(x)
  }
)

base_df$rna_margin <- apply(
  do.call(cbind, lapply(rna_details, function(x) x$rna_margin)),
  1,
  function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[!is.na(x)]
    if (length(x) == 0) NA_real_ else max(x)
  }
)

base_df$final_decision <- ifelse(
  is.na(base_df$dna_status) | base_df$dna_status == "NO_DATA",
  ifelse(base_df$rna_status == "MATCH_B25", "RESEQ", ifelse(base_df$rna_status == "MATCH_Z23", "REVIEW", "NO_DATA")),
  ifelse(
    base_df$dna_status == "MATCH",
    ifelse(base_df$rna_status == "MATCH_B25", "REVIEW", "KEEP"),
    ifelse(base_df$rna_status == "MATCH_B25", "RESEQ", ifelse(base_df$rna_status == "MATCH_Z23", "REVIEW", "REMOVE"))
  )
)

base_df$action <- ifelse(
  base_df$final_decision == "KEEP", "keep",
  ifelse(base_df$final_decision == "RESEQ", "keep_for_reseq",
    ifelse(base_df$final_decision == "REVIEW", "manual_review",
      ifelse(base_df$final_decision == "REMOVE", "remove", "manual_review")
    )
  )
)

base_df$comment <- ifelse(
  base_df$final_decision == "KEEP", "dna_match_primary",
  ifelse(base_df$final_decision == "RESEQ", "dna_issue_rna_supports_b25",
    ifelse(base_df$final_decision == "REVIEW", "check_margin_or_conflicting_support",
      ifelse(base_df$final_decision == "REMOVE", "dna_and_rna_do_not_support_retention", "insufficient_data")
    )
  )
)

final_cols <- c(
  "sample_id", "ploidy", "dataset", "expected_z23", "expected_b25", "expected_tc", "expected_sc", "expected_fc",
  "dna_best_match", "dna_ibs_expected", "dna_ibs_best", "dna_second_best", "dna_margin", "dna_status",
  "rna_source", "rna_best_match", "rna_ibs_expected", "rna_margin", "rna_status",
  "final_decision", "action", "comment"
)
final_df <- base_df[, final_cols, drop = FALSE]

outdir <- opt[["outdir"]]
prefix <- opt[["prefix"]]

write.table(final_df, file = file.path(outdir, paste0(prefix, "_final_decision.tsv")), quote = FALSE, sep = "\t", row.names = FALSE)
write.csv(final_df, file = file.path(outdir, paste0(prefix, "_final_decision.csv")), row.names = FALSE, quote = TRUE)

summary_df <- data.frame(
  dataset = prefix,
  ploidy = opt[["ploidy"]],
  total_samples = nrow(final_df),
  keep_count = sum(final_df$final_decision == "KEEP"),
  reseq_count = sum(final_df$final_decision == "RESEQ"),
  review_count = sum(final_df$final_decision == "REVIEW"),
  remove_count = sum(final_df$final_decision == "REMOVE"),
  no_data_count = sum(final_df$final_decision == "NO_DATA"),
  final_retained_count = sum(final_df$final_decision %in% c("KEEP", "RESEQ")),
  stringsAsFactors = FALSE
)
write.table(summary_df, file = file.path(outdir, paste0(prefix, "_final_summary.tsv")), quote = FALSE, sep = "\t", row.names = FALSE)

write.table(final_df[final_df$final_decision == "REMOVE", , drop = FALSE], file = file.path(outdir, paste0(prefix, "_remove_candidates.tsv")), quote = FALSE, sep = "\t", row.names = FALSE)
write.table(final_df[final_df$final_decision == "RESEQ", , drop = FALSE], file = file.path(outdir, paste0(prefix, "_reseq_candidates.tsv")), quote = FALSE, sep = "\t", row.names = FALSE)
write.table(final_df[final_df$final_decision == "REVIEW", , drop = FALSE], file = file.path(outdir, paste0(prefix, "_review_candidates.tsv")), quote = FALSE, sep = "\t", row.names = FALSE)

write_html_table <- function(df, summary_df, file) {
  row_color <- function(x) {
    switch(
      x,
      KEEP = "#d9ead3",
      RESEQ = "#fff2cc",
      REVIEW = "#fce5cd",
      REMOVE = "#f4cccc",
      NO_DATA = "#d9d9d9",
      "#ffffff"
    )
  }

  con <- file(file, "w")
  on.exit(close(con), add = TRUE)
  writeLines("<html><head><meta charset='utf-8'><style>body{font-family:Arial,sans-serif} table{border-collapse:collapse;font-size:12px} th,td{border:1px solid #999;padding:4px 6px} th{background:#f0f0f0;position:sticky;top:0} .summary{margin-bottom:20px}</style></head><body>", con)
  writeLines("<h2>Final Decision Summary</h2>", con)
  writeLines("<table class='summary'>", con)
  writeLines("<tr>" %+% paste(sprintf("<th>%s</th>", names(summary_df)), collapse = "") %+% "</tr>", con)
  writeLines("<tr>" %+% paste(sprintf("<td>%s</td>", summary_df[1, ]), collapse = "") %+% "</tr>", con)
  writeLines("</table>", con)
  writeLines("<h2>Final Decision Table</h2>", con)
  writeLines("<table>", con)
  writeLines("<tr>" %+% paste(sprintf("<th>%s</th>", names(df)), collapse = "") %+% "</tr>", con)
  for (i in seq_len(nrow(df))) {
    color <- row_color(df$final_decision[i])
    writeLines(
      "<tr style='background:" %+% color %+% "'>" %+%
        paste(sprintf("<td>%s</td>", ifelse(is.na(df[i, ]), "", as.character(df[i, ]))), collapse = "") %+%
        "</tr>",
      con
    )
  }
  writeLines("</table></body></html>", con)
}

`%+%` <- function(a, b) paste0(a, b)
write_html_table(final_df, summary_df, file.path(outdir, paste0(prefix, "_final_decision.html")))

if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "final_decision")
  openxlsx::writeData(wb, "final_decision", final_df)
  styles <- list(
    KEEP = openxlsx::createStyle(fgFill = "#d9ead3"),
    RESEQ = openxlsx::createStyle(fgFill = "#fff2cc"),
    REVIEW = openxlsx::createStyle(fgFill = "#fce5cd"),
    REMOVE = openxlsx::createStyle(fgFill = "#f4cccc"),
    NO_DATA = openxlsx::createStyle(fgFill = "#d9d9d9")
  )
  for (i in seq_len(nrow(final_df))) {
    st <- styles[[final_df$final_decision[i]]]
    if (!is.null(st)) {
      openxlsx::addStyle(wb, "final_decision", st, rows = i + 1, cols = seq_len(ncol(final_df)), gridExpand = TRUE, stack = TRUE)
    }
  }
  openxlsx::saveWorkbook(wb, file.path(outdir, paste0(prefix, "_final_decision.xlsx")), overwrite = TRUE)
}

plot_color <- c(KEEP = "#5ab769", RESEQ = "#d8b500", REVIEW = "#e68a2e", REMOVE = "#d64545", NO_DATA = "#7f7f7f")

pdf(file.path(outdir, paste0(prefix, "_dna_expected_best_scatter.pdf")), width = 7, height = 6)
plot(
  final_df$dna_ibs_expected,
  final_df$dna_ibs_best,
  pch = 19,
  col = plot_color[final_df$final_decision],
  xlab = "DNA IBS to expected pair",
  ylab = "DNA IBS to best match",
  main = paste(prefix, "DNA Expected vs Best")
)
abline(0, 1, lty = 2, col = "grey40")
legend("bottomright", legend = names(plot_color), col = plot_color, pch = 19, bty = "n")
dev.off()

box_df <- data.frame(
  value = c(final_df$dna_ibs_expected, final_df$dna_ibs_best, final_df$rna_ibs_expected),
  metric = c(
    rep("DNA_expected", nrow(final_df)),
    rep("DNA_best", nrow(final_df)),
    rep("RNA_expected", nrow(final_df))
  )
)
box_df <- box_df[!is.na(box_df$value), , drop = FALSE]
if (nrow(box_df) > 0) {
  pdf(file.path(outdir, paste0(prefix, "_ibs_boxplot.pdf")), width = 7, height = 5)
  boxplot(value ~ metric, data = box_df, col = c("#abd9e9", "#74add1", "#fdae61"), ylab = "IBS", main = paste(prefix, "IBS Distribution"))
  dev.off()
}

margin_df <- final_df[!is.na(final_df$dna_margin), , drop = FALSE]
if (nrow(margin_df) > 0) {
  margin_df <- margin_df[order(margin_df$dna_margin), , drop = FALSE]
  pdf(file.path(outdir, paste0(prefix, "_dna_margin_rank.pdf")), width = 10, height = 5)
  barplot(
    margin_df$dna_margin,
    names.arg = margin_df$sample_id,
    las = 2,
    cex.names = 0.6,
    col = plot_color[margin_df$final_decision],
    ylab = "DNA margin",
    main = paste(prefix, "DNA Margin Rank")
  )
  abline(h = 0.01, col = "red", lty = 2)
  dev.off()
}

if (nzchar(opt[["dna-ibs-matrix"]]) && file.exists(opt[["dna-ibs-matrix"]])) {
  ibs_mat_df <- read.table(opt[["dna-ibs-matrix"]], header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  if (ncol(ibs_mat_df) > 1) {
    row_ids <- ibs_mat_df[[1]]
    ibs_mat <- as.matrix(ibs_mat_df[, -1, drop = FALSE])
    mode(ibs_mat) <- "numeric"
    rownames(ibs_mat) <- row_ids
    colnames(ibs_mat) <- colnames(ibs_mat_df)[-1]
    common_ids <- intersect(final_df$expected_b25, row_ids)
    if (length(common_ids) >= 3) {
      dist_mat <- 1 - ibs_mat[common_ids, common_ids, drop = FALSE]
      dist_mat[is.na(dist_mat)] <- max(dist_mat, na.rm = TRUE)
      mds <- cmdscale(as.dist(dist_mat), k = 2)
      plot_df <- data.frame(
        sample = rownames(mds),
        x = mds[, 1],
        y = mds[, 2],
        stringsAsFactors = FALSE
      )
      plot_df <- merge(plot_df, final_df[, c("expected_b25", "sample_id", "final_decision")], by.x = "sample", by.y = "expected_b25", all.x = TRUE)
      pdf(file.path(outdir, paste0(prefix, "_dna_mds.pdf")), width = 7, height = 6)
      plot(plot_df$x, plot_df$y, pch = 19, col = plot_color[plot_df$final_decision], xlab = "MDS1", ylab = "MDS2", main = paste(prefix, "DNA MDS (1-IBS)"))
      text(plot_df$x, plot_df$y, labels = ifelse(plot_df$final_decision %in% c("RESEQ", "REVIEW", "REMOVE", "NO_DATA"), plot_df$sample_id, ""), pos = 3, cex = 0.7)
      legend("topright", legend = names(plot_color), col = plot_color, pch = 19, bty = "n")
      dev.off()
    }
  }
}

cat("Final decision outputs written to:", outdir, "\n")
