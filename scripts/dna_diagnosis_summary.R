#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
if (is.na(script_arg)) stop("Cannot determine script path for dna_diagnosis_summary.R")
script_dirname <- dirname(normalizePath(sub("^--file=", "", script_arg)))
source(file.path(script_dirname, "ibs_common.R"))

opt <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("input-dir", "outdir", "prefix")
missing <- required[!required %in% names(opt)]
if (length(missing) > 0) stop("Missing arguments: ", paste(missing, collapse = ", "))

dir.create(opt[["outdir"]], recursive = TRUE, showWarnings = FALSE)

files <- list.files(opt[["input-dir"]], pattern = "_full_matching_table\\.tsv$|_expected_pair_ibs\\.tsv$|_matching_summary\\.tsv$", recursive = TRUE, full.names = TRUE)
if (length(files) == 0) stop("No diagnosis TSV files found under: ", opt[["input-dir"]])

tables <- lapply(files, function(path) {
  df <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  df$source_file <- path
  df
})

all_cols <- unique(unlist(lapply(tables, names)))
tables <- lapply(tables, function(df) {
  missing_cols <- setdiff(all_cols, names(df))
  for (col in missing_cols) df[[col]] <- NA
  df[, all_cols, drop = FALSE]
})
combined <- do.call(rbind, tables)
write.table(combined, file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_combined_diagnosis_rows.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

label_col <- if ("match_type" %in% names(combined)) {
  "match_type"
} else if ("diagnosis_label" %in% names(combined)) {
  "diagnosis_label"
} else {
  NA_character_
}

if (!is.na(label_col)) {
  summary_df <- as.data.frame(table(label = combined[[label_col]], useNA = "ifany"), stringsAsFactors = FALSE)
  names(summary_df) <- c("label", "count")
  write.table(summary_df, file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_diagnosis_label_summary.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

  taxa_col <- if ("Taxa" %in% names(combined)) {
    "Taxa"
  } else if ("taxa" %in% names(combined)) {
    "taxa"
  } else {
    NA_character_
  }

  if (!is.na(taxa_col)) {
    taxa_df <- combined[, c(taxa_col, label_col), drop = FALSE]
    names(taxa_df) <- c("taxa", "label")
    taxa_df$taxa <- standardize_id(taxa_df$taxa)
    taxa_df$label <- standardize_id(taxa_df$label)
    taxa_df$taxa[is.na(taxa_df$taxa)] <- "Unknown"
    taxa_df$label[is.na(taxa_df$label)] <- "NA"

    taxa_long <- as.data.frame(table(taxa = taxa_df$taxa, label = taxa_df$label), stringsAsFactors = FALSE)
    names(taxa_long) <- c("taxa", "label", "count")
    taxa_long <- taxa_long[taxa_long$count > 0, , drop = FALSE]
    write.table(taxa_long, file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_taxa_diagnosis_count_long.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

    taxa_wide <- reshape(taxa_long, idvar = "taxa", timevar = "label", direction = "wide")
    names(taxa_wide) <- sub("^count\\.", "", names(taxa_wide))
    taxa_wide[is.na(taxa_wide)] <- 0
    taxa_wide$total_n <- rowSums(taxa_wide[, setdiff(names(taxa_wide), "taxa"), drop = FALSE])
    taxa_wide <- taxa_wide[, c("taxa", "total_n", setdiff(names(taxa_wide), c("taxa", "total_n"))), drop = FALSE]
    write.table(taxa_wide, file.path(opt[["outdir"]], paste0(opt[["prefix"]], "_taxa_diagnosis_count_table.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

cat("Summarized", length(files), "diagnosis files\n")
