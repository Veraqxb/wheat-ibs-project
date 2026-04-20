#!/usr/bin/env Rscript

`%||%` <- function(x, y) {
  if (is.null(x) || identical(x, "")) y else x
}

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

script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)[1]
  if (is.na(file_arg)) return(getwd())
  dirname(normalizePath(sub("^--file=", "", file_arg)))
}

standardize_id <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "-", "NULL")] <- NA_character_
  x
}

extract_sample_number <- function(x) {
  x <- standardize_id(x)
  out <- rep(NA_character_, length(x))
  ok <- !is.na(x)
  m <- regexpr("([0-9]+)$", x[ok], perl = TRUE)
  out[ok] <- ifelse(m > 0, regmatches(x[ok], m), NA_character_)
  out
}

read_sample_map <- function(path) {
  map_df <- read.table(
    path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    fill = TRUE,
    na.strings = c("", "NA", " ", "-", "NULL")
  )
  if (ncol(map_df) == 0) stop("Sample map has no columns: ", path)
  for (col in names(map_df)) {
    map_df[[col]] <- standardize_id(map_df[[col]])
  }
  map_df
}

read_mibs_matrix <- function(mibs_path, id_path) {
  ids_raw <- read.table(id_path, colClasses = "character", stringsAsFactors = FALSE)
  if (ncol(ids_raw) < 2) {
    stop("IBS id file format is invalid: expected at least 2 columns in ", id_path)
  }
  ids <- standardize_id(ids_raw[[2]])
  n <- length(ids)
  if (n == 0) stop("No sample IDs found in ", id_path)

  mibs_vec <- scan(mibs_path, quiet = TRUE)
  expected_tri <- n * (n + 1) / 2
  expected_square <- n * n

  if (length(mibs_vec) == expected_square) {
    mat <- matrix(mibs_vec, nrow = n, ncol = n, byrow = TRUE)
  } else if (length(mibs_vec) == expected_tri) {
    mat <- matrix(NA_real_, n, n)
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
      "IBS matrix size mismatch: expected either ", expected_square,
      " (square) or ", expected_tri, " (lower triangle) values, got ",
      length(mibs_vec), " in ", mibs_path
    )
  }

  rownames(mat) <- ids
  colnames(mat) <- ids
  mat
}

write_matrix_tsv <- function(mat, path) {
  write.table(
    cbind(sample = rownames(mat), mat),
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

get_upper_triangle_values <- function(m) {
  if (is.null(m) || nrow(m) < 2 || ncol(m) < 2) return(numeric(0))
  vals <- m[lower.tri(m, diag = FALSE)]
  vals[!is.na(vals)]
}

build_similarity_clusters <- function(sample_ids, ibs_mat, threshold = 0.99) {
  ids <- unique(standardize_id(sample_ids))
  ids <- ids[!is.na(ids) & ids %in% rownames(ibs_mat)]
  if (length(ids) == 0) {
    return(data.frame(
      sample_id = character(0),
      cluster_id = character(0),
      cluster_size = integer(0),
      cluster_members = character(0),
      neighbor_ids = character(0),
      max_neighbor_ibs = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  sub_mat <- ibs_mat[ids, ids, drop = FALSE]
  adj <- setNames(vector("list", length(ids)), ids)
  hi_idx <- which(upper.tri(sub_mat) & !is.na(sub_mat) & sub_mat > threshold, arr.ind = TRUE)
  if (nrow(hi_idx) > 0) {
    for (k in seq_len(nrow(hi_idx))) {
      a <- rownames(sub_mat)[hi_idx[k, 1]]
      b <- colnames(sub_mat)[hi_idx[k, 2]]
      adj[[a]] <- unique(c(adj[[a]], b))
      adj[[b]] <- unique(c(adj[[b]], a))
    }
  }

  visited <- setNames(rep(FALSE, length(ids)), ids)
  clusters <- list()
  cluster_idx <- 0L
  for (id in ids) {
    if (visited[[id]]) next
    cluster_idx <- cluster_idx + 1L
    stack <- id
    members <- character(0)
    while (length(stack) > 0) {
      cur <- stack[[1]]
      stack <- stack[-1]
      if (visited[[cur]]) next
      visited[[cur]] <- TRUE
      members <- c(members, cur)
      nbrs <- adj[[cur]]
      if (length(nbrs) > 0) {
        stack <- c(stack, nbrs[!visited[nbrs]])
      }
    }
    clusters[[sprintf("cluster_%03d", cluster_idx)]] <- sort(unique(members))
  }

  do.call(
    rbind,
    lapply(names(clusters), function(cid) {
      members <- clusters[[cid]]
      do.call(
        rbind,
        lapply(members, function(sample_id) {
          neighbors <- setdiff(members, sample_id)
          max_ibs <- if (length(neighbors) > 0) suppressWarnings(max(sub_mat[sample_id, neighbors], na.rm = TRUE)) else NA_real_
          if (!is.finite(max_ibs)) max_ibs <- NA_real_
          data.frame(
            sample_id = sample_id,
            cluster_id = if (length(members) > 1) cid else NA_character_,
            cluster_size = length(members),
            cluster_members = if (length(members) > 1) paste(members, collapse = ";") else sample_id,
            neighbor_ids = if (length(neighbors) > 0) paste(neighbors, collapse = ";") else "",
            max_neighbor_ibs = max_ibs,
            stringsAsFactors = FALSE
          )
        })
      )
    })
  )
}

draw_heatmap <- function(sub_mat, file, title, zlim = c(0.7, 1.0), show_values = TRUE) {
  if (is.null(sub_mat) || nrow(sub_mat) == 0 || ncol(sub_mat) == 0) return(invisible(NULL))

  cols <- colorRampPalette(rev(c(
    "#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090",
    "#FFFFBF", "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4", "#313695"
  )))(100)

  nx <- ncol(sub_mat)
  ny <- nrow(sub_mat)
  pdf(file, width = max(8, nx * 0.45 + 3), height = max(8, ny * 0.45 + 3))
  par(mar = c(12, 12, 4, 2))
  image_mat <- t(sub_mat)[, ny:1, drop = FALSE]
  image(seq_len(nx), seq_len(ny), image_mat, col = cols, axes = FALSE, zlim = zlim, main = title, xlab = "", ylab = "")
  axis(1, at = seq_len(nx), labels = colnames(sub_mat), las = 2, cex.axis = 0.7)
  axis(2, at = seq_len(ny), labels = rev(rownames(sub_mat)), las = 1, cex.axis = 0.7)
  abline(h = seq(0.5, ny + 0.5, by = 1), col = "grey80")
  abline(v = seq(0.5, nx + 0.5, by = 1), col = "grey80")
  box(col = "grey50")
  if (show_values) {
    for (i in seq_len(nx)) {
      for (j in seq_len(ny)) {
        val <- sub_mat[j, i]
        if (is.na(val)) next
        y <- ny - j + 1
        col_txt <- ifelse(val > 0.94 | val < 0.78, "white", "black")
        text(i, y, sprintf("%.2f", val), cex = ifelse(val >= 0.99, 0.72, 0.58), col = col_txt, font = ifelse(val >= 0.99, 2, 1))
      }
    }
  }
  dev.off()
}

draw_density_plot <- function(internal_vals, pair_vals, out_file, title, internal_label, pair_label, threshold) {
  internal_vals <- internal_vals[!is.na(internal_vals)]
  pair_vals <- pair_vals[!is.na(pair_vals)]
  if (length(internal_vals) + length(pair_vals) == 0) return(invisible(NULL))

  df <- data.frame(
    IBS = c(internal_vals, pair_vals),
    Category = c(rep(internal_label, length(internal_vals)), rep(pair_label, length(pair_vals))),
    stringsAsFactors = FALSE
  )
  df$Category <- factor(df$Category, levels = c(internal_label, pair_label))

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    library(ggplot2)
    cols <- c("#4575b4", "#d73027")
    names(cols) <- levels(df$Category)
    pdf(out_file, width = 8.5, height = 5.5)
    print(
      ggplot(df, aes(x = IBS, color = Category, fill = Category)) +
        geom_density(alpha = 0.18, linewidth = 1.15, adjust = 1.1) +
        geom_vline(xintercept = threshold, linetype = "dashed", color = "grey35", linewidth = 0.7) +
        scale_color_manual(values = cols) +
        scale_fill_manual(values = cols) +
        labs(title = title, x = "IBS", y = "Density") +
        theme_bw(base_size = 12) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5), panel.grid.minor = element_blank())
    )
    dev.off()
  } else {
    pdf(out_file, width = 8.5, height = 5.5)
    plot(NA, xlim = range(df$IBS, na.rm = TRUE), ylim = c(0, 1), xlab = "IBS", ylab = "Density", main = title)
    if (length(internal_vals) >= 2 && diff(range(internal_vals)) > 0) {
      lines(density(internal_vals), col = "#4575b4", lwd = 2)
    }
    if (length(pair_vals) >= 2 && diff(range(pair_vals)) > 0) {
      lines(density(pair_vals), col = "#d73027", lwd = 2)
    }
    abline(v = threshold, lty = 2, col = "grey40")
    legend("topright", legend = c(internal_label, pair_label), col = c("#4575b4", "#d73027"), lwd = 2, bty = "n")
    dev.off()
  }
}

draw_category_barplot <- function(summary_df, out_file, title, group_col = "group_name", category_col = "match_type", count_col = "count") {
  if (!requireNamespace("ggplot2", quietly = TRUE) || nrow(summary_df) == 0) return(invisible(NULL))
  library(ggplot2)
  pdf(out_file, width = 9, height = 5.5)
  print(
    ggplot(summary_df, aes_string(x = group_col, y = count_col, fill = category_col)) +
      geom_col(position = "stack") +
      theme_bw(base_size = 12) +
      labs(title = title, x = "Group", y = "Sample count") +
      theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.text.x = element_text(angle = 30, hjust = 1))
  )
  dev.off()
}

pick_threshold <- function(group_name, rna_groups = c("TC", "SC", "FC"), dna_threshold = 0.99, rna_threshold = 0.90) {
  if (group_name %in% rna_groups) rna_threshold else dna_threshold
}

collect_neighbor_hits <- function(neighbor_ids, candidates) {
  if (is.na(neighbor_ids) || !nzchar(neighbor_ids)) return("")
  neighbors <- unlist(strsplit(neighbor_ids, ";", fixed = TRUE))
  neighbors <- neighbors[nzchar(neighbors)]
  candidates <- unique(candidates[!is.na(candidates) & candidates != ""])
  hits <- intersect(neighbors, candidates)
  if (length(hits) == 0) "" else paste(hits, collapse = ";")
}
