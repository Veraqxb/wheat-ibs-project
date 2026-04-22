#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
if (is.na(script_arg)) stop("Cannot determine script path for samples.R")
script_dirname <- dirname(normalizePath(sub("^--file=", "", script_arg)))
source(file.path(script_dirname, "ibs_common.R"))

opt <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("config")
missing <- required[!required %in% names(opt)]
if (length(missing) > 0) stop("Missing arguments: ", paste(missing, collapse = ", "))

pipeline <- file.path(script_dirname, "wheat_ibs_modular_pipeline.sh")
if (!file.exists(pipeline)) stop("Cannot find modular pipeline: ", pipeline)

message("Running CAMP DNA/RNA sample identity workflow with config: ", opt[["config"]])
message("Decision rules are defined in the modular step scripts; VCF/PLINK/IBS generation is not changed by this wrapper.")
status <- system2("bash", c(pipeline, opt[["config"]]))
if (!identical(status, 0L)) stop("CAMP DNA/RNA sample workflow failed for config: ", opt[["config"]])

cat("Completed CAMP DNA/RNA sample workflow\n")
