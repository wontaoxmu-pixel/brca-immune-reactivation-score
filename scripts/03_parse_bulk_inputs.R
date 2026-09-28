#!/usr/bin/env Rscript

set.seed(123)

project_root <- normalizePath(getwd(), mustWork = TRUE)
if (basename(project_root) == "scripts") {
  project_root <- normalizePath(file.path(project_root, ".."), mustWork = TRUE)
}

project_lib <- file.path(project_root, ".Rlib")
if (dir.exists(project_lib)) {
  .libPaths(c(project_lib, .libPaths()))
}

suppressPackageStartupMessages({
  library(data.table)
})

raw_dir <- file.path(project_root, "data", "raw")
processed_dir <- file.path(project_root, "data", "processed")
table_dir <- file.path(project_root, "results", "tables")
log_dir <- file.path(project_root, "logs")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

clean_token <- function(x) {
  x <- gsub('^"|"$', "", x)
  x <- trimws(x)
  x[x %in% c("", "NA", "na", "N/A", "n/a")] <- NA_character_
  x
}

clean_name <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
  ifelse(nchar(x) == 0, "field", x)
}

read_gz_lines <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  readLines(con, warn = FALSE)
}

split_matrix_line <- function(line) {
  clean_token(strsplit(line, "\t", fixed = TRUE)[[1]])
}

parse_series_sample_metadata <- function(path, dataset, platform_hint = NA_character_) {
  lines <- read_gz_lines(path)
  sample_lines <- grep("^!Sample_", lines, value = TRUE)
  titles <- sample_lines[grep("^!Sample_title\t", sample_lines)]
  if (length(titles) == 0) {
    stop("No !Sample_title row found in ", basename(path))
  }
  sample_title <- split_matrix_line(titles[1])[-1]
  n_samples <- length(sample_title)

  meta <- data.table(
    dataset = dataset,
    source_file = basename(path),
    platform = platform_hint,
    sample_title = sample_title
  )

  characteristic_seen <- list()

  for (line in sample_lines) {
    parts <- split_matrix_line(line)
    field <- sub("^!", "", parts[1])
    values <- parts[-1]
    if (length(values) != n_samples) {
      next
    }

    if (field == "Sample_title") {
      next
    }

    if (field == "Sample_characteristics_ch1") {
      first_value <- values[which(!is.na(values))[1]]
      if (!is.na(first_value) && grepl(":", first_value, fixed = TRUE)) {
        key <- sub(":.*$", "", first_value)
        col <- clean_name(key)
        parsed <- sub("^[^:]+:\\s*", "", values)
      } else {
        col <- paste0("characteristics_ch1_", length(characteristic_seen) + 1)
        parsed <- values
      }
      characteristic_seen[[length(characteristic_seen) + 1]] <- col
    } else {
      col <- clean_name(sub("^Sample_", "", field))
      parsed <- values
    }

    if (col %in% names(meta)) {
      col <- make.unique(c(names(meta), col), sep = "_")[length(names(meta)) + 1]
    }
    meta[, (col) := clean_token(parsed)]
  }

  meta[]
}

extract_series_expression <- function(path) {
  lines <- read_gz_lines(path)
  begin <- grep("^!series_matrix_table_begin", lines)
  end <- grep("^!series_matrix_table_end", lines)
  if (length(begin) != 1 || length(end) != 1 || end <= begin + 1) {
    return(NULL)
  }

  table_lines <- lines[(begin + 1):(end - 1)]
  if (length(table_lines) < 2) {
    return(NULL)
  }

  expr <- fread(text = paste(table_lines, collapse = "\n"))
  setnames(expr, 1, "feature_id")
  expr[]
}

write_table <- function(x, path) {
  fwrite(x, path, sep = "\t", quote = FALSE, na = "NA")
}

series_inputs <- data.table(
  dataset = c("GSE58812", "GSE96058", "GSE96058"),
  platform = c("GPL570", "GPL11154", "GPL18573"),
  path = file.path(raw_dir, c(
    "GSE58812_series_matrix.txt.gz",
    "GSE96058-GPL11154_series_matrix.txt.gz",
    "GSE96058-GPL18573_series_matrix.txt.gz"
  ))
)

sample_metadata_list <- list()
expression_outputs <- list()
summary_rows <- list()

for (i in seq_len(nrow(series_inputs))) {
  input <- series_inputs[i]
  if (!file.exists(input$path)) {
    summary_rows[[length(summary_rows) + 1]] <- data.table(
      dataset = input$dataset,
      source_file = basename(input$path),
      data_type = "series_matrix",
      samples = NA_integer_,
      features = NA_integer_,
      output_file = NA_character_,
      status = "missing raw file"
    )
    next
  }

  meta <- parse_series_sample_metadata(input$path, input$dataset, input$platform)
  sample_metadata_list[[length(sample_metadata_list) + 1]] <- meta

  expr <- extract_series_expression(input$path)
  expr_out <- NA_character_
  feature_count <- NA_integer_
  status <- "metadata parsed"
  if (!is.null(expr) && nrow(expr) > 0 && ncol(expr) > 1) {
    expr_out <- file.path(
      processed_dir,
      paste0(input$dataset, "_", input$platform, "_series_expression.tsv.gz")
    )
    write_table(expr, expr_out)
    expression_outputs[[length(expression_outputs) + 1]] <- expr_out
    feature_count <- nrow(expr)
    status <- "metadata and expression parsed"
  }

  summary_rows[[length(summary_rows) + 1]] <- data.table(
    dataset = input$dataset,
    source_file = basename(input$path),
    data_type = "series_matrix",
    samples = nrow(meta),
    features = feature_count,
    output_file = ifelse(is.na(expr_out), NA_character_, basename(expr_out)),
    status = status
  )
}

all_series_metadata <- rbindlist(sample_metadata_list, fill = TRUE)
write_table(all_series_metadata, file.path(processed_dir, "series_sample_metadata.tsv.gz"))

gse176078_path <- file.path(raw_dir, "GSE176078_bulkRNAseq_raw_counts.txt.gz")
if (file.exists(gse176078_path)) {
  bulk_counts <- fread(
    cmd = paste("tar -xOzf", shQuote(gse176078_path)),
    check.names = FALSE
  )
  setnames(bulk_counts, 1, "gene_symbol")
  counts_out <- file.path(processed_dir, "GSE176078_bulkRNAseq_raw_counts.tsv.gz")
  write_table(bulk_counts, counts_out)

  bulk_sample_metadata <- data.table(
    dataset = "GSE176078",
    sample_id = setdiff(names(bulk_counts), "gene_symbol"),
    source_file = basename(gse176078_path),
    assay = "matched bulk RNA-seq raw counts"
  )
  write_table(bulk_sample_metadata, file.path(processed_dir, "GSE176078_bulk_sample_metadata.tsv.gz"))

  summary_rows[[length(summary_rows) + 1]] <- data.table(
    dataset = "GSE176078",
    source_file = basename(gse176078_path),
    data_type = "bulk_raw_counts",
    samples = ncol(bulk_counts) - 1,
    features = nrow(bulk_counts),
    output_file = basename(counts_out),
    status = "expression parsed"
  )
}

parse_summary <- rbindlist(summary_rows, fill = TRUE)
write.csv(parse_summary, file.path(table_dir, "parsed_bulk_data_summary.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(log_dir, "sessionInfo_03_parse_bulk_inputs.txt"))

cat("Bulk input parsing complete\n")
cat("Series metadata rows:", nrow(all_series_metadata), "\n")
cat("Expression outputs:", length(expression_outputs) + as.integer(file.exists(file.path(processed_dir, "GSE176078_bulkRNAseq_raw_counts.tsv.gz"))), "\n")
print(parse_summary)
