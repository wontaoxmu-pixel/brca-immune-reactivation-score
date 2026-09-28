#!/usr/bin/env Rscript

set.seed(123)

project_root <- normalizePath(getwd(), mustWork = TRUE)
if (basename(project_root) == "scripts") {
  project_root <- normalizePath(file.path(project_root, ".."), mustWork = TRUE)
}

raw_dir <- file.path(project_root, "data", "raw")
out_dir <- file.path(project_root, "results", "tables")
log_dir <- file.path(project_root, "logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

raw_files <- list.files(raw_dir, full.names = TRUE)
file_info <- data.frame(
  file = basename(raw_files),
  size_bytes = file.info(raw_files)$size,
  stringsAsFactors = FALSE
)

peek_gz <- function(path, n = 20) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  readLines(con, n = n, warn = FALSE)
}

extract_series_fields <- function(path) {
  lines <- peek_gz(path, n = 200)
  clean_field <- function(pattern) {
    values <- grep(pattern, lines, value = TRUE)
    values <- sub(paste0("^", pattern, "\\t"), "", values)
    values <- gsub('^"|"$', "", values)
    values
  }
  fields <- list(
    series_title = clean_field("!Series_title")[1],
    geo_accession = clean_field("!Series_geo_accession")[1],
    pubmed_id = paste(clean_field("!Series_pubmed_id"), collapse = ";"),
    platform = paste(clean_field("!Series_platform_id"), collapse = ";")
  )
  as.data.frame(fields, stringsAsFactors = FALSE)
}

series_files <- raw_files[grepl("series_matrix\\.txt\\.gz$", raw_files)]
series_summary <- do.call(rbind, lapply(series_files, extract_series_fields))

gse176078_tar_listing <- character()
gse176078_path <- file.path(raw_dir, "GSE176078_bulkRNAseq_raw_counts.txt.gz")
if (file.exists(gse176078_path)) {
  gse176078_tar_listing <- system2("tar", c("-tzf", gse176078_path), stdout = TRUE)
}

write.csv(file_info, file.path(out_dir, "raw_file_inventory.csv"), row.names = FALSE)
write.csv(series_summary, file.path(out_dir, "series_matrix_summary.csv"), row.names = FALSE)
writeLines(gse176078_tar_listing, file.path(log_dir, "GSE176078_bulk_tar_listing.txt"))
capture.output(sessionInfo(), file = file.path(log_dir, "sessionInfo_02_raw_data_probe.txt"))

cat("Raw data probe complete\n")
cat("Files inventoried:", nrow(file_info), "\n")
cat("Series matrix files summarized:", nrow(series_summary), "\n")
