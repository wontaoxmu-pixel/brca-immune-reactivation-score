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

processed_dir <- file.path(project_root, "data", "processed")
raw_dir <- file.path(project_root, "data", "raw")
table_dir <- file.path(project_root, "results", "tables")
log_dir <- file.path(project_root, "logs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

metadata_path <- file.path(processed_dir, "series_sample_metadata.tsv.gz")
if (!file.exists(metadata_path)) {
  stop("Missing parsed series metadata: ", metadata_path)
}

meta <- fread(metadata_path)

to_num <- function(x) suppressWarnings(as.numeric(x))
to_int <- function(x) suppressWarnings(as.integer(x))
known <- function(x) !is.na(x) & x != ""

gse58812 <- copy(meta[dataset == "GSE58812"])
gse58812[, `:=`(
  age_at_diag_num = to_num(age_at_diag),
  meta_event = to_int(meta),
  mfs_days_num = to_num(mfs_days),
  death_event = to_int(death),
  os_days_num = to_num(os_days),
  er_ihc_num = to_int(er_ihc),
  pr_ihc_num = to_int(pr_ihc),
  her2_ihc_num = to_int(her2_ihc)
)]
gse58812[, tnbc_pathology := er_ihc_num == 0 & pr_ihc_num == 0 & her2_ihc_num == 0]

gse96058 <- copy(meta[dataset == "GSE96058"])
gse96058[, `:=`(
  is_technical_replicate = grepl("repl$", sample_title),
  primary_sample_id = sub("repl$", "", sample_title),
  er_status_num = to_int(er_status),
  pgr_status_num = to_int(pgr_status),
  her2_status_num = to_int(her2_status),
  os_days_num = to_num(overall_survival_days),
  os_event = to_int(overall_survival_event),
  tumor_size_num = to_num(tumor_size),
  age_at_diagnosis_num = to_num(age_at_diagnosis)
)]
gse96058[, pathology_complete := !is.na(er_status_num) & !is.na(pgr_status_num) & !is.na(her2_status_num)]
gse96058[, tnbc_pathology := pathology_complete & er_status_num == 0 & pgr_status_num == 0 & her2_status_num == 0]

primary_gse96058 <- gse96058[is_technical_replicate == FALSE]
gse96058_expression_path <- file.path(raw_dir, "GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz")
gse96058_expression_available <- file.exists(gse96058_expression_path)

cohort_summary <- rbindlist(list(
  data.table(
    dataset = "GSE58812",
    platform = "GPL570",
    gsm_rows = nrow(gse58812),
    primary_samples = nrow(gse58812),
    technical_replicates = 0L,
    tnbc_pathology = sum(gse58812$tnbc_pathology, na.rm = TRUE),
    pathology_complete = sum(!is.na(gse58812$er_ihc_num) & !is.na(gse58812$pr_ihc_num) & !is.na(gse58812$her2_ihc_num)),
    os_available = sum(!is.na(gse58812$os_days_num) & !is.na(gse58812$death_event)),
    os_events = sum(gse58812$death_event == 1, na.rm = TRUE),
    mfs_available = sum(!is.na(gse58812$mfs_days_num) & !is.na(gse58812$meta_event)),
    mfs_events = sum(gse58812$meta_event == 1, na.rm = TRUE),
    expression_available = TRUE,
    recommended_use = "TNBC expression validation with OS/MFS"
  ),
  data.table(
    dataset = "GSE96058",
    platform = "GPL11154+GPL18573",
    gsm_rows = nrow(gse96058),
    primary_samples = nrow(primary_gse96058),
    technical_replicates = sum(gse96058$is_technical_replicate),
    tnbc_pathology = sum(primary_gse96058$tnbc_pathology, na.rm = TRUE),
    pathology_complete = sum(primary_gse96058$pathology_complete, na.rm = TRUE),
    os_available = sum(!is.na(primary_gse96058$os_days_num) & !is.na(primary_gse96058$os_event)),
    os_events = sum(primary_gse96058$os_event == 1, na.rm = TRUE),
    mfs_available = NA_integer_,
    mfs_events = NA_integer_,
    expression_available = gse96058_expression_available,
    recommended_use = if (gse96058_expression_available) {
      "Large OS validation; remove technical replicates"
    } else {
      "Large OS validation after expression file download; remove technical replicates"
    }
  )
), fill = TRUE)

pam50_summary <- primary_gse96058[, .N, by = .(pam50_subtype)][order(-N)]
setnames(pam50_summary, "N", "primary_sample_count")

gse96058_tnbc_candidates <- primary_gse96058[
  tnbc_pathology == TRUE,
  .(
    dataset,
    platform,
    sample_title,
    geo_accession,
    scan_b_external_id,
    er_status,
    pgr_status,
    her2_status,
    pam50_subtype,
    overall_survival_days,
    overall_survival_event,
    tumor_size,
    lymph_node_status,
    lymph_node_group,
    endocrine_treated,
    chemo_treated
  )
]

gse58812_survival <- gse58812[
  ,
  .(
    dataset,
    platform,
    sample_title,
    geo_accession,
    diagnosis,
    age_at_diag,
    er_ihc,
    pr_ihc,
    her2_ihc,
    meta,
    mfs_days,
    death,
    os_days
  )
]

write.csv(cohort_summary, file.path(table_dir, "bulk_cohort_metadata_summary.csv"), row.names = FALSE)
write.csv(pam50_summary, file.path(table_dir, "GSE96058_pam50_primary_sample_summary.csv"), row.names = FALSE)
write.csv(gse96058_tnbc_candidates, file.path(table_dir, "GSE96058_pathology_TNBC_primary_samples.csv"), row.names = FALSE)
write.csv(gse58812_survival, file.path(table_dir, "GSE58812_TNBC_survival_metadata.csv"), row.names = FALSE)

capture.output(sessionInfo(), file = file.path(log_dir, "sessionInfo_05_summarize_bulk_cohort_metadata.txt"))

cat("Bulk cohort metadata summary complete\n")
print(cohort_summary)
cat("\nGSE96058 PAM50 primary samples:\n")
print(pam50_summary)
