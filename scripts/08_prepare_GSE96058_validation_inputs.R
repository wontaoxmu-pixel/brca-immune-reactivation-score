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

expression_file <- file.path(raw_dir, "GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz")

to_num <- function(x) suppressWarnings(as.numeric(x))
to_int <- function(x) suppressWarnings(as.integer(x))
not_missing <- function(x) !is.na(x) & x != ""

meta <- fread(metadata_path)
gse96058 <- copy(meta[dataset == "GSE96058"])
if (nrow(gse96058) == 0) {
  stop("No GSE96058 rows found in metadata: ", metadata_path)
}

gse96058[, `:=`(
  is_technical_replicate = grepl("repl$", sample_title),
  primary_sample_id = sub("repl$", "", sample_title),
  age_at_diagnosis_num = to_num(age_at_diagnosis),
  tumor_size_num = to_num(tumor_size),
  er_status_num = to_int(er_status),
  pgr_status_num = to_int(pgr_status),
  her2_status_num = to_int(her2_status),
  os_days_num = to_num(overall_survival_days),
  os_event = to_int(overall_survival_event),
  endocrine_treated_num = to_int(endocrine_treated),
  chemo_treated_num = to_int(chemo_treated)
)]

gse96058[, pathology_complete := !is.na(er_status_num) & !is.na(pgr_status_num) & !is.na(her2_status_num)]
gse96058[, tnbc_pathology := pathology_complete & er_status_num == 0 & pgr_status_num == 0 & her2_status_num == 0]
gse96058[, basal_pam50 := not_missing(pam50_subtype) & pam50_subtype == "Basal"]
gse96058[, tnbc_or_basal := tnbc_pathology | basal_pam50]
gse96058[, os_complete := !is.na(os_days_num) & os_days_num > 0 & os_event %in% c(0L, 1L)]

primary_samples <- gse96058[is_technical_replicate == FALSE]
setorder(primary_samples, primary_sample_id, platform, geo_accession)

validation_sample_set <- primary_samples[
  ,
  .(
    dataset,
    platform,
    primary_sample_id,
    sample_title,
    geo_accession,
    scan_b_external_id,
    age_at_diagnosis = age_at_diagnosis_num,
    tumor_size = tumor_size_num,
    lymph_node_status,
    lymph_node_group,
    er_status = er_status_num,
    pgr_status = pgr_status_num,
    her2_status = her2_status_num,
    pathology_complete,
    tnbc_pathology,
    pam50_subtype,
    basal_pam50,
    tnbc_or_basal,
    os_days = os_days_num,
    os_event,
    os_complete,
    endocrine_treated = endocrine_treated_num,
    chemo_treated = chemo_treated_num
  )
]

make_design_row <- function(dt, label) {
  data.table(
    analysis_set = label,
    n = nrow(dt),
    os_complete = sum(dt$os_complete, na.rm = TRUE),
    os_events = sum(dt$os_event == 1 & dt$os_complete, na.rm = TRUE),
    pathology_complete = sum(dt$pathology_complete, na.rm = TRUE),
    tnbc_pathology = sum(dt$tnbc_pathology, na.rm = TRUE),
    basal_pam50 = sum(dt$basal_pam50, na.rm = TRUE),
    median_age = median(dt$age_at_diagnosis, na.rm = TRUE),
    age_available = sum(!is.na(dt$age_at_diagnosis)),
    tumor_size_available = sum(!is.na(dt$tumor_size)),
    lymph_node_available = sum(not_missing(dt$lymph_node_status)),
    endocrine_treated_available = sum(!is.na(dt$endocrine_treated)),
    chemo_treated_available = sum(!is.na(dt$chemo_treated))
  )
}

design_summary <- rbindlist(list(
  make_design_row(validation_sample_set, "all_primary"),
  make_design_row(validation_sample_set[tnbc_pathology == TRUE], "pathology_tnbc"),
  make_design_row(validation_sample_set[basal_pam50 == TRUE], "pam50_basal"),
  make_design_row(validation_sample_set[tnbc_pathology == TRUE & basal_pam50 == TRUE], "pathology_tnbc_and_pam50_basal"),
  make_design_row(validation_sample_set[tnbc_or_basal == TRUE], "pathology_tnbc_or_pam50_basal")
), fill = TRUE)

replicate_summary <- gse96058[
  ,
  .(
    gsm_rows = .N,
    technical_replicates = sum(is_technical_replicate),
    primary_samples = uniqueN(primary_sample_id),
    duplicated_primary_ids = sum(duplicated(primary_sample_id))
  )
]

expression_status <- data.table(
  expected_expression_file = expression_file,
  expression_file_available = file.exists(expression_file),
  expression_file_size_mb = if (file.exists(expression_file)) {
    round(file.info(expression_file)$size / 1024^2, 2)
  } else {
    NA_real_
  },
  note = if (file.exists(expression_file)) {
    "Ready for GSE96058 signature scoring."
  } else {
    "Expression matrix is deferred because the GEO file is >500 MB; run scoring only after explicit download confirmation."
  }
)

write.csv(validation_sample_set, file.path(table_dir, "GSE96058_validation_sample_set.csv"), row.names = FALSE)
write.csv(design_summary, file.path(table_dir, "GSE96058_validation_design_summary.csv"), row.names = FALSE)
write.csv(replicate_summary, file.path(table_dir, "GSE96058_replicate_summary.csv"), row.names = FALSE)
write.csv(expression_status, file.path(table_dir, "GSE96058_expression_file_status.csv"), row.names = FALSE)

capture.output(sessionInfo(), file = file.path(log_dir, "sessionInfo_08_prepare_GSE96058_validation_inputs.txt"))

cat("GSE96058 validation input preparation complete\n")
print(replicate_summary)
cat("\nDesign summary:\n")
print(design_summary)
cat("\nExpression file status:\n")
print(expression_status)
