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

table_dir <- file.path(project_root, "results", "tables")
manuscript_table_dir <- file.path(project_root, "manuscript", "tables")
log_dir <- file.path(project_root, "logs")
dir.create(manuscript_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

format_p <- function(x) {
  out <- rep(NA_character_, length(x))
  out[!is.na(x) & x < 0.001] <- "<0.001"
  out[!is.na(x) & x >= 0.001] <- sprintf("%.3f", x[!is.na(x) & x >= 0.001])
  out
}

write_md_table <- function(dt, path, title, note) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(paste0("# ", title), con)
  writeLines("", con)
  writeLines(note, con)
  writeLines("", con)
  header <- paste(names(dt), collapse = " | ")
  sep <- paste(rep("---", ncol(dt)), collapse = " | ")
  writeLines(paste0("| ", header, " |"), con)
  writeLines(paste0("| ", sep, " |"), con)
  for (i in seq_len(nrow(dt))) {
    vals <- vapply(dt[i], as.character, character(1))
    writeLines(paste0("| ", paste(vals, collapse = " | "), " |"), con)
  }
}

cox <- fread(file.path(table_dir, "estimate_adjusted_cox.csv"))
cor <- fread(file.path(table_dir, "estimate_score_correlations.csv"))

model_labels <- c(
  age_plus_estimate_immune = "Age + ESTIMATE ImmuneScore",
  age_plus_estimate_stromal = "Age + ESTIMATE StromalScore",
  age_plus_estimate_score = "Age + ESTIMATEScore",
  age_plus_estimate_purity = "Age + ESTIMATE-derived TumorPurity",
  age_plus_estimate_immune_stromal_purity = "Age + ImmuneScore + StromalScore + TumorPurity"
)

analysis_labels <- c(
  tnbc = "GSE58812 TNBC",
  all_primary = "GSE96058 all primary",
  pam50_basal = "GSE96058 PAM50 Basal",
  pathology_tnbc = "GSE96058 pathology TNBC",
  pathology_tnbc_or_pam50_basal = "GSE96058 pathology TNBC or PAM50 Basal"
)

endpoint_labels <- c(
  overall_survival = "OS",
  metastasis_free_survival = "MFS"
)

selected <- cox[
  model %in% names(model_labels) &
    (
      (dataset == "GSE58812" & endpoint %in% c("overall_survival", "metastasis_free_survival")) |
        (dataset == "GSE96058" & analysis_set %in% c(
          "all_primary",
          "pam50_basal",
          "pathology_tnbc",
          "pathology_tnbc_or_pam50_basal"
        ))
    )
]

selected[, `:=`(
  cohort = analysis_labels[analysis_set],
  endpoint_label = endpoint_labels[endpoint],
  model_label = model_labels[model],
  p_value = format_p(p_value),
  hr_ci = sprintf("%.3f (%.3f-%.3f)", hr, ci_lower, ci_upper)
)]

selected <- selected[
  ,
  .(
    Cohort = cohort,
    Endpoint = endpoint_label,
    Model = model_label,
    `n` = n,
    Events = events,
    `HR (95% CI)` = hr_ci,
    `P value` = p_value
  )
]

selected[, cohort_order := match(Cohort, c(
  "GSE58812 TNBC",
  "GSE96058 all primary",
  "GSE96058 PAM50 Basal",
  "GSE96058 pathology TNBC",
  "GSE96058 pathology TNBC or PAM50 Basal"
))]
selected[, endpoint_order := match(Endpoint, c("OS", "MFS"))]
selected[, model_order := match(Model, model_labels)]
setorder(selected, cohort_order, endpoint_order, model_order)
selected[, c("cohort_order", "endpoint_order", "model_order") := NULL]

variable_labels <- c(
  StromalScore_z = "ESTIMATE StromalScore",
  ImmuneScore_z = "ESTIMATE ImmuneScore",
  ESTIMATEScore_z = "ESTIMATEScore",
  TumorPurity_z = "ESTIMATE-derived TumorPurity"
)

summary_cor <- cor[variable %in% names(variable_labels)]
summary_cor[, `:=`(
  cohort = analysis_labels[analysis_set],
  variable_label = variable_labels[variable],
  rho = sprintf("%.3f", spearman_rho),
  p_value = format_p(p_value)
)]
summary_cor <- summary_cor[
  ,
  .(
    Cohort = cohort,
    Variable = variable_label,
    `n` = n,
    `Spearman rho` = rho,
    `P value` = p_value
  )
]
summary_cor[, cohort_order := match(Cohort, c(
  "GSE58812 TNBC",
  "GSE96058 all primary",
  "GSE96058 PAM50 Basal",
  "GSE96058 pathology TNBC",
  "GSE96058 pathology TNBC or PAM50 Basal"
))]
summary_cor[, variable_order := match(Variable, variable_labels)]
setorder(summary_cor, cohort_order, variable_order)
summary_cor[, c("cohort_order", "variable_order") := NULL]

write_md_table(
  selected,
  file.path(manuscript_table_dir, "Table_6_ESTIMATE_purity_sensitivity.md"),
  "Table 6. ESTIMATE package-level immune/stromal/purity sensitivity Cox models",
  "The immune reactivation score HR is reported per 1 SD increment. ESTIMATE scores were generated with the `estimate` R package. GSE96058 TumorPurity was derived from ESTIMATEScore using the published ESTIMATE formula because the non-Affymetrix ESTIMATE run did not emit TumorPurity directly; it should not be interpreted as histology-measured tumor purity."
)

write_md_table(
  summary_cor,
  file.path(manuscript_table_dir, "Table_S6_ESTIMATE_correlations.md"),
  "Table S6. Correlations between the immune reactivation score and ESTIMATE scores",
  "Spearman correlations were calculated between the standardized immune reactivation score and standardized ESTIMATE StromalScore, ImmuneScore, ESTIMATEScore, and ESTIMATE-derived TumorPurity."
)

writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_19_make_estimate_manuscript_table.txt"))
cat("ESTIMATE manuscript tables written\n")
