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

cox <- fread(file.path(table_dir, "signature_deconvolution_adjusted_cox.csv"))
cor <- fread(file.path(table_dir, "signature_deconvolution_score_correlations.csv"))

model_labels <- c(
  age_plus_purity_surrogate = "Age + purity surrogate",
  age_plus_deconv_immune_summary = "Age + immune-lineage summary",
  age_plus_stromal_plus_purity = "Age + stromal summary + purity surrogate",
  age_plus_immune_stromal_purity = "Age + immune-lineage summary + stromal summary + purity surrogate"
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
  model %in% c(
    "age_plus_purity_surrogate",
    "age_plus_deconv_immune_summary",
    "age_plus_immune_stromal_purity"
  ) &
    (
      (dataset == "GSE58812" & endpoint %in% c("overall_survival", "metastasis_free_survival")) |
        (dataset == "GSE96058" & analysis_set %in% c("all_primary", "pam50_basal", "pathology_tnbc"))
    )
]

selected[, `:=`(
  cohort = fifelse(
    dataset == "GSE58812",
    analysis_labels[analysis_set],
    analysis_labels[analysis_set]
  ),
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

setorder(selected, Cohort, Endpoint, Model)

summary_cor <- cor[
  deconvolution_variable %in% c(
    "deconv_immune_summary_z",
    "deconv_stromal_summary_z",
    "deconv_epithelial_summary_z",
    "purity_surrogate_z"
  )
]
summary_cor[, `:=`(
  cohort = paste(dataset, analysis_set, sep = " "),
  variable = fifelse(
    deconvolution_variable == "deconv_immune_summary_z", "Immune-lineage summary",
    fifelse(
      deconvolution_variable == "deconv_stromal_summary_z", "Stromal summary",
      fifelse(
        deconvolution_variable == "deconv_epithelial_summary_z", "Epithelial summary",
        "Purity surrogate"
      )
    )
  ),
  rho = sprintf("%.3f", spearman_rho),
  p_value = format_p(p_value)
)]
summary_cor <- summary_cor[
  ,
  .(
    Cohort = cohort,
    Variable = variable,
    `n` = n,
    `Spearman rho` = rho,
    `P value` = p_value
  )
]

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

write_md_table(
  selected,
  file.path(manuscript_table_dir, "Table_5_Signature_deconvolution_purity_sensitivity.md"),
  "Table 5. Lineage-signature and purity-surrogate sensitivity Cox models",
  "The immune reactivation score HR is reported per 1 SD increment. These are marker/signature-based sensitivity models, not package-level ESTIMATE, MCP-counter, xCell, CIBERSORT, or measured tumor-purity outputs."
)

write_md_table(
  summary_cor,
  file.path(manuscript_table_dir, "Table_S5_Signature_deconvolution_correlations.md"),
  "Table S5. Correlations between the immune reactivation score and lineage-signature/purity surrogates",
  "Spearman correlations were calculated between the standardized immune reactivation score and signature-based immune, stromal, epithelial, and purity-surrogate scores."
)

writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_16_make_deconvolution_manuscript_table.txt"))
cat("Deconvolution manuscript tables written\n")
