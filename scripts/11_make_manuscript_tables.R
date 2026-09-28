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

required_files <- file.path(
  table_dir,
  c(
    "parsed_bulk_data_summary.csv",
    "GSE96058_validation_design_summary.csv",
    "GSE176078_signature_gene_coverage.csv",
    "GSE58812_signature_gene_coverage.csv",
    "GSE96058_signature_gene_coverage.csv",
    "bulk_validation_continuous_cox_publication_table.csv",
    "microenvironment_proxy_gene_coverage.csv",
    "microenvironment_adjusted_cox.csv",
    "microenvironment_proxy_correlations.csv"
  )
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required input files:\n", paste(missing_files, collapse = "\n"))
}

read_csv <- function(name) {
  fread(file.path(table_dir, name))
}

write_markdown_table <- function(dt, path, title, note = NULL) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(paste0("# ", title), con)
  writeLines("", con)
  if (!is.null(note) && nzchar(note)) {
    writeLines(note, con)
    writeLines("", con)
  }
  headers <- names(dt)
  writeLines(paste0("| ", paste(headers, collapse = " | "), " |"), con)
  writeLines(paste0("| ", paste(rep("---", length(headers)), collapse = " | "), " |"), con)
  for (i in seq_len(nrow(dt))) {
    vals <- vapply(dt[i], as.character, character(1))
    vals[is.na(vals)] <- ""
    vals <- gsub("\\|", "/", vals)
    writeLines(paste0("| ", paste(vals, collapse = " | "), " |"), con)
  }
}

parsed <- read_csv("parsed_bulk_data_summary.csv")
gse96058_design <- read_csv("GSE96058_validation_design_summary.csv")

table1 <- data.table(
  Dataset = c("GSE176078", "GSE58812", "GSE96058"),
  Role = c(
    "Feasibility scoring",
    "TNBC external survival validation",
    "Large RNA-seq OS validation"
  ),
  Platform_or_input = c(
    "Matched bulk RNA-seq supplementary matrix",
    "GPL570 microarray series matrix",
    "GPL11154/GPL18573 metadata plus transformed RNA-seq matrix"
  ),
  Analyzed_samples = c(
    parsed[dataset == "GSE176078", as.character(samples)][1],
    "107",
    gse96058_design[analysis_set == "all_primary", as.character(n)][1]
  ),
  Subtype_focus = c(
    "Breast cancer matched bulk feasibility set",
    "Pathology TNBC",
    "All primary; pathology TNBC; PAM50 Basal; TNBC/Basal combined subsets"
  ),
  Endpoint = c(
    "Not used for survival",
    "OS and MFS",
    "OS"
  ),
  Events = c(
    "NA",
    "OS 29; MFS 31",
    paste0(
      "All primary ",
      gse96058_design[analysis_set == "all_primary", os_events][1],
      "; pathology TNBC ",
      gse96058_design[analysis_set == "pathology_tnbc", os_events][1],
      "; PAM50 Basal ",
      gse96058_design[analysis_set == "pam50_basal", os_events][1],
      "; pathology TNBC + PAM50 Basal ",
      gse96058_design[analysis_set == "pathology_tnbc_and_pam50_basal", os_events][1],
      "; pathology TNBC or PAM50 Basal ",
      gse96058_design[analysis_set == "pathology_tnbc_or_pam50_basal", os_events][1]
    )
  )
)

write_markdown_table(
  table1,
  file.path(manuscript_table_dir, "Table_1_Dataset_inventory.md"),
  "Table 1. Dataset inventory and validation design",
  "NA indicates not applicable. TNBC, triple-negative breast cancer; OS, overall survival; MFS, metastasis-free survival."
)

coverage_files <- c(
  GSE176078 = "GSE176078_signature_gene_coverage.csv",
  GSE58812 = "GSE58812_signature_gene_coverage.csv",
  GSE96058 = "GSE96058_signature_gene_coverage.csv"
)
coverage <- rbindlist(
  lapply(names(coverage_files), function(dataset) {
    dt <- read_csv(coverage_files[[dataset]])
    dt[, Dataset := dataset]
    dt
  }),
  fill = TRUE
)

signature_labels <- c(
  exhausted_cd8_t_cell = "Exhausted/cytotoxic CD8+ T-cell",
  antigen_presentation = "Antigen presentation",
  interferon_response = "Interferon response"
)
coverage[, Component := unname(signature_labels[signature])]
coverage[is.na(Component), Component := signature]
coverage[, Coverage := paste0(present_genes, "/", requested_genes)]
coverage[, missing_genes := as.character(missing_genes)]
coverage[, Missing_genes := fifelse(is.na(missing_genes) | missing_genes == "", "None", missing_genes)]
table2 <- coverage[
  ,
  .(Dataset, Component, Coverage, Missing_genes, Present_genes = present_gene_list)
]
setorder(table2, Dataset, Component)

write_markdown_table(
  table2,
  file.path(manuscript_table_dir, "Table_2_Signature_gene_coverage.md"),
  "Table 2. Signature gene coverage across datasets",
  "Coverage is reported as present/requested genes after dataset-specific gene-symbol mapping and collapse."
)

cox_publication <- read_csv("bulk_validation_continuous_cox_publication_table.csv")
setnames(
  cox_publication,
  old = c("cohort", "endpoint", "model", "n", "events", "HR (95% CI)", "P value", "concordance", "direction"),
  new = c("Cohort", "Endpoint", "Model", "n", "Events", "HR_95_CI", "P_value", "Concordance", "Direction")
)
cox_publication[, Concordance := sprintf("%.3f", as.numeric(Concordance))]
cox_publication[, P_value := {
  p_num <- suppressWarnings(as.numeric(P_value))
  fifelse(is.na(p_num), as.character(P_value), sprintf("%.3f", p_num))
}]

write_markdown_table(
  cox_publication,
  file.path(manuscript_table_dir, "Table_3_Continuous_Cox_validation.md"),
  "Table 3. Continuous Cox validation results",
  "Hazard ratios are reported per 1 SD increase in the immune reactivation score. Age-adjusted models were fitted when age was available."
)

proxy_coverage <- read_csv("microenvironment_proxy_gene_coverage.csv")
proxy_cox <- read_csv("microenvironment_adjusted_cox.csv")
proxy_cor <- read_csv("microenvironment_proxy_correlations.csv")

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
proxy_labels <- c(
  leukocyte_proxy = "Leukocyte proxy",
  stromal_proxy = "Stromal proxy",
  epithelial_proxy = "Epithelial proxy"
)

proxy_cox[, proxy := sub("^age_plus_", "", model)]
proxy_cox[, proxy := sub("_proxy$", "_proxy", proxy)]
proxy_cox[, Analysis_set := unname(analysis_labels[analysis_set])]
proxy_cox[is.na(Analysis_set), Analysis_set := analysis_set]
proxy_cox[, Endpoint := unname(endpoint_labels[endpoint])]
proxy_cox[is.na(Endpoint), Endpoint := endpoint]
proxy_cox[, Proxy := unname(proxy_labels[proxy])]
proxy_cox[is.na(Proxy), Proxy := proxy]
proxy_cox[, HR_95_CI := hr_95_ci]
proxy_cox[, P_value := p_value_formatted]

proxy_cor[, endpoint := fifelse(
  analysis_set == "tnbc_os_rows",
  "overall_survival",
  fifelse(analysis_set == "tnbc_mfs_rows", "metastasis_free_survival", "overall_survival")
)]
proxy_cor[, analysis_set_join := fifelse(
  analysis_set %chin% c("tnbc_os_rows", "tnbc_mfs_rows"),
  "tnbc",
  analysis_set
)]
proxy_cor[, proxy_label := unname(proxy_labels[proxy])]
proxy_cor[is.na(proxy_label), proxy_label := proxy]
proxy_cor[, Spearman_rho := fifelse(
  is.na(spearman_rho),
  NA_character_,
  sprintf("%.3f", spearman_rho)
)]
proxy_cor_small <- proxy_cor[
  ,
  .(
    analysis_set = analysis_set_join,
    endpoint,
    Proxy = proxy_label,
    Spearman_rho,
    Correlation_P = p_value_formatted,
    Correlation_status = status
  )
]

table4 <- merge(
  proxy_cox[
    ,
    .(
      analysis_set,
      endpoint,
      Analysis_set,
      Endpoint,
      Proxy,
      n,
      Events = events,
      HR_95_CI,
      P_value
    )
  ],
  proxy_cor_small,
  by = c("analysis_set", "endpoint", "Proxy"),
  all.x = TRUE
)
table4 <- table4[
  ,
  .(
    Analysis_set,
    Endpoint,
    Proxy,
    n,
    Events,
    Spearman_rho,
    Correlation_P,
    HR_95_CI,
    P_value
  )
]
setorder(table4, Analysis_set, Endpoint, Proxy)

coverage_note <- paste(
  proxy_coverage[
    ,
    paste0(dataset, " ", proxy, " ", present_genes, "/", requested_genes)
  ],
  collapse = "; "
)
write_markdown_table(
  table4,
  file.path(manuscript_table_dir, "Table_4_Microenvironment_proxy_sensitivity.md"),
  "Table 4. Marker-based microenvironment proxy sensitivity analysis",
  paste0(
    "Spearman rho reports correlation between the immune reactivation score and each marker-based proxy. ",
    "Hazard ratios are reported per 1 SD increase in the immune reactivation score after age plus proxy adjustment. ",
    "Marker coverage: ", coverage_note, "."
  )
)

writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_11_make_manuscript_tables.txt"))

message("Wrote manuscript tables to: ", manuscript_table_dir)
