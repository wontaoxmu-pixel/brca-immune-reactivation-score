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
  library(survival)
  library(survminer)
  library(ggplot2)
})

raw_dir <- file.path(project_root, "data", "raw")
processed_dir <- file.path(project_root, "data", "processed")
table_dir <- file.path(project_root, "results", "tables")
figure_dir <- file.path(project_root, "results", "figures")
log_dir <- file.path(project_root, "logs")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

expr_path <- file.path(raw_dir, "GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz")
sample_path <- file.path(table_dir, "GSE96058_validation_sample_set.csv")
if (!file.exists(expr_path)) {
  stop("Missing GSE96058 expression file: ", expr_path)
}
if (!file.exists(sample_path)) {
  stop("Missing GSE96058 validation sample set. Run scripts/08_prepare_GSE96058_validation_inputs.R first.")
}

signature_sets <- list(
  exhausted_cd8_t_cell = c("CD8A", "CD8B", "PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX", "CXCL13", "GZMB", "PRF1", "NKG7"),
  antigen_presentation = c("HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "TAPBP", "PSMB8", "PSMB9", "NLRC5"),
  interferon_response = c("IFNG", "STAT1", "IRF1", "CXCL9", "CXCL10", "GBP1", "GBP5", "ISG15", "IFIT1", "MX1")
)
signature_genes <- unique(toupper(unlist(signature_sets, use.names = FALSE)))

sample_set <- fread(sample_path)
required_sample_cols <- c(
  "sample_title",
  "primary_sample_id",
  "os_days",
  "os_event",
  "os_complete",
  "age_at_diagnosis",
  "tumor_size",
  "tnbc_pathology",
  "basal_pam50",
  "tnbc_or_basal",
  "pam50_subtype",
  "lymph_node_status",
  "endocrine_treated",
  "chemo_treated"
)
missing_sample_cols <- setdiff(required_sample_cols, names(sample_set))
if (length(missing_sample_cols) > 0) {
  stop("Validation sample set missing columns: ", paste(missing_sample_cols, collapse = ", "))
}

expr_all <- fread(expr_path, check.names = FALSE)
setnames(expr_all, 1, "gene_symbol")
expr_all[, gene_symbol := toupper(gene_symbol)]
expr_sig <- expr_all[gene_symbol %in% signature_genes]
rm(expr_all)

available_samples <- intersect(sample_set$sample_title, setdiff(names(expr_sig), "gene_symbol"))
if (length(available_samples) == 0) {
  stop("No validation samples from GSE96058 metadata were found in the expression matrix columns.")
}

for (col in available_samples) {
  set(expr_sig, j = col, value = as.numeric(expr_sig[[col]]))
}

gene_expr <- expr_sig[
  ,
  lapply(.SD, mean, na.rm = TRUE),
  by = gene_symbol,
  .SDcols = available_samples
]

mat <- as.matrix(gene_expr[, ..available_samples])
mode(mat) <- "numeric"
rownames(mat) <- gene_expr$gene_symbol
z_mat <- t(scale(t(mat)))
z_mat[is.na(z_mat)] <- 0

score_dt <- data.table(sample_title = colnames(z_mat))
coverage_rows <- list()

for (set_name in names(signature_sets)) {
  genes <- unique(toupper(signature_sets[[set_name]]))
  present <- intersect(genes, rownames(z_mat))
  missing <- setdiff(genes, rownames(z_mat))
  coverage_rows[[length(coverage_rows) + 1]] <- data.table(
    signature = set_name,
    requested_genes = length(genes),
    present_genes = length(present),
    missing_genes = paste(missing, collapse = ";"),
    present_gene_list = paste(present, collapse = ";")
  )
  score_dt[, (set_name) := colMeans(z_mat[present, , drop = FALSE], na.rm = TRUE)]
}

score_cols <- names(signature_sets)
score_dt[, immune_reactivation_score := rowMeans(.SD, na.rm = TRUE), .SDcols = score_cols]
score_dt[, immune_reactivation_score_z := as.numeric(scale(immune_reactivation_score))]

coverage_dt <- rbindlist(coverage_rows)
mapping_summary <- data.table(
  dataset = "GSE96058",
  expression_file = basename(expr_path),
  rows_read = nrow(expr_sig),
  unique_signature_gene_rows = uniqueN(expr_sig$gene_symbol),
  gene_symbols_after_collapse = nrow(gene_expr),
  expression_samples_available = length(setdiff(names(expr_sig), "gene_symbol")),
  validation_primary_samples = nrow(sample_set),
  validation_samples_in_expression = length(available_samples),
  unmatched_validation_samples = sum(!sample_set$sample_title %in% available_samples)
)

analysis_dt <- merge(sample_set, score_dt, by = "sample_title", all = FALSE)
analysis_dt[, immune_reactivation_group := fifelse(
  immune_reactivation_score >= median(immune_reactivation_score, na.rm = TRUE),
  "High",
  "Low"
)]
analysis_dt[, immune_reactivation_group := factor(immune_reactivation_group, levels = c("Low", "High"))]

run_cox <- function(dt, analysis_set_label) {
  endpoint_dt <- dt[
    os_complete == TRUE &
      !is.na(os_days) &
      os_days > 0 &
      os_event %in% c(0L, 1L) &
      !is.na(immune_reactivation_score_z)
  ]
  if (nrow(endpoint_dt) < 10 || sum(endpoint_dt$os_event == 1, na.rm = TRUE) < 5) {
    return(list(
      cox = data.table(),
      zph = data.table(),
      logrank = data.table(),
      group_summary = data.table()
    ))
  }

  endpoint_dt[, `:=`(
    endpoint_time = os_days,
    endpoint_event = os_event
  )]

  models <- list(
    cox_continuous_per_1sd = coxph(Surv(endpoint_time, endpoint_event) ~ immune_reactivation_score_z, data = endpoint_dt),
    cox_median_high_vs_low = coxph(Surv(endpoint_time, endpoint_event) ~ immune_reactivation_group, data = endpoint_dt)
  )

  if (all(!is.na(endpoint_dt$age_at_diagnosis))) {
    models$cox_continuous_per_1sd_age_adjusted <- coxph(
      Surv(endpoint_time, endpoint_event) ~ immune_reactivation_score_z + age_at_diagnosis,
      data = endpoint_dt
    )
    models$cox_median_high_vs_low_age_adjusted <- coxph(
      Surv(endpoint_time, endpoint_event) ~ immune_reactivation_group + age_at_diagnosis,
      data = endpoint_dt
    )
  }

  extract_model <- function(model_name, model) {
    model_summary <- summary(model)
    conf <- model_summary$conf.int
    coefs <- model_summary$coefficients
    data.table(
      analysis_set = analysis_set_label,
      endpoint = "overall_survival",
      model = model_name,
      term = rownames(coefs),
      n = model_summary$n,
      events = model_summary$nevent,
      hr = unname(conf[, "exp(coef)"]),
      ci_lower = unname(conf[, "lower .95"]),
      ci_upper = unname(conf[, "upper .95"]),
      p_value = unname(coefs[, "Pr(>|z|)"]),
      concordance = unname(model_summary$concordance[1])
    )
  }

  cox_rows <- rbindlist(Map(extract_model, names(models), models), fill = TRUE)

  zph_rows <- rbindlist(lapply(names(models), function(model_name) {
    zph <- cox.zph(models[[model_name]])
    data.table(
      analysis_set = analysis_set_label,
      endpoint = "overall_survival",
      model = model_name,
      term = rownames(zph$table),
      chisq = unname(zph$table[, "chisq"]),
      p_value = unname(zph$table[, "p"])
    )
  }), fill = TRUE)

  logrank_fit <- survdiff(Surv(endpoint_time, endpoint_event) ~ immune_reactivation_group, data = endpoint_dt)
  logrank_row <- data.table(
    analysis_set = analysis_set_label,
    endpoint = "overall_survival",
    test = "logrank_median_high_vs_low",
    n = nrow(endpoint_dt),
    events = sum(endpoint_dt$endpoint_event == 1, na.rm = TRUE),
    chisq = unname(logrank_fit$chisq),
    p_value = 1 - pchisq(logrank_fit$chisq, length(logrank_fit$n) - 1)
  )

  group_summary <- endpoint_dt[
    ,
    .(
      n = .N,
      events = as.integer(sum(endpoint_event == 1, na.rm = TRUE)),
      median_os_days = as.numeric(median(endpoint_time, na.rm = TRUE)),
      median_score = as.numeric(median(immune_reactivation_score, na.rm = TRUE))
    ),
    by = .(immune_reactivation_group)
  ]
  group_summary[, analysis_set := analysis_set_label]
  setcolorder(group_summary, c("analysis_set", setdiff(names(group_summary), "analysis_set")))

  km_fit <- survfit(Surv(endpoint_time, endpoint_event) ~ immune_reactivation_group, data = endpoint_dt)
  km_plot <- ggsurvplot(
    km_fit,
    data = endpoint_dt,
    pval = TRUE,
    risk.table = TRUE,
    conf.int = FALSE,
    palette = c("#3b6fb6", "#c43b3b"),
    legend.title = "Immune reactivation",
    legend.labs = c("Low", "High"),
    xlab = "Days",
    ylab = "Overall survival probability",
    risk.table.height = 0.28,
    ggtheme = theme_classic(base_size = 10)
  )

  ggsave(
    filename = file.path(figure_dir, paste0("GSE96058_", analysis_set_label, "_overall_survival_KM_median_score.png")),
    plot = km_plot$plot,
    width = 5.5,
    height = 4.2,
    dpi = 300
  )
  ggsave(
    filename = file.path(figure_dir, paste0("GSE96058_", analysis_set_label, "_overall_survival_KM_median_score_risktable.png")),
    plot = arrange_ggsurvplots(list(km_plot), print = FALSE),
    width = 5.8,
    height = 5.6,
    dpi = 300
  )

  list(cox = cox_rows, zph = zph_rows, logrank = logrank_row, group_summary = group_summary)
}

analysis_sets <- list(
  all_primary = analysis_dt,
  pathology_tnbc = analysis_dt[tnbc_pathology == TRUE],
  pam50_basal = analysis_dt[basal_pam50 == TRUE],
  pathology_tnbc_and_pam50_basal = analysis_dt[tnbc_pathology == TRUE & basal_pam50 == TRUE],
  pathology_tnbc_or_pam50_basal = analysis_dt[tnbc_or_basal == TRUE]
)

survival_results <- lapply(names(analysis_sets), function(label) run_cox(copy(analysis_sets[[label]]), label))
names(survival_results) <- names(analysis_sets)

cox_results <- rbindlist(lapply(survival_results, `[[`, "cox"), fill = TRUE)
zph_results <- rbindlist(lapply(survival_results, `[[`, "zph"), fill = TRUE)
logrank_results <- rbindlist(lapply(survival_results, `[[`, "logrank"), fill = TRUE)
group_summary <- rbindlist(lapply(survival_results, `[[`, "group_summary"), fill = TRUE)

write.csv(score_dt, file.path(table_dir, "GSE96058_immune_reactivation_scores.csv"), row.names = FALSE)
write.csv(coverage_dt, file.path(table_dir, "GSE96058_signature_gene_coverage.csv"), row.names = FALSE)
write.csv(mapping_summary, file.path(table_dir, "GSE96058_signature_mapping_summary.csv"), row.names = FALSE)
write.csv(analysis_dt, file.path(table_dir, "GSE96058_signature_survival_analysis_dataset.csv"), row.names = FALSE)
write.csv(cox_results, file.path(table_dir, "GSE96058_survival_cox_results.csv"), row.names = FALSE)
write.csv(zph_results, file.path(table_dir, "GSE96058_survival_cox_zph_results.csv"), row.names = FALSE)
write.csv(logrank_results, file.path(table_dir, "GSE96058_survival_logrank_results.csv"), row.names = FALSE)
write.csv(group_summary, file.path(table_dir, "GSE96058_survival_group_summary.csv"), row.names = FALSE)

capture.output(sessionInfo(), file = file.path(log_dir, "sessionInfo_09_score_GSE96058_signature_survival.txt"))

cat("GSE96058 signature scoring and survival analysis complete\n")
print(mapping_summary)
cat("\nCoverage:\n")
print(coverage_dt)
cat("\nCox results:\n")
print(cox_results)
