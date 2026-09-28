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
  library(ggplot2)
})

raw_dir <- file.path(project_root, "data", "raw")
processed_dir <- file.path(project_root, "data", "processed")
table_dir <- file.path(project_root, "results", "tables")
figure_dir <- file.path(project_root, "results", "figures")
log_dir <- file.path(project_root, "logs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

format_p <- function(x) {
  out <- rep(NA_character_, length(x))
  out[!is.na(x) & x < 0.001] <- "<0.001"
  out[!is.na(x) & x >= 0.001] <- sprintf("%.3f", x[!is.na(x) & x >= 0.001])
  out
}

zscore <- function(x) {
  out <- as.numeric(scale(x))
  out[is.na(out)] <- 0
  out
}

mean_score <- function(z_mat, genes) {
  genes <- intersect(toupper(genes), rownames(z_mat))
  if (length(genes) == 0) {
    return(rep(NA_real_, ncol(z_mat)))
  }
  colMeans(z_mat[genes, , drop = FALSE], na.rm = TRUE)
}

signature_sets <- list(
  mcp_like_t_cells = c("CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "TRBC2"),
  mcp_like_cd8_t_cells = c("CD8A", "CD8B", "GZMK", "GZMB", "NKG7", "PRF1", "CCL5", "CXCL13"),
  mcp_like_cytotoxic_lymphocytes = c("NKG7", "GNLY", "PRF1", "GZMA", "GZMB", "GZMH", "KLRD1", "CCL5"),
  mcp_like_nk_cells = c("KLRD1", "NCR1", "FCGR3A", "TYROBP", "GNLY", "NKG7"),
  mcp_like_b_lineage = c("MS4A1", "CD79A", "CD79B", "BANK1", "CD19", "JCHAIN", "MZB1"),
  mcp_like_myeloid = c("LYZ", "LST1", "AIF1", "TYROBP", "FCER1G", "C1QA", "C1QB", "C1QC", "CD68"),
  mcp_like_dendritic_cells = c("ITGAX", "FCER1A", "CLEC10A", "LILRA4", "IRF8", "BATF3", "CLEC9A"),
  mcp_like_neutrophils = c("S100A8", "S100A9", "FCGR3B", "CSF3R", "CXCR2", "MMP9"),
  stromal_fibroblast = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FAP", "PDGFRA", "PDGFRB", "ACTA2", "TAGLN"),
  stromal_endothelial = c("PECAM1", "VWF", "ENG", "KDR", "ESAM", "CLDN5", "RAMP2", "ACKR1"),
  epithelial_tumor = c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT5", "KRT14", "KRT17", "MUC1", "CDH1", "ERBB2"),
  checkpoint_exhaustion = c("PDCD1", "CD274", "PDCD1LG2", "LAG3", "HAVCR2", "TIGIT", "CTLA4", "TOX"),
  antigen_presentation = c("HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "TAPBP", "PSMB8", "PSMB9", "NLRC5"),
  interferon_response = c("IFNG", "STAT1", "IRF1", "CXCL9", "CXCL10", "GBP1", "GBP5", "ISG15", "IFIT1", "MX1")
)

immune_signature_names <- c(
  "mcp_like_t_cells",
  "mcp_like_cd8_t_cells",
  "mcp_like_cytotoxic_lymphocytes",
  "mcp_like_nk_cells",
  "mcp_like_b_lineage",
  "mcp_like_myeloid",
  "mcp_like_dendritic_cells",
  "mcp_like_neutrophils"
)
stromal_signature_names <- c("stromal_fibroblast", "stromal_endothelial")
epithelial_signature_names <- "epithelial_tumor"

target_genes <- unique(toupper(unlist(signature_sets, use.names = FALSE)))

prepare_expression <- function(expr_dt, gene_col, sample_cols, dataset, transform = c("none", "log_cpm")) {
  transform <- match.arg(transform)
  expr_dt <- copy(expr_dt)
  sample_cols <- setdiff(sample_cols, gene_col)
  setnames(expr_dt, gene_col, "gene_symbol")
  expr_dt[, gene_symbol := toupper(gene_symbol)]
  expr_dt <- expr_dt[gene_symbol %in% target_genes]
  if (nrow(expr_dt) == 0) {
    stop("No deconvolution marker genes found in ", dataset)
  }
  for (col in sample_cols) {
    set(expr_dt, j = col, value = suppressWarnings(as.numeric(expr_dt[[col]])))
  }
  collapsed <- expr_dt[
    ,
    lapply(.SD, mean, na.rm = TRUE),
    by = gene_symbol,
    .SDcols = sample_cols
  ]
  mat <- as.matrix(collapsed[, ..sample_cols])
  mode(mat) <- "numeric"
  rownames(mat) <- collapsed$gene_symbol
  if (transform == "log_cpm") {
    library_sizes <- colSums(mat, na.rm = TRUE)
    library_sizes[library_sizes <= 0 | is.na(library_sizes)] <- 1
    mat <- log2(t(t(mat) / library_sizes * 1e6) + 1)
  }
  z_mat <- t(scale(t(mat)))
  z_mat[is.na(z_mat)] <- 0
  list(mat = mat, z_mat = z_mat)
}

score_deconvolution <- function(prepared, dataset) {
  z_mat <- prepared$z_mat
  scores <- data.table(sample_id = colnames(z_mat), dataset = dataset)
  coverage <- rbindlist(lapply(names(signature_sets), function(set_name) {
    requested <- unique(toupper(signature_sets[[set_name]]))
    present <- intersect(requested, rownames(z_mat))
    missing <- setdiff(requested, rownames(z_mat))
    scores[, (set_name) := mean_score(z_mat, present)]
    scores[, (paste0(set_name, "_z")) := zscore(get(set_name))]
    data.table(
      dataset = dataset,
      signature = set_name,
      requested_genes = length(requested),
      present_genes = length(present),
      missing_genes = paste(missing, collapse = ";"),
      present_gene_list = paste(present, collapse = ";")
    )
  }), fill = TRUE)

  scores[, deconv_immune_summary := rowMeans(.SD, na.rm = TRUE), .SDcols = immune_signature_names]
  scores[, deconv_stromal_summary := rowMeans(.SD, na.rm = TRUE), .SDcols = stromal_signature_names]
  scores[, deconv_epithelial_summary := rowMeans(.SD, na.rm = TRUE), .SDcols = epithelial_signature_names]
  scores[, deconv_immune_summary_z := zscore(deconv_immune_summary)]
  scores[, deconv_stromal_summary_z := zscore(deconv_stromal_summary)]
  scores[, deconv_epithelial_summary_z := zscore(deconv_epithelial_summary)]
  scores[, purity_surrogate := deconv_epithelial_summary_z - rowMeans(.SD, na.rm = TRUE), .SDcols = c("deconv_immune_summary_z", "deconv_stromal_summary_z")]
  scores[, purity_surrogate_z := zscore(purity_surrogate)]
  scores[, immune_to_purity_axis := deconv_immune_summary_z - purity_surrogate_z]
  scores[, immune_to_purity_axis_z := zscore(immune_to_purity_axis)]
  setcolorder(scores, c("dataset", "sample_id", setdiff(names(scores), c("dataset", "sample_id"))))
  list(scores = scores, coverage = coverage)
}

correlate_with_score <- function(dt, dataset, analysis_set) {
  vars <- c(
    paste0(names(signature_sets), "_z"),
    "deconv_immune_summary_z",
    "deconv_stromal_summary_z",
    "deconv_epithelial_summary_z",
    "purity_surrogate_z",
    "immune_to_purity_axis_z"
  )
  rbindlist(lapply(vars, function(var) {
    if (!var %in% names(dt)) {
      return(NULL)
    }
    keep <- complete.cases(dt[, .(immune_reactivation_score_z, get(var))])
    if (sum(keep) < 3) {
      return(data.table(
        dataset = dataset,
        analysis_set = analysis_set,
        deconvolution_variable = var,
        n = sum(keep),
        spearman_rho = NA_real_,
        p_value = NA_real_,
        p_value_formatted = NA_character_
      ))
    }
    test <- suppressWarnings(cor.test(dt$immune_reactivation_score_z[keep], dt[[var]][keep], method = "spearman"))
    data.table(
      dataset = dataset,
      analysis_set = analysis_set,
      deconvolution_variable = var,
      n = sum(keep),
      spearman_rho = unname(test$estimate),
      p_value = test$p.value,
      p_value_formatted = format_p(test$p.value)
    )
  }), fill = TRUE)
}

extract_score_term <- function(model_name, fit, dataset, analysis_set, endpoint) {
  model_summary <- summary(fit)
  conf <- model_summary$conf.int
  coefs <- model_summary$coefficients
  out <- data.table(
    dataset = dataset,
    analysis_set = analysis_set,
    endpoint = endpoint,
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
  out[term == "immune_reactivation_score_z"]
}

safe_cox <- function(dt, rhs, dataset, analysis_set, endpoint, model_name, min_events = 25) {
  needed <- all.vars(as.formula(paste("~", rhs)))
  model_vars <- unique(c("endpoint_time", "endpoint_event", needed))
  model_dt <- dt[complete.cases(dt[, ..model_vars])]
  event_count <- sum(model_dt$endpoint_event == 1L, na.rm = TRUE)
  if (nrow(model_dt) < 50 || event_count < min_events) {
    return(data.table(
      dataset = dataset,
      analysis_set = analysis_set,
      endpoint = endpoint,
      model = model_name,
      term = "MODEL_NOT_FITTED_LOW_EVENTS",
      n = nrow(model_dt),
      events = event_count,
      hr = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      p_value = NA_real_,
      concordance = NA_real_
    ))
  }
  fit <- coxph(as.formula(paste("Surv(endpoint_time, endpoint_event) ~", rhs)), data = model_dt)
  extract_score_term(model_name, fit, dataset, analysis_set, endpoint)
}

gse58812_expr <- fread(file.path(processed_dir, "GSE58812_GPL570_gene_symbol_expression.tsv.gz"), check.names = FALSE)
gse58812_prepared <- prepare_expression(
  expr_dt = gse58812_expr,
  gene_col = "SYMBOL",
  sample_cols = setdiff(names(gse58812_expr), "SYMBOL"),
  dataset = "GSE58812",
  transform = "none"
)
rm(gse58812_expr)
gse58812_deconv <- score_deconvolution(gse58812_prepared, "GSE58812")
rm(gse58812_prepared)

gse96058_expr <- fread(file.path(raw_dir, "GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz"), check.names = FALSE)
setnames(gse96058_expr, 1, "gene_symbol")
gse96058_sample_set <- fread(file.path(table_dir, "GSE96058_validation_sample_set.csv"))
gse96058_sample_cols <- intersect(gse96058_sample_set$sample_title, setdiff(names(gse96058_expr), "gene_symbol"))
gse96058_prepared <- prepare_expression(
  expr_dt = gse96058_expr,
  gene_col = "gene_symbol",
  sample_cols = gse96058_sample_cols,
  dataset = "GSE96058",
  transform = "none"
)
rm(gse96058_expr)
gse96058_deconv <- score_deconvolution(gse96058_prepared, "GSE96058")
rm(gse96058_prepared)

gse176078_expr <- fread(file.path(processed_dir, "GSE176078_bulkRNAseq_raw_counts.tsv.gz"), check.names = FALSE)
gse176078_prepared <- prepare_expression(
  expr_dt = gse176078_expr,
  gene_col = "gene_symbol",
  sample_cols = setdiff(names(gse176078_expr), "gene_symbol"),
  dataset = "GSE176078",
  transform = "log_cpm"
)
rm(gse176078_expr)
gse176078_deconv <- score_deconvolution(gse176078_prepared, "GSE176078")
rm(gse176078_prepared)

deconv_scores <- rbindlist(
  list(gse176078_deconv$scores, gse58812_deconv$scores, gse96058_deconv$scores),
  use.names = TRUE,
  fill = TRUE
)
deconv_coverage <- rbindlist(
  list(gse176078_deconv$coverage, gse58812_deconv$coverage, gse96058_deconv$coverage),
  use.names = TRUE,
  fill = TRUE
)

gse58812_rows <- fread(file.path(table_dir, "GSE58812_survival_analysis_rows.tsv.gz"))
gse58812_age <- fread(file.path(table_dir, "GSE58812_TNBC_survival_metadata.csv"))[
  ,
  .(geo_accession, age_at_diag_num = suppressWarnings(as.numeric(age_at_diag)))
]
gse58812_deconv_dt <- gse58812_deconv$scores[dataset == "GSE58812"]
gse58812_dt <- merge(gse58812_rows, gse58812_deconv_dt, by.x = "geo_accession", by.y = "sample_id", all.x = TRUE)
gse58812_dt <- merge(gse58812_dt, gse58812_age, by = "geo_accession", all.x = TRUE)
gse58812_dt[, `:=`(
  endpoint_time = as.numeric(time_days),
  endpoint_event = as.integer(event)
)]

gse96058_analysis <- fread(file.path(table_dir, "GSE96058_signature_survival_analysis_dataset.csv"))
gse96058_deconv_dt <- gse96058_deconv$scores[dataset == "GSE96058"]
setnames(gse96058_deconv_dt, "sample_id", "sample_title")
gse96058_dt <- merge(gse96058_analysis, gse96058_deconv_dt, by = "sample_title", all.x = TRUE)
gse96058_dt <- gse96058_dt[
  os_complete == TRUE &
    !is.na(os_days) &
    os_days > 0 &
    os_event %in% c(0L, 1L) &
    !is.na(immune_reactivation_score_z)
]
gse96058_dt[, `:=`(
  endpoint_time = as.numeric(os_days),
  endpoint_event = as.integer(os_event),
  age_num = suppressWarnings(as.numeric(age_at_diagnosis))
)]

cox_rows <- list()
for (endpoint_name in unique(gse58812_dt$endpoint)) {
  endpoint_dt <- gse58812_dt[endpoint == endpoint_name]
  cox_rows[[length(cox_rows) + 1]] <- safe_cox(
    endpoint_dt,
    "immune_reactivation_score_z + age_at_diag_num + purity_surrogate_z",
    "GSE58812",
    "tnbc",
    endpoint_name,
    "age_plus_purity_surrogate"
  )
  cox_rows[[length(cox_rows) + 1]] <- safe_cox(
    endpoint_dt,
    "immune_reactivation_score_z + age_at_diag_num + deconv_immune_summary_z",
    "GSE58812",
    "tnbc",
    endpoint_name,
    "age_plus_deconv_immune_summary"
  )
  cox_rows[[length(cox_rows) + 1]] <- safe_cox(
    endpoint_dt,
    "immune_reactivation_score_z + age_at_diag_num + deconv_stromal_summary_z + purity_surrogate_z",
    "GSE58812",
    "tnbc",
    endpoint_name,
    "age_plus_stromal_plus_purity"
  )
  cox_rows[[length(cox_rows) + 1]] <- safe_cox(
    endpoint_dt,
    "immune_reactivation_score_z + age_at_diag_num + deconv_immune_summary_z + deconv_stromal_summary_z + purity_surrogate_z",
    "GSE58812",
    "tnbc",
    endpoint_name,
    "age_plus_immune_stromal_purity"
  )
}

gse96058_sets <- list(
  all_primary = gse96058_dt,
  pam50_basal = gse96058_dt[pam50_subtype == "Basal"],
  pathology_tnbc = gse96058_dt[tnbc_pathology == TRUE],
  pathology_tnbc_or_pam50_basal = gse96058_dt[tnbc_or_basal == TRUE]
)
for (set_name in names(gse96058_sets)) {
  set_dt <- gse96058_sets[[set_name]]
  cox_rows[[length(cox_rows) + 1]] <- safe_cox(
    set_dt,
    "immune_reactivation_score_z + age_num + purity_surrogate_z",
    "GSE96058",
    set_name,
    "overall_survival",
    "age_plus_purity_surrogate"
  )
  cox_rows[[length(cox_rows) + 1]] <- safe_cox(
    set_dt,
    "immune_reactivation_score_z + age_num + deconv_immune_summary_z",
    "GSE96058",
    set_name,
    "overall_survival",
    "age_plus_deconv_immune_summary"
  )
  cox_rows[[length(cox_rows) + 1]] <- safe_cox(
    set_dt,
    "immune_reactivation_score_z + age_num + deconv_stromal_summary_z + purity_surrogate_z",
    "GSE96058",
    set_name,
    "overall_survival",
    "age_plus_stromal_plus_purity"
  )
  cox_rows[[length(cox_rows) + 1]] <- safe_cox(
    set_dt,
    "immune_reactivation_score_z + age_num + deconv_immune_summary_z + deconv_stromal_summary_z + purity_surrogate_z",
    "GSE96058",
    set_name,
    "overall_survival",
    "age_plus_immune_stromal_purity"
  )
}
cox_results <- rbindlist(cox_rows, fill = TRUE)
cox_results[, `:=`(
  p_value_formatted = format_p(p_value),
  hr_95_ci = fifelse(is.na(hr), NA_character_, sprintf("%.3f (%.3f-%.3f)", hr, ci_lower, ci_upper))
)]

cor_rows <- rbindlist(list(
  correlate_with_score(gse58812_dt[endpoint == "overall_survival"], "GSE58812", "tnbc_os_rows"),
  correlate_with_score(gse58812_dt[endpoint == "metastasis_free_survival"], "GSE58812", "tnbc_mfs_rows"),
  correlate_with_score(gse96058_dt, "GSE96058", "all_primary"),
  correlate_with_score(gse96058_dt[pam50_subtype == "Basal"], "GSE96058", "pam50_basal"),
  correlate_with_score(gse96058_dt[tnbc_pathology == TRUE], "GSE96058", "pathology_tnbc"),
  correlate_with_score(gse96058_dt[tnbc_or_basal == TRUE], "GSE96058", "pathology_tnbc_or_pam50_basal")
), fill = TRUE)

summary_vars <- c(
  "deconv_immune_summary_z",
  "deconv_stromal_summary_z",
  "deconv_epithelial_summary_z",
  "purity_surrogate_z",
  "immune_to_purity_axis_z"
)
summary_cor <- cor_rows[deconvolution_variable %in% summary_vars]
summary_cor[, label := sprintf("%.2f", spearman_rho)]

heatmap_plot <- ggplot(
  summary_cor,
  aes(x = deconvolution_variable, y = paste(dataset, analysis_set, sep = ": "), fill = spearman_rho)
) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = label), size = 2.6) +
  scale_fill_gradient2(low = "#3b6fb6", mid = "white", high = "#c43b3b", midpoint = 0, limits = c(-1, 1)) +
  labs(x = NULL, y = NULL, fill = "Spearman rho") +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "right"
  )

ggsave(
  filename = file.path(figure_dir, "deconvolution_score_correlation_heatmap.png"),
  plot = heatmap_plot,
  width = 8,
  height = 4.5,
  dpi = 300
)

write.csv(deconv_scores, file.path(table_dir, "signature_deconvolution_scores.csv"), row.names = FALSE)
write.csv(deconv_coverage, file.path(table_dir, "signature_deconvolution_gene_coverage.csv"), row.names = FALSE)
write.csv(cor_rows, file.path(table_dir, "signature_deconvolution_score_correlations.csv"), row.names = FALSE)
write.csv(cox_results, file.path(table_dir, "signature_deconvolution_adjusted_cox.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_15_signature_deconvolution_purity_sensitivity.txt"))

cat("Signature-based deconvolution and purity sensitivity complete\n")
cat("\nCoverage:\n")
print(deconv_coverage)
cat("\nAdjusted Cox rows:\n")
print(cox_results)
cat("\nScore-deconvolution correlations:\n")
print(cor_rows)
