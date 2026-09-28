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
})

raw_dir <- file.path(project_root, "data", "raw")
processed_dir <- file.path(project_root, "data", "processed")
table_dir <- file.path(project_root, "results", "tables")
log_dir <- file.path(project_root, "logs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

marker_sets <- list(
  leukocyte_proxy = c("PTPRC", "CD3D", "CD3E", "CD4", "MS4A1", "CD79A", "CD68", "LST1", "ITGAM", "FCGR3A", "MZB1", "JCHAIN"),
  stromal_proxy = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FAP", "ACTA2", "TAGLN", "VIM", "PDGFRB", "PECAM1", "VWF", "ENG"),
  epithelial_proxy = c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT5", "KRT14", "KRT17", "MUC1", "CDH1", "KRT7")
)
marker_genes <- unique(toupper(unlist(marker_sets, use.names = FALSE)))

format_p <- function(x) {
  out <- rep(NA_character_, length(x))
  out[!is.na(x) & x < 0.001] <- "<0.001"
  out[!is.na(x) & x >= 0.001] <- sprintf("%.3f", x[!is.na(x) & x >= 0.001])
  out
}

zscore_vector <- function(x) {
  out <- as.numeric(scale(x))
  out[is.na(out)] <- 0
  out
}

score_marker_sets <- function(expr_dt, gene_col, sample_cols, dataset) {
  sample_cols <- setdiff(sample_cols, gene_col)
  expr_dt <- copy(expr_dt)
  setnames(expr_dt, gene_col, "gene_symbol")
  expr_dt[, gene_symbol := toupper(gene_symbol)]
  expr_dt <- expr_dt[gene_symbol %in% marker_genes]
  if (nrow(expr_dt) == 0) {
    stop("No marker genes found for ", dataset)
  }
  for (col in sample_cols) {
    set(expr_dt, j = col, value = suppressWarnings(as.numeric(expr_dt[[col]])))
  }
  gene_expr <- expr_dt[
    ,
    lapply(.SD, mean, na.rm = TRUE),
    by = gene_symbol,
    .SDcols = sample_cols
  ]
  mat <- as.matrix(gene_expr[, ..sample_cols])
  mode(mat) <- "numeric"
  rownames(mat) <- gene_expr$gene_symbol
  z_mat <- t(scale(t(mat)))
  z_mat[is.na(z_mat)] <- 0

  score_dt <- data.table(sample_id = colnames(z_mat))
  coverage_rows <- vector("list", length(marker_sets))
  names(coverage_rows) <- names(marker_sets)
  for (set_name in names(marker_sets)) {
    genes <- unique(toupper(marker_sets[[set_name]]))
    present <- intersect(genes, rownames(z_mat))
    missing <- setdiff(genes, rownames(z_mat))
    coverage_rows[[set_name]] <- data.table(
      dataset = dataset,
      proxy = set_name,
      requested_genes = length(genes),
      present_genes = length(present),
      missing_genes = paste(missing, collapse = ";"),
      present_gene_list = paste(present, collapse = ";")
    )
    if (length(present) == 0) {
      score_dt[, (set_name) := NA_real_]
      score_dt[, (paste0(set_name, "_z")) := NA_real_]
    } else {
      score_dt[, (set_name) := colMeans(z_mat[present, , drop = FALSE], na.rm = TRUE)]
      score_dt[, (paste0(set_name, "_z")) := zscore_vector(get(set_name))]
    }
  }
  score_dt[, dataset := dataset]
  setcolorder(score_dt, c("dataset", "sample_id", setdiff(names(score_dt), c("dataset", "sample_id"))))
  list(scores = score_dt, coverage = rbindlist(coverage_rows, use.names = TRUE, fill = TRUE))
}

extract_term <- function(model_name, model, dataset, analysis_set, endpoint, term_pattern = "^immune_reactivation_score_z$") {
  model_summary <- summary(model)
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
  out[grepl(term_pattern, term)]
}

safe_cox <- function(dt, formula_rhs, dataset, analysis_set, endpoint, model_name, min_events = 25) {
  needed <- all.vars(as.formula(paste("~", formula_rhs)))
  proxy_vars <- grep("_proxy_z$", needed, value = TRUE)
  unavailable_proxy <- proxy_vars[
    !proxy_vars %in% names(dt) |
      vapply(proxy_vars, function(var) {
        if (!var %in% names(dt)) {
          return(TRUE)
        }
        all(is.na(dt[[var]]))
      }, logical(1))
  ]
  if (length(unavailable_proxy) > 0) {
    return(data.table(
      dataset = dataset,
      analysis_set = analysis_set,
      endpoint = endpoint,
      model = model_name,
      term = "MODEL_NOT_FITTED_NO_MARKERS",
      n = sum(complete.cases(dt[, .(endpoint_time, endpoint_event)])),
      events = sum(dt$endpoint_event == 1L, na.rm = TRUE),
      hr = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      p_value = NA_real_,
      concordance = NA_real_
    ))
  }
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
  fit <- coxph(as.formula(paste("Surv(endpoint_time, endpoint_event) ~", formula_rhs)), data = model_dt)
  extract_term(model_name, fit, dataset, analysis_set, endpoint)
}

make_correlations <- function(dt, dataset, analysis_set) {
  proxy_cols <- paste0(names(marker_sets), "_z")
  rbindlist(lapply(proxy_cols, function(proxy_col) {
    proxy_name <- sub("_z$", "", proxy_col)
    if (!proxy_col %in% names(dt) || all(is.na(dt[[proxy_col]]))) {
      return(data.table(
        dataset = dataset,
        analysis_set = analysis_set,
        proxy = proxy_name,
        n = sum(!is.na(dt$immune_reactivation_score_z)),
        spearman_rho = NA_real_,
        p_value = NA_real_,
        p_value_formatted = NA_character_,
        status = "NO_MARKERS_AVAILABLE"
      ))
    }
    keep <- complete.cases(dt[, .(immune_reactivation_score_z, get(proxy_col))])
    if (sum(keep) < 3) {
      return(data.table(
        dataset = dataset,
        analysis_set = analysis_set,
        proxy = proxy_name,
        n = sum(keep),
        spearman_rho = NA_real_,
        p_value = NA_real_,
        p_value_formatted = NA_character_,
        status = "NOT_FITTED_TOO_FEW_COMPLETE_CASES"
      ))
    }
    test <- suppressWarnings(cor.test(dt$immune_reactivation_score_z[keep], dt[[proxy_col]][keep], method = "spearman"))
    data.table(
      dataset = dataset,
      analysis_set = analysis_set,
      proxy = proxy_name,
      n = sum(keep),
      spearman_rho = unname(test$estimate),
      p_value = test$p.value,
      p_value_formatted = format_p(test$p.value),
      status = "OK"
    )
  }), fill = TRUE)
}

gse58812_expr_path <- file.path(processed_dir, "GSE58812_GPL570_gene_symbol_expression.tsv.gz")
gse96058_expr_path <- file.path(raw_dir, "GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz")
if (!file.exists(gse58812_expr_path)) {
  stop("Missing GSE58812 gene-symbol expression matrix: ", gse58812_expr_path)
}
if (!file.exists(gse96058_expr_path)) {
  stop("Missing GSE96058 expression matrix: ", gse96058_expr_path)
}

gse58812_expr <- fread(gse58812_expr_path, check.names = FALSE)
gse58812_scores <- score_marker_sets(
  expr_dt = gse58812_expr,
  gene_col = "SYMBOL",
  sample_cols = setdiff(names(gse58812_expr), "SYMBOL"),
  dataset = "GSE58812"
)
rm(gse58812_expr)

gse96058_expr <- fread(gse96058_expr_path, check.names = FALSE)
setnames(gse96058_expr, 1, "gene_symbol")
gse96058_sample_set <- fread(file.path(table_dir, "GSE96058_validation_sample_set.csv"))
gse96058_sample_cols <- intersect(gse96058_sample_set$sample_title, setdiff(names(gse96058_expr), "gene_symbol"))
gse96058_scores <- score_marker_sets(
  expr_dt = gse96058_expr,
  gene_col = "gene_symbol",
  sample_cols = gse96058_sample_cols,
  dataset = "GSE96058"
)
rm(gse96058_expr)

proxy_scores <- rbindlist(list(gse58812_scores$scores, gse96058_scores$scores), fill = TRUE)
proxy_coverage <- rbindlist(list(gse58812_scores$coverage, gse96058_scores$coverage), fill = TRUE)

gse58812_endpoint_rows <- fread(file.path(table_dir, "GSE58812_survival_analysis_rows.tsv.gz"))
gse58812_surv_meta <- fread(file.path(table_dir, "GSE58812_TNBC_survival_metadata.csv"))
gse58812_age <- gse58812_surv_meta[, .(geo_accession, age_at_diag_num = suppressWarnings(as.numeric(age_at_diag)))]
gse58812_proxy <- gse58812_scores$scores[, .(sample_id, leukocyte_proxy_z, stromal_proxy_z, epithelial_proxy_z)]
gse58812_dt <- merge(gse58812_endpoint_rows, gse58812_proxy, by.x = "geo_accession", by.y = "sample_id", all.x = TRUE)
gse58812_dt <- merge(gse58812_dt, gse58812_age, by = "geo_accession", all.x = TRUE)
gse58812_dt[, `:=`(
  endpoint_time = as.numeric(time_days),
  endpoint_event = as.integer(event)
)]

gse96058_analysis <- fread(file.path(table_dir, "GSE96058_signature_survival_analysis_dataset.csv"))
gse96058_proxy <- gse96058_scores$scores[
  ,
  .(
    sample_id,
    leukocyte_proxy_z,
    stromal_proxy_z,
    epithelial_proxy_z
  )
]
setnames(gse96058_proxy, "sample_id", "sample_title")
gse96058_dt <- merge(gse96058_analysis, gse96058_proxy, by = "sample_title", all.x = TRUE)
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
  age_num = suppressWarnings(as.numeric(age_at_diagnosis)),
  tumor_size_num = suppressWarnings(as.numeric(tumor_size)),
  lymph_node_factor = factor(lymph_node_group),
  chemo_treated_factor = factor(chemo_treated)
)]

cox_rows <- list()
for (endpoint_name in unique(gse58812_dt$endpoint)) {
  endpoint_dt <- gse58812_dt[endpoint == endpoint_name]
  cox_rows[[length(cox_rows) + 1]] <- safe_cox(
    endpoint_dt,
    "immune_reactivation_score_z + age_at_diag_num + leukocyte_proxy_z",
    "GSE58812",
    "tnbc",
    endpoint_name,
    "age_plus_leukocyte_proxy"
  )
  cox_rows[[length(cox_rows) + 1]] <- safe_cox(
    endpoint_dt,
    "immune_reactivation_score_z + age_at_diag_num + stromal_proxy_z",
    "GSE58812",
    "tnbc",
    endpoint_name,
    "age_plus_stromal_proxy"
  )
  cox_rows[[length(cox_rows) + 1]] <- safe_cox(
    endpoint_dt,
    "immune_reactivation_score_z + age_at_diag_num + epithelial_proxy_z",
    "GSE58812",
    "tnbc",
    endpoint_name,
    "age_plus_epithelial_proxy"
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
    "immune_reactivation_score_z + age_num + leukocyte_proxy_z",
    "GSE96058",
    set_name,
    "overall_survival",
    "age_plus_leukocyte_proxy"
  )
  cox_rows[[length(cox_rows) + 1]] <- safe_cox(
    set_dt,
    "immune_reactivation_score_z + age_num + stromal_proxy_z",
    "GSE96058",
    set_name,
    "overall_survival",
    "age_plus_stromal_proxy"
  )
  cox_rows[[length(cox_rows) + 1]] <- safe_cox(
    set_dt,
    "immune_reactivation_score_z + age_num + epithelial_proxy_z",
    "GSE96058",
    set_name,
    "overall_survival",
    "age_plus_epithelial_proxy"
  )
}
cox_results <- rbindlist(cox_rows, fill = TRUE)
cox_results[, `:=`(
  p_value_formatted = format_p(p_value),
  hr_95_ci = fifelse(
    is.na(hr),
    NA_character_,
    sprintf("%.3f (%.3f-%.3f)", hr, ci_lower, ci_upper)
  )
)]

cor_rows <- rbindlist(list(
  make_correlations(gse58812_dt[endpoint == "overall_survival"], "GSE58812", "tnbc_os_rows"),
  make_correlations(gse58812_dt[endpoint == "metastasis_free_survival"], "GSE58812", "tnbc_mfs_rows"),
  make_correlations(gse96058_dt, "GSE96058", "all_primary"),
  make_correlations(gse96058_dt[pam50_subtype == "Basal"], "GSE96058", "pam50_basal"),
  make_correlations(gse96058_dt[tnbc_pathology == TRUE], "GSE96058", "pathology_tnbc"),
  make_correlations(gse96058_dt[tnbc_or_basal == TRUE], "GSE96058", "pathology_tnbc_or_pam50_basal")
), fill = TRUE)

write.csv(proxy_scores, file.path(table_dir, "microenvironment_proxy_scores.csv"), row.names = FALSE)
write.csv(proxy_coverage, file.path(table_dir, "microenvironment_proxy_gene_coverage.csv"), row.names = FALSE)
write.csv(cox_results, file.path(table_dir, "microenvironment_adjusted_cox.csv"), row.names = FALSE)
write.csv(cor_rows, file.path(table_dir, "microenvironment_proxy_correlations.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_13_microenvironment_adjustment_sensitivity.txt"))

cat("Microenvironment adjustment sensitivity complete\n")
cat("\nMarker coverage:\n")
print(proxy_coverage)
cat("\nProxy-adjusted Cox rows:\n")
print(cox_results)
cat("\nScore-proxy Spearman correlations:\n")
print(cor_rows)
