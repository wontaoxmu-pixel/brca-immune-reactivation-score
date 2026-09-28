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
  library(estimate)
})

raw_dir <- file.path(project_root, "data", "raw")
processed_dir <- file.path(project_root, "data", "processed")
table_dir <- file.path(project_root, "results", "tables")
estimate_dir <- file.path(processed_dir, "estimate")
log_dir <- file.path(project_root, "logs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(estimate_dir, recursive = TRUE, showWarnings = FALSE)
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

read_estimate_gct <- function(path, dataset) {
  score_dt <- fread(path, skip = 2)
  setnames(score_dt, 1, "metric")
  score_dt[, Description := NULL]
  long <- melt(
    score_dt,
    id.vars = "metric",
    variable.name = "sample_id",
    value.name = "score"
  )
  wide <- dcast(long, sample_id ~ metric, value.var = "score")
  wide[, dataset := dataset]
  if (!"TumorPurity" %in% names(wide) && "ESTIMATEScore" %in% names(wide)) {
    wide[, TumorPurity := cos(0.6049872018 + 0.0001467884 * ESTIMATEScore)]
    wide[TumorPurity < 0, TumorPurity := NA_real_]
  }
  for (col in intersect(c("StromalScore", "ImmuneScore", "ESTIMATEScore", "TumorPurity"), names(wide))) {
    wide[, (paste0(col, "_z")) := zscore(get(col))]
  }
  setcolorder(wide, c("dataset", "sample_id", setdiff(names(wide), c("dataset", "sample_id"))))
  wide
}

write_expression_txt <- function(expr_dt, gene_col, sample_cols, out_txt, transform = c("none", "log2_plus_1")) {
  transform <- match.arg(transform)
  expr_dt <- copy(expr_dt)
  setnames(expr_dt, gene_col, "gene_symbol")
  expr_dt[, gene_symbol := toupper(gene_symbol)]
  expr_dt <- expr_dt[!is.na(gene_symbol) & nzchar(gene_symbol)]
  sample_cols <- intersect(sample_cols, names(expr_dt))
  for (col in sample_cols) {
    set(expr_dt, j = col, value = suppressWarnings(as.numeric(expr_dt[[col]])))
  }
  collapsed <- expr_dt[
    ,
    lapply(.SD, mean, na.rm = TRUE),
    by = gene_symbol,
    .SDcols = sample_cols
  ]
  if (transform == "log2_plus_1") {
    for (col in sample_cols) {
      vals <- collapsed[[col]]
      vals[is.na(vals)] <- 0
      min_val <- min(vals, na.rm = TRUE)
      if (is.finite(min_val) && min_val < 0) {
        vals <- vals - min_val
      }
      collapsed[[col]] <- log2(vals + 1)
    }
  }
  fwrite(collapsed, out_txt, sep = "\t")
  data.table(
    input_file = basename(out_txt),
    genes_written = nrow(collapsed),
    samples_written = length(sample_cols),
    transform = transform
  )
}

run_estimate <- function(expr_txt, dataset, platform) {
  filtered_gct <- file.path(estimate_dir, paste0(dataset, "_estimate_common_genes.gct"))
  score_gct <- file.path(estimate_dir, paste0(dataset, "_estimate_scores.gct"))
  if (!file.exists(score_gct)) {
    if (!file.exists(filtered_gct)) {
      filterCommonGenes(input.f = expr_txt, output.f = filtered_gct, id = "GeneSymbol")
    }
    estimateScore(filtered_gct, score_gct, platform = platform)
  }
  read_estimate_gct(score_gct, dataset)
}

summarize_existing_input <- function(path, transform) {
  if (!file.exists(path)) {
    return(data.table(
      input_file = basename(path),
      genes_written = NA_integer_,
      samples_written = NA_integer_,
      transform = transform
    ))
  }
  header <- names(fread(path, nrows = 0))
  wc_out <- system2("wc", c("-l", path), stdout = TRUE)
  line_count <- suppressWarnings(as.integer(strsplit(trimws(wc_out[1]), "\\s+")[[1]][1]))
  data.table(
    input_file = basename(path),
    genes_written = ifelse(!is.na(line_count), line_count - 1L, NA_integer_),
    samples_written = max(length(header) - 1L, 0L),
    transform = transform
  )
}

correlate_estimate <- function(dt, dataset, analysis_set) {
  vars <- intersect(
    c("StromalScore_z", "ImmuneScore_z", "ESTIMATEScore_z", "TumorPurity_z"),
    names(dt)
  )
  rbindlist(lapply(vars, function(var) {
    keep <- complete.cases(dt[, .(immune_reactivation_score_z, get(var))])
    if (sum(keep) < 3) {
      return(data.table(
        dataset = dataset,
        analysis_set = analysis_set,
        variable = var,
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
      variable = var,
      n = sum(keep),
      spearman_rho = unname(test$estimate),
      p_value = test$p.value,
      p_value_formatted = format_p(test$p.value)
    )
  }), fill = TRUE)
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

run_endpoint_models <- function(dt, dataset, analysis_set, endpoint, time_col, event_col, age_col) {
  model_dt <- copy(dt)
  model_dt[, endpoint_time := get(time_col)]
  model_dt[, endpoint_event := get(event_col)]
  model_dt[, age_adjustment := get(age_col)]
  models <- list(
    age_plus_estimate_immune = "immune_reactivation_score_z + age_adjustment + ImmuneScore_z",
    age_plus_estimate_stromal = "immune_reactivation_score_z + age_adjustment + StromalScore_z",
    age_plus_estimate_score = "immune_reactivation_score_z + age_adjustment + ESTIMATEScore_z",
    age_plus_estimate_purity = "immune_reactivation_score_z + age_adjustment + TumorPurity_z",
    age_plus_estimate_immune_stromal_purity = "immune_reactivation_score_z + age_adjustment + ImmuneScore_z + StromalScore_z + TumorPurity_z"
  )
  rbindlist(lapply(names(models), function(model_name) {
    safe_cox(model_dt, models[[model_name]], dataset, analysis_set, endpoint, model_name)
  }), fill = TRUE)
}

cat("Preparing GSE58812 ESTIMATE input\n")
gse58812_score_gct <- file.path(estimate_dir, "GSE58812_estimate_scores.gct")
if (file.exists(gse58812_score_gct)) {
  gse58812_input_summary <- summarize_existing_input(
    file.path(estimate_dir, "GSE58812_estimate_input.txt"),
    transform = "none_cached_score_gct"
  )
  gse58812_estimate <- read_estimate_gct(gse58812_score_gct, "GSE58812")
} else {
  gse58812_expr <- fread(file.path(processed_dir, "GSE58812_GPL570_gene_symbol_expression.tsv.gz"))
  gse58812_sample_cols <- setdiff(names(gse58812_expr), "SYMBOL")
  gse58812_input <- file.path(estimate_dir, "GSE58812_estimate_input.txt")
  gse58812_input_summary <- write_expression_txt(gse58812_expr, "SYMBOL", gse58812_sample_cols, gse58812_input, transform = "none")
  gse58812_estimate <- run_estimate(gse58812_input, "GSE58812", platform = "affymetrix")
}

cat("Preparing GSE96058 ESTIMATE input\n")
gse96058_score_gct <- file.path(estimate_dir, "GSE96058_estimate_scores.gct")
if (file.exists(gse96058_score_gct)) {
  gse96058_input_summary <- summarize_existing_input(
    file.path(estimate_dir, "GSE96058_estimate_input.txt"),
    transform = "none_cached_score_gct"
  )
  gse96058_estimate <- read_estimate_gct(gse96058_score_gct, "GSE96058")
} else {
  gse96058_samples <- fread(file.path(table_dir, "GSE96058_validation_sample_set.csv"))
  gse96058_expr <- fread(file.path(raw_dir, "GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz"))
  setnames(gse96058_expr, 1, "gene_symbol")
  gse96058_sample_cols <- intersect(gse96058_samples$primary_sample_id, names(gse96058_expr))
  gse96058_input <- file.path(estimate_dir, "GSE96058_estimate_input.txt")
  gse96058_input_summary <- write_expression_txt(gse96058_expr, "gene_symbol", gse96058_sample_cols, gse96058_input, transform = "none")
  gse96058_estimate <- run_estimate(gse96058_input, "GSE96058", platform = "illumina")
}

input_summary <- rbindlist(list(
  cbind(dataset = "GSE58812", gse58812_input_summary),
  cbind(dataset = "GSE96058", gse96058_input_summary)
), fill = TRUE)

fwrite(input_summary, file.path(table_dir, "estimate_input_summary.csv"))
estimate_scores <- rbindlist(list(gse58812_estimate, gse96058_estimate), fill = TRUE)
fwrite(estimate_scores, file.path(table_dir, "estimate_scores.csv"))

gse58812_scores <- fread(file.path(table_dir, "GSE58812_immune_reactivation_scores.csv"))
if (!"immune_reactivation_score_z" %in% names(gse58812_scores)) {
  gse58812_scores[, immune_reactivation_score_z := zscore(immune_reactivation_score)]
}
gse58812_surv <- fread(file.path(table_dir, "GSE58812_TNBC_survival_metadata.csv"))
gse58812_dt <- merge(gse58812_surv, gse58812_estimate, by.x = "geo_accession", by.y = "sample_id", all.x = TRUE)
gse58812_dt[, age_at_diag_num := suppressWarnings(as.numeric(age_at_diag))]
gse58812_dt <- merge(
  gse58812_dt,
  gse58812_scores[, .(sample_id, immune_reactivation_score_z)],
  by.x = "geo_accession",
  by.y = "sample_id",
  all.x = TRUE,
  suffixes = c("", "_score")
)
if ("immune_reactivation_score_z_score" %in% names(gse58812_dt)) {
  gse58812_dt[is.na(immune_reactivation_score_z), immune_reactivation_score_z := immune_reactivation_score_z_score]
  gse58812_dt[, immune_reactivation_score_z_score := NULL]
}

gse96058_dt <- fread(file.path(table_dir, "GSE96058_signature_survival_analysis_dataset.csv"))
gse96058_dt <- merge(
  gse96058_dt,
  gse96058_estimate,
  by.x = "primary_sample_id",
  by.y = "sample_id",
  all.x = TRUE
)
gse96058_dt[, age_num := suppressWarnings(as.numeric(age_at_diagnosis))]

cor_rows <- rbindlist(list(
  correlate_estimate(gse58812_dt, "GSE58812", "tnbc"),
  correlate_estimate(gse96058_dt, "GSE96058", "all_primary"),
  correlate_estimate(gse96058_dt[basal_pam50 == TRUE], "GSE96058", "pam50_basal"),
  correlate_estimate(gse96058_dt[tnbc_pathology == TRUE], "GSE96058", "pathology_tnbc"),
  correlate_estimate(gse96058_dt[tnbc_or_basal == TRUE], "GSE96058", "pathology_tnbc_or_pam50_basal")
), fill = TRUE)
fwrite(cor_rows, file.path(table_dir, "estimate_score_correlations.csv"))

cox_rows <- rbindlist(list(
  run_endpoint_models(gse58812_dt, "GSE58812", "tnbc", "overall_survival", "os_days", "death", "age_at_diag_num"),
  run_endpoint_models(gse58812_dt, "GSE58812", "tnbc", "metastasis_free_survival", "mfs_days", "meta", "age_at_diag_num"),
  run_endpoint_models(gse96058_dt, "GSE96058", "all_primary", "overall_survival", "os_days", "os_event", "age_num"),
  run_endpoint_models(gse96058_dt[basal_pam50 == TRUE], "GSE96058", "pam50_basal", "overall_survival", "os_days", "os_event", "age_num"),
  run_endpoint_models(gse96058_dt[tnbc_pathology == TRUE], "GSE96058", "pathology_tnbc", "overall_survival", "os_days", "os_event", "age_num"),
  run_endpoint_models(gse96058_dt[tnbc_or_basal == TRUE], "GSE96058", "pathology_tnbc_or_pam50_basal", "overall_survival", "os_days", "os_event", "age_num")
), fill = TRUE)
cox_rows[, `:=`(
  p_value_formatted = format_p(p_value),
  hr_ci = fifelse(
    is.na(hr),
    NA_character_,
    sprintf("%.3f (%.3f-%.3f)", hr, ci_lower, ci_upper)
  )
)]
fwrite(cox_rows, file.path(table_dir, "estimate_adjusted_cox.csv"))

writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_18_estimate_purity_sensitivity.txt"))
cat("ESTIMATE purity sensitivity complete\n")
