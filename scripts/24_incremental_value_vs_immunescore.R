#!/usr/bin/env Rscript

# Round-9 review-driven analysis:
# Symmetric / incremental-value test of the immune reactivation score versus
# ESTIMATE ImmuneScore. Earlier scripts only showed that the score loses
# significance after adjusting for ImmuneScore. This script adds the reverse
# direction and model-fit comparison so that "redundant" vs "adds value" is
# tested, not asserted:
#   - M_age   : Surv ~ age
#   - M_score : Surv ~ score_z + age
#   - M_immune: Surv ~ ImmuneScore_z + age
#   - M_both  : Surv ~ score_z + ImmuneScore_z + age
# Likelihood-ratio tests:
#   score added to immune  = anova(M_immune, M_both)
#   immune added to score  = anova(M_score , M_both)
# plus Harrell C-index and AIC for each model on identical complete-case rows.

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

table_dir <- file.path(project_root, "results", "tables")
log_dir <- file.path(project_root, "logs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
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

# ---- Load per-sample scores and ESTIMATE scores (same sources as script 18) ----
estimate_scores <- fread(file.path(table_dir, "estimate_scores.csv"))
gse58812_estimate <- estimate_scores[dataset == "GSE58812"]
gse96058_estimate <- estimate_scores[dataset == "GSE96058"]

# ---- GSE58812 merge (key: geo_accession <-> sample_id) ----
gse58812_scores <- fread(file.path(table_dir, "GSE58812_immune_reactivation_scores.csv"))
if (!"immune_reactivation_score_z" %in% names(gse58812_scores)) {
  gse58812_scores[, immune_reactivation_score_z := zscore(immune_reactivation_score)]
}
gse58812_surv <- fread(file.path(table_dir, "GSE58812_TNBC_survival_metadata.csv"))
gse58812_dt <- merge(gse58812_surv, gse58812_estimate,
                     by.x = "geo_accession", by.y = "sample_id", all.x = TRUE)
gse58812_dt <- merge(
  gse58812_dt,
  gse58812_scores[, .(sample_id, immune_reactivation_score_z)],
  by.x = "geo_accession", by.y = "sample_id", all.x = TRUE
)
gse58812_dt[, age_num := suppressWarnings(as.numeric(age_at_diag))]

# ---- GSE96058 merge (key: primary_sample_id <-> sample_id) ----
gse96058_dt <- fread(file.path(table_dir, "GSE96058_signature_survival_analysis_dataset.csv"))
gse96058_dt <- merge(gse96058_dt, gse96058_estimate,
                     by.x = "primary_sample_id", by.y = "sample_id", all.x = TRUE)
gse96058_dt[, age_num := suppressWarnings(as.numeric(age_at_diagnosis))]

# ---- Core comparison routine ----
# Fits the four nested age-adjusted models on identical complete-case rows and
# returns one summary row plus the focal-term HR rows.
compare_models <- function(dt, dataset, analysis_set, endpoint,
                           time_col, event_col, age_col, min_events = 20) {
  d <- copy(dt)
  d[, endpoint_time := get(time_col)]
  d[, endpoint_event := as.integer(get(event_col))]
  d[, age_adjustment := get(age_col)]
  vars <- c("endpoint_time", "endpoint_event", "age_adjustment",
            "immune_reactivation_score_z", "ImmuneScore_z")
  d <- d[complete.cases(d[, ..vars]) & endpoint_time > 0]
  n <- nrow(d)
  events <- sum(d$endpoint_event == 1L, na.rm = TRUE)

  empty_summary <- data.table(
    dataset = dataset, analysis_set = analysis_set, endpoint = endpoint,
    n = n, events = events,
    cindex_age = NA_real_, cindex_score = NA_real_,
    cindex_immune = NA_real_, cindex_both = NA_real_,
    aic_age = NA_real_, aic_score = NA_real_,
    aic_immune = NA_real_, aic_both = NA_real_,
    lrt_score_added_to_immune_p = NA_real_,
    lrt_immune_added_to_score_p = NA_real_,
    hr_score_in_both = NA_real_, hr_immune_in_both = NA_real_,
    p_score_in_both = NA_real_, p_immune_in_both = NA_real_
  )
  if (n < 40 || events < min_events) {
    empty_summary[, note := "MODEL_NOT_FITTED_LOW_EVENTS"]
    return(list(summary = empty_summary, terms = data.table()))
  }

  m_age <- coxph(Surv(endpoint_time, endpoint_event) ~ age_adjustment, data = d)
  m_score <- coxph(Surv(endpoint_time, endpoint_event) ~ immune_reactivation_score_z + age_adjustment, data = d)
  m_immune <- coxph(Surv(endpoint_time, endpoint_event) ~ ImmuneScore_z + age_adjustment, data = d)
  m_both <- coxph(Surv(endpoint_time, endpoint_event) ~ immune_reactivation_score_z + ImmuneScore_z + age_adjustment, data = d)

  ci <- function(m) unname(summary(m)$concordance[1])
  lrt_p <- function(reduced, full) {
    a <- anova(reduced, full)
    pcol <- grep("^P", names(a), value = TRUE)[1]
    a[[pcol]][2]
  }
  both_s <- summary(m_both)
  both_conf <- both_s$conf.int
  both_coef <- both_s$coefficients

  summary_row <- data.table(
    dataset = dataset, analysis_set = analysis_set, endpoint = endpoint,
    n = n, events = events,
    cindex_age = ci(m_age), cindex_score = ci(m_score),
    cindex_immune = ci(m_immune), cindex_both = ci(m_both),
    aic_age = AIC(m_age), aic_score = AIC(m_score),
    aic_immune = AIC(m_immune), aic_both = AIC(m_both),
    lrt_score_added_to_immune_p = lrt_p(m_immune, m_both),
    lrt_immune_added_to_score_p = lrt_p(m_score, m_both),
    hr_score_in_both = unname(both_conf["immune_reactivation_score_z", "exp(coef)"]),
    hr_immune_in_both = unname(both_conf["ImmuneScore_z", "exp(coef)"]),
    p_score_in_both = unname(both_coef["immune_reactivation_score_z", "Pr(>|z|)"]),
    p_immune_in_both = unname(both_coef["ImmuneScore_z", "Pr(>|z|)"]),
    note = "fitted"
  )

  term_rows <- rbindlist(lapply(
    list(score_only = m_score, immune_only = m_immune, both = m_both),
    function(m) {
      s <- summary(m); conf <- s$conf.int; coef <- s$coefficients
      data.table(
        term = rownames(coef),
        hr = unname(conf[, "exp(coef)"]),
        ci_lower = unname(conf[, "lower .95"]),
        ci_upper = unname(conf[, "upper .95"]),
        p_value = unname(coef[, "Pr(>|z|)"]),
        cindex = ci(m), aic = AIC(m)
      )
    }
  ), idcol = "model")
  term_rows[, `:=`(dataset = dataset, analysis_set = analysis_set, endpoint = endpoint)]
  list(summary = summary_row, terms = term_rows[term %in% c("immune_reactivation_score_z", "ImmuneScore_z")])
}

specs <- list(
  list(gse58812_dt, "GSE58812", "tnbc", "overall_survival", "os_days", "death", "age_num"),
  list(gse58812_dt, "GSE58812", "tnbc", "metastasis_free_survival", "mfs_days", "meta", "age_num"),
  list(gse96058_dt, "GSE96058", "all_primary", "overall_survival", "os_days", "os_event", "age_num"),
  list(gse96058_dt[basal_pam50 == TRUE], "GSE96058", "pam50_basal", "overall_survival", "os_days", "os_event", "age_num"),
  list(gse96058_dt[tnbc_pathology == TRUE], "GSE96058", "pathology_tnbc", "overall_survival", "os_days", "os_event", "age_num")
)

results <- lapply(specs, function(s) do.call(compare_models, s))
summary_dt <- rbindlist(lapply(results, `[[`, "summary"), fill = TRUE)
terms_dt <- rbindlist(lapply(results, `[[`, "terms"), fill = TRUE)

# Readable formatting
summary_dt[, `:=`(
  delta_cindex_both_vs_immune = cindex_both - cindex_immune,
  delta_cindex_both_vs_score = cindex_both - cindex_score,
  delta_aic_both_vs_immune = aic_both - aic_immune,
  lrt_score_added_to_immune_p_fmt = format_p(lrt_score_added_to_immune_p),
  lrt_immune_added_to_score_p_fmt = format_p(lrt_immune_added_to_score_p)
)]
if (nrow(terms_dt) > 0) {
  terms_dt[, `:=`(
    hr_ci = sprintf("%.3f (%.3f-%.3f)", hr, ci_lower, ci_upper),
    p_value_fmt = format_p(p_value)
  )]
}

fwrite(summary_dt, file.path(table_dir, "incremental_value_model_comparison.csv"))
fwrite(terms_dt, file.path(table_dir, "incremental_value_cox_terms.csv"))

cat("\n==== Incremental value: model comparison ====\n")
print(summary_dt[, .(dataset, analysis_set, endpoint, n, events,
                     cindex_score, cindex_immune, cindex_both,
                     lrt_score_added_to_immune_p_fmt,
                     lrt_immune_added_to_score_p_fmt,
                     hr_score_in_both, p_score_in_both,
                     hr_immune_in_both, p_immune_in_both)])

writeLines(capture.output(sessionInfo()),
           file.path(log_dir, "sessionInfo_24_incremental_value_vs_immunescore.txt"))
cat("\nIncremental value analysis complete\n")
