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

table_dir <- file.path(project_root, "results", "tables")
log_dir <- file.path(project_root, "logs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

analysis_path <- file.path(table_dir, "GSE96058_signature_survival_analysis_dataset.csv")
cox_path <- file.path(table_dir, "bulk_validation_continuous_cox_summary.csv")
if (!file.exists(analysis_path)) {
  stop("Missing GSE96058 analysis dataset: ", analysis_path)
}
if (!file.exists(cox_path)) {
  stop("Missing continuous Cox summary: ", cox_path)
}

analysis_dt <- fread(analysis_path)
cox_dt <- fread(cox_path)

analysis_dt <- analysis_dt[
  os_complete == TRUE &
    !is.na(os_days) &
    os_days > 0 &
    os_event %in% c(0L, 1L) &
    !is.na(immune_reactivation_score_z)
]
analysis_dt[, `:=`(
  endpoint_time = as.numeric(os_days),
  endpoint_event = as.integer(os_event),
  basal_pam50_factor = factor(ifelse(basal_pam50 == TRUE, "Basal", "Non-basal"), levels = c("Non-basal", "Basal")),
  pam50_factor = factor(pam50_subtype),
  lymph_node_factor = factor(lymph_node_group),
  endocrine_treated_factor = factor(endocrine_treated),
  chemo_treated_factor = factor(chemo_treated),
  tumor_size_num = suppressWarnings(as.numeric(tumor_size)),
  age_num = suppressWarnings(as.numeric(age_at_diagnosis))
)]

extract_term <- function(model_name, model, term_pattern = NULL) {
  model_summary <- summary(model)
  conf <- model_summary$conf.int
  coefs <- model_summary$coefficients
  out <- data.table(
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
  if (!is.null(term_pattern)) {
    out <- out[grepl(term_pattern, term)]
  }
  out
}

format_p <- function(x) {
  out <- rep(NA_character_, length(x))
  out[!is.na(x) & x < 0.001] <- "<0.001"
  out[!is.na(x) & x >= 0.001] <- sprintf("%.3f", x[!is.na(x) & x >= 0.001])
  out
}

# FDR summaries for the already reported continuous-score Cox results.
score_rows <- cox_dt[
  term == "immune_reactivation_score_z" &
    grepl("continuous", model)
]
score_rows[, fdr_scope := fifelse(grepl("age_adjusted", model), "age_adjusted_continuous_models", "unadjusted_continuous_models")]
score_rows[, q_value := p.adjust(p_value, method = "BH"), by = fdr_scope]
score_rows[, `:=`(
  p_value_formatted = format_p(p_value),
  q_value_formatted = format_p(q_value)
)]
fdr_table <- score_rows[
  ,
  .(
    dataset,
    analysis_set,
    analysis_set_label,
    endpoint,
    endpoint_label,
    model,
    model_label,
    n,
    events,
    hr,
    ci_lower,
    ci_upper,
    p_value,
    q_value,
    p_value_formatted,
    q_value_formatted,
    fdr_scope
  )
]

# Interaction models in GSE96058 all-primary samples.
all_primary <- copy(analysis_dt)
interaction_models <- list(
  age_score_basal = coxph(Surv(endpoint_time, endpoint_event) ~ immune_reactivation_score_z * basal_pam50_factor + age_num, data = all_primary),
  age_score_pam50 = coxph(Surv(endpoint_time, endpoint_event) ~ immune_reactivation_score_z * pam50_factor + age_num, data = all_primary)
)
basal_no_interaction <- coxph(Surv(endpoint_time, endpoint_event) ~ immune_reactivation_score_z + basal_pam50_factor + age_num, data = all_primary)
pam50_no_interaction <- coxph(Surv(endpoint_time, endpoint_event) ~ immune_reactivation_score_z + pam50_factor + age_num, data = all_primary)

interaction_terms <- rbindlist(list(
  extract_term("age_score_basal", interaction_models$age_score_basal, ":"),
  extract_term("age_score_pam50", interaction_models$age_score_pam50, ":")
), fill = TRUE)

interaction_lrt <- rbindlist(list(
  data.table(
    comparison = "score_by_basal_binary",
    n = interaction_models$age_score_basal$n,
    events = interaction_models$age_score_basal$nevent,
    loglik_without_interaction = as.numeric(logLik(basal_no_interaction)),
    loglik_with_interaction = as.numeric(logLik(interaction_models$age_score_basal)),
    df = attr(logLik(interaction_models$age_score_basal), "df") - attr(logLik(basal_no_interaction), "df")
  ),
  data.table(
    comparison = "score_by_pam50_factor",
    n = interaction_models$age_score_pam50$n,
    events = interaction_models$age_score_pam50$nevent,
    loglik_without_interaction = as.numeric(logLik(pam50_no_interaction)),
    loglik_with_interaction = as.numeric(logLik(interaction_models$age_score_pam50)),
    df = attr(logLik(interaction_models$age_score_pam50), "df") - attr(logLik(pam50_no_interaction), "df")
  )
), fill = TRUE)
interaction_lrt[, `:=`(
  chisq = 2 * (loglik_with_interaction - loglik_without_interaction),
  p_value = pchisq(2 * (loglik_with_interaction - loglik_without_interaction), df = df, lower.tail = FALSE)
)]

# Multivariable sensitivity models. Keep small-event subgroups out of high-dimensional models.
run_multivariable <- function(dt, label, formula_rhs) {
  model_dt <- dt[complete.cases(dt[, all.vars(as.formula(paste("~", formula_rhs))), with = FALSE])]
  if (nrow(model_dt) < 50 || sum(model_dt$endpoint_event == 1L) < 30) {
    return(data.table(
      analysis_set = label,
      model = paste0("multivariable: ", formula_rhs),
      term = "MODEL_NOT_FITTED_LOW_EVENTS",
      n = nrow(model_dt),
      events = sum(model_dt$endpoint_event == 1L),
      hr = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      p_value = NA_real_,
      concordance = NA_real_
    ))
  }
  model <- coxph(as.formula(paste("Surv(endpoint_time, endpoint_event) ~", formula_rhs)), data = model_dt)
  extract_term(paste0("multivariable: ", formula_rhs), model, "^immune_reactivation_score_z$")[
    ,
    analysis_set := label
  ][
    ,
    setcolorder(.SD, c("analysis_set", setdiff(names(.SD), "analysis_set")))
  ]
}

multivariable_rows <- rbindlist(list(
  run_multivariable(
    all_primary,
    "all_primary",
    "immune_reactivation_score_z + age_num + pam50_factor + tumor_size_num + lymph_node_factor + endocrine_treated_factor + chemo_treated_factor"
  ),
  run_multivariable(
    all_primary[pam50_subtype == "Basal"],
    "pam50_basal",
    "immune_reactivation_score_z + age_num + tumor_size_num + lymph_node_factor + chemo_treated_factor"
  ),
  run_multivariable(
    all_primary[tnbc_pathology == TRUE],
    "pathology_tnbc",
    "immune_reactivation_score_z + age_num + tumor_size_num + lymph_node_factor + chemo_treated_factor"
  ),
  run_multivariable(
    all_primary[tnbc_or_basal == TRUE],
    "pathology_tnbc_or_pam50_basal",
    "immune_reactivation_score_z + age_num + tumor_size_num + lymph_node_factor + chemo_treated_factor"
  )
), fill = TRUE)

write.csv(fdr_table, file.path(table_dir, "bulk_validation_continuous_cox_fdr_sensitivity.csv"), row.names = FALSE)
write.csv(interaction_terms, file.path(table_dir, "GSE96058_score_subtype_interaction_terms.csv"), row.names = FALSE)
write.csv(interaction_lrt, file.path(table_dir, "GSE96058_score_subtype_interaction_lrt.csv"), row.names = FALSE)
write.csv(multivariable_rows, file.path(table_dir, "GSE96058_multivariable_sensitivity_cox.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_12_GSE96058_sensitivity_models.txt"))

cat("GSE96058 sensitivity models complete\n")
cat("\nFDR table:\n")
print(fdr_table)
cat("\nInteraction likelihood-ratio tests:\n")
print(interaction_lrt)
cat("\nMultivariable sensitivity rows:\n")
print(multivariable_rows)
