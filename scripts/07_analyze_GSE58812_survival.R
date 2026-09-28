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

table_dir <- file.path(project_root, "results", "tables")
figure_dir <- file.path(project_root, "results", "figures")
log_dir <- file.path(project_root, "logs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

score_path <- file.path(table_dir, "GSE58812_immune_reactivation_scores.csv")
survival_path <- file.path(table_dir, "GSE58812_TNBC_survival_metadata.csv")
if (!file.exists(score_path)) {
  stop("Missing score file: ", score_path)
}
if (!file.exists(survival_path)) {
  stop("Missing survival metadata file: ", survival_path)
}

scores <- fread(score_path)
survival_meta <- fread(survival_path)

required_score_cols <- c("sample_id", "immune_reactivation_score")
required_surv_cols <- c("geo_accession", "age_at_diag", "meta", "mfs_days", "death", "os_days")
missing_score <- setdiff(required_score_cols, names(scores))
missing_surv <- setdiff(required_surv_cols, names(survival_meta))
if (length(missing_score) > 0) {
  stop("Score file missing columns: ", paste(missing_score, collapse = ", "))
}
if (length(missing_surv) > 0) {
  stop("Survival metadata file missing columns: ", paste(missing_surv, collapse = ", "))
}
if (anyDuplicated(scores$sample_id) > 0) {
  stop("Duplicated sample_id values in score file")
}
if (anyDuplicated(survival_meta$geo_accession) > 0) {
  stop("Duplicated geo_accession values in survival metadata")
}

analysis_dt <- merge(
  survival_meta,
  scores,
  by.x = "geo_accession",
  by.y = "sample_id",
  all = FALSE
)

to_num <- function(x) suppressWarnings(as.numeric(x))
to_int <- function(x) suppressWarnings(as.integer(x))
analysis_dt[, `:=`(
  age_at_diag_num = to_num(age_at_diag),
  os_days_num = to_num(os_days),
  death_event = to_int(death),
  mfs_days_num = to_num(mfs_days),
  meta_event = to_int(meta),
  immune_reactivation_score_z = as.numeric(scale(immune_reactivation_score))
)]

median_score <- median(analysis_dt$immune_reactivation_score, na.rm = TRUE)
analysis_dt[, immune_reactivation_group := fifelse(
  immune_reactivation_score >= median_score,
  "High",
  "Low"
)]
analysis_dt[, immune_reactivation_group := factor(immune_reactivation_group, levels = c("Low", "High"))]

valid_event <- function(event) {
  !is.na(event) & event %in% c(0L, 1L)
}

endpoint_specs <- list(
  overall_survival = list(time = "os_days_num", event = "death_event"),
  metastasis_free_survival = list(time = "mfs_days_num", event = "meta_event")
)

run_endpoint <- function(dt, endpoint_name, time_col, event_col) {
  endpoint_dt <- dt[
    !is.na(get(time_col)) &
      get(time_col) > 0 &
      valid_event(get(event_col)) &
      !is.na(immune_reactivation_score_z) &
      !is.na(immune_reactivation_group)
  ]

  if (nrow(endpoint_dt) == 0) {
    stop("No analyzable rows for endpoint: ", endpoint_name)
  }

  endpoint_dt[, `:=`(
    endpoint_time = get(time_col),
    endpoint_event = get(event_col)
  )]

  cox_continuous <- coxph(
    Surv(endpoint_time, endpoint_event) ~ immune_reactivation_score_z,
    data = endpoint_dt
  )
  cox_group <- coxph(
    Surv(endpoint_time, endpoint_event) ~ immune_reactivation_group,
    data = endpoint_dt
  )
  age_adjusted_dt <- endpoint_dt[!is.na(age_at_diag_num)]
  run_age_adjusted <- nrow(age_adjusted_dt) == nrow(endpoint_dt)
  if (run_age_adjusted) {
    cox_continuous_age <- coxph(
      Surv(endpoint_time, endpoint_event) ~ immune_reactivation_score_z + age_at_diag_num,
      data = age_adjusted_dt
    )
    cox_group_age <- coxph(
      Surv(endpoint_time, endpoint_event) ~ immune_reactivation_group + age_at_diag_num,
      data = age_adjusted_dt
    )
  }

  extract_cox_rows <- function(model, model_name) {
    model_summary <- summary(model)
    conf <- model_summary$conf.int
    coefs <- model_summary$coefficients
    data.table(
      endpoint = endpoint_name,
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

  continuous_summary <- summary(cox_continuous)
  continuous_zph <- cox.zph(cox_continuous)
  group_summary <- summary(cox_group)
  group_zph <- cox.zph(cox_group)
  cox_rows <- list(
    extract_cox_rows(cox_continuous, "cox_continuous_per_1sd"),
    extract_cox_rows(cox_group, "cox_median_high_vs_low")
  )
  zph_row_list <- list(
    data.table(
      endpoint = endpoint_name,
      model = "cox_continuous_per_1sd",
      term = rownames(continuous_zph$table),
      chisq = unname(continuous_zph$table[, "chisq"]),
      p_value = unname(continuous_zph$table[, "p"])
    ),
    data.table(
      endpoint = endpoint_name,
      model = "cox_median_high_vs_low",
      term = rownames(group_zph$table),
      chisq = unname(group_zph$table[, "chisq"]),
      p_value = unname(group_zph$table[, "p"])
    )
  )
  if (run_age_adjusted) {
    continuous_age_zph <- cox.zph(cox_continuous_age)
    group_age_zph <- cox.zph(cox_group_age)
    cox_rows <- c(
      cox_rows,
      list(
        extract_cox_rows(cox_continuous_age, "cox_continuous_per_1sd_age_adjusted"),
        extract_cox_rows(cox_group_age, "cox_median_high_vs_low_age_adjusted")
      )
    )
    zph_row_list <- c(
      zph_row_list,
      list(
        data.table(
          endpoint = endpoint_name,
          model = "cox_continuous_per_1sd_age_adjusted",
          term = rownames(continuous_age_zph$table),
          chisq = unname(continuous_age_zph$table[, "chisq"]),
          p_value = unname(continuous_age_zph$table[, "p"])
        ),
        data.table(
          endpoint = endpoint_name,
          model = "cox_median_high_vs_low_age_adjusted",
          term = rownames(group_age_zph$table),
          chisq = unname(group_age_zph$table[, "chisq"]),
          p_value = unname(group_age_zph$table[, "p"])
        )
      )
    )
  }

  logrank_fit <- survdiff(
    Surv(endpoint_time, endpoint_event) ~ immune_reactivation_group,
    data = endpoint_dt
  )
  logrank_p <- 1 - pchisq(logrank_fit$chisq, length(logrank_fit$n) - 1)
  logrank_row <- data.table(
    endpoint = endpoint_name,
    test = "logrank_median_high_vs_low",
    n = nrow(endpoint_dt),
    events = sum(endpoint_dt[[event_col]] == 1, na.rm = TRUE),
    chisq = unname(logrank_fit$chisq),
    p_value = logrank_p
  )

  zph_rows <- rbindlist(zph_row_list, fill = TRUE)

  km_fit <- survfit(
    Surv(endpoint_time, endpoint_event) ~ immune_reactivation_group,
    data = endpoint_dt
  )
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
    ylab = ifelse(endpoint_name == "overall_survival", "Overall survival probability", "MFS probability"),
    risk.table.height = 0.28,
    ggtheme = theme_classic(base_size = 10)
  )

  ggsave(
    filename = file.path(figure_dir, paste0("GSE58812_", endpoint_name, "_KM_median_score.png")),
    plot = km_plot$plot,
    width = 5.5,
    height = 4.2,
    dpi = 300
  )
  ggsave(
    filename = file.path(figure_dir, paste0("GSE58812_", endpoint_name, "_KM_median_score_risktable.png")),
    plot = arrange_ggsurvplots(list(km_plot), print = FALSE),
    width = 5.8,
    height = 5.6,
    dpi = 300
  )

  list(
    cox = rbindlist(cox_rows, fill = TRUE),
    logrank = logrank_row,
    zph = zph_rows,
    analysis_rows = endpoint_dt[
      ,
      .(
        geo_accession,
        endpoint = endpoint_name,
        time_days = get(time_col),
        event = get(event_col),
        immune_reactivation_score,
        immune_reactivation_score_z,
        immune_reactivation_group
      )
    ]
  )
}

endpoint_results <- lapply(names(endpoint_specs), function(endpoint_name) {
  spec <- endpoint_specs[[endpoint_name]]
  run_endpoint(analysis_dt, endpoint_name, spec$time, spec$event)
})
names(endpoint_results) <- names(endpoint_specs)

cox_results <- rbindlist(lapply(endpoint_results, `[[`, "cox"), fill = TRUE)
logrank_results <- rbindlist(lapply(endpoint_results, `[[`, "logrank"), fill = TRUE)
zph_results <- rbindlist(lapply(endpoint_results, `[[`, "zph"), fill = TRUE)
endpoint_rows <- rbindlist(lapply(endpoint_results, `[[`, "analysis_rows"), fill = TRUE)

group_summary <- endpoint_rows[
  ,
  .(
    n = .N,
    events = sum(event == 1, na.rm = TRUE),
    median_time_days = median(time_days, na.rm = TRUE),
    median_score = median(immune_reactivation_score, na.rm = TRUE)
  ),
  by = .(endpoint, immune_reactivation_group)
]

merge_qc <- data.table(
  score_rows = nrow(scores),
  survival_rows = nrow(survival_meta),
  merged_rows = nrow(analysis_dt),
  unmatched_score_rows = sum(!scores$sample_id %in% survival_meta$geo_accession),
  unmatched_survival_rows = sum(!survival_meta$geo_accession %in% scores$sample_id),
  median_score_cutoff = median_score,
  os_analyzable = sum(!is.na(analysis_dt$os_days_num) & analysis_dt$os_days_num > 0 & valid_event(analysis_dt$death_event)),
  os_events = sum(analysis_dt$death_event == 1, na.rm = TRUE),
  mfs_analyzable = sum(!is.na(analysis_dt$mfs_days_num) & analysis_dt$mfs_days_num > 0 & valid_event(analysis_dt$meta_event)),
  mfs_events = sum(analysis_dt$meta_event == 1, na.rm = TRUE)
)

write.csv(cox_results, file.path(table_dir, "GSE58812_survival_cox_results.csv"), row.names = FALSE)
write.csv(logrank_results, file.path(table_dir, "GSE58812_survival_logrank_results.csv"), row.names = FALSE)
write.csv(zph_results, file.path(table_dir, "GSE58812_survival_cox_zph_results.csv"), row.names = FALSE)
write.csv(group_summary, file.path(table_dir, "GSE58812_survival_group_summary.csv"), row.names = FALSE)
write.csv(merge_qc, file.path(table_dir, "GSE58812_survival_merge_qc.csv"), row.names = FALSE)
fwrite(endpoint_rows, file.path(table_dir, "GSE58812_survival_analysis_rows.tsv.gz"), sep = "\t", quote = FALSE, na = "NA")

capture.output(sessionInfo(), file = file.path(log_dir, "sessionInfo_07_analyze_GSE58812_survival.txt"))

cat("GSE58812 survival prototype analysis complete\n")
print(merge_qc)
cat("\nCox results:\n")
print(cox_results)
cat("\nLog-rank results:\n")
print(logrank_results)
cat("\nCox proportional hazards tests:\n")
print(zph_results)
