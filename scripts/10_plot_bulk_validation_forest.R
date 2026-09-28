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
  library(ggplot2)
})

table_dir <- file.path(project_root, "results", "tables")
figure_dir <- file.path(project_root, "results", "figures")
log_dir <- file.path(project_root, "logs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

gse58812_path <- file.path(table_dir, "GSE58812_survival_cox_results.csv")
gse96058_path <- file.path(table_dir, "GSE96058_survival_cox_results.csv")
if (!file.exists(gse58812_path)) {
  stop("Missing GSE58812 Cox results: ", gse58812_path)
}
if (!file.exists(gse96058_path)) {
  stop("Missing GSE96058 Cox results: ", gse96058_path)
}

format_p <- function(p) {
  out <- rep(NA_character_, length(p))
  idx <- !is.na(p)
  out[idx] <- ifelse(p[idx] < 0.001, "<0.001", sprintf("%.3f", p[idx]))
  out
}

format_hr_ci <- function(hr, lo, hi) {
  sprintf("%.3f (%.3f-%.3f)", hr, lo, hi)
}

endpoint_label <- function(endpoint) {
  fifelse(
    endpoint == "overall_survival", "OS",
    fifelse(endpoint == "metastasis_free_survival", "MFS", endpoint)
  )
}

gse58812 <- fread(gse58812_path)
gse58812[, `:=`(
  dataset = "GSE58812",
  analysis_set = "tnbc",
  analysis_set_label = "TNBC"
)]

gse96058 <- fread(gse96058_path)
gse96058[, dataset := "GSE96058"]
analysis_set_labels <- c(
  all_primary = "All primary",
  pathology_tnbc = "Pathology TNBC",
  pam50_basal = "PAM50 Basal",
  pathology_tnbc_and_pam50_basal = "Pathology TNBC + PAM50 Basal",
  pathology_tnbc_or_pam50_basal = "Pathology TNBC or PAM50 Basal"
)
gse96058[, analysis_set_label := unname(analysis_set_labels[analysis_set])]
gse96058[is.na(analysis_set_label), analysis_set_label := analysis_set]

cox_all <- rbindlist(list(gse58812, gse96058), fill = TRUE)
score_terms <- c("immune_reactivation_score_z")
continuous_models <- c("cox_continuous_per_1sd", "cox_continuous_per_1sd_age_adjusted")

summary_dt <- cox_all[
  term %in% score_terms &
    model %in% continuous_models
]

summary_dt[, `:=`(
  endpoint_label = endpoint_label(endpoint),
  hr_ci = format_hr_ci(hr, ci_lower, ci_upper),
  p_formatted = format_p(p_value),
  direction = fcase(
    ci_upper < 1, "lower_risk",
    ci_lower > 1, "higher_risk",
    default = "not_significant"
  )
)]
model_labels <- c(
  cox_continuous_per_1sd = "Unadjusted continuous score per 1 SD",
  cox_continuous_per_1sd_age_adjusted = "Age-adjusted continuous score per 1 SD"
)
summary_dt[, model_label := unname(model_labels[model])]
summary_dt[is.na(model_label), model_label := model]

setcolorder(
  summary_dt,
  c(
    "dataset", "analysis_set", "analysis_set_label", "endpoint", "endpoint_label",
    "model", "model_label", "n", "events", "hr", "ci_lower", "ci_upper",
    "hr_ci", "p_value", "p_formatted", "concordance", "direction"
  )
)
setorder(summary_dt, dataset, analysis_set, endpoint, model)
fwrite(summary_dt, file.path(table_dir, "bulk_validation_continuous_cox_summary.csv"))

publication_dt <- copy(summary_dt)
publication_dt[, `:=`(
  cohort = paste(dataset, analysis_set_label, sep = ": "),
  endpoint = endpoint_label,
  model = model_label,
  `HR (95% CI)` = hr_ci,
  `P value` = p_formatted
)]
publication_dt <- publication_dt[
  ,
  .(cohort, endpoint, model, n, events, `HR (95% CI)`, `P value`, concordance, direction)
]
fwrite(publication_dt, file.path(table_dir, "bulk_validation_continuous_cox_publication_table.csv"))

forest_dt <- summary_dt[model == "cox_continuous_per_1sd_age_adjusted"]
forest_dt[, cohort_endpoint := paste0(dataset, ": ", analysis_set_label, " (", endpoint_label, ")")]
forest_order <- c(
  "GSE58812: TNBC (OS)",
  "GSE58812: TNBC (MFS)",
  "GSE96058: All primary (OS)",
  "GSE96058: Pathology TNBC (OS)",
  "GSE96058: PAM50 Basal (OS)",
  "GSE96058: Pathology TNBC + PAM50 Basal (OS)",
  "GSE96058: Pathology TNBC or PAM50 Basal (OS)"
)
forest_dt <- forest_dt[cohort_endpoint %in% forest_order]
forest_dt[, cohort_endpoint := factor(cohort_endpoint, levels = rev(forest_order))]
forest_dt[, plot_label := sprintf("HR %s; P %s", hr_ci, p_formatted)]
forest_dt[, signal := fcase(
  direction == "lower_risk", "Lower risk",
  direction == "higher_risk", "Higher risk",
  default = "Not significant"
)]
forest_dt[, signal := factor(signal, levels = c("Lower risk", "Higher risk", "Not significant"))]

if (nrow(forest_dt) == 0) {
  stop("No age-adjusted continuous Cox rows found for forest plot.")
}

x_min <- max(0.35, min(forest_dt$ci_lower, na.rm = TRUE) * 0.80)
label_x <- max(1.60, max(forest_dt$ci_upper, na.rm = TRUE) * 1.10)
x_max <- max(3.00, label_x * 1.75)

forest_plot <- ggplot(forest_dt, aes(x = hr, y = cohort_endpoint)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.45, color = "grey45") +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper, color = signal), width = 0.20, linewidth = 0.85, orientation = "y") +
  geom_point(aes(fill = signal), shape = 21, size = 3.4, color = "black", stroke = 0.25) +
  geom_text(aes(x = label_x, label = plot_label), hjust = 0, size = 3.1, color = "grey15") +
  scale_x_log10(
    limits = c(x_min, x_max),
    breaks = c(0.5, 0.75, 1.0, 1.25, 1.5, 2.0),
    labels = c("0.50", "0.75", "1.00", "1.25", "1.50", "2.00")
  ) +
  scale_color_manual(values = c("Lower risk" = "#2B8C6E", "Higher risk" = "#B84A4A", "Not significant" = "#6B7280")) +
  scale_fill_manual(values = c("Lower risk" = "#2B8C6E", "Higher risk" = "#B84A4A", "Not significant" = "#6B7280")) +
  labs(
    x = "Hazard ratio per 1 SD immune reactivation score",
    y = NULL,
    color = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    axis.text.y = element_text(size = 10, color = "grey10"),
    axis.text.x = element_text(color = "grey10"),
    axis.title.x = element_text(size = 11, margin = margin(t = 8)),
    plot.margin = margin(t = 8, r = 12, b = 8, l = 8)
  )

png_path <- file.path(figure_dir, "bulk_validation_continuous_cox_forest_age_adjusted.png")
pdf_path <- file.path(figure_dir, "bulk_validation_continuous_cox_forest_age_adjusted.pdf")
ggsave(png_path, forest_plot, width = 10.5, height = 4.8, dpi = 300, bg = "white")
ggsave(pdf_path, forest_plot, width = 10.5, height = 4.8, bg = "white")

writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_10_plot_bulk_validation_forest.txt"))

message("Wrote: ", file.path(table_dir, "bulk_validation_continuous_cox_summary.csv"))
message("Wrote: ", file.path(table_dir, "bulk_validation_continuous_cox_publication_table.csv"))
message("Wrote: ", png_path)
message("Wrote: ", pdf_path)
