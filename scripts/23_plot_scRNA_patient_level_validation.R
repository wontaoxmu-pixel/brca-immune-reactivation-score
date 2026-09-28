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
manuscript_table_dir <- file.path(project_root, "manuscript", "tables")
log_dir <- file.path(project_root, "logs")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manuscript_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

paired <- fread(file.path(table_dir, "scRNA_TNBC_patient_Tcell_vs_malignant_paired_scores.csv"))
summary_dt <- fread(file.path(table_dir, "scRNA_TNBC_patient_level_robustness_summary.csv"))

plot_dt <- rbindlist(list(
  paired[, .(
    patient_id,
    compartment = "T-cells",
    n_cells = n_cells_t_cells,
    composite_score = immune_reactivation_score_t_cells
  )],
  paired[, .(
    patient_id,
    compartment = "Cancer Epithelial",
    n_cells = n_cells_malignant,
    composite_score = immune_reactivation_score_malignant
  )]
))
plot_dt[, compartment := factor(compartment, levels = c("Cancer Epithelial", "T-cells"))]

summary_row <- summary_dt[1]
subtitle <- sprintf(
  "%d paired tumors; T-cells > Cancer Epithelial in %d/%d; one-sided Wilcoxon P = %.3f",
  summary_row$paired_patients_t_cells_and_malignant,
  summary_row$patients_t_cells_composite_greater_than_malignant,
  summary_row$paired_patients_t_cells_and_malignant,
  summary_row$wilcoxon_p_greater_t_cells_gt_malignant
)

p <- ggplot(plot_dt, aes(x = compartment, y = composite_score, group = patient_id)) +
  geom_line(color = "#7A7A7A", linewidth = 0.45, alpha = 0.75) +
  geom_point(aes(fill = compartment), shape = 21, color = "black", size = 3, stroke = 0.35) +
  scale_fill_manual(values = c("Cancer Epithelial" = "#4C78A8", "T-cells" = "#D55E00")) +
  labs(
    x = NULL,
    y = "Composite score, log2(CPM + 1)",
    title = "Patient-level scRNA-seq robustness of TNBC compartment scores",
    subtitle = subtitle
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9),
    axis.text.x = element_text(size = 10)
  )

ggsave(file.path(figure_dir, "scRNA_TNBC_patient_Tcell_vs_malignant_paired_composite.png"), p, width = 7.2, height = 4.8, dpi = 300)
ggsave(file.path(figure_dir, "scRNA_TNBC_patient_Tcell_vs_malignant_paired_composite.pdf"), p, width = 7.2, height = 4.8)

format_p <- function(x) {
  ifelse(is.na(x), "NA", ifelse(x < 0.001, "<0.001", sprintf("%.3f", x)))
}

summary_table <- data.table(
  Metric = c(
    "TNBC patients assessed",
    "Patients with T-cells, >=20 cells",
    "Patients with Cancer Epithelial cells, >=20 cells",
    "Paired patients with both compartments",
    "Pairs with T-cell composite score > Cancer Epithelial score",
    "Median T-cell composite score",
    "Median Cancer Epithelial composite score",
    "Median paired difference, T-cells minus Cancer Epithelial",
    "One-sided paired Wilcoxon P value",
    "Two-sided paired Wilcoxon P value",
    "Patients where T-cells ranked first across compartments",
    "Median Cancer Epithelial composite rank"
  ),
  Value = c(
    as.character(summary_row$tnbc_patients_total),
    as.character(summary_row$patients_with_t_cells_min_cells),
    as.character(summary_row$patients_with_malignant_min_cells),
    as.character(summary_row$paired_patients_t_cells_and_malignant),
    sprintf(
      "%d/%d (%.1f%%)",
      summary_row$patients_t_cells_composite_greater_than_malignant,
      summary_row$paired_patients_t_cells_and_malignant,
      100 * summary_row$fraction_t_cells_composite_greater_than_malignant
    ),
    sprintf("%.3f", summary_row$median_t_cells_composite),
    sprintf("%.3f", summary_row$median_malignant_composite),
    sprintf("%.3f", summary_row$median_composite_difference_t_minus_malignant),
    format_p(summary_row$wilcoxon_p_greater_t_cells_gt_malignant),
    format_p(summary_row$wilcoxon_p_two_sided),
    sprintf("%d/%d", summary_row$patients_t_cells_rank_1, summary_row$tnbc_patients_total),
    sprintf("%.1f", summary_row$median_malignant_composite_rank)
  )
)

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
  summary_table,
  file.path(manuscript_table_dir, "Table_8_scRNA_patient_level_robustness.md"),
  "Table 8. Patient-level scRNA-seq robustness of TNBC compartment scores",
  "Scores were calculated from GSE176078 TNBC author-annotated compartments using pseudobulked signature-gene counts. The paired comparison includes samples with at least 20 cells in both T-cell and Cancer Epithelial compartments."
)

writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_23_plot_scRNA_patient_level_validation.txt"))
cat("scRNA patient-level validation plot and table written\n")
