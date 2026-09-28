#!/usr/bin/env Rscript

# Round-9 review-driven figure:
# Figure 9, two panels.
#   Panel A: nested-model C-index comparison (score+age, ImmuneScore+age, both)
#            per analysis set, annotated with the directional likelihood-ratio
#            test P value for "score added to ImmuneScore".
#   Panel B: per-tumor antigen-presentation component rank of the malignant
#            (Cancer Epithelial) compartment versus other compartments in
#            GSE176078 TNBC tumors (lower rank = higher AP).
# Panels are combined with base grid viewports (no patchwork/gridExtra needed).

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
  library(grid)
})

table_dir <- file.path(project_root, "results", "tables")
figure_dir <- file.path(project_root, "results", "figures")
log_dir <- file.path(project_root, "logs")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------- Panel A: incremental-value C-index ----------------
mc <- fread(file.path(table_dir, "incremental_value_model_comparison.csv"))

set_labels <- c(
  "GSE58812|tnbc|overall_survival" = "GSE58812 TNBC\nOS",
  "GSE58812|tnbc|metastasis_free_survival" = "GSE58812 TNBC\nMFS",
  "GSE96058|all_primary|overall_survival" = "GSE96058 all-primary\nOS",
  "GSE96058|pam50_basal|overall_survival" = "GSE96058 PAM50 Basal\nOS",
  "GSE96058|pathology_tnbc|overall_survival" = "GSE96058 pathology TNBC\nOS"
)
mc[, set_key := paste(dataset, analysis_set, endpoint, sep = "|")]
mc[, set_label := set_labels[set_key]]
mc[, set_label := factor(set_label, levels = set_labels)]

ci_long <- melt(
  mc[, .(set_label, lrt_score_added_to_immune_p,
         `score + age` = cindex_score,
         `ImmuneScore + age` = cindex_immune,
         `score + ImmuneScore + age` = cindex_both)],
  id.vars = c("set_label", "lrt_score_added_to_immune_p"),
  variable.name = "model", value.name = "cindex"
)
ci_long[, model := factor(model,
  levels = c("score + age", "ImmuneScore + age", "score + ImmuneScore + age"))]

# annotation: LRT P for score added to immune (placed above the "both" bar)
lrt_lab <- mc[, .(
  set_label,
  ymax = pmax(cindex_score, cindex_immune, cindex_both),
  label = sprintf("LRT P (score added) = %s",
                  ifelse(lrt_score_added_to_immune_p < 0.001, "<0.001",
                         sprintf("%.3f", lrt_score_added_to_immune_p)))
)]

panel_a <- ggplot(ci_long, aes(x = set_label, y = cindex, fill = model)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_text(aes(label = sprintf("%.3f", cindex)),
            position = position_dodge(width = 0.8),
            vjust = -0.4, size = 2.6) +
  geom_text(data = lrt_lab, inherit.aes = FALSE,
            aes(x = set_label, y = ymax + 0.045, label = label),
            size = 2.7, fontface = "italic") +
  scale_fill_manual(values = c("score + age" = "#C8102E",
                               "ImmuneScore + age" = "#1A3D7C",
                               "score + ImmuneScore + age" = "#B08D57")) +
  coord_cartesian(ylim = c(0.5, 0.85)) +
  labs(title = "A  Discrimination (Harrell C-index): the score adds nothing over ImmuneScore where protective",
       x = NULL, y = "C-index", fill = "Age-adjusted Cox model") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", size = 10.5),
        axis.text.x = element_text(size = 8.5))

# ---------------- Panel B: per-tumor malignant AP rank ----------------
ap <- fread(file.path(table_dir, "scRNA_malignant_AP_per_tumor_rank.csv"))
ap[, compartment_lab := fcase(
  broad_compartment == "malignant_epithelial", "Malignant (Cancer Epithelial)",
  broad_compartment == "t_nk", "T/NK",
  broad_compartment == "myeloid", "Myeloid",
  broad_compartment == "b_plasma", "B/plasma",
  broad_compartment == "stromal", "Stromal",
  broad_compartment == "normal_epithelial", "Normal epithelial",
  default = "Other"
)]
ap[, is_malignant := broad_compartment == "malignant_epithelial"]
# order tumors by malignant AP rank for readability
mal_order <- ap[is_malignant == TRUE][order(ap_rank), patient_id]
ap[, patient_id := factor(patient_id, levels = mal_order)]
ap_plot <- ap[!is.na(patient_id)]

panel_b <- ggplot(ap_plot, aes(x = patient_id, y = ap_rank)) +
  geom_point(aes(color = is_malignant, size = is_malignant), alpha = 0.85) +
  geom_text(data = ap_plot[is_malignant == TRUE],
            aes(label = ap_rank), color = "white", size = 2.5) +
  scale_y_reverse(breaks = 1:9) +
  scale_color_manual(values = c("FALSE" = "#9AA3B2", "TRUE" = "#C8102E"),
                     labels = c("Other compartments", "Malignant (Cancer Epithelial)"),
                     name = NULL) +
  scale_size_manual(values = c("FALSE" = 2.2, "TRUE" = 5), guide = "none") +
  labs(title = "B  Antigen-presentation rank per TNBC tumor: malignant cells are never the top AP source (rank 1 = highest AP)",
       x = "TNBC tumor (GSE176078)", y = "AP component rank\n(1 = highest)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", size = 10.5),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

# ---------------- Combine with base grid viewports ----------------
out_png <- file.path(figure_dir, "Figure_9_incremental_value_and_malignant_AP.png")
out_pdf <- file.path(figure_dir, "Figure_9_incremental_value_and_malignant_AP.pdf")

draw_both <- function() {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 1, heights = unit(c(1, 1), "null"))))
  print(panel_a, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(panel_b, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  popViewport()
}

png(out_png, width = 2400, height = 2200, res = 200)
draw_both()
dev.off()

pdf(out_pdf, width = 12, height = 11)
draw_both()
dev.off()

cat("Figure 9 written:\n  ", out_png, "\n  ", out_pdf, "\n")
writeLines(capture.output(sessionInfo()),
           file.path(log_dir, "sessionInfo_26_plot_incremental_value_and_malignant_AP.txt"))
cat("Figure 9 plotting complete\n")
