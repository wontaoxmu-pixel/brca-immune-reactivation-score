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

scores <- fread(file.path(table_dir, "scRNA_compartment_signature_scores.csv"))
coverage <- fread(file.path(table_dir, "scRNA_signature_gene_coverage.csv"))

component_labels <- c(
  exhausted_cd8_t_cell = "Exhausted/cytotoxic CD8+ T-cell",
  antigen_presentation = "Antigen presentation",
  interferon_response = "Interferon response",
  immune_reactivation_score = "Composite score"
)

plot_dt <- melt(
  scores[subtype == "TNBC"],
  id.vars = c("subtype", "broad_compartment", "celltype_major", "n_cells"),
  measure.vars = names(component_labels),
  variable.name = "component",
  value.name = "score"
)
compartment_levels <- scores[
  subtype == "TNBC"
][order(-immune_reactivation_score), paste0(broad_compartment, "\n", celltype_major, "\n", "n=", n_cells)]
plot_dt[, component_label := component_labels[as.character(component)]]
plot_dt[, compartment_label := paste0(broad_compartment, "\n", celltype_major, "\n", "n=", n_cells)]
plot_dt[, compartment_label := factor(compartment_label, levels = rev(compartment_levels))]
plot_dt[, component_label := factor(component_label, levels = component_labels)]

p <- ggplot(plot_dt, aes(x = component_label, y = compartment_label, fill = score)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", score)), size = 3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = median(plot_dt$score, na.rm = TRUE), name = "log2 CPM+1\nscore") +
  labs(x = NULL, y = NULL, title = "GSE176078 TNBC scRNA-seq compartment-level signature scores") +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid = element_blank()
  )

ggsave(file.path(figure_dir, "scRNA_TNBC_compartment_signature_heatmap.png"), p, width = 8.5, height = 5.8, dpi = 300)
ggsave(file.path(figure_dir, "scRNA_TNBC_compartment_signature_heatmap.pdf"), p, width = 8.5, height = 5.8)

summary_table <- scores[
  subtype == "TNBC",
  .(
    Compartment = broad_compartment,
    `Major cell type` = celltype_major,
    `n cells` = n_cells,
    `Exhausted/cytotoxic CD8+ T-cell score` = sprintf("%.3f", exhausted_cd8_t_cell),
    `Antigen-presentation score` = sprintf("%.3f", antigen_presentation),
    `Interferon-response score` = sprintf("%.3f", interferon_response),
    `Composite score` = sprintf("%.3f", immune_reactivation_score)
  )
]
summary_table[, composite_score_numeric := as.numeric(`Composite score`)]
setorder(summary_table, -composite_score_numeric)
summary_table[, composite_score_numeric := NULL]

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

coverage_summary <- coverage[, .(
  requested_genes = .N,
  present_genes = sum(present),
  missing_genes = paste(gene[!present], collapse = ";")
), by = component]
coverage_summary[, missing_genes := fifelse(missing_genes == "", "None", missing_genes)]
fwrite(coverage_summary, file.path(table_dir, "scRNA_signature_gene_coverage_summary.csv"))

write_md_table(
  summary_table,
  file.path(manuscript_table_dir, "Table_7_scRNA_compartment_validation.md"),
  "Table 7. Lightweight scRNA-seq compartment validation in GSE176078 TNBC samples",
  "Scores were calculated from pseudobulked raw counts within author-annotated broad compartments. Values are mean log2(CPM + 1) scores for each predefined component and should be interpreted as compartment-level support rather than new cell-state discovery."
)

writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_21_plot_scRNA_compartment_validation.txt"))
cat("scRNA compartment validation plots and table written\n")
