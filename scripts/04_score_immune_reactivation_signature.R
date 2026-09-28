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

processed_dir <- file.path(project_root, "data", "processed")
table_dir <- file.path(project_root, "results", "tables")
figure_dir <- file.path(project_root, "results", "figures")
log_dir <- file.path(project_root, "logs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

counts_path <- file.path(processed_dir, "GSE176078_bulkRNAseq_raw_counts.tsv.gz")
if (!file.exists(counts_path)) {
  stop("Missing parsed count matrix: ", counts_path)
}

signature_sets <- list(
  exhausted_cd8_t_cell = c("CD8A", "CD8B", "PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX", "CXCL13", "GZMB", "PRF1", "NKG7"),
  antigen_presentation = c("HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "TAPBP", "PSMB8", "PSMB9", "NLRC5"),
  interferon_response = c("IFNG", "STAT1", "IRF1", "CXCL9", "CXCL10", "GBP1", "GBP5", "ISG15", "IFIT1", "MX1")
)

counts <- fread(counts_path, check.names = FALSE)
stopifnot("gene_symbol" %in% names(counts))

sample_cols <- setdiff(names(counts), "gene_symbol")
counts[, gene_symbol := toupper(gene_symbol)]
counts <- counts[!is.na(gene_symbol) & gene_symbol != ""]

collapsed <- counts[, lapply(.SD, sum, na.rm = TRUE), by = gene_symbol, .SDcols = sample_cols]
mat <- as.matrix(collapsed[, ..sample_cols])
mode(mat) <- "numeric"
rownames(mat) <- collapsed$gene_symbol

library_sizes <- colSums(mat, na.rm = TRUE)
log_cpm <- log2(t(t(mat) / library_sizes * 1e6) + 1)

z_mat <- t(scale(t(log_cpm)))
z_mat[is.na(z_mat)] <- 0

score_dt <- data.table(sample_id = colnames(z_mat))
coverage_rows <- list()

for (set_name in names(signature_sets)) {
  genes <- unique(toupper(signature_sets[[set_name]]))
  present <- intersect(genes, rownames(z_mat))
  missing <- setdiff(genes, rownames(z_mat))

  coverage_rows[[length(coverage_rows) + 1]] <- data.table(
    signature = set_name,
    requested_genes = length(genes),
    present_genes = length(present),
    missing_genes = paste(missing, collapse = ";"),
    present_gene_list = paste(present, collapse = ";")
  )

  if (length(present) == 0) {
    score_dt[, (set_name) := NA_real_]
  } else {
    score_dt[, (set_name) := colMeans(z_mat[present, , drop = FALSE], na.rm = TRUE)]
  }
}

score_cols <- names(signature_sets)
score_dt[, immune_reactivation_score := rowMeans(.SD, na.rm = TRUE), .SDcols = score_cols]
setorder(score_dt, -immune_reactivation_score)
score_dt[, immune_reactivation_rank := seq_len(.N)]

coverage_dt <- rbindlist(coverage_rows)

write.csv(score_dt, file.path(table_dir, "GSE176078_immune_reactivation_scores.csv"), row.names = FALSE)
write.csv(coverage_dt, file.path(table_dir, "GSE176078_signature_gene_coverage.csv"), row.names = FALSE)

plot_dt <- melt(
  score_dt,
  id.vars = c("sample_id", "immune_reactivation_rank"),
  measure.vars = c(score_cols, "immune_reactivation_score"),
  variable.name = "signature",
  value.name = "score"
)
plot_dt[, sample_id := factor(sample_id, levels = score_dt$sample_id)]
plot_dt[, signature := factor(
  signature,
  levels = c(score_cols, "immune_reactivation_score"),
  labels = c("Exhausted CD8 T cell", "Antigen presentation", "IFN response", "Composite")
)]

heatmap_plot <- ggplot(plot_dt, aes(x = sample_id, y = signature, fill = score)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient2(low = "#3b6fb6", mid = "white", high = "#c43b3b", midpoint = 0) +
  labs(x = NULL, y = NULL, fill = "Score") +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 7),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank(),
    legend.position = "right"
  )

ggsave(
  filename = file.path(figure_dir, "GSE176078_immune_reactivation_scores_heatmap.png"),
  plot = heatmap_plot,
  width = 8,
  height = 2.6,
  dpi = 300
)

capture.output(sessionInfo(), file = file.path(log_dir, "sessionInfo_04_score_immune_reactivation_signature.txt"))

cat("Immune reactivation signature scoring complete\n")
cat("Samples scored:", nrow(score_dt), "\n")
print(coverage_dt)
