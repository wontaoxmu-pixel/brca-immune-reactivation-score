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
  library(AnnotationDbi)
  library(hgu133plus2.db)
  library(ggplot2)
})

processed_dir <- file.path(project_root, "data", "processed")
table_dir <- file.path(project_root, "results", "tables")
figure_dir <- file.path(project_root, "results", "figures")
log_dir <- file.path(project_root, "logs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

expr_path <- file.path(processed_dir, "GSE58812_GPL570_series_expression.tsv.gz")
if (!file.exists(expr_path)) {
  stop("Missing GSE58812 expression matrix: ", expr_path)
}

signature_sets <- list(
  exhausted_cd8_t_cell = c("CD8A", "CD8B", "PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX", "CXCL13", "GZMB", "PRF1", "NKG7"),
  antigen_presentation = c("HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "TAPBP", "PSMB8", "PSMB9", "NLRC5"),
  interferon_response = c("IFNG", "STAT1", "IRF1", "CXCL9", "CXCL10", "GBP1", "GBP5", "ISG15", "IFIT1", "MX1")
)

expr <- fread(expr_path, check.names = FALSE)
stopifnot("feature_id" %in% names(expr))
sample_cols <- setdiff(names(expr), "feature_id")

probe_map <- AnnotationDbi::select(
  hgu133plus2.db,
  keys = unique(expr$feature_id),
  columns = c("SYMBOL"),
  keytype = "PROBEID"
)
probe_map <- as.data.table(probe_map)
probe_map <- probe_map[!is.na(SYMBOL) & SYMBOL != ""]
probe_map[, SYMBOL := toupper(SYMBOL)]
probe_map <- unique(probe_map, by = c("PROBEID", "SYMBOL"))

expr_mapped <- merge(expr, probe_map, by.x = "feature_id", by.y = "PROBEID", allow.cartesian = TRUE)
for (col in sample_cols) {
  set(expr_mapped, j = col, value = as.numeric(expr_mapped[[col]]))
}

gene_expr <- expr_mapped[
  ,
  lapply(.SD, mean, na.rm = TRUE),
  by = SYMBOL,
  .SDcols = sample_cols
]

mat <- as.matrix(gene_expr[, ..sample_cols])
mode(mat) <- "numeric"
rownames(mat) <- gene_expr$SYMBOL

z_mat <- t(scale(t(mat)))
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

  score_dt[, (set_name) := colMeans(z_mat[present, , drop = FALSE], na.rm = TRUE)]
}

score_cols <- names(signature_sets)
score_dt[, immune_reactivation_score := rowMeans(.SD, na.rm = TRUE), .SDcols = score_cols]
setorder(score_dt, -immune_reactivation_score)
score_dt[, immune_reactivation_rank := seq_len(.N)]

coverage_dt <- rbindlist(coverage_rows)

mapping_summary <- data.table(
  dataset = "GSE58812",
  platform = "GPL570",
  probes_total = nrow(expr),
  probes_with_symbol = uniqueN(probe_map$PROBEID),
  symbols_mapped = uniqueN(probe_map$SYMBOL),
  gene_symbols_after_collapse = nrow(gene_expr),
  samples = length(sample_cols)
)

write.csv(score_dt, file.path(table_dir, "GSE58812_immune_reactivation_scores.csv"), row.names = FALSE)
write.csv(coverage_dt, file.path(table_dir, "GSE58812_signature_gene_coverage.csv"), row.names = FALSE)
write.csv(mapping_summary, file.path(table_dir, "GSE58812_probe_mapping_summary.csv"), row.names = FALSE)
fwrite(gene_expr, file.path(processed_dir, "GSE58812_GPL570_gene_symbol_expression.tsv.gz"), sep = "\t", quote = FALSE, na = "NA")

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
  geom_tile(color = "white", linewidth = 0.15) +
  scale_fill_gradient2(low = "#3b6fb6", mid = "white", high = "#c43b3b", midpoint = 0) +
  labs(x = NULL, y = NULL, fill = "Score") +
  theme_minimal(base_size = 8) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank(),
    legend.position = "right"
  )

ggsave(
  filename = file.path(figure_dir, "GSE58812_immune_reactivation_scores_heatmap.png"),
  plot = heatmap_plot,
  width = 8,
  height = 2.4,
  dpi = 300
)

capture.output(sessionInfo(), file = file.path(log_dir, "sessionInfo_06_score_GSE58812_signature.txt"))

cat("GSE58812 signature scoring complete\n")
print(mapping_summary)
print(coverage_dt)
