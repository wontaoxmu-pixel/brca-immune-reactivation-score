#!/usr/bin/env Rscript

set.seed(123)

project_root <- normalizePath(getwd(), mustWork = TRUE)
if (basename(project_root) == "scripts") {
  project_root <- normalizePath(file.path(project_root, ".."), mustWork = TRUE)
}

dir.create(file.path(project_root, "logs"), recursive = TRUE, showWarnings = FALSE)
project_lib <- file.path(project_root, ".Rlib")
if (dir.exists(project_lib)) {
  .libPaths(c(project_lib, .libPaths()))
}

core_packages <- c(
  "tidyverse",
  "data.table",
  "survival",
  "survminer",
  "glmnet",
  "pROC",
  "ComplexHeatmap",
  "clusterProfiler",
  "GSVA",
  "GEOquery",
  "UCSCXenaTools",
  "SingleCellExperiment",
  "scater",
  "scran",
  "SingleR",
  "celldex",
  "limma",
  "edgeR",
  "DESeq2"
)

optional_packages <- c(
  "Seurat",
  "SeuratObject",
  "harmony",
  "CellChat",
  "NicheNet",
  "timeROC",
  "rms",
  "estimate",
  "immunedeconv",
  "BiocManager"
)

required_packages <- c(core_packages, optional_packages)

check_package <- function(pkg) {
  code <- sprintf(
    ".libPaths(c(%s, .libPaths())); suppressPackageStartupMessages(ok <- requireNamespace(%s, quietly = TRUE)); quit(status = ifelse(ok, 0, 1))",
    shQuote(project_lib),
    shQuote(pkg)
  )
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    args = c("--vanilla", "-e", shQuote(code)),
    stdout = FALSE,
    stderr = FALSE
  )
  identical(status, 0L)
}

pkg_status <- data.frame(
  package = required_packages,
  group = c(rep("core", length(core_packages)), rep("optional", length(optional_packages))),
  installed = vapply(required_packages, check_package, logical(1)),
  stringsAsFactors = FALSE
)

write.csv(pkg_status, file.path(project_root, "logs", "r_package_status.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(project_root, "logs", "sessionInfo.txt"))

cat("R environment check complete\n")
cat("Core packages installed:", sum(pkg_status$installed[pkg_status$group == "core"]), "/", length(core_packages), "\n")
cat("Optional packages installed:", sum(pkg_status$installed[pkg_status$group == "optional"]), "/", length(optional_packages), "\n")
cat("Missing core packages:\n")
print(pkg_status$package[pkg_status$group == "core" & !pkg_status$installed])
cat("Missing optional packages:\n")
print(pkg_status$package[pkg_status$group == "optional" & !pkg_status$installed])
