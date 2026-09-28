#!/usr/bin/env Rscript

set.seed(123)

project_root <- normalizePath(getwd(), mustWork = TRUE)
if (basename(project_root) == "scripts") {
  project_root <- normalizePath(file.path(project_root, ".."), mustWork = TRUE)
}

project_lib <- file.path(project_root, ".Rlib")
dir.create(project_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(project_lib, .libPaths()))

log_dir <- file.path(project_root, "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(log_dir, "install_deconvolution_dependencies.log")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

cat("Project root:", project_root, "\n")
cat("Library paths:\n")
print(.libPaths())

cran_repo <- "https://mirrors.tuna.tsinghua.edu.cn/CRAN"
options(repos = c(CRAN = cran_repo))

ensure_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, lib = project_lib, repos = cran_repo, type = "binary")
  }
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package still unavailable after install attempt: ", pkg)
  }
}

ensure_package("remotes")
ensure_package("curl")

if (!requireNamespace("MCPcounter", quietly = TRUE)) {
  cat("Installing MCPcounter from GitHub ebecht/MCPcounter, subdir Source\n")
  github_result <- try(
    remotes::install_github(
      "ebecht/MCPcounter",
      subdir = "Source",
      lib = project_lib,
      upgrade = "never",
      dependencies = c("Depends", "Imports")
    ),
    silent = TRUE
  )
  if (inherits(github_result, "try-error") && !requireNamespace("MCPcounter", quietly = TRUE)) {
    cat("GitHub install failed; installing MCPcounter Zenodo source archive\n")
    mcp_tar <- file.path(tempdir(), "MCPcounter_1.1.0.tar.gz")
    utils::download.file(
      "https://zenodo.org/records/61372/files/MCPcounter_1.1.0.tar.gz?download=1",
      destfile = mcp_tar,
      mode = "wb",
      quiet = FALSE
    )
    install.packages(
      mcp_tar,
      lib = project_lib,
      repos = NULL,
      type = "source"
    )
  }
}

if (!requireNamespace("estimate", quietly = TRUE)) {
  cat("Installing estimate from R-Forge source\n")
  install.packages(
    "estimate",
    lib = project_lib,
    repos = "http://r-forge.r-project.org",
    type = "source",
    dependencies = c("Depends", "Imports")
  )
}

status <- data.frame(
  package = c("MCPcounter", "estimate"),
  installed = vapply(c("MCPcounter", "estimate"), requireNamespace, logical(1), quietly = TRUE),
  stringsAsFactors = FALSE
)
print(status)

if (!all(status$installed)) {
  stop("One or more deconvolution dependencies failed to install")
}

writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_17_install_deconvolution_dependencies.txt"))
cat("Deconvolution dependencies installed successfully\n")
