#!/usr/bin/env Rscript

set.seed(123)

project_root <- normalizePath(file.path(getwd()), mustWork = TRUE)
if (basename(project_root) == "scripts") {
  project_root <- normalizePath(file.path(project_root, ".."), mustWork = TRUE)
}

dir.create(file.path(project_root, "results", "tables"), recursive = TRUE, showWarnings = FALSE)

manifest <- data.frame(
  dataset = c("GSE176078", "GSE161529", "TCGA-BRCA", "METABRIC", "GEO_TNBC_bulk", "TNBC_immunotherapy"),
  source = c("NCBI GEO", "NCBI GEO", "GDC/UCSC Xena", "cBioPortal/original METABRIC", "NCBI GEO", "GEO/clinical trial repositories"),
  role = c("Primary scRNA-seq discovery", "scRNA-seq confirmation", "Bulk validation", "External bulk validation", "External bulk validation", "Optional treatment-response validation"),
  status = c("confirmed", "confirmed", "to_download", "to_download", "to_identify", "to_identify"),
  key_reference = c(
    "PMID 34493872; DOI 10.1038/s41588-021-00911-1",
    "GEO accession confirmed; publication metadata to finalize",
    "PMID 23000897; DOI 10.1038/nature11412",
    "PMID 22522925; DOI 10.1038/nature10983",
    "pending",
    "pending"
  ),
  planned_use = c(
    "Define TNBC immune reactivation phenotype",
    "Confirm immune cell states and marker stability",
    "Project signature and test survival/immune associations",
    "Independent survival validation",
    "Additional TNBC-specific validation",
    "Optional future immune-checkpoint response validation; not included in the current analysis"
  ),
  stringsAsFactors = FALSE
)

out_file <- file.path(project_root, "results", "tables", "dataset_manifest.csv")
write.csv(manifest, out_file, row.names = FALSE)

cat("Dataset manifest written to:", out_file, "\n")
capture.output(sessionInfo(), file = file.path(project_root, "logs", "sessionInfo_01_build_dataset_manifest.txt"))
