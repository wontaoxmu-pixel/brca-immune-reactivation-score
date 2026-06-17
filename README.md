  ## Public datasets

  Raw datasets are available from the Gene Expression Omnibus:

  - GSE176078
  - GSE58812
  - GSE96058

  No newly generated sequencing data are included in this repository.

  ## Repository contents

  - `scripts/`: analysis scripts used for data processing, scoring, survival
  analysis, sensitivity analyses, figure/table generation, and TCRT submission-
  package preparation.
  - `tcrt_display_item_manifest.tsv`: mapping of selected main and supplementary
  display items for the TCRT submission package.

  ## Reproducibility

  Scripts were run from the project root in the order documented in `scripts/
  README.md`.

  The original analyses used R 4.2.0 with a project-local R package library.
  Each R script sets `set.seed(123)` where applicable and writes session
  information logs.

  ## Data availability

  Raw public datasets are available from GEO under GSE176078, GSE58812, and
  GSE96058.

  Derived non-identifiable result files and analysis scripts are made available
  for reproducibility.

  ## License

  This repository is released under the MIT License.
