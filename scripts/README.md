# Analysis scripts and reproduction guide

This folder is the versioned source code for manuscript PONE-D-26-40686. The revision analyses are in `revision_R1/`; root-level scripts `00`–`26` are retained historical/provenance scripts. The R1 scripts implement the analyses described in the revised manuscript. Sample-level intermediate inputs and result tables are not stored in this GitHub repository; obtain the manuscript Supporting Information files S1 and S2 and extract them into the same project root before reproducing the reported revision analyses. Public source datasets remain available from GEO/cBioPortal as documented in the manifests. Controlled-access raw GSE176078 reads are not included.

## Reproduce the reported R1 results from the supplied supporting files

Use R 4.2 or later with `data.table` and `survival`, and Python 3.10 or later with `numpy`, `pandas`, `matplotlib` and `Pillow`. Package versions used by the analyses are recorded in `logs/revision_R1/*sessionInfo*` in S2 File. Start from the project root after extracting S1 and S2 File, preserving their directory structure.

```sh
Rscript scripts/revision_R1/03_survival_uncertainty.R
Rscript scripts/revision_R1/10_verify_core_results.R
python scripts/revision_R1/05_make_figures.py
```

The first command refits the four-cohort survival and bootstrap analyses (2,000 paired draws; it may take time). The verification script compares the regenerated core estimates with the retained analysis records. The figure script reads the supporting-file result tables and generates the manuscript figures. `10_verify_core_results.R` is a consistency check on the regenerated outputs; it does not replace the analysis step.

Optional analyses and source rebuilding:

- `revision_R1/00_fetch_sources.py` retrieves public expression/clinical sources and official MCP-counter/singscore resources listed in the source manifest.
- `revision_R1/01_score_cohorts.R` rebuilds scores when source expression inputs are available; existing cached analysis RDS files are skipped.
- `revision_R1/02_singlecell.py` processes the public GSE176078 processed matrix; its raw sequencing reads are access-controlled by the source repository and are not part of this package.
- `revision_R1/09_audit_records.py`, `12_exploratory_KM.R` and `14_time_effect_details.R` produce supplementary audit, exploratory KM and time-dependent sensitivity outputs.
- `revision_R1/17_check_figure_exports.py` checks the exported TIFF color mode, compression and resolution metadata.

Root-level scripts `00`–`26` document the earlier exploratory workflow and are retained for provenance. Their filenames and historical labels should not be interpreted as prospective validation. Scripts containing `install` are dependency-installation helpers: inspect them and install only the required packages for your environment; they are not part of the reproduction command sequence. Third-party methods retain their own terms and citations. The project code is distributed under the repository MIT license; third-party data and software remain subject to their source terms.

## Data and reporting boundaries

The repository contains author-generated analysis code and instructions. Sample-level inputs and detailed result tables are supplied with the manuscript as Supporting Information S1/S2, rather than duplicated here. GEO and cBioPortal sources and retrieval details are listed in `raw_input_manifest.csv` within S2 File. The analyses are retrospective and exploratory; reproducing them does not convert associations into causal or treatment-response claims. See the manuscript for endpoint definitions, exclusions, model specifications and limitations.
