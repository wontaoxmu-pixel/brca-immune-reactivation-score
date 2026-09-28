# Bulk-projected immune reactivation score in breast cancer

Analysis code for the PLOS ONE revision **“Immune-cell abundance coupling and incremental prognostic information of a bulk-projected immune reactivation score in breast cancer”** (PONE-D-26-40686).

## Repository contents

- [`scripts/`](scripts/README.md): complete author-generated analysis and figure scripts. The revised analyses are under `scripts/revision_R1/`; root-level scripts `00`–`26` are retained for historical provenance.
- [`tcrt_display_item_manifest.tsv`](tcrt_display_item_manifest.tsv): legacy manifest retained from an earlier submission workflow; it is not the display-item manifest for the PLOS ONE revision.
- [`LICENSE`](LICENSE): MIT license for the author-generated code. Third-party data and software retain their source terms.

## Reproduction

Follow [`scripts/README.md`](scripts/README.md). Reproduction of the reported R1 analyses requires the S1 and S2 Supporting Information files accompanying the manuscript; sample-level inputs and derived result tables are provided there and are not duplicated in this code repository. Public input sources and retrieval details are documented in the S2 File manifest. No controlled-access raw GSE176078 sequencing reads are redistributed.

The code provides retrospective analysis, not a clinical prediction tool. See the manuscript for the analysis design, endpoints, interpretation and limitations.
