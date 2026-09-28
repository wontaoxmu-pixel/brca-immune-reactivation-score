#!/usr/bin/env python3

# Round-9 review-driven analysis:
# Malignant-cell-intrinsic antigen-presentation (AP) exploration in GSE176078
# TNBC tumors. The Introduction motivates whether malignant cells retain
# antigen-presentation machinery, but earlier drafts never scored AP *within*
# the Cancer Epithelial compartment. This script builds that focused,
# descriptive analysis from the already-verified patient x compartment
# pseudobulk table (scRNA_TNBC_patient_compartment_signature_scores.csv) so it
# does not re-walk the 100k-cell matrix.
#
# Questions, all descriptive (no causal / prognostic claim):
#   1. Is the AP component detectable and how variable is it across TNBC tumors
#      within the malignant (Cancer Epithelial) compartment?
#   2. Within each tumor, where does the malignant compartment rank for the AP
#      component relative to immune/stromal compartments?
#   3. Is malignant-cell AP decoupled from the exhausted/cytotoxic CD8 and IFN
#      components (i.e. AP can be retained without local T-cell cytotoxicity)?

import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
TABLE_DIR = PROJECT_ROOT / "results" / "tables"
LOG_DIR = PROJECT_ROOT / "logs"
TABLE_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)

SCORES = TABLE_DIR / "scRNA_TNBC_patient_compartment_signature_scores.csv"
MIN_CELLS = 20


def main() -> None:
    if not SCORES.exists():
        raise FileNotFoundError(SCORES)
    df = pd.read_csv(SCORES)

    # Collapse to one row per patient x major cell type with >= MIN_CELLS.
    df = df[df["n_cells"] >= MIN_CELLS].copy()

    # ---- Q1: malignant-compartment AP distribution across TNBC tumors ----
    mal = df[df["broad_compartment"] == "malignant_epithelial"].copy()
    ap = mal["antigen_presentation"].dropna()
    q1 = pd.DataFrame([{
        "compartment": "malignant_epithelial (Cancer Epithelial)",
        "n_tumors": int(mal["patient_id"].nunique()),
        "ap_detected_tumors": int((ap > 0).sum()),
        "ap_median": float(np.median(ap)) if len(ap) else np.nan,
        "ap_min": float(ap.min()) if len(ap) else np.nan,
        "ap_max": float(ap.max()) if len(ap) else np.nan,
        "ap_iqr_low": float(np.percentile(ap, 25)) if len(ap) else np.nan,
        "ap_iqr_high": float(np.percentile(ap, 75)) if len(ap) else np.nan,
        "exhausted_cd8_median": float(mal["exhausted_cd8_t_cell"].median()),
        "interferon_median": float(mal["interferon_response"].median()),
    }])
    q1.to_csv(TABLE_DIR / "scRNA_malignant_AP_distribution.csv", index=False)

    # ---- Q2: per-tumor AP rank of the malignant compartment ----
    rank_rows = []
    for pid, sub in df.groupby("patient_id"):
        ranked = sub.sort_values("antigen_presentation", ascending=False).reset_index(drop=True)
        n_comp = len(ranked)
        for rank, (_, row) in enumerate(ranked.iterrows(), start=1):
            rank_rows.append({
                "patient_id": pid,
                "ap_rank": rank,
                "n_compartments": n_comp,
                "broad_compartment": row["broad_compartment"],
                "celltype_major": row["celltype_major"],
                "antigen_presentation": float(row["antigen_presentation"]),
            })
    rank_dt = pd.DataFrame(rank_rows)
    rank_dt.to_csv(TABLE_DIR / "scRNA_malignant_AP_per_tumor_rank.csv", index=False)
    mal_rank = rank_dt[rank_dt["broad_compartment"] == "malignant_epithelial"]

    # ---- Q3: per-tumor decoupling of malignant AP from exhausted-CD8 / IFN ----
    mal_decoupling = mal[[
        "patient_id", "n_cells",
        "antigen_presentation", "exhausted_cd8_t_cell", "interferon_response",
    ]].copy()
    mal_decoupling["ap_minus_exhausted_cd8"] = (
        mal_decoupling["antigen_presentation"] - mal_decoupling["exhausted_cd8_t_cell"]
    )
    mal_decoupling.to_csv(TABLE_DIR / "scRNA_malignant_AP_component_decoupling.csv", index=False)

    summary = pd.DataFrame([{
        "tumors_with_malignant_compartment": int(mal["patient_id"].nunique()),
        "ap_detected_in_all_malignant": bool((ap > 0).all()),
        "malignant_ap_median": float(np.median(ap)),
        "malignant_ap_range": f"{ap.min():.3f}-{ap.max():.3f}",
        "malignant_exhausted_cd8_median": float(mal["exhausted_cd8_t_cell"].median()),
        "tumors_malignant_ap_rank1": int((mal_rank["ap_rank"] == 1).sum()),
        "tumors_malignant_ap_top3": int((mal_rank["ap_rank"] <= 3).sum()),
        "median_malignant_ap_rank": float(mal_rank["ap_rank"].median()),
        "median_n_compartments_per_tumor": float(mal_rank["n_compartments"].median()),
        "tumors_ap_exceeds_exhausted_cd8_in_malignant": int((mal_decoupling["ap_minus_exhausted_cd8"] > 0).sum()),
    }])
    summary.to_csv(TABLE_DIR / "scRNA_malignant_AP_summary.csv", index=False)

    print("==== Malignant-compartment AP distribution ====")
    print(q1.to_string(index=False))
    print("\n==== Malignant AP per-tumor rank (AP component) ====")
    print(mal_rank[["patient_id", "ap_rank", "n_compartments", "antigen_presentation"]].to_string(index=False))
    print("\n==== Summary ====")
    print(summary.to_string(index=False))

    with (LOG_DIR / "sessionInfo_25_scRNA_malignant_AP_exploration.txt").open("wt") as log:
        log.write(f"python_version={sys.version}\n")
        log.write(f"numpy={np.__version__}\n")
        log.write(f"pandas={pd.__version__}\n")
        log.write(f"input={SCORES}\n")
        log.write(f"min_cells={MIN_CELLS}\n")

    print("\nMalignant-compartment AP exploration complete")


if __name__ == "__main__":
    main()
