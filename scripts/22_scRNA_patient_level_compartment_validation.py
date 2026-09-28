#!/usr/bin/env python3

import math
import random
import sys
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from scipy.stats import wilcoxon
except Exception:  # pragma: no cover
    wilcoxon = None

random.seed(123)
np.random.seed(123)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed" / "scRNA"
TABLE_DIR = PROJECT_ROOT / "results" / "tables"
LOG_DIR = PROJECT_ROOT / "logs"
TABLE_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)

MTX = PROCESSED_DIR / "count_matrix_sparse.mtx"
GENES = PROCESSED_DIR / "count_matrix_genes.tsv"
BARCODES = PROCESSED_DIR / "count_matrix_barcodes.tsv"
METADATA = PROCESSED_DIR / "metadata.csv"

SIGNATURES = {
    "exhausted_cd8_t_cell": [
        "CD8A", "CD8B", "PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX",
        "CXCL13", "GZMB", "PRF1", "NKG7",
    ],
    "antigen_presentation": [
        "HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "TAPBP",
        "PSMB8", "PSMB9", "NLRC5",
    ],
    "interferon_response": [
        "IFNG", "STAT1", "IRF1", "CXCL9", "CXCL10", "GBP1", "GBP5",
        "ISG15", "IFIT1", "MX1",
    ],
}

MIN_CELLS_FOR_PAIR = 20


def broad_compartment(row: pd.Series) -> str:
    major = str(row["celltype_major"])
    minor = str(row["celltype_minor"])
    if major == "Cancer Epithelial":
        return "malignant_epithelial"
    if major == "Normal Epithelial":
        return "normal_epithelial"
    if major == "T-cells" or "NK" in minor:
        return "t_nk"
    if major in {"B-cells", "Plasmablasts"}:
        return "b_plasma"
    if major == "Myeloid":
        return "myeloid"
    if major in {"CAFs", "PVL", "Endothelial"}:
        return "stromal"
    return "other"


def read_lines(path: Path) -> list[str]:
    with path.open("rt") as handle:
        return [line.rstrip("\n") for line in handle]


def main() -> None:
    for path in [MTX, GENES, BARCODES, METADATA]:
        if not path.exists():
            raise FileNotFoundError(path)

    metadata = pd.read_csv(METADATA)
    genes = [g.upper() for g in read_lines(GENES)]
    barcodes = read_lines(BARCODES)

    cell_id_col = "Unnamed: 0"
    if len(barcodes) != metadata.shape[0]:
        raise ValueError(f"Barcode/metadata row mismatch: {len(barcodes)} vs {metadata.shape[0]}")
    if not np.array_equal(metadata[cell_id_col].astype(str).values, np.array(barcodes, dtype=str)):
        raise ValueError("Barcode order does not match metadata cell order")

    metadata["broad_compartment"] = metadata.apply(broad_compartment, axis=1)
    metadata["is_tnbc"] = metadata["subtype"].astype(str).eq("TNBC")
    metadata["group"] = (
        metadata["orig.ident"].astype(str)
        + "|"
        + metadata["broad_compartment"].astype(str)
        + "|"
        + metadata["celltype_major"].astype(str)
    )

    tnbc_metadata = metadata[metadata["is_tnbc"]].copy()
    group_levels = sorted(tnbc_metadata["group"].unique())
    group_to_idx = {g: i for i, g in enumerate(group_levels)}
    cell_group = metadata["group"].map(group_to_idx).fillna(-1).to_numpy(np.int32)

    cell_counts = tnbc_metadata.groupby(
        ["orig.ident", "broad_compartment", "celltype_major"], dropna=False
    ).agg(
        n_cells=(cell_id_col, "size"),
        total_umi=("nCount_RNA", "sum"),
        median_features=("nFeature_RNA", "median"),
        median_percent_mito=("percent.mito", "median"),
    ).reset_index()
    cell_counts.to_csv(TABLE_DIR / "scRNA_TNBC_patient_compartment_cell_counts.csv", index=False)

    requested = []
    for component, component_genes in SIGNATURES.items():
        for gene in component_genes:
            requested.append((component, gene))
    gene_to_row = {}
    for idx, gene in enumerate(genes, start=1):
        gene_to_row.setdefault(gene, idx)
    present = [(component, gene, gene_to_row.get(gene)) for component, gene in requested]
    selected = [(component, gene, row_idx) for component, gene, row_idx in present if row_idx is not None]
    selected_genes = [gene for _, gene, _ in selected]
    row_to_gene_idx = {row_idx: i for i, (_, _, row_idx) in enumerate(selected)}
    row_to_gene_idx_bytes = {
        str(row_idx).encode(): gene_idx
        for row_idx, gene_idx in row_to_gene_idx.items()
    }

    n_groups = len(group_levels)
    n_genes = len(selected_genes)
    count_sum = np.zeros((n_groups, n_genes), dtype=np.float64)
    detected_cells = np.zeros((n_groups, n_genes), dtype=np.int64)
    n_cells_by_group = tnbc_metadata.groupby("group").size().reindex(group_levels).to_numpy(np.int64)
    total_umi_by_group = tnbc_metadata.groupby("group")["nCount_RNA"].sum().reindex(group_levels).to_numpy(np.float64)

    nnz_seen = 0
    with MTX.open("rb") as handle:
        for raw in handle:
            if raw.startswith(b"%"):
                continue
            dims = raw.split()
            if len(dims) == 3:
                n_rows, n_cols, n_nnz = map(int, dims)
                break
        if n_rows != len(genes) or n_cols != len(barcodes):
            raise ValueError(f"Matrix dimension mismatch: {n_rows}x{n_cols} vs {len(genes)}x{len(barcodes)}")

        for raw in handle:
            nnz_seen += 1
            parts = raw.split()
            gene_pos = row_to_gene_idx_bytes.get(parts[0])
            if gene_pos is None:
                continue
            cell_idx = int(parts[1]) - 1
            group_idx = cell_group[cell_idx]
            if group_idx < 0:
                continue
            val = float(parts[2])
            count_sum[group_idx, gene_pos] += val
            detected_cells[group_idx, gene_pos] += 1
            if nnz_seen % 20_000_000 == 0:
                print(f"processed {nnz_seen}/{n_nnz} nonzero entries", flush=True)

    group_rows = []
    for group, group_idx in group_to_idx.items():
        patient_id, broad, major = group.split("|", 2)
        for gene_idx, (component, gene, _) in enumerate(selected):
            total_umi = total_umi_by_group[group_idx]
            cpm = (count_sum[group_idx, gene_idx] / total_umi * 1e6) if total_umi > 0 else np.nan
            group_rows.append({
                "patient_id": patient_id,
                "broad_compartment": broad,
                "celltype_major": major,
                "component": component,
                "gene": gene,
                "n_cells": int(n_cells_by_group[group_idx]),
                "total_umi": float(total_umi),
                "gene_counts": float(count_sum[group_idx, gene_idx]),
                "detected_cells": int(detected_cells[group_idx, gene_idx]),
                "detection_fraction": float(detected_cells[group_idx, gene_idx] / n_cells_by_group[group_idx]),
                "cpm": float(cpm),
                "log2_cpm_plus1": float(math.log2(cpm + 1)) if not math.isnan(cpm) else np.nan,
            })
    gene_expr = pd.DataFrame(group_rows)
    gene_expr.to_csv(TABLE_DIR / "scRNA_TNBC_patient_compartment_signature_gene_expression.csv", index=False)

    score = gene_expr.groupby(
        ["patient_id", "broad_compartment", "celltype_major", "component"],
        dropna=False,
    ).agg(
        n_cells=("n_cells", "first"),
        total_umi=("total_umi", "first"),
        genes_present=("gene", "nunique"),
        mean_detection_fraction=("detection_fraction", "mean"),
        component_log2_cpm_score=("log2_cpm_plus1", "mean"),
    ).reset_index()
    wide = score.pivot_table(
        index=["patient_id", "broad_compartment", "celltype_major", "n_cells", "total_umi"],
        columns="component",
        values="component_log2_cpm_score",
        aggfunc="first",
    ).reset_index()
    component_cols = list(SIGNATURES.keys())
    for col in component_cols:
        if col not in wide.columns:
            wide[col] = np.nan
    wide["immune_reactivation_score"] = wide[component_cols].mean(axis=1)
    wide = wide.sort_values(["patient_id", "broad_compartment", "celltype_major"])
    wide.to_csv(TABLE_DIR / "scRNA_TNBC_patient_compartment_signature_scores.csv", index=False)

    t_cells = wide[
        (wide["broad_compartment"] == "t_nk")
        & (wide["celltype_major"] == "T-cells")
        & (wide["n_cells"] >= MIN_CELLS_FOR_PAIR)
    ][["patient_id", "n_cells", "immune_reactivation_score", "exhausted_cd8_t_cell", "antigen_presentation", "interferon_response"]]
    malignant = wide[
        (wide["broad_compartment"] == "malignant_epithelial")
        & (wide["celltype_major"] == "Cancer Epithelial")
        & (wide["n_cells"] >= MIN_CELLS_FOR_PAIR)
    ][["patient_id", "n_cells", "immune_reactivation_score", "exhausted_cd8_t_cell", "antigen_presentation", "interferon_response"]]
    paired = t_cells.merge(malignant, on="patient_id", suffixes=("_t_cells", "_malignant"))
    paired["composite_difference_t_minus_malignant"] = (
        paired["immune_reactivation_score_t_cells"] - paired["immune_reactivation_score_malignant"]
    )
    paired["antigen_presentation_difference_t_minus_malignant"] = (
        paired["antigen_presentation_t_cells"] - paired["antigen_presentation_malignant"]
    )
    paired["interferon_difference_t_minus_malignant"] = (
        paired["interferon_response_t_cells"] - paired["interferon_response_malignant"]
    )
    paired.to_csv(TABLE_DIR / "scRNA_TNBC_patient_Tcell_vs_malignant_paired_scores.csv", index=False)

    ranks = []
    for patient_id, patient_dt in wide.groupby("patient_id"):
        ranked = patient_dt.sort_values("immune_reactivation_score", ascending=False).reset_index(drop=True)
        for rank, (_, row) in enumerate(ranked.iterrows(), start=1):
            ranks.append({
                "patient_id": patient_id,
                "rank": rank,
                "broad_compartment": row["broad_compartment"],
                "celltype_major": row["celltype_major"],
                "n_cells": int(row["n_cells"]),
                "immune_reactivation_score": float(row["immune_reactivation_score"]),
            })
    rank_dt = pd.DataFrame(ranks)
    rank_dt.to_csv(TABLE_DIR / "scRNA_TNBC_patient_compartment_composite_rankings.csv", index=False)

    diffs = paired["composite_difference_t_minus_malignant"].dropna().to_numpy()
    wilcoxon_p_greater = np.nan
    wilcoxon_p_two_sided = np.nan
    if wilcoxon is not None and len(diffs) > 0 and np.any(diffs != 0):
        wilcoxon_p_greater = float(wilcoxon(diffs, alternative="greater", zero_method="wilcox").pvalue)
        wilcoxon_p_two_sided = float(wilcoxon(diffs, alternative="two-sided", zero_method="wilcox").pvalue)

    t_rank = rank_dt[(rank_dt["broad_compartment"] == "t_nk") & (rank_dt["celltype_major"] == "T-cells")]
    m_rank = rank_dt[(rank_dt["broad_compartment"] == "malignant_epithelial") & (rank_dt["celltype_major"] == "Cancer Epithelial")]
    summary = pd.DataFrame([{
        "tnbc_patients_total": int(tnbc_metadata["orig.ident"].nunique()),
        "patients_with_t_cells_min_cells": int(t_cells["patient_id"].nunique()),
        "patients_with_malignant_min_cells": int(malignant["patient_id"].nunique()),
        "paired_patients_t_cells_and_malignant": int(paired["patient_id"].nunique()),
        "min_cells_threshold": MIN_CELLS_FOR_PAIR,
        "patients_t_cells_composite_greater_than_malignant": int((diffs > 0).sum()),
        "fraction_t_cells_composite_greater_than_malignant": float((diffs > 0).mean()) if len(diffs) else np.nan,
        "median_t_cells_composite": float(paired["immune_reactivation_score_t_cells"].median()) if len(paired) else np.nan,
        "median_malignant_composite": float(paired["immune_reactivation_score_malignant"].median()) if len(paired) else np.nan,
        "median_composite_difference_t_minus_malignant": float(np.median(diffs)) if len(diffs) else np.nan,
        "wilcoxon_p_greater_t_cells_gt_malignant": wilcoxon_p_greater,
        "wilcoxon_p_two_sided": wilcoxon_p_two_sided,
        "patients_t_cells_rank_1": int((t_rank["rank"] == 1).sum()),
        "patients_t_cells_rank_top2": int((t_rank["rank"] <= 2).sum()),
        "median_t_cells_composite_rank": float(t_rank["rank"].median()) if len(t_rank) else np.nan,
        "median_malignant_composite_rank": float(m_rank["rank"].median()) if len(m_rank) else np.nan,
    }])
    summary.to_csv(TABLE_DIR / "scRNA_TNBC_patient_level_robustness_summary.csv", index=False)

    with (LOG_DIR / "sessionInfo_22_scRNA_patient_level_compartment_validation.txt").open("wt") as log:
        log.write(f"python_version={sys.version}\n")
        log.write(f"numpy={np.__version__}\n")
        log.write(f"pandas={pd.__version__}\n")
        log.write(f"scipy_wilcoxon_available={wilcoxon is not None}\n")
        log.write(f"matrix={MTX}\n")
        log.write(f"selected_signature_genes={n_genes}\n")
        log.write(f"tnbc_patients={tnbc_metadata['orig.ident'].nunique()}\n")
        log.write(f"patient_compartment_groups={n_groups}\n")

    print("scRNA patient-level compartment validation complete")


if __name__ == "__main__":
    main()
