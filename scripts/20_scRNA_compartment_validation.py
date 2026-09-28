#!/usr/bin/env python3

import csv
import math
import random
import tarfile
import time
from pathlib import Path

import numpy as np
import pandas as pd

random.seed(123)
np.random.seed(123)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = PROJECT_ROOT / "data" / "raw"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed" / "scRNA"
TABLE_DIR = PROJECT_ROOT / "results" / "tables"
LOG_DIR = PROJECT_ROOT / "logs"
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
TABLE_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)

ARCHIVE = RAW_DIR / "GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz"
INNER_PREFIX = "Wu_etal_2021_BRCA_scRNASeq"
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


def extract_member_if_missing(member_name: str, out_path: Path) -> None:
    if out_path.exists() and out_path.stat().st_size > 0:
        return
    if not ARCHIVE.exists():
        raise FileNotFoundError(f"Missing archive: {ARCHIVE}")
    with tarfile.open(ARCHIVE, "r:gz") as tar:
        member = tar.getmember(f"{INNER_PREFIX}/{member_name}")
        src = tar.extractfile(member)
        if src is None:
            raise FileNotFoundError(member_name)
        tmp = out_path.with_suffix(out_path.suffix + ".tmp")
        with tmp.open("wb") as dst:
            while True:
                chunk = src.read(1024 * 1024)
                if not chunk:
                    break
                dst.write(chunk)
        tmp.replace(out_path)


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
    for member, out in [
        ("metadata.csv", METADATA),
        ("count_matrix_genes.tsv", GENES),
        ("count_matrix_barcodes.tsv", BARCODES),
        ("count_matrix_sparse.mtx", MTX),
    ]:
        extract_member_if_missing(member, out)

    metadata = pd.read_csv(METADATA)
    genes = [g.upper() for g in read_lines(GENES)]
    barcodes = read_lines(BARCODES)

    cell_id_col = "Unnamed: 0"
    if len(barcodes) != metadata.shape[0]:
        raise ValueError(f"Barcode/metadata row mismatch: {len(barcodes)} vs {metadata.shape[0]}")
    if not np.array_equal(metadata[cell_id_col].astype(str).values, np.array(barcodes, dtype=str)):
        raise ValueError("Barcode order does not match metadata cell order")

    metadata["broad_compartment"] = metadata.apply(broad_compartment, axis=1)
    metadata["group"] = (
        metadata["subtype"].astype(str)
        + "|"
        + metadata["broad_compartment"].astype(str)
        + "|"
        + metadata["celltype_major"].astype(str)
    )
    group_levels = sorted(metadata["group"].unique())
    group_to_idx = {g: i for i, g in enumerate(group_levels)}
    cell_group = metadata["group"].map(group_to_idx).to_numpy(np.int32)

    cell_counts = metadata.groupby(["subtype", "broad_compartment", "celltype_major"], dropna=False).agg(
        n_cells=(cell_id_col, "size"),
        total_umi=("nCount_RNA", "sum"),
        median_features=("nFeature_RNA", "median"),
        median_percent_mito=("percent.mito", "median"),
    ).reset_index()
    cell_counts.to_csv(TABLE_DIR / "scRNA_compartment_cell_counts.csv", index=False)

    requested = []
    for component, component_genes in SIGNATURES.items():
        for gene in component_genes:
            requested.append((component, gene))
    gene_to_row = {}
    for idx, gene in enumerate(genes, start=1):
        gene_to_row.setdefault(gene, idx)
    present = [(component, gene, gene_to_row.get(gene)) for component, gene in requested]
    coverage = pd.DataFrame([
        {
            "component": component,
            "gene": gene,
            "present": row_idx is not None,
            "matrix_row": row_idx if row_idx is not None else "",
        }
        for component, gene, row_idx in present
    ])
    coverage.to_csv(TABLE_DIR / "scRNA_signature_gene_coverage.csv", index=False)

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
    n_cells_by_group = metadata.groupby("group").size().reindex(group_levels).to_numpy(np.int64)
    total_umi_by_group = metadata.groupby("group")["nCount_RNA"].sum().reindex(group_levels).to_numpy(np.float64)

    started = time.time()
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
            val = float(parts[2])
            group_idx = cell_group[cell_idx]
            count_sum[group_idx, gene_pos] += val
            detected_cells[group_idx, gene_pos] += 1
            if nnz_seen % 20_000_000 == 0:
                elapsed = time.time() - started
                print(f"processed {nnz_seen}/{n_nnz} nonzero entries in {elapsed:.1f}s", flush=True)

    group_rows = []
    for group, group_idx in group_to_idx.items():
        subtype, broad, major = group.split("|", 2)
        for gene_idx, (component, gene, _) in enumerate(selected):
            total_umi = total_umi_by_group[group_idx]
            cpm = (count_sum[group_idx, gene_idx] / total_umi * 1e6) if total_umi > 0 else np.nan
            group_rows.append({
                "subtype": subtype,
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
    gene_expr.to_csv(TABLE_DIR / "scRNA_compartment_signature_gene_expression.csv", index=False)

    score = gene_expr.groupby(
        ["subtype", "broad_compartment", "celltype_major", "component"],
        dropna=False,
    ).agg(
        n_cells=("n_cells", "first"),
        total_umi=("total_umi", "first"),
        genes_present=("gene", "nunique"),
        mean_detection_fraction=("detection_fraction", "mean"),
        component_log2_cpm_score=("log2_cpm_plus1", "mean"),
    ).reset_index()
    wide = score.pivot_table(
        index=["subtype", "broad_compartment", "celltype_major", "n_cells", "total_umi"],
        columns="component",
        values="component_log2_cpm_score",
        aggfunc="first",
    ).reset_index()
    component_cols = list(SIGNATURES.keys())
    for col in component_cols:
        if col not in wide.columns:
            wide[col] = np.nan
    wide["immune_reactivation_score"] = wide[component_cols].mean(axis=1)
    wide = wide.sort_values(["subtype", "broad_compartment", "celltype_major"])
    wide.to_csv(TABLE_DIR / "scRNA_compartment_signature_scores.csv", index=False)

    tnbc = wide[wide["subtype"] == "TNBC"].copy()
    rank_rows = []
    for col in component_cols + ["immune_reactivation_score"]:
        tmp = tnbc.sort_values(col, ascending=False)
        for rank, (_, row) in enumerate(tmp.iterrows(), start=1):
            rank_rows.append({
                "subtype": row["subtype"],
                "component": col,
                "rank": rank,
                "broad_compartment": row["broad_compartment"],
                "celltype_major": row["celltype_major"],
                "n_cells": int(row["n_cells"]),
                "score": float(row[col]),
            })
    pd.DataFrame(rank_rows).to_csv(TABLE_DIR / "scRNA_TNBC_compartment_signature_rankings.csv", index=False)

    with (LOG_DIR / "sessionInfo_20_scRNA_compartment_validation.txt").open("wt") as log:
        log.write(f"python_version={__import__('sys').version}\n")
        log.write(f"numpy={np.__version__}\n")
        log.write(f"pandas={pd.__version__}\n")
        log.write(f"archive={ARCHIVE}\n")
        log.write(f"matrix={MTX}\n")
        log.write(f"matrix_rows={len(genes)}\n")
        log.write(f"matrix_cols={len(barcodes)}\n")
        log.write(f"selected_signature_genes={n_genes}\n")
        log.write(f"groups={n_groups}\n")

    print("scRNA compartment validation aggregation complete")


if __name__ == "__main__":
    main()
