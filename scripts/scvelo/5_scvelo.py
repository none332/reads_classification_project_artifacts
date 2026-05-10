#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pathlib import Path

import scvelo as scv

##################################
# This script runs scVelo analysis using processed count matrices of the method and extracted matrices of velocyto. In addition, this script provide a comparison between our baseline method data and precursor-ambiguous merged data.
##################################

# Construct function to write information in the log file.
def log(message, log_file):
    with open(log_file, "a") as f:
        f.write(str(message) + "\n")


# Construct function to check whether the input file exists.
def require_file(path):
    path = Path(path)

    if not path.is_file():
        raise FileNotFoundError(f"Required file not found: {path}")


# Construct function to read the matrix (rows = genes, columns = cells).
def read_matrix(path, name, log_file):
    path = Path(path)
    require_file(path)

    df = pd.read_csv(path, sep="\t", index_col=0)
    # Convert the matrix to numeric values
    df = df.apply(pd.to_numeric, errors="coerce")

    log(f"Read: {name}: shape={df.shape}", log_file)

    return df


# Construct function to clean the count matrix. 
def sanitize_matrix(df, name, log_file):
    df = df.copy()

    # Convert row and column names into strings.
    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)

    # Sum duplicated genes and rows.
    if df.index.duplicated().any():
        log(f"Warn: {name}: duplicated genes found, summing duplicated rows.", log_file)
        df = df.groupby(df.index).sum()

    if df.columns.duplicated().any():
        log(f"Warn: {name}: duplicated cells found, summing duplicated columns.", log_file)
        df = df.T.groupby(level=0).sum().T

    # Replace NA, inf, and negative values to 0.
    if df.isna().values.any():
        log(f"Warn: {name}: NA values found, filling with 0.", log_file)
        df = df.fillna(0)

    df = df.replace([np.inf, -np.inf], 0)

    if (df.values < 0).any():
        log(f"Warn: {name}: negative values found, setting them to 0.", log_file)
        df[df < 0] = 0

    return df


# Construct function to read embedding coordinates (prepared in advance, can be founded in root/data/scvelo directories).
# If import external embbeding file, it is expected to have cells as rows and have two columns.
def read_embedding(path, log_file):
    path = Path(path)
    require_file(path)

    embedding = pd.read_csv(path, index_col=0)
    embedding.index = embedding.index.astype(str)

    if embedding.shape[1] < 2:
        raise ValueError("Embedding file must contain at least two columns.")

    embedding = embedding.iloc[:, :2].copy()
    embedding = embedding.apply(pd.to_numeric, errors="coerce")

    if embedding.isna().values.any():
        raise ValueError("Embedding file contains NA or non-numeric coordinates.")

    return embedding

# Construct function to read cell color information (prepared in advance, can be founded in root/data/scvelo directories).
# If the file does not exist, all cells will be assigned to one group.
def read_cell_colors(path, log_file):
    path = Path(path)

    if not path.is_file():
        log("Info: cell color file not found; all cells will be assigned to one group.", log_file)
        return None

    cell_colors = pd.read_csv(path, index_col=0)
    cell_colors.index = cell_colors.index.astype(str)

    if cell_colors.shape[1] == 0:
        log("Warn: cell color file has no columns; ignored.", log_file)
        return None

    return cell_colors

# Construct function to build AnnData from spliced and unspliced matrices.
def build_adata(spliced_df, unspliced_df, embedding_df, cell_colors_df=None):
    # Input: gene-by-cell, AnnData expects: cell-by-gene -> matrices need to be transposed.
    spliced = spliced_df.T.copy()
    unspliced = unspliced_df.T.copy()

    adata = ad.AnnData(X=spliced.values.astype(np.float32))

    adata.layers["spliced"] = spliced.values.astype(np.float32)
    adata.layers["unspliced"] = unspliced.values.astype(np.float32)

    adata.obs_names = spliced.index.astype(str)
    adata.var_names = spliced.columns.astype(str)

    adata.obsm["X_embedding"] = embedding_df.loc[adata.obs_names].values.astype(np.float32)

    if cell_colors_df is not None:
        group_col = cell_colors_df.columns[0]
        groups = cell_colors_df.reindex(adata.obs_names)[group_col].fillna("unknown")
        adata.obs["group"] = groups.astype(str).astype("category")
    else:
        adata.obs["group"] = "all_cells"
        adata.obs["group"] = adata.obs["group"].astype("category")

    return adata

# Construct function to save embedding coordinates together with velocity vectors.
def save_embedding_and_velocity(adata, prefix, outdir):
    emb = pd.DataFrame(
        adata.obsm["X_embedding"],
        index=adata.obs_names,
        columns=["embedding1", "embedding2"]
    )

    if "velocity_embedding" in adata.obsm:
        vel = pd.DataFrame(
            adata.obsm["velocity_embedding"],
            index=adata.obs_names,
            columns=["velo_embedding1", "velo_embedding2"]
        )
        out = pd.concat([emb, vel], axis=1)
    else:
        out = emb
    # when no calculated velocity embedding

    out.to_csv(outdir / f"{prefix}_embedding_with_velocity.csv")


# Construct function to remove invalid (NA, inf) numerical values after normalization.
def clean_adata_values(adata):
    X = np.asarray(adata.X, dtype=np.float32)
    X[np.isnan(X)] = 0
    X[np.isinf(X)] = 0
    adata.X = X

    for layer in ["spliced", "unspliced"]:
        arr = np.asarray(adata.layers[layer], dtype=np.float32)
        arr[np.isnan(arr)] = 0
        arr[np.isinf(arr)] = 0
        adata.layers[layer] = arr

    return adata


# Function to generate plots:
# (1) embedding scatter by group
# (2) velocity stream by group
# (3) velocity confidence stream
def save_basic_plots(adata, prefix, outdir, log_file):
    try:
        scv.pl.scatter(
            adata,
            basis="embedding",
            color="group",
            size=80,
            legend_loc="right",
            show=False
        )
        plt.savefig(outdir / f"{prefix}_embedding_by_group.png", dpi=300, bbox_inches="tight")
        plt.close()
    except Exception as e:
        if log_file is not None:
            log(f"Warn: {prefix}: embedding group plot failed: {e}", log_file)

    try:
        scv.pl.velocity_embedding_stream(
            adata,
            basis="embedding",
            color="group",
            size=80,
            legend_loc="right",
            figsize=(7, 7),
            show=False
        )
        plt.savefig(outdir / f"{prefix}_velocity_stream_by_group.png", dpi=300, bbox_inches="tight")
        plt.close()
    except Exception as e:
        if log_file is not None:
            log(f"Warn: {prefix}: velocity stream group plot failed: {e}", log_file)

    if "velocity_confidence" in adata.obs.columns:
        try:
            scv.pl.velocity_embedding_stream(
                adata,
                basis="embedding",
                color="velocity_confidence",
                cmap="viridis",
                size=80,
                figsize=(7, 7),
                show=False
            )
            plt.savefig(outdir / f"{prefix}_velocity_stream_by_confidence.png", dpi=300, bbox_inches="tight")
            plt.close()
        except Exception as e:
            if log_file is not None:
                log(f"Warn: {prefix}: velocity confidence plot failed: {e}", log_file)


# Construct function to run scVelo for one dataset.
# workflow:
#  (1) remove zero, filter genes by min_shared_counts, normalize, calculate moments
#  (2) estimate RNA velocity, build graph and calculate confidence
#  (3) save AnnData, summary, embedding, and plots
def run_one(
    adata,
    prefix,
    outdir,
    log_file,
    velocity_mode="deterministic",
    min_shared_counts=20,
    n_pcs=30,
    n_neighbors=15
):
    log(f"Run: {prefix}", log_file)

    if adata.n_obs < 3:
        raise ValueError(f"{prefix}: too few cells for scVelo analysis.")

    gene_mask = np.array(
        adata.layers["spliced"].sum(axis=0) + adata.layers["unspliced"].sum(axis=0)
    ).flatten() > 0

    adata = adata[:, gene_mask].copy()
    n_genes_after_nonzero_filter = adata.n_vars

    if adata.n_vars == 0:
        raise ValueError(f"{prefix}: no genes left after nonzero-gene filtering.")

    scv.pp.filter_genes(adata, min_shared_counts=min_shared_counts)
    n_genes_after_min_shared_counts = adata.n_vars

    if adata.n_vars == 0:
        raise ValueError(f"{prefix}: no genes left after min_shared_counts filtering.")

    scv.pp.normalize_per_cell(adata)
    adata = clean_adata_values(adata)

    n_pcs_eff = min(n_pcs, adata.n_obs - 1, adata.n_vars - 1)
    n_neighbors_eff = min(n_neighbors, adata.n_obs - 1)

    if n_pcs_eff < 2:
        n_pcs_eff = 2

    if n_neighbors_eff < 2:
        n_neighbors_eff = 2

    scv.pp.moments(
        adata,
        n_pcs=n_pcs_eff,
        n_neighbors=n_neighbors_eff
    )

    scv.tl.velocity(adata, mode=velocity_mode)
    scv.tl.velocity_graph(adata)
    scv.tl.velocity_confidence(adata)

    # calculate velocity vectors on the external embedding
    scv.tl.velocity_embedding(adata, basis="embedding")

    summary_row = {
        "dataset": prefix,
        "n_cells": adata.n_obs,
        "n_genes_after_nonzero_filter": n_genes_after_nonzero_filter,
        "n_genes_after_min_shared_counts": n_genes_after_min_shared_counts,
        "spliced_total_after_filter": float(np.sum(adata.layers["spliced"])),
        "unspliced_total_after_filter": float(np.sum(adata.layers["unspliced"])),
        "n_pcs_used": n_pcs_eff,
        "n_neighbors_used": n_neighbors_eff
    }

    if "velocity_confidence" in adata.obs.columns:
        summary_row["mean_velocity_confidence"] = float(adata.obs["velocity_confidence"].mean())
        summary_row["median_velocity_confidence"] = float(adata.obs["velocity_confidence"].median())

    if "velocity_length" in adata.obs.columns:
        summary_row["mean_velocity_length"] = float(adata.obs["velocity_length"].mean())
        summary_row["median_velocity_length"] = float(adata.obs["velocity_length"].median())

    adata.obs.to_csv(outdir / f"{prefix}_obs.csv")

    emb = pd.DataFrame(
        adata.obsm["X_embedding"],
        index=adata.obs_names,
        columns=["embedding1", "embedding2"]
    )
    emb.to_csv(outdir / f"{prefix}_embedding.csv")

    save_embedding_and_velocity(adata, prefix, outdir)

    save_basic_plots(
        adata=adata,
        prefix=prefix,
        outdir=outdir,
        log_file=log_file
    )

    adata.write(outdir / f"{prefix}.h5ad")

    log(
        f"Done: {prefix}: cells={adata.n_obs}, genes={adata.n_vars}, "
        f"mean_confidence={summary_row.get('mean_velocity_confidence', 'NA')}",
        log_file
    )

    return adata, summary_row

# Main function.
def main():

    # Get the location of this script.
    script_dir = Path(__file__).resolve().parent
    # Back to the root of the project.
    project_root = (script_dir / "../../").resolve()

    parser = argparse.ArgumentParser(
        description="Run scVelo analysis using processed method and velocyto matrices."
    )

    parser.add_argument(
        "--matrix-root",
        default=str(project_root / "results/scvelo/3_processed_matrices"),
        help="Directory containing processed matrices."
    )

    parser.add_argument(
        "--outdir",
        default=str(project_root / "results/scvelo/4_scvelo_final"),
        help="Output directory for final scVelo results."
    )

    parser.add_argument(
        "--velocity-mode",
        default="deterministic",
        choices=["deterministic", "stochastic"],
        help="Velocity mode used by scv.tl.velocity."
    )

    parser.add_argument(
        "--min-shared-counts",
        type=int,
        default=20,
        help="Minimum shared counts used by scv.pp.filter_genes."
    )

    parser.add_argument(
        "--n-pcs",
        type=int,
        default=30,
        help="Number of PCs used by scv.pp.moments."
    )

    parser.add_argument(
        "--n-neighbors",
        type=int,
        default=15,
        help="Number of neighbors used by scv.pp.moments."
    )

    args = parser.parse_args()

    matrix_root = Path(args.matrix_root).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    log_file = outdir / "run.log"
    log_file.write_text("Start scVelo analysis\n")

    scv.settings.verbosity = 1

    # input files from processed matrix directory (both the method and velocyto)
    method_spliced_file = matrix_root / "method_spliced.gene_name.tsv"
    method_unspliced_file = matrix_root / "method_unspliced.gene_name.tsv"
    method_ambiguous_file = matrix_root / "method_ambiguous.gene_name.tsv"

    velo_spliced_file = matrix_root / "velocyto_merged.spliced.tsv"
    velo_unspliced_file = matrix_root / "velocyto_merged.unspliced.tsv"
    velo_ambiguous_file = matrix_root / "velocyto_merged.ambiguous.tsv"

    embedding_file = project_root / "data/scvelo/embedding_data.csv"
    cell_colors_file = project_root / "data/scvelo/cell_colors_data.csv"

    log(f"Matrix directory: {matrix_root}", log_file)
    log(f"Output directory: {outdir}", log_file)

    # read method matrices
    m_s = sanitize_matrix(
        read_matrix(method_spliced_file, "method_spliced", log_file),
        "method_spliced",
        log_file
    )

    m_u = sanitize_matrix(
        read_matrix(method_unspliced_file, "method_unspliced", log_file),
        "method_unspliced",
        log_file
    )

    m_a = sanitize_matrix(
        read_matrix(method_ambiguous_file, "method_ambiguous", log_file),
        "method_ambiguous",
        log_file
    )

    # read and clean velocyto matrices
    v_s = sanitize_matrix(
        read_matrix(velo_spliced_file, "velocyto_spliced", log_file),
        "velocyto_spliced",
        log_file
    )

    v_u = sanitize_matrix(
        read_matrix(velo_unspliced_file, "velocyto_unspliced", log_file),
        "velocyto_unspliced",
        log_file
    )

    v_a = sanitize_matrix(
        read_matrix(velo_ambiguous_file, "velocyto_ambiguous", log_file),
        "velocyto_ambiguous",
        log_file
    )

    embedding = read_embedding(embedding_file, log_file)
    cell_colors = read_cell_colors(cell_colors_file, log_file)

    # Keep only cells and genes shared by all matrices to make method and velocyto results comparable.
    common_cells = sorted(
        set(m_s.columns)
        & set(m_u.columns)
        & set(m_a.columns)
        & set(v_s.columns)
        & set(v_u.columns)
        & set(v_a.columns)
        & set(embedding.index)
    )

    common_genes = sorted(
        set(m_s.index)
        & set(m_u.index)
        & set(m_a.index)
        & set(v_s.index)
        & set(v_u.index)
        & set(v_a.index)
    )

    if len(common_cells) == 0:
        raise ValueError("No common cells found across matrices and embedding.")

    if len(common_genes) == 0:
        raise ValueError("No common genes found across method and velocyto matrices.")

    log(f"Common cells: {len(common_cells)}", log_file)
    log(f"Common genes: {len(common_genes)}", log_file)

    # subset all matrices to the common cell and gene group
    m_s = m_s.loc[common_genes, common_cells]
    m_u = m_u.loc[common_genes, common_cells]
    m_a = m_a.loc[common_genes, common_cells]

    v_s = v_s.loc[common_genes, common_cells]
    v_u = v_u.loc[common_genes, common_cells]
    v_a = v_a.loc[common_genes, common_cells]

    embedding = embedding.loc[common_cells]

    if cell_colors is not None:
        cell_colors = cell_colors.reindex(common_cells)

    # save the final cells and genes used for comparison
    pd.Series(common_cells, name="cell").to_csv(outdir / "common_cells_used.csv", index=False)
    pd.Series(common_genes, name="gene").to_csv(outdir / "common_genes_used.csv", index=False)

    intersection_summary = pd.DataFrame([{
        "n_common_cells": len(common_cells),
        "n_common_genes": len(common_genes),
        "method_spliced_total_after_intersection": float(m_s.values.sum()),
        "method_unspliced_total_after_intersection": float(m_u.values.sum()),
        "method_ambiguous_total_after_intersection": float(m_a.values.sum()),
        "velocyto_spliced_total_after_intersection": float(v_s.values.sum()),
        "velocyto_unspliced_total_after_intersection": float(v_u.values.sum()),
        "velocyto_ambiguous_total_after_intersection": float(v_a.values.sum()),
    }])

    intersection_summary.to_csv(outdir / "intersection_summary.csv", index=False)

    # build AnnData objects
    adata_method = build_adata(m_s, m_u, embedding, cell_colors)
    adata_method_plus_ambiguous = build_adata(m_s, m_u + m_a, embedding, cell_colors)
    adata_velocyto = build_adata(v_s, v_u, embedding, cell_colors)

    # run scvelo for each method
    summary_rows = []

    _, row = run_one(
        adata=adata_method,
        prefix="method_baseline",
        outdir=outdir,
        log_file=log_file,
        velocity_mode=args.velocity_mode,
        min_shared_counts=args.min_shared_counts,
        n_pcs=args.n_pcs,
        n_neighbors=args.n_neighbors
    )
    summary_rows.append(row)

    # This is the combination of precursor and ambiguous data for comparison with baseline method.
    _, row = run_one(
        adata=adata_method_plus_ambiguous,
        prefix="method_plus_ambiguous_to_unspliced",
        outdir=outdir,
        log_file=log_file,
        velocity_mode=args.velocity_mode,
        min_shared_counts=args.min_shared_counts,
        n_pcs=args.n_pcs,
        n_neighbors=args.n_neighbors
    )
    summary_rows.append(row)

    _, row = run_one(
        adata=adata_velocyto,
        prefix="velocyto_baseline",
        outdir=outdir,
        log_file=log_file,
        velocity_mode=args.velocity_mode,
        min_shared_counts=args.min_shared_counts,
        n_pcs=args.n_pcs,
        n_neighbors=args.n_neighbors
    )
    summary_rows.append(row)

    summary_df = pd.DataFrame(summary_rows)
    summary_file = outdir / "scvelo_compare_summary.csv"
    summary_df.to_csv(summary_file, index=False)

    log(f"Summary file: {summary_file}", log_file)
    log("All analyses finished.", log_file)

    print(f"Finished scVelo analysis. Summary: {summary_file}")


if __name__ == "__main__":
    main()