"""Shared Python utilities for logging and data export."""

import logging
import sys
import traceback
import numpy as np

# Configure timestamped logging to stdout
logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] - %(message)s",
    datefmt="%d-%m-%Y %H:%M:%S",
    handlers=[logging.StreamHandler(sys.stdout)],
)


def log_message(message):
    logging.info(message)


def handle_error(e):
    error_message = "---------------------------------------\nERROR: " + str(e)
    traceback_message = traceback.format_exc()
    system_info = (
        f"System Information:\n"
        f"  Python version: {sys.version}\n"
        f"  Platform: {sys.platform}\n"
        f"---------------------------------------\n"
    )
    full_message = f"{error_message}\nTraceback:\n{traceback_message}\n{system_info}"
    log_message(full_message)
    raise e


def beautify_array(data, col_names, row_names):
    data = np.column_stack((row_names, data))
    data = np.row_stack((col_names, data))
    return data


def save_pca_umap_data(adata, suffix, path):
    """Export PCA and UMAP coordinates and loadings to TSV files."""
    cell_names = adata.obs_names.tolist()
    gene_names = adata.var_names.tolist()

    pca_scores = adata.obsm["X_pca"].copy()
    pca_loadings = adata.varm["PCs"].copy()
    umap_scores = adata.obsm["X_umap"].copy()

    PC_names = ["Index"] + [f"PC_{i + 1}" for i in range(pca_scores.shape[1])]
    UMAP_names = ["Index"] + [f"UMAP_{i + 1}" for i in range(umap_scores.shape[1])]

    # Save PCA scores
    pca_scores_beautified = beautify_array(pca_scores, PC_names, cell_names)
    np.savetxt(
        f"{path}{suffix}_PCA_scores.tsv",
        pca_scores_beautified,
        delimiter="\t",
        fmt="%s",
    )

    # Save PCA loadings
    pca_loadings_beautified = beautify_array(pca_loadings, PC_names, gene_names)
    np.savetxt(
        f"{path}{suffix}_PCA_loadings.tsv",
        pca_loadings_beautified,
        delimiter="\t",
        fmt="%s",
    )

    # Save PCA variance and variance ratio
    variance = np.column_stack((PC_names[1:], adata.uns["pca"]["variance"]))
    np.savetxt(f"{path}{suffix}_PCA_variance.tsv", variance, delimiter="\t", fmt="%s")

    variance_ratio = np.column_stack((PC_names[1:], adata.uns["pca"]["variance_ratio"]))
    np.savetxt(
        f"{path}{suffix}_PCA_variance_ratio.tsv",
        variance_ratio,
        delimiter="\t",
        fmt="%s",
    )

    # Save UMAP scores
    umap_scores_beautified = beautify_array(umap_scores, UMAP_names, cell_names)
    np.savetxt(
        f"{path}{suffix}_UMAP_scores.tsv",
        umap_scores_beautified,
        delimiter="\t",
        fmt="%s",
    )
