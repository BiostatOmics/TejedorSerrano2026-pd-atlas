"""
Harmony batch correction for the reference dataset, followed by re-clustering and tSNE.
Generates diagnostic tSNE plots before and after correction.

Usage:
    python 07_batch_correction_reference.py <input_ref> <output_ref> <save_path>
"""

import argparse
import numpy as np
import anndata as ad
import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "utils")))
from my_utils import log_message

# CLI
parser = argparse.ArgumentParser(description="Harmony batch correction for reference dataset")
parser.add_argument("input_ref", type=str, help="Path to preprocessed reference h5ad")
parser.add_argument("output_ref", type=str, help="Output path for batch-corrected reference h5ad")
parser.add_argument("save_path", type=str, help="Directory prefix for diagnostic plots")
args = parser.parse_args()

sc.settings.verbosity = 4
sc.settings.n_jobs = 8
sc.settings.autoshow = False
sc.set_figure_params(dpi=80, dpi_save=300, figsize=(14, 10), frameon=True, color_map="magma")

xkcd_colors = list(matplotlib.colors.XKCD_COLORS.values())
plt.rcParams.update({
    "axes.titlesize": 24, "axes.labelsize": 18,
    "xtick.labelsize": 14, "ytick.labelsize": 14,
    "legend.fontsize": 16, "font.size": 16,
})
size = 0.5

# LOAD DATA
log_message(f"Loading reference from {args.input_ref}...")
ref = ad.read_h5ad(args.input_ref)
print(ref)

# HELPER: GENERATE tSNE PLOTS
def plot_tsne_panels(adata, base_file_name, variables):
    np.random.seed(0)
    idx = np.random.permutation(adata.shape[0])
    adata_rand = adata[idx, :]
    for color in variables:
        ncat = len(adata.obs[color].cat.categories)
        palette = sns.husl_palette(ncat) if ncat <= 30 else xkcd_colors
        ncol = 1 if ncat <= 30 else (2 if ncat <= 60 else (3 if ncat <= 90 else (4 if ncat <= 150 else 6)))
        fig, ax = plt.subplots()
        sc.pl.tsne(adata_rand, color=color, ax=ax, show=False, s=size, alpha=1,
                   palette=palette, frameon=True)
        plt.legend(ncol=ncol, bbox_to_anchor=(1, 1), frameon=False,
                   fontsize=(6 if ncol > 4 else 16))
        plt.savefig(f"{base_file_name}_{color}.png", bbox_inches="tight")
        plt.close()

# PRE-CORRECTION PLOTS
to_plot = ["donor_id", "sample_id", "leiden_0_10", "dataset", "cell_type", "supercluster_term"]
plot_tsne_panels(ref, f"{args.save_path}ref_tsne_precorrected", to_plot)

# HARMONY BATCH CORRECTION
log_message("Running Harmony batch correction (by sample_id)...")
sc.external.pp.harmony_integrate(ref, key="sample_id")

# Replace default PCA with Harmony-corrected PCA, preserve original
ref.obsm["X_pca_uncorrected"] = ref.obsm["X_pca"]
ref.obsm["X_pca"] = ref.obsm["X_pca_harmony"]

# Re-run clustering pipeline on corrected embeddings
log_message("Re-running nearest neighbors and clustering on corrected embeddings...")
sc.pp.neighbors(ref, n_neighbors=30, n_pcs=50, metric="manhattan")
sc.tl.leiden(ref, resolution=0.1, flavor="leidenalg", key_added="leiden_0_10_corrected")
sc.tl.leiden(ref, resolution=0.25, flavor="leidenalg", key_added="leiden_0_25_corrected")
sc.tl.leiden(ref, resolution=1, flavor="leidenalg", key_added="leiden_1_00_corrected")

# Re-run tSNE on corrected embeddings
ref.obsm["X_tsne_uncorrected"] = ref.obsm["X_tsne"]
sc.tl.tsne(ref, n_pcs=50, perplexity=500, early_exaggeration=20)

# SAVE
log_message(f"Saving batch-corrected reference to {args.output_ref}...")
ref.write_h5ad(args.output_ref, compression="gzip")

# POST-CORRECTION PLOTS
plot_tsne_panels(ref, f"{args.save_path}ref_tsne_postcorrected", to_plot)
log_message("Done!")
