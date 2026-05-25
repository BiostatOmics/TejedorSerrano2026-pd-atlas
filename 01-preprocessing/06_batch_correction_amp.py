"""
Batch correction of AMP-PD dataset using Harmony (on PCA embeddings),
followed by re-clustering and tSNE.

Usage:
    python 05_batch_correction_amp.py <input_amp> <output_amp> <save_path>
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
parser = argparse.ArgumentParser(description="Harmony batch correction for AMP-PD dataset")
parser.add_argument("input_amp", type=str, help="Path to input h5ad")
parser.add_argument("output_amp", type=str, help="Output path for batch-corrected h5ad")
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
log_message(f"Loading AMP dataset from {args.input_amp}...")
amp = ad.read_h5ad(args.input_amp)

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
        sc.pl.tsne(adata_rand, color=color, ax=ax, show=False, s=size, alpha=1, palette=palette, frameon=True)
        plt.legend(ncol=ncol, bbox_to_anchor=(1, 1), frameon=False, fontsize=(6 if ncol > 4 else 16))
        plt.savefig(f"{base_file_name}_{color}.png", bbox_inches="tight")
        plt.close()

# PRE-CORRECTION PLOTS
to_plot = ["leiden_0_10", "participant_id", "sample_id", "set", "case_control", "brain_region"]
plot_tsne_panels(amp, f"{args.save_path}amppd_tsne_precorrected", to_plot)

# HARMONY BATCH CORRECTION
log_message("Running Harmony batch correction (by participant_id)...")
sc.external.pp.harmony_integrate(amp, key="participant_id")

# Replace default PCA with Harmony-corrected PCA, preserve original
amp.obsm["X_pca_uncorrected"] = amp.obsm["X_pca"]
amp.obsm["X_pca"] = amp.obsm["X_pca_harmony"]

# Re-run clustering pipeline on corrected embeddings
log_message("Re-running nearest neighbors and clustering on corrected embeddings...")
sc.pp.neighbors(amp, n_neighbors=30, n_pcs=50, metric="manhattan")
sc.tl.leiden(amp, resolution=0.1, flavor="leidenalg", key_added="leiden_0_10_corrected")
sc.tl.leiden(amp, resolution=0.25, flavor="leidenalg", key_added="leiden_0_25_corrected")
sc.tl.leiden(amp, resolution=1, flavor="leidenalg", key_added="leiden_1_00_corrected")

# Re-run tSNE on corrected embeddings
amp.obsm["X_tsne_uncorrected"] = amp.obsm["X_tsne"]
sc.tl.tsne(amp, n_pcs=50, perplexity=500, early_exaggeration=20)

# SAVE
log_message(f"Saving batch-corrected object to {args.output_amp}...")
amp.write_h5ad(args.output_amp, compression="gzip")

# POST-CORRECTION PLOTS
plot_tsne_panels(amp, f"{args.save_path}amppd_tsne_postcorrected", to_plot)
log_message("Done!")
