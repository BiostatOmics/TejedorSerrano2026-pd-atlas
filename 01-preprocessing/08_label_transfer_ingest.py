"""
Cell type label transfer from reference to query using scanpy ingest.

Computes UMAPs for both datasets, identifies common genes, and transfers
cell type labels from reference to query via ingest.

Usage:
    python 07_label_transfer_ingest.py <input_query> <input_ref> <output_query> <save_path>
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
parser = argparse.ArgumentParser(description="Ingest label transfer")
parser.add_argument("input_query", type=str, help="Path to query h5ad")
parser.add_argument("input_ref", type=str, help="Path to reference h5ad")
parser.add_argument("output_query", type=str, help="Output path for annotated query h5ad")
parser.add_argument("save_path", type=str, help="Directory prefix for diagnostic plots")
args = parser.parse_args()

sc.settings.verbosity = 4
sc.settings.n_jobs = 8
sc.settings.autoshow = False
sc.set_figure_params(dpi=80, dpi_save=150, figsize=(14, 10), frameon=True, color_map="magma")

xkcd_colors = list(matplotlib.colors.XKCD_COLORS.values())
plt.rcParams.update({
    "axes.titlesize": 24, "axes.labelsize": 18,
    "xtick.labelsize": 14, "ytick.labelsize": 14,
    "legend.fontsize": 16, "font.size": 16,
})
size = 0.5

# LOAD DATA
log_message(f"Loading query from {args.input_query}...")
query = ad.read_h5ad(args.input_query)

log_message(f"Loading reference from {args.input_ref}...")
ref = ad.read_h5ad(args.input_ref)

# UMAP
log_message("Computing UMAPs...")
sc.tl.umap(query, min_dist=0.1)
sc.tl.umap(ref, min_dist=0.1)

# LABEL TRANSFER (via INGEST)
log_message("Identifying common genes between query and reference...")
common_genes = ref.var_names.intersection(query.var_names)
query_sub = query[:, common_genes]
ref_sub = ref[:, common_genes]

obs_to_transfer = ["ROIGroupFine", "supercluster_term", "cell_type"]
log_message(f"Transferring labels: {obs_to_transfer}...")
sc.tl.ingest(query_sub, ref_sub, obs=obs_to_transfer)

# Transfer annotations back to the full query object
for obs in obs_to_transfer:
    query.obs[obs] = query_sub.obs[obs]

# SAVE
log_message(f"Saving annotated query to {args.output_query}...")
query.write_h5ad(args.output_query, compression="gzip")
log_message("Done!")
