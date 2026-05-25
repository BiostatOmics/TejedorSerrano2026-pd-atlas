"""
QC, normalization, HVG selection, PCA, clustering, and tSNE for the reference dataset.

Usage:
    python 04_preprocess_reference.py <input_ref> <output_ref>
"""

import argparse
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "utils")))
from my_utils import log_message

# CLI
parser = argparse.ArgumentParser(description="Preprocess reference dataset")
parser.add_argument("input_ref", type=str, help="Path to merged reference h5ad")
parser.add_argument("output_ref", type=str, help="Output path for preprocessed reference h5ad")
args = parser.parse_args()

sc.settings.verbosity = 4
sc.settings.n_jobs = 8

# LOAD DATA
log_message(f"Loading reference from {args.input_ref}...")
ref = ad.read_h5ad(args.input_ref)
print(ref)

# QUALITY CONTROL
log_message("Applying QC filters...")
sc.pp.filter_cells(ref, min_counts=1500)
sc.pp.filter_cells(ref, max_counts=110000)
sc.pp.filter_cells(ref, min_genes=1100)
sc.pp.filter_cells(ref, max_genes=12500)
ref = ref[ref.obs["fraction_mitochondrial"] < 0.02]
print(ref)

# DOUBLET DETECTION
log_message("Detecting doublets with Scrublet...")
try:
    ref = sc.pp.scrublet(ref, n_prin_comps=50, batch_key="sample_id", copy=True, verbose=True)
    log_message("Doublet detection completed (batch-corrected by sample).")
except Exception as e:
    log_message(f"Batch-aware Scrublet failed ({e}); running without batch correction.")
    ref = sc.pp.scrublet(ref, n_prin_comps=50, copy=True, verbose=True)

# NORMALIZATION
log_message("Normalizing (shifted-log: size factors + log1p)...")
counts = ref.X.copy()
sc.pp.normalize_total(ref, target_sum=None, inplace=True)
sc.pp.log1p(ref)
ref.layers["counts"] = counts
print(ref)

# HVG SELECTION
log_message("Selecting HVGs (top 20%, seurat_v3, batch-aware)...")
top_genes = np.floor(0.2 * ref.shape[1]).astype(int)
try:
    sc.pp.highly_variable_genes(
        ref, layer="counts", n_top_genes=top_genes,
        batch_key="sample_id", flavor="seurat_v3", span=0.3
    )
    log_message("HVG selection completed (batch-aware).")
except Exception as e:
    log_message(f"Batch-aware HVG selection failed ({e}); running without batch correction.")
    sc.pp.highly_variable_genes(
        ref, layer="counts", n_top_genes=top_genes,
        flavor="seurat_v3", span=0.3
    )

# PCA
log_message("Computing PCA (50 components)...")
sc.pp.pca(ref, n_comps=50)

# Checkpoint: save before clustering (clustering is expensive)
log_message(f"Checkpoint 1: saving to {args.output_ref}...")
ref.write_h5ad(args.output_ref, compression="gzip")

# CLUSTERING
log_message("Computing nearest-neighbor graph...")
sc.pp.neighbors(ref, n_neighbors=30, n_pcs=50, metric="manhattan")

log_message("Running Leiden clustering (resolution=0.1)...")
sc.tl.leiden(ref, resolution=0.1, flavor="leidenalg", key_added="leiden_0_10")

log_message(f"Checkpoint 2: saving to {args.output_ref}...")
ref.write_h5ad(args.output_ref, compression="gzip")

# tSNE
log_message("Running tSNE...")
sc.tl.tsne(ref, n_pcs=50, perplexity=500, early_exaggeration=20)

# SAVE
log_message(f"Saving preprocessed reference to {args.output_ref}...")
ref.write_h5ad(args.output_ref, compression="gzip")
log_message("Done!")
