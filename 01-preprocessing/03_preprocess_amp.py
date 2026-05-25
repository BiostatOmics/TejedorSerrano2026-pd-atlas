"""
Preprocess AMP-PD h5ad: gene metadata transfer, normalization, HVG selection, PCA,
clustering, and tSNE.

Usage:
    python 03_preprocess_amp.py <input_amp> <input_gene_metadata> <output_amp>
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
parser = argparse.ArgumentParser(description="Preprocess AMP-PD dataset")
parser.add_argument("input_amp", type=str, help="Path to input h5ad")
parser.add_argument("input_gene_metadata", type=str, help="Path to gene metadata TSV")
parser.add_argument("output_amp", type=str, help="Output path for processed h5ad")
args = parser.parse_args()

sc.settings.verbosity = 4
sc.settings.n_jobs = 8

# LOAD DATA
log_message(f"Loading AMP dataset from {args.input_amp}...")
amp = ad.read_h5ad(args.input_amp)
print(amp)

log_message(f"Loading gene metadata from {args.input_gene_metadata}...")
gene_metadata = pd.read_csv(args.input_gene_metadata, sep="\t")
print(gene_metadata)

# GENE METADATA TRANSFER
log_message("Transferring gene metadata...")
gene_metadata.index = gene_metadata["original_ensembl_id"]
gene_metadata = gene_metadata[["hgnc_symbol", "description"]]
amp.var = amp.var.merge(
    gene_metadata, "inner", validate="one_to_one", left_index=True, right_index=True
)
print(amp.var)

# DOUBLET DETECTION
log_message("Detecting doublets with scrublet...")
try:
    amp = sc.pp.scrublet(amp, n_prin_comps=50, batch_key="participant_id", copy=True, verbose=True)
    log_message("Doublet detection completed (batch-corrected by participant)")
except Exception as e1:
    log_message(f"Batch-corrected scrublet failed: {e1}. Running without batch correction...")
    amp = sc.pp.scrublet(amp, n_prin_comps=50, copy=True, verbose=True)
    log_message("Doublet detection completed (no batch correction)")

# QC METRICS
log_message("Computing QC metrics...")
sc.pp.calculate_qc_metrics(amp, percent_top=[20], inplace=True)

# NORMALIZATION
log_message("Normalizing...")
counts = amp.X.copy()  # preserve raw counts before normalization

sc.pp.normalize_total(amp, target_sum=None, inplace=True)
cpm = amp.X.copy()  # size-factor corrected counts

sc.pp.log1p(amp)

amp.layers["counts"] = counts
amp.layers["cpm"] = cpm
log_message("Normalization complete.")

# HVG SELECTION
log_message("Selecting highly variable genes (top 20%)...")
top_genes = np.floor(0.2 * amp.shape[1]).astype(int)
try:
    # seurat_v3: selects HVGs consistent across donors, suitable for case-control designs
    sc.pp.highly_variable_genes(
        amp, layer="counts", n_top_genes=top_genes, batch_key="participant_id",
        flavor="seurat_v3", span=0.3
    )
    log_message("HVG selection completed (batch-corrected by participant)")
except Exception as e2:
    log_message(f"Batch-corrected HVG selection failed: {e2}. Running without batch correction...")
    sc.pp.highly_variable_genes(
        amp, layer="counts", n_top_genes=top_genes, flavor="seurat_v3", span=0.3
    )
    log_message("HVG selection completed (no batch correction)")

# PCA
log_message("Running PCA (50 components)...")
sc.pp.pca(amp, n_comps=50)

# Checkpoint: save before clustering (allows re-running from this point)
log_message(f"CHECKPOINT: Saving processed object to {args.output_amp}...")
amp.write_h5ad(args.output_amp, compression="gzip")

# CLUSTERING
log_message("Computing nearest neighbors (k=30, manhattan metric)...")
sc.pp.neighbors(amp, n_neighbors=30, n_pcs=50, metric="manhattan")

log_message("Running Leiden clustering at resolutions 0.1, 0.25, 1.0...")
sc.tl.leiden(amp, resolution=0.1, flavor="leidenalg", key_added="leiden_0_10")
sc.tl.leiden(amp, resolution=0.25, flavor="leidenalg", key_added="leiden_0_25")
sc.tl.leiden(amp, resolution=1, flavor="leidenalg", key_added="leiden_1_00")

# Checkpoint: save after clustering
log_message(f"CHECKPOINT: Saving processed object to {args.output_amp}...")
amp.write_h5ad(args.output_amp, compression="gzip")

# tSNE
log_message("Running tSNE...")
sc.tl.tsne(amp, n_pcs=50, perplexity=500, early_exaggeration=20)

# SAVE FINAL OBJECT
log_message(f"Saving final processed object to {args.output_amp}...")
amp.write_h5ad(args.output_amp, compression="gzip")
log_message("Done!")
