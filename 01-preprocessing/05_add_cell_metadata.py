"""
Transfer cell-level metadata from an external CSV to the AMP-PD h5ad object.

Usage:
    python 04_add_cell_metadata.py <input_amp> <cell_metadata_path> <output_amp>
"""

import argparse
import pandas as pd
import anndata as ad
import scanpy as sc
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "utils")))
from my_utils import log_message

# CLI
parser = argparse.ArgumentParser(description="Transfer cell metadata to AMP h5ad")
parser.add_argument("input_amp", type=str, help="Path to input h5ad")
parser.add_argument("cell_metadata_path", type=str, help="Path to cell metadata CSV")
parser.add_argument("output_amp", type=str, help="Output path for updated h5ad")
args = parser.parse_args()

sc.settings.verbosity = 4
sc.settings.n_jobs = 8

# LOAD DATA
log_message(f"Loading AMP dataset from {args.input_amp}...")
amp = ad.read_h5ad(args.input_amp)

log_message(f"Loading cell metadata from {args.cell_metadata_path}...")
cell_metadata = pd.read_csv(args.cell_metadata_path, sep=",", header=0, index_col=0)

# METADATA TRANSFER
log_message("Transferring cell metadata...")
# Identify columns not already present in the AnnData object
columns_to_transfer = cell_metadata.columns.difference(amp.obs.columns)
cell_metadata_to_transfer = cell_metadata[columns_to_transfer]
# Align cells between datasets
cell_metadata_to_transfer = cell_metadata_to_transfer[
    cell_metadata_to_transfer.index.isin(amp.obs.index)
]
amp.obs = pd.concat([amp.obs, cell_metadata_to_transfer], axis=1)

# Convert ordinal/staging columns to categorical type
for col in ["hoehn_yahr_stage", "path_braak_lb", "path_braak_nft", "barcodekey"]:
    amp.obs[col] = amp.obs[col].astype("category")

log_message("Metadata transfer complete.")

# SAVE
log_message(f"Saving updated object to {args.output_amp}...")
amp.write_h5ad(args.output_amp, compression="gzip")
log_message("Done!")
