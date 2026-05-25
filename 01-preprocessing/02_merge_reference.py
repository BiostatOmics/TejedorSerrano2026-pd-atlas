"""
Merge neuronal and non-neuronal reference h5ad datasets into a single object.

Usage:
    python 02_merge_reference.py <input_h5ad_1> <input_h5ad_2> <output_h5ad>
"""

import argparse
import anndata as ad
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "utils")))
from my_utils import log_message

# CLI
parser = argparse.ArgumentParser(description="Merge two reference h5ad datasets")
parser.add_argument("input_h5ad_1", type=str, help="Path to neuronal h5ad")
parser.add_argument("input_h5ad_2", type=str, help="Path to non-neuronal h5ad")
parser.add_argument("output_h5ad", type=str, help="Output path for merged h5ad")
args = parser.parse_args()

# MAIN

log_message(f"Loading h5ad_1 from {args.input_h5ad_1}...")
h5ad_1 = ad.read_h5ad(args.input_h5ad_1)

log_message(f"Loading h5ad_2 from {args.input_h5ad_2}...")
h5ad_2 = ad.read_h5ad(args.input_h5ad_2)

log_message("Merging h5ad_1 and h5ad_2...")
reference = ad.concat(
    [h5ad_1, h5ad_2],
    label="dataset_title",
    keys=["neuronal", "non_neuronal"],
    join="outer",
    merge="unique",
    uns_merge="unique",
)

# Store per-dataset metadata that could not be auto-merged
log_message("Storing per-dataset metadata in uns...")
reference.uns["citation_neurons"] = h5ad_1.uns["citation"]
reference.uns["title_neurons"] = h5ad_1.uns["title"]
reference.uns["citation_non_neurons"] = h5ad_2.uns["citation"]
reference.uns["title_non_neurons"] = h5ad_2.uns["title"]
reference.uns["X_UMAP_neurons"] = h5ad_1.obsm["X_UMAP"]
reference.uns["X_tSNE_neurons"] = h5ad_1.obsm["X_tSNE"]
reference.uns["X_UMAP_non_neurons"] = h5ad_2.obsm["X_UMAP"]
reference.uns["X_tSNE_non_neurons"] = h5ad_2.obsm["X_tSNE"]

log_message(f"Saving merged h5ad to {args.output_h5ad}...")
reference.write_h5ad(args.output_h5ad, compression="gzip")
log_message("Done!")
