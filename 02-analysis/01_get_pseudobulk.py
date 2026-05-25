# Creation date: 2025-09-11
# Description: Generate pseudobulk profiles aggregated by donor x brain region x cell subtype.

import scanpy as sc
import pandas as pd
import decoupler as dc
import argparse
from sys import path
import os

path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "utils")))
from my_utils import log_message

# %% Parse CLI arguments
parser = argparse.ArgumentParser()
parser.add_argument("input_file", type=str)  # annotated h5ad
parser.add_argument("annotations_file", type=str)  # cell type annotation TSV
parser.add_argument("output_file", type=str)  # pseudobulk h5ad output
args = parser.parse_args()

# %% Load input files
log_message(f"Loading '{args.input_file}'...")
amp = sc.read_h5ad(args.input_file)

log_message(f"Loading '{args.annotations_file}'...")
annotations = pd.read_csv(args.annotations_file, sep="\t")

# %% Merge cell annotations
log_message("Merging cell annotation metadata...")
amp.obs = amp.obs.drop(columns=["cluster_id"])
amp.obs = amp.obs.astype({"leiden_0_25_corrected": "string"})
annotations = annotations.rename(columns={"cluster_id": "leiden_0_25_corrected"})
annotations = annotations.astype({"leiden_0_25_corrected": "string"})

original_index = amp.obs.index
amp.obs = amp.obs.merge(annotations, on="leiden_0_25_corrected", how="left")
amp.obs.index = original_index

# %% Remove doublets and low-quality clusters
log_message("Removing doublets and low-quality clusters (18 and 24)...")
amp = amp[~amp.obs["predicted_doublet"]]
amp = amp[~amp.obs["leiden_0_25_corrected"].isin(["18", "24"])]

# %% Aggregate cell-level QC metrics at participant level
log_message("Aggregating cell-level QC metrics at the participant level...")
columns_to_aggregate = [
    "n_genes",
    "total_counts",
    "pct_counts_in_top_20_genes",
    "doublet_score",
    "S.Score",
    "G2M.Score",
    "percent_apop",
    "percent_dna_repair",
    "percent_ieg",
    "percent_mito",
    "percent_mito_ribo",
    "percent_oxphos",
    "percent_ribo",
]
rename_dict = {key: f"{key}_median" for key in columns_to_aggregate}
aggregated_columns = (
    amp.obs.groupby("participant_id", observed=True)[columns_to_aggregate]
    .median()
    .rename(columns=rename_dict)
)

# %% Generate pseudobulk profiles
log_message("Computing pseudobulk profiles (donor x region x cell subtype)...")
pseudobulk = dc.pp.pseudobulk(
    amp,
    sample_col="participant_id",
    groups_col=["cell_subtype", "brain_region"],
    layer="counts",
    mode="sum",
    empty=True,
    verbose=True,
)

log_message("Removing empty pseudobulk profiles...")
dc.pp.filter_samples(pseudobulk, min_cells=1, min_counts=0, inplace=True)

# %% Attach participant-level aggregated metrics
log_message("Attaching participant-aggregated QC metrics...")
pseudobulk.obs = pseudobulk.obs.join(aggregated_columns, on="participant_id")

# %% Save
log_message(f"Saving pseudobulk h5ad to '{args.output_file}'...")
pseudobulk.write_h5ad(args.output_file, compression="gzip")
log_message("Done!")
