"""
Apply manual cell type annotations based on Leiden cluster assignments.

Assigns cell type labels from a predefined cluster-to-annotation mapping,
extracts abbreviations from full names, and removes low-quality clusters.

Usage:
    python 09_manual_annotation.py <input_amp> <output_amp> <output_tsv>
"""

import argparse
import pandas as pd
import anndata as ad
import scanpy as sc
import re
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "utils")))
from my_utils import log_message

# CLI
parser = argparse.ArgumentParser(description="Manual cell type annotation")
parser.add_argument("input_amp", type=str, help="Path to input h5ad")
parser.add_argument("output_amp", type=str, help="Output path for annotated h5ad")
parser.add_argument("output_tsv", type=str, help="Output path for annotation TSV")
args = parser.parse_args()

sc.settings.verbosity = 4
sc.settings.n_jobs = 8

# LOAD DATA
log_message(f"Loading AMP dataset from {args.input_amp}...")
amp = ad.read_h5ad(args.input_amp)

# PREPARE OBJECT
# Convert ordinal columns to categorical
for col in ["hoehn_yahr_stage", "path_braak_lb", "path_braak_nft"]:
    amp.obs[col] = amp.obs[col].astype("category")

# Remove CPM layer to reduce file size
del amp.layers["cpm"]

# CLUSTER-TO-ANNOTATION MAPPING
# Annotations derived from Leiden clustering at resolution 0.10 and 0.25
# Format: cluster_id -> "Full Name (Abbreviation)"

mapping_0_10 = {
    0: "Oligodendrocytes (OLGs)",
    2: "Astrocytes (ASTs)",
    3: "Microglia (MG)",
    4: "Oligodendrocyte Precursor Cells (OPCs)",
    6: "MGE-derived Interneurons (MGE-derived Interneurons)",
    7: "Deep-layer Intratelencephalic PCP4+ Neurons (DL-ITC PCP4+ Neurons)",
    10: "Deep-layer Corticothalamic and Layer 6b TLE4+ Neurons (DL-CTC-6b TLE4+ Neurons)",
    11: "Endothelial Cells (ECs)",
    12: "Deep-layer Near-projecting CRYM+ Neurons (DL-NPC CRYM+ Neurons)",
    13: "Chandelier Cells (ChCs)",
    14: "Upper-layer Intratelencephalic CUX2- Neurons (UL-ITC CUX2- Neurons)",
}

mapping_0_25 = {
    2: "Upper-layer Intratelencephalic CUX2+ Neurons (UL-ITC CUX2+ Neurons)",
    7: "Deep-layer Intratelencephalic PCP4- Neurons (DL-ITC PCP4- Neurons)",
    15: "LAMP5+ Interneurons (LAMP5+ Interneurons)",
    10: "CGE-derived Interneurons (CGE-derived Interneurons)",
    14: "Neurovascular Support Cells (NVSCs)",
    17: "Cerebrospinal Fluid-related Cells of the Dorsal Motor Nucleus of the Vagus (CSFCs-DMNX)",
    18: "**REMOVE**",
    16: "DRD2+ Neurons (DRD2+ Neurons)",
}

# ASSIGN ANNOTATIONS
df = amp.obs[["leiden_0_10_corrected", "leiden_0_25_corrected"]].copy()
df = df.apply(pd.to_numeric)

df_0_10 = df[["leiden_0_10_corrected"]].copy()
df_0_10["annotations"] = df_0_10["leiden_0_10_corrected"].map(mapping_0_10)

df_0_25 = df[["leiden_0_25_corrected"]].copy()
df_0_25["annotations"] = df_0_25["leiden_0_25_corrected"].map(mapping_0_25)

# Merge and resolve: prefer resolution-0.25 annotations when available
df_combined = df_0_10.join(df_0_25, lsuffix="_0_10", rsuffix="_0_25")

def combine_annotations(row):
    for col in ["annotations_0_10", "annotations_0_25"]:
        if col in row and row[col] != "Unknown" and pd.notna(row[col]):
            return row[col]
    return "Unknown"

df_combined["final_annotations_desc"] = df_combined.apply(combine_annotations, axis=1)

# Extract abbreviation from parentheses: "Full Name (Abbrev)" -> "Abbrev"
def extract_abbreviation(annotation):
    match = re.search(r"\(([^)]+)\)", annotation)
    return match.group(1) if match else annotation

df_combined["annotations"] = df_combined["final_annotations_desc"].apply(extract_abbreviation)

# Transfer annotations to AnnData object
amp.obs["final_annotations_desc"] = df_combined["final_annotations_desc"]
amp.obs["annotations"] = df_combined["annotations"]

# Remove low-quality clusters flagged during annotation
log_message("Removing low-quality clusters flagged during annotation...")
amp = amp[amp.obs["annotations"] != "**REMOVE**"]

# SAVE
log_message(f"Saving annotated object to {args.output_amp}...")
amp.write_h5ad(args.output_amp, compression="gzip")

log_message(f"Saving annotation table to {args.output_tsv}...")
df_final = amp.obs[[
    "leiden_0_10_corrected", "leiden_0_25_corrected", "leiden_1_00_corrected",
    "final_annotations_desc", "annotations", "supercluster_term"
]]
df_final.to_csv(args.output_tsv, sep="\t", index_label="cell_id")
log_message("Done!")
