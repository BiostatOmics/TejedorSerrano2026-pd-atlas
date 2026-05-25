# Creation date: 2025-07-01
# Description: Differential cell type abundance analysis using tascCODA (pertpy).
#   Runs two models per brain region: raw (case/control only) and full (with covariates).
#   Also performs credibility analysis by rotating the reference cell type across all subtypes.

import pertpy as pt
from pertpy.tools._coda._base_coda import df2newick
import pandas as pd
import anndata as ad
import toytree as tt
import matplotlib.pyplot as plt
import arviz as az
import mudata as md
import pickle as pkl
from pathlib import Path
import argparse
import os

os.environ["QT_QPA_PLATFORM"] = "offscreen"  # prevent GUI crashes on HPC

import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "utils")))
from my_utils import log_message

# %% Parse CLI arguments
parser = argparse.ArgumentParser()
parser.add_argument("input_file", type=str, help="Cell metadata CSV")
parser.add_argument("annotations_file", type=str, help="Cell type annotations TSV")
parser.add_argument("n_burn", type=str, help="Number of NUTS warmup steps")
parser.add_argument("n_iter", type=str, help="Number of NUTS sampling steps")
parser.add_argument(
    "cred_iter_fraction",
    type=str,
    help="Fraction of iterations for credibility analysis",
)
parser.add_argument("phi", type=str, help="Spike-and-slab regularization (phi)")
parser.add_argument("fdr", type=str, help="FDR threshold for credible effects")
parser.add_argument(
    "reference_cell_type", type=str, help="Reference cell type (or 'automatic')"
)
parser.add_argument("full_formula", type=str, help="Patsy formula for the full model")
parser.add_argument("output_dir", type=str, help="Output directory")
args = parser.parse_args()

NUM_WARMUP = int(args.n_burn)
NUM_SAMPLES = int(args.n_iter)
CREDIBILITY_ITER_FRACTION = float(args.cred_iter_fraction)
PHI = int(args.phi)
FDR = float(args.fdr)
REFERENCE_CT = args.reference_cell_type
FULL_FORMULA = args.full_formula
output_dir = Path(args.output_dir)

# %% Load input files
log_message(f"Loading cell metadata from '{args.input_file}'...")
amp_obs = pd.read_csv(args.input_file)

log_message(f"Loading annotations from '{args.annotations_file}'...")
annotations = pd.read_csv(args.annotations_file, sep="\t")

output_dir.mkdir(parents=True, exist_ok=True)

# %% Clean inputs
log_message("Cleaning inputs...")
amp_obs = amp_obs.drop(
    columns=[
        "ROIGroupFine",
        "supercluster_term",
        "cell_type",
        "cluster_id",
        "leiden_0_10_ingested",
        "final_annotations_desc",
        "annotations",
    ]
)
amp_obs = amp_obs[~(amp_obs["predicted_doublet"])]
amp_obs = amp_obs[~(amp_obs["case_control"] == "Other")]
amp_obs = amp_obs[~(amp_obs["leiden_0_25_corrected"].isin([18, 24]))]

# Scale age per donor
donor_ages = amp_obs.drop_duplicates("participant_id").set_index("participant_id")[
    ["age_at_baseline"]
]
donor_ages["age_at_baseline_z_scaled"] = donor_ages["age_at_baseline"].transform(
    lambda x: (x - x.mean()) / x.std()
)
donor_ages["age_at_baseline_mm_scaled"] = donor_ages["age_at_baseline"].transform(
    lambda x: (x - x.min()) / (x.max() - x.min())
)
amp_obs = amp_obs.merge(
    donor_ages[["age_at_baseline_mm_scaled", "age_at_baseline_z_scaled"]],
    on="participant_id",
)

# Clean annotations
unnamed_cols = annotations.columns[annotations.columns.str.startswith("Unn")]
annotations = annotations.drop(columns=unnamed_cols)
annotations = annotations.rename(columns={"cluster_id": "leiden_0_25_corrected"})
amp_obs = amp_obs.merge(annotations, on="leiden_0_25_corrected")

# %% Filter rare cell subtype x region pairs (< 5% of subtype total)
log_message("Filtering rare cell subtype x region combinations...")
MIN_PROP = 0.05
amp_obs["count_subtype_in_region"] = amp_obs.groupby(
    ["cell_subtype_short", "brain_region"]
)["participant_id"].transform("size")
amp_obs["count_subtype_total"] = amp_obs.groupby(["cell_subtype_short"])[
    "participant_id"
].transform("size")
amp_obs["prop_subtype_in_region"] = (
    amp_obs["count_subtype_in_region"] / amp_obs["count_subtype_total"]
)
amp_obs = amp_obs[amp_obs["prop_subtype_in_region"] >= MIN_PROP].copy()

# %% Build per-region CoDA objects
log_message("Building per-region CoDA objects...")
coda = pt.tl.Tasccoda()
amp_coda = {}
regions = amp_obs["brain_region"].unique().tolist()

for r in regions + ["ALL"]:
    region_subset = amp_obs if r == "ALL" else amp_obs[amp_obs["brain_region"] == r]

    counts = pd.crosstab(
        [
            region_subset["case_control"],
            region_subset["cohort"],
            region_subset["participant_id"],
            region_subset["sample_id"],
            region_subset["age_at_baseline_mm_scaled"],
            region_subset["sex"],
            region_subset["brain_region"],
        ],
        region_subset["cell_subtype_short"],
    )
    obs = pd.DataFrame(index=counts.index).reset_index()
    var = region_subset.groupby("cell_subtype_short")[
        ["cell_type_short", "cell_supertype_short"]
    ].first()

    tree_levels = ["cell_supertype_short", "cell_type_short", "cell_subtype_short"]
    newick = df2newick(region_subset[tree_levels].reset_index(drop=True), tree_levels)
    tree = tt.tree(newick)

    region_adata = ad.AnnData(
        X=counts.reset_index(drop=True),
        var=var,
        obs=obs,
        uns={"newick": newick, "toytree": tree},
    )
    region_coda = coda.load(
        region_adata,
        type="sample_level",
        levels_agg=tree_levels,
        key_added="tree",
        add_level_name=True,
    )
    del region_coda.mod["rna"], region_coda.obsm["rna"], region_coda.varm["rna"]
    amp_coda[r] = region_coda

# %% Run tascCODA models
log_message("Running tascCODA models...")
amp_coda_processed_raw = {}
amp_coda_processed_full = {}
formulas = [
    (amp_coda_processed_raw, "raw", "C(case_control, Treatment('Control'))"),
    (amp_coda_processed_full, "full", FULL_FORMULA),
]

all_arviz = {}
for r, original in amp_coda.items():
    arviz_dict = {}
    for processed_dict, suffix, formula in formulas:
        log_message(f"Running {suffix}-model for region: {r}")
        c = original.copy()
        c = coda.prepare(
            c,
            modality_key="coda",
            tree_key="tree",
            reference_cell_type=REFERENCE_CT,
            formula=formula,
            pen_args={"phi": PHI},
        )
        coda.run_nuts(
            c, modality_key="coda", num_warmup=NUM_WARMUP, num_samples=NUM_SAMPLES
        )
        processed_dict[r] = c

        intercept_df, effect_df, node_df = coda.summary_prepare(c["coda"], est_fdr=FDR)

        # Export arviz diagnostics
        arviz = coda.make_arviz(c)
        arviz.to_netcdf(output_dir / f"arviz-{r}-{suffix}.nc")
        arviz_dict[suffix] = arviz

        # Export tree
        c["coda"].uns["tree"].write(output_dir / f"tree-{r}-{suffix}.newick")

        # Export results to CSV
        intercept_df.to_csv(output_dir / f"intercept-{r}-{suffix}.csv")
        effect_df.to_csv(output_dir / f"effects-{r}-{suffix}.csv")
        node_df.to_csv(output_dir / f"node-{r}-{suffix}.csv")
        pd.DataFrame(
            [
                [
                    c["coda"].uns["scCODA_params"]["reference_cell_type"],
                    c["coda"].uns["scCODA_params"]["formula"],
                    c["coda"].X.shape,
                ]
            ],
            index=[r],
            columns=["reference_cell_type", "formula", "samples x cell types"],
        ).to_csv(output_dir / f"model-{r}-{suffix}.csv")

        # Export model params
        params = c["coda"].uns["scCODA_params"]
        with open(output_dir / f"tascCODA-params-{r}-{suffix}.pkl", "wb") as f:
            pkl.dump(params, f)

        save_compatible = c.copy()
        del save_compatible["coda"].uns["tree"]
        del save_compatible["coda"].uns["toytree"]
        del save_compatible["coda"].uns["scCODA_params"]
        md.write_h5mu(
            output_dir / f"tascCODA-object-{r}-{suffix}.h5md", save_compatible
        )

    all_arviz[r] = arviz_dict

# %% Credibility analysis: rotate reference cell type
log_message("Running credibility analysis (rotating reference cell type)...")
for r, model in amp_coda.items():
    credible_dict = {}
    cell_types = model["coda"].var.index
    for ct in cell_types:
        log_message(f"Region: {r} | Reference: {ct}")
        model_copy = model.copy()
        model_copy = coda.prepare(
            model_copy,
            modality_key="coda",
            tree_key="tree",
            reference_cell_type=ct,
            formula=FULL_FORMULA,
            pen_args={"phi": PHI},
        )
        coda.run_nuts(
            model_copy,
            modality_key="coda",
            num_warmup=int(NUM_WARMUP * CREDIBILITY_ITER_FRACTION),
            num_samples=int(NUM_SAMPLES * CREDIBILITY_ITER_FRACTION),
        )
        credible = coda.get_node_df(model_copy, modality_key="coda")["Final Parameter"]
        credible = credible.reset_index()
        credible["reference"] = ct
        credible_dict[f"reference_{ct}"] = credible.set_index(["Covariate", "Node"])[
            "Final Parameter"
        ]

    credible_df = pd.concat(credible_dict, axis=1)
    ref_cols = credible_df.columns
    credible_df["times_credible"] = (credible_df[ref_cols] != 0).sum(axis=1)
    credible_df["n_tests"] = len(ref_cols)
    credible_df["pct_credible"] = credible_df["times_credible"] / credible_df["n_tests"]
    credible_df["is_credible"] = credible_df["pct_credible"] >= 0.5
    credible_df["region"] = r
    credible_df.reset_index().to_csv(
        output_dir / f"all-references-credibility-{r}-full.csv", index=False
    )

log_message("Differential abundance analysis complete!")
