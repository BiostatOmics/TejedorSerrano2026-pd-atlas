# Creation date: 2025-09-25
# Description: Preprocess pseudobulk matrix before limma DEA.
#   - Loads pseudobulk h5ad as DGEList
#   - Removes "Other" disease donors and donors with missing PMI
#   - Removes rare cell type x region combinations
#   - Filters pseudobulk profiles with fewer than N_CELLS_THRESHOLD cells
#   - Scales covariates for downstream modeling

# %% Load libraries
suppressPackageStartupMessages({
  library(zellkonverter)
  library(SingleCellExperiment)
  library(conflicted)
  library(tidyverse)
  library(edgeR)
  library(futile.logger)
})

# %% Parse CLI args
args <- commandArgs(trailingOnly = TRUE)
INPUT_FILE <- args[1] # pseudobulk h5ad
CLINICAL_METADATA <- args[2] # CSV with PMI per donor
OUTPUT_DGE_FILE <- args[3] # output RDS path
N_CELLS_THRESHOLD <- as.numeric(args[4]) # minimum cells per pseudobulk profile (default: 10)

flog.info("Arguments: %s", paste(args, collapse = " | "))

# %% Load input files
flog.info("Loading input files...")
sce <- readH5AD(INPUT_FILE, reader = "R")
metadata <- read_csv(CLINICAL_METADATA, col_types = "f")

# %% Clean metadata
flog.info("Cleaning metadata...")
names(assays(sce))[1] <- "counts"

clean_metadata <- colData(sce) %>%
  as.data.frame() %>%
  rownames_to_column("psbulk_profile") %>%
  select(!starts_with("Unnamed..")) %>%
  select(!c(orig.ident, predicted_doublet)) %>%
  mutate(across(where(is.factor), as.character)) %>%
  type_convert() %>%
  left_join(metadata %>% select(participant_id, PMI), by = "participant_id")

# %% Build DGEList
flog.info("Building DGEList...")
dge <- DGEList(counts = counts(sce))

dge$samples <- dge$samples %>%
  rownames_to_column("psbulk_profile") %>%
  left_join(clean_metadata, by = "psbulk_profile") %>%
  column_to_rownames("psbulk_profile") %>%
  mutate(
    rowname = colnames(dge$counts),
    group = factor(paste(
      cell_subtype_short %>%
        str_replace_all(" ", "_") %>%
        str_replace_all("/", "__"),
      brain_region,
      case_control,
      sep = "."
    ))
  )

# %% Remove "Other" disease donors and donors with missing PMI
flog.info("Removing 'Other' donors and missing PMI...")
dge <- dge[, dge$samples$case_control != "Other"]
dge <- dge[, !is.na(dge$samples$PMI)]

# %% Remove rare cell type x region combinations (< 5% of total cells for that subtype)
dge$samples <- dge$samples %>%
  left_join(
    dge$samples %>%
      group_by(cell_subtype_short, brain_region) %>%
      summarise(n_cells = sum(psbulk_cells), .groups = "drop_last") %>%
      mutate(prop = n_cells / sum(n_cells), rare_ct_reg = prop < 0.05) %>%
      select(-c(n_cells, prop))
  )
dge <- dge[, !dge$samples$rare_ct_reg]

# %% Filter by minimum cell number
flog.info(
  "Filtering pseudobulk profiles by minimum cell count (%d)...",
  N_CELLS_THRESHOLD
)
dge_preprocessed <- dge[, dge$samples$psbulk_cells >= N_CELLS_THRESHOLD]

# %% Scale covariates
dge_preprocessed$samples <- dge_preprocessed$samples %>%
  mutate(
    age_scaled = as.numeric(scale(age_at_baseline)),
    pmi_scaled = as.numeric(scale(PMI)),
    sex = factor(sex),
    cohort = factor(cohort),
    participant_id = factor(participant_id),
    log_libsize_scaled = as.numeric(scale(log1p(lib.size)))
  )

# %% Save
flog.info("Saving preprocessed DGEList to: %s", OUTPUT_DGE_FILE)
saveRDS(dge_preprocessed, file = OUTPUT_DGE_FILE)
flog.info("Preprocessing complete.")
