# Creation date: 2025-09-25
# Description: Apply limma DEA between PD cases and controls independently for each cell subtype.
#   Per cell type: gene filtering, TMM normalization, voomLmFit with donor blocking.
#   Covariates: age, PMI, library size, sex, cohort.
# Reference: https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf

# %% Load libraries
suppressPackageStartupMessages({
  library(conflicted)
  library(tidyverse)
  library(edgeR)
  library(limma)
  library(futile.logger)
  library(glue)
})

# %% Parse CLI args
args <- commandArgs(trailingOnly = TRUE)
PREPROCESSED_DGE_FILE <- args[1] # preprocessed DGEList RDS
OUTPUT_DIR <- args[2] # output directory
FINAL_DGE_LIST_FILE <- args[3] # filename for filtered/normalized DGE list RDS
VOOMFIT_LIST_FILE <- args[4] # filename for voomLmFit list RDS
MIN_COUNT <- as.numeric(args[5]) # minimum count for filterByExpr (default: 1)
MIN_TOTAL_COUNT <- as.numeric(args[6]) # minimum total count (default: 50)

flog.info("Arguments: %s", paste(args, collapse = " | "))

# %% Load input
flog.info("Loading preprocessed DGEList...")
dge <- read_rds(PREPROCESSED_DGE_FILE)

# %% Split DGE by cell type
flog.info("Splitting DGE by cell subtype...")
dge_per_cell_type <- unique(dge$samples$cell_subtype_short) %>%
  set_names() %>%
  map(~ dge[, dge$samples$cell_subtype_short == .x])

# %% Gene filtering per cell type
flog.info("Filtering lowly expressed genes per cell type...")
dge_per_cell_type_gene_filtered <- dge_per_cell_type %>%
  imap(function(dge, ct) {
    genes_to_keep <- filterByExpr(
      dge,
      group = dge$samples$group,
      min.total.count = MIN_TOTAL_COUNT,
      min.count = MIN_COUNT,
      min.prop = 0.85
    )
    flog.info(
      "%s: keeping %d / %d genes",
      ct,
      sum(genes_to_keep),
      length(genes_to_keep)
    )
    dge[genes_to_keep, ]
  })

# %% TMM normalization per cell type
flog.info("Computing TMM normalization factors...")
dge_per_cell_type_normalized <- dge_per_cell_type_gene_filtered %>%
  map(~ calcNormFactors(.x, method = "TMMwsp"))

# %% Save filtered/normalized DGE list
flog.info(
  "Saving filtered DGE list to: %s",
  file.path(OUTPUT_DIR, FINAL_DGE_LIST_FILE)
)
saveRDS(
  dge_per_cell_type_normalized,
  file = file.path(OUTPUT_DIR, FINAL_DGE_LIST_FILE)
)

# %% Fit voomLmFit per cell type
flog.info("Fitting voomLmFit models...")
voomfit_per_cell_type <- dge_per_cell_type_normalized %>%
  imap(function(dge, ct) {
    design <- model.matrix(
      ~ 0 + group + age_scaled + pmi_scaled + log_libsize_scaled + sex + cohort,
      data = dge$samples
    )
    colnames(design) <- colnames(design) %>% str_remove_all("group")

    flog.info("Fitting '%s': %d genes x %d samples", ct, nrow(dge), ncol(dge))

    voomLmFit(
      dge,
      design = design,
      sample.weights = TRUE,
      block = dge$samples$participant_id,
      save.plot = TRUE
    )
  })

# %% Save voomfit list
flog.info(
  "Saving voomfit list to: %s",
  file.path(OUTPUT_DIR, VOOMFIT_LIST_FILE)
)
saveRDS(voomfit_per_cell_type, file = file.path(OUTPUT_DIR, VOOMFIT_LIST_FILE))
flog.info("DEA complete.")
