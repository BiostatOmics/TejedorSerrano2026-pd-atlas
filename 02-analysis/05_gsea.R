# Creation date: 2025-11-04
# Description: GSEA functional enrichment for each cell subtype x brain region contrast.
#   Uses MSigDB gene sets (excluding positional, cancer, oncogenic, and cell type collections).
#   Input: DEA results ranked by t-statistic. Output: combined GSEA results RDS.

# %% Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(logger)
  library(glue)
  library(org.Hs.eg.db)
  library(msigdbr)
  library(clusterProfiler)
  library(conflicted)
})
conflicts_prefer(dplyr::filter, dplyr::select, dplyr::rename)

# %% Parse CLI args
args <- commandArgs(trailingOnly = TRUE)
INPUT_FILE <- args[1] # DEA combined results RDS
OUTPUT_DIR <- args[2] # output directory
FINAL_GSEA_FILE <- args[3] # output filename

log_info("Arguments: {paste(args, collapse = ' | ')}")

# %% Load DEA results
log_info("Loading DEA results...")
de_results <- read_rds(INPUT_FILE)

# %% Load MSigDB gene sets
log_info("Loading MSigDB gene sets...")
MSIG <- msigdbr(species = "Homo sapiens")
MSIG_OF_INTEREST <- MSIG %>%
  filter(
    !gs_collection %>% str_to_lower() %in% c("c1", "c4", "c6", "c8"),
    !gs_collection_name %>% str_to_lower() %>% str_detect("_legacy"),
    !gs_subcollection %>% str_to_lower() %>% str_detect("vax"),
  )
MSIG_TERM2GENE <- MSIG_OF_INTEREST %>% select(gs_id, ensembl_gene)
MSIG_TERM2NAME <- MSIG_OF_INTEREST %>% select(gs_id, gs_name)
MSIG_METADATA <- MSIG_OF_INTEREST %>%
  distinct(
    gs_id,
    gs_collection,
    gs_subcollection,
    gs_collection_name,
    gs_description
  )

# %% Build ranked gene lists per contrast
log_info("Building ranked gene lists per contrast...")
extract_gene_list <- function(de_results_df, gene_id_type, stat = "t_stat") {
  de_results_df %>%
    dplyr::select(all_of(c(gene_id_type, stat))) %>%
    arrange(desc(.data[[stat]])) %>%
    filter(!is.na(.data[[gene_id_type]])) %>%
    distinct(.data[[gene_id_type]], .keep_all = TRUE) %>%
    deframe()
}

contrasts <- unique(de_results$contrast)
ranked_gene_list_per_contrast <- contrasts %>%
  set_names() %>%
  map(function(c) {
    extract_gene_list(de_results %>% filter(contrast == c), "ensembl_gene_id")
  })

# %% Run GSEA for each contrast
log_info("Running GSEA for {length(contrasts)} contrasts...")
gsea_all_contrasts <- ranked_gene_list_per_contrast %>%
  imap_dfr(function(gene_list, contrast) {
    log_info("Running GSEA for '{contrast}' ({length(gene_list)} genes)")
    GSEA(
      gene_list,
      minGSSize = 10,
      maxGSSize = 500,
      TERM2GENE = MSIG_TERM2GENE,
      TERM2NAME = MSIG_TERM2NAME,
      pvalueCutoff = 1,
      eps = 1e-50,
      pAdjustMethod = "fdr",
      verbose = FALSE,
      seed = TRUE
    ) %>%
      as_tibble() %>%
      mutate(contrast = contrast)
  }) %>%
  rename(gs_id = ID, gs_name = Description, gs_size = setSize) %>%
  left_join(MSIG_METADATA, by = "gs_id") %>%
  relocate(starts_with("gs"), .before = enrichmentScore) %>%
  mutate(
    cell_subtype_short = contrast %>%
      str_split_i("\\.", 1) %>%
      str_replace_all("__", "/") %>%
      str_replace_all("_", " "),
    brain_region = contrast %>% str_split_i("\\.", 2)
  )

# %% Save
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}
output_path <- file.path(OUTPUT_DIR, FINAL_GSEA_FILE)
log_info("Saving GSEA results to '{output_path}'...")
saveRDS(gsea_all_contrasts, file = output_path)
log_info("GSEA complete.")
