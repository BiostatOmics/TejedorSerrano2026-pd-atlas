# Creation date: 2025-09-25
# Description: Post-DEA processing: build contrasts, apply eBayes, extract topTable results.
#   Builds region-level (cell subtype x region) and cell-type-level (global) contrasts.
#   Outputs a combined DEA results dataframe annotated with gene metadata.

# %% Load libraries
suppressPackageStartupMessages({
  library(conflicted)
  library(tidyverse)
  library(limma)
  library(futile.logger)
})
conflicts_prefer(dplyr::filter)

# %% Input paths (relative to project root)
VOOMFIT_LIST_FILE <- "04-analysis/03-02-dea-disease/pseudobulk-dea-voomLmFit-per-cell-type.rds"
DGE_LIST_FILE <- "04-analysis/03-02-dea-disease/pseudobulk-dea-final-dge-per-cell-type.rds"
GENE_METADATA_FILE <- "03-data/02-metadata/gene-metadata.csv"
OUTPUT_DIR <- "04-analysis/03-02-dea-disease/"

flog.info("Loading input files...")
voomfit_per_cell_type <- readRDS(VOOMFIT_LIST_FILE)
dge_per_cell_type <- readRDS(DGE_LIST_FILE)
gene_metadata <- read_csv(
  GENE_METADATA_FILE,
  col_types = cols(entrezgene_id = "c")
) %>%
  filter(transcript_is_canonical == 1) %>%
  select(
    ensembl_gene_id,
    hgnc_symbol,
    entrezgene_id,
    gene_biotype,
    chromosome_name,
    percentage_gene_gc_content
  ) %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

# %% Build contrast matrices per cell type
flog.info("Building contrasts...")
contrast_voomfit_per_cell_type <- voomfit_per_cell_type %>%
  map(function(fit) {
    # Identify cell_type.region groups with both Case and Control samples
    celltype_region_group <- colnames(fit$design) %>%
      str_subset("\\.Case$|\\.Control$") %>%
      str_remove("\\.Case$|\\.Control$")

    celltype_group <- celltype_region_group %>%
      str_split_i("\\.", i = 1) %>%
      unique()

    groups_to_keep <- !(table(celltype_region_group) == 1)
    celltype_region_group <- names(table(celltype_region_group))[groups_to_keep]

    # Region-level contrasts: Case vs Control for each cell subtype x region
    contrasts_list <- map(
      celltype_region_group,
      ~ paste0(.x, ".Case - ", .x, ".Control")
    )
    names(contrasts_list) <- paste0(
      celltype_region_group,
      ".CasevsControl.region_level"
    )

    # Cell-type-level (global) contrast: average across regions
    if (length(celltype_region_group) > 1) {
      formula <- map_chr(
        celltype_region_group,
        ~ paste0(.x, ".Case - ", .x, ".Control")
      ) %>%
        str_flatten(collapse = " + ")
      formula <- str_c("(", formula, ")/", length(celltype_region_group))
      contrasts_list[[paste0(
        celltype_group,
        ".GLOBAL.CasevsControl.celltype_level"
      )]] <- formula
    }

    contrasts_list$levels <- fit$design
    contrasts_matrix <- do.call(makeContrasts, contrasts_list)
    contrasts.fit(fit, contrasts_matrix)
  })

# %% Apply eBayes
flog.info("Applying eBayes...")
bayes_voomfit_per_cell_type <- contrast_voomfit_per_cell_type %>%
  map(~ eBayes(.x, robust = TRUE))

# %% Extract topTable results
flog.info("Extracting topTable results...")

extract_subtype <- function(contrast) {
  contrast %>%
    str_split_i("\\.", i = 1) %>%
    str_replace_all("__", "/") %>%
    str_replace_all("_", " ")
}
extract_region <- function(contrast) contrast %>% str_split_i("\\.", i = 2)

de_results_df_combined <- bayes_voomfit_per_cell_type %>%
  map(function(fit) {
    colnames(fit$contrasts) %>%
      map_dfr(function(contrast) {
        topTable(
          fit,
          coef = contrast,
          number = Inf,
          adjust.method = "fdr",
          confint = TRUE
        ) %>%
          rownames_to_column("ensembl_gene_id") %>%
          rename(
            logfc = logFC,
            ci_lower = CI.L,
            ci_upper = CI.R,
            avg_expr = AveExpr,
            t_stat = t,
            p_value = P.Value,
            fdr = adj.P.Val,
            b_stat = B
          ) %>%
          mutate(
            contrast = contrast,
            contrast_type = contrast %>%
              str_extract("(region_level|celltype_level)$"),
            cell_subtype_short = extract_subtype(contrast),
            brain_region = extract_region(contrast)
          )
      })
  }) %>%
  bind_rows() %>%
  left_join(gene_metadata, by = join_by(ensembl_gene_id)) %>%
  relocate(colnames(gene_metadata), .before = 1) %>%
  as_tibble()

# %% Save results
flog.info("Saving combined DEA results...")
saveRDS(
  de_results_df_combined %>% filter(contrast_type == "region_level"),
  file = file.path(
    OUTPUT_DIR,
    "pseudobulk-dea-combined-results-ebayes-region-level.rds"
  )
)
saveRDS(
  de_results_df_combined,
  file = file.path(
    OUTPUT_DIR,
    "pseudobulk-dea-combined-results-ebayes-fc-1_00.rds"
  )
)
flog.info("Downstream DEA processing complete.")
