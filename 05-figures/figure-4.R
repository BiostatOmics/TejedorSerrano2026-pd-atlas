# Description: Generate Figure 4 panels:
#   panel a — logFC Spearman correlation heatmap
#   panels b, c — PCA of logFCs (PC1 vs PC2)
#   panel d — PCA of logFCs highlighting DMNX (PC5 vs PC9)
#   panel e — DMNX-axis gene loadings barplot
#   panel f — DMNX-axis GSEA dotplot

# %% Setup
library(tidyverse)
library(patchwork)
library(glue)
library(clusterProfiler)
library(msigdbr)
library(ComplexHeatmap)
library(circlize)
library(vegan)
source("utils/my_utils.R")


# %% Load data
COMBINED_DEA_RESULTS <- "04-analysis/private/pseudobulk-dea-combined-results-ebayes-fc-1_00.rds"
ANNOTATIONS <- "03-data/public/annotations_7_11_agg.txt"

dea <- readRDS(COMBINED_DEA_RESULTS)
dea_wo_global <- dea %>%
  filter(brain_region != "GLOBAL") %>%
  mutate(symbol = if_else(is.na(hgnc_symbol), ensembl_gene_id, hgnc_symbol))

annotations <- read_tsv(ANNOTATIONS) %>%
  select(
    cell_subtype_short,
    cell_subtype_color_vibrant,
    cell_type_short,
    cell_type_color_vibrant,
    cell_supertype_short,
    cell_supertype_color
  ) %>%
  distinct(cell_subtype_short, .keep_all = TRUE)


# %% Figure 4 panels b-e: PCA of logFC matrix across contexts

# Build gene x context logFC matrix; fill absent context-gene pairs with 0
logfc_matrix <- dea_wo_global %>%
  mutate(context = glue("{cell_subtype_short} - {brain_region}")) %>%
  select(ensembl_gene_id, context, logfc) %>%
  pivot_wider(names_from = context, values_from = logfc, values_fill = 0) %>%
  column_to_rownames("ensembl_gene_id") %>%
  scale() %>%
  t()

pca <- prcomp(logfc_matrix)
var_explained <- summary(pca)$importance[2, ]

pca_df <- as_tibble(pca$x, rownames = "context") %>%
  separate(
    context,
    into = c("cell_subtype_short", "brain_region"),
    sep = " - ",
    remove = FALSE
  )


# Shared layer definitions for PCA plots
common_theme <- list(theme_publication(5), no_legend)

layers_ellipse <- list(
  stat_ellipse(aes(color = cell_type_short), level = .8, linewidth = .1),
  scale_color_manual(values = correct_cell_type_colors)
)

layers_points <- list(
  geom_point(aes(shape = brain_region, color = cell_subtype_short), alpha = .9),
  scale_color_manual(values = correct_cell_subtype_colors),
  scale_shape_manual(values = region_shapes)
)

layers_p2 <- list(
  geom_point(aes(shape = brain_region, color = brain_region), alpha = .9),
  stat_ellipse(aes(color = brain_region), level = .9, linewidth = .1),
  scale_color_manual(values = region_colors),
  scale_shape_manual(values = region_shapes),
  scale_x_continuous(breaks = scales::breaks_pretty(3)),
  scale_y_continuous(breaks = scales::breaks_pretty(3))
)

# %% Figure 4 panels b, c: PC1 vs PC2
p1 <- pca_df %>%
  left_join(annotations) %>%
  ggplot(aes(x = PC1, y = PC2)) +
  layers_ellipse +
  ggnewscale::new_scale_color() +
  layers_points +
  labs(
    x = paste0("PC1 (", round(var_explained[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2] * 100, 1), "%)")
  ) +
  common_theme

p2 <- pca_df %>%
  left_join(annotations) %>%
  filter(cell_type_short %in% c("Exc", "Inh")) %>%
  ggplot(aes(x = PC1, y = PC2)) +
  layers_p2 +
  common_theme

p_pca_subtypes <- p1 / (plot_spacer() + p2) + plot_layout(heights = c(2, 1))
ggsave_pdf_svg(
  base_name = "05-results/02-plots/figure-3-pca-logfc-subtypes",
  plot = p_pca_subtypes,
  width = 45,
  height = 60,
  units = "mm"
)

# %% Figure 4 panel d: PC5 vs PC9 — DMNX regional signature
p1 <- pca_df %>%
  left_join(annotations) %>%
  ggplot(aes(x = PC5, y = PC9)) +
  stat_ellipse(aes(color = brain_region), level = .8, linewidth = .1) +
  scale_color_manual(values = region_colors) +
  ggnewscale::new_scale_color() +
  geom_vline(xintercept = 0, linewidth = .1, linetype = "dotted") +
  geom_hline(yintercept = 0, linewidth = .1, linetype = "dotted") +
  geom_point(aes(shape = brain_region, color = brain_region), alpha = .9) +
  scale_color_manual(values = region_colors) +
  scale_shape_manual(values = region_shapes) +
  labs(
    x = paste0("PC5 (", round(var_explained[5] * 100, 1), "%)"),
    y = paste0("PC9 (", round(var_explained[9] * 100, 1), "%)")
  ) +
  common_theme

p2 <- pca_df %>%
  left_join(annotations) %>%
  filter(brain_region == "DMNX") %>%
  ggplot(aes(x = PC5, y = PC9)) +
  layers_p2 +
  common_theme

p_pca_dmnx <- p1 / (plot_spacer() + p2) + plot_layout(heights = c(2, 1))
ggsave_pdf_svg(
  base_name = "05-results/02-plots/figure-3-pca-dmnx-logfc-subtypes",
  plot = p_pca_dmnx,
  width = 45,
  height = 60,
  units = "mm"
)


# %% Figure 4 panel e: DMNX-axis gene loadings barplot
pc_cols <- grep("^PC", colnames(pca_df), value = TRUE)
pc_cor_test <- function(pc, region_vector) {
  test <- cor.test(pca_df[[pc]], as.numeric(region_vector == "DMNX"))
  c(cor = test$estimate, p = test$p.value)
}
cor_results <- sapply(pc_cols, pc_cor_test, region_vector = pca_df$brain_region)

pc_dmnx <- cor_results %>%
  t() %>%
  as_tibble(rownames = "PC") %>%
  filter(p < 0.01) %>%
  arrange(desc(abs(cor.cor))) %>%
  pull(PC)

gene_map <- dea_wo_global %>% distinct(ensembl_gene_id, symbol)

loadings <- as_tibble(pca$rotation, rownames = "ensembl_gene_id") %>%
  left_join(gene_map, by = "ensembl_gene_id") %>%
  mutate(
    regional_PD_loading = rowMeans(across(all_of(pc_dmnx)), na.rm = TRUE)
  ) %>%
  arrange(desc(abs(regional_PD_loading)))

top_n_genes <- 30
top_genes <- loadings %>%
  slice_max(regional_PD_loading, n = top_n_genes) %>%
  pull(symbol)

dmnx_loadings_plot <- ggplot(
  loadings %>% filter(symbol %in% top_genes),
  aes(y = fct_reorder(symbol, regional_PD_loading), x = regional_PD_loading)
) +
  geom_col() +
  theme_publication(5) +
  labs(x = "Contribution\nto DMNX-axis") +
  theme(
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(face = "italic")
  )

ggsave_pdf_svg(
  base_name = "05-results/02-plots/figure-3-pca-dmnx-signature-loadings",
  width = 30,
  height = 60,
  units = "mm",
  plot = dmnx_loadings_plot
)


# %% Figure 4 panel f: DMNX-axis GSEA dotplot
MSIG <- msigdbr(species = "Homo sapiens")
MSIG_OF_INTEREST <- MSIG %>%
  filter(
    !gs_collection %>% str_to_lower() %in% c("c1", "c4", "c6", "c8"),
    !gs_collection_name %>% str_to_lower() %>% str_detect("_legacy"),
    !gs_subcollection %>% str_to_lower() %>% str_detect("vax")
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

gene_list_vector <- loadings$regional_PD_loading
names(gene_list_vector) <- loadings$ensembl_gene_id
gene_list_vector <- sort(gene_list_vector, decreasing = TRUE)

gsea_result <- GSEA(
  gene_list_vector,
  minGSSize = 10,
  maxGSSize = 500,
  TERM2GENE = MSIG_TERM2GENE,
  TERM2NAME = MSIG_TERM2NAME,
  pvalueCutoff = 1,
  eps = 1e-50,
  pAdjustMethod = "fdr",
  verbose = FALSE,
  seed = TRUE
)

gsea_df <- gsea_result %>%
  as_tibble() %>%
  rename(gs_id = ID, gs_name = Description, gs_size = setSize) %>%
  left_join(MSIG_METADATA, by = "gs_id") %>%
  relocate(starts_with("gs"), .before = enrichmentScore)

dmnx_gsea_ids <- c(
  "M17305",
  "M11131",
  "M37223",
  "M27446",
  "M37960",
  "M37959",
  "M26004",
  "M13015",
  "M47021",
  "M11274",
  "M10361",
  "M39508",
  "M42695",
  "M898",
  "M13812",
  "M39826",
  "M17230",
  "M49327",
  "M11847",
  "M22563"
)

dmnx_gsea_plot <- gsea_df %>%
  mutate(
    count = str_count(core_enrichment, "/") + 1,
    gene_ratio = count / gs_size,
    gs_name_pretty = prettify_genesets_names(gs_name, width = 100),
    gs_name_pretty = str_glue("{gs_name_pretty} [{gs_subcollection}]"),
    group = case_when(
      str_detect(gs_name, "MITO") ~ "mito",
      str_detect(gs_name, "INFLA|TYROBP|ANTI|CYTO|IMMU") ~ "immu",
      str_detect(gs_name, "NEURO") ~ "neu",
      str_detect(gs_name, "CILI") ~ "cil",
      str_detect(gs_name, "GLYC") ~ "glyc",
      TRUE ~ "other"
    )
  ) %>%
  filter(
    NES > 0,
    p.adjust < 0.01,
    !gs_subcollection %in%
      c(
        "IMMUNESIGDB",
        "TFT:GTRD",
        "MIR:MIRDB",
        "CGP",
        "CP:BIOCARTA",
        "CP:KEGG_MEDICUS",
        "CP:PID",
        "CP"
      ),
    group != "other",
    gs_id %in% dmnx_gsea_ids
  ) %>%
  ggplot(aes(x = gene_ratio, y = reorder(gs_name_pretty, gene_ratio))) +
  geom_point(aes(size = gs_size, color = -log10(p.adjust))) +
  facet_grid(
    forcats::fct_recode(
      forcats::fct_relevel(group, "mito", "immu", "neu", "glyc", "cil"),
      "Mitochondrial\nfunction" = "mito",
      "Immune\nresponse" = "immu",
      "Neuronal\nfunction" = "neu",
      "GAG\nsynthesis" = "glyc",
      "Cilia" = "cil"
    ) ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_size_continuous(range = c(.5, 3), name = "Gene\nset size") +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0))) +
  scale_color_viridis_c(
    option = "viridis",
    limits = c(2, 7),
    breaks = c(2, 3, 4, 5, 6, 7),
    direction = 1,
    name = expression(-log[10]("FDR"))
  ) +
  theme_publication(5) +
  coord_cartesian(clip = "off") +
  vertical_scientific_colormap +
  labs(x = "Gene ratio") +
  theme(
    panel.grid.major.y = element_line(color = "grey80", linewidth = .1),
    axis.title.y = element_blank(),
    strip.clip = "off"
  )

ggsave_pdf_svg(
  base_name = "05-results/02-plots/figure-3-pca-dmnx-signature-gsea",
  width = 120,
  height = 60,
  units = "mm",
  plot = dmnx_gsea_plot
)


# %% Figure 4 panel a: logFC Spearman correlation heatmap (ComplexHeatmap)
cor_matrix <- cor(
  dea_wo_global %>%
    mutate(context = glue("{cell_subtype_short} - {brain_region}")) %>%
    select(ensembl_gene_id, context, logfc) %>%
    pivot_wider(names_from = context, values_from = logfc) %>%
    column_to_rownames("ensembl_gene_id"),
  use = "pairwise.complete.obs",
  method = "spearman"
)

annotation_df <- data.frame(group = colnames(cor_matrix)) %>%
  tidyr::separate(
    group,
    into = c("cell_subtype_short", "brain_region"),
    sep = " - ",
    remove = FALSE
  ) %>%
  left_join(annotations, by = "cell_subtype_short") %>%
  tibble::column_to_rownames("group")

base_gp <- gpar(fontsize = 5, fontfamily = "Helvetica")

col_annotations <- HeatmapAnnotation(
  df = annotation_df %>%
    select(
      `Cell supertype` = cell_supertype_short,
      `Cell type` = cell_type_short,
      `Cell subtype` = cell_subtype_short,
      `Brain region` = brain_region
    ),
  col = list(
    `Cell supertype` = correct_cell_supertype_colors,
    `Cell type` = correct_cell_type_colors,
    `Cell subtype` = correct_cell_subtype_colors,
    `Brain region` = region_colors
  ),
  annotation_name_gp = base_gp,
  simple_anno_size = unit(1.5, "mm"),
  gp = gpar(lwd = 0.1, col = "black")
)

no_diagonal <- cor_matrix[upper.tri(cor_matrix)]
correlation_color <- colorRamp2(
  breaks = seq(min(no_diagonal), max(no_diagonal), length = 10),
  colors = viridisLite::inferno(10)
)

heatmap <- Heatmap(
  cor_matrix,
  name = "logFC Spearman correlation",
  col = correlation_color,
  clustering_distance_rows = function(x) as.dist(1 - x),
  clustering_distance_columns = function(x) as.dist(1 - x),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  column_split = 10,
  row_split = 10,
  top_annotation = col_annotations,
  show_column_names = FALSE,
  row_title = NULL,
  row_names_gp = base_gp,
  column_dend_height = unit(3, "mm"),
  column_dend_gp = gpar(lwd = .5),
  row_dend_width = unit(3, "mm"),
  row_dend_gp = gpar(lwd = .5),
  row_gap = unit(.5, "mm"),
  column_gap = unit(.5, "mm"),
  row_title_gp = base_gp,
  column_title_gp = base_gp,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 5, fontfamily = "Helvetica", fontface = "plain"),
    labels_gp = base_gp,
    direction = "horizontal",
    title_position = "topcenter",
    border = "black",
    grid_height = unit(2, "mm"),
    legend_gp = gpar(lwd = .5),
    tick_length = unit(-.25, "mm")
  )
)

pdf(
  file = "05-results/02-plots/figure-3-heatmap-logfc-correlation.pdf",
  width = 150 / 25.4,
  height = 150 / 25.4
)
draw(heatmap, show_annotation_legend = FALSE, heatmap_legend_side = "bottom")
dev.off()


# %% Stats related to figure 4
# PERMANOVA: quantify variance in logFC profiles explained by cell subtype
# identity vs brain region. Cell subtype expected to explain more, except for
# DMNX where regional context is also prominent.
dist_mat <- as.dist(1 - cor_matrix)
metadata_permanova <- annotation_df[labels(dist_mat), ]
adonis_res <- adonis2(
  dist_mat ~ cell_subtype_short + brain_region,
  data = metadata_permanova,
  by = "margin"
)
print(adonis_res)
