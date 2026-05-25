# Description: Generate Figure 3 DEA panels:
#   panel a — raw DEG count barplot
#   panel b — DEG count vs cell count scatterplot (Michaelis-Menten fit)
#   panel c — power-corrected DEG heatmap
#   panel d — shared DEGs across cell subtypes barplot
#   panel e — divergence summary barplot
#   panel f — shared DEGs across brain regions barplots
#   panel g — top DEGs per cell subtype heatmap

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
COMBINED_DEA_RESULTS <- "04-analysis/pseudobulk-dea-combined-results-ebayes-fc-1_00.rds"
VOOMFIT_LIST <- "04-analysis/pseudobulk-dea-voomLmFit-per-cell-type.rds"
ANNOTATIONS <- "03-data/annotations_7_11_agg.txt"

dea <- readRDS(COMBINED_DEA_RESULTS)
dea_wo_global <- dea %>%
  filter(brain_region != "GLOBAL") %>%
  mutate(symbol = if_else(is.na(hgnc_symbol), ensembl_gene_id, hgnc_symbol))

voomfits <- readRDS(VOOMFIT_LIST)

pseudobulk_metadata <- map_dfr(voomfits, ~ .x$targets %>% as_tibble())

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


# %% Compute DEG counts vs pseudobulk cell counts per contrast
n_degs_vs_n_cells <- dea_wo_global %>%
  group_by(cell_subtype_short, brain_region) %>%
  summarise(
    n_degs = sum(fdr < 0.05),
    n_degs_up = sum(fdr < 0.05 & logfc > 0),
    n_degs_down = sum(fdr < 0.05 & logfc < 0),
    .groups = "drop"
  ) %>%
  left_join(
    pseudobulk_metadata %>%
      group_by(cell_subtype_short, brain_region) %>%
      summarise(median_psbulk_cells = median(psbulk_cells), .groups = "drop")
  )

# Fit Michaelis-Menten model: DEGs ~ f(median pseudobulk cell count)
mm_model <- nls(
  n_degs ~ (DEGs_max * median_psbulk_cells^n) / (Km^n + median_psbulk_cells^n),
  data = n_degs_vs_n_cells,
  start = list(
    DEGs_max = max(n_degs_vs_n_cells$n_degs),
    Km = median(n_degs_vs_n_cells$median_psbulk_cells),
    n = 1
  )
)
summary(mm_model)

n_degs_vs_n_cells <- n_degs_vs_n_cells %>%
  mutate(
    predicted_n_degs = predict(mm_model),
    corrected_n_degs = n_degs - predicted_n_degs,
    corrected_n_degs_up = if_else(
      n_degs == 0,
      corrected_n_degs * 0.5,
      corrected_n_degs * n_degs_up / n_degs
    ),
    corrected_n_degs_down = if_else(
      n_degs == 0,
      corrected_n_degs * 0.5,
      corrected_n_degs * n_degs_down / n_degs
    )
  )

# Summary statistics for main text
dea_wo_global %>% filter(fdr < 0.05) %>% distinct(ensembl_gene_id) %>% nrow()
cor.test(
  n_degs_vs_n_cells$n_degs,
  n_degs_vs_n_cells$median_psbulk_cells,
  method = "spearman"
)


# %% Figure 3 panel b: DEG count vs pseudobulk cell count scatterplot
ndegs_vs_ncells_plot <- n_degs_vs_n_cells %>%
  ggplot(aes(x = median_psbulk_cells)) +
  geom_segment(
    data = n_degs_vs_n_cells %>% filter(corrected_n_degs > 0),
    aes(y = predicted_n_degs, yend = n_degs),
    color = cool_red,
    linewidth = 0.1,
    alpha = .7
  ) +
  geom_segment(
    data = n_degs_vs_n_cells %>% filter(corrected_n_degs < 0),
    aes(y = predicted_n_degs, yend = n_degs),
    color = cool_blue,
    linewidth = 0.1,
    alpha = .7
  ) +
  geom_point(aes(
    y = n_degs,
    color = cell_subtype_short,
    shape = brain_region
  )) +
  geom_line(aes(y = predicted_n_degs)) +
  scale_x_log10(guide = "axis_logticks") +
  scale_colour_manual(
    values = correct_cell_subtype_colors,
    name = "Cell subtypes:",
    guide = guide_legend(ncol = 2)
  ) +
  scale_shape_manual(
    name = "Brain region of origin:",
    values = region_shapes,
    guide = guide_legend(ncol = 2)
  ) +
  labs(
    x = "Median pseudobulk cell count per contrast\n(log scale)",
    y = "Number of DEGs\n(FDR < 0.05)"
  ) +
  theme_publication(base_size = 5) +
  theme(legend.position = "bottom", legend.direction = "vertical")

ggsave_pdf_svg(
  base_name = "05-results/02-plots/figure-3-ndegs-vs-ncells-scatterplot",
  width = 60,
  height = 60,
  units = "mm",
  plot = ndegs_vs_ncells_plot
)


# %% Figure 3 panel c: power-corrected DEG heatmap with marginal barplots
signed_sqrt_trans <- scales::new_transform(
  name = "signed_sqrt",
  transform = function(x) sign(x) * sqrt(abs(x)),
  inverse = function(x) sign(x) * x^2
)

common_scale_fill <- fill_distiller_gradient(
  name = "Power-corrected\nnumber of DEGs",
  transform = signed_sqrt_trans,
  limits = range(n_degs_vs_n_cells$corrected_n_degs),
  breaks = c(-2000, -1000, -500, 0, 500, 1000, 2000),
)

heatmap <- n_degs_vs_n_cells %>%
  ggplot() +
  geom_tile(aes(
    x = brain_region,
    y = factor(cell_subtype_short, rev(correct_cell_subtype_order)),
    fill = corrected_n_degs
  )) +
  common_scale_fill +
  theme_publication(5) +
  vertical_scientific_colormap +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
  )

barplots_template <- list(
  geom_col(width = 1, col = "white", linewidth = .1),
  common_scale_fill,
  theme_publication(5),
  no_legend,
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())
)

right_barplot <- n_degs_vs_n_cells %>%
  group_by(cell_subtype_short) %>%
  summarise(
    mean_corrected_n_degs = mean(corrected_n_degs),
    .groups = "drop"
  ) %>%
  ggplot(aes(
    x = mean_corrected_n_degs,
    y = factor(cell_subtype_short, rev(correct_cell_subtype_order)),
    fill = mean_corrected_n_degs
  )) +
  barplots_template +
  geom_vline(xintercept = 0, linewidth = .1) +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
  )

top_barplot <- n_degs_vs_n_cells %>%
  group_by(brain_region) %>%
  summarise(
    mean_corrected_n_degs = mean(corrected_n_degs),
    .groups = "drop"
  ) %>%
  ggplot(aes(
    x = brain_region,
    y = mean_corrected_n_degs,
    fill = mean_corrected_n_degs
  )) +
  barplots_template +
  geom_hline(yintercept = 0, linewidth = .1) +
  theme(
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
  )

legend <- cowplot::get_legend(heatmap)

design <- "
AB
CD
"
combined <- top_barplot +
  legend +
  heatmap +
  no_legend +
  right_barplot +
  plot_layout(design = design, widths = c(2, 2), heights = c(1, 2.5))

ggsave_pdf_svg(
  base_name = "05-results/02-plots/figure-3-ndegs-heatmap",
  width = 40,
  height = 70,
  units = "mm",
  plot = combined
)


# %% Figure 3 panel a: raw DEG count barplot
raw_n_degs_plot <- n_degs_vs_n_cells %>%
  ggplot(aes(
    y = fct_reorder(cell_subtype_short, n_degs, .fun = sum, .desc = TRUE),
    fill = factor(brain_region, levels = rev(correct_brain_region_order)),
    x = n_degs
  )) +
  geom_col(color = "white", linewidth = .1, width = 1) +
  geom_text(
    data = n_degs_vs_n_cells %>%
      group_by(cell_subtype_short) %>%
      summarise(n_degs_total = sum(n_degs), .groups = "drop"),
    aes(
      label = n_degs_total,
      y = fct_reorder(cell_subtype_short, n_degs_total, .desc = TRUE),
      x = n_degs_total
    ),
    inherit.aes = FALSE,
    hjust = -0.1
  ) +
  scale_fill_manual(values = region_colors, name = "Brain region of origin:") +
  labs(x = "Number of DEGs\n(FDR < 0.05)") +
  theme_publication(5) +
  no_legend +
  coord_cartesian(clip = "off") +
  theme(
    axis.title.y = element_blank(),
    plot.margin = margin(r = 4, unit = "mm")
  )

ggsave_pdf_svg(
  base_name = "05-results/02-plots/figure-3-raw-ndegs-barplot",
  width = 45,
  height = 40,
  units = "mm",
  plot = raw_n_degs_plot
)


# %% Figure 3 panel f: shared DEGs across brain regions (summary, compact)
shared_degs_across_regions_df <- dea_wo_global %>%
  filter(fdr < 0.05) %>%
  group_by(cell_subtype_short, ensembl_gene_id) %>%
  summarise(n_regions_where_is_deg = n(), .groups = "drop") %>%
  mutate(n_regions_where_is_deg = factor(n_regions_where_is_deg)) %>%
  group_by(n_regions_where_is_deg, cell_subtype_short) %>%
  summarise(n_degs = n(), .groups = "drop") %>%
  mutate(
    sharing_category = if_else(
      n_regions_where_is_deg == 1,
      "Specific",
      "Shared"
    ),
    sharing_category = factor(
      sharing_category,
      levels = c("Specific", "Shared")
    )
  ) %>%
  group_by(cell_subtype_short) %>%
  mutate(
    subtype_total_degs = sum(n_degs * as.numeric(n_regions_where_is_deg)),
    subtype_unique_degs = sum(n_degs),
    prop_degs = n_degs / subtype_unique_degs,
  )

summary_shared_degs_across_regions_plot <- ggplot(
  shared_degs_across_regions_df,
  aes(x = sharing_category, y = n_degs, fill = n_regions_where_is_deg)
) +
  facet_wrap(~cell_subtype_short, scales = "free_y", nrow = 2) +
  geom_col(col = "white", linewidth = .1, width = 1) +
  geom_text(
    data = . %>%
      group_by(cell_subtype_short, sharing_category) %>%
      summarise(y_pos = sum(n_degs), .groups = "drop_last") %>%
      mutate(prop = y_pos / sum(y_pos)) %>%
      ungroup(),
    aes(
      x = sharing_category,
      y = y_pos,
      label = scales::percent(prop, accuracy = 1)
    ),
    inherit.aes = FALSE,
    vjust = -0.1
  ) +
  scale_fill_manual(
    name = "Number of brain regions\nin which a gene is DE",
    values = c(
      "1" = "#420a68",
      "2" = "#dd513a",
      "3" = "#f3771a",
      "4" = "#fca50a",
      "5" = "#f6d645"
    )
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(3),
    expand = expansion(mult = c(0, 0.2))
  ) +
  coord_cartesian(clip = "off") +
  labs(y = "Number of DEGs\n(FDR < 0.05)") +
  theme_publication(5) +
  theme(
    strip.clip = "off",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
    axis.title.x = element_blank(),
    legend.position = "right",
    legend.direction = "vertical",
    legend.title.position = "top"
  )

ggsave_pdf_svg(
  base_name = "05-results/02-plots/figure-3-summary-shared-degs-across-regions-barplots",
  width = 180,
  height = 30,
  units = "mm",
  plot = summary_shared_degs_across_regions_plot
)

# Statistics for main text: proportion of region-specific vs shared DEGs
shared_degs_across_regions_df %>%
  group_by(sharing_category) %>%
  summarise(sum_n_degs = sum(n_degs)) %>%
  mutate(prop = sum_n_degs / sum(sum_n_degs))


# %% Figure 3 panel d: shared DEGs across cell subtypes barplot
shared_degs_across_cts_df <- dea_wo_global %>%
  filter(fdr < 0.05) %>%
  distinct(ensembl_gene_id, cell_subtype_short) %>%
  group_by(ensembl_gene_id) %>%
  summarise(n_cts_where_is_deg = n(), .groups = "drop") %>%
  mutate(
    n_cts_where_is_deg = case_when(
      n_cts_where_is_deg %in% 1:4 ~ as.character(n_cts_where_is_deg),
      n_cts_where_is_deg >= 5 & n_cts_where_is_deg <= 10 ~ "5-10",
      n_cts_where_is_deg > 10 ~ "10+"
    ),
    n_cts_where_is_deg = factor(
      n_cts_where_is_deg,
      levels = c("1", "2", "3", "4", "5-10", "10+")
    )
  ) %>%
  group_by(n_cts_where_is_deg) %>%
  summarise(n_degs = n(), .groups = "drop") %>%
  mutate(
    sharing_category = if_else(n_cts_where_is_deg == 1, "Specific", "Shared"),
    sharing_category = factor(
      sharing_category,
      levels = c("Specific", "Shared")
    ),
    total_degs = sum(n_degs * as.numeric(n_cts_where_is_deg)),
    unique_degs = sum(n_degs),
    prop_degs = n_degs / unique_degs,
  )

# Summary statistics for main text
shared_degs_across_cts_df

shared_degs_across_cts_plot <- shared_degs_across_cts_df %>%
  ggplot(aes(x = n_cts_where_is_deg, y = n_degs, fill = n_cts_where_is_deg)) +
  geom_col(col = "white", linewidth = .1, width = 1) +
  geom_text(
    mapping = aes(label = scales::percent(prop_degs, accuracy = 1)),
    vjust = -0.1
  ) +
  scale_fill_manual(
    guide = "none",
    values = c(
      "1" = "#420a68",
      "2" = "#dd513a",
      "3" = "#f3771a",
      "4" = "#fca50a",
      "5-10" = "#f6d645",
      "10+" = "#ffffb2"
    )
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(5),
    expand = expansion(mult = c(0, 0.2))
  ) +
  labs(
    x = "Number of cell subtypes\nin which a gene is DE",
    y = "Number of DEGs\n(FDR < 0.05)"
  ) +
  theme_publication(5)

ggsave_pdf_svg(
  base_name = "05-results/02-plots/figure-3-shared-degs-across-cts-barplots",
  width = 45,
  height = 40,
  units = "mm",
  plot = shared_degs_across_cts_plot
)


# %% Figure 3 panel e: divergence summary barplot
divergence_across_regions <- dea_wo_global %>%
  filter(fdr < 0.05) %>%
  group_by(cell_subtype_short, ensembl_gene_id) %>%
  filter(n() > 1) %>%
  group_by(cell_subtype_short, ensembl_gene_id) %>%
  summarise(divergent = n_distinct(sign(logfc)) > 1, .groups = "drop") %>%
  count(divergent) %>%
  mutate(comparison = "region")

divergence_across_cts <- dea_wo_global %>%
  filter(fdr < 0.05) %>%
  group_by(ensembl_gene_id, cell_subtype_short) %>%
  summarise(
    mean_logfc = mean(logfc),
    sign = sign(mean_logfc),
    divergent_regions = n_distinct(sign(logfc)) > 1,
    n_regions = n(),
    .groups = "drop"
  ) %>%
  group_by(ensembl_gene_id) %>%
  filter(n() > 1) %>%
  summarise(
    divergent = n_distinct(sign(mean_logfc)) > 1,
    n_cts_where_is_deg = n(),
    .groups = "drop"
  ) %>%
  count(divergent) %>%
  mutate(comparison = "cell_type")

divergence_summary_df <- bind_rows(
  divergence_across_regions,
  divergence_across_cts
)

# Summary statistics for main text
divergence_summary_df %>%
  group_by(comparison) %>%
  mutate(pct = n / sum(n) * 100)

divergence_summary_plot <- divergence_summary_df %>%
  ggplot(aes(
    x = n,
    y = factor(comparison, levels = rev(c("cell_type", "region"))),
    fill = divergent
  )) +
  geom_col(position = "fill") +
  scale_fill_manual(
    name = "Regulatory\ndirection:",
    labels = c("TRUE" = "Divergent", "FALSE" = "Consistent"),
    values = c("TRUE" = "grey20", "FALSE" = "grey80")
  ) +
  scale_x_continuous(labels = scales::label_percent()) +
  scale_y_discrete(
    labels = c(
      "cell_type" = "Across cell subtypes\n(all regions)",
      "region" = "Across regions\n(within cell subtypes)"
    )
  ) +
  ggforce::facet_zoom(xlim = c(0, 0.01), zoom.size = .75) +
  labs(x = "Proportion of shared DEGs") +
  theme_publication(5) +
  theme(
    axis.title.y = element_blank(),
    strip.background = element_rect(
      fill = "grey80",
      color = "grey95",
      linewidth = .05
    )
  )

ggsave_pdf_svg(
  base_name = "05-results/02-plots/figure-3-divergence-summary-barplot",
  width = 55,
  height = 30,
  units = "mm",
  plot = divergence_summary_plot
)


# %% Figure 3 panel g: top 10 DEGs per cell subtype horizontal heatmap

# Define shared horizontal heatmap template
horizontal_heatmap_template <- list(
  aes(
    x = brain_region,
    y = forcats::fct_reorder(hgnc_symbol, logfc, median),
    label = label_sig,
    fill = logfc
  ),
  facet_grid(
    cols = vars(factor(cell_subtype_short, correct_cell_subtype_order)),
    scales = "free_x",
    space = "free_x",
    labeller = label_wrap_gen(width = 10)
  ),
  geom_tile(),
  geom_text(vjust = .75, size = 2.5),
  scale_fill_distiller(
    palette = "RdBu",
    limits = c(-1, 1),
    oob = scales::oob_squish,
    name = "log<sub>2</sub> Fold Change<br>(PD:Control)"
  ),
  coord_cartesian(clip = "off"),
  scale_x_discrete(expand = c(0, 0)),
  scale_y_discrete(position = "left", expand = c(0, 0)),
  theme_publication(5),
  horizontal_scientific_colormap,
  theme(
    axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1),
    axis.text.y = element_text(face = "italic"),
    strip.placement = "outside",
    axis.title = element_blank(),
    strip.text.x = element_text(hjust = 0, angle = 90),
    strip.clip = "off",
    panel.spacing.x = unit(1, "mm"),
    panel.background = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.2
    ),
    legend.title = ggtext::element_markdown(),
    legend.position = "bottom"
  )
)

top_genes_per_ct <- dea_wo_global %>%
  filter(fdr < 0.05, !is.na(hgnc_symbol)) %>%
  group_by(cell_subtype_short, hgnc_symbol) %>%
  summarise(max_abs_logfc = max(abs(logfc)), .groups = "drop") %>%
  group_by(cell_subtype_short) %>%
  slice_max(order_by = max_abs_logfc, n = 10) %>%
  ungroup()

top_10_degs_per_ct_df <- dea_wo_global %>%
  inner_join(top_genes_per_ct, by = c("cell_subtype_short", "hgnc_symbol")) %>%
  mutate(
    hgnc_symbol_ordered = tidytext::reorder_within(
      hgnc_symbol,
      logfc,
      cell_subtype_short,
      fun = max
    ),
    label_sig = case_when(fdr < 0.05 ~ "*", TRUE ~ "")
  )

top_10_degs_per_ct_plot <- top_10_degs_per_ct_df %>%
  ggplot() +
  horizontal_heatmap_template +
  aes(y = hgnc_symbol_ordered) +
  facet_wrap(
    vars(factor(cell_subtype_short, levels = correct_cell_subtype_order)),
    nrow = 3,
    scales = "free"
  ) +
  tidytext::scale_y_reordered(expand = c(0, 0)) +
  scale_fill_distiller(
    palette = "RdBu",
    limits = c(-2, 2),
    oob = scales::oob_squish,
    name = "log<sub>2</sub> Fold Change<br>(PD:Control)"
  ) +
  coord_cartesian(clip = "off") +
  theme(
    strip.text.x = element_text(hjust = 0.5, angle = 0),
    plot.margin = margin(r = 2, unit = "mm"),
    legend.position = c(0.9, 0.1)
  )

ggsave_pdf_svg(
  base_name = "05-results/02-plots/figure-3-top-degs-per-ct-heatmap",
  width = 180,
  height = 80,
  units = "mm",
  plot = top_10_degs_per_ct_plot
)
