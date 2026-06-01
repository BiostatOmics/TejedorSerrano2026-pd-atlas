# Description: Generate Supplementary Figure 5 panels:
#   panel a — up/down power-corrected DEG barplot
#   panel b — top shared DEGs across cell subtypes heatmap
#   panel c — top divergent DEGs across cell subtypes heatmap
#   panel d — top divergent DEGs across regions within cell subtype heatmap

# %% Setup
library(tidyverse)
library(patchwork)
library(glue)
library(ComplexHeatmap)
library(circlize)
source("utils/my_utils.R")


# %% Load data
COMBINED_DEA_RESULTS <- "04-analysis/pseudobulk-dea-combined-results-ebayes-fc-1_00.rds"
VOOMFIT_LIST <- "04-analysis/pseudobulk-dea-voomLmFit-per-cell-type.rds"

dea <- readRDS(COMBINED_DEA_RESULTS)
dea_wo_global <- dea %>%
  filter(brain_region != "GLOBAL") %>%
  mutate(symbol = if_else(is.na(hgnc_symbol), ensembl_gene_id, hgnc_symbol))

voomfits <- readRDS(VOOMFIT_LIST)

pseudobulk_metadata <- map_dfr(voomfits, ~ .x$targets %>% as_tibble())


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

mm_model <- nls(
  n_degs ~ (DEGs_max * median_psbulk_cells^n) / (Km^n + median_psbulk_cells^n),
  data = n_degs_vs_n_cells,
  start = list(
    DEGs_max = max(n_degs_vs_n_cells$n_degs),
    Km = median(n_degs_vs_n_cells$median_psbulk_cells),
    n = 1
  )
)

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


# %% Supplementary figure 3 panel a: up/down power-corrected DEG barplot
corrected_updown_degs_barplot <- n_degs_vs_n_cells %>%
  select(
    cell_subtype_short,
    brain_region,
    corrected_n_degs_up,
    corrected_n_degs_down
  ) %>%
  pivot_longer(
    cols = c(corrected_n_degs_up, corrected_n_degs_down),
    names_to = "direction",
    values_to = "corrected_n_degs"
  ) %>%
  ggplot(aes(
    y = corrected_n_degs,
    x = factor(brain_region, levels = c("DMNX", "GPI", "PFC", "PMC", "PVC")),
    fill = direction
  )) +
  geom_col() +
  geom_hline(yintercept = 0, linewidth = .1) +
  facet_wrap(
    ~ factor(cell_subtype_short, levels = correct_cell_subtype_order),
    nrow = 2,
    scales = "free_x"
  ) +
  scale_y_continuous(
    breaks = scales::breaks_width(1000),
    limits = c(-3000, 3000),
    labels = scales::label_number(scale = 1e-3),
  ) +
  scale_fill_manual(
    breaks = c("corrected_n_degs_up", "corrected_n_degs_down"),
    labels = c(
      "corrected_n_degs_up" = "Up (FC > 1)",
      "corrected_n_degs_down" = "Down (FC < 1)"
    ),
    name = "Direction of regulation:",
    values = c(
      "corrected_n_degs_up" = cool_red,
      "corrected_n_degs_down" = cool_blue
    )
  ) +
  labs(y = "Number of power-corrected DEGs (\u00d710\u00b3)") +
  theme_publication(5) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
    axis.title.y = ggtext::element_markdown(),
    panel.grid.major.y = element_line(color = "grey80", linewidth = 0.05),
    panel.grid.minor.y = element_line(color = "grey95", linewidth = 0.05),
    legend.position = "right",
    legend.direction = "vertical",
    strip.clip = "off"
  )

ggsave_pdf_svg(
  base_name = "05-results/02-plots/figure-3-updown-corrected-ndegs-barplot",
  width = 150,
  height = 50,
  units = "mm",
  plot = corrected_updown_degs_barplot
)


# %% Compute top shared DEGs across cell subtypes (manual curation of gene list)
pan_pd_DEGs <- c(
  "RELA",
  "ANK3",
  "RBPJ",
  "CIRBP",
  "TARBP1",
  "HMGN1",
  "DNAJC10",
  "EEF1D",
  "PIK3CA",
  "IRS2",
  "IGF1R",
  "NEDD4L",
  "UBE2D3"
)


# %% Define shared horizontal heatmap template
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


# %% Supplementary figure 3 panel b: top shared DEGs across cell subtypes heatmap
top_shared_degs_across_cts_df <- dea_wo_global %>%
  filter(symbol %in% pan_pd_DEGs, !is.na(hgnc_symbol)) %>%
  mutate(label_sig = case_when(fdr < 0.05 ~ "*", TRUE ~ ""))

top_shared_degs_across_cts_plot <- top_shared_degs_across_cts_df %>%
  ggplot() +
  horizontal_heatmap_template

ggsave_pdf_svg(
  base_name = "05-results/02-plots/figure-3-top-shared-degs-across-cts-heatmap",
  width = 150,
  height = 50,
  units = "mm",
  plot = top_shared_degs_across_cts_plot
)


# %% Supplementary figure 3 panel c: top divergent DEGs across cell subtypes heatmap
divergent_DEGs_across_cts <- dea_wo_global %>%
  filter(fdr < 0.05) %>%
  distinct(
    ensembl_gene_id,
    cell_subtype_short,
    brain_region,
    .keep_all = TRUE
  ) %>%
  group_by(ensembl_gene_id) %>%
  summarise(
    n_cts_where_deg_is_divergent = n_distinct(sign(logfc)),
    n_cts_where_is_deg = n(),
    sd = sd(logfc) * n_cts_where_is_deg,
    .groups = "drop"
  ) %>%
  slice_max(
    order_by = tibble(n_cts_where_deg_is_divergent, n_cts_where_is_deg),
    n = 10
  ) %>%
  pull(ensembl_gene_id)

top_divergent_degs_across_cts_df <- dea_wo_global %>%
  filter(
    ensembl_gene_id %in% divergent_DEGs_across_cts,
    !is.na(hgnc_symbol)
  ) %>%
  mutate(label_sig = case_when(fdr < 0.05 ~ "*", TRUE ~ ""))

top_divergent_degs_across_cts_plot <- top_divergent_degs_across_cts_df %>%
  group_by(cell_subtype_short, brain_region) %>%
  filter(any(fdr < 0.05)) %>%
  ungroup() %>%
  ggplot(aes(y = forcats::fct_reorder(hgnc_symbol, logfc, sd))) +
  horizontal_heatmap_template

ggsave_pdf_svg(
  base_name = "05-results/02-plots/figure-3-top-divergent-degs-across-cts-heatmap",
  width = 110,
  height = 50,
  units = "mm",
  plot = top_divergent_degs_across_cts_plot
)


# %% Supplementary figure 3 panel d: top divergent DEGs across regions within cell subtype
divergent_DEGs_bt_reg_wt_ct <- dea_wo_global %>%
  filter(fdr < 0.05) %>%
  group_by(cell_subtype_short, ensembl_gene_id) %>%
  filter(n() > 1) %>%
  summarise(divergent = n_distinct(sign(logfc)) > 1, .groups = "drop") %>%
  filter(divergent) %>%
  pull(ensembl_gene_id)

top_divergent_degs_bt_reg_wt_ct_df <- dea_wo_global %>%
  filter(ensembl_gene_id %in% divergent_DEGs_bt_reg_wt_ct) %>%
  mutate(label_sig = case_when(fdr < 0.05 ~ "*", TRUE ~ ""))

top_divergent_degs_bt_reg_wt_ct_plot <- top_divergent_degs_bt_reg_wt_ct_df %>%
  group_by(cell_subtype_short, brain_region) %>%
  filter(any(fdr < 0.05)) %>%
  ungroup() %>%
  group_by(cell_subtype_short, ensembl_gene_id) %>%
  filter(fdr < 0.05, n_distinct(sign(logfc)) > 1) %>%
  mutate(hgnc_symbol = symbol) %>%
  ggplot() +
  horizontal_heatmap_template +
  ggforce::facet_col(
    vars(factor(cell_subtype_short, correct_cell_subtype_order)),
    scales = "free",
    space = "free"
  ) +
  theme(strip.text.x = element_text(hjust = 0.5, angle = 0)) +
  no_legend

ggsave_pdf_svg(
  base_name = "05-results/02-plots/figure-3-top-divergent-degs-across-regions-within-celltype-heatmap",
  width = 26,
  height = 35,
  units = "mm",
  plot = top_divergent_degs_bt_reg_wt_ct_plot
)
