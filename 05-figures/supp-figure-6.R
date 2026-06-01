# Description: Generate Supplementary Figure 6 panels:
#   panel a — enrichment profiles of manually selected pathways across cell subtypes and brain regions
#   panel b — top dysregulated pathway per cell subtype and brain region

# %% Setup
library(tidyverse)
library(glue)
library(ggtext)
source("utils/my_utils.R")


# %% Load data
gsea <- readRDS("04-analysis/03-gsea-disease/gsea-disease-combined-results.rds") %>%
  filter(brain_region != "GLOBAL")


# %% Supplementary figure 6 panel a: manually curated pathways heatmap
gss <- c(
  "KEGG_RIBOSOME",
  "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION",
  "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",
  "KEGG_OXIDATIVE_PHOSPHORYLATION",
  "GOBP_CHAPERONE_MEDIATED_PROTEIN_FOLDING",
  "GOMF_ATP_DEPENDENT_PROTEIN_FOLDING_CHAPERONE",
  "REACTOME_HSF1_DEPENDENT_TRANSACTIVATION"
)

relevant_plot <- gsea %>%
  filter(gs_name %in% gss) %>%
  mutate(
    gs_name = factor(gs_name, levels = rev(gss)),
    gs_name_pretty = prettify_genesets_names(gs_name, width = 50, truncate = TRUE),
    gs_subcollection = if_else(gs_subcollection == "", "HALLMARK", gs_subcollection),
    gs_name_pretty = str_glue("{gs_name_pretty} [{gs_subcollection}]"),
    gs_name_pretty = factor(
      gs_name_pretty,
      levels = rev(c(
        "CHAPERONE MEDIATED PROTEIN FOLDING [GO:BP]",
        "ATP DEPENDENT PROTEIN FOLDING CHAPERONE [GO:MF]",
        "HSF1 DEPENDENT TRANSACTIVATION [CP:REACTOME]",
        "RESPIRATORY ELECTRON TRANSPORT [CP:REACTOME]",
        "OXIDATIVE PHOSPHORYLATION [CP:KEGG_LEGACY]",
        "EUKARYOTIC TRANSLATION ELONGATION [CP:REACTOME]",
        "RIBOSOME [CP:KEGG_LEGACY]"
      ))
    ),
    label_sig = if_else(p.adjust < 0.01 / 100, "*", "")
  ) %>%
  ggplot(aes(
    x = factor(brain_region, levels = c("DMNX", "GPI", "PFC", "PMC", "PVC")),
    y = gs_name_pretty,
    label = label_sig,
    fill = NES
  )) +
  geom_tile() +
  geom_text(vjust = .75, size = 2.5, color = "black") +
  facet_wrap(
    ~ factor(cell_subtype_short, correct_cell_subtype_order),
    ncol = 12,
    scales = "free_x"
  ) +
  scale_fill_distiller(palette = "RdBu", name = "NES\n(PD:Control)") +
  coord_cartesian(clip = "off") +
  tidytext::scale_y_reordered(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  theme_publication(5) +
  horizontal_scientific_colormap +
  theme(
    strip.text.x = element_text(hjust = 0.5, angle = 0),
    plot.margin = margin(r = 2, unit = "mm"),
    axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1),
    strip.clip = "off",
    panel.spacing.x = unit(2.5, "mm"),
    panel.background = element_rect(color = "black", fill = NA, linewidth = 0.2),
    strip.placement = "outside",
    axis.title = element_blank(),
    legend.key.width = unit(3, "mm"),
    legend.position = c(0.93, 0.2)
  )

ggsave_pdf_svg(
  base_name = "05-results/02-plots/figure-3-relevant-terms-per-ct-heatmap",
  width = 180,
  height = 45,
  units = "mm",
  plot = relevant_plot
)


# %% Supplementary figure 6 panel b: top pathway per cell subtype and region
x_labels_formatted <- setNames(
  glue("<span style='color:{region_colors};'>{names(region_colors)}</span>"),
  names(region_colors)
)

top_terms_summary <- gsea %>%
  filter(
    !gs_subcollection %in% c(
      "IMMUNESIGDB", "CGP", "HPO", "CP", "CP:BIOCARTA",
      "CP:PID", "CP:KEGG_MEDICUS", "TFT:GTRD", "MIR:MIRDB"
    ),
    p.adjust < 0.01 / 100
  ) %>%
  group_by(cell_subtype_short, brain_region) %>%
  slice_max(order_by = abs(NES), n = 1) %>%
  ungroup() %>%
  mutate(brain_region = factor(brain_region, levels = correct_brain_region_order)) %>%
  arrange(cell_subtype_short, brain_region) %>%
  group_by(cell_subtype_short, gs_name) %>%
  summarise(
    regions_where_top = paste(
      glue("<span style='color:{region_colors[as.character(brain_region)]};'>{brain_region}</span>"),
      collapse = ", "
    ),
    primary_region = first(brain_region),
    .groups = "drop"
  )

plot_df <- gsea %>%
  inner_join(top_terms_summary, by = c("cell_subtype_short", "gs_name")) %>%
  mutate(
    gs_name_pretty = prettify_genesets_names(gs_name, width = 50, truncate = TRUE),
    gs_subcollection = case_when(
      gs_subcollection == "" ~ "HALLMARK",
      gs_subcollection == "CP:REACTOME" ~ "REACT.",
      gs_subcollection == "CP:KEGG_LEGACY" ~ "KEGG",
      gs_subcollection == "CP:KEGG_MEDICUS" ~ "KEGG_MED.",
      gs_subcollection == "CP:WIKIPATHWAYS" ~ "WP",
      TRUE ~ gs_subcollection
    ),
    gs_name_pretty = str_glue("{gs_name_pretty} [{gs_subcollection}]"),
    gs_name_colored = glue(
      "{gs_name_pretty} <span style='font-size:5pt;'>({regions_where_top})</span>"
    ),
    region_rank = match(as.character(primary_region), rev(correct_brain_region_order)),
    gs_name_ordered = tidytext::reorder_within(gs_name_colored, region_rank, cell_subtype_short)
  )

gsea_plot <- plot_df %>%
  ggplot(aes(
    x = factor(brain_region, levels = correct_brain_region_order),
    y = gs_name_ordered,
    fill = NES
  )) +
  geom_tile() +
  geom_text(
    aes(label = if_else(p.adjust < 0.01 / 100, "*", "")),
    vjust = 0.75,
    size = 2.5
  ) +
  facet_wrap(
    vars(factor(cell_subtype_short, levels = correct_cell_subtype_order)),
    ncol = 2,
    scales = "free"
  ) +
  tidytext::scale_y_reordered(expand = c(0, 0)) +
  scale_x_discrete(labels = x_labels_formatted, expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  scale_fill_distiller(palette = "RdBu", name = "NES\n(PD:Control)") +
  theme_publication(5) +
  vertical_scientific_colormap +
  theme(
    axis.text.y = element_markdown(lineheight = 0.9),
    axis.text.x = element_markdown(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    strip.clip = "off",
    legend.position = "right"
  )

ggsave_pdf_svg(
  base_name = "05-results/02-plots/figure-3-top-terms-per-ct-and-region-heatmap",
  width = 170,
  height = 200,
  units = "mm",
  plot = gsea_plot
)