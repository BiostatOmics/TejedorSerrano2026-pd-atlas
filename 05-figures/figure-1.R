# Description: Generate all Figure 1 main panels (a through j):
#   panel a   — scheme icons
#   panels b-e — cohort description (barplots and density plots)
#   panel f   — cell subtype UMAP
#   panels g-h — marker gene dotplots
#   panels i-j — cell composition barplots

# %% Libraries
library(tidyverse)
library(glue)
library(svglite)
library(systemfonts)
library(sysfonts)
library(ggrastr)
library(ggrepel)
library(scales)
library(ggtext)
library(patchwork)
library(cowplot)
library(legendry)
library(edgeR)
library(conflicted)
source("utils/my_utils.R")
source("utils/figure_helpers.R")
conflicts_prefer(dplyr::filter, dplyr::select)
set.seed(42)


# Panel a: Scheme icons

set.seed(1)

# %% Compositional bar chart icon
tibble(
  x = rep(c("x", "y"), each = 3),
  color = rep(c("a", "b", "c"), 2),
  y = c(0.33, 0.33, 0.33, 0.56, 0.1, 0.33)
) %>%
  ggplot(aes(x = x, fill = color, y = y)) +
  geom_bar(
    stat = "identity",
    position = "stack",
    color = "black",
    linewidth = .1,
    width = .5
  ) +
  scale_y_continuous(labels = scales::label_percent()) +
  scale_x_discrete(labels = c("x" = "Healthy", "y" = "Disease")) +
  scale_fill_manual(
    values = c("a" = cool_red, "b" = cool_blue, "c" = cool_grey)
  ) +
  theme_publication(base_size = 4) +
  no_legend +
  theme(axis.title = element_blank())

ggsave_pdf_svg(
  base_name = "05-results/02-plots/figure-1-scheme-compositional",
  width = 20,
  height = 15,
  units = "mm"
)


# %% Volcano plot icon
n_genes <- 500

tibble(
  x = c(
    rnorm(n_genes * 0.8, mean = 0, sd = 0.5),
    rnorm(n_genes * 0.1, mean = 2, sd = 0.8),
    rnorm(n_genes * 0.1, mean = -2, sd = 0.8)
  )
) %>%
  mutate(
    y = abs(x) * runif(n_genes, 0.5, 3) + rexp(n_genes, rate = 5),
    status = case_when(
      x > 1 & y > -log10(0.05) ~ "up",
      x < -1 & y > -log10(0.05) ~ "down",
      TRUE ~ "ns"
    )
  ) %>%
  ggplot(aes(x = x, y = y, color = status)) +
  geom_point(alpha = .9, size = .4, shape = 16) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dotted",
    color = "black",
    linewidth = .1
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dotted",
    color = "black",
    linewidth = .1
  ) +
  scale_color_manual(
    values = c("up" = cool_red, "down" = cool_blue, "ns" = cool_grey)
  ) +
  labs(
    x = expression(Log[2] ~ "Fold Change"),
    y = expression(-Log[10] ~ "(FDR)")
  ) +
  theme_publication(base_size = 4) +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

ggsave_pdf_svg(
  base_name = "05-results/02-plots/figure-1-scheme-volcano",
  width = 15,
  height = 15,
  units = "mm"
)


# %% GSEA lollipop icon
tibble(
  pathway = factor(paste0("Pathway", 1:5)),
  nes = c(2.8, 2.2, 1.8, -1.5, -1.2),
  gene_count = c(85, 60, 45, 25, 15),
  significance = c(5.0, 3.5, 2.2, 1.5, 1.2)
) %>%
  mutate(pathway = fct_reorder(pathway, significance, .desc = FALSE)) %>%
  ggplot(aes(x = significance, y = pathway, size = gene_count, color = nes)) +
  geom_segment(
    aes(x = 0, xend = significance, y = pathway, yend = pathway),
    color = "grey50",
    linewidth = 0.15
  ) +
  geom_point() +
  scale_size_continuous(range = c(0.5, 2)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_color_gradient(low = cool_blue, high = cool_red) +
  labs(x = -Log[10] ~ "(FDR)", y = "Pathways") +
  theme_publication(base_size = 4) +
  theme(
    legend.position = "none",
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  )

ggsave_pdf_svg(
  base_name = "05-results/02-plots/figure-1-scheme-gsea",
  width = 16,
  height = 15,
  units = "mm"
)


# %% UMAP icon
n1 <- 60
n2 <- 45
n3 <- 35

df_umap <- tibble(
  umap1 = c(rnorm(n1, 0, 0.4), rnorm(n2, 2.5, 0.5), rnorm(n3, 3.5, 0.4)),
  umap2 = c(rnorm(n1, 0, 0.4), rnorm(n2, 2.5, 0.3), rnorm(n3, -0.5, 0.4)),
  cluster = rep(c("C1", "C2", "C3"), c(n1, n2, n3))
)

min_x <- min(df_umap$umap1) - 0.5
min_y <- min(df_umap$umap2) - 0.5
range_x <- diff(range(df_umap$umap1))
range_y <- diff(range(df_umap$umap2))
len_x <- range_x * 0.75
len_y <- range_y * 0.75

p_umap <- df_umap %>%
  ggplot(aes(x = umap1, y = umap2)) +
  geom_point(aes(color = cluster), size = 0.2, alpha = 0.8, shape = 16) +
  scale_color_manual(
    values = c("C1" = "#a21082", "C2" = "#79929e", "C3" = "#107161")
  ) +
  geom_segment(
    aes(x = min_x, y = min_y, xend = min_x + len_x, yend = min_y),
    arrow = arrow(length = unit(0.5, "mm"), type = "closed"),
    color = "black",
    linewidth = 0.1
  ) +
  geom_segment(
    aes(x = min_x, y = min_y, xend = min_x, yend = min_y + len_y),
    arrow = arrow(length = unit(0.5, "mm"), type = "closed"),
    color = "black",
    linewidth = 0.1
  ) +
  annotate(
    "text",
    x = min_x + len_x / 2,
    y = min_y - range_y * 0.1,
    label = "UMAP 1",
    size = 3.5 * 0.3528
  ) +
  annotate(
    "text",
    x = min_x - range_x * 0.1,
    y = min_y + len_y / 2,
    label = "UMAP 2",
    size = 3.5 * 0.35285,
    angle = 90
  ) +
  theme_void() +
  theme(legend.position = "none")

ggsave_pdf_svg(
  base_name = "05-results/02-plots/figure-1-scheme-umap",
  plot = p_umap,
  width = 10,
  height = 10,
  units = "mm"
)


# Panels b-e: Cohort description

# %% Load data
amp_metadata <- read_csv("03-data/private/00-AMP_V1_metadata.csv")
clinical_metadata <- read_csv(
  "03-data/private/00-AMP_V1_clinical_metadata.csv"
)

sample_metadata <- amp_metadata %>%
  filter(predicted_doublet == FALSE) %>%
  filter(!(leiden_0_25_corrected %in% c(18, 24))) %>%
  filter(case_control != "Other") %>%
  distinct(
    sample_id,
    participant_id,
    brain_region,
    age_at_baseline,
    case_control,
    cohort,
    education_level_years,
    ethnicity,
    hoehn_yahr_stage,
    path_braak_lb,
    path_braak_nft,
    race,
    sex,
    set
  ) %>%
  mutate(
    brain_region_long = recode(brain_region, !!!short_to_long_brain_region),
    case_control_long = factor(
      recode(case_control, !!!short_to_long_case_control),
      levels = c(
        "Parkinson's\nDisease",
        "Control",
        "Other disorder\n(Parkinsonism or\nLB Dementia)"
      )
    )
  ) %>%
  left_join(
    clinical_metadata %>% select(participant_id, PMI, RIN_PVC),
    by = "participant_id"
  )

participant_metadata <- sample_metadata %>%
  distinct(participant_id, .keep_all = TRUE)

# %% Figure 1 panel d: age of death density plot by disease status
disease_densityplot_template <- list(
  geom_density(alpha = 0.1),
  geom_rug(alpha = 0.7, linewidth = .1),
  scale_fill_manual(
    values = case_control_colors,
    name = "Disease status:",
    labels = c("Parkinson's Disease", "Control")
  ),
  scale_color_manual(
    values = case_control_colors,
    name = "Disease status:",
    labels = c("Parkinson's Disease", "Control")
  ),
  theme_publication(base_size = 5, base_family = font),
  no_legend
)

disease_by_age_densityplot <- participant_metadata %>%
  ggplot(aes(x = age_at_baseline, fill = case_control, color = case_control)) +
  disease_densityplot_template +
  geom_vline(
    data = participant_metadata %>%
      group_by(case_control) %>%
      summarise(median_value = median(age_at_baseline)),
    mapping = aes(xintercept = median_value, color = case_control),
    linetype = "dashed"
  ) +
  labs(x = "Age of death\n(years)", y = "Density")

ggsave_pdf_svg(
  base_name = "05-results/02-plots/figure-1-disease-by-age-densityplot",
  width = 25,
  height = 20,
  units = "mm",
  plot = disease_by_age_densityplot
)

# %% Figure 1 panel e: PMI density plot by disease status
disease_by_pmi_densityplot <- participant_metadata %>%
  ggplot(aes(x = PMI, fill = case_control, color = case_control)) +
  disease_densityplot_template +
  geom_vline(
    data = participant_metadata %>%
      group_by(case_control) %>%
      summarise(median_value = median(PMI, na.rm = TRUE)),
    mapping = aes(xintercept = median_value, color = case_control),
    linetype = "dashed"
  ) +
  labs(x = "Post-mortem interval\n(hours)", y = "Density")

ggsave_pdf_svg(
  base_name = "05-results/02-plots/figure-1-disease-by-pmi-densityplot",
  width = 25,
  height = 20,
  units = "mm",
  plot = disease_by_pmi_densityplot
)

# %% Figure 1 panel b: categorical cohort covariates barplot (sex x bank x disease)
categorical_cohort_covariates_barplot <- participant_metadata %>%
  add_count(sex, case_control, cohort, name = "n_donors") %>%
  add_count(cohort, name = "n_cohort_total") %>%
  mutate(cohort_n = glue("{cohort} (n={n_cohort_total})")) %>%
  ggplot(aes(x = sex, fill = case_control)) +
  geom_bar(color = "white", linewidth = .1, width = 1) +
  geom_text(
    aes(label = glue("n={n_donors}")),
    stat = "count",
    position = position_stack(vjust = .5),
    size = 1.5
  ) +
  facet_wrap(~cohort_n) +
  scale_x_discrete(labels = c("Female" = "F", "Male" = "M")) +
  scale_fill_manual(
    values = case_control_colors,
    name = "Disease status:",
    labels = c("Parkinson's Disease", "Control")
  ) +
  theme_publication(base_size = 5, base_family = font) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1)
  ) +
  no_legend

ggsave_pdf_svg(
  base_name = "05-results/02-plots/figure-1-categorical-cohort-covariates-barplot",
  width = 25,
  height = 25,
  units = "mm",
  plot = categorical_cohort_covariates_barplot
)

# %% Figure 1 panel c: sample counts by brain region and cohort
categorical_samples_covariates_barplot <- sample_metadata %>%
  add_count(case_control, cohort, brain_region, name = "n_samples") %>%
  add_count(brain_region, name = "n_brain_region_total") %>%
  mutate(
    brain_region = factor(
      brain_region,
      levels = rev(correct_brain_region_order)
    ),
    brain_region_n = glue("{brain_region}\n(n={n_brain_region_total})"),
    brain_region_n = factor(
      brain_region_n,
      levels = unique(brain_region_n[order(brain_region)])
    )
  ) %>%
  ggplot(aes(x = cohort, fill = case_control)) +
  geom_bar(color = "white", linewidth = .1, width = 1) +
  geom_text(
    aes(label = glue("n={n_samples}")),
    stat = "count",
    position = position_stack(vjust = .5),
    size = 1.5
  ) +
  facet_wrap(~brain_region_n, ncol = 1) +
  theme_publication(5, font) +
  scale_fill_manual(
    values = case_control_colors,
    name = "Disease status:",
    labels = c("Parkinson's Disease", "Control")
  ) +
  theme(axis.title = element_blank(), strip.placement = "outer") +
  no_legend

ggsave_pdf_svg(
  base_name = "05-results/02-plots/figure-1-categorical-samples-covariates-barplot",
  width = 30,
  height = 25 + 25 + 20,
  units = "mm",
  plot = categorical_samples_covariates_barplot
)

# %% Summary statistics for main text
cat("Case/Control counts:\n")
print(table(participant_metadata$case_control))
pd_control_ratio <- table(participant_metadata$case_control) /
  (table(participant_metadata$case_control))[2]
print(pd_control_ratio)

cat("\nProportion of nuclei by brain region:\n")
print(prop.table(table(amp_metadata$brain_region)) * 100)

cat("\nFisher test: sex x disease status:\n")
print(fisher.test(table(
  participant_metadata$sex,
  participant_metadata$case_control
)))

cat("\nFisher test: brain bank x disease status:\n")
print(fisher.test(
  table(participant_metadata$cohort, participant_metadata$case_control),
  simulate.p.value = TRUE
))

cat("\nFisher test: region x disease status:\n")
print(fisher.test(
  table(sample_metadata$brain_region, sample_metadata$case_control),
  simulate.p.value = TRUE
))

cat("\nWilcoxon test: age at death (PD vs Control):\n")
print(wilcox.test(age_at_baseline ~ case_control, participant_metadata))
participant_metadata %>%
  group_by(case_control) %>%
  summarise(
    median_age = median(age_at_baseline),
    iqr_age = IQR(age_at_baseline)
  ) %>%
  print()

cat("\nWilcoxon test: PMI (PD vs Control):\n")
print(wilcox.test(PMI ~ case_control, participant_metadata))
participant_metadata %>%
  group_by(case_control) %>%
  summarise(
    median_pmi = median(PMI, na.rm = TRUE),
    iqr_pmi = IQR(PMI, na.rm = TRUE)
  ) %>%
  print()


# Panel f: Cell subtype UMAP

# %% Load data
amp_metadata <- read_csv("03-data/private/00-AMP_V1_metadata.csv")
cell_annotations <- read_tsv("03-data/public/annotations_7_11_agg.txt")
clinical_metadata <- read_csv(
  "03-data/private/00-AMP_V1_clinical_metadata.csv"
)

amp_metadata <- amp_metadata %>%
  filter(predicted_doublet == FALSE) %>%
  filter(!(leiden_0_25_corrected %in% c(18, 24))) %>%
  filter(case_control != "Other") %>%
  rename(
    predicted_cluster_id = cluster_id,
    predicted_cell_type = cell_type,
    cluster_id = leiden_0_25_corrected
  ) %>%
  mutate(
    brain_region_long = recode(brain_region, !!!short_to_long_brain_region),
    case_control_long = factor(
      recode(case_control, !!!short_to_long_case_control),
      levels = correct_case_control_order
    )
  ) %>%
  left_join(cell_annotations, by = "cluster_id") %>%
  left_join(
    clinical_metadata %>% select(participant_id, PMI, RIN_PVC),
    by = "participant_id"
  ) %>%
  slice_sample(n = nrow(.))


# (UMAP helpers loaded from utils/figure_helpers.R)

# %% Figure 1 panel f: cell subtype UMAP (vibrant coloring)
cell_subtype_colors_vibrant <- amp_metadata %>%
  distinct(cell_subtype_short, cell_subtype_color_vibrant) %>%
  deframe()

panel_subtype <- ggplot(data = amp_metadata) +
  rasterise(
    geom_umap_points_mapped_color(amp_metadata, cell_subtype_short),
    dpi = dpi
  ) +
  cluster_annotation(
    cluster_centroid(amp_metadata, "cell_subtype_short"),
    mapping = aes(x = X, y = Y, label = label, color = color),
    bg.color = "white",
    bg.r = 0.10,
    fontface = "bold",
    size = 7
  ) +
  scale_color_manual(values = cell_subtype_colors_vibrant) +
  top_label_annotation(
    "Cell subtypes",
    mean(range(amp_metadata$X_umap)),
    size = 8
  ) +
  nuclei_number_annotation(nrow(amp_metadata), size = 8) +
  theme_umap_plot()

ggsave_pdf_svg(
  base_name = "05-results/02-plots/figure-1-umap-cell-subtype_vibrant",
  plot = panel_subtype,
  width = fig_width + 100,
  height = fig_height + 100,
  units = "mm"
)


# Panels g-h: Marker gene dotplots

# %% Load data
dge <- read_rds("04-analysis/private/pseudobulk-dge-cell-types.rds")

# Arrange samples by cell type for consistent ordering
metadata <- as_tibble(dge$samples, rownames = "sample_name") %>%
  arrange(cell_type_short)
dge$counts <- dge$counts[, metadata$sample_name]
dge$samples <- dge$samples[metadata$sample_name, ]
lcpm <- cpm(dge, log = TRUE, prior.count = 0.25)

qlf <- read_rds("04-analysis/private/pseudobulk-dea-cell-types.rds")
qlf_neurons_vs_other_neurons <- read_rds(
  "04-analysis/private/pseudobulk-dea-neurons-vs-other-neurons.rds"
)
qlf_neurons_vs_other <- read_rds(
  "04-analysis/private/pseudobulk-dea-cell-supertypes.rds"
)

cell_annotations <- read_tsv("03-data/public/annotations_7_11_agg.txt")
gene_annotations <- read_csv("03-data/public/gene-metadata.csv") %>%
  filter(transcript_is_canonical == 1) %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

gene_pct_cells_cell_type <- read_csv(
  "03-data/public/gene-pct-cells-expressing-cell_type_short.csv"
)
gene_pct_long_cell_type <- gene_pct_cells_cell_type %>%
  pivot_longer(
    cols = starts_with("pct_detected"),
    names_to = "cell_type_short",
    values_to = "pct_detected"
  ) %>%
  mutate(cell_type_short = str_extract(cell_type_short, "[a-zA-Z-/]+$"))

gene_pct_cells_cell_supertype <- read_csv(
  "03-data/public/gene-pct-cells-expressing-cell-supertype-short.csv"
)
gene_pct_long_cell_supertype <- gene_pct_cells_cell_supertype %>%
  pivot_longer(
    cols = starts_with("pct_detected"),
    names_to = "cell_supertype_short",
    values_to = "pct_detected"
  ) %>%
  mutate(
    cell_supertype_short = str_extract(cell_supertype_short, "[a-zA-Z-/]+$")
  )


# %% Figure 1 panel g: dotplot of data-driven marker genes

# Select top 5 DEGs per cell type (logFC >= 2, detected in >= 50% of cells)
top <- 5

extract_top_markers <- function(
  qlf_obj,
  pct_data,
  pct_join_col,
  filter_types = NULL
) {
  qlf_obj %>%
    map(~ topTags(.x, n = nrow(.x), adjust.method = "bonferroni")) %>%
    map("table") %>%
    map(~ as_tibble(.x, rownames = "ensembl_gene_id")) %>%
    map(~ select(.x, -c(hgnc_symbol, description))) %>%
    map(~ left_join(.x, gene_annotations, by = "ensembl_gene_id")) %>%
    imap(function(x, idx) {
      x %>%
        mutate(!!pct_join_col := idx) %>%
        left_join(
          pct_data %>% filter(.data[[pct_join_col]] == idx),
          by = c("ensembl_gene_id", pct_join_col)
        )
    }) %>%
    map(~ drop_na(.x, hgnc_symbol)) %>%
    map(~ dplyr::filter(.x, logFC >= 2)) %>%
    map(~ dplyr::filter(.x, pct_detected >= 0.5)) %>%
    map(~ slice_head(.x, n = top)) %>%
    map(~ pull(.x, ensembl_gene_id)) %>%
    enframe(name = pct_join_col, value = "ensembl_gene_id") %>%
    unnest(ensembl_gene_id) %>%
    {
      if (!is.null(filter_types)) {
        filter(., .data[[pct_join_col]] %in% filter_types)
      } else {
        .
      }
    }
}

marker_genes_general <- extract_top_markers(
  qlf,
  gene_pct_long_cell_type,
  "cell_type_short",
  filter_types = setdiff(correct_cell_type_order, c("Exc", "Inh", "DMNX/GPI"))
)

# Invert filter: keep only Exc and Inh types
marker_genes_neurons_vs_other_neurons <- extract_top_markers(
  qlf_neurons_vs_other_neurons,
  gene_pct_long_cell_type,
  "cell_type_short"
)

marker_genes_neurons_vs_other <- extract_top_markers(
  qlf_neurons_vs_other,
  gene_pct_long_cell_supertype,
  "cell_type_short",
  filter_types = "Neu"
)
# Fix: cell_type_short here refers to cell_supertype_short in the original; rename to match
marker_genes_neurons_vs_other <- marker_genes_neurons_vs_other %>%
  rename(cell_type_short = cell_type_short)

all_markers <- bind_rows(
  marker_genes_general,
  marker_genes_neurons_vs_other_neurons,
  marker_genes_neurons_vs_other
) %>%
  distinct(ensembl_gene_id, .keep_all = TRUE) %>%
  left_join(
    gene_annotations %>% select(ensembl_gene_id, hgnc_symbol),
    by = "ensembl_gene_id"
  )

correct_gene_order <- all_markers %>%
  mutate(
    cell_type_short = factor(cell_type_short, levels = correct_cell_type_order)
  ) %>%
  arrange(cell_type_short) %>%
  pull(hgnc_symbol) %>%
  rev()

top_markers_lcpm <- lcpm[all_markers$ensembl_gene_id, ]
rownames(top_markers_lcpm) <- deframe(gene_annotations[, 1:2])[rownames(
  top_markers_lcpm
)]

top_markers_lcpm_plot <- top_markers_lcpm %>%
  as.data.frame() %>%
  rownames_to_column("hgnc_symbol") %>%
  pivot_longer(-hgnc_symbol, names_to = "sample", values_to = "expression") %>%
  mutate(cell_type_short = str_extract(sample, "[a-zA-Z-/]+$")) %>%
  group_by(hgnc_symbol, cell_type_short) %>%
  summarise(mean_expr = mean(expression), .groups = "drop") %>%
  left_join(
    gene_annotations %>% select(ensembl_gene_id, hgnc_symbol),
    by = "hgnc_symbol"
  ) %>%
  left_join(
    gene_pct_long_cell_type,
    by = c("ensembl_gene_id", "cell_type_short")
  ) %>%
  group_by(hgnc_symbol) %>%
  mutate(scaled_expr = scale(mean_expr)[, 1]) %>%
  ungroup() %>%
  left_join(
    cell_annotations %>%
      select(starts_with("cell_type")) %>%
      distinct(cell_type_short, .keep_all = TRUE),
    by = "cell_type_short"
  ) %>%
  mutate(
    hgnc_symbol = factor(hgnc_symbol, levels = correct_gene_order),
    cell_type_short = factor(cell_type_short, levels = correct_cell_type_order)
  )

# Correct cell type and gene axis colors
correct_cell_colors_data <- top_markers_lcpm_plot %>%
  distinct(cell_type_short, cell_type_color_vibrant) %>%
  deframe()
correct_cell_colors_data <- correct_cell_colors_data[
  correct_cell_type_order[correct_cell_type_order != "Neu"]
]

correct_gene_colors_data <- all_markers %>%
  left_join(
    top_markers_lcpm_plot %>%
      distinct(cell_type_short, cell_type_color_vibrant),
    by = "cell_type_short"
  ) %>%
  replace_na(list(cell_type_color_vibrant = "#da8bc3")) %>%
  select(hgnc_symbol, cell_type_color_vibrant) %>%
  deframe()
correct_gene_colors_data <- correct_gene_colors_data[correct_gene_order]

ggplot(top_markers_lcpm_plot, aes(x = cell_type_short, y = hgnc_symbol)) +
  geom_point(aes(size = pct_detected, color = scaled_expr)) +
  scale_color_distiller(
    palette = "RdBu",
    limits = c(-2, 2),
    oob = scales::squish,
    name = "Average scaled\nexpression",
    guide = guide_colorbar(theme = theme(legend.title.position = "top"))
  ) +
  scale_size(
    range = c(0.1, 1.6),
    name = "Proportion of\ncells expressing",
    labels = scales::percent_format(accuracy = 1),
    guide = guide_legend(
      theme = theme(
        legend.title.position = "top",
        legend.text.position = "bottom"
      )
    )
  ) +
  theme_publication(base_size = 5, base_family = font) +
  horizontal_scientific_colormap +
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      colour = correct_cell_colors_data
    ),
    axis.text.y = element_text(colour = correct_gene_colors_data),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(hjust = 0.5),
    legend.key.width = unit(.5, "lines"),
    legend.box = "vertical"
  )

ggsave_pdf_svg(
  base_name = "05-results/02-plots/figure-1-dotplot-data-driven-cell-types-marker-genes-color",
  width = 30,
  height = 120,
  units = "mm"
)


# %% Figure 1 panel h: dotplot of literature-driven marker genes
gene_map <- setNames(
  gene_annotations$ensembl_gene_id,
  gene_annotations$hgnc_symbol
)

markers_dict_hgnc <- c(
  # Astrocytes
  "GFAP" = "Astro",
  "ALDH1L1" = "Astro",
  "AQP4" = "Astro",
  "SLC1A2" = "Astro",
  # OPCs
  "PDGFRA" = "OPC",
  "CSPG4" = "OPC",
  "VCAN" = "OPC",
  "PTPRZ1" = "OPC",
  # Oligodendrocytes
  "MBP" = "Oligo",
  "PLP1" = "Oligo",
  "MOG" = "Oligo",
  "OPALIN" = "Oligo",
  "CNP" = "Oligo",
  # Ependymal
  "FOXJ1" = "Epen",
  "CFAP126" = "Epen",
  # Microglia
  "AIF1" = "Micro",
  "ITGAM" = "Micro",
  "P2RY12" = "Micro",
  "TMEM119" = "Micro",
  "C1QA" = "Micro",
  # CTL (cytotoxic T lymphocytes)
  "CD3D" = "CTL",
  "CD8A" = "CTL",
  "GZMB" = "CTL",
  "PRF1" = "CTL",
  "NKG7" = "CTL",
  # Endothelial
  "PECAM1" = "Endo",
  "CLDN5" = "Endo",
  "FLT1" = "Endo",
  # Pericytes
  "PDGFRB" = "Peri",
  "RGS5" = "Peri",
  "KCNJ8" = "Peri",
  "COL1A1" = "Peri",
  "LUM" = "Peri",
  "DCN" = "Peri",
  # Pan-neuronal
  "RBFOX3" = "Neu",
  "SYT1" = "Neu",
  "SNAP25" = "Neu",
  # Excitatory neurons
  "SLC17A7" = "Exc",
  "CAMK2A" = "Exc",
  "SATB2" = "Exc",
  # Inhibitory neurons
  "GAD1" = "Inh",
  "GAD2" = "Inh",
  "SLC32A1" = "Inh"
)

markers_list_ensembl <- gene_map[names(markers_dict_hgnc)]
markers_lcpm <- lcpm[markers_list_ensembl, ]
rownames(markers_lcpm) <- deframe(gene_annotations[, 1:2])[rownames(
  markers_lcpm
)]

correct_gene_order_lit <- rev(unique(names(markers_dict_hgnc)))

plot_df <- markers_lcpm %>%
  as.data.frame() %>%
  rownames_to_column("hgnc_symbol") %>%
  pivot_longer(-hgnc_symbol, names_to = "sample", values_to = "expression") %>%
  mutate(cell_type_short = str_extract(sample, "[a-zA-Z-/]+$")) %>%
  group_by(hgnc_symbol, cell_type_short) %>%
  summarise(mean_expr = mean(expression), .groups = "drop") %>%
  left_join(
    gene_annotations %>% select(ensembl_gene_id, hgnc_symbol),
    by = "hgnc_symbol"
  ) %>%
  left_join(
    gene_pct_long_cell_type,
    by = c("ensembl_gene_id", "cell_type_short")
  ) %>%
  group_by(hgnc_symbol) %>%
  mutate(scaled_expr = scale(mean_expr)[, 1]) %>%
  ungroup() %>%
  left_join(
    cell_annotations %>% select(starts_with("cell_type")),
    by = "cell_type_short"
  ) %>%
  mutate(
    hgnc_symbol = factor(hgnc_symbol, levels = correct_gene_order_lit),
    cell_type_short = factor(cell_type_short, levels = correct_cell_type_order)
  )

correct_cell_colors_lit <- plot_df %>%
  distinct(cell_type_short, cell_type_color_vibrant) %>%
  deframe()
correct_cell_colors_lit <- correct_cell_colors_lit[
  correct_cell_type_order[correct_cell_type_order != "Neu"]
]

correct_gene_colors_lit <- tibble(
  hgnc_symbol = names(markers_dict_hgnc),
  cell_type_short = unname(markers_dict_hgnc)
) %>%
  left_join(
    plot_df %>% distinct(cell_type_short, cell_type_color_vibrant),
    by = "cell_type_short"
  ) %>%
  replace_na(list(cell_type_color_vibrant = "#da8bc3")) %>%
  select(hgnc_symbol, cell_type_color_vibrant) %>%
  deframe()
correct_gene_colors_lit <- correct_gene_colors_lit[correct_gene_order_lit]

ggplot(plot_df, aes(x = cell_type_short, y = hgnc_symbol)) +
  geom_point(aes(size = pct_detected, color = scaled_expr)) +
  scale_color_distiller(
    palette = "RdBu",
    limits = c(-2, 2),
    oob = scales::squish,
    name = "Average scaled\nexpression",
    guide = guide_colorbar(theme = theme(legend.title.position = "top"))
  ) +
  scale_size(
    range = c(0.1, 1.6),
    name = "Proportion of\ncells expressing",
    labels = scales::percent_format(accuracy = 1),
    guide = guide_legend(
      theme = theme(
        legend.title.position = "top",
        legend.text.position = "bottom"
      )
    )
  ) +
  theme_publication(base_size = 5, base_family = font) +
  horizontal_scientific_colormap +
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      colour = correct_cell_colors_lit
    ),
    axis.text.y = element_text(colour = correct_gene_colors_lit),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(hjust = 0.5),
    legend.key.width = unit(.5, "lines"),
    legend.box = "vertical"
  )

ggsave_pdf_svg(
  base_name = "05-results/02-plots/figure-1-dotplot-literature-driven-cell-types-marker-genes-color",
  width = 30,
  height = 120,
  units = "mm"
)


# Panels i-j: Cell composition barplots

# %% Load data
amp_metadata <- read_csv("03-data/private/00-AMP_V1_metadata.csv")
clinical_metadata <- read_csv(
  "03-data/private/00-AMP_V1_clinical_metadata.csv"
)
cell_annotations <- read_tsv("03-data/public/annotations_7_11_agg.txt")

amp_metadata <- amp_metadata %>%
  filter(predicted_doublet == FALSE) %>%
  filter(!(leiden_0_25_corrected %in% c(18, 24))) %>%
  filter(case_control != "Other") %>%
  rename(
    predicted_cluster_id = cluster_id,
    predicted_cell_type = cell_type,
    cluster_id = leiden_0_25_corrected
  ) %>%
  mutate(
    brain_region_long = recode(brain_region, !!!short_to_long_brain_region)
  ) %>%
  left_join(cell_annotations, by = "cluster_id") %>%
  left_join(
    clinical_metadata %>% select(participant_id, PMI, RIN_PVC),
    by = "participant_id"
  ) %>%
  slice_sample(n = nrow(.))


# %% Figure 1 panel j: horizontal stacked barplot — cell type by brain region
region_props_type <- amp_metadata %>%
  count(brain_region, cell_type_short, name = "n_cells") %>%
  group_by(brain_region) %>%
  mutate(
    total_cells_region = sum(n_cells),
    prop = n_cells / total_cells_region
  ) %>%
  ungroup()

global_props_type <- amp_metadata %>%
  count(cell_type_short, name = "n_cells") %>%
  mutate(
    brain_region = "All",
    total_cells_region = sum(n_cells),
    prop = n_cells / total_cells_region
  )

threshold_type <- 0.03

final_props_type <- bind_rows(region_props_type, global_props_type) %>%
  mutate(
    cell_type_grouped = if_else(prop < threshold_type, "Other", cell_type_short)
  ) %>%
  group_by(brain_region, cell_type_grouped) %>%
  summarise(prop = sum(prop), n_cells = sum(n_cells), .groups = "drop") %>%
  mutate(
    brain_region = factor(
      brain_region,
      levels = c("DMNX", "GPI", "PFC", "PMC", "PVC", "", "All")
    ),
    cell_type_grouped = factor(
      cell_type_grouped,
      levels = rev(c(correct_cell_type_order, "Other"))
    )
  )

all_cell_type_colors <- c(correct_cell_type_colors, Other = cool_grey)

ggplot(
  final_props_type,
  aes(x = prop, y = brain_region, fill = cell_type_grouped)
) +
  geom_col(position = "stack", width = 1, color = "white", linewidth = 0.01) +
  geom_text(
    aes(
      label = if_else(
        prop > threshold_type,
        scales::percent(prop, accuracy = 1),
        ""
      )
    ),
    position = position_stack(vjust = .5),
    size = 5 * 5 / 14
  ) +
  scale_fill_manual(
    values = all_cell_type_colors,
    name = "Cell type:",
    guide = guide_legend(reverse = TRUE)
  ) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0)),
    limits = c(0, 1)
  ) +
  scale_y_discrete(drop = FALSE, breaks = function(levels) {
    levels[levels != ""]
  }) +
  labs(x = "Cell proportion", y = "Brain region") +
  theme_publication(base_size = 5, base_family = font)

ggsave_pdf_svg(
  base_name = "05-results/02-plots/figure-1-horizontal-stacked-barplot-cell-type-regions",
  width = 90,
  height = 25,
  units = "mm"
)


# %% Figure 1 panel i: multiple horizontal stacked barplots (donor-normalized)

stacked_barplots_df <- amp_metadata %>%
  filter(case_control != "Other") %>%
  mutate(
    cell_subtype_short = factor(
      cell_subtype_short,
      levels = rev(correct_cell_subtype_order)
    ),
    brain_region = factor(
      brain_region,
      levels = rev(correct_brain_region_order)
    )
  )

# Nested y-axis key for cell supertype groupings
nested_keys <- key_range_manual(
  start = c(
    "Astro 1",
    "Oligo 1",
    "Micro 1",
    "Exc L2/3 IT",
    "Inh LAMP5",
    "Astro 1",
    "Micro 1",
    "Endo",
    "Exc L2/3 IT"
  ),
  end = c(
    "Astro 3",
    "Oligo 2",
    "Micro 2",
    "Exc L6/6b CT",
    "Inh SST",
    "Epen",
    "CTL",
    "Peri",
    "DMNX/GPI"
  ),
  name = c(
    "Astrocyte",
    "Oligodendrocyte",
    "Microglia",
    "Excitatory neuron\n(glutamatergic)",
    "Inhibitory neuron\n(GABAergic)",
    "Non-neuron",
    "Immune",
    "Vascular",
    "Neuron"
  ),
  level = c(1, 1, 1, 1, 1, 2, 2, 2, 2)
)

base_stacked_barplot <- ggplot(
  stacked_barplots_df,
  aes(y = cell_subtype_short)
) +
  geom_bar(position = "fill", width = 1, color = "white", linewidth = 0.1) +
  scale_x_continuous(
    breaks = c(0, .5, 1),
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0)),
    limits = c(0, 1)
  ) +
  theme_publication(base_size = 5, base_family = font) +
  theme(legend.position = "bottom")

nuclei_by_region_barplot <- base_stacked_barplot +
  aes(fill = brain_region) +
  scale_y_discrete(guide = guide_axis_nested(key = nested_keys)) +
  labs(x = "Nuclei by\nbrain region") +
  scale_fill_manual(values = region_colors, name = "Brain region of origin:") +
  guides(
    fill = guide_legend(
      ncol = 1,
      direction = "horizontal",
      title.position = "top"
    )
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
    legend.key.spacing = unit(.7, "mm"),
    legend.background = element_rect(fill = "white", color = NA),
    legend.box.margin = margin(0, 0, 0, 0)
  ) +
  theme_guide(bracket = element_line(linewidth = .1))

# Donor-normalized disease barplot
# (1. mean nuclei counts per donor within each subtype x disease group;
#  2. relative abundance of means — corrects for donor count imbalance between groups)
nuclei_by_disease_barplot_norm <- ggplot(
  stacked_barplots_df %>%
    group_by(cell_subtype_short, case_control, participant_id) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(cell_subtype_short, case_control) %>%
    summarise(mean_n = mean(n), .groups = "drop") %>%
    group_by(cell_subtype_short) %>%
    mutate(
      prop = mean_n / sum(mean_n),
      cell_subtype_short = factor(
        cell_subtype_short,
        levels = rev(correct_cell_subtype_order)
      )
    ),
  aes(y = cell_subtype_short, x = prop, fill = case_control)
) +
  geom_col(width = 1, color = "white", linewidth = 0.1) +
  scale_x_continuous(
    breaks = c(0, .5, 1),
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0)),
    limits = c(0, 1)
  ) +
  scale_fill_manual(
    values = case_control_colors,
    name = "Disease status:",
    labels = c("Parkinson's Disease", "Control")
  ) +
  labs(x = "Nuclei by\ndisease status") +
  guides(
    fill = guide_legend(
      ncol = 1,
      direction = "horizontal",
      title.position = "top"
    )
  ) +
  theme_publication(base_size = 5, base_family = font) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
    plot.margin = margin(0, 5, 0, 5)
  )

nuclei_by_subtype_barplot <- ggplot(
  stacked_barplots_df,
  aes(y = cell_subtype_short, fill = cell_subtype_short)
) +
  geom_bar(width = 1, color = "white", linewidth = 0.1) +
  scale_x_continuous(labels = scales::label_number(scale = 1 / 1000)) +
  labs(x = "Total number of<br>nuclei (\u00d710\u00b3)") +
  scale_fill_manual(
    values = correct_cell_subtype_colors,
    name = "Cell subtypes: "
  ) +
  guides(
    fill = guide_legend(
      ncol = 2,
      direction = "horizontal",
      title.position = "top",
      reverse = TRUE
    )
  ) +
  theme_publication(base_size = 5, base_family = font) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = ggtext::element_markdown(),
    legend.key.spacing = unit(.7, "mm"),
    legend.background = element_rect(fill = "white", color = NA),
    legend.box.margin = margin(0, 0, 0, 0)
  )

composed_horizontal_stacked_barplots_norm <-
  (nuclei_by_region_barplot + no_legend) +
  (nuclei_by_disease_barplot_norm + no_legend) +
  (nuclei_by_subtype_barplot + no_legend) +
  plot_layout(widths = c(.9, .75, .75))

ggsave_pdf_svg(
  base_name = "05-results/02-plots/figure-1-multiple-horizontal-stacked-barplots-normalized",
  width = 90,
  height = 55,
  units = "mm",
  plot = composed_horizontal_stacked_barplots_norm
)
