# Description: Generate Supplementary Figure 1 panels:
#   panel a — cohort covariates heatmap (ComplexHeatmap)
#   panel b — sequencing batch UMAP
#   panel c — brain bank UMAP
#   panel d — PMI UMAP
#   panel e — donor UMAP
#   panel f — disease status UMAP
#   panel g — brain region UMAP

# %% Setup
library(tidyverse)
library(glue)
library(svglite)
library(systemfonts)
library(ComplexHeatmap)
library(sysfonts)
library(ggrastr)
library(ggrepel)
library(scales)
library(conflicted)
source("utils/my_utils.R")
source("utils/figure_helpers.R")
conflicts_prefer(dplyr::filter, dplyr::select)
set.seed(42)


# %% Load data
amp_metadata <- read_csv("03-data/private/00-AMP_V1_metadata.csv")
cell_annotations <- read_tsv("03-data/public/annotations_7_11_agg.txt")
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


# %% Supplementary figure 1 panel a: cohort metadata heatmap
target_size <- 5
target_lwd <- 0.5
cell_anno_size <- unit(1.5, "mm")

text_gp <- gpar(fontfamily = font, fontsize = target_size)
border_gp <- gpar(lwd = target_lwd, col = "black")

common_legend_param <- list(
  title_gp = text_gp,
  labels_gp = text_gp,
  grid_height = unit(2, "mm"),
  grid_width = unit(2, "mm"),
  title_position = "topcenter",
  legend_gp = gpar(lwd = target_lwd),
  tick_length = unit(-.25, "mm")
)

annotation_df <- participant_metadata %>%
  arrange(
    desc(case_control),
    !is.na(path_braak_lb),
    path_braak_lb,
    !is.na(path_braak_nft),
    path_braak_nft,
    !is.na(hoehn_yahr_stage),
    hoehn_yahr_stage,
    age_at_baseline
  )

col_fun_lb <- circlize::colorRamp2(
  breaks = seq(
    min(annotation_df$path_braak_lb, na.rm = TRUE),
    max(annotation_df$path_braak_lb, na.rm = TRUE),
    length.out = length(unique(na.omit(annotation_df$path_braak_lb)))
  ),
  colors = RColorBrewer::brewer.pal(
    length(unique(na.omit(annotation_df$path_braak_lb))),
    "Reds"
  )
)
col_fun_nft <- circlize::colorRamp2(
  breaks = seq(
    min(annotation_df$path_braak_nft, na.rm = TRUE),
    max(annotation_df$path_braak_nft, na.rm = TRUE),
    length.out = length(unique(na.omit(annotation_df$path_braak_nft)))
  ),
  colors = RColorBrewer::brewer.pal(
    length(unique(na.omit(annotation_df$path_braak_nft))),
    "Blues"
  )
)
col_fun_hy <- circlize::colorRamp2(
  breaks = seq(
    min(annotation_df$hoehn_yahr_stage, na.rm = TRUE),
    max(annotation_df$hoehn_yahr_stage, na.rm = TRUE),
    length.out = length(unique(na.omit(annotation_df$hoehn_yahr_stage)))
  ),
  colors = RColorBrewer::brewer.pal(
    length(unique(na.omit(annotation_df$hoehn_yahr_stage))),
    "Purples"
  )
)

hoehn_yahr_pch <- ifelse(is.na(annotation_df$hoehn_yahr_stage), "?", NA)
lb_braak_pch <- ifelse(is.na(annotation_df$path_braak_lb), "?", NA)
nft_braak_pch <- ifelse(is.na(annotation_df$path_braak_nft), "?", NA)

participant_metadata_annotations <- HeatmapAnnotation(
  `Disease status` = annotation_df$case_control,
  `LB Braak stage` = anno_simple(
    annotation_df$path_braak_lb,
    col = col_fun_lb,
    gp = border_gp,
    simple_anno_size = cell_anno_size,
    na_col = cool_grey,
    pch = lb_braak_pch,
    pt_size = unit(3, "pt")
  ),
  `NFT Braak stage` = anno_simple(
    annotation_df$path_braak_nft,
    col = col_fun_nft,
    gp = border_gp,
    simple_anno_size = cell_anno_size,
    na_col = cool_grey,
    pch = nft_braak_pch,
    pt_size = unit(3, "pt")
  ),
  `Hoehn-Yahr stage` = anno_simple(
    annotation_df$hoehn_yahr_stage,
    col = col_fun_hy,
    gp = border_gp,
    simple_anno_size = cell_anno_size,
    na_col = cool_grey,
    pch = hoehn_yahr_pch,
    pt_size = unit(3, "pt")
  ),
  `Age of death` = annotation_df$age_at_baseline,
  `Post-mortem interval` = annotation_df$PMI,
  `Sex` = annotation_df$sex,
  `Brain bank` = annotation_df$cohort,
  col = list(
    `Disease status` = case_control_colors,
    `Age of death` = circlize::colorRamp2(
      seq(
        min(annotation_df$age_at_baseline, na.rm = TRUE),
        max(annotation_df$age_at_baseline, na.rm = TRUE),
        length.out = 9
      ),
      RColorBrewer::brewer.pal(9, "Greys")
    ),
    `Post-mortem interval` = circlize::colorRamp2(
      seq(
        min(annotation_df$PMI, na.rm = TRUE),
        max(annotation_df$PMI, na.rm = TRUE),
        length.out = 9
      ),
      rev(RColorBrewer::brewer.pal(9, "RdYlGn"))
    ),
    `Sex` = sex_colors,
    `Brain bank` = brain_bank_colors
  ),
  annotation_legend_param = list(
    `Post-mortem interval` = c(
      common_legend_param,
      list(
        title = "Post-mortem interval (hours)",
        direction = "horizontal",
        border = "black",
        legend_width = unit(20, "mm")
      )
    ),
    `Age of death` = c(
      common_legend_param,
      list(
        title = "Age of death (years)",
        direction = "horizontal",
        border = "black",
        legend_width = unit(20, "mm")
      )
    ),
    `Disease status` = c(
      common_legend_param,
      list(
        title = "Disease status",
        labels = c("Parkinson's Disease", "Control")
      )
    ),
    `Sex` = c(common_legend_param, list(title = "Sex")),
    `Brain bank` = c(common_legend_param, list(title = "Brain bank", nrow = 2))
  ),
  na_col = cool_grey,
  gp = border_gp,
  annotation_name_gp = text_gp,
  simple_anno_size = cell_anno_size
)

make_styled_legend <- function(title, col_fun) {
  Legend(
    border = "black",
    direction = "horizontal",
    title = title,
    title_gp = text_gp,
    labels_gp = text_gp,
    title_position = "topcenter",
    legend_width = unit(15, "mm"),
    grid_height = unit(2, "mm"),
    grid_width = unit(2, "mm"),
    legend_gp = gpar(lwd = target_lwd),
    tick_length = unit(-.25, "mm"),
    col_fun = col_fun
  )
}
manual_legends <- list(
  make_styled_legend("LB Braak stage", col_fun_lb),
  make_styled_legend("NFT Braak stage", col_fun_nft),
  make_styled_legend("Hoehn-Yahr stage", col_fun_hy)
)

participant_metadata_heatmap <- Heatmap(
  matrix(nrow = 0, ncol = nrow(annotation_df)),
  top_annotation = participant_metadata_annotations,
  show_heatmap_legend = FALSE
)

svglite::svglite(
  filename = "05-results/02-plots/figure-1-cohort-covariates-heatmap.svg",
  width = 120 / 25.4,
  height = 40 / 25.4,
  fix_text_size = FALSE
)
draw(
  participant_metadata_heatmap,
  annotation_legend_list = manual_legends,
  annotation_legend_side = "bottom"
)
dev.off()


# (UMAP helpers loaded from utils/figure_helpers.R)

# %% Supplementary figure 1 panel g: brain region UMAP
umap_colored_regions <- ggplot(amp_metadata, aes(X_umap, Y_umap)) +
  rasterise(
    geom_point(
      aes(color = brain_region),
      shape = 16,
      size = 0.1,
      stroke = 0,
      alpha = 0.3
    ),
    dpi = dpi
  ) +
  scale_color_manual(values = region_colors) +
  guides(
    color = guide_legend(
      title = "Brain region of origin:",
      override.aes = list(alpha = 1, shape = 15, size = 3),
      nrow = 2
    )
  ) +
  theme_publication(6, font) +
  theme(
    legend.position = "bottom",
    legend.direction = "vertical",
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

ggsave(
  filename = "05-results/02-plots/figure-1-umap-brain-region.pdf",
  plot = umap_colored_regions,
  width = 90,
  height = 90,
  units = "mm"
)


# %% Supplementary figure 1 panel f: disease status UMAP
umap_colored_case_control <- ggplot(amp_metadata, aes(X_umap, Y_umap)) +
  rasterise(
    geom_point(
      aes(color = case_control),
      shape = 16,
      size = 0.1,
      stroke = 0,
      alpha = 0.3
    ),
    dpi = dpi
  ) +
  scale_color_manual(
    values = case_control_colors,
    labels = c("Case" = "Parkinson's Disease")
  ) +
  guides(
    color = guide_legend(
      title = "Disease status:",
      override.aes = list(alpha = 1, shape = 15, size = 3)
    )
  ) +
  theme_publication(6, font) +
  theme(
    legend.position = "bottom",
    legend.direction = "vertical",
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

ggsave(
  filename = "05-results/02-plots/figure-1-umap-case-control.pdf",
  plot = umap_colored_case_control,
  width = 90,
  height = 90,
  units = "mm"
)


# %% Supplementary figure 1 panel c: brain bank UMAP
umap_colored_brain_bank <- ggplot(amp_metadata, aes(X_umap, Y_umap)) +
  rasterise(
    geom_point(
      aes(color = cohort),
      shape = 16,
      size = 0.1,
      stroke = 0,
      alpha = 0.3
    ),
    dpi = dpi
  ) +
  scale_color_manual(values = brain_bank_colors) +
  guides(
    color = guide_legend(
      title = "Brain bank:",
      override.aes = list(alpha = 1, shape = 15, size = 3),
      ncol = 2
    )
  ) +
  theme_publication(6, font) +
  theme(
    legend.position = "bottom",
    legend.direction = "vertical",
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

ggsave(
  filename = "05-results/02-plots/figure-1-umap-brain-bank.pdf",
  plot = umap_colored_brain_bank,
  width = 90,
  height = 90,
  units = "mm"
)


# %% Supplementary figure 1 panel e: donor UMAP
set.seed(123)
donor_palette <- Polychrome::createPalette(
  length(unique(amp_metadata$participant_id)),
  seedcolors = "#FF0D45"
)
names(donor_palette) <- NULL

umap_colored_donor <- ggplot(
  amp_metadata %>%
    mutate(
      p_id = paste0("P", as.integer(factor(amp_metadata$participant_id))),
      p_id = fct_reorder(p_id, as.numeric(str_extract(p_id, "\\d+")))
    ),
  aes(X_umap, Y_umap)
) +
  rasterise(
    geom_point(
      aes(color = p_id),
      shape = 16,
      size = 0.1,
      stroke = 0,
      alpha = 0.3
    ),
    dpi = dpi
  ) +
  scale_color_manual(values = donor_palette) +
  guides(
    color = guide_legend(
      title = "Participant ID:",
      override.aes = list(alpha = 1, shape = 15, size = 2),
      ncol = 3
    )
  ) +
  theme_publication(6, font) +
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

ggsave(
  filename = "05-results/02-plots/figure-1-umap-donor.pdf",
  plot = umap_colored_donor,
  width = 120,
  height = 90,
  units = "mm"
)


# %% Supplementary figure 1 panel b: sequencing batch UMAP
set.seed(123)
seq_batch_palette <- Polychrome::createPalette(
  length(unique(amp_metadata$set)),
  seedcolors = "#FF0D45"
)
names(seq_batch_palette) <- NULL

umap_colored_seq_batch <- ggplot(
  amp_metadata %>%
    mutate(set = fct_reorder(set, as.numeric(str_extract(set, "\\d+")))),
  aes(X_umap, Y_umap)
) +
  rasterise(
    geom_point(
      aes(color = set),
      shape = 16,
      size = 0.1,
      stroke = 0,
      alpha = 0.3
    ),
    dpi = dpi
  ) +
  scale_color_manual(values = seq_batch_palette) +
  guides(
    color = guide_legend(
      title = "Sequencing batch:",
      override.aes = list(alpha = 1, shape = 15, size = 2),
      ncol = 3
    )
  ) +
  theme_publication(6, font) +
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

ggsave(
  filename = "05-results/02-plots/figure-1-umap-seq-batch.pdf",
  plot = umap_colored_seq_batch,
  width = 120,
  height = 90,
  units = "mm"
)


# %% Supplementary figure 1 panel d: PMI UMAP
umap_colored_PMI <- ggplot(amp_metadata, aes(X_umap, Y_umap)) +
  rasterise(
    geom_point(
      aes(color = PMI),
      shape = 16,
      size = 0.1,
      stroke = 0,
      alpha = 0.3
    ),
    dpi = dpi
  ) +
  scale_color_viridis_c(direction = -1) +
  theme_publication(6, font) +
  horizontal_scientific_colormap +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

ggsave(
  filename = "05-results/02-plots/figure-1-umap-PMI.pdf",
  plot = umap_colored_PMI,
  width = 90,
  height = 90,
  units = "mm"
)
