# Description: Generate Supplementary Figure 2 panels:
#   panels a-e — jitter + logFC bars per brain region
#   panel f    — faceted cell proportion boxplot
#   panels a-e (credibility) — credibility heatmaps per brain region

# %% Load libraries
library(tidyverse)
library(patchwork)
library(ggnewscale)
library(ggtree)
library(glue)
library(aplot)
library(ggtext)
source("utils/my_utils.R")

# %% Load data
coda_results_path <- "04-analysis/public/02-differential-abundance-phi5-itermore-peri-final"

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
  )

subtype_region_props <- amp_metadata %>%
  group_by(cell_subtype_short, brain_region) %>%
  summarise(total_cells_in_subtype_region = n(), .groups = "drop") %>%
  group_by(cell_subtype_short) %>%
  mutate(
    subtype_regional_distribution = total_cells_in_subtype_region /
      sum(total_cells_in_subtype_region)
  ) %>%
  ungroup()

real_counts <- amp_metadata %>%
  count(sample_id, cell_subtype_short, name = "cells_count")

sample_metadata <- amp_metadata %>%
  distinct(
    sample_id,
    participant_id,
    case_control,
    brain_region,
    brain_region_long,
    cohort,
    age_at_baseline,
    sex
  )

MIN_PROP <- 0.05
cell_counts_per_sample <- sample_metadata %>%
  crossing(cell_subtype_short = unique(amp_metadata$cell_subtype_short)) %>%
  left_join(real_counts, by = c("sample_id", "cell_subtype_short")) %>%
  mutate(cells_count = replace_na(cells_count, 0)) %>%
  left_join(
    subtype_region_props,
    by = c("cell_subtype_short", "brain_region")
  ) %>%
  filter(subtype_regional_distribution >= MIN_PROP) %>%
  group_by(sample_id) %>%
  mutate(prop_ct_in_sample = cells_count / sum(cells_count)) %>%
  ungroup() %>%
  group_by(cell_subtype_short, brain_region) %>%
  mutate(donor_contribution_to_ct_pool = cells_count / sum(cells_count)) %>%
  ungroup()

regions_map <- cell_counts_per_sample %>%
  ungroup() %>%
  distinct(brain_region, brain_region_long)

effect_files <- list.files(
  coda_results_path,
  pattern = ".*effects.*\\.csv$",
  full.names = TRUE
)
node_files <- list.files(
  coda_results_path,
  pattern = ".*node.*\\.csv$",
  full.names = TRUE
)
credibility_files <- list.files(
  coda_results_path,
  pattern = ".*credibility.*\\.csv$",
  full.names = TRUE
)
model_files <- list.files(
  coda_results_path,
  pattern = ".*model.*full\\.csv$",
  full.names = TRUE
)

model_n_by_region <- map(model_files, function(path) {
  read_csv(path, show_col_types = FALSE) %>%
    rename(brain_region = ...1) %>%
    left_join(
      tibble(
        brain_region = names(short_to_long_brain_region),
        brain_region_long = short_to_long_brain_region
      ),
      by = "brain_region"
    ) %>%
    mutate(
      n_samples = as.integer(str_extract(`samples x cell types`, "\\d+")),
      brain_region_model_n = as.character(glue(
        "{brain_region}\n(n = {n_samples})"
      )),
      brain_region_long_model_n = as.character(glue(
        "{brain_region_long}\n(n = {n_samples})"
      ))
    ) %>%
    select(brain_region, brain_region_model_n, brain_region_long_model_n)
}) %>%
  list_rbind() %>%
  filter(brain_region != "ALL")

read_and_annotate_effects <- function(file_path) {
  file_name <- basename(file_path)
  match <- str_match(file_name, "effects-([^-]+)-([^.]+)\\.csv")
  brain_region <- match[2]
  model <- match[3]

  read_csv(file_path, show_col_types = FALSE) %>%
    rename_with(~ str_to_lower(str_replace_all(.x, "\\s+", "_"))) %>%
    mutate(
      brain_region = brain_region,
      model = model,
      is_father_or_children_credible = effect != 0
    )
}

all_coda_effects <- map_dfr(effect_files, read_and_annotate_effects) %>%
  rename(cell_subtype_short = cell_type, logFC = `log2-fold_change`)

full_disease_coda_effects <- all_coda_effects %>%
  filter(
    model == "full",
    covariate == "C(case_control, Treatment('Control'))[T.Case]"
  ) %>%
  mutate(
    logFC_color = case_when(
      !is_father_or_children_credible ~ "non_credible",
      logFC > 0 ~ "up",
      logFC < 0 ~ "down"
    ),
    effect_color = ifelse(effect > 0, "up", "down")
  )


# %% Supplementary figure 2, panel f: faceted cell proportion boxplot
boxplot_faceted <- ggplot(
  cell_counts_per_sample,
  aes(y = brain_region, x = prop_ct_in_sample)
) +
  geom_boxplot(
    aes(fill = case_control),
    outliers = FALSE,
    position = position_dodge(width = 0.8),
    show.legend = FALSE
  ) +
  geom_jitter(
    mapping = aes(group = case_control),
    size = 0.1,
    alpha = .5,
    shape = 16,
    color = "black",
    stroke = 0.1,
    position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8)
  ) +
  facet_wrap(
    ~ factor(cell_subtype_short, levels = correct_cell_subtype_order),
    scales = "free",
    ncol = 2
  ) +
  scale_x_continuous(
    breaks = scales::pretty_breaks(4),
    labels = scales::percent_format()
  ) +
  scale_fill_manual(
    name = "Disease status:",
    values = case_control_colors,
    labels = c("Case" = "Parkinson's Disease")
  ) +
  xlab("Cell subtype proportion") +
  theme_publication(5, font) +
  no_legend +
  theme(axis.title.y = element_blank())

ggsave_pdf_svg(
  base_name = "05-results/02-plots/figure-2-coda-boxplot-cell-subtype-faceted",
  plot = boxplot_faceted,
  width = 80,
  height = 220,
  units = "mm"
)


# %% Supplementary figure 2, panels a-e: jitter + logFC bars per region
to_plot <- cell_counts_per_sample %>%
  left_join(
    full_disease_coda_effects,
    by = c("cell_subtype_short", "brain_region")
  ) %>%
  mutate(
    cell_subtype_short = factor(
      cell_subtype_short,
      levels = rev(correct_cell_subtype_order)
    )
  )

plot_coda_jitter_counts <- function(df) {
  ggplot(
    df,
    aes(y = cell_subtype_short, x = cells_count, color = case_control)
  ) +
    scale_x_continuous(
      transform = "log1p",
      breaks = c(0, 1, 10, 100, 1000),
      limits = c(0, 10000),
      labels = scales::label_comma()
    ) +
    guides(
      x = guide_axis_logticks(negative.small = 1),
      color = guide_legend(
        title = "Disease status:",
        override.aes = list(size = 1, alpha = 1)
      )
    ) +
    geom_boxplot(
      mapping = aes(group = interaction(cell_subtype_short, case_control)),
      color = "#333333",
      outliers = FALSE
    ) +
    geom_jitter(
      position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.8),
      size = 0.3,
      alpha = .4,
      shape = 16
    ) +
    scale_color_manual(
      values = case_control_colors,
      labels = c("Case" = "Parkinson's Disease")
    ) +
    theme_publication(5, font) +
    theme(
      legend.position = "bottom",
      legend.direction = "vertical",
      panel.grid.major.y = element_line(linewidth = 0.1, color = "grey90"),
      axis.title.y = element_blank()
    ) +
    labs(x = "Cell counts\n(log scale)")
}

plot_coda_bar <- function(df, metric, metric_color, limit) {
  ggplot(
    df %>% distinct(cell_subtype_short, .keep_all = TRUE),
    aes(
      y = cell_subtype_short,
      x = .data[[metric]],
      fill = .data[[metric_color]]
    )
  ) +
    geom_col(width = 0.5) +
    scale_fill_manual(
      values = c(
        "up" = "#b31b2c",
        "down" = "#2166ac",
        "non_credible" = "grey15"
      ),
      labels = c(
        "up" = "Positive effect",
        "down" = "Negative effect",
        "non_credible" = "Non-credible effect"
      ),
      name = element_blank()
    ) +
    scale_x_continuous(
      breaks = scales::breaks_extended(n = 5),
      limits = c(-limit, limit)
    ) +
    geom_vline(xintercept = 0, linetype = "solid", linewidth = 0.2) +
    theme_publication(5, font) +
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
      panel.grid.major.x = element_line(color = "grey80", linewidth = 0.1),
      panel.grid.minor.x = element_line(color = "grey90", linewidth = 0.1),
      legend.position = "bottom",
      legend.direction = "vertical",
      plot.margin = margin(l = 2, unit = "mm")
    )
}

plot_coda_combined <- function(df, r, metric, metric_color, limit) {
  jitter_plot <- plot_coda_jitter_counts(df)
  bar_plot <- plot_coda_bar(df, metric, metric_color, limit)

  if (metric == "logFC") {
    bar_plot <- bar_plot +
      labs(x = "Estimated<br>log<sub>2</sub>fold-change") +
      theme(axis.title.x = element_markdown())
  } else {
    bar_plot <- bar_plot + labs(x = "Estimated effect size")
  }

  jitter_plot +
    no_legend +
    bar_plot +
    no_legend +
    plot_layout(widths = c(1, 0.4)) +
    plot_annotation(
      title = str_replace_all(r, "\\n", " "),
      theme = theme(
        plot.title = element_text(
          size = 5,
          family = font,
          hjust = 0.5,
          face = "plain",
          margin = margin(b = 0)
        )
      )
    )
}

regions <- unique(amp_metadata$brain_region_long)
for (r in regions) {
  subset <- to_plot %>% filter(brain_region_long == r)
  r_compatible <- str_replace_all(r, "\\s", "_")

  combined_logFC <- plot_coda_combined(
    subset,
    r,
    "logFC",
    "logFC_color",
    max(abs(to_plot$logFC))
  )
  ggsave_pdf_svg(
    base_name = glue(
      "05-results/02-plots/figure-2-coda-jitter-logFC-{r_compatible}"
    ),
    plot = combined_logFC,
    width = 50,
    height = 75,
    units = "mm"
  )
}


# %% Supplementary figure 2b: credibility heatmaps per region
clean_credibility_file <- function(path) {
  read_csv(path, show_col_types = FALSE) %>%
    mutate(
      Covariate = recode(
        Covariate,
        `C(case_control, Treatment('Control'))[T.Case]_node` = "Disease status [PD]",
        `sex[T.Male]_node` = "Sex [male]",
        `age_at_baseline_mm_scaled_node` = "Age at death"
      ),
      Node = str_remove(Node, "^cell_(supertype|type)_short_")
    ) %>%
    mutate(
      Covariate = factor(
        Covariate,
        levels = c("Disease status [PD]", "Sex [male]", "Age at death")
      )
    ) %>%
    filter(
      Covariate %in% c("Disease status [PD]", "Sex [male]", "Age at death")
    ) %>%
    pivot_longer(
      cols = starts_with("reference"),
      names_to = "reference_cell_type",
      values_to = "final_parameter"
    ) %>%
    mutate(
      reference_cell_type = str_remove(reference_cell_type, "^reference_"),
      reference_cell_type = factor(
        reference_cell_type,
        levels = correct_cell_subtype_order
      ),
      Node = factor(Node, levels = rev(correct_node_order))
    )
}

get_heatmap_plot <- function(heatmap_df) {
  diag_tiles <- heatmap_df %>%
    filter(as.character(Node) == as.character(reference_cell_type))

  reference_parents <- cell_annotations %>%
    transmute(
      reference_cell_type = cell_subtype_short,
      parent_1 = cell_type_short,
      parent_2 = cell_supertype_short
    ) %>%
    pivot_longer(
      cols = starts_with("parent"),
      names_to = NULL,
      values_to = "Node"
    ) %>%
    filter(reference_cell_type != Node)

  parent_diag_tiles <- heatmap_df %>%
    inner_join(
      reference_parents,
      by = c("reference_cell_type", "Node"),
      relationship = "many-to-many"
    )

  heatmap_df %>%
    ggplot(aes(x = reference_cell_type, y = Node, fill = final_parameter)) +
    geom_tile(color = "gray80", linewidth = 0.01) +
    facet_wrap(~Covariate, strip.position = "left", ncol = 1) +
    fill_distiller_gradient(name = "Effect size") +
    geom_tile(data = diag_tiles, fill = "black") +
    geom_tile(data = parent_diag_tiles, fill = "black") +
    theme_publication(4, font) +
    horizontal_scientific_colormap +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
      axis.title.y = element_blank(),
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text = element_text(angle = 90, hjust = 0.5, face = "plain"),
      panel.spacing = unit(0.2, "lines")
    ) +
    labs(x = "Reference")
}

type_subtype_map <- cell_annotations %>%
  select(parent = cell_type_short, child = cell_subtype_short)
supertype_subtype_map <- cell_annotations %>%
  select(parent = cell_supertype_short, child = cell_subtype_short)
supertype_type_map <- cell_annotations %>%
  select(parent = cell_supertype_short, child = cell_type_short)

subtype_parents <- bind_rows(type_subtype_map, supertype_subtype_map) %>%
  filter(parent != child) %>%
  distinct()
all_parents <- bind_rows(
  type_subtype_map,
  supertype_subtype_map,
  supertype_type_map
) %>%
  filter(parent != child) %>%
  distinct()

get_family <- function(node, hierarchy_df) {
  ancestors <- c()
  stack <- c(node)
  visited <- c()
  while (length(stack) > 0) {
    current_node <- stack[1]
    stack <- stack[-1]
    if (current_node %in% visited) {
      next
    }
    visited <- c(visited, current_node)
    ancestors <- c(ancestors, current_node)
    parents <- hierarchy_df %>% filter(child == current_node) %>% pull(parent)
    for (p in parents) {
      if (!(p %in% visited)) stack <- c(stack, p)
    }
  }
  unique(ancestors)
}

get_bar_df <- function(heatmap_df, hierarchy_map) {
  base_stats <- heatmap_df %>%
    group_by(Covariate, Node) %>%
    summarise(
      n_tests = length(unique(reference_cell_type)),
      times_credible = sum(final_parameter != 0),
      .groups = "drop"
    )

  specific_hierarchy_map <- hierarchy_map %>%
    filter(child %in% unique(heatmap_df$Node))
  node_exclusions <- specific_hierarchy_map %>%
    group_by(parent) %>%
    summarise(n_children = n()) %>%
    rename(Node = parent)

  bar_df <- base_stats %>%
    left_join(node_exclusions, by = "Node") %>%
    mutate(
      n_children = ifelse(is.na(n_children), 1, n_children),
      n_possible_tests = n_tests - n_children,
      pct_credible = times_credible / n_possible_tests,
      is_credible = pct_credible > 0.5
    )

  bar_df_propagated <- map(unique(bar_df$Node), function(current_node) {
    bar_df %>%
      group_by(Covariate) %>%
      filter(Node %in% get_family(current_node, all_parents)) %>%
      mutate(propagated_pct_credible = max(pct_credible)) %>%
      filter(Node == current_node) %>%
      select(Covariate, Node, propagated_pct_credible)
  }) %>%
    list_rbind()

  bar_df %>%
    left_join(bar_df_propagated, by = c("Covariate", "Node")) %>%
    mutate(
      Covariate = factor(
        Covariate,
        levels = c("Disease status [PD]", "Sex [male]", "Age at death")
      ),
      Node = factor(Node, levels = rev(correct_node_order))
    )
}

get_bar_plot <- function(bar_df) {
  bar_df %>%
    ggplot(aes(x = pct_credible, y = Node, fill = is_credible)) +
    geom_col() +
    facet_wrap(~Covariate, strip.position = "left", ncol = 1) +
    scale_x_continuous(limits = c(0, 1), labels = scales::label_percent()) +
    scale_fill_manual(values = c("FALSE" = cool_grey, "TRUE" = cool_green)) +
    geom_vline(
      xintercept = .5,
      colour = "black",
      alpha = .5,
      linetype = "dashed",
      linewidth = 0.2
    ) +
    theme_publication(4, font) +
    xlab("Consistency") +
    theme(
      panel.grid.major.x = element_line(linewidth = 0.10, color = "grey85"),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
      legend.position = "None",
      strip.background = element_blank(),
      strip.placement = "None",
      strip.text = element_blank(),
      panel.spacing = unit(0.2, "lines")
    )
}

for (r in credibility_files) {
  short_region_name <- str_extract(
    r,
    ".*credibility-(.*)-full\\.csv$",
    group = 1
  )
  region_name <- short_to_long_brain_region[short_region_name]
  region_compatible <- str_replace_all(region_name, "\\s", "_")
  region <- clean_credibility_file(r)

  credibility_heatmap <- get_heatmap_plot(region) +
    get_bar_plot(get_bar_df(region, subtype_parents)) +
    plot_layout(widths = c(4, 1)) +
    plot_annotation(
      title = str_replace_all(region_name, "\n", " "),
      theme = theme(
        plot.title = element_text(
          size = 5,
          hjust = .5,
          face = "plain",
          margin = margin(b = 2, 0, 0, 0)
        )
      )
    )
  ggsave_pdf_svg(
    base_name = glue(
      "05-results/02-plots/figure-2-coda-multicovariate-credibility-heatmap-{region_compatible}"
    ),
    plot = credibility_heatmap,
    width = 60,
    height = 120,
    units = "mm"
  )
}
