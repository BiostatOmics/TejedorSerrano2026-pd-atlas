# Description: Generate Supplementary Figure 3 panels:
#   panels a-e — credibility heatmaps per brain region

# %% Setup
library(tidyverse)
library(patchwork)
library(glue)
source("utils/my_utils.R")


# %% Load data
coda_results_path <- "04-analysis/02-differential-abundance-phi5-itermore-peri-final"

credibility_files <- list.files(
  coda_results_path,
  pattern = ".*credibility.*\\.csv$",
  full.names = TRUE
)


# %% Helper functions
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


# %% Supplementary figure 3, panels a-e: credibility heatmaps per region
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
