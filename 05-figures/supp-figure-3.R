# Description: Generate Supplementary Figure 3 panels:
#   panels a-e — multi-covariate effect dotplot (disease status, sex, age)

# %% Setup
library(tidyverse)
library(glue)
source("utils/my_utils.R")


# %% Load data
coda_results_path <- "04-analysis/02-differential-abundance-phi5-itermore-peri-final"

effect_files <- list.files(
  coda_results_path,
  pattern = ".*effects.*\\.csv$",
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
      ))
    ) %>%
    select(brain_region, brain_region_model_n)
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


# %% Helper functions for credibility consistency
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
    if (current_node %in% visited) next
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

node_consistency <- credibility_files %>%
  map(function(path) {
    region_name <- str_match(path, ".*credibility-(.*)-full\\.csv$")[, 2]
    clean_credibility_file(path) %>%
      get_bar_df(subtype_parents) %>%
      mutate(brain_region = region_name)
  }) %>%
  list_rbind() %>%
  filter(brain_region != "ALL") %>%
  rename(
    covariate = Covariate,
    cell_subtype_short = Node,
    consistency = propagated_pct_credible
  )


# %% Supplementary figure 3, panels a-e: multi-covariate effect dotplot
multicovariate_heatmap <- all_coda_effects %>%
  mutate(
    covariate = recode(
      covariate,
      `C(case_control, Treatment('Control'))[T.Case]` = "Disease status [PD]",
      `sex[T.Male]` = "Sex [male]",
      `age_at_baseline_mm_scaled` = "Age at death"
    )
  ) %>%
  mutate(
    covariate = factor(
      covariate,
      levels = rev(c("Disease status [PD]", "Sex [male]", "Age at death"))
    ),
    fill_category = case_when(
      median == 0 & sd == 0 ~ "reference",
      effect == 0 ~ "non_credible",
      TRUE ~ "gradient"
    ),
    cell_subtype_short = factor(
      cell_subtype_short,
      levels = correct_cell_subtype_order
    )
  ) %>%
  filter(
    model == "full",
    brain_region != "ALL",
    !str_detect(covariate, "cohort")
  ) %>%
  left_join(
    node_consistency,
    by = c("covariate", "cell_subtype_short", "brain_region")
  )

multicovariate_plot <- ggplot() +
  facet_wrap(
    ~ factor(
      brain_region,
      levels = rev(c("DMNX", "GPI", "PFC", "PMC", "PVC"))
    ),
    nrow = 5,
    strip.position = "left",
    scales = "free_x",
    labeller = as_labeller(setNames(
      model_n_by_region$brain_region_model_n,
      model_n_by_region$brain_region
    ))
  ) +
  geom_point(
    data = multicovariate_heatmap,
    aes(
      x = cell_subtype_short,
      y = covariate,
      color = effect,
      size = consistency
    )
  ) +
  scale_radius(
    name = "Consistency",
    range = c(-.25, 2),
    labels = scales::label_percent()
  ) +
  color_distiller_gradient(name = "Propagated\neffect size:") +
  theme_publication(5, font) +
  vertical_scientific_colormap +
  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(angle = 90, hjust = 0.5, face = "plain"),
    panel.spacing = unit(0, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    legend.position = "right",
    legend.direction = "vertical",
    strip.clip = "off"
  )

ggsave_pdf_svg(
  base_name = "05-results/02-plots/figure-2-coda-multicovariate-heatmap",
  plot = multicovariate_plot,
  width = 80,
  height = 80,
  units = "mm"
)
