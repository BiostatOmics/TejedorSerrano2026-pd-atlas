# Description: Generate Figure 2 panels:
#   panels a-e — forest + tree plots per brain region

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
# Path to tascCODA results directory
coda_results_path <- "04-analysis/02-differential-abundance-phi5-itermore-peri-final"

amp_metadata <- read_csv("03-data/00-AMP_V1_metadata.csv")
cell_annotations <- read_tsv("03-data/annotations_7_11_agg.txt")
clinical_metadata <- read_csv(
  "03-data/00-AMP_V1_clinical_metadata.csv"
)

# %% Format metadata
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

# Calculate proportion of each cell subtype in each brain region
subtype_region_props <- amp_metadata %>%
  group_by(cell_subtype_short, brain_region) %>%
  summarise(total_cells_in_subtype_region = n(), .groups = "drop") %>%
  group_by(cell_subtype_short) %>%
  mutate(
    subtype_regional_distribution = total_cells_in_subtype_region /
      sum(total_cells_in_subtype_region)
  ) %>%
  ungroup()

# Cell counts per sample x subtype (only existing combinations)
real_counts <- amp_metadata %>%
  count(sample_id, cell_subtype_short, name = "cells_count")

# Sample-level metadata
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
  # Create all sample x subtype combinations (including those with 0 cells)
  crossing(cell_subtype_short = unique(amp_metadata$cell_subtype_short)) %>%
  left_join(real_counts, by = c("sample_id", "cell_subtype_short")) %>%
  mutate(cells_count = replace_na(cells_count, 0)) %>%
  left_join(
    subtype_region_props,
    by = c("cell_subtype_short", "brain_region")
  ) %>%
  filter(subtype_regional_distribution >= MIN_PROP) %>%
  # Relative abundance within each sample (composition profile)
  group_by(sample_id) %>%
  mutate(prop_ct_in_sample = cells_count / sum(cells_count)) %>%
  ungroup() %>%
  # Donor contribution to each cell type's pool
  group_by(cell_subtype_short, brain_region) %>%
  mutate(donor_contribution_to_ct_pool = cells_count / sum(cells_count)) %>%
  ungroup()

regions_map <- cell_counts_per_sample %>%
  ungroup() %>%
  distinct(brain_region, brain_region_long)

# %% Read CoDA results files
effect_files <- list.files(
  coda_results_path,
  pattern = ".*effects.*\\.csv$",
  full.names = TRUE
)
intercept_files <- list.files(
  coda_results_path,
  pattern = ".*intercept.*\\.csv$",
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

# Helper: read and annotate an effects CSV file
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

# Helper: read and annotate a node CSV file
read_and_annotate_node <- function(file_path) {
  file_name <- basename(file_path)
  match <- str_match(file_name, "node-([^-]+)-([^.]+)\\.csv")
  brain_region <- match[2]
  model <- match[3]

  read_csv(file_path, show_col_types = FALSE) %>%
    rename_with(~ str_to_lower(str_replace_all(.x, "\\s+", "_"))) %>%
    mutate(brain_region = brain_region, model = model)
}

all_coda_effects <- map_dfr(effect_files, read_and_annotate_effects) %>%
  rename(cell_subtype_short = cell_type, logFC = `log2-fold_change`)
all_coda_node <- map_dfr(node_files, read_and_annotate_node)

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

full_disease_coda_node <- all_coda_node %>%
  left_join(regions_map, by = "brain_region") %>%
  filter(
    model == "full",
    covariate == "C(case_control, Treatment('Control'))[T.Case]_node",
    brain_region != "ALL"
  ) %>%
  rename(hdi3 = `hdi_3%`, hdi97 = `hdi_97%`) %>%
  mutate(
    label = str_remove_all(node, "cell.*_"),
    abs_final_parameter = abs(final_parameter)
  ) %>%
  select(-node)

# %% Figure 2, panels a-e: forest + tree plots per region
{
  get_forest_df <- function(coda_effects_df) {
    coda_effects_df %>%
      left_join(regions_map, by = "brain_region") %>%
      rename(hdi3 = `hdi_3%`, hdi97 = `hdi_97%`) %>%
      filter(brain_region != "ALL") %>%
      select(
        cell_subtype_short,
        effect,
        median,
        hdi3,
        hdi97,
        brain_region_long
      ) %>%
      mutate(
        cell_subtype_short = factor(
          cell_subtype_short,
          levels = rev(correct_cell_subtype_order)
        ),
        effect_type = case_when(
          median == 0 ~ "Reference",
          effect == 0 ~ "Non-credible effect",
          effect > 0 ~ "Positive effect",
          effect < 0 ~ "Negative effect"
        )
      )
  }

  get_forest_plot <- function(forest_df) {
    forest_df %>%
      ggplot(aes(x = median, y = cell_subtype_short, color = effect_type)) +
      geom_vline(
        xintercept = 0,
        linetype = "dashed",
        linewidth = 0.2,
        color = "grey50"
      ) +
      geom_point(size = 1.5, shape = 15) +
      geom_linerange(
        aes(xmin = hdi3, xmax = hdi97),
        linewidth = 0.2,
        alpha = 0.5
      ) +
      scale_color_manual(
        values = c(
          "Positive effect" = "#b31b2c",
          "Negative effect" = "#2166ac",
          "Non-credible effect" = "grey15",
          "Reference" = "#29ac20"
        )
      ) +
      scale_x_continuous(breaks = scales::extended_breaks(n = 7)) +
      labs(x = "Propagated effect size", y = "Cell subtype") +
      theme_publication(6, font) +
      theme(
        plot.title = element_text(size = 6, face = "plain"),
        legend.title = element_blank()
      )
  }

  get_tree_dfs <- function(coda_node_df) {
    hierarchy <- cell_annotations %>%
      filter(cell_subtype_short %in% coda_node_df$label) %>%
      select(cell_supertype_short, cell_type_short, cell_subtype_short) %>%
      unique() %>%
      mutate(
        path_string = case_when(
          cell_type_short == cell_subtype_short ~
            paste("all", cell_supertype_short, cell_type_short, sep = "->"),
          TRUE ~
            paste(
              "all",
              cell_supertype_short,
              cell_type_short,
              cell_subtype_short,
              sep = "->"
            )
        )
      )

    data_tree <- data.tree::FromDataFrameTable(
      hierarchy,
      pathName = "path_string",
      pathDelimiter = "->"
    )
    phylo_tree <- data.tree::as.phylo.Node(data_tree)
    phylo_tree$tip.label <- str_replace_all(phylo_tree$tip.label, "_", " ")
    phylo_tree$node.label <- str_replace_all(phylo_tree$node.label, "_", " ")
    phylo_tree$node.label[1] <- ""

    tree <- ape::rotateConstr(phylo_tree, rev(correct_cell_subtype_order))
    tree_metadata <- as_tibble(phylo_tree) %>%
      left_join(coda_node_df, by = "label")

    list(tree = tree, tree_metadata = tree_metadata)
  }

  get_tree_plot <- function(tree, tree_metadata) {
    ggtree(tree, ladderize = FALSE, linewidth = 0.1) %<+%
      tree_metadata +
      geom_nodelab(hjust = 1.2, vjust = -0.7, family = font, size = 1.7) +
      geom_nodepoint(
        data = . %>% filter(final_parameter != 0),
        mapping = aes(size = abs(final_parameter), color = final_parameter)
      ) +
      geom_tippoint(
        data = . %>% filter(final_parameter != 0),
        mapping = aes(size = abs(final_parameter), color = final_parameter)
      ) +
      scale_color_gradientn(
        colours = scales::pal_brewer(palette = "RdBu", direction = -1)(7),
        limits = c(-0.4, 0.4),
        name = "Effect size",
        breaks = c(-0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4),
        guide = guide_legend(
          override.aes = list(
            size = scales::rescale(
              abs(c(-0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4)),
              from = c(0, 0.4),
              to = c(0, 3)
            )
          )
        )
      ) +
      scale_radius(
        limits = c(0, 0.4),
        range = c(0, 3),
        name = "Effect size",
        guide = "none"
      ) +
      theme(
        text = element_text(family = font, size = 6),
        legend.key.size = unit(3, "mm")
      )
  }

  forest_df <- get_forest_df(full_disease_coda_effects)
  for (r in unique(forest_df$brain_region_long)) {
    r_compatible <- str_replace_all(r, "\\s", "_")

    forest_df_subset <- forest_df %>% filter(brain_region_long == r)
    forest_plot <- get_forest_plot(forest_df_subset) +
      labs(title = str_replace_all(r, "\\n", " "))

    node_subset <- full_disease_coda_node %>% filter(brain_region_long == r)
    tree_objects <- get_tree_dfs(node_subset)
    tree_plot <- get_tree_plot(tree_objects$tree, tree_objects$tree_metadata)

    tree_forest_plot <- ((forest_plot +
      no_legend +
      theme(
        axis.title.y = element_blank(),
        axis.text.y = element_text(hjust = 0)
      )) %>%
      insert_left(tree_plot, width = 0.8))

    ggsave_pdf_svg(
      base_name = glue(
        "05-results/02-plots/figure-2-coda-forest&treeplot-{r_compatible}"
      ),
      plot = tree_forest_plot,
      width = 90,
      height = 60,
      units = "mm"
    )
  }
}
