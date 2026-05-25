# Shared helper functions for figure scripts.
# Source this file in scripts that need UMAP plotting helpers or gene set clustering.

# UMAP plotting helpers
# (used in figure-1.R and supp-figure-1.R)

nuclei_number_annotation <- function(n, size = 6) {
  annotate(
    "text",
    x = -Inf,
    y = -Inf,
    label = glue("{format(n, big.mark = ',')} nuclei"),
    hjust = -0.05,
    vjust = -0.8,
    size = size * scale_factor / .pt,
    family = font
  )
}

top_label_annotation <- function(label, x_position, size = 6) {
  annotate(
    "label",
    x = x_position,
    y = Inf,
    label = label,
    size = size * scale_factor / .pt,
    family = font,
    vjust = "inward",
    hjust = "inward"
  )
}

geom_umap_points_mapped_color <- function(data, color = NULL) {
  geom_point(
    data = data,
    aes(x = X_umap, y = Y_umap, color = {{ color }}),
    size = 0.1 / 1000,
    alpha = 0.3,
    shape = 16
  )
}

geom_umap_points_fixed_color <- function(data, color = NULL) {
  geom_point(
    data = data,
    aes(x = X_umap, y = Y_umap),
    color = color,
    size = 0.1 / 1000,
    alpha = 0.3,
    shape = 16
  )
}

cluster_annotation <- function(
  centroids_df,
  mapping,
  bg.color = "white",
  bg.r = 0.1,
  fontface = "bold",
  size = 4.5
) {
  geom_text_repel(
    data = centroids_df,
    mapping = mapping,
    size = size * scale_factor / .pt,
    fontface = fontface,
    box.padding = 0.5,
    force = 1.25,
    segment.color = "black",
    min.segment.length = 1.5,
    bg.color = bg.color,
    bg.r = bg.r,
    lineheight = 0.8
  )
}

theme_umap_plot <- function(size = 8) {
  theme_publication(base_size = size * scale_factor) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      legend.position = "none",
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.25 * scale_factor
      ),
      plot.margin = unit(rep(0.7 * scale_factor, 4), "mm")
    )
}

cluster_centroid <- function(data, annotation) {
  annotation_sym <- sym(annotation)
  total <- nrow(data)
  props <- data %>%
    dplyr::count(!!annotation_sym, name = "count") %>%
    mutate(prop = count / total)
  centroids <- data %>%
    group_by(!!annotation_sym) %>%
    summarise(X = median(X_umap), Y = median(Y_umap), .groups = "drop")
  left_join(props, centroids, by = annotation) %>%
    mutate(
      label = paste(
        !!annotation_sym,
        glue("({percent(prop, accuracy = 0.1)})"),
        sep = "\n"
      )
    ) %>%
    dplyr::rename(color = !!annotation_sym)
}


# Gene set clustering helper
# (used in figure-5.R and supp-figure-5.R)
# Requires MSIG (msigdbr tibble) to be loaded in the calling environment.

get_genesets_clusters <- function(gsea_results_df, n_clusters = 5) {
  target_genesets <- gsea_results_df$gs_name

  genesets_and_its_genes <- MSIG %>%
    filter(gs_name %in% target_genesets) %>%
    distinct(gs_name, gene_symbol)

  freq_table <- table(
    genesets_and_its_genes$gs_name,
    genesets_and_its_genes$gene_symbol
  )

  intersection <- tcrossprod(freq_table)
  sizes <- rowSums(freq_table)
  min_sizes <- outer(sizes, sizes, FUN = pmin)
  overlap_matrix <- intersection / min_sizes

  dist_matrix <- as.dist(1 - overlap_matrix)
  clustering <- hclust(dist_matrix, method = "average")
  clusters <- cutree(clustering, k = n_clusters)

  enframe(
    factor(clusters),
    name = "gs_name",
    value = "gs_cluster"
  )
}
