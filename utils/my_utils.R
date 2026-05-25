# Description: Shared R utilities for publication figure generation.
#   Includes logging helpers, ggplot2 theme, color palettes, ordering vectors,
#   custom scales, gene set name prettifier, and data loaders.

library(cli)
library(glue)
library(ggplot2)
library(readr)
library(dplyr)


# %% Logging helpers

log_info <- function(msg) {
  cli::cli_alert_info(paste0(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    " - ",
    msg
  ))
}

log_success <- function(msg) {
  cli::cli_alert_success(paste0(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    " - ",
    msg
  ))
}

log_warning <- function(msg) {
  cli::cli_alert_warning(paste0(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    " - ",
    msg
  ))
}

log_danger <- function(msg) {
  cli::cli_alert_danger(paste0(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    " - ",
    msg
  ))
}


# %% Save helper: write both PDF and SVG from a single call
ggsave_pdf_svg <- function(base_name, device = "pdf", ...) {
  ggplot2::ggsave(
    filename = glue::glue("{base_name}.pdf"),
    device = device,
    ...
  )
  ggplot2::ggsave(filename = glue::glue("{base_name}.svg"), ...)
}


# %% Publication theme and colormap helpers

font <- "Helvetica" # Nature journals font
dpi <- 1000
fig_width <- 90   # one-column width in mm
fig_height <- 70
scale_factor <- 188 / fig_width  # used by figure_helpers.R UMAP annotation helpers

theme_publication <- function(base_size = 8, base_family = font) {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(
        family = base_family,
        size = base_size,
        color = "black"
      ),
      axis.title = element_text(size = base_size, color = "black"),
      axis.text = element_text(size = base_size, color = "black"),
      axis.line = element_line(linewidth = 0.1),
      axis.ticks = element_line(linewidth = 0.1),
      plot.title = element_text(size = base_size, face = "plain", hjust = 0.5),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size),
      legend.key.size = unit(2, "mm"),
      legend.position = "right",
      legend.background = element_blank(),
      legend.key = element_blank(),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "plain"),
      plot.margin = margin(0, 0, 0, 0)
    )
}

no_legend <- theme(legend.position = "none")

# Colormap legend styling: adds ticks and frame to continuous color scales
base_scientific_colormap <- theme(
  legend.ticks.length = unit(c(-.1, 0), "lines"),
  legend.ticks = element_line(color = "black", linewidth = 0.1),
  legend.frame = element_rect(fill = NA, color = "black", linewidth = 0.1),
)

horizontal_scientific_colormap <- base_scientific_colormap +
  theme(
    legend.key.width = unit(1, "lines"),
    legend.key.height = unit(0.3, "lines"),
    legend.box = "horizontal",
    legend.justification = "center",
    legend.title.position = "top",
    legend.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    legend.direction = "horizontal",
  )

vertical_scientific_colormap <- base_scientific_colormap +
  theme(
    legend.key.width = unit(0.3, "lines"),
    legend.key.height = unit(0.5, "lines"),
    legend.box = "vertical",
    legend.justification = "center",
    legend.title.position = "top",
    legend.title = element_text(hjust = 0.5),
    legend.position = "right",
    legend.direction = "vertical",
  )


# %% Cell type annotations and color/order vectors
# Loaded at source() time; scripts must be run from the project root.
cell_annotations <- readr::read_tsv(
  "03-data/annotations_7_11_agg.txt"
) %>%
  dplyr::select(!dplyr::starts_with("..."))

correct_cell_subtype_colors <- setNames(
  cell_annotations$cell_subtype_color_vibrant,
  cell_annotations$cell_subtype_short
)

correct_cell_type_colors <- setNames(
  unique(cell_annotations$cell_type_color_vibrant),
  unique(cell_annotations$cell_type_short)
)

correct_cell_supertype_colors <- setNames(
  unique(cell_annotations$cell_supertype_color),
  unique(cell_annotations$cell_supertype_short)
)

correct_cell_subtype_order <- c(
  "Astro 1",
  "Astro 2",
  "Astro 3",
  "OPC",
  "Oligo 1",
  "Oligo 2",
  "Epen",
  "Micro 1",
  "Micro 2",
  "CTL",
  "Endo",
  "Peri",
  "Exc L2/3 IT",
  "Exc L4/5 IT",
  "Exc L5/6 NP",
  "Exc L6 IT CA3",
  "Exc L6/6b CT",
  "Inh LAMP5",
  "Inh VIP",
  "Inh PVALB",
  "Inh SST",
  "DMNX/GPI"
)

correct_cell_type_order <- c(
  "Astro",
  "OPC",
  "Oligo",
  "Epen",
  "Micro",
  "CTL",
  "Endo",
  "Peri",
  "Exc",
  "Inh",
  "DMNX/GPI"
)

correct_node_order <- c(
  # subtypes
  "Astro 1",
  "Astro 2",
  "Astro 3",
  "OPC",
  "Oligo 1",
  "Oligo 2",
  "Epen",
  "Micro 1",
  "Micro 2",
  "CTL",
  "Endo",
  "Peri",
  "Exc L2/3 IT",
  "Exc L4/5 IT",
  "Exc L5/6 NP",
  "Exc L6 IT CA3",
  "Exc L6/6b CT",
  "Inh LAMP5",
  "Inh VIP",
  "Inh PVALB",
  "Inh SST",
  "DMNX/GPI",
  # types
  "Astro",
  "Oligo",
  "Micro",
  "Exc",
  "Inh",
  # supertypes
  "NonNeu",
  "Immu",
  "Neu",
  "Vasc"
)

correct_brain_region_order <- c("DMNX", "GPI", "PMC", "PFC", "PVC")

correct_case_control_order <- c(
  "Parkinson's\nDisease",
  "Control",
  "Other disorder\n(Parkinsonism or\nLB Dementia)"
)

short_to_long_brain_region <- c(
  "DMNX" = "Dorsal Motor Nucleus\nof the Xth Nerve",
  "GPI" = "Globus Pallidus Interna",
  "PMC" = "Primary Motor Cortex",
  "PFC" = "Dorsolateral Prefrontal Cortex",
  "PVC" = "Primary Visual Cortex"
)

short_to_long_case_control <- c(
  "Case" = "Parkinson's\nDisease",
  "Control" = "Control",
  "Other" = "Other disorder\n(Parkinsonism or\nLB Dementia)"
)


# %% Color palettes

case_control_colors <- c("Case" = "#A5243D", "Control" = "#90abb0")

region_colors <- c(
  "DMNX" = "#dbb481",
  "GPI" = "#bb5c55",
  "PFC" = "#977faf",
  "PMC" = "#d06eac",
  "PVC" = "#749b86"
)

region_shapes <- c(DMNX = 19, GPI = 17, PFC = 15, PMC = 18, PVC = 3)

brain_bank_colors <- c(
  "UM" = "#da9021",
  "UD" = "#219fda",
  "MS" = "#7259cc",
  "HA" = "#7acc59"
)

sex_colors <- c("Male" = "#56B4E9", "Female" = "#D55E00")
binary_colors <- c("#56B4E9", "#D55E00")

cool_grey <- "#D2D4C8"
cool_blue <- "#2166AC"
cool_red <- "#b2182b"
cool_green <- "#248232"


# %% Custom divergent fill/color scales (RdBu, midpoint anchored at 0)

fill_distiller_gradient <- function(...) {
  scale_fill_gradientn(
    colours = scales::pal_brewer(palette = "RdBu", direction = -1)(7),
    rescaler = function(x, to = c(0, 1), from = NULL, mid = 0) {
      scales::rescale_mid(x, to = to, from = from, mid = mid)
    },
    ...
  )
}

color_distiller_gradient <- function(...) {
  scale_color_gradientn(
    colours = scales::pal_brewer(palette = "RdBu", direction = -1)(7),
    rescaler = ~ scales::rescale_mid(.x, mid = 0),
    ...
  )
}

neg_log10_transform <- scales::new_transform(
  name = "-log10",
  transform = function(x) -log10(x),
  inverse = function(x) 10^-x
)


# %% Gene set name prettifier

prettify_genesets_names <- function(x, width = 30, truncate = FALSE) {
  x <- x %>%
    # KEGG MEDICUS: move type qualifier to end in parentheses
    stringr::str_replace(
      "^KEGG_MEDICUS_([^_]+)_(.*)",
      function(m) {
        full <- stringr::str_match(m, "^KEGG_MEDICUS_([^_]+)_(.*)")
        tipo <- full[, 2] |> stringr::str_sub(1, 3) |> stringr::str_to_upper()
        nombre <- full[, 3]
        glue::glue("{nombre} ({tipo}.)")
      }
    ) %>%
    # Remove standard database prefixes
    stringr::str_remove(
      "^(HALLMARK|REACTOME|BIOCARTA|PID|WP|KEGG|HP|GO(BP|CC|MF))_"
    ) %>%
    # miRNA strand suffixes
    stringr::str_replace("_([35]P)$", " (\\1)") %>%
    # Replace remaining underscores with spaces
    stringr::str_replace_all("_", " ") %>%
    stringr::str_trim()

  if (truncate) stringr::str_trunc(x, width) else stringr::str_wrap(x, width)
}


# %% Data loaders

load_pseudobulk_dge <- function() {
  readRDS("04-analysis/pseudobulk-dge-preprocessed.rds")
}

load_gene_metadata <- function() {
  readr::read_csv(
    "03-data/gene-metadata.csv",
    col_types = readr::cols(entrezgene_id = "c")
  ) %>%
    dplyr::filter(transcript_is_canonical == 1) %>%
    dplyr::select(
      ensembl_gene_id,
      hgnc_symbol,
      entrezgene_id,
      gene_biotype,
      chromosome_name,
      percentage_gene_gc_content,
      transcript_length
    ) %>%
    dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)
}
