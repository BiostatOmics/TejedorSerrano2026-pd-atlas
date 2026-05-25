# Description: Generate Supplementary Figure 5 panels:
#   panel a — NES scatter comparing PD response in Oligo 1 vs Oligo 2 (GPI)
#   panel b — PD response GSEA dotplot for Oligo 1 and Oligo 2 (GPI)

# %% Setup
library(tidyverse)
library(glue)
library(ggrastr)
library(patchwork)
library(msigdbr)
library(clusterProfiler)
source("utils/my_utils.R")
source("utils/figure_helpers.R")


# %% Load MSigDB
MSIG <- msigdbr(species = "Homo sapiens")
MSIG_OF_INTEREST <- MSIG %>%
  filter(
    !gs_collection %>% str_to_lower() %in% c("c1", "c4", "c6", "c8"),
    !gs_collection_name %>% str_to_lower() %>% str_detect("_legacy"),
    !gs_subcollection %>% str_to_lower() %>% str_detect("vax"),
  )
MSIG_TERM2GENE <- MSIG_OF_INTEREST %>% select(gs_id, ensembl_gene)
MSIG_TERM2NAME <- MSIG_OF_INTEREST %>% select(gs_id, gs_name)
MSIG_METADATA <- MSIG_OF_INTEREST %>%
  distinct(
    gs_id,
    gs_collection,
    gs_subcollection,
    gs_collection_name,
    gs_description
  )


# %% Load pre-computed GSEA results
gsea_pd_effect_oligo2_in_gpi <- read_rds(
  "04-analysis/public/04-gseas-of-interest/gsea-pd-effect-oligo2-in-gpi.rds"
)

gsea_pd_effect_oligo1_in_gpi <- read_rds(
  "04-analysis/public/04-gseas-of-interest/gsea-pd-effect-oligo1-in-gpi.rds"
)


# Supplementary figure 5 panel b: PD GSEA dotplot Oligo 1 & Oligo 2
oxphos_resp <- "ELECTRON_TRANSPORT_CHAIN_OXPHOS|RESPIRATORY_CHAIN_COMPLEX|GOBP_OXIDATIVE_PHOSPHORYLATION|ATP_SYNTHESIS|AEROBIC_RESPIRATION_AND_RESP"
folding <- "ATP_DEPENDENT_PROTEIN|UNFOLDED_PROTEIN_BINDING"
transl <- "GOCC_RIBOSOME|EUKARYOTIC_TRANSLATION_(ELONGATION|INITIATION)|REACTOME_TRANSLATION"
stress <- "REGULATION_OF_INTRINSIC_APOPTOTIC|APOPTOSIS_MODULATION|HYPOXIA|UV_RESPONSE|ARYLSULF"
neurodeg <- "PARKINSON"
cilium <- "AXONEME|CILIARY_PLASM|MICROTUBULE_BUNDLE|CILIPATHIES|CILIUM_ORGANIZATION"
gss <- paste(oxphos_resp, folding, transl, stress, neurodeg, cilium, sep = "|")

top_n_each <- 150

df <- bind_rows(
  gsea_pd_effect_oligo1_in_gpi %>% mutate(cell_subtype_short = "Oligo 1"),
  gsea_pd_effect_oligo2_in_gpi %>% mutate(cell_subtype_short = "Oligo 2")
) %>%
  filter(
    !gs_subcollection %in%
      c(
        "IMMUNESIGDB",
        "TFT:GTRD",
        "MIR:MIRDB",
        "CGP",
        "CP:BIOCARTA",
        "HPO",
        "CP:PID"
      ),
    p.adjust < 0.01
  ) %>%
  mutate(
    sign = if_else(NES > 0, "up", "down"),
    count = str_count(core_enrichment, "/") + 1,
    gene_ratio = count / gs_size,
    gs_name_pretty = prettify_genesets_names(
      gs_name,
      width = 65,
      truncate = TRUE
    ),
    gs_subcollection = if_else(
      gs_collection_name == "Hallmark",
      "HALLMARK",
      gs_subcollection
    ),
    gs_name_pretty = str_glue("{gs_name_pretty} [{gs_subcollection}]"),
  )

df_plot <- df %>%
  group_by(sign, gs_name) %>%
  filter(n() > 1) %>%
  summarise(best_p = min(p.adjust), .groups = "drop") %>%
  slice_min(best_p, n = top_n_each, by = sign) %>%
  select(sign, gs_name) %>%
  left_join(df, by = c("sign", "gs_name")) %>%
  filter(str_detect(gs_name, gss)) %>%
  left_join(get_genesets_clusters(., n_clusters = 7), by = "gs_name") %>%
  group_by(gs_cluster) %>%
  mutate(mean_gene_ratio = mean(gene_ratio, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(gs_name = fct_reorder(gs_name, mean_gene_ratio, .desc = TRUE)) %>%
  mutate(
    gs_cluster = forcats::fct_collapse(
      forcats::fct_relevel(gs_cluster, "4", "5", "1", "6", "3", "2", "7"),
      "Translation" = "4",
      "Folding by\nchaperones" = "5",
      "Aerobic\nrespiration" = "1",
      "Apoptosis &\nstress-related" = c("3"),
      "Cilium" = "2",
      " " = "6",
      "  " = "7",
    )
  )

scales <- list(
  x = scale_x_continuous(limits = range(df_plot$gene_ratio)),
  color = scale_color_viridis_c(
    name = expression(-log[10]("FDR")),
    limits = range(-log10(df_plot$p.adjust)),
    option = "viridis",
    direction = 1
  ),
  size = scale_size_continuous(
    name = "Gene\nset size",
    limits = range(df_plot$gs_size),
    range = c(.5, 2),
  ),
  shape = scale_shape_manual(
    name = "Cell type",
    values = c("Oligo 1" = 1, "Oligo 2" = 13),
    labels = c("Oligo 2" = "Oligo 2\n(vuln.)"),
    guide = guide_legend(override.aes = list(size = 1))
  )
)

make_plot <- function(
  sign_label,
  title = sign_label,
  show_legend = TRUE,
  hide_x = FALSE
) {
  ggplot(
    df_plot %>% filter(sign == sign_label),
    aes(gene_ratio, reorder(gs_name_pretty, gene_ratio))
  ) +
    geom_point(
      aes(size = gs_size, color = -log10(p.adjust), shape = cell_subtype_short),
      alpha = 1,
      stroke = .3
    ) +
    facet_grid(gs_cluster ~ ., scales = "free_y", space = "free_y") +
    scales +
    labs(title = title, x = if (hide_x) "" else "Gene ratio", y = NULL) +
    theme_publication(6) +
    coord_cartesian(clip = "off") +
    vertical_scientific_colormap +
    theme(
      axis.text.x = if (hide_x) element_blank() else element_text(),
      axis.line.x = if (hide_x) element_blank() else element_line(),
      axis.ticks.x = if (hide_x) element_blank() else element_line(),
      legend.position = if (show_legend) "right" else "none",
      panel.grid.major.y = element_line(color = "grey85", linewidth = 0.1),
      strip.clip = "off",
      plot.title = ggtext::element_textbox(
        hjust = 0.5,
        halign = 0.5,
        box.color = "black",
        linetype = 1,
        linewidth = 0.2,
        color = "black",
        padding = margin(2, 2, 2, 2),
        margin = margin(b = 0, t = 1),
        width = unit(1, "npc")
      ),
      plot.title.position = "panel",
    )
}

gseaplot <- make_plot(
  "up",
  title = "Upregulated in PD",
  show_legend = TRUE,
  hide_x = TRUE
) /
  make_plot(
    "down",
    title = "Downregulated in PD",
    show_legend = FALSE,
    hide_x = FALSE
  ) +
  plot_layout(heights = c(1, .3))

ggsave_pdf_svg(
  base_name = "05-results/02-plots/figure-5-pd-oligo1-and-oligo2-gpi-gseaplot",
  plot = gseaplot,
  width = 180,
  height = 90,
  units = "mm",
)


# Supplementary figure 5 panel a: NES scatter Oligo 1 vs Oligo 2
df_scatter <- gsea_pd_effect_oligo2_in_gpi %>%
  select(gs_name, gs_subcollection, NES, p.adjust) %>%
  rename(NES_oligo2 = NES, padj_oligo2 = p.adjust) %>%
  inner_join(
    gsea_pd_effect_oligo1_in_gpi %>%
      select(gs_name, NES, p.adjust) %>%
      rename(NES_oligo1 = NES, padj_oligo1 = p.adjust),
    by = "gs_name"
  ) %>%
  filter(
    !gs_subcollection %in%
      c(
        "IMMUNESIGDB",
        "TFT:GTRD",
        "MIR:MIRDB",
        "CGP",
        "CP:BIOCARTA",
        "HPO",
        "CP:PID"
      ),
  ) %>%
  mutate(
    sig = if_else(padj_oligo1 < 0.01 | padj_oligo2 < 0.01, "sign", "nosign")
  )

rho <- cor(df_scatter$NES_oligo1, df_scatter$NES_oligo2, method = "spearman")

# Stats related to supp-figure 5
# Spearman correlation between Oligo 1 and Oligo 2 NES values (all gene sets)
cor.test(df_scatter$NES_oligo1, df_scatter$NES_oligo2, method = "spearman")

p_scatter <- ggplot(df_scatter, aes(NES_oligo1, NES_oligo2, color = sig)) +
  rasterise(geom_point(alpha = 0.5, stroke = 0), dpi = dpi) +
  geom_abline(slope = 1, intercept = 0, linewidth = .1, linetype = "dashed") +
  geom_vline(xintercept = 0, linewidth = .1) +
  geom_hline(yintercept = 0, linewidth = .1) +
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    label = paste0("italic(rho)==", round(rho, 2)),
    parse = TRUE,
    hjust = 1.1,
    vjust = -0.5
  ) +
  scale_color_manual(
    values = c("nosign" = "grey80", "sign" = "#d62728"),
    labels = c("sign" = "FDR < 0.01", "nosign" = "FDR > 0.01"),
    guide = guide_legend(override.aes = list(alpha = 1, size = 1))
  ) +
  theme_publication(6) +
  coord_fixed() +
  labs(x = "NES Oligo 1", y = "NES Oligo 2", color = NULL) +
  theme(legend.position = "bottom")

ggsave_pdf_svg(
  base_name = "05-results/02-plots/figure-5-pd-oligo1-and-oligo2-gpi-scatterplot",
  plot = p_scatter,
  width = 45,
  height = 45,
  units = "mm",
)
