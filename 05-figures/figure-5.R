# Description: Generate Figure 5 panels:
#   figure-5a panel a — biological difference Oligo2 vs Oligo1 volcano (GPI)
#   figure-5a panels b, c — biological difference Oligo2 vs Oligo1 combined GSEA
#   figure-5b panel a — PD differential response Oligo2 vs Oligo1 volcano (GPI)
#   figure-5b panels b, c — PD differential response Oligo2 vs Oligo1 combined GSEA

# %% Setup
library(tidyverse)
library(glue)
library(ggrepel)
library(ggrastr)
library(patchwork)
library(edgeR)
library(limma)
library(msigdbr)
library(clusterProfiler)
library(enrichplot)
source("utils/my_utils.R")
source("utils/figure_helpers.R")


# %% Load data
gene_metadata <- load_gene_metadata()
dge <- load_pseudobulk_dge()

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

go_dbs <- c("GO:BP", "GO:MF", "GO:CC")
go_ids_df <- MSIG_OF_INTEREST %>%
  filter(gs_subcollection %in% go_dbs) %>%
  distinct(gs_id, gs_exact_source)


# %% Pre-processing: filter to Oligo 1 and Oligo 2 samples
oligo_samples <- dge$samples$cell_subtype_short %in%
  c("Oligo 1", "Oligo 2") &
  dge$samples$case_control %in% c("Case", "Control")

dge_oligos <- dge[, oligo_samples]

dge_oligos$samples <- dge_oligos$samples %>%
  mutate(
    rowname = colnames(dge$dge_oligos),
    group = factor(paste(
      cell_subtype_short %>%
        str_replace_all(" ", "_") %>%
        str_replace_all("/", "__"),
      brain_region,
      case_control,
      sep = "."
    ))
  )

keep_genes <- filterByExpr(
  dge_oligos,
  group = dge_oligos$samples$group,
  min.count = 1,
  min.total.count = 50,
  min.prop = 0.85
)
dge_oligos <- dge_oligos[keep_genes, ]

dge_oligos <- calcNormFactors(dge_oligos, method = "TMMwsp")

design <- model.matrix(
  ~ 0 +
    group +
    age_scaled +
    pmi_scaled +
    log_libsize_scaled +
    sex +
    cohort,
  data = dge_oligos$samples
)

colnames(design) <- str_remove_all(colnames(design), "group")


# %% Load pre-computed voomLmFit model
oligo_fit <- readRDS(
  "04-analysis/pseudobulk-dea-voomLmFit-oligos.rds"
)


# %% Contrasts definition and eBayes
contrasts <- makeContrasts(
  PDvsCTRL_Oligo1_GPI = Oligo_1.GPI.Case - Oligo_1.GPI.Control,
  PDvsCTRL_Oligo2_GPI = Oligo_2.GPI.Case - Oligo_2.GPI.Control,
  Oligo2_vs_Oligo1_CTRL_GPI = Oligo_2.GPI.Control - Oligo_1.GPI.Control,
  Oligo2_vs_Oligo1_GPI = (Oligo_2.GPI.Case + Oligo_2.GPI.Control) /
    2 -
    (Oligo_1.GPI.Case + Oligo_1.GPI.Control) / 2,
  Oligo2_vs_Oligo1 = ((Oligo_2.DMNX.Case + Oligo_2.DMNX.Control) +
    (Oligo_2.GPI.Case + Oligo_2.GPI.Control) +
    (Oligo_2.PFC.Case + Oligo_2.PFC.Control) +
    (Oligo_2.PMC.Case + Oligo_2.PMC.Control) +
    (Oligo_2.PVC.Case + Oligo_2.PVC.Control)) /
    10 -
    ((Oligo_1.DMNX.Case + Oligo_1.DMNX.Control) +
      (Oligo_1.GPI.Case + Oligo_1.GPI.Control) +
      (Oligo_1.PFC.Case + Oligo_1.PFC.Control) +
      (Oligo_1.PMC.Case + Oligo_1.PMC.Control) +
      (Oligo_1.PVC.Case + Oligo_1.PVC.Control)) /
      10,
  PD_effect_Oligo2_vs_Oligo1_GPI = (Oligo_2.GPI.Case - Oligo_2.GPI.Control) -
    (Oligo_1.GPI.Case - Oligo_1.GPI.Control),
  PD_effect_GPI_vs_rest_Oligo2 = (Oligo_2.GPI.Case - Oligo_2.GPI.Control) -
    ((Oligo_2.DMNX.Case - Oligo_2.DMNX.Control) +
      (Oligo_2.PFC.Case - Oligo_2.PFC.Control) +
      (Oligo_2.PMC.Case - Oligo_2.PMC.Control) +
      (Oligo_2.PVC.Case - Oligo_2.PVC.Control)) /
      4,
  PD_effect_GPI_vs_DMNX_Oligo2 = (Oligo_2.GPI.Case - Oligo_2.GPI.Control) -
    (Oligo_2.DMNX.Case - Oligo_2.DMNX.Control),
  PD_effect_GPI_vs_PFC_Oligo2 = (Oligo_2.GPI.Case - Oligo_2.GPI.Control) -
    (Oligo_2.PFC.Case - Oligo_2.PFC.Control),
  PD_effect_GPI_vs_PMC_Oligo2 = (Oligo_2.GPI.Case - Oligo_2.GPI.Control) -
    (Oligo_2.PMC.Case - Oligo_2.PMC.Control),
  PD_effect_GPI_vs_PVC_Oligo2 = (Oligo_2.GPI.Case - Oligo_2.GPI.Control) -
    (Oligo_2.PVC.Case - Oligo_2.PVC.Control),
  levels = design
)

oligo_fit_contrasts <- contrasts.fit(oligo_fit, contrasts)
oligo_efit <- eBayes(oligo_fit_contrasts, robust = TRUE)

de_oligo_df <- colnames(oligo_efit$contrasts) %>%
  map_dfr(function(contrast) {
    topTable(
      oligo_efit,
      coef = contrast,
      number = Inf,
      adjust.method = "fdr",
      confint = TRUE
    ) %>%
      rownames_to_column("ensembl_gene_id") %>%
      rename(
        logfc = logFC,
        ci_lower = CI.L,
        ci_upper = CI.R,
        avg_expr = AveExpr,
        t_stat = t,
        p_value = P.Value,
        fdr = adj.P.Val,
        b_stat = B
      ) %>%
      mutate(contrast = contrast)
  }) %>%
  mutate(
    ensembl_gene_id = stringr::str_replace(ensembl_gene_id, "^MT-", "")
  ) %>%
  left_join(gene_metadata, by = dplyr::join_by(ensembl_gene_id)) %>%
  relocate(colnames(gene_metadata), .before = 1) %>%
  tibble::as_tibble()


# %% GSEA helpers
create_ranked_gene_list <- function(df) {
  gene_list_vector <- df$t_stat
  names(gene_list_vector) <- df$ensembl_gene_id
  sort(gene_list_vector, decreasing = TRUE)
}

run_gsea <- function(ranked_list) {
  GSEA(
    ranked_list,
    minGSSize = 10,
    maxGSSize = 500,
    TERM2GENE = MSIG_TERM2GENE,
    TERM2NAME = MSIG_TERM2NAME,
    pvalueCutoff = 1,
    eps = 1e-50,
    pAdjustMethod = "fdr",
    verbose = FALSE,
    seed = TRUE
  ) %>%
    as_tibble() %>%
    rename(gs_id = ID, gs_name = Description, gs_size = setSize) %>%
    left_join(MSIG_METADATA, by = "gs_id") %>%
    relocate(starts_with("gs"), .before = enrichmentScore)
}

# %% Load pre-computed GSEA results
gsea_pd_effect_oligo2_vs_oligo1_in_gpi <- read_rds(
  "04-analysis/04-gseas-of-interest/gsea-pd-effect-oligo2-vs-oligo1-in-gpi.rds"
)

gsea_bio_oligo2_vs_oligo1_in_gpi <- read_rds(
  "04-analysis/04-gseas-of-interest/gsea-bio-oligo2-vs-oligo1-in-gpi.rds"
)


# Figure 5b panel a: PD differential response volcano (GPI)
fdr_threshold <- 0.1
df <- de_oligo_df %>%
  filter(contrast == "PD_effect_Oligo2_vs_Oligo1_GPI") %>%
  mutate(
    sig = case_when(
      fdr > 0.1 ~ "ns",
      fdr <= 0.1 & logfc > 0 ~ "up",
      fdr <= 0.1 & logfc < 0 ~ "down",
      TRUE ~ "ns"
    ),
    symbol = if_else(is.na(hgnc_symbol), ensembl_gene_id, hgnc_symbol)
  )

p_pd_volcano <- df %>%
  ggplot(aes(x = logfc, y = fdr, fill = sig, shape = sig)) +
  scale_y_continuous(
    transform = neg_log10_transform,
    breaks = c(0.05, 0.1, 0.2, 0.5, 1)
  ) +
  scale_x_continuous(
    limits = c(-2.6, 2.6),
    breaks = seq(-2, 2, .5)
  ) +
  rasterise(
    geom_point(
      data = . %>% filter(sig == "ns"),
      size = .7, stroke = 0, alpha = .75
    ),
    dpi = dpi
  ) +
  geom_point(
    data = . %>% filter(sig != "ns"),
    size = 1, stroke = 0, alpha = .75
  ) +
  scale_fill_manual(
    values = c(ns = cool_grey, up = cool_red, down = cool_blue)
  ) +
  scale_shape_manual(values = c(ns = 21, up = 24, down = 25)) +
  geom_hline(
    yintercept = fdr_threshold,
    linetype = "dashed", linewidth = .1, alpha = 0.5
  ) +
  geom_vline(
    xintercept = c(-log2(1.5), log2(1.5)),
    linetype = "dashed", linewidth = .1, alpha = 0.5
  ) +
  ggrepel::geom_text_repel(
    data = df %>% filter(sig != "ns", logfc < 0),
    aes(label = symbol),
    seed = 1, segment.size = 0.1, min.segment.length = 0,
    size = 5 * 0.3528, color = "black", bg.color = "white", bg.r = 0.05,
    direction = "y",
    nudge_x = df %>%
      filter(sig != "ns", logfc < 0) %>%
      mutate(nudge = -1.5 - logfc) %>%
      pull(nudge),
    ylim = c(-log10(0.1), -log10(0.045)),
    hjust = 1, segment.curvature = 1e-20,
  ) +
  ggrepel::geom_text_repel(
    data = df %>% filter(sig != "ns", logfc > 0),
    aes(label = symbol),
    seed = 1, segment.size = 0.1, min.segment.length = 0,
    size = 5 * 0.3528, color = "black", bg.color = "white", bg.r = 0.05,
    direction = "y",
    nudge_x = df %>%
      filter(sig != "ns", logfc > 0) %>%
      mutate(nudge = 2 - logfc) %>%
      pull(nudge),
    ylim = c(-log10(0.1), -log10(0.045)),
    hjust = 0, segment.curvature = -1e-20,
  ) +
  labs(
    x = "log<sub>2</sub> Fold Change<br>(&Delta;PD Oligo2:Oligo1 in GPI)",
    y = "False Discovery Rate\n(FDR)",
  ) +
  theme_publication(5) +
  theme(axis.title.x = ggtext::element_markdown()) +
  no_legend

ggsave_pdf_svg(
  base_name = "05-results/02-plots/figure-5-pd-oligo2-vs-oligo1-gpi-volcanoplot",
  plot = p_pd_volcano, width = 80, height = 40, units = "mm",
)


# Figure 5b panels b, c: PD differential response combined GSEA
gsea_oligo2_up_terms <- gsea_pd_effect_oligo2_vs_oligo1_in_gpi %>%
  filter(
    !gs_subcollection %in%
      c("IMMUNESIGDB", "CGP", "CP:BIOCARTA", "HPO", "CP:PID")
  ) %>%
  filter(NES > 0, p.adjust < 0.1) %>%
  mutate(
    sign = "up",
    count = str_count(core_enrichment, "/") + 1,
    gene_ratio = count / gs_size,
    gs_cluster = case_when(
      str_detect(gs_name, "PROTEIN") ~ "Unfolded\nprotein",
      str_detect(gs_name, "FATTY|LIPID") ~ "Lipid\nmetabolism",
      str_detect(gs_name, "PH|ACIDIFICATI") ~ "Vesicular\nacidification",
      TRUE ~ "Other"
    ) %>%
      factor(
        levels = c(
          "Vesicular\nacidification",
          "Unfolded\nprotein",
          "Lipid\nmetabolism",
          "Other"
        )
      )
  )

gsea_oligo2_down_filtered <- gsea_pd_effect_oligo2_vs_oligo1_in_gpi %>%
  dplyr::filter(
    NES < 0,
    p.adjust < 0.05,
    !gs_subcollection %in%
      c("IMMUNESIGDB", "TFT:GTRD", "MIR:MIRDB", "CGP", "CP:BIOCARTA", "HPO", "CP:PID")
  ) %>%
  filter(!str_detect(gs_name, "(UV_RESPONSE)|(CELL_CORTEX)|TEMPERATURE|COLD")) %>%
  slice_min(p.adjust, n = 24)

down_clusters <- get_genesets_clusters(gsea_oligo2_down_filtered, n_clusters = 5)

gsea_oligo2_down_terms <- gsea_oligo2_down_filtered %>%
  left_join(down_clusters, by = "gs_name") %>%
  mutate(
    sign = "down",
    count = str_count(core_enrichment, "/") + 1,
    gene_ratio = count / gs_size,
    gs_cluster = forcats::fct_collapse(
      forcats::fct_relevel(gs_cluster, "3", "5", "4", "1", "2"),
      "Progenitor state\n(dedifferentiation)" = "1",
      "Neuronal\nsupport" = "2",
      "RNA metabolism" = "3",
      " " = "4",
      "Transcriptional\nregulation" = "5"
    )
  )

df_plot <- bind_rows(gsea_oligo2_up_terms, gsea_oligo2_down_terms) %>%
  mutate(
    gs_name_pretty = prettify_genesets_names(gs_name, width = 40, truncate = TRUE),
    gs_name_pretty = str_glue("{gs_name_pretty} [{gs_subcollection}]"),
  )

shared_scales <- list(
  scale_x_continuous(limits = range(df_plot$gene_ratio) + c(0, 0.1)),
  scale_color_viridis_c(
    option = "viridis",
    limits = range(-log10(df_plot$p.adjust)),
    direction = 1,
    name = expression(-log[10]("FDR"))
  ),
  scale_size_continuous(
    range = c(.5, 3),
    limits = range(df_plot$gs_size),
    name = "Gene\nset size"
  )
)

make_plot <- function(sign_label, title, show_legend = TRUE, hide_x = FALSE, box.color = "black") {
  ggplot(
    df_plot %>% filter(sign == sign_label),
    aes(x = gene_ratio, y = reorder(gs_name_pretty, gene_ratio))
  ) +
    shared_scales +
    geom_point(aes(size = gs_size, color = -log10(p.adjust))) +
    facet_grid(gs_cluster ~ ., scales = "free_y", space = "free_y") +
    geom_text(
      aes(label = scales::label_pvalue(accuracy = 0.01)(p.adjust)),
      hjust = 0, nudge_x = 0.04, size = 5 * 0.3528
    ) +
    labs(title = title, x = if (hide_x) "" else "Gene ratio", y = NULL) +
    theme_publication(5) +
    coord_cartesian(clip = "off") +
    vertical_scientific_colormap +
    theme(
      axis.text.x  = if (hide_x) element_blank() else element_text(),
      axis.line.x  = if (hide_x) element_blank() else element_line(),
      axis.ticks.x = if (hide_x) element_blank() else element_line(),
      legend.position = if (show_legend) "bottom" else "none",
      legend.box = "horizontal",
      panel.grid.major.y = element_line(color = "grey80", linewidth = .1),
      strip.clip = "off",
      plot.title = ggtext::element_textbox(
        hjust = 0.5, halign = 0.5, box.color = box.color,
        linetype = 1, linewidth = 0.2, color = "black",
        padding = margin(2, 2, 2, 2), margin = margin(b = 2, t = 1),
        width = unit(1, "npc")
      ),
      plot.title.position = "panel",
    )
}

p_up <- make_plot("up", title = "Enriched in Oligo 2", show_legend = FALSE, hide_x = TRUE,
                  box.color = correct_cell_subtype_colors["Oligo 2"])
p_down <- make_plot("down", title = "Enriched in Oligo 1", show_legend = TRUE, hide_x = FALSE,
                    box.color = correct_cell_subtype_colors["Oligo 1"])

n_up   <- nrow(gsea_oligo2_up_terms)
n_down <- nrow(gsea_oligo2_down_terms)
final_plot <- p_up / p_down + plot_layout(heights = c(n_up, n_down))

ggsave_pdf_svg(
  base_name = "05-results/02-plots/figure-5-pd-oligo2-vs-oligo1-combined-gseaplot",
  plot = final_plot, width = 90, height = 145, units = "mm"
)


# Figure 5a panel a: biological difference Oligo2 vs Oligo1 volcano (GPI)
fdr_threshold <- 1e-130
df <- de_oligo_df %>%
  filter(contrast == "Oligo2_vs_Oligo1_GPI") %>%
  mutate(
    sig = case_when(
      fdr > 1e-130 ~ "ns",
      fdr <= 1e-130 & logfc > 3 ~ "up",
      fdr <= 1e-130 & logfc < -3 ~ "down",
      TRUE ~ "ns"
    ),
    symbol = if_else(is.na(hgnc_symbol), ensembl_gene_id, hgnc_symbol)
  )

x_limit <- max(abs(df$logfc))

p_bio_volcano <- df %>%
  ggplot(aes(x = logfc, y = fdr, fill = sig, shape = sig)) +
  scale_y_continuous(
    transform = neg_log10_transform,
    breaks = c(1e-250, 1e-200, 1e-150, 1e-100, 1e-50, 1)
  ) +
  scale_x_continuous(limits = c(-6.5, 7.5), breaks = seq(-6, 6, 2)) +
  rasterise(
    geom_point(data = . %>% filter(sig == "ns"), size = .7, stroke = 0, alpha = .75),
    dpi = dpi
  ) +
  geom_point(data = . %>% filter(sig != "ns"), size = 1, stroke = 0, alpha = .75) +
  scale_fill_manual(values = c(ns = cool_grey, up = cool_red, down = cool_blue)) +
  scale_shape_manual(values = c(ns = 21, up = 24, down = 25)) +
  geom_hline(yintercept = fdr_threshold, linetype = "dashed", linewidth = .1, alpha = 0.5) +
  geom_vline(xintercept = c(-log2(8), log2(8)), linetype = "dashed", linewidth = .1, alpha = 0.5) +
  ggrepel::geom_text_repel(
    data = df %>% filter(sig != "ns", logfc > 0) %>% group_by(sig) %>% slice_head(n = 9),
    aes(label = symbol),
    seed = 1, segment.size = 0.1, min.segment.length = 0,
    size = 5 * 0.3528, color = "black", bg.color = "white", bg.r = 0.05,
    direction = "y",
    nudge_x = df %>%
      filter(sig != "ns", logfc > 0) %>%
      group_by(sig) %>%
      slice_head(n = 9) %>%
      ungroup() %>%
      mutate(nudge = 5.9 - logfc) %>%
      pull(nudge),
    ylim = c(100, 250), hjust = 0, segment.curvature = -1e-20,
  ) +
  ggrepel::geom_text_repel(
    data = df %>% filter(sig != "ns", logfc < 0) %>% group_by(sig) %>% slice_head(n = 10),
    aes(label = symbol),
    seed = 1, segment.size = 0.1, min.segment.length = 0,
    size = 5 * 0.3528, color = "black", bg.color = "white", bg.r = 0.05,
    direction = "y",
    nudge_x = df %>%
      filter(sig != "ns", logfc < 0) %>%
      group_by(sig) %>%
      slice_head(n = 10) %>%
      ungroup() %>%
      mutate(nudge = -5 - logfc) %>%
      pull(nudge),
    ylim = c(100, 250), hjust = 1, segment.curvature = 1e-20,
  ) +
  labs(
    y = "False Discovery Rate\n(FDR)",
    x = "log<sub>2</sub> Fold Change<br>(Oligo2:Oligo1 in GPI)"
  ) +
  theme_publication(5) +
  theme(axis.title.x = ggtext::element_markdown()) +
  no_legend

ggsave_pdf_svg(
  base_name = "05-results/02-plots/figure-5-bio-oligo2-vs-oligo1-gpi-volcanoplot",
  plot = p_bio_volcano, width = 80, height = 40, units = "mm",
)


# Figure 5a panels b, c: biological difference combined GSEA
n_clusters_up <- 6

gsea_bio_up_filtered <- gsea_bio_oligo2_vs_oligo1_in_gpi %>%
  filter(
    NES > 0, p.adjust < 0.05,
    !gs_subcollection %in% c("IMMUNESIGDB", "TFT:GTRD", "MIR:MIRDB", "CGP", "CP:BIOCARTA", "HPO", "CP:PID")
  ) %>%
  filter(!str_detect(gs_name, "SKIN|KERATIN")) %>%
  slice_min(p.adjust, n = 40) %>%
  filter(!str_detect(gs_name, "ABERRANT|ESCRT_III|BRUSH_BORDER|PROTEIN_HORMONE|BASED_MOVEMENT|MOTILE|LINEAR_AMIDINES|SULFUR_COMPOUND|PHOSPHOLIPID_METABOLISM|PROTEINOGENIC|ALPHA_AMINO_ACID|MOVEMENT|MOTILITY|AXONAL|CATSPER|MICROTUBULE_ASSOCIATED_COMPLEX"))

up_clusters <- get_genesets_clusters(gsea_bio_up_filtered, n_clusters_up)

gsea_bio_up_terms <- gsea_bio_up_filtered %>%
  left_join(up_clusters, by = "gs_name") %>%
  mutate(
    sign = "up",
    count = str_count(core_enrichment, "/") + 1,
    gene_ratio = count / gs_size,
    gs_cluster = forcats::fct_collapse(
      forcats::fct_relevel(gs_cluster, "1", "2", "3", "6", "4", "5"),
      "Lipid\nmetabolism" = "1",
      "Met\nmetabolism" = "2",
      "Cytoskeleton\n& transport" = c("3", "6"),
      "SNARE" = "4",
      " " = "5"
    )
  )

n_clusters_down <- 4

term_to_keep <- gsea_bio_oligo2_vs_oligo1_in_gpi %>%
  filter(gs_name == "GOBP_NEGATIVE_REGULATION_OF_OLIGODENDROCYTE_DIFFERENTIATION")

gsea_bio_down_filtered <- gsea_bio_oligo2_vs_oligo1_in_gpi %>%
  filter(
    NES < 0, p.adjust < 0.05,
    !gs_subcollection %in% c("IMMUNESIGDB", "TFT:GTRD", "MIR:MIRDB", "CGP", "CP:BIOCARTA", "HPO", "CP:PID")
  ) %>%
  filter(!str_detect(gs_name, "INFLUENZA|SECOND_MESSENGERS|COTRANSLATIONAL_PROTEIN|LYMPHO|TNFA|DNA_BINDING|SYNAPTIC|EIF2AK4|STARVATION")) %>%
  slice_min(p.adjust, n = 32) %>%
  filter(!str_detect(gs_name, "RRNA_PROCESSING|REACTOME_TRANSLATION|RRNA_METABOLIC_PROCESS|PROTEIN_RNA_COMPLEX|MRNA_METABOLIC_PROCESS|MRNA_BINDING|GOBP_MRNA_PROCESSING|REGULATION_OF_TRANSLATION|NUCLEAR_SPECK")) %>%
  bind_rows(term_to_keep) %>%
  distinct(gs_name, .keep_all = TRUE)

down_clusters <- get_genesets_clusters(gsea_bio_down_filtered, n_clusters_down)

gsea_bio_down_terms <- gsea_bio_down_filtered %>%
  left_join(down_clusters, by = "gs_name") %>%
  mutate(
    sign = "down",
    count = str_count(core_enrichment, "/") + 1,
    gene_ratio = count / gs_size,
    gs_cluster = forcats::fct_collapse(
      forcats::fct_relevel(gs_cluster, "2", "4", "1", "3"),
      "Cell\ngrowth" = "1",
      "Translation" = "2",
      "Glial diff." = "3",
      "RNA splicing" = "4"
    ),
  )

df_plot_bio <- bind_rows(gsea_bio_up_terms, gsea_bio_down_terms) %>%
  mutate(
    gs_name_pretty = prettify_genesets_names(gs_name, width = 40, truncate = TRUE),
    gs_name_pretty = str_glue("{gs_name_pretty} [{gs_subcollection}]"),
  )

max_log_fdr <- ceiling(max(-log10(df_plot_bio$p.adjust)))

shared_scales_bio <- list(
  scale_x_continuous(
    limits = range(df_plot_bio$gene_ratio) + c(0, 0.2),
    expand = expansion(mult = c(0.05, 0)),
    breaks = c(0.25, 0.5, 0.75, 1)
  ),
  scale_color_viridis_c(
    option = "viridis",
    limits = c(1, max_log_fdr),
    breaks = seq(1, max_log_fdr, by = 3),
    direction = 1,
    name = expression(-log[10]("FDR"))
  ),
  scale_size_continuous(
    range = c(.5, 3),
    limits = range(df_plot_bio$gs_size),
    name = "Gene\nset size"
  )
)

make_plot_bio <- function(sign_label, title, show_legend = TRUE, hide_x = FALSE, box.color = "black") {
  ggplot(
    df_plot_bio %>% filter(sign == sign_label),
    aes(x = gene_ratio, y = reorder(gs_name_pretty, gene_ratio))
  ) +
    shared_scales_bio +
    geom_point(aes(size = gs_size, color = -log10(p.adjust))) +
    facet_grid(gs_cluster ~ ., scales = "free_y", space = "free_y") +
    geom_text(
      aes(label = scales::label_pvalue(accuracy = 0.01)(p.adjust)),
      hjust = 0, nudge_x = 0.05, size = 5 * 0.3528
    ) +
    labs(title = title, x = if (hide_x) "" else "Gene ratio", y = NULL) +
    theme_publication(5) +
    coord_cartesian(clip = "off") +
    vertical_scientific_colormap +
    theme(
      axis.text.x  = if (hide_x) element_blank() else element_text(),
      axis.line.x  = if (hide_x) element_blank() else element_line(),
      axis.ticks.x = if (hide_x) element_blank() else element_line(),
      legend.position = if (show_legend) "bottom" else "none",
      legend.box = "horizontal",
      panel.grid.major.y = element_line(color = "grey80", linewidth = .1),
      strip.clip = "off",
      plot.title = ggtext::element_textbox(
        hjust = 0.5, halign = 0.5, box.color = box.color,
        linetype = 1, linewidth = 0.2, color = "black",
        padding = margin(2, 2, 2, 2), margin = margin(b = 2, t = 1),
        width = unit(1, "npc")
      ),
      plot.title.position = "panel",
    )
}

p_up_bio <- make_plot_bio("up", title = "Enriched in Oligo 2", show_legend = FALSE, hide_x = TRUE,
                          box.color = correct_cell_subtype_colors["Oligo 2"])
p_down_bio <- make_plot_bio("down", title = "Enriched in Oligo 1", show_legend = TRUE, hide_x = FALSE,
                            box.color = correct_cell_subtype_colors["Oligo 1"])

n_up_bio   <- nrow(gsea_bio_up_terms)
n_down_bio <- nrow(gsea_bio_down_terms)
final_plot_bio <- p_up_bio / p_down_bio + plot_layout(heights = c(n_up_bio, n_down_bio))

ggsave_pdf_svg(
  base_name = "05-results/02-plots/figure-5-bio-oligo2-vs-oligo1-combined-gseaplot",
  plot = final_plot_bio, width = 90, height = 145, units = "mm"
)
