#!/usr/bin/env Rscript
# ==============================================================================
# BENG 628: Lecture 2 — Nature 2025 Figure 1 Replication
# Panels: 1D (PCA), 1E (Venn), 1F (Heatmap), 1H (Volcano), 1I (SynGO Lollipop)
#
# Paper : https://www.nature.com/articles/s41586-025-09987-9
# Data  : PRIDE Archive PXD057020
#         P1_MS2.raw   → Cam;PheRS*  (Fig 1D, F, H, I)
#         Met_MS2.raw  → Cam;MetRS*  (Fig 1D, E, F)
#         Tyr_MS2.raw  → Cam;TyrRS*  (negative control)
#
# NOTE: Values below are seeded simulations that approximate the published
#       distributions. Replace with real MaxQuant/Spectronaut output to
#       reproduce the exact paper figures.
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------------------------
required_cran <- c("tidyverse", "ggpubr", "viridis", "ggrepel", "ggvenn", "ggplotify")
new_pkg <- required_cran[!(required_cran %in% installed.packages()[, "Package"])]
if (length(new_pkg)) install.packages(new_pkg, repos = "https://cloud.r-project.org")

if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("ComplexHeatmap", ask = FALSE)
}
if (!requireNamespace("circlize", quietly = TRUE)) {
  BiocManager::install("circlize", ask = FALSE)
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(viridis)
  library(ggrepel)
  library(ggvenn)
  library(ComplexHeatmap)
  library(circlize)
  library(ggplotify)
  library(grid)
})

# Output directory — change to your preferred path
OUT_DIR <- "."
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

MODEL_COLORS <- c(
  "Cam;PheRS*" = "#56B4E9",
  "Cam;MetRS*" = "#E69F00",
  "Cam;TyrRS*" = "#009E73"
)

cat("=== BENG 628 Lecture 2: Figure 1 Panels ===\n\n")


# ==============================================================================
# Fig 1D — PCA
# Three-replicate PCA of LFQ intensities across the three synthetase models.
# Each model clusters tightly; PheRS* separates strongly along PC1.
# ==============================================================================
cat("Rendering Fig 1D: PCA...\n")

set.seed(123)
pca_df <- data.frame(
  PC1   = c(rnorm(3,  18, 1.5), rnorm(3, -14, 1.5), rnorm(3, -16, 1.0)),
  PC2   = c(rnorm(3,   0, 1.5), rnorm(3,   8, 1.5), rnorm(3,  -8, 1.5)),
  Model = factor(rep(c("Cam;PheRS*", "Cam;MetRS*", "Cam;TyrRS*"), each = 3))
)

# Draw ellipses only when ≥ 3 points exist per group (stat_ellipse default)
p1d <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Model)) +
  stat_ellipse(aes(fill = Model), geom = "polygon",
               level = 0.90, alpha = 0.08, show.legend = FALSE) +
  geom_point(size = 5, alpha = 0.9) +
  scale_color_manual(values = MODEL_COLORS) +
  scale_fill_manual(values  = MODEL_COLORS) +
  labs(
    title    = "Principal Component Analysis",
    subtitle = "LFQ intensities — triplicates per model",
    x        = "PC1 (81.0%)",
    y        = "PC2 (17.2%)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold"),
    legend.title  = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(OUT_DIR, "Fig1D_PCA.pdf"),
       p1d, width = 5, height = 4, useDingbats = FALSE)


# ==============================================================================
# Fig 1E — Venn Diagram
# Protein identifications shared/unique across the three models.
# Cam;TyrRS* (negative control) contributes almost nothing unique.
# ==============================================================================
cat("Rendering Fig 1E: Venn diagram...\n")

# Protein IDs approximating published overlap numbers
venn_list <- list(
  `Cam;PheRS*` = paste0("P", 1:3787),
  `Cam;MetRS*` = paste0("P", 1824:4143),   # ~1 964 overlap with PheRS*
  `Cam;TyrRS*` = paste0("P", c(1, 10, 100, 1000))  # 4 proteins (negative ctrl)
)

p1e <- ggvenn(
  venn_list,
  fill_color   = unname(MODEL_COLORS),
  stroke_size  = 0.5,
  set_name_size = 3.8,
  text_size     = 3.5
) +
  labs(title = "Protein Identification Overlap") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12))

ggsave(file.path(OUT_DIR, "Fig1E_Venn.pdf"),
       p1e, width = 5, height = 4.5, useDingbats = FALSE)


# ==============================================================================
# Fig 1F — Clustered Heatmap
# -log10(p-value) enrichment matrix across the three models (rows = proteins).
# Cam;PheRS* shows highest and broadest enrichment; TyrRS* baseline noise.
# ==============================================================================
cat("Rendering Fig 1F: Clustered heatmap...\n")

set.seed(42)
n_proteins <- 100  # increase to ~3 000 for full-scale figure

# Simulate enrichment scores: PheRS* high, MetRS* moderate, TyrRS* near-zero
heat_mat <- matrix(
  c(
    rnorm(n_proteins, mean = 5.5, sd = 2),   # Cam;PheRS*
    rnorm(n_proteins, mean = 3.0, sd = 2),   # Cam;MetRS*
    rnorm(n_proteins, mean = 0.5, sd = 0.8)  # Cam;TyrRS*
  ),
  ncol = 3
)
heat_mat <- pmax(heat_mat, 0)   # -log10 p cannot be negative
colnames(heat_mat) <- c("Cam;PheRS*", "Cam;MetRS*", "Cam;TyrRS*")

col_fun <- colorRamp2(
  breaks = c(0, 5, 10),
  colors = c("white", "orange", "darkred")
)

col_anno <- HeatmapAnnotation(
  Model = colnames(heat_mat),
  col   = list(Model = MODEL_COLORS),
  show_legend = FALSE,
  show_annotation_name = FALSE
)

hm_obj <- Heatmap(
  heat_mat,
  name             = "-log10[P]",
  col              = col_fun,
  top_annotation   = col_anno,
  cluster_columns  = TRUE,
  cluster_rows     = TRUE,
  show_row_names   = FALSE,
  column_title     = "Enrichment Heatmap (Fig 1F)",
  column_title_gp  = gpar(fontsize = 12, fontface = "bold"),
  column_names_gp  = gpar(fontsize = 10),
  heatmap_legend_param = list(title_gp = gpar(fontsize = 9),
                              labels_gp = gpar(fontsize = 8))
)

# Capture ComplexHeatmap as a ggplot-compatible grob
p1f <- as.ggplot(grid.grabExpr(draw(hm_obj), wrap = TRUE))

ggsave(file.path(OUT_DIR, "Fig1F_Heatmap.pdf"),
       p1f, width = 4, height = 6, useDingbats = FALSE)


# ==============================================================================
# Fig 1H — Volcano Plot
# PheRS*-enriched proteins plotted by log2 signal vs. -log10 p-value.
# Key synaptic proteins (SNCA, REEP1, LY6H) are labelled.
# ==============================================================================
cat("Rendering Fig 1H: Volcano plot...\n")

set.seed(7)

# Key labelled proteins
key_proteins <- data.frame(
  gene    = c("SNCA",  "REEP1", "SACS", "LY6H", "ARF3"),
  log2FC  = c(5.2,     4.8,     0.8,    6.1,    4.5),
  log10p  = c(6.5,     7.2,     4.1,    8.9,    7.8),
  Type    = c("Neuron","Neuron","Glia","Neuron","Neuron")
)

# Background proteome
bg_proteins <- data.frame(
  gene   = paste0("P", seq_len(1000)),
  log2FC = rnorm(1000, mean = 3.0, sd = 1.5),
  log10p = pmax(0, rnorm(1000, mean = 5.0, sd = 2.0)),
  Type   = "Other"
)

volcano_df <- bind_rows(bg_proteins, key_proteins)

# Significance thresholds
FC_THRESH  <- 2.0   # log2 fold-change cutoff
P_THRESH   <- 3.0   # -log10 p-value cutoff

p1h <- ggplot(volcano_df, aes(x = log2FC, y = log10p, color = Type)) +
  # Threshold lines
  geom_hline(yintercept = P_THRESH, linetype = "dashed",
             color = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = FC_THRESH, linetype = "dashed",
             color = "grey50", linewidth = 0.4) +
  # Points
  geom_point(data = filter(volcano_df, Type == "Other"),
             alpha = 0.4, size = 1.5) +
  geom_point(data = filter(volcano_df, Type != "Other"),
             alpha = 0.9, size = 3.0) +
  # Labels for key proteins
  geom_text_repel(
    data          = filter(volcano_df, gene %in% c("SNCA", "REEP1", "LY6H", "SACS")),
    aes(label = gene),
    color         = "black",
    fontface      = "bold",
    size          = 3.5,
    box.padding   = 0.5,
    point.padding = 0.3,
    max.overlaps  = 20,
    segment.color = "grey60"
  ) +
  scale_color_manual(
    values = c("Neuron" = "#D73027", "Glia" = "#E69F00", "Other" = "grey75"),
    guide  = guide_legend(override.aes = list(size = 3, alpha = 1))
  ) +
  labs(
    title    = "PheRS* Synaptic Protein Enrichment",
    subtitle = "Cam;PheRS* vs. input (n = 3 replicates)",
    x        = "log\u2082 [Signal Ratio]",
    y        = "-log\u2081\u2080 [P-value]",
    color    = "Cell type"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold"),
    legend.position = c(0.02, 0.97),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "white", colour = "grey80", linewidth = 0.3)
  )

ggsave(file.path(OUT_DIR, "Fig1H_Volcano.pdf"),
       p1h, width = 6, height = 5, useDingbats = FALSE)


# ==============================================================================
# Fig 1I — SynGO Compartment Enrichment (Lollipop / Dot Plot)
# Dot size = fold-enrichment; colour = enrichment magnitude (viridis magma).
# Ordered by -log10[FDR] ascending so the strongest hit sits at the top.
# ==============================================================================
cat("Rendering Fig 1I: SynGO lollipop...\n")

syngo_df <- data.frame(
  Compartment = c("Postsynapse", "Presynapse", "Dendrite",
                  "Myelin sheath", "Axon", "Synapse"),
  FDR         = c(28, 32, 38, 52, 42, 65),     # -log10[FDR]
  Enrichment  = c(2.5, 2.3, 2.1, 1.8, 2.4, 2.8)
) %>%
  mutate(Compartment = fct_reorder(Compartment, FDR))

p1i <- ggplot(syngo_df, aes(x = FDR, y = Compartment)) +
  geom_segment(aes(xend = 0, yend = Compartment),
               color = "grey70", linewidth = 0.8) +
  geom_point(aes(size = Enrichment, color = Enrichment)) +
  scale_color_viridis_c(
    option = "magma",
    name   = "Fold\nenrichment",
    limits = c(1.5, 3.0),
    guide  = guide_colorbar(barwidth = 0.6, barheight = 5,
                            title.position = "top")
  ) +
  scale_size_continuous(
    range  = c(4, 10),
    name   = "Fold\nenrichment",
    limits = c(1.5, 3.0),
    guide  = guide_legend(title.position = "top")
  ) +
  # Merge the two legends into one
  guides(
    color = guide_colorbar(barwidth = 0.6, barheight = 5,
                           title.position = "top"),
    size  = "none"
  ) +
  labs(
    title    = "SynGO Cellular Compartment Enrichment",
    subtitle = "Cam;PheRS*-enriched synaptic proteome",
    x        = "-log\u2081\u2080 [FDR]",
    y        = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "right"
  )

ggsave(file.path(OUT_DIR, "Fig1I_SynGO.pdf"),
       p1i, width = 6, height = 4, useDingbats = FALSE)


# ==============================================================================
# Combined layout (matches paper Fig 1 arrangement: D | E | F  then  H | I)
# ==============================================================================
cat("Rendering combined Figure 1 layout...\n")

top_row    <- ggarrange(p1d, p1e, p1f,  ncol = 3, labels = c("d", "e", "f"))
bottom_row <- ggarrange(p1h, p1i,       ncol = 2, labels = c("h", "i"),
                        widths = c(1.3, 1))
fig1_full  <- ggarrange(top_row, bottom_row, nrow = 2)

ggsave(file.path(OUT_DIR, "Fig1_Combined.pdf"),
       fig1_full, width = 14, height = 10, useDingbats = FALSE)

cat("\nDone. Output files saved to:", normalizePath(OUT_DIR), "\n")
cat("  Fig1D_PCA.pdf\n")
cat("  Fig1E_Venn.pdf\n")
cat("  Fig1F_Heatmap.pdf\n")
cat("  Fig1H_Volcano.pdf\n")
cat("  Fig1I_SynGO.pdf\n")
cat("  Fig1_Combined.pdf\n")
