library(Seurat)
library(dplyr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(harmony)

rm(list = ls()); gc()
ORIGINAL_DIR <- ""
setwd(ORIGINAL_DIR)

scdata <- readRDS(file = "scdata.RDS")
tam <- subset(scdata, subset = TAM_group == "T_TAM-LYVE1-like TAM")
tam <- NormalizeData(
  tam,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

tam <- FindVariableFeatures(
  tam,
  selection.method = "vst",
  nfeatures = 2000
)

tam <- ScaleData(
  tam,
  features = VariableFeatures(tam)
)

tam <- RunPCA(
  tam,
  features = VariableFeatures(tam),
  npcs = 50
)

ElbowPlot(tam, ndims = 50)

tam <- RunHarmony(
  tam,
  group.by.vars = "sample",
  theta = 4,
  lambda = 2
)

tam <- FindNeighbors(
  tam,
  reduction = "harmony",
  dims = 1:40
)

tam <- FindClusters(
  tam,
  resolution = 0.5
)

tam <- RunUMAP(
  tam,
  reduction = "harmony",
  dims = 1:40,
  n.neighbors = 10,
  min.dist = 0.5,
  spread = 0.5,
  seed.use = 123
)

tam$cluster <- Idents(tam)

DimPlot(tam, group.by = "cluster")+
  ggtitle("UMAP of Samples") +             
  theme(
    plot.title = element_text(size = 16),   
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 16),
    legend.text = element_text(size = 15)
  )

FeaturePlot(
  tam,
  features = "TCN2",
  reduction = "umap",
  order = TRUE
) +
  ggtitle("TCN2 Expression") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15)
  )

DotPlot(tam, features = "TCN2", group.by = "cluster")+
  RotatedAxis() +
  scale_x_discrete(labels = "TCN2")+
  scale_color_gradientn(colors = c("blue", "white", "red"))+
  theme(
    axis.text.x = element_text(size = 14, angle = 45),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )+
  xlab("Marker Genes") +
  ylab("Cell Type")

table(Idents(tam))
new.cluster.ids <- c("TCN2-high","TCN2-low","TCN2-low","TCN2-low","TCN2-high","TCN2-low",
                     "TCN2-low","TCN2-high","TCN2-high")
names(new.cluster.ids) <- levels(tam)
tam <- RenameIdents(tam, new.cluster.ids)
tam$tam_group <- Idents(tam)
DimPlot(tam, group.by = "tam_group")+
  ggtitle("UMAP of Samples") +             
  theme(
    plot.title = element_text(size = 16),   
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 16),
    legend.text = element_text(size = 15)
  )
summary(tam$tam_group)

tam_group <- data.frame(
  barcode = colnames(tam),
  tam_group = tam$tam_group
)
save(tam_group, file = "TCN2_group_result_60samples.Rdata")
saveRDS(tam, file = "tam_group_60samples.RDS")

#--------S5A--------
tam <- readRDS(file = tam_group_60samples.RDS")
tam <- RenameIdents(
  tam,
  "TCN2-low" = "TCN2-low TAM_LYVE1",
  "TCN2-high" = "TCN2-high TAM_LYVE1"
)

my_colors <- c(
  "TCN2-low TAM_LYVE1" = "#A6CEE3",
  "TCN2-high TAM_LYVE1" = "#8B5A9E"
)

p1 <- DimPlot(
  tam,
  cols = my_colors,
  label.size = 5
) +
  labs(
    x = "UMAP1",
    y = "UMAP2"
  ) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    plot.title = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal"
  ) +
  guides(
    color = guide_legend(
      ncol = 1,
      override.aes = list(size = 5)
    )
  )

p1

ggsave(
  file.path("FigS5A-UMAP.png"),
  p1,
  width = 3.5,
  height = 4,
  dpi = 300,
  bg = "white"
)

#--------S5B--------
#通路分析
tam <- readRDS(file = "tam_group_60samples.RDS")
options(stringsAsFactors = FALSE)
set.seed(20260713)

input_gmt <- "c2.cp.reactome.v2026.1.Hs.symbols.gmt"

tam_group <- data.frame(
  barcode = colnames(tam),
  tam_group = tam$tam_group
)
tam_group <- tam_group[colnames(tam),]
tam$TCN2_group <- tam_group$tam_group
tam$TCN2_group <- factor(tam$TCN2_group, levels = c("TCN2-low", "TCN2-high"))
group_counts <- as.data.frame(table(tam$TCN2_group))
colnames(group_counts) <- c("group", "n_cells")

expr <- GetAssayData(tam, assay = "RNA", layer = "data")
de <- presto::wilcoxauc(expr, y = tam$TCN2_group)
de_high <- de[de$group == "TCN2-high", , drop = FALSE]
de_high$avg_log2FC <- log2((de_high$avgExpr + 1e-9) /
                             (de[de$group == "TCN2-low", "avgExpr"] + 1e-9))

rank_keep <- de_high$pct_in >= 1 | de_high$pct_out >= 1
rank_stats <- 2 * (de_high$auc[rank_keep] - 0.5) +
  de_high$logFC[rank_keep] * 1e-7

names(rank_stats) <- de_high$feature[rank_keep]
rank_stats <- sort(rank_stats[is.finite(rank_stats)], decreasing = TRUE)

pathways <- fgsea::gmtPathways(input_gmt)
gsea <- fgsea::fgseaMultilevel(
  pathways = pathways,
  stats = rank_stats,
  minSize = 10,
  maxSize = 500,
  eps = 0,
  nPermSimple = 10000
)
gsea <- as.data.frame(gsea)
if (anyNA(gsea$pval)) {
  na_paths <- gsea$pathway[is.na(gsea$pval)]
  fallback <- suppressWarnings(fgsea::fgseaSimple(
    pathways = pathways[na_paths], stats = rank_stats,
    nperm = 10000, minSize = 10, maxSize = 500
  ))
  fallback <- as.data.frame(fallback)
  replace_cols <- intersect(c("pval", "padj", "ES", "NES", "size", "leadingEdge"),
                            colnames(fallback))
  for (p in na_paths) {
    i <- match(p, gsea$pathway)
    j <- match(p, fallback$pathway)
    if (!is.na(j)) gsea[i, replace_cols] <- fallback[j, replace_cols]
  }
  gsea$padj <- p.adjust(gsea$pval, method = "BH")
}
gsea$enriched_in <- "TCN2-high"
gsea$leadingEdge <- vapply(gsea$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
gsea_high <- gsea[order(gsea$padj, -abs(gsea$NES)), ]
gsea_sig <- gsea[!is.na(gsea$padj) & gsea$padj < 0.05, ]
write.csv(gsea_sig, file = gsea_result_60samples.csv", row.names = F)


#计算low组
options(stringsAsFactors = FALSE)
set.seed(20260713)

input_gmt <- "c2.cp.reactome.v2026.1.Hs.symbols.gmt"

tam_group <- data.frame(
  barcode = colnames(tam),
  tam_group = tam$tam_group
)
tam_group <- tam_group[colnames(tam),]
tam$TCN2_group <- tam_group$tam_group
tam$TCN2_group <- factor(tam$TCN2_group, levels = c("TCN2-low", "TCN2-high"))
group_counts <- as.data.frame(table(tam$TCN2_group))
colnames(group_counts) <- c("group", "n_cells")

expr <- GetAssayData(tam, assay = "RNA", layer = "data")
de <- presto::wilcoxauc(expr, y = tam$TCN2_group)
de_low <- de[de$group == "TCN2-low", , drop = FALSE]
de_low$avg_log2FC <- log2((de_low$avgExpr + 1e-9) /
                             (de[de$group == "TCN2-low", "avgExpr"] + 1e-9))

rank_keep <- de_low$pct_in >= 1 | de_low$pct_out >= 1
rank_stats <- 2 * (de_low$auc[rank_keep] - 0.5) +
  de_low$logFC[rank_keep] * 1e-7

names(rank_stats) <- de_low$feature[rank_keep]
rank_stats <- sort(rank_stats[is.finite(rank_stats)], decreasing = TRUE)

pathways <- fgsea::gmtPathways(input_gmt)
gsea <- fgsea::fgseaMultilevel(
  pathways = pathways,
  stats = rank_stats,
  minSize = 10,
  maxSize = 500,
  eps = 0,
  nPermSimple = 10000
)
gsea <- as.data.frame(gsea)
if (anyNA(gsea$pval)) {
  na_paths <- gsea$pathway[is.na(gsea$pval)]
  fallback <- suppressWarnings(fgsea::fgseaSimple(
    pathways = pathways[na_paths], stats = rank_stats,
    nperm = 10000, minSize = 10, maxSize = 500
  ))
  fallback <- as.data.frame(fallback)
  replace_cols <- intersect(c("pval", "padj", "ES", "NES", "size", "leadingEdge"),
                            colnames(fallback))
  for (p in na_paths) {
    i <- match(p, gsea$pathway)
    j <- match(p, fallback$pathway)
    if (!is.na(j)) gsea[i, replace_cols] <- fallback[j, replace_cols]
  }
  gsea$padj <- p.adjust(gsea$pval, method = "BH")
}
gsea$enriched_in <- "TCN2-low"
gsea$leadingEdge <- vapply(gsea$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
gsea_low <- gsea[order(gsea$padj, -abs(gsea$NES)), ]
gsea_sig_low <- gsea[!is.na(gsea$padj) & gsea$padj < 0.05, ]

pathway_keep <- c(
  "REACTOME_INTERFERON_GAMMA_SIGNALING",
  "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING",
  "REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION",
  "REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION",
  "REACTOME_ANTIGEN_PRESENTATION_FOLDING_ASSEMBLY_AND_PEPTIDE_LOADING_OF_CLASS_I_MHC",
  "REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL",
  "REACTOME_INTERLEUKIN_10_SIGNALING"
)


pathway_order <- c(
  "REACTOME_INTERFERON_GAMMA_SIGNALING",
  "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING",
  "REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION",
  "REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION",
  "REACTOME_ANTIGEN_PRESENTATION_FOLDING_ASSEMBLY_AND_PEPTIDE_LOADING_OF_CLASS_I_MHC",
  "REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL",
  "REACTOME_INTERLEUKIN_10_SIGNALING"
)

pathway_labels <- c(
  "REACTOME_INTERFERON_GAMMA_SIGNALING" =
    "IFN-gamma signaling",
  
  "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING" =
    "IFN-alpha/beta signaling",
  
  "REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION" =
    "Cross-presentation",
  
  "REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION" =
    "MHC-II antigen presentation",
  
  "REACTOME_ANTIGEN_PRESENTATION_FOLDING_ASSEMBLY_AND_PEPTIDE_LOADING_OF_CLASS_I_MHC" =
    "MHC-I antigen presentation",
  
  "REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL" =
    "Lymphoid-myeloid immunoregulation",
  
  "REACTOME_INTERLEUKIN_10_SIGNALING" =
    "IL-10 signaling"
)
gsea <- rbind(gsea_high, gsea_low)
plot_df <- gsea %>%
  filter(pathway %in% pathway_keep) %>%
  dplyr::select(
    pathway,
    enriched_in,
    padj,
    NES
  ) %>%
  mutate(
    log10P = -log10(pmax(padj, 1e-300)),
    pathway = factor(
      pathway,
      levels = pathway_order
    ),
    enriched_in = factor(
      enriched_in,
      levels = c("TCN2-low", "TCN2-high")
    )
  )
print(plot_df)

p1 <- ggplot(
  plot_df,
  aes(
    x = pathway,
    y = enriched_in
  )
) +
  geom_point(
    aes(
      size = log10P,
      fill = NES
    ),
    shape = 21,
    color = "black",
    stroke = 0.6
  ) +
  scale_fill_gradient2(
    name = "NES",
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    breaks = pretty_breaks(5)
  ) +
  scale_size_continuous(
    name = expression(-log[10](P[adj])),
    range = c(3, 10),
    breaks = pretty_breaks(4)
  ) +
  scale_x_discrete(
    position = "top",
    labels = pathway_labels,
    drop = FALSE
  ) +
  scale_y_discrete(
    position = "right",
    drop = FALSE
  ) +
  coord_cartesian(clip = "off") +
  theme_bw() +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme(
    panel.border = element_rect(
      fill = NA,
      color = "black",
      linewidth = 1
    ),
    panel.grid.major = element_line(
      color = "grey90",
      linewidth = 0.4
    ),
    panel.grid.minor = element_blank(),
    
    axis.text.x = element_text(
      size = 12,
      angle = 55,
      colour = "black",
      hjust = 0,
      vjust = 0
    ),
    axis.text.y = element_text(
      size = 13,
      colour = "black"
    ),
    
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    legend.position = "bottom",
    legend.box = "horizontal",
    
    plot.margin = margin(
      t = 5,
      r = 5,
      b = 2,
      l = 5,
      unit = "pt"
    )
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      order = 1,
      barwidth = unit(4, "cm")
    ),
    size = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      order = 2
    )
  )

p1
ggsave(file.path("FigS5B-heatmap.pdf"), p1, width =7, height = 4.3, bg = "white",
       device = cairo_pdf)
