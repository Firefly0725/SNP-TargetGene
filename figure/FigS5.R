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

# -------------------------- TAM_LYVE1 preprocessing --------------------------
scdata <- readRDS(file = "scdata.RDS")
tam <- subset(scdata, subset = TAM_group == "T_TAM-LYVE1-like TAM")
tam <- NormalizeData(tam, normalization.method = "LogNormalize", scale.factor = 10000)
tam <- FindVariableFeatures(tam, selection.method = "vst", nfeatures = 2000)
tam <- ScaleData(tam, features = VariableFeatures(tam))
tam <- RunPCA(tam, features = VariableFeatures(tam), npcs = 50)
ElbowPlot(tam, ndims = 50)
tam <- RunHarmony(tam, group.by.vars = "sample", theta = 4, lambda = 2)
tam <- FindNeighbors(tam, reduction = "harmony", dims = 1:40)
tam <- FindClusters(tam, resolution = 0.5)
tam <- RunUMAP(tam, reduction = "harmony", dims = 1:40, n.neighbors = 10,
               min.dist = 0.5, spread = 0.5, seed.use = 123)
tam$cluster <- Idents(tam)

DimPlot(tam, group.by = "cluster") + ggtitle("UMAP of Samples") +
  theme(plot.title = element_text(size = 16), axis.title = element_text(size = 16),
        axis.text = element_text(size = 16), legend.text = element_text(size = 15))

FeaturePlot(tam, features = "TCN2", reduction = "umap", order = TRUE) +
  ggtitle("TCN2 Expression") + theme_classic() +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        axis.title = element_text(size = 16), axis.text = element_text(size = 16),
        legend.text = element_text(size = 15), legend.title = element_text(size = 15))

DotPlot(tam, features = "TCN2", group.by = "cluster") + RotatedAxis() +
  scale_x_discrete(labels = "TCN2") +
  scale_color_gradientn(colors = c("blue", "white", "red")) +
  theme(axis.text.x = element_text(size = 14, angle = 45),
        axis.text.y = element_text(size = 16), axis.title = element_text(size = 16),
        legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
  labs(x = "Marker Genes", y = "Cell Type")

table(Idents(tam))
new.cluster.ids <- c("TCN2-high", "TCN2-low", "TCN2-low", "TCN2-low",
                     "TCN2-high", "TCN2-low", "TCN2-low", "TCN2-high", "TCN2-high")
names(new.cluster.ids) <- levels(tam)
tam <- RenameIdents(tam, new.cluster.ids)
tam$tam_group <- Idents(tam)

tam_group <- data.frame(barcode = colnames(tam), tam_group = tam$tam_group)
save(tam_group, file = "TCN2_group_result_60samples.Rdata")
saveRDS(tam, file = "tam_group_60samples.RDS")

# Common files and appearance used by Fig. S5A-C
rds_file <- "tam_group_60samples.RDS"
input_gmt <- "c2.cp.reactome.v2026.1.Hs.symbols.gmt"
group_levels <- c("TCN2-low", "TCN2-high")
group_labels <- paste(group_levels, "TAM_LYVE1")
group_colors <- setNames(c("#A6CEE3", "#8B5A9E"), group_labels)
common_theme <- theme(
  axis.title = element_text(size = 14), axis.text = element_text(size = 14),
  legend.title = element_text(size = 14), legend.text = element_text(size = 14),
  legend.position = "bottom", legend.direction = "horizontal",
  plot.title = element_blank()
)

# Shared pathway set and plotting style for Fig. S5B-C
pathways <- fgsea::gmtPathways(input_gmt)
run_fgsea <- function(stats) {
  ans <- as.data.frame(fgsea::fgseaMultilevel(
    pathways, stats, minSize = 10, maxSize = 500, eps = 0,
    nPermSimple = 10000, nproc = 1
  ))
  if (anyNA(ans$pval)) {
    bad <- ans$pathway[is.na(ans$pval)]
    fb <- as.data.frame(suppressWarnings(fgsea::fgseaSimple(
      pathways[bad], stats, nperm = 10000, minSize = 10, maxSize = 500,
      nproc = 1
    )))
    cols <- intersect(c("pval", "ES", "NES", "size", "leadingEdge"), names(fb))
    for (p in bad) ans[match(p, ans$pathway), cols] <- fb[match(p, fb$pathway), cols]
    ans$padj <- p.adjust(ans$pval, method = "BH")
  }
  ans
}
pathway_order <- c(
  "REACTOME_INTERFERON_GAMMA_SIGNALING",
  "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING",
  "REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION",
  "REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION",
  "REACTOME_ANTIGEN_PRESENTATION_FOLDING_ASSEMBLY_AND_PEPTIDE_LOADING_OF_CLASS_I_MHC",
  "REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL",
  "REACTOME_INTERLEUKIN_10_SIGNALING"
)
pathway_labels <- setNames(c(
  "IFN-gamma signaling", "IFN-alpha/beta signaling", "Cross-presentation",
  "MHC-II antigen presentation", "MHC-I antigen presentation",
  "Lymphoid-myeloid immunoregulation", "IL-10 signaling"
), pathway_order)

format_gsea_data <- function(x) {
  x %>% filter(pathway %in% pathway_order) %>%
    dplyr::select(pathway, enriched_in, padj, NES) %>%
    mutate(log10P = -log10(pmax(padj, 1e-300)),
           pathway = factor(pathway, levels = pathway_order),
           enriched_in = factor(enriched_in, levels = group_levels))
}

gsea_dotplot <- function(x) {
  ggplot(x, aes(pathway, enriched_in)) +
    geom_point(aes(size = log10P, fill = NES), shape = 21,
               color = "black", stroke = 0.6) +
    scale_fill_gradient2(name = "NES", low = "#2166AC", mid = "white",
                         high = "#B2182B", midpoint = 0,
                         breaks = scales::pretty_breaks(5)) +
    scale_size_continuous(name = expression(-log[10](P[adj])), range = c(3, 10),
                          breaks = scales::pretty_breaks(4)) +
    scale_x_discrete(position = "top", labels = pathway_labels, drop = FALSE) +
    scale_y_discrete(position = "right", drop = FALSE) +
    coord_cartesian(clip = "off") + labs(x = NULL, y = NULL) + theme_bw() +
    theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
          panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = 12, angle = 55, colour = "black",
                                     hjust = 0, vjust = 0),
          axis.text.y = element_text(size = 13, colour = "black"),
          legend.title = element_text(size = 13), legend.text = element_text(size = 11),
          legend.position = "bottom", legend.box = "horizontal",
          plot.margin = margin(5, 5, 2, 5, unit = "pt")) +
    guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5,
                                 order = 1, barwidth = grid::unit(4, "cm")),
           size = guide_legend(title.position = "top", title.hjust = 0.5, order = 2))
}

# ---------------------------------- Fig. S5A ---------------------------------
tam <- readRDS(file = rds_file)
tam <- RenameIdents(tam, "TCN2-low" = group_labels[1], "TCN2-high" = group_labels[2])

pA <- DimPlot(tam, cols = group_colors, label.size = 5) +
  labs(x = "UMAP1", y = "UMAP2") + common_theme +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 5)))

ggsave("FigS5A-UMAP.png", pA, width = 3.5, height = 4, dpi = 300, bg = "white")

# ---------------------------------- Fig. S5B ---------------------------------
options(stringsAsFactors = FALSE)
set.seed(20260713)
tam <- readRDS(file = rds_file)
tam$TCN2_group <- factor(as.character(tam$tam_group), levels = group_levels)
expr <- GetAssayData(tam, assay = "RNA", layer = "data")
cell_de <- presto::wilcoxauc(expr, y = tam$TCN2_group)

run_cell_gsea <- function(label) {
  z <- cell_de[cell_de$group == label, , drop = FALSE]
  z <- z[z$pct_in >= 1 | z$pct_out >= 1, ]
  stats <- sort(setNames(2 * (z$auc - 0.5) + z$logFC * 1e-7, z$feature),
                decreasing = TRUE)
  ans <- run_fgsea(stats)
  ans$enriched_in <- label
  ans
}

gsea_B <- bind_rows(run_cell_gsea("TCN2-high"), run_cell_gsea("TCN2-low"))
gsea_B_out <- gsea_B
gsea_B_out$leadingEdge <- vapply(gsea_B_out$leadingEdge, paste, collapse = ";",
                                 FUN.VALUE = character(1))
write.csv(gsea_B_out[!is.na(gsea_B_out$padj) & gsea_B_out$padj < 0.05, ],
          "gsea_result_60samples.csv", row.names = FALSE)

plot_B <- format_gsea_data(gsea_B)
pB <- gsea_dotplot(plot_B)
ggsave("FigS5B-heatmap.pdf", pB, width = 7, height = 4.3,
       bg = "white", device = cairo_pdf)
ggsave("FigS5B-heatmap.png", pB, width = 7, height = 4.3,
       dpi = 300, bg = "white")

# ---------------------------------- Fig. S5C ---------------------------------
options(stringsAsFactors = FALSE)
set.seed(20260713)
tam <- readRDS(file = rds_file)
meta <- tam[[]]
meta$donor <- as.character(meta$sample)
meta$TCN2_group <- factor(as.character(meta$tam_group), levels = group_levels)
meta$pb_id <- paste(meta$donor, meta$TCN2_group, sep = "__")

# Sum raw counts within each patient x TCN2 group.
counts <- GetAssayData(tam, assay = "RNA", layer = "counts")
meta <- meta[colnames(counts), , drop = FALSE]
pb_levels <- unique(meta$pb_id)
cell_to_pb <- Matrix::sparseMatrix(
  i = seq_len(nrow(meta)), j = match(meta$pb_id, pb_levels), x = 1,
  dims = c(nrow(meta), length(pb_levels)), dimnames = list(rownames(meta), pb_levels)
)
pb_counts <- counts %*% cell_to_pb
pb_meta <- unique(meta[, c("pb_id", "donor", "TCN2_group")])
rownames(pb_meta) <- pb_meta$pb_id
pb_meta <- pb_meta[colnames(pb_counts), , drop = FALSE]
pb_meta$donor <- factor(pb_meta$donor)

# limma-voom with patient as the duplicateCorrelation block.
design <- model.matrix(~ 0 + TCN2_group, pb_meta)
colnames(design) <- c("TCN2_low", "TCN2_high")
y <- edgeR::DGEList(pb_counts)
y <- y[edgeR::filterByExpr(y, design), , keep.lib.sizes = FALSE]
y <- edgeR::calcNormFactors(y)
v0 <- limma::voom(y, design, plot = FALSE)
cor0 <- limma::duplicateCorrelation(v0, design, block = pb_meta$donor)
v <- limma::voom(y, design, plot = FALSE, block = pb_meta$donor,
                 correlation = cor0$consensus.correlation)
corfit <- limma::duplicateCorrelation(v, design, block = pb_meta$donor)
fit <- limma::lmFit(v, design, block = pb_meta$donor,
                    correlation = corfit$consensus.correlation)
fit <- limma::contrasts.fit(
  fit, limma::makeContrasts(TCN2_high_vs_low = TCN2_high - TCN2_low,
                            levels = design)
)
fit <- limma::eBayes(fit, robust = TRUE)

de <- limma::topTable(fit, coef = "TCN2_high_vs_low", number = Inf, sort.by = "P")
de$gene <- rownames(de)
write.csv(de[, c("gene", setdiff(colnames(de), "gene"))],
          "limma_DE_TCN2-high_vs_low.csv", row.names = FALSE)

# GSEA ranked by the donor-aware moderated t statistic.
rank_stats <- sort(setNames(fit$t[, 1], rownames(fit$t)), decreasing = TRUE)
gsea <- run_fgsea(rank_stats)
gsea <- gsea[order(gsea$padj, -abs(gsea$NES)), ]
gsea_out <- gsea
gsea_out$leadingEdge <- vapply(gsea_out$leadingEdge, paste, collapse = ";",
                               FUN.VALUE = character(1))
write.csv(gsea_out, "GSEA_Reactome_pseudobulk.csv", row.names = FALSE)

# Mirror the same contrast: low-vs-high has the opposite NES but the same FDR.
plot_high <- gsea[gsea$pathway %in% pathway_order,
                  c("pathway", "pval", "padj", "NES", "size")]
plot_high$enriched_in <- "TCN2-high"
plot_low <- plot_high
plot_low$enriched_in <- "TCN2-low"
plot_low$NES <- -plot_low$NES
plot_C <- format_gsea_data(bind_rows(plot_high, plot_low))
write.csv(plot_C, "FigS5C_selected_pathways.csv", row.names = FALSE)

pC <- gsea_dotplot(plot_C)

ggsave("FigS5C-heatmap.pdf", pC, width = 7, height = 4.3,
       bg = "white", device = cairo_pdf)
ggsave("FigS5C-heatmap.png", pC, width = 7, height = 4.3,
       dpi = 300, bg = "white")
print(pC)
