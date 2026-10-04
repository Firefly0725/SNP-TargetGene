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


#--------S5C--------
rds_file <- "tam_group_60samples.RDS"
input_gmt <- "c2.cp.reactome.v2026.1.Hs.symbols.gmt"
output_dir <- Sys.getenv(
  "FIGS5C_OUTPUT_DIR",
  unset = "FigS5C_pseudobulk_results"
)

donor_col <- "sample"
group_col <- "tam_group"
low_label <- "TCN2-low"
high_label <- "TCN2-high"

# Keep 1 to analyze the entire 60-patient cohort. For a stricter sensitivity
# analysis, change this to 5 or 10 and optionally set paired_donors_only=TRUE.
min_cells_per_pseudobulk <- as.integer(Sys.getenv(
  "FIGS5C_MIN_CELLS",
  unset = "1"
))
paired_donors_only <- tolower(Sys.getenv(
  "FIGS5C_PAIRED_ONLY",
  unset = "false"
)) %in% c("true", "t", "1", "yes", "y")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------- input checks --------------------------------
stopifnot(file.exists(rds_file), file.exists(input_gmt))
tam <- readRDS(rds_file)
meta <- tam[[]]

required_cols <- c(donor_col, group_col)
missing_cols <- setdiff(required_cols, colnames(meta))
if (length(missing_cols) > 0L) {
  stop("Missing metadata column(s): ", paste(missing_cols, collapse = ", "))
}

meta$donor <- as.character(meta[[donor_col]])
meta$TCN2_group <- factor(
  as.character(meta[[group_col]]),
  levels = c(low_label, high_label)
)
if (anyNA(meta$donor) || anyNA(meta$TCN2_group)) {
  stop("Donor or TCN2 group contains missing/unrecognized values.")
}

counts <- GetAssayData(tam, assay = "RNA", layer = "counts")
if (!identical(colnames(counts), rownames(meta))) {
  meta <- meta[colnames(counts), , drop = FALSE]
  stopifnot(identical(colnames(counts), rownames(meta)))
}

# ------------------------ patient-level pseudobulking -------------------------
meta$pseudobulk_id <- paste(meta$donor, meta$TCN2_group, sep = "__")
pseudobulk_levels <- unique(meta$pseudobulk_id)

# Sparse one-hot cell-to-pseudobulk matrix; summing raw UMI counts creates one
# expression profile per patient and TCN2 group.
cell_to_pb <- sparseMatrix(
  i = seq_len(nrow(meta)),
  j = match(meta$pseudobulk_id, pseudobulk_levels),
  x = 1,
  dims = c(nrow(meta), length(pseudobulk_levels)),
  dimnames = list(rownames(meta), pseudobulk_levels)
)
pb_counts <- counts %*% cell_to_pb

pb_meta <- unique(meta[, c("pseudobulk_id", "donor", "TCN2_group")])
rownames(pb_meta) <- pb_meta$pseudobulk_id
pb_meta <- pb_meta[colnames(pb_counts), , drop = FALSE]
pb_meta$n_cells <- as.integer(table(meta$pseudobulk_id)[pb_meta$pseudobulk_id])
pb_meta$library_size <- as.numeric(colSums(pb_counts))

keep_pb <- pb_meta$n_cells >= min_cells_per_pseudobulk
if (paired_donors_only) {
  donor_group_n <- table(pb_meta$donor[keep_pb], pb_meta$TCN2_group[keep_pb])
  paired_donors <- rownames(donor_group_n)[rowSums(donor_group_n > 0) == 2L]
  keep_pb <- keep_pb & pb_meta$donor %in% paired_donors
}

pb_counts <- pb_counts[, keep_pb, drop = FALSE]
pb_meta <- pb_meta[keep_pb, , drop = FALSE]
pb_meta$TCN2_group <- droplevels(pb_meta$TCN2_group)
pb_meta$donor <- factor(pb_meta$donor)

if (nlevels(pb_meta$TCN2_group) != 2L) {
  stop("Both TCN2 groups are required after pseudobulk QC.")
}

write.csv(
  pb_meta,
  file.path(output_dir, "pseudobulk_sample_QC.csv"),
  row.names = FALSE
)

# -------------------- limma-voom + duplicateCorrelation ----------------------
dge <- DGEList(counts = pb_counts)
design <- model.matrix(~ 0 + TCN2_group, data = pb_meta)
colnames(design) <- c("TCN2_low", "TCN2_high")

keep_gene <- filterByExpr(dge, design = design)
dge <- dge[keep_gene, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")

# First pass estimates within-patient correlation. The second voom pass uses
# that correlation when estimating observation-level precision weights.
v0 <- voom(dge, design = design, plot = FALSE)
corfit0 <- duplicateCorrelation(
  v0,
  design = design,
  block = pb_meta$donor
)
v <- voom(
  dge,
  design = design,
  plot = FALSE,
  block = pb_meta$donor,
  correlation = corfit0$consensus.correlation
)
corfit <- duplicateCorrelation(
  v,
  design = design,
  block = pb_meta$donor
)

fit <- lmFit(
  v,
  design = design,
  block = pb_meta$donor,
  correlation = corfit$consensus.correlation
)
contrast_matrix <- makeContrasts(
  TCN2_high_vs_low = TCN2_high - TCN2_low,
  levels = design
)
fit <- contrasts.fit(fit, contrast_matrix)
fit <- eBayes(fit, robust = TRUE)

de <- topTable(
  fit,
  coef = "TCN2_high_vs_low",
  number = Inf,
  sort.by = "P"
)
de$gene <- rownames(de)
de <- de[, c("gene", setdiff(colnames(de), "gene"))]
write.csv(
  de,
  file.path(output_dir, "limma_DE_TCN2-high_vs_low.csv"),
  row.names = FALSE
)

# ------------------------------- Reactome GSEA -------------------------------
# Positive t / NES: enriched in TCN2-high; negative t / NES: TCN2-low.
rank_stats <- fit$t[, "TCN2_high_vs_low"]
names(rank_stats) <- rownames(fit$t)
rank_stats <- sort(rank_stats[is.finite(rank_stats)], decreasing = TRUE)

if (anyDuplicated(names(rank_stats))) {
  rank_stats <- tapply(rank_stats, names(rank_stats), function(z) z[which.max(abs(z))])
  rank_stats <- sort(rank_stats, decreasing = TRUE)
}

pathways <- gmtPathways(input_gmt)
gsea <- fgseaMultilevel(
  pathways = pathways,
  stats = rank_stats,
  minSize = 10,
  maxSize = 500,
  eps = 0,
  nPermSimple = 10000,
  nproc = 1
)
gsea <- as.data.frame(gsea)
gsea$enriched_in <- ifelse(gsea$NES >= 0, high_label, low_label)
gsea <- gsea[order(gsea$padj, -abs(gsea$NES)), ]

gsea_export <- gsea
gsea_export$leadingEdge <- vapply(
  gsea_export$leadingEdge,
  paste,
  collapse = ";",
  FUN.VALUE = character(1)
)
write.csv(
  gsea_export,
  file.path(output_dir, "GSEA_Reactome_all_pathways.csv"),
  row.names = FALSE
)
write.csv(
  gsea_export[!is.na(gsea_export$padj) & gsea_export$padj < 0.05, ],
  file.path(output_dir, "GSEA_Reactome_FDR_lt_0.05.csv"),
  row.names = FALSE
)

# --------------------------------- Fig. S5C ----------------------------------
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
  "REACTOME_INTERFERON_GAMMA_SIGNALING" = "IFN-gamma signaling",
  "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING" = "IFN-alpha/beta signaling",
  "REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION" = "Cross-presentation",
  "REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION" = "MHC-II antigen presentation",
  "REACTOME_ANTIGEN_PRESENTATION_FOLDING_ASSEMBLY_AND_PEPTIDE_LOADING_OF_CLASS_I_MHC" = "MHC-I antigen presentation",
  "REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL" = "Lymphoid-myeloid immunoregulation",
  "REACTOME_INTERLEUKIN_10_SIGNALING" = "IL-10 signaling"
)

missing_pathways <- setdiff(pathway_order, gsea$pathway)
if (length(missing_pathways) > 0L) {
  warning("Selected pathway(s) absent after size filtering: ",
          paste(missing_pathways, collapse = ", "))
}

# Mirror display requested for Fig. S5C: the low-vs-high row is the exact
# sign-reversed view of the high-vs-low contrast. P values and FDR are shared
# because these are two orientations of the same statistical test.
plot_df_high <- gsea[gsea$pathway %in% pathway_order, c(
  "pathway", "enriched_in", "pval", "padj", "NES", "size"
)]
plot_df_high$enriched_in <- high_label
plot_df_low <- plot_df_high
plot_df_low$enriched_in <- low_label
plot_df_low$NES <- -plot_df_low$NES
plot_df <- rbind(plot_df_high, plot_df_low)
plot_df$log10FDR <- -log10(pmax(plot_df$padj, 1e-300))
plot_df$pathway <- factor(plot_df$pathway, levels = pathway_order)
plot_df$enriched_in <- factor(
  plot_df$enriched_in,
  levels = c(low_label, high_label)
)
write.csv(
  plot_df,
  file.path(output_dir, "FigS5C_selected_pathways.csv"),
  row.names = FALSE
)

p1 <- ggplot(plot_df, aes(x = pathway, y = enriched_in)) +
  geom_point(
    aes(size = log10FDR, fill = NES),
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
    name = expression(-log[10](FDR)),
    range = c(3, 10),
    breaks = pretty_breaks(4)
  ) +
  scale_x_discrete(
    position = "top",
    labels = pathway_labels,
    drop = FALSE
  ) +
  scale_y_discrete(position = "right", drop = FALSE) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 12) +
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      size = 12, angle = 55, colour = "black", hjust = 0, vjust = 0
    ),
    axis.text.y = element_text(size = 13, colour = "black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    legend.position = "bottom",
    legend.box = "horizontal",
    plot.margin = margin(t = 5, r = 5, b = 2, l = 5, unit = "pt")
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "top", title.hjust = 0.5, order = 1,
      barwidth = grid::unit(4, "cm")
    ),
    size = guide_legend(
      title.position = "top", title.hjust = 0.5, order = 2
    )
  )

ggsave(
  file.path(output_dir, "FigS5C-pseudobulk-GSEA.pdf"),
  p1, width = 7, height = 4.3, bg = "white"
)
ggsave(
  file.path(output_dir, "FigS5C-pseudobulk-GSEA.png"),
  p1, width = 7, height = 4.3, dpi = 600, bg = "white"
)

# ---------------------------- reproducibility log -----------------------------
n_donors <- nlevels(pb_meta$donor)
n_paired <- sum(rowSums(table(pb_meta$donor, pb_meta$TCN2_group) > 0) == 2L)
summary_lines <- c(
  paste0("Cells in Seurat object: ", ncol(tam)),
  paste0("Pseudobulks analyzed: ", ncol(pb_counts)),
  paste0("Donors analyzed: ", n_donors),
  paste0("Donors represented in both groups: ", n_paired),
  paste0("Genes after filterByExpr: ", nrow(dge)),
  paste0("Consensus within-donor correlation: ",
         signif(corfit$consensus.correlation, 5)),
  paste0("Reactome pathways tested: ", nrow(gsea)),
  paste0("Reactome pathways at FDR < 0.05: ",
         sum(!is.na(gsea$padj) & gsea$padj < 0.05))
)
writeLines(summary_lines, file.path(output_dir, "analysis_summary.txt"))
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))

message(paste(summary_lines, collapse = "\n"))
print(p1)

