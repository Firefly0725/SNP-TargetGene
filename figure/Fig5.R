#!/usr/bin/env Rscript
# Figure 5 revision: paired differential expression and TCGA survival analysis.

#--------SETUP AND INPUT PATHS--------
project_dir <- "/Users/L/Desktop/ESCC/figure5_revision"
data_dir <- file.path(project_dir, "data")
results_dir <- file.path(project_dir, "results")
derived_dir <- file.path(project_dir, "derived")
qa_dir <- file.path(project_dir, "qa")
for (directory in c(results_dir, derived_dir, qa_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

use_original_96_samples <- TRUE
run_supplementary_heatmap <- FALSE
suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(survival)
  library(survminer)
  library(ggtext)
})
write_tsv <- function(x, filename, directory = derived_dir) {
  write.table(x, file.path(directory, filename), sep = "\t",
              quote = FALSE, row.names = FALSE, na = "NA")
}
read_clinical <- function(filename) {
  read.delim(file.path(data_dir, filename), comment.char = "#", check.names = FALSE,
             na.strings = c("NA", "", "[Not Available]", "[Not Applicable]",
                            "[Unknown]", "unknown"), stringsAsFactors = FALSE)
}
pdf_device <- function(filename, width, height, bg = "white", ...) {
  if (capabilities("aqua")) {
    grDevices::quartz(type = "pdf", file = filename, width = width,
                     height = height, bg = bg, family = "Arial")
  } else {
    grDevices::pdf(file = filename, width = width, height = height, bg = bg)
  }
}

#--------5A--------
gse_file <- file.path(data_dir, "GSE53625_clean.RData")

de_adjP_cutoff <- 0.05
de_logFC_cutoff <- 0.5
plot_base_size <- 14

loaded_names <- load(gse_file)
gse <- GSE53625_clean

expr_raw <- as.matrix(gse$expr)          
clinical <- as.data.frame(gse$clinical, stringsAsFactors = FALSE)

common <- intersect(colnames(expr_raw), rownames(clinical))
expr <- expr_raw[, common, drop = FALSE]
clinical <- clinical[common, , drop = FALSE]

group_raw <- trimws(tolower(clinical$group))
group <- ifelse(group_raw %in% c("tumor","tumour","cancer","escc","case"), "cancer",
                ifelse(group_raw %in% c("normal","control","adjacent normal","adjacent_normal"), "normal", NA))

n <- ncol(expr)
cancer_pos <- seq(1, n, by = 2)
normal_pos <- seq(2, n, by = 2)

pair_id <- factor(rep(1:(n/2), each = 2))
tissue <- factor(group, levels = c("normal", "cancer"))

# Preserve the original adjacent tumor/normal pairing, but verify its ordering.
stopifnot(n %% 2 == 0, !anyNA(tissue),
          all(tissue[cancer_pos] == "cancer"), all(tissue[normal_pos] == "normal"))
pair_info <- data.frame(sample_id = colnames(expr), pair_id = pair_id, tissue = tissue)

design <- model.matrix(~ pair_id + tissue, data = pair_info)
fit <- lmFit(expr, design)
fit <- eBayes(fit)

coef_name <- "tissuecancer"
de_all <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P", adjust.method = "BH")
de_all$gene <- rownames(de_all)
rownames(de_all) <- NULL

de_sig <- de_all[de_all$adj.P.Val < de_adjP_cutoff & abs(de_all$logFC) >= de_logFC_cutoff, ]


de_all$neglog10FDR <- -log10(pmax(de_all$adj.P.Val, 1e-300))
de_all$DE_status <- "NS"
de_all$DE_status[de_all$adj.P.Val < de_adjP_cutoff & de_all$logFC >  de_logFC_cutoff] <- "Up"
de_all$DE_status[de_all$adj.P.Val < de_adjP_cutoff & de_all$logFC < -de_logFC_cutoff] <- "Down"


load(file.path(data_dir, "filtered_biomarker.Rdata"))
genes_to_label <- biomarker
de_all_sig <- de_all %>% filter(DE_status %in% c("Up","Down"))

label_data <- de_all_sig[de_all_sig$gene %in% genes_to_label, ]

volcano_plot <- ggplot(de_all, aes(x = logFC, y = neglog10FDR)) +
  geom_point(aes(color = DE_status), size = 1.2, alpha = 0.75) +
  geom_vline(xintercept = c(-de_logFC_cutoff, de_logFC_cutoff), linetype = 2, color = "grey60") +
  geom_hline(yintercept = -log10(de_adjP_cutoff), linetype = 2, color = "grey60") +
  scale_color_manual(values = c(NS = "grey70", Up = "#D55E00", Down = "#0072B2")) +
  geom_point(data = label_data, 
             aes(x = logFC, y = neglog10FDR),
             color = "black", fill = "gold", shape = 21, size = 3, stroke = 0.8) +
  geom_text_repel(data = label_data,
                  aes(x = logFC, y = neglog10FDR, label = gene),
                  size = 14 / .pt,   
                  box.padding = 0.35, point.padding = 0.2,
                  segment.color = "black", segment.size = 0.5) +
  labs(x = "logFC (Tumor - Normal)", y = expression(-log[10]("FDR")), color = NULL) +
  theme_classic(base_size = plot_base_size) +  
  theme(plot.title = element_blank(), 
        panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 14),   
        axis.title = element_text(color = "black", size = 14))
set.seed(2026)  # Reproducible label placement; statistical analysis is unchanged.
volcano_plot
ggsave(file.path(results_dir, "Fig5A-paired_DE_volcano.png"), volcano_plot, width = 5.5, height = 5, dpi = 300)
ggsave(file.path(results_dir, "Fig5A-paired_DE_volcano.pdf"), volcano_plot, width = 5.5, height = 5, bg = "white")



write_tsv(de_all, "Fig5A_all_gene_differential_expression.tsv")
write_tsv(pair_info, "Fig5A_pair_assignments.tsv")

if (run_supplementary_heatmap) {
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
#--------S1A--------
font_size <- 14
row_var <- apply(expr, 1, var, na.rm = TRUE)
top_genes <- names(sort(row_var, decreasing = TRUE))[1:min(30, nrow(expr))]
mat <- expr[top_genes, , drop = FALSE]
mat_z <- t(scale(t(mat)))
mat_z[!is.finite(mat_z)] <- 0

tissue_colors <- c("cancer" = "#D55E00", "normal" = "#0072B2")
pair_levels <- levels(pair_info$pair_id)
n_pair <- length(pair_levels)
pair_colors <- setNames(rainbow(n_pair), pair_levels)  # Original pair-color palette.

annotation_col <- data.frame(
  Tissue = pair_info$tissue,
  Pair = pair_info$pair_id
)
rownames(annotation_col) <- pair_info$sample_id

ha <- HeatmapAnnotation(
  df = annotation_col,
  col = list(
    Tissue = tissue_colors,
    Pair = pair_colors
  ),
  show_legend = c(Tissue = FALSE, Pair = FALSE)
)

all_genes <- rownames(mat_z)
row_labels <- ifelse(all_genes %in% label_data$gene, all_genes, "")

ht <- Heatmap(
  mat_z,
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  top_annotation = ha,
  show_column_names = FALSE,
  row_labels = row_labels,
  row_names_gp = gpar(fontsize = font_size),
  column_title_gp = gpar(fontsize = font_size),
  heatmap_legend_param = list(
    title = "Z-score",
    at = c(-2, 0, 2),
    title_gp = gpar(fontsize = font_size),
    labels_gp = gpar(fontsize = font_size)
  )
)

tissue_legend <- Legend(
  title = "Tissue",
  at = c("cancer", "normal"),
  legend_gp = gpar(fill = tissue_colors),
  labels = c("Cancer", "Normal"),
  title_gp = gpar(fontsize = font_size),
  labels_gp = gpar(fontsize = font_size)
)

draw(ht,
     annotation_legend_list = list(tissue_legend),
     annotation_legend_side = "right",
     merge_legend = FALSE)


png(file.path(results_dir, "SFig1A-heatmap_with_tissue_legend.png"), width = 2000, height = 1600, res = 300)
draw(ht,
     annotation_legend_list = list(tissue_legend),
     annotation_legend_side = "right",
     merge_legend = FALSE)
dev.off()

pdf(file.path(results_dir, "SFig1A-heatmap_with_tissue_legend.pdf"), width = 7, height = 5)
draw(ht,
     annotation_legend_list = list(tissue_legend),
     annotation_legend_side = "right",
     merge_legend = FALSE)
dev.off()




}

#--------TCGA DATA PREPARATION--------
genes <- c("GBA", "CASP10", "PPP1R14A", "MUC1", "EFNA1", "TCN2",
           "ZNF7", "CLDN4", "WBSCR27", "TAPBPL", "TTC28")
patients <- read_clinical("TCGA_legacy_clinical_patient.txt")
samples <- read_clinical("TCGA_legacy_clinical_sample.txt")
rsem <- read_clinical("TCGA_legacy_mrna_seq_v2_rsem.txt")
stopifnot(all(vapply(genes, function(g) sum(rsem$Hugo_Symbol == g, na.rm = TRUE), integer(1)) == 1))
clinical_tcga <- merge(samples, patients, by = "PATIENT_ID")
clinical_tcga <- clinical_tcga[
  !is.na(clinical_tcga$ONCOTREE_CODE) & clinical_tcga$ONCOTREE_CODE == "ESCC" &
    clinical_tcga$SAMPLE_ID %in% names(rsem), ]
write_tsv(clinical_tcga, "TCGA_all_ESCC_sample_audit.tsv")
clinical_tcga <- clinical_tcga[clinical_tcga$SAMPLE_TYPE == "Primary", ]
clinical_tcga <- clinical_tcga[order(clinical_tcga$PATIENT_ID, clinical_tcga$SAMPLE_ID), ]
stopifnot(!anyDuplicated(clinical_tcga$PATIENT_ID))
stage <- toupper(sub("^Stage ", "", clinical_tcga$AJCC_PATHOLOGIC_TUMOR_STAGE))
stage <- sub("[ABC]$", "", stage)
tcga <- data.frame(
  sample_id = clinical_tcga$SAMPLE_ID, patient_id = clinical_tcga$PATIENT_ID,
  OS_months = as.numeric(clinical_tcga$OS_MONTHS),
  OS_event = match(clinical_tcga$OS_STATUS, c("0:LIVING", "1:DECEASED")) - 1L,
  age10 = as.numeric(clinical_tcga$AGE) / 10,
  sex = factor(tolower(clinical_tcga$SEX), levels = c("female", "male")),
  stage = factor(stage, levels = c("I", "II", "III", "IV")),
  # GX is unassessable grade; discrepant smoking codes are also missing.
  grade = factor(clinical_tcga$GRADE, levels = c("G1", "G2", "G3")),
  smoking = factor(clinical_tcga$TOBACCO_SMOKING_HISTORY_INDICATOR,
                   levels = c("1", "2", "3", "4")),
  alcohol = factor(tolower(clinical_tcga$ALCOHOL_HISTORY_DOCUMENTED), levels = c("no", "yes"))
)
tcga <- tcga[is.finite(tcga$OS_months) & tcga$OS_months > 0 & !is.na(tcga$OS_event), ]
stopifnot(nrow(tcga) == 95, sum(tcga$OS_event) == 32)
expression <- log2(as.matrix(rsem[match(genes, rsem$Hugo_Symbol), tcga$sample_id]) + 1)
rownames(expression) <- genes
stopifnot(all(is.finite(expression)))
# Standardization uses the same full 95-patient cohort for every adjusted model.
for (gene in genes) {
  tcga[[gene]] <- as.numeric(expression[gene, ])
  tcga[[paste0(gene, "_z")]] <- as.numeric(scale(tcga[[gene]]))
}
tcga$OS_days <- tcga$OS_months * 30.4
write_tsv(tcga, "TCGA_95_patient_model_inputs.tsv")

#--------5B--------

fits <- model_rows <- ph_rows <- inclusion_rows <- list()
for (gene in genes) {
  gene_term <- paste0(gene, "_z")
  formula <- as.formula(paste("Surv(OS_months, OS_event) ~", gene_term,
                              "+ age10 + sex + stage + grade + smoking + alcohol"))
  variables <- all.vars(formula)
  included <- complete.cases(tcga[, variables])
  analysis_data <- droplevels(tcga[included, ])
  stopifnot(nrow(analysis_data) == 79, sum(analysis_data$OS_event) == 26)
  fit_warnings <- character()
  fit <- withCallingHandlers(
    coxph(formula, data = analysis_data, ties = "efron", x = TRUE, y = TRUE, model = TRUE),
    warning = function(w) {
      fit_warnings <<- c(fit_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  fits[[gene]] <- fit
  sm <- summary(fit)
  ph <- cox.zph(fit)
  model_rows[[gene]] <- data.frame(
    gene, term = rownames(sm$coefficients), HR = sm$conf.int[, "exp(coef)"],
    CI_low = sm$conf.int[, "lower .95"], CI_high = sm$conf.int[, "upper .95"],
    p = sm$coefficients[, "Pr(>|z|)"], n = fit$n, events = fit$nevent,
    parameters = length(coef(fit)), events_per_parameter = fit$nevent / length(coef(fit)),
    warning = paste(unique(fit_warnings), collapse = " | ")
  )
  ph_rows[[gene]] <- data.frame(gene, term = rownames(ph$table), ph$table)
  inclusion_rows[[gene]] <- data.frame(
    gene, sample_id = tcga$sample_id, patient_id = tcga$patient_id, included,
    missing_variables = apply(is.na(tcga[, variables]), 1, function(x) paste(variables[x], collapse = ";"))
  )
}
all_terms <- do.call(rbind, model_rows)
forest_results <- all_terms[all_terms$term == paste0(all_terms$gene, "_z"), ]
forest_results <- forest_results[match(genes, forest_results$gene), ]
forest_results$FDR_BH_11 <- p.adjust(forest_results$p, method = "BH")
write_tsv(all_terms, "Fig5B_all_model_coefficients.tsv")
write_tsv(forest_results, "Fig5B_11genes_multivariable_cox.tsv")
write_tsv(do.call(rbind, ph_rows), "Fig5B_proportional_hazards_tests.tsv")
write_tsv(do.call(rbind, inclusion_rows), "Fig5B_patient_inclusion.tsv")
saveRDS(fits, file.path(derived_dir, "Fig5B_multivariable_models.rds"))

# Plot from the exported table. Preserve the original order and highlight only TCN2.
plot_df <- read.delim(file.path(derived_dir, "Fig5B_11genes_multivariable_cox.tsv"))
plot_df$gene <- factor(plot_df$gene, levels = rev(genes))
pB <- ggplot(plot_df, aes(x = HR, y = gene)) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey60", linewidth = 0.8) +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high), orientation = "y",
                width = 0.20, color = "grey9", linewidth = 0.8) +
  geom_point(size = 3, color = "#D55E00") + scale_x_log10() +
  labs(x = "Hazard ratio (per 1 SD increase, log scale)", y = NULL) +
  scale_y_discrete(labels = function(x) ifelse(x == "TCN2", paste0("**", x, "**"), x)) +
  theme_classic(base_size = 15) +
  theme(plot.title = element_blank(), panel.grid = element_blank(),
        axis.title = element_text(size = 15, face = "plain"),
        axis.text.y = ggtext::element_markdown(size = 15, face = "plain", color = "black"),
        axis.text.x = element_text(size = 15, face = "plain", color = "black"))
ggsave(file.path(results_dir, "Fig5B_multivariable_forest.png"), pB,
       width = 7, height = 8, dpi = 300, bg = "white")
ggsave(file.path(results_dir, "Fig5B_multivariable_forest.pdf"), pB,
       width = 7, height = 8, bg = "white", device = pdf_device)

#--------SHARED SURVIVAL PREDICTION AND PLOTTING FUNCTIONS--------
# Low/High groups are for display and observed risk counts only. The Cox term is
# continuous TCN2_z. Each curve fixes expression at its group median and averages
# predictions over the same model patients (marginal standardization).
# The p-value is the continuous-expression Wald test, not a log-rank test.
plot_cox_panel <- function(fit, analysis_data, reference_data, heading, prefix) {
  cutoff <- median(reference_data$TCN2)
  reference_data$group <- factor(ifelse(reference_data$TCN2 > cutoff, "High", "Low"),
                                 levels = c("Low", "High"))
  analysis_data$group <- factor(ifelse(analysis_data$TCN2 > cutoff, "High", "Low"),
                                levels = c("Low", "High"))
  analysis_data$OS_days <- analysis_data$OS_months * 30.4
  representatives <- aggregate(cbind(TCN2, TCN2_z) ~ group, reference_data, median)
  representatives$cutoff <- cutoff
  write_tsv(representatives, paste0(prefix, "_representative_expression.tsv"))
  baseline <- basehaz(fit, centered = FALSE)
  curve_rows <- validation <- list()
  for (group in c("Low", "High")) {
    newdata <- analysis_data
    newdata$TCN2_z <- representatives$TCN2_z[representatives$group == group]
    lp <- predict(fit, newdata = newdata, type = "lp", reference = "zero")
    survival <- rowMeans(exp(-outer(baseline$hazard, exp(lp))))
    independent <- survfit(fit, newdata = newdata, se.fit = FALSE)
    expected <- if (is.matrix(independent$surv)) rowMeans(independent$surv) else independent$surv
    index <- match(baseline$time, independent$time)
    stopifnot(!anyNA(index))
    error <- max(abs(survival - expected[index]))
    stopifnot(error < 1e-10, all(survival >= 0 & survival <= 1), all(diff(survival) <= 1e-12))
    validation[[group]] <- data.frame(group, max_abs_difference_vs_survfit = error)
    # Match the previous figures: no curve beyond its display group's follow-up.
    keep <- baseline$time <= max(analysis_data$OS_months[analysis_data$group == group])
    curve_rows[[group]] <- data.frame(
      time = c(0, baseline$time[keep] * 30.4), survival = c(1, survival[keep]), group
    )
  }
  curves <- do.call(rbind, curve_rows)
  curves$group <- factor(curves$group, levels = c("Low", "High"))
  write_tsv(curves, paste0(prefix, "_predicted_survival.tsv"))
  write_tsv(do.call(rbind, validation), paste0(prefix, "_curve_validation.tsv"), qa_dir)
  sm <- summary(fit)
  write_tsv(data.frame(
    term = "TCN2_z", HR = sm$conf.int["TCN2_z", "exp(coef)"],
    CI_low = sm$conf.int["TCN2_z", "lower .95"], CI_high = sm$conf.int["TCN2_z", "upper .95"],
    p = sm$coefficients["TCN2_z", "Pr(>|z|)"], n = fit$n,
    unique_patients = length(unique(analysis_data$patient_id)), events = fit$nevent
  ), paste0(prefix, "_model_summary.tsv"))
  p_label <- sprintf("p=%.3f", sm$coefficients["TCN2_z", "Pr(>|z|)"])

  # Original Fig5 D/E styling: 15-point text, blue/orange curves, risk table.
  my_theme <- theme_bw() +
    theme(panel.grid = element_blank(), axis.title = element_text(size = 15),
          axis.text = element_text(size = 15), legend.title = element_text(size = 15),
          legend.text = element_text(size = 15), plot.title = element_text(size = 15))
  my_table_theme <- theme_cleantable() +
    theme(axis.text = element_text(size = 15), text = element_text(size = 15))
  observed <- survfit(Surv(OS_days, OS_event) ~ group, data = analysis_data)
  p <- ggsurvplot(
    observed, data = analysis_data, risk.table = TRUE, risk.table.col = "strata",
    pval = FALSE, pval.method = FALSE, conf.int = FALSE,
    palette = c("#0072B2", "#D55E00"), xlab = "Overall survival time",
    ylab = "Survival probability", legend.title = heading, legend.labs = c("Low", "High"),
    risk.table.y.text.col = TRUE, risk.table.y.text = FALSE,
    tables.theme = my_table_theme, ggtheme = my_theme,
    font.x = 15, font.y = 15, font.tickslab = 15, font.legend = 15,
    font.pval = 15, font.risk.table = 15
  )
  time_range <- range(c(0, analysis_data$OS_days))
  breaks <- pretty(time_range)
  censor <- analysis_data[analysis_data$OS_event == 0, c("OS_days", "group")]
  names(censor) <- c("time", "group")
  censor$survival <- NA_real_
  for (group in c("Low", "High")) {
    cc <- curves[curves$group == group, ]
    index <- which(censor$group == group)
    step <- findInterval(censor$time[index], cc$time)
    censor$survival[index] <- cc$survival[pmax(1, step)]
  }
  # Replace the observed KM curve with predictions from the requested Cox model.
  p$plot <- ggplot(curves, aes(x = time, y = survival, color = group)) +
    geom_step(linewidth = 1) +
    geom_point(data = censor, shape = 3, size = 4.5, show.legend = FALSE) +
    scale_color_manual(values = c(Low = "#0072B2", High = "#D55E00"), name = heading) +
    scale_x_continuous(breaks = breaks, limits = time_range, expand = expansion(mult = 0.05)) +
    scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1),
                       labels = function(x) sprintf("%.2f", x), expand = expansion(mult = 0.05)) +
    labs(x = "Overall survival time", y = "Survival probability") +
    annotate("text", x = time_range[2] * 0.02, y = 0.2, label = p_label, hjust = 0, size = 5) +
    guides(color = guide_legend(title.position = "top", title.hjust = 0.5)) +
    my_theme + theme(legend.position = "top")
  p$table <- p$table + xlab("Time") +
    scale_x_continuous(breaks = breaks, limits = time_range) +
    theme(axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),
          legend.position = "none")
  risk <- summary(observed, times = breaks[breaks >= 0 & breaks <= max(analysis_data$OS_days)], extend = TRUE)
  risk_table <- data.frame(group = as.character(risk$strata), time_days = risk$time,
                            n_risk = risk$n.risk, n_event = risk$n.event, n_censor = risk$n.censor)
  stopifnot(sum(risk_table$n_risk[risk_table$time_days == 0]) == fit$n)
  write_tsv(risk_table, paste0(prefix, "_observed_risk_table.tsv"))
  curve_grob <- ggplotGrob(p$plot)
  risk_grob <- ggplotGrob(p$table)
  stopifnot(length(curve_grob$widths) == length(risk_grob$widths))
  widths <- grid::unit.pmax(curve_grob$widths, risk_grob$widths)
  curve_grob$widths <- risk_grob$widths <- widths
  combined <- gridExtra::arrangeGrob(curve_grob, risk_grob, ncol = 1, heights = c(0.75, 0.25))
  # Original PNG: 1600 x 1600 at 300 dpi. Original PDF: 6 x 6 inches.
  png(file.path(results_dir, paste0(prefix, ".png")), width = 1600, height = 1600, res = 300)
  grid::grid.newpage()
  grid::grid.draw(combined)
  dev.off()
  pdf(file.path(results_dir, paste0(prefix, ".pdf")), width = 6, height = 6)
  grid::grid.newpage()
  grid::grid.draw(combined)
  dev.off()
  invisible(p)
}

#--------5C--------

legacy_env <- new.env()
load(file.path(data_dir, "ESCC_two_tcga_extracted_clean.RData"), envir = legacy_env)
legacy <- legacy_env$ESCC_tcga_legacy
colnames(legacy$expr) <- gsub("\\.", "-", colnames(legacy$expr))
legacy_clinical <- legacy$clinical
uni_data <- data.frame(
  sample_id = legacy_clinical$SAMPLE_ID, patient_id = legacy_clinical$PATIENT_ID,
  OS_months = as.numeric(legacy_clinical$OS_MONTHS),
  OS_event = match(legacy_clinical$OS_STATUS, c("0:LIVING", "1:DECEASED")) - 1L,
  TCN2 = log2(as.numeric(legacy$expr["TCN2", legacy_clinical$SAMPLE_ID]) + 1)
)
uni_data$duplicate_patient <- duplicated(uni_data$patient_id) |
  duplicated(uni_data$patient_id, fromLast = TRUE)
if (!use_original_96_samples) {
  uni_data <- uni_data[substr(uni_data$sample_id, 14, 15) == "01", ]
  stopifnot(!anyDuplicated(uni_data$patient_id))
}
stopifnot(all(is.finite(uni_data$TCN2)), all(uni_data$OS_months > 0), !anyNA(uni_data$OS_event))
uni_data$TCN2_z <- as.numeric(scale(uni_data$TCN2))
write_tsv(uni_data, "Fig5C_univariable_input_audit.tsv")
fit_uni <- coxph(Surv(OS_months, OS_event) ~ TCN2_z, data = uni_data,
                 ties = "efron", x = TRUE, y = TRUE, model = TRUE)
if (use_original_96_samples) stopifnot(fit_uni$n == 96, fit_uni$nevent == 33)
pC <- plot_cox_panel(fit_uni, uni_data, uni_data,
                     "Univariable Cox: TCGA legacy ESCC", "Fig5C_univariable_survival")

#--------5D--------
fit_multi <- fits[["TCN2"]]
multi_complete <- complete.cases(tcga[, all.vars(formula(fit_multi))])
multi_data <- droplevels(tcga[multi_complete, ])
pD <- plot_cox_panel(fit_multi, multi_data, tcga,
                     "Multivariable Cox: TCGA legacy ESCC", "Fig5D_multivariable_survival")
saveRDS(list(univariable = fit_uni, multivariable = fit_multi),
        file.path(derived_dir, "Fig5CD_TCN2_models.rds"))
saveRDS(list(A = volcano_plot, B = pB, C = pC, D = pD),
        file.path(derived_dir, "Fig5_plot_objects.rds"))
capture.output(sessionInfo(), file = file.path(qa_dir, "sessionInfo.txt"))
message("Figure 5 complete. Outputs: ", results_dir)

#--------PREVIOUS FIGURE 5: COMMENTED ARCHIVE (NOT EXECUTED)--------
# Previous Figure 5 workflow, archived as comments only.
# library(limma)
# library(ggplot2)
# library(pheatmap)
# library(ggrepel)
# library(ComplexHeatmap)
# library(circlize)
# library(survival)
# library(ggtext)
# library(survminer)
# 
# rm(list = ls()); gc()
# ORIGINAL_DIR <- ""
# setwd(ORIGINAL_DIR)
# 
# #--------5A--------
# data_dir <- "E:/ESCC_SNP"
# gse_file <- file.path(data_dir, "GSE53625_clean.RData")
# 
# de_adjP_cutoff <- 0.05
# de_logFC_cutoff <- 0.5
# plot_base_size <- 14
# 
# loaded_names <- load(gse_file)
# gse <- GSE53625_clean
# 
# expr_raw <- as.matrix(gse$expr)          
# clinical <- as.data.frame(gse$clinical, stringsAsFactors = FALSE)
# 
# common <- intersect(colnames(expr_raw), rownames(clinical))
# expr <- expr_raw[, common, drop = FALSE]
# clinical <- clinical[common, , drop = FALSE]
# 
# group_raw <- trimws(tolower(clinical$group))
# group <- ifelse(group_raw %in% c("tumor","tumour","cancer","escc","case"), "cancer",
#                 ifelse(group_raw %in% c("normal","control","adjacent normal","adjacent_normal"), "normal", NA))
# 
# n <- ncol(expr)
# cancer_pos <- seq(1, n, by = 2)
# normal_pos <- seq(2, n, by = 2)
# 
# pair_id <- factor(rep(1:(n/2), each = 2))
# tissue <- factor(group, levels = c("normal", "cancer"))
# 
# pair_info <- data.frame(sample_id = colnames(expr), pair_id = pair_id, tissue = tissue)
# 
# design <- model.matrix(~ pair_id + tissue, data = pair_info)
# fit <- lmFit(expr, design)
# fit <- eBayes(fit)
# 
# coef_name <- "tissuecancer"
# de_all <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P", adjust.method = "BH")
# de_all$gene <- rownames(de_all)
# rownames(de_all) <- NULL
# 
# de_sig <- de_all[de_all$adj.P.Val < de_adjP_cutoff & abs(de_all$logFC) >= de_logFC_cutoff, ]
# 
# 
# de_all$neglog10FDR <- -log10(pmax(de_all$adj.P.Val, 1e-300))
# de_all$DE_status <- "NS"
# de_all$DE_status[de_all$adj.P.Val < de_adjP_cutoff & de_all$logFC >  de_logFC_cutoff] <- "Up"
# de_all$DE_status[de_all$adj.P.Val < de_adjP_cutoff & de_all$logFC < -de_logFC_cutoff] <- "Down"
# 
# 
# load("filtered_biomarker.Rdata")
# genes_to_label <- biomarker
# de_all_sig <- de_all %>% filter(DE_status %in% c("Up","Down"))
# 
# label_data <- de_all_sig[de_all_sig$gene %in% genes_to_label, ]
# 
# volcano_plot <- ggplot(de_all, aes(x = logFC, y = neglog10FDR)) +
#   geom_point(aes(color = DE_status), size = 1.2, alpha = 0.75) +
#   geom_vline(xintercept = c(-de_logFC_cutoff, de_logFC_cutoff), linetype = 2, color = "grey60") +
#   geom_hline(yintercept = -log10(de_adjP_cutoff), linetype = 2, color = "grey60") +
#   scale_color_manual(values = c(NS = "grey70", Up = "#D55E00", Down = "#0072B2")) +
#   geom_point(data = label_data, 
#              aes(x = logFC, y = neglog10FDR),
#              color = "black", fill = "gold", shape = 21, size = 3, stroke = 0.8) +
#   geom_text_repel(data = label_data,
#                   aes(x = logFC, y = neglog10FDR, label = gene),
#                   size = 14 / .pt,   
#                   box.padding = 0.35, point.padding = 0.2,
#                   segment.color = "black", segment.size = 0.5) +
#   labs(x = "logFC (Tumor - Normal)", y = expression(-log[10]("FDR")), color = NULL) +
#   theme_classic(base_size = plot_base_size) +  
#   theme(plot.title = element_blank(), 
#         panel.grid = element_blank(),
#         axis.text = element_text(color = "black", size = 14),   
#         axis.title = element_text(color = "black", size = 14))
# volcano_plot
# ggsave(file.path("Fig5A-paired_DE_volcano.png"), volcano_plot, width = 5.5, height = 5, dpi = 300)
# ggsave(file.path("Fig5A-paired_DE_volcano.pdf"), volcano_plot, width = 5.5, height = 5, bg = "white")
# 
# 
# #--------S1A--------
# font_size <- 14
# row_var <- apply(expr, 1, var, na.rm = TRUE)
# top_genes <- names(sort(row_var, decreasing = TRUE))[1:min(30, nrow(expr))]
# mat <- expr[top_genes, , drop = FALSE]
# mat_z <- t(scale(t(mat)))
# mat_z[!is.finite(mat_z)] <- 0
# 
# tissue_colors <- c("cancer" = "#D55E00", "normal" = "#0072B2")
# pair_levels <- levels(pair_info$pair_id)
# n_pair <- length(pair_levels)
# pair_colors <- setNames(rainbow(n_pair), pair_levels)  # 或使用 viridis
# 
# annotation_col <- data.frame(
#   Tissue = pair_info$tissue,
#   Pair = pair_info$pair_id
# )
# rownames(annotation_col) <- pair_info$sample_id
# 
# ha <- HeatmapAnnotation(
#   df = annotation_col,
#   col = list(
#     Tissue = tissue_colors,
#     Pair = pair_colors
#   ),
#   show_legend = c(Tissue = FALSE, Pair = FALSE)
# )
# 
# all_genes <- rownames(mat_z)
# row_labels <- ifelse(all_genes %in% label_data$gene, all_genes, "")
# 
# ht <- Heatmap(
#   mat_z,
#   name = "Z-score",
#   col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
#   top_annotation = ha,
#   show_column_names = FALSE,
#   row_labels = row_labels,
#   row_names_gp = gpar(fontsize = font_size),
#   column_title_gp = gpar(fontsize = font_size),
#   heatmap_legend_param = list(
#     title = "Z-score",
#     at = c(-2, 0, 2),
#     title_gp = gpar(fontsize = font_size),
#     labels_gp = gpar(fontsize = font_size)
#   )
# )
# 
# tissue_legend <- Legend(
#   title = "Tissue",
#   at = c("cancer", "normal"),
#   legend_gp = gpar(fill = tissue_colors),
#   labels = c("Cancer", "Normal"),
#   title_gp = gpar(fontsize = font_size),
#   labels_gp = gpar(fontsize = font_size)
# )
# 
# draw(ht,
#      annotation_legend_list = list(tissue_legend),
#      annotation_legend_side = "right",
#      merge_legend = FALSE)
# 
# 
# png("SFig1A-heatmap_with_tissue_legend.png", width = 2000, height = 1600, res = 300)
# draw(ht,
#      annotation_legend_list = list(tissue_legend),
#      annotation_legend_side = "right",
#      merge_legend = FALSE)
# dev.off()
# 
# pdf("SFig1A-heatmap_with_tissue_legend.pdf", width = 7, height = 5)
# draw(ht,
#      annotation_legend_list = list(tissue_legend),
#      annotation_legend_side = "right",
#      merge_legend = FALSE)
# dev.off()
# 
# 
# 
# #--------5B--------
# data_dir <- "E:/ESCC_SNP"
# discovery_file   <- file.path("GSE53625_clean_cancer.RData")
# replication_file <- file.path("ESCC_two_tcga_extracted_clean.RData")
# replication_object_name <- "ESCC_tcga_legacy"
# biomarker_file <- file.path(data_dir, "filtered_gene.Rdata")
# 
# discovery_p_cutoff <- 0.05
# use_discovery_fdr  <- TRUE
# discovery_fdr_cutoff <- 0.10
# require_ph_ok <- TRUE
# ph_p_cutoff <- 0.05
# forest_top_n <- 25
# plot_base_size <- 15
# 
# load(file = "Fig5_dis_and_val_expr_clin.Rdata")
# load(file = "discovery_univariate_cox_all.Rdata")
# 
# plot_df <- head(discovery_res, forest_top_n)
# plot_df$gene <- factor(plot_df$gene, levels = rev(plot_df$gene))
#   
# theme_pub <- theme_classic(base_size = plot_base_size) +
#     theme(plot.title = element_blank(),
#           panel.grid = element_blank(),
#           axis.title = element_text(size = plot_base_size, face = "plain"),
#           axis.text  = element_text(size = plot_base_size, face = "plain"),
#           legend.title = element_text(size = plot_base_size, face = "plain"),
#           legend.text  = element_text(size = plot_base_size, face = "plain"))
#   
# p <- ggplot(plot_df, aes(x = HR, y = gene)) +
#   geom_vline(xintercept = 1, linetype = 2, color = "grey60", linewidth = 0.8) +     
#   geom_errorbarh(aes(xmin = HR_L, xmax = HR_H), height = 0.20, color = "grey9", linewidth = 0.8) +  
#   geom_point(size = 3, color = "#D55E00") +
#   scale_x_log10() +
#   labs(x = "Hazard ratio (per 1 SD increase, log scale)", y = NULL) +
#   scale_y_discrete(labels = function(x) ifelse(x == "TCN2", paste0("**", x, "**"), x)) +
#   theme_pub +
#   theme(axis.text.y = ggtext::element_markdown(size = plot_base_size, face = "plain", color = "black"),
#         axis.text.x = element_text(size = plot_base_size, face = "plain", color = "black"))
# p 
# ggsave(file.path("Fig5B-discovery_top_forestplot.png"), p, width = 7,
#          height = 8, dpi = 300)
# ggsave(file.path("Fig5B-discovery_top_forestplot.pdf"), p, width = 7,
#        height = 8, bg = "white")
# 
# 
# 
# #--------5C--------
# plot_base_size = 15
# load(file = "Fig5C-rep_df.Rdata")
# long_df <- rbind(
#   data.frame(
#     gene = rep_df$gene,
#     cohort = "Discovery",
#     HR = rep_df$HR_discovery,
#     HR_L = rep_df$HR_L_discovery,
#     HR_H = rep_df$HR_H_discovery
#   ),
#   data.frame(
#     gene = rep_df$gene,
#     cohort = "Replication",
#     HR = rep_df$HR_replication,
#     HR_L = rep_df$HR_L_replication,
#     HR_H = rep_df$HR_H_replication
#   )
# )
# long_df$gene <- factor(long_df$gene, levels = rev(unique(rep_df$gene)))
# long_df$cohort <- factor(long_df$cohort, levels = c("Discovery", "Replication"))
#   
# theme_pub <- theme_classic(base_size = plot_base_size) +
#   theme(plot.title = element_blank(),
#           panel.grid = element_blank(),
#           axis.title = element_text(size = plot_base_size, face = "plain"),
#           axis.text = element_text(size = plot_base_size, face = "plain"),
#           legend.title = element_text(size = plot_base_size, face = "plain"),
#           legend.text = element_text(size = plot_base_size, face = "plain"))
# p <- ggplot(long_df, aes(x = HR, y = gene, color = cohort)) +
#   geom_vline(xintercept = 1, linetype = 2, color = "grey60", linewidth = 1) +
#   geom_errorbarh(aes(xmin = HR_L, xmax = HR_H), height = 0.20, 
#                  position = position_dodge(width = 0.55), linewidth = 1) +
#   geom_point(position = position_dodge(width = 0.55), size = 2.5) +
#   scale_x_log10() +
#   scale_color_manual(values = c("Discovery" = "#D55E00",   
#                                 "Replication" = "#0072B2"), 
#                      name = "Cohort") +
#   labs(x = "Hazard ratio (per 1 SD increase, log scale)", y = NULL, color = NULL) +
#   theme_pub +
#   theme(axis.text = element_text(color = "black")) 
# p
# ggsave(file.path("Fig3C-replicated_gene_forestplot.png"), p, width = 5.5,
#        height = 2.5, dpi = 300)
# ggsave(file.path("Fig3C-replicated_gene_forestplot.pdf"), p, width = 5.5,
#        height = 2.5)
# 
# #--------5D and 5E--------
# #data_dir <- "E:/ESCC_SNP"X
# data_dir <- "/Users/L/Desktop/ESCC/ESCC代码和数据/data"
# discovery_file  <- file.path("GSE53625_clean_cancer.RData")
# validation_file <- file.path("ESCC_two_tcga_extracted_clean.RData")
# discovery_label  <- "Discovery: GSE53625"
# validation_label <- "Validation: TCGA legacy ESCC"
# load(file = "Fig5C-rep_df.Rdata")
# merged_results_file <- rep_df
# gene_to_plot <- "TCN2"   
# split_method <- "median"
# min_group_size <- 5
# cex_base <- 1.3
# 
# 
# merged_df <- merged_results_file
# if (!is.null(gene_to_plot) && gene_to_plot %in% merged_df$gene) {
#   gene_selected <- gene_to_plot
# } else if (any(merged_df$replicated, na.rm = TRUE)) {
#   tmp <- merged_df[merged_df$replicated, ]
#   gene_selected <- tmp$gene[order(tmp$replication_pval, tmp$discovery_pval)][1]
# } else if (any(merged_df$direction_consistent, na.rm = TRUE)) {
#   tmp <- merged_df[merged_df$direction_consistent, ]
#   gene_selected <- tmp$gene[order(tmp$replication_pval, tmp$discovery_pval)][1]
# } else {
#   gene_selected <- merged_df$gene[order(merged_df$discovery_pval, merged_df$replication_pval)][1]
# }
# message("Selected gene: ", gene_selected)
# 
# load( file = "Fig5_dis_and_val_expr_clin.Rdata")
# 
# plot_km_with_risk <- function(expr_mat, clin_df, gene, cohort_label) {
#   if (!gene %in% rownames(expr_mat)) return(NULL)
#   x <- as.numeric(expr_mat[gene, ])
#   keep <- is.finite(x) & is.finite(clin_df$OS.time) & is.finite(clin_df$OS.event)
#   x <- x[keep]; time <- clin_df$OS.time[keep]; event <- clin_df$OS.event[keep]
#   if (length(unique(x)) < 2) return(NULL)
#   cutv <- median(x, na.rm = TRUE)
#   group <- factor(ifelse(x > cutv, "High", "Low"), levels = c("Low", "High"))
#   if (any(table(group) < 5)) return(NULL)
#   fit <- survfit(Surv(time, event) ~ group)
#   
#   my_theme <- theme_bw() +
#     theme(panel.grid = element_blank(),
#           axis.title = element_text(size = 15),
#           axis.text = element_text(size = 15),
#           legend.title = element_text(size = 15),
#           legend.text = element_text(size = 15),
#           plot.title = element_text(size = 15))
#   my_table_theme <- theme_cleantable() +
#     theme(axis.text = element_text(size = 15),
#           text = element_text(size = 15))
#   p <- ggsurvplot(fit, 
#                   data = data.frame(time = time, event = event, group = group),
#                   risk.table = TRUE,
#                   risk.table.col = "strata",
#                   pval = TRUE,
#                   pval.method = FALSE,
#                   conf.int = FALSE,
#                   palette = c("#0072B2", "#D55E00"),
#                   xlab = "Overall survival time",
#                   ylab = "Survival probability",
#                   legend.title = cohort_label,
#                   legend.labs = c("Low", "High"),
#                   risk.table.y.text.col = TRUE,
#                   risk.table.y.text = FALSE,
#                   tables.theme = my_table_theme,
#                   ggtheme = my_theme,
#                   font.x = 15,
#                   font.y = 15,
#                   font.tickslab = 15,
#                   font.legend = 15,
#                   font.pval = 15,
#                   font.risk.table = 15)
#   if (!is.null(p$table)) {
#     time_range <- range(c(0, time), na.rm = TRUE)
#     breaks <- pretty(time_range)
#     p$table <- p$table +
#       xlab("Time") +
#       scale_x_continuous(breaks = breaks, limits = time_range) +
#       theme(axis.text.x = element_text(size = 15),
#             axis.title.x = element_text(size = 15))
#   }
#   return(p)
# }
# 
# png(file.path("Fig5D-KM-dis.png"), width=1600, height=1600, res=300)
# plot_km_with_risk(expr_dis, clin_dis, gene_selected, discovery_label)
# dev.off()
# 
# png(file.path("Fig5E-KM-val.png"), width=1600, height=1600, res=300)
# plot_km_with_risk(expr_val, clin_val, gene_selected, validation_label)
# dev.off()
# 
# pdf(file.path("Fig5D-KM-dis.pdf"), width=6, height=6)
# plot_km_with_risk(expr_dis, clin_dis, gene_selected, discovery_label)
# dev.off()
# 
# pdf(file.path("Fig5E-KM-val.pdf"), width=6, height=6)
# plot_km_with_risk(expr_val, clin_val, gene_selected, validation_label)
# dev.off()
# 
# 
