# ==============================================================================
# Script 25: Cell-Type Deconvolution (EpiDISH)
# Project: OPMD Methylation Analysis (GSE204943)
# Author: Michael E. Troka
#
# Purpose:
#   Estimate epithelial/fibroblast/immune cell proportions using EpiDISH.
#   Test if cell composition differs across groups, confounds key findings,
#   and whether ZNF471/TFAP2A methylation changes survive cell-type adjustment.
#
# Inputs:
#   - analysis_objects.rds (beta matrix, pd, Mval.bc, Mset.noob)
#
# Outputs:
#   - results/cell_type_proportions.csv
#   - results/DMP_sensitivity_celltype.csv
#   - figures/CellType_by_group.pdf
#   - figures/CellType_confounding_scatter.pdf
#   - results/source_data/SourceData_ED_CellType.xlsx
#
# Dependencies:
#   EpiDISH, ggplot2, dplyr, limma, openxlsx
#
# Upstream: Script 01
# ==============================================================================

suppressPackageStartupMessages({
  library(EpiDISH)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(limma)
  library(minfi)
  library(openxlsx)
  library(ggpubr)
})

# ============================================================
# 0. Paths
# ============================================================
results_dir <- "GSE204943_data/GSE204943/results"
fig_dir     <- "GSE204943_data/GSE204943/figures"
source_dir  <- file.path(results_dir, "source_data")
dir.create(source_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. Load data
# ============================================================
cat("Loading analysis objects...\n")
obj <- readRDS(file.path(results_dir, "analysis_objects.rds"))
beta   <- obj$beta
pd     <- obj$pd

# Fix empty/NULL colnames on beta: assign Basename IDs so sample matching works
sample_ids <- basename(pd$Basename)
colnames(beta) <- sample_ids
cat("Assigned sample IDs as beta colnames (", ncol(beta), "samples)\n")

group <- factor(pd$Sample_Group, levels = c("3_Normal", "2_OPL", "1_Tumor"))
group_label <- factor(group,
                       levels = c("3_Normal", "2_OPL", "1_Tumor"),
                       labels = c("Normal", "OPL", "Tumor"))

# ============================================================
# 2. Run EpiDISH deconvolution
# ============================================================
cat("\n--- Running EpiDISH RPC deconvolution ---\n")

# Use centEpiFibIC reference (epithelial / fibroblast / immune)
data(centEpiFibIC.m)

epidish_result <- epidish(
  beta.m = beta,
  ref.m = centEpiFibIC.m,
  method = "RPC"
)

cell_props <- as.data.frame(epidish_result$estF)
cell_props$Sample <- rownames(cell_props)
colnames(cell_props)[1:3] <- c("Epithelial", "Fibroblast", "Immune")

# Add group info — match Sample (Basename IDs) to pd$Basename
cell_props <- cell_props %>%
  mutate(
    Group = as.character(group[match(Sample, sample_ids)]),
    Group_Label = factor(Group,
                         levels = c("3_Normal", "2_OPL", "1_Tumor"),
                         labels = c("Normal", "OPL", "Tumor"))
  )

cat("Cell-type proportions computed for", nrow(cell_props), "samples\n")
cat("\nMean proportions by group:\n")
cell_props %>%
  group_by(Group_Label) %>%
  summarise(across(c(Epithelial, Fibroblast, Immune), mean), .groups = "drop") %>%
  print()

# ============================================================
# 3. Statistical tests: cell proportions across groups
# ============================================================
cat("\n--- Statistical tests ---\n")

cell_tests <- list()
for (ct in c("Epithelial", "Fibroblast", "Immune")) {
  # Kruskal-Wallis overall
  kw <- kruskal.test(cell_props[[ct]] ~ cell_props$Group)

  # Pairwise Wilcoxon
  pw <- pairwise.wilcox.test(cell_props[[ct]], cell_props$Group,
                              p.adjust.method = "BH")

  cat("\n", ct, ":\n")
  cat("  Kruskal-Wallis p =", format.pval(kw$p.value), "\n")
  cat("  Pairwise p-values:\n")
  print(pw$p.value)

  cell_tests[[ct]] <- list(kruskal = kw, pairwise = pw)
}

# ============================================================
# 4. Visualization: cell proportions by group
# ============================================================
cat("\n--- Generating cell-type plots ---\n")

cell_long <- cell_props %>%
  pivot_longer(cols = c(Epithelial, Fibroblast, Immune),
               names_to = "Cell_Type", values_to = "Proportion")

p_celltype <- ggplot(cell_long, aes(x = Group_Label, y = Proportion, fill = Group_Label)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 0.8, alpha = 0.4) +
  facet_wrap(~ Cell_Type, scales = "free_y") +
  scale_fill_manual(values = c("Normal" = "#3C5488", "OPL" = "#00A087", "Tumor" = "#E64B35")) +
  stat_compare_means(method = "kruskal.test", label = "p.format", size = 3.5) +
  labs(
    title = "Cell-Type Composition Across Disease Groups",
    subtitle = "EpiDISH RPC deconvolution (epithelial/fibroblast/immune reference)",
    y = "Estimated Proportion",
    x = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"))

ggsave(file.path(fig_dir, "CellType_by_group.pdf"), p_celltype,
       width = 10, height = 4, dpi = 300)
ggsave(file.path(fig_dir, "CellType_by_group.png"), p_celltype,
       width = 10, height = 4, dpi = 300)

# ============================================================
# 5. Confounding analysis: cell props vs PMD load & tumor prob
# ============================================================
cat("\n--- Confounding analysis ---\n")

# Load PMD load data
pmd_file <- file.path(results_dir, "PMD_load_per_sample.csv")
if (file.exists(pmd_file)) {
  pmd_load <- read.csv(pmd_file, stringsAsFactors = FALSE)

  cell_pmd <- cell_props %>%
    left_join(pmd_load %>% select(Sample, PMD_load_Mb), by = "Sample")

  # Spearman correlations
  for (ct in c("Epithelial", "Fibroblast", "Immune")) {
    ct_vals <- cell_pmd[[ct]]
    pmd_vals <- cell_pmd$PMD_load_Mb
    valid <- !is.na(ct_vals) & !is.na(pmd_vals)
    if (sum(valid) > 5) {
      cor_test <- cor.test(ct_vals[valid], pmd_vals[valid], method = "spearman")
      cat(ct, "vs PMD_load: rho =", round(cor_test$estimate, 3),
          ", p =", format.pval(cor_test$p.value), "\n")
    }
  }

  # Scatter plot
  scatter_data <- cell_pmd %>%
    pivot_longer(cols = c(Epithelial, Fibroblast, Immune),
                 names_to = "Cell_Type", values_to = "Proportion") %>%
    filter(!is.na(PMD_load_Mb))

  p_confound <- ggplot(scatter_data, aes(x = Proportion, y = PMD_load_Mb, color = Group_Label)) +
    geom_point(size = 1.5, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "black", linewidth = 0.5) +
    facet_wrap(~ Cell_Type, scales = "free_x") +
    scale_color_manual(values = c("Normal" = "#3C5488", "OPL" = "#00A087", "Tumor" = "#E64B35")) +
    labs(
      title = "Cell-Type Proportions vs PMD Erosion Load",
      x = "Cell-Type Proportion",
      y = "PMD Load (Mb)",
      color = "Group"
    ) +
    theme_classic(base_size = 11) +
    theme(legend.position = "bottom")

  ggsave(file.path(fig_dir, "CellType_confounding_scatter.pdf"), p_confound,
         width = 12, height = 4, dpi = 300)
  ggsave(file.path(fig_dir, "CellType_confounding_scatter.png"), p_confound,
         width = 12, height = 4, dpi = 300)
}

# ============================================================
# 6. Sensitivity: Re-run limma with cell-type covariates
# ============================================================
cat("\n--- Sensitivity analysis: limma with cell-type covariates ---\n")

# Re-extract Mval fresh from obj (earlier operations may corrupt the variable)
Mval <- obj$Mval.bc
cat("Mval re-loaded: dim =", dim(Mval), "\n")

# Align cell proportions to pd/Mval sample order (same as sample_ids)
cell_props_aligned <- cell_props[match(sample_ids, cell_props$Sample), ]
stopifnot(nrow(cell_props_aligned) == ncol(Mval))

# Design matrix with group + cell-type proportions
design_ct <- model.matrix(~ 0 + group + cell_props_aligned$Epithelial + cell_props_aligned$Fibroblast)
colnames(design_ct) <- c("Normal", "OPL", "Tumor", "Epithelial", "Fibroblast")

# Contrast: Tumor vs Normal
contrast_mat <- makeContrasts(
  TvN = Tumor - Normal,
  OvN = OPL - Normal,
  levels = design_ct
)

# Fit
fit_ct <- lmFit(Mval, design_ct)
fit_ct2 <- contrasts.fit(fit_ct, contrast_mat)
fit_ct2 <- eBayes(fit_ct2)

# Extract results
dmps_ct_TvN <- topTable(fit_ct2, coef = "TvN", number = Inf)
dmps_ct_OvN <- topTable(fit_ct2, coef = "OvN", number = Inf)

sig_ct_TvN <- dmps_ct_TvN %>% filter(adj.P.Val < 0.05, abs(logFC) > 1)
sig_ct_OvN <- dmps_ct_OvN %>% filter(adj.P.Val < 0.05, abs(logFC) > 1)

# Compare with original (without cell-type adjustment)
sig_orig_TvN <- obj$sig$TvN
sig_orig_OvN <- obj$sig$OvN

sensitivity_summary <- data.frame(
  Comparison = c("Tumor vs Normal", "OPL vs Normal"),
  Original_DMPs = c(nrow(sig_orig_TvN), nrow(sig_orig_OvN)),
  CellType_Adjusted_DMPs = c(nrow(sig_ct_TvN), nrow(sig_ct_OvN)),
  Overlap = c(
    length(intersect(rownames(sig_orig_TvN), rownames(sig_ct_TvN))),
    length(intersect(rownames(sig_orig_OvN), rownames(sig_ct_OvN)))
  ),
  stringsAsFactors = FALSE
)
sensitivity_summary$Pct_Retained <- round(100 * sensitivity_summary$Overlap /
                                           sensitivity_summary$Original_DMPs, 1)

cat("\nDMP sensitivity to cell-type adjustment:\n")
print(sensitivity_summary)

# Check ZNF471 and TFAP2A specifically
ann_obj <- as.data.frame(getAnnotation(obj$Mset.noob))
znf471_cpgs <- rownames(ann_obj)[base::grepl("ZNF471", ann_obj$UCSC_RefGene_Name)]
tfap2a_cpgs <- rownames(ann_obj)[base::grepl("TFAP2A", ann_obj$UCSC_RefGene_Name)]

cat("\n--- ZNF471 CpGs after cell-type adjustment (TvN) ---\n")
znf471_ct <- dmps_ct_TvN[intersect(znf471_cpgs, rownames(dmps_ct_TvN)), ]
cat("Significant ZNF471 CpGs (FDR<0.05):",
    sum(znf471_ct$adj.P.Val < 0.05), "of", nrow(znf471_ct), "\n")

cat("\n--- TFAP2A CpGs after cell-type adjustment (TvN) ---\n")
tfap2a_ct <- dmps_ct_TvN[intersect(tfap2a_cpgs, rownames(dmps_ct_TvN)), ]
cat("Significant TFAP2A CpGs (FDR<0.05):",
    sum(tfap2a_ct$adj.P.Val < 0.05), "of", nrow(tfap2a_ct), "\n")

# ============================================================
# 7. Save results
# ============================================================
cat("\n--- Saving results ---\n")

write.csv(cell_props, file.path(results_dir, "cell_type_proportions.csv"),
          row.names = FALSE)
write.csv(sensitivity_summary, file.path(results_dir, "DMP_sensitivity_celltype.csv"),
          row.names = FALSE)

# ============================================================
# 7b. Correlate immune fractions with key methylation markers (Gap E)
# ============================================================
cat("\n--- Correlating immune fractions with key genes ---\n")

# Load PMD erosion data which has ZNF471/TFAP2A methylation
pmd_erosion_file <- file.path(results_dir, "PMD_erosion_score_per_sample.csv")
immune_corr_results <- data.frame(
  Cell_Type = character(), Methylation_Feature = character(),
  Spearman_rho = numeric(), p_value = numeric(), n = integer(),
  stringsAsFactors = FALSE
)

if (file.exists(pmd_erosion_file)) {
  pmd_erosion <- read.csv(pmd_erosion_file, stringsAsFactors = FALSE)
  pmd_erosion <- pmd_erosion %>% mutate(across(where(is.character), trimws))

  cell_meth <- cell_props %>%
    left_join(pmd_erosion %>% select(Sample, ZNF471_prom_beta, TFAP2A_prom_beta,
                                      TFAP2A_body_beta, PMD_load_Mb,
                                      TumorProb, PMD_erosion_score),
              by = "Sample")

  meth_features <- c("ZNF471_prom_beta", "TFAP2A_prom_beta", "TFAP2A_body_beta",
                       "PMD_load_Mb", "TumorProb", "PMD_erosion_score")
  cell_types <- c("Epithelial", "Fibroblast", "Immune")

  for (ct in cell_types) {
    for (mf in meth_features) {
      if (mf %in% names(cell_meth)) {
        valid <- !is.na(cell_meth[[ct]]) & !is.na(cell_meth[[mf]])
        if (sum(valid) >= 10) {
          ct_test <- cor.test(cell_meth[[ct]][valid], cell_meth[[mf]][valid],
                               method = "spearman", exact = FALSE)
          immune_corr_results <- rbind(immune_corr_results, data.frame(
            Cell_Type = ct, Methylation_Feature = mf,
            Spearman_rho = ct_test$estimate,
            p_value = ct_test$p.value,
            n = sum(valid),
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }

  if (nrow(immune_corr_results) > 0) {
    immune_corr_results$FDR <- p.adjust(immune_corr_results$p_value, method = "BH")
    immune_corr_results <- immune_corr_results %>% arrange(p_value)
    cat("\n--- Immune/cell-type correlations with methylation features ---\n")
    print(immune_corr_results)
    write.csv(immune_corr_results,
              file.path(results_dir, "celltype_methylation_correlations.csv"),
              row.names = FALSE)
  }
}

# ============================================================
# 8. Source data XLSX
# ============================================================
wb <- createWorkbook()

addWorksheet(wb, "README")
writeData(wb, "README", data.frame(
  Field = c("Analysis", "Method", "Reference", "Key Question", "Date"),
  Value = c("Cell-Type Deconvolution",
            "EpiDISH RPC (Robust Partial Correlation)",
            "centEpiFibIC.m (Epithelial/Fibroblast/Immune)",
            "Do cell-type shifts confound methylation findings?",
            format(Sys.time(), "%Y-%m-%d"))
))

addWorksheet(wb, "EDFig_CellType_Boxplot")
ct_summary <- cell_long %>%
  group_by(Cell_Type, Group_Label) %>%
  summarise(N = n(), Mean = mean(Proportion), SD = sd(Proportion),
            Median = median(Proportion), .groups = "drop")
writeData(wb, "EDFig_CellType_Boxplot", ct_summary)

addWorksheet(wb, "EDFig_CellType_Scatter")
if (exists("scatter_data")) {
  scatter_summary <- scatter_data %>%
    group_by(Cell_Type, Group_Label) %>%
    summarise(
      N = n(),
      Mean_Proportion = mean(Proportion),
      Mean_PMD = mean(PMD_load_Mb, na.rm = TRUE),
      .groups = "drop"
    )
  writeData(wb, "EDFig_CellType_Scatter", scatter_summary)
}

addWorksheet(wb, "Cell_Proportions")
writeData(wb, "Cell_Proportions", cell_props)

addWorksheet(wb, "Sensitivity_Summary")
writeData(wb, "Sensitivity_Summary", sensitivity_summary)

addWorksheet(wb, "Statistics")
stat_rows <- data.frame(
  Test = character(), Cell_Type = character(),
  Statistic = numeric(), P_value = numeric(),
  stringsAsFactors = FALSE
)
for (ct in names(cell_tests)) {
  stat_rows <- rbind(stat_rows, data.frame(
    Test = "Kruskal-Wallis",
    Cell_Type = ct,
    Statistic = cell_tests[[ct]]$kruskal$statistic,
    P_value = cell_tests[[ct]]$kruskal$p.value,
    stringsAsFactors = FALSE
  ))
}
writeData(wb, "Statistics", stat_rows)

if (nrow(immune_corr_results) > 0) {
  addWorksheet(wb, "CellType_Meth_Correlations")
  writeData(wb, "CellType_Meth_Correlations", immune_corr_results)
}

saveWorkbook(wb, file.path(source_dir, "SourceData_ED_CellType.xlsx"),
             overwrite = TRUE)

cat("\n=== Script 25 complete ===\n")
