# ==============================================================================
# Script 07: Locus-Specific TFAP2A and ZNF471 Methylation Analysis
# Project: OPMD Methylation Analysis (GSE204943)
# Author: Michael E. Troka
#
# Purpose:
#   Test the central hypothesis: ZNF471 promoter hypermethylation (silencing)
#   and TFAP2A derepression/activation occur progressively Normal->OPL->Tumor.
#   Extracts all CpGs at target loci, stratifies by promoter vs gene body,
#   performs pairwise statistical tests, and generates the epigenetic switch
#   scatter plot
#
# Inputs:
#   - analysis_objects.rds (from Script 01)
#
# Outputs:
#   - TFAP2A_ZNF471_summary.csv
#   - CpG detail CSVs
#   - gene boxplots
#   - trajectory plots
#   - epigenetic switch scatter
#
# Dependencies:
#   minfi, ggplot2, ggbeeswarm, ggsignif, patchwork, dplyr
#
# Upstream: Script 01
# ==============================================================================

############################################################
# TFAP2A and ZNF471 Locus-Specific Methylation Analysis
#
# Purpose: Directly test the central hypothesis that:
#   1. ZNF471 promoter becomes hypermethylated (silenced)
#   2. TFAP2A is derepressed and activated
#   3. These changes occur progressively: Normal → OPL → Tumor
#
# This analysis extracts all CpGs at these loci and examines
# promoter vs gene body methylation patterns.
############################################################

########################
## 0. Load packages
########################

suppressPackageStartupMessages({
  library(minfi)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(ggbeeswarm)
  library(ggsignif)
  library(RColorBrewer)
  library(scales)
})

########################
## 1. Load saved analysis objects
########################

cat("========================================\n")
cat("TFAP2A/ZNF471 LOCUS-SPECIFIC ANALYSIS\n")
cat("========================================\n\n")

# Adjust paths as needed
results_dir <- "GSE204943_data/GSE204943/results"
fig_dir <- "GSE204943_data/GSE204943/figures"

# Create output subdirectory
locus_fig_dir <- file.path(fig_dir, "locus_analysis")
dir.create(locus_fig_dir, showWarnings = FALSE, recursive = TRUE)

# Load analysis objects
analysis <- readRDS(file.path(results_dir, "analysis_objects.rds"))

beta     <- analysis$beta
pd       <- analysis$pd
Mset.noob <- analysis$Mset.noob

# Get annotation
ann <- getAnnotation(Mset.noob)

# Set up groups
group <- factor(pd$Sample_Group)
group <- relevel(group, ref = "3_Normal")

cat("Data loaded successfully!\n")
cat("Samples:", ncol(beta), "\n")
cat("CpGs:", nrow(beta), "\n\n")

########################
## 2. Define target genes and regions
########################

# Key genes from the hypothesis
target_genes <- c(
  "TFAP2A",   # Central transcription factor - should show activation pattern

  "ZNF471",   # Tumor suppressor - should show promoter hypermethylation
  "ZNF582",   # Another ZNF on 19q13.43 cluster - often co-methylated
  "MYC",      # Downstream of ZNF471 repression
  "HMGA1",    # TFAP2A target, metabolic regulator
  "CDKN2A",   # p16, key tumor suppressor
  "DNMT3B",   # De novo methyltransferase
  "VEGFA",    # Angiogenic target of TFAP2A/HIF1A
  "ZEB1",     # EMT factor regulated by TFAP2A
  "HIF1A"     # Hypoxia master regulator, cooperates with TFAP2A
)

cat("Target genes for analysis:\n")
cat(paste(target_genes, collapse = ", "), "\n\n")

########################
## 3. Extract CpGs for each gene
########################

# Function to extract CpGs for a gene
extract_gene_cpgs <- function(gene_symbol, annotation) {
  # Get CpGs associated with the gene
  # UCSC_RefGene_Name can contain multiple genes separated by ";"
  gene_pattern <- paste0("(^|;)", gene_symbol, "(;|$)")

  gene_cpgs <- annotation %>%
    as.data.frame() %>%
    filter(grepl(gene_pattern, UCSC_RefGene_Name)) %>%
    mutate(
      Gene = gene_symbol,
      # Extract the specific region for this gene
      Region = sapply(seq_len(n()), function(i) {
        names <- strsplit(UCSC_RefGene_Name[i], ";")[[1]]
        groups <- strsplit(UCSC_RefGene_Group[i], ";")[[1]]
        idx <- which(names == gene_symbol)[1]
        if (!is.na(idx) && idx <= length(groups)) groups[idx] else groups[1]
      })
    ) %>%
    select(Name, Gene, Region, chr, pos, strand,
           Relation_to_Island, UCSC_RefGene_Name, UCSC_RefGene_Group)

  return(gene_cpgs)
}

# Extract CpGs for all target genes
all_gene_cpgs <- lapply(target_genes, function(g) {
  cpgs <- extract_gene_cpgs(g, ann)
  if (nrow(cpgs) > 0) {
    cat(sprintf("  %s: %d CpGs\n", g, nrow(cpgs)))
  } else {
    cat(sprintf("  %s: No CpGs found\n", g))
  }
  return(cpgs)
})
names(all_gene_cpgs) <- target_genes

# Combine into single dataframe
gene_cpg_df <- bind_rows(all_gene_cpgs)
cat("\nTotal gene-associated CpGs:", nrow(gene_cpg_df), "\n\n")

########################
## 4. Define promoter vs gene body regions
########################

# Classify regions into functional categories
classify_region <- function(region) {
  case_when(
    region %in% c("TSS200", "TSS1500", "5'UTR", "1stExon") ~ "Promoter",
    region %in% c("Body", "ExonBnd") ~ "Gene Body",
    region == "3'UTR" ~ "3'UTR",
    TRUE ~ "Other"
  )
}

gene_cpg_df <- gene_cpg_df %>%
  mutate(
    Region_Class = classify_region(Region),
    Region_Class = factor(Region_Class, levels = c("Promoter", "Gene Body", "3'UTR", "Other"))
  )

# Summary of regions per gene
cat("CpGs by gene and region:\n")
region_summary <- gene_cpg_df %>%
  group_by(Gene, Region_Class) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = Region_Class, values_from = n, values_fill = 0)
print(region_summary)
cat("\n")

########################
## 5. Calculate mean methylation per sample
########################

sample_ids <- as.character(pd$Sample_Name)   # or pd$Basename

stopifnot(length(sample_ids) == ncol(beta))
if (anyDuplicated(sample_ids)) stop("Duplicate Sample_Name values in pd")

colnames(beta) <- sample_ids

# Create sample_info from the full beta matrix
sample_info <- data.frame(
  Sample = colnames(beta),
  Group  = as.character(group),
  stringsAsFactors = FALSE
)

# Safely subset CpGs (prevents empty/NA indexing)
gene_cpgs <- intersect(
  na.omit(gene_cpg_df$Name),
  rownames(beta)
)

if (length(gene_cpgs) == 0)
  stop("No matching CpGs found between gene_cpg_df and beta rownames")

# Get beta values for gene CpGs
gene_beta <- beta[gene_cpgs, , drop = FALSE]

# Calculate mean beta per gene, region, and sample
calculate_gene_methylation <- function(gene_name, region_class = NULL) {

  if (is.null(region_class)) {
    cpg_subset <- gene_cpg_df %>% filter(Gene == gene_name)
  } else {
    cpg_subset <- gene_cpg_df %>% filter(Gene == gene_name, Region_Class == region_class)
  }

  if (nrow(cpg_subset) == 0) return(NULL)

  cpg_names <- cpg_subset$Name
  beta_subset <- gene_beta[cpg_names, , drop = FALSE]

  # Mean across CpGs for each sample
  mean_beta <- colMeans(beta_subset, na.rm = TRUE)

  data.frame(
    Sample = names(mean_beta),
    Gene = gene_name,
    Region = ifelse(is.null(region_class), "All", region_class),
    Mean_Beta = mean_beta,
    stringsAsFactors = FALSE
  ) %>%
    left_join(sample_info, by = "Sample")
}

# Calculate for key genes with promoter and gene body separately
key_genes <- c("TFAP2A", "ZNF471", "ZNF582", "CDKN2A", "MYC", "HMGA1")

methylation_data <- bind_rows(
  lapply(key_genes, function(g) {
    bind_rows(
      calculate_gene_methylation(g, "Promoter"),
      calculate_gene_methylation(g, "Gene Body")
    )
  })
) %>%
  filter(!is.na(Mean_Beta))

cat("Methylation data calculated for", length(unique(methylation_data$Gene)), "genes\n\n")

########################
## 6. Statistical testing
########################

# Function to perform pairwise comparisons
perform_comparisons <- function(data, gene_name, region) {
  subset_data <- data %>%
    filter(Gene == gene_name, Region == region)

  if (nrow(subset_data) < 10) return(NULL)

  # Pairwise Wilcoxon tests
  comparisons <- list(
    c("3_Normal", "2_OPL"),
    c("2_OPL", "1_Tumor"),
    c("3_Normal", "1_Tumor")
  )

  results <- lapply(comparisons, function(comp) {
    g1 <- subset_data$Mean_Beta[subset_data$Group == comp[1]]
    g2 <- subset_data$Mean_Beta[subset_data$Group == comp[2]]

    if (length(g1) < 3 || length(g2) < 3) return(NULL)

    test <- wilcox.test(g1, g2)

    data.frame(
      Gene = gene_name,
      Region = region,
      Comparison = paste(comp, collapse = " vs "),
      Group1 = comp[1],
      Group2 = comp[2],
      Mean1 = mean(g1, na.rm = TRUE),
      Mean2 = mean(g2, na.rm = TRUE),
      Delta = mean(g2, na.rm = TRUE) - mean(g1, na.rm = TRUE),
      p_value = test$p.value,
      stringsAsFactors = FALSE
    )
  })

  bind_rows(results)
}

# Run comparisons for all genes and regions
stat_results <- bind_rows(
  lapply(key_genes, function(g) {
    bind_rows(
      perform_comparisons(methylation_data, g, "Promoter"),
      perform_comparisons(methylation_data, g, "Gene Body")
    )
  })
)

# Adjust p-values
stat_results <- stat_results %>%
  mutate(
    FDR = p.adjust(p_value, method = "BH"),
    Significance = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01 ~ "**",
      FDR < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

cat("Statistical comparisons:\n")
print(stat_results %>% select(Gene, Region, Comparison, Delta, p_value, FDR, Significance))
cat("\n")

# Save statistics
write.csv(stat_results, file.path(results_dir, "TFAP2A_ZNF471_statistics.csv"), row.names = FALSE)

########################
## 7. Visualization: Main hypothesis figure
########################

# Color palette
group_colors <- c(
  "3_Normal" = "#2166AC",  # Blue
  "2_OPL" = "#F6AE2D",     # Yellow/Orange
  "1_Tumor" = "#C41E3A"    # Red
)

# Clean group labels
clean_labels <- c(
  "3_Normal" = "Normal",
  "2_OPL" = "OPL",
  "1_Tumor" = "Tumor"
)

# Function to create gene methylation boxplot
create_gene_boxplot <- function(data, gene_name, stats_df, title_suffix = "") {

  plot_data <- data %>%
    filter(Gene == gene_name) %>%
    mutate(
      Group_Label = clean_labels[as.character(Group)],
      Group_Label = factor(Group_Label, levels = c("Normal", "OPL", "Tumor"))
    )

  if (nrow(plot_data) == 0) return(NULL)

  # Get significance for this gene
  gene_stats <- stats_df %>%
    filter(Gene == gene_name, Comparison == "3_Normal vs 1_Tumor")

  p <- ggplot(plot_data, aes(x = Group_Label, y = Mean_Beta, fill = Group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
    geom_beeswarm(size = 1.5, alpha = 0.6, cex = 2) +
    facet_wrap(~Region, ncol = 2) +
    scale_fill_manual(values = group_colors, guide = "none") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(
      title = paste0(gene_name, title_suffix),
      x = NULL,
      y = "Mean β-value"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      strip.text = element_text(face = "bold", size = 10),
      strip.background = element_rect(fill = "gray95", color = NA),
      axis.text = element_text(color = "black"),
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3)
    )

  # Add significance annotations
  for (region in unique(plot_data$Region)) {
    region_stats <- gene_stats %>% filter(Region == region)
    if (nrow(region_stats) > 0 && region_stats$FDR[1] < 0.05) {
      # Would add significance bars here
    }
  }

  return(p)
}

########################
## 8. Create main hypothesis figure: TFAP2A vs ZNF471
########################

cat("Creating main hypothesis figure...\n")

# TFAP2A plot
p_tfap2a <- create_gene_boxplot(methylation_data, "TFAP2A", stat_results)

# ZNF471 plot
p_znf471 <- create_gene_boxplot(methylation_data, "ZNF471", stat_results)

# Combine
main_fig <- p_znf471 / p_tfap2a +
  plot_annotation(
    title = "Locus-Specific Methylation: ZNF471 Silencing and TFAP2A Regulation",
    subtitle = "Testing the epigenetic switch hypothesis in oral carcinogenesis",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray30")
    )
  )

ggsave(file.path(locus_fig_dir, "TFAP2A_ZNF471_main_hypothesis.pdf"),
       main_fig, width = 8, height = 10, dpi = 300)
ggsave(file.path(locus_fig_dir, "TFAP2A_ZNF471_main_hypothesis.png"),
       main_fig, width = 8, height = 10, dpi = 300)

########################
## 9. Extended figure: All key genes
########################

cat("Creating extended figure for all key genes...\n")

# Create individual plots
gene_plots <- lapply(key_genes, function(g) {
  create_gene_boxplot(methylation_data, g, stat_results)
})
names(gene_plots) <- key_genes

# Remove NULL plots
gene_plots <- gene_plots[!sapply(gene_plots, is.null)]

# Combine into multi-panel figure
if (length(gene_plots) >= 4) {
  extended_fig <- wrap_plots(gene_plots, ncol = 2) +
    plot_annotation(
      title = "Methylation at Key Loci in the TFAP2A-Metabolic Axis",
      subtitle = "Promoter hypermethylation (silencing) vs Gene body patterns",
      theme = theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray30")
      )
    )

  ggsave(file.path(locus_fig_dir, "Key_genes_methylation_extended.pdf"),
         extended_fig, width = 12, height = 14, dpi = 300)
  ggsave(file.path(locus_fig_dir, "Key_genes_methylation_extended.png"),
         extended_fig, width = 12, height = 14, dpi = 300)
}

########################
## 10. Trajectory plot: Progressive changes
########################

cat("Creating trajectory plots...\n")

# Calculate group means for trajectory
trajectory_data <- methylation_data %>%
  group_by(Gene, Region, Group) %>%
  summarise(
    Mean = mean(Mean_Beta, na.rm = TRUE),
    SD = sd(Mean_Beta, na.rm = TRUE),
    SE = SD / sqrt(n()),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    Group_Numeric = case_when(
      Group == "3_Normal" ~ 1,
      Group == "2_OPL" ~ 2,
      Group == "1_Tumor" ~ 3
    ),
    Group_Label = clean_labels[as.character(Group)]
  )

# Focus on key hypothesis genes
hypothesis_genes <- c("ZNF471", "TFAP2A", "CDKN2A", "MYC")

trajectory_subset <- trajectory_data %>%
  filter(Gene %in% hypothesis_genes)

# Create trajectory plot
p_trajectory <- ggplot(trajectory_subset,
                       aes(x = Group_Numeric, y = Mean, color = Region, group = Region)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.1, linewidth = 0.8) +
  facet_wrap(~Gene, ncol = 2, scales = "free_y") +
  scale_x_continuous(
    breaks = 1:3,
    labels = c("Normal", "OPL", "Tumor")
  ) +
  scale_color_manual(
    values = c("Promoter" = "#D73027", "Gene Body" = "#4575B4"),
    name = "Region"
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Progressive Methylation Changes at Key Loci",
    subtitle = "Normal → OPL → Tumor trajectory",
    x = "Disease Stage",
    y = "Mean β-value (± SE)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40"),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "gray95", color = NA),
    legend.position = "bottom",
    axis.text = element_text(color = "black"),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3)
  )

ggsave(file.path(locus_fig_dir, "Trajectory_key_genes.pdf"),
       p_trajectory, width = 9, height = 8, dpi = 300)
ggsave(file.path(locus_fig_dir, "Trajectory_key_genes.png"),
       p_trajectory, width = 9, height = 8, dpi = 300)

########################
## 11. The "Epigenetic Switch" scatter plot
########################

cat("Creating epigenetic switch scatter plot...\n")

# Get per-sample data for ZNF471 promoter and TFAP2A (any region)
znf471_promoter <- methylation_data %>%
  filter(Gene == "ZNF471", Region == "Promoter") %>%
  select(Sample, Group, ZNF471_Promoter = Mean_Beta)

tfap2a_body <- methylation_data %>%
  filter(Gene == "TFAP2A", Region == "Gene Body") %>%
  select(Sample, TFAP2A_Body = Mean_Beta)

tfap2a_promoter <- methylation_data %>%
  filter(Gene == "TFAP2A", Region == "Promoter") %>%
  select(Sample, TFAP2A_Promoter = Mean_Beta)

# Merge
switch_data <- znf471_promoter %>%
  left_join(tfap2a_body, by = "Sample") %>%
  left_join(tfap2a_promoter, by = "Sample") %>%
  filter(!is.na(ZNF471_Promoter))

# If we have both promoter and body data for TFAP2A
if (sum(!is.na(switch_data$TFAP2A_Body)) > 10) {

  p_switch <- ggplot(switch_data, aes(x = ZNF471_Promoter, y = TFAP2A_Body, color = Group)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed", aes(group = 1), color = "gray40") +
    scale_color_manual(
      values = group_colors,
      labels = clean_labels,
      name = "Group"
    ) +
    labs(
      title = "The Epigenetic Switch: ZNF471 Silencing ↔ TFAP2A Activation",
      subtitle = "ZNF471 promoter methylation vs TFAP2A gene body methylation",
      x = "ZNF471 Promoter Methylation (β)",
      y = "TFAP2A Gene Body Methylation (β)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40"),
      legend.position = "right"
    ) +
    annotate("text", x = 0.15, y = 0.85, label = "Normal\n(ZNF471 active,\nTFAP2A repressed)",
             color = "#2166AC", size = 3, fontface = "italic", hjust = 0) +
    annotate("text", x = 0.7, y = 0.3, label = "Tumor\n(ZNF471 silenced,\nTFAP2A active)",
             color = "#C41E3A", size = 3, fontface = "italic", hjust = 1)

  ggsave(file.path(locus_fig_dir, "Epigenetic_switch_ZNF471_TFAP2A.pdf"),
         p_switch, width = 8, height = 7, dpi = 300)
  ggsave(file.path(locus_fig_dir, "Epigenetic_switch_ZNF471_TFAP2A.png"),
         p_switch, width = 8, height = 7, dpi = 300)

  # Calculate correlation
  cor_test <- cor.test(switch_data$ZNF471_Promoter, switch_data$TFAP2A_Body,
                       use = "complete.obs", method = "spearman")
  cat("\nCorrelation between ZNF471 promoter and TFAP2A body methylation:\n")
  cat(sprintf("  Spearman rho = %.3f, p = %.2e\n", cor_test$estimate, cor_test$p.value))
}

########################
## 12. Summary statistics table
########################

cat("\n========================================\n")
cat("SUMMARY OF KEY FINDINGS\n")
cat("========================================\n\n")

# Create summary table
summary_table <- stat_results %>%
  filter(Comparison == "3_Normal vs 1_Tumor") %>%
  select(Gene, Region, Delta, FDR, Significance) %>%
  mutate(
    Direction = ifelse(Delta > 0, "↑ Hyper", "↓ Hypo"),
    Interpretation = case_when(
      Gene == "ZNF471" & Region == "Promoter" & Delta > 0 ~ "SILENCED (supports hypothesis)",
      Gene == "TFAP2A" & Region == "Gene Body" & Delta > 0 ~ "ACTIVE (supports hypothesis)",
      Gene == "CDKN2A" & Region == "Promoter" & Delta > 0 ~ "SILENCED (tumor suppressor loss)",
      Gene == "MYC" & Delta != 0 ~ "Potential activation",
      TRUE ~ ""
    )
  ) %>%
  arrange(Gene, Region)

cat("Methylation changes (Normal → Tumor):\n\n")
print(summary_table)

# Save summary
write.csv(summary_table, file.path(results_dir, "TFAP2A_ZNF471_summary.csv"), row.names = FALSE)

########################
## 13. Individual CpG analysis for ZNF471 and TFAP2A
########################

cat("\n\nDetailed CpG-level analysis for key genes...\n")

# Function to create CpG-level heatmap data
create_cpg_detail <- function(gene_name) {

  gene_cpgs <- gene_cpg_df %>% filter(Gene == gene_name)

  if (nrow(gene_cpgs) == 0) return(NULL)

  # Get beta values
  cpg_beta <- beta[gene_cpgs$Name, , drop = FALSE]

  # Calculate mean by group
  group_means <- sapply(levels(group), function(g) {
    rowMeans(cpg_beta[, group == g, drop = FALSE], na.rm = TRUE)
  })

  # Create detailed table
  detail_df <- gene_cpgs %>%
    mutate(
      Normal = group_means[Name, "3_Normal"],
      OPL = group_means[Name, "2_OPL"],
      Tumor = group_means[Name, "1_Tumor"],
      Delta_TvN = Tumor - Normal,
      Progressive = sign(OPL - Normal) == sign(Tumor - Normal)
    ) %>%
    arrange(Region, pos)

  return(detail_df)
}

# Get details for ZNF471 and TFAP2A
znf471_detail <- create_cpg_detail("ZNF471")
tfap2a_detail <- create_cpg_detail("TFAP2A")

if (!is.null(znf471_detail)) {
  cat("\nZNF471 CpG details:\n")
  print(znf471_detail %>% select(Name, Region, Normal, OPL, Tumor, Delta_TvN, Progressive))
  write.csv(znf471_detail, file.path(results_dir, "ZNF471_CpG_details.csv"), row.names = FALSE)
}

if (!is.null(tfap2a_detail)) {
  cat("\nTFAP2A CpG details:\n
")
  print(tfap2a_detail %>% select(Name, Region, Normal, OPL, Tumor, Delta_TvN, Progressive))
  write.csv(tfap2a_detail, file.path(results_dir, "TFAP2A_CpG_details.csv"), row.names = FALSE)
}

########################
## 14. Final output summary
########################

cat("\n========================================\n")
cat("ANALYSIS COMPLETE\n")
cat("========================================\n\n")

cat("Figures saved to:", locus_fig_dir, "\n")
cat("  - TFAP2A_ZNF471_main_hypothesis.pdf/png\n")
cat("  - Key_genes_methylation_extended.pdf/png\n")
cat("  - Trajectory_key_genes.pdf/png\n")
cat("  - Epigenetic_switch_ZNF471_TFAP2A.pdf/png\n")

cat("\nData saved to:", results_dir, "\n")
cat("  - TFAP2A_ZNF471_statistics.csv\n")
cat("  - TFAP2A_ZNF471_summary.csv\n")
cat("  - ZNF471_CpG_details.csv\n")
cat("  - TFAP2A_CpG_details.csv\n")

cat("\n========================================\n")
cat("KEY HYPOTHESIS TEST RESULTS\n")
cat("========================================\n\n")

# Print key findings
znf471_promoter_stat <- stat_results %>%
  filter(Gene == "ZNF471", Region == "Promoter", Comparison == "3_Normal vs 1_Tumor")

tfap2a_body_stat <- stat_results %>%
  filter(Gene == "TFAP2A", Region == "Gene Body", Comparison == "3_Normal vs 1_Tumor")

if (nrow(znf471_promoter_stat) > 0) {
  cat("ZNF471 Promoter (Normal → Tumor):\n")
  cat(sprintf("  Δβ = %.3f, FDR = %.2e %s\n",
              znf471_promoter_stat$Delta,
              znf471_promoter_stat$FDR,
              znf471_promoter_stat$Significance))
  if (znf471_promoter_stat$Delta > 0 && znf471_promoter_stat$FDR < 0.05) {
    cat("  → SUPPORTS HYPOTHESIS: ZNF471 promoter hypermethylation = silencing\n")
  }
}

if (nrow(tfap2a_body_stat) > 0) {
  cat("\nTFAP2A Gene Body (Normal → Tumor):\n")
  cat(sprintf("  Δβ = %.3f, FDR = %.2e %s\n",
              tfap2a_body_stat$Delta,
              tfap2a_body_stat$FDR,
              tfap2a_body_stat$Significance))
  if (tfap2a_body_stat$Delta > 0 && tfap2a_body_stat$FDR < 0.05) {
    cat("  → SUPPORTS HYPOTHESIS: TFAP2A gene body hypermethylation = active transcription\n")
  }
}

cat("\n")
