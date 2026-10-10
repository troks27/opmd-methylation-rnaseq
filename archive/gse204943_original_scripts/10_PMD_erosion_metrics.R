# ==============================================================================
# Script 10: PMD Erosion Quantification
# Project: OPMD Methylation Analysis (GSE204943)
# Author: Michael E. Troka
#
# Purpose:
#   Identify PMD-like hypomethylated blocks using blockFinder (case vs control),
#   quantify per-sample PMD erosion metrics (load in Mb, inside vs outside
#   methylation delta), generate publication-quality PMD distribution and
#   metric plots
#
# Inputs:
#   - Mset.noob, beta matrix, pd (from Script 01)
#   - block coordinates
#
# Outputs:
#   - hypo_PMD_blocks.csv/.rds
#   - PMD_load_per_sample.csv
#   - PMD_inside_outside_stats_per_sample.csv
#   - PMD erosion figures (Fig1-Fig5)
#
# Dependencies:
#   minfi, limma, bumphunter, GenomicRanges, ggplot2, pheatmap, dplyr,
#   matrixStats
#
# Upstream: Script 09
# ==============================================================================

############################################################
# 14. PMD erosion analysis (OPMD vs Normal)
# Purpose:
#   - Identify PMD-like hypomethylated blocks (case vs control)
#   - Quantify per-sample PMD erosion metrics
#   - Generate publication-quality plots
############################################################

suppressPackageStartupMessages({
  library(minfi)
  library(GenomicRanges)
  library(IRanges)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
})

########################
## 14a. Parameters
########################
pmd_params <- list(
  # blockFinder settings
  type          = "OpenSea",
  cutoff_delta  = 0.10,   # min abs meanDiff to call a block
  lambda        = 1000,   # smoothing width (bp)
  C            = 2,       # cluster parameter

  # per-sample erosion call settings
  sample_beta_cutoff      = 0.50,  # mean beta within block below this => eroded in sample
  sample_delta_from_out   = 0.05,  # must be at least this much lower than outside-PMD mean beta

  # plotting
  palette_heatmap = colorRampPalette(c("navy","white","firebrick3"))(100)
)

########################
## 14b. Make a GenomicRatioSet (robustly)
########################
# blockFinder requires a GenomicRatioSet/GenomicMethylSet with genomic coordinates.
# We map Mset.noob to genome once.
message("Building GenomicRatioSet for blockFinder...")

if (!inherits(Mset.noob, "GenomicMethylSet")) {
  # Mset.noob from preprocessNoob is usually a MethylSet; map to genome.
  grSet <- mapToGenome(Mset.noob)
} else {
  grSet <- Mset.noob
}

# Ensure we have M-values in grSet
if (!"M" %in% assayNames(grSet)) {
  grSet <- ratioConvert(grSet, what = "M", keepCN = TRUE)
}

########################
## 14c. Define case vs control for OPMD erosion
########################
# Controls = Normal only
# Cases = OPL + Tumor (OPMD spectrum)
pd2 <- as.data.frame(colData(grSet))
pd2$case_status <- ifelse(pd2$Sample_Group %in% c("1_Tumor","2_OPL"), "Case", "Control")
pd2$case_status <- factor(pd2$case_status, levels = c("Control","Case"))
colData(grSet)$case_status <- pd2$case_status

design_pmd <- model.matrix(~ case_status, data = pd2)

########################
## 14d (PATCH v3). Use blockgroup as clusters + restore names
########################
message("Preparing OpenSea clusters with cpgCollapse (older minfi)...")

`%||%` <- function(a, b) if (!is.null(a)) a else b

cc <- cpgCollapse(grSet, returnBlockInfo = FALSE)  # your version ignores blockInfo anyway
# cc is a list: use the object slot
grSet_collapsed <- cc$object

# restore column names (helpful later)
colnames(grSet_collapsed) <- colnames(grSet)

# clusters live on rowRanges metadata as blockgroup
cluster_vec <- rowRanges(grSet_collapsed)$blockgroup

message("Clusters built: ", length(unique(cluster_vec[!is.na(cluster_vec)])))

keep_objs <- c(
  "grSet_collapsed", "design_pmd", "cluster_vec", "pmd_params",
  "results_dir", "fig_dir", "params"
)

rm(list = setdiff(ls(envir = .GlobalEnv), keep_objs), envir = .GlobalEnv)
gc()

cl <- makeCluster(5)
registerDoParallel(cl)

blocks_case_control <- blockFinder(
  object  = grSet_collapsed,
  design  = design_pmd,
  coef    = 2,
  cluster = cluster_vec,
  cutoff  = pmd_params$cutoff_delta,
  pickCutoff = FALSE,
  smooth  = TRUE,
  B       = 100,
  verbose = TRUE
)

stopCluster(cl)

str(blocks_case_control)

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
})

# ---- 14d-1) bumps -> data.frame -> GRanges ----
blocks_df <- blocks_case_control$table

blocks_gr <- GRanges(
  seqnames = blocks_df$chr,
  ranges   = IRanges(start = blocks_df$start, end = blocks_df$end),
  value    = blocks_df$value,
  area     = blocks_df$area,
  L        = blocks_df$L,
  p.value  = blocks_df$p.value,
  fwer     = blocks_df$fwer,
  cluster  = blocks_df$cluster
)

# ---- 14d-2) keep ONLY hypomethylated PMD-like blocks ----
# "value" is smoothed Case - Control; negative = hypomethylated in cases
hypo_pmd_gr <- blocks_gr[
  mcols(blocks_gr)$value < -pmd_params$cutoff_delta
]

message("Hypo-PMD blocks retained: ", length(hypo_pmd_gr))

saveRDS(hypo_pmd_gr, file.path(results_dir, "hypo_PMD_blocks.rds"))
write.csv(as.data.frame(hypo_pmd_gr),
          file.path(results_dir, "hypo_PMD_blocks.csv"),
          row.names = FALSE)

########################
## 14e. Compute inside- vs outside-PMD methylation per sample
########################

suppressPackageStartupMessages({
  library(minfi)
  library(GenomicRanges)
  library(dplyr)
})

# --- 0) Build GRanges-aware GenomicRatioSet from obj$Mset.noob ---
# obj$Mset.noob is a MethylSet, so convert -> mapToGenome
grSet <- ratioConvert(obj$Mset.noob, what = "Beta")
grSet <- mapToGenome(grSet)

# preserve sample names if they were dropped upstream
if (is.null(colnames(grSet)) || all(is.na(colnames(grSet)))) {
  colnames(grSet) <- obj$pd$Basename %||% obj$pd$Sample_Name %||% paste0("S", seq_len(ncol(grSet)))
}

# --- 1) Collapse CpGs into contiguous clusters/blocks (for blockFinder) ---
cc <- cpgCollapse(
  object = grSet,
  returnBlockInfo = TRUE
)
grSet_collapsed <- cc$object

# clusters live on rowRanges metadata as blockgroup
cluster_vec <- rowRanges(grSet_collapsed)$blockgroup

# sanity checks
message("Clusters built: ", length(unique(cluster_vec[!is.na(cluster_vec)])))
message("Collapsed dim: ", paste(dim(grSet_collapsed), collapse=" x "))

# --- 2) Recreate hypo_blocks_gr from your blockFinder output ---
# blocks_case_control is a 'bumps' object
blocks_tbl <- blocks_case_control$table

# hypo blocks = negative Case-Control delta beta
hypo_tbl <- subset(blocks_tbl, value < 0)

# (optional) apply significance filter if desired
# hypo_tbl <- subset(blocks_tbl, value < 0 & fwerArea < 0.05)

hypo_blocks_gr <- GRanges(
  seqnames = hypo_tbl$chr,
  ranges   = IRanges(start = hypo_tbl$start, end = hypo_tbl$end),
  value    = hypo_tbl$value,
  area     = hypo_tbl$area,
  cluster  = hypo_tbl$cluster,
  L        = hypo_tbl$L
)

# ----- FIX: restore sample names on grSet BEFORE getBeta -----
# pick the best available name column
sample_names <- obj$pd$Basename
if (is.null(sample_names) || all(is.na(sample_names)) || all(sample_names == "")) {
  sample_names <- obj$pd$Sample_Name
}
if (is.null(sample_names) || length(sample_names) != ncol(grSet)) {
  sample_names <- paste0("S", seq_len(ncol(grSet)))
}

colnames(grSet) <- sample_names
# ------------------------------------------------------------

beta_gr  <- getBeta(grSet)
probe_gr <- rowRanges(grSet)

# --- 3) Get beta matrix aligned to probes and probe GRanges ---
beta_gr  <- getBeta(grSet)
probe_gr <- rowRanges(grSet)

# --- 4) Overlaps of probes with hypo-PMD blocks ---
ol <- findOverlaps(probe_gr, hypo_blocks_gr)

in_block_probes  <- unique(queryHits(ol))
out_block_probes <- setdiff(seq_along(probe_gr), in_block_probes)

# guard against empty sets
stopifnot(length(in_block_probes)  > 0)
stopifnot(length(out_block_probes) > 0)

# --- 5) Per-sample means inside vs outside PMDs ---
mean_in  <- colMeans(beta_gr[in_block_probes,  , drop = FALSE], na.rm = TRUE)
mean_out <- colMeans(beta_gr[out_block_probes, , drop = FALSE], na.rm = TRUE)
delta_io <- mean_in - mean_out

pmd_sample_stats <- data.frame(
  Sample              = colnames(beta_gr),
  mean_beta_in_PMD    = mean_in,
  mean_beta_out_PMD   = mean_out,
  delta_in_minus_out  = delta_io,
  stringsAsFactors = FALSE
)

# --- 6) Join to phenotype (use obj$pd unless you have pd2 explicitly) ---
pd_for_join <- obj$pd %>%
  mutate(Sample = Basename %||% Sample_Name)

pmd_sample_stats <- pmd_sample_stats %>%
  left_join(pd_for_join, by = "Sample")

write.csv(
  pmd_sample_stats,
  file.path(results_dir, "PMD_inside_outside_stats_per_sample.csv"),
  row.names = FALSE
)

########################
## 14f. Per-sample PMD load and number of eroded blocks
########################
library(dplyr)
library(GenomicRanges)

# EXPECTS from 14e:
# - probe_gr (GRanges of probes)
# - hypo_blocks_gr (GRanges of hypo-PMD blocks)
# - beta_gr (beta matrix aligned to probe_gr; rows=probes, cols=samples)
# - mean_out (named vector: mean beta outside-PMD per sample)
# - pmd_params$sample_beta_cutoff
# - pmd_params$sample_delta_from_out
# - results_dir
# - obj$pd (phenotype table)

stopifnot(
  exists("probe_gr"),
  exists("hypo_blocks_gr"),
  exists("beta_gr"),
  exists("mean_out")
)

# --- 1) Build list of probe indices per block (efficiently, one overlap call) ---
ol_blocks <- findOverlaps(probe_gr, hypo_blocks_gr)
block_probe_list <- split(queryHits(ol_blocks), subjectHits(ol_blocks))

# Ensure list is in block order (fill empty blocks, if any)
block_probe_list <- lapply(seq_along(hypo_blocks_gr), function(i) block_probe_list[[as.character(i)]] %||% integer(0))

# --- 2) Block-by-sample mean beta matrix ---
block_beta_mat <- sapply(block_probe_list, function(idx) {
  if (length(idx) == 0) {
    rep(NA_real_, ncol(beta_gr))
  } else {
    colMeans(beta_gr[idx, , drop = FALSE], na.rm = TRUE)
  }
})

block_beta_mat <- t(block_beta_mat)  # blocks in rows
rownames(block_beta_mat) <- paste0("PMD_", seq_len(nrow(block_beta_mat)))
colnames(block_beta_mat) <- colnames(beta_gr)

# --- 3) Define erosions per sample ---
# eroded if:
# (a) block mean <= mean_out - sample_delta_from_out
# (b) block mean < sample_beta_cutoff
mean_out_vec <- mean_out[colnames(beta_gr)]  # align order

eroded_bool <- sweep(block_beta_mat, 2, mean_out_vec, FUN = "<=") &
  (block_beta_mat < pmd_params$sample_beta_cutoff) &
  (sweep(block_beta_mat, 2, mean_out_vec, FUN = "-") <= -pmd_params$sample_delta_from_out)

# --- 4) Compute PMD load per sample ---
block_widths <- width(hypo_blocks_gr)

pmd_load_bp <- colSums(eroded_bool * block_widths, na.rm = TRUE)
n_eroded_blocks <- colSums(eroded_bool, na.rm = TRUE)

# --- 5) Join phenotype + save ---
pd_join <- obj$pd %>% mutate(Sample = Basename)

pmd_load_df <- data.frame(
  Sample = names(pmd_load_bp),
  PMD_load_bp = as.numeric(pmd_load_bp),
  PMD_load_Mb = as.numeric(pmd_load_bp) / 1e6,
  n_eroded_blocks = as.numeric(n_eroded_blocks),
  stringsAsFactors = FALSE
) %>%
  left_join(pd_join, by = "Sample") %>%
  arrange(desc(PMD_load_bp))

write.csv(
  pmd_load_df,
  file.path(results_dir, "PMD_load_per_sample.csv"),
  row.names = FALSE
)

head(pmd_load_df)

########################
## 14g. Plots (publication-quality)
########################

# 1) Cumulative PMD erosion plot (samples sorted by PMD load)
library(ggplot2)

p_cum <- ggplot(
  pmd_load_df,
  aes(
    x = reorder(Sample, PMD_load_bp),
    y = PMD_load_Mb,
    fill = Sample_Group
  )
) +
  geom_col(width = 0.9) +
  coord_flip() +
  labs(
    title = "Cumulative hypo-PMD erosion per sample",
    x = NULL,
    y = "Hypo-PMD coverage (Mb)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(),
    legend.title = element_blank()
  )

ggsave(file.path(fig_dir, "PMD_cumulative_erosion_per_sample.pdf"),
       p_cum, width = 7, height = 8)

# 2) Heatmap of inside vs outside PMD methylation (delta)
library(dplyr)
library(ggplot2)
library(pheatmap)

# --- Ensure pmd_sample_stats has needed annotation columns ---
# (Assumes pmd_sample_stats already left_join'ed with pd/pd2 earlier)
pmd_sample_stats <- pmd_sample_stats %>%
  mutate(
    # if Sample_Group exists, derive Case/Control from it
    case_status = ifelse(Sample_Group %in% c("1_Tumor", "2_OPL"),
                         "Case", "Control")
  )

# --- 2) Heatmap of inside vs outside PMD methylation (delta) ---
delta_mat <- matrix(
  pmd_sample_stats$delta_in_minus_out,
  nrow = 1,
  dimnames = list("Δβ inside−outside", pmd_sample_stats$Sample)
)

ann_col <- data.frame(
  Group      = pmd_sample_stats$Sample_Group,
  CaseStatus = pmd_sample_stats$case_status,
  row.names  = pmd_sample_stats$Sample,
  stringsAsFactors = FALSE
)

# optional: sort samples by delta for nicer viz
ord <- order(pmd_sample_stats$delta_in_minus_out, decreasing = TRUE, na.last = TRUE)
delta_mat <- delta_mat[, ord, drop = FALSE]
ann_col   <- ann_col[ord, , drop = FALSE]

pheatmap(
  delta_mat,
  annotation_col = ann_col,
  cluster_cols = TRUE,
  cluster_rows = FALSE,
  show_colnames = FALSE,
  fontsize = 10,
  main = "Δβ (inside hypo-PMD − outside) per sample"
)

ggsave(
  file.path(fig_dir, "PMD_inside_minus_outside_heatmap.pdf"),
  width = 10, height = 2.5
)

# 3) Boxplot of PMD sizes per group (block-level property, shown by case/control calling)
library(dplyr)
library(ggplot2)
library(GenomicRanges)

# --- Build hypo_blocks_df from GRanges ---
hypo_blocks_df <- as.data.frame(hypo_blocks_gr) %>%
  transmute(
    chr   = as.character(seqnames),
    start = start,
    end   = end,
    width_bp = width,
    width_kb = width / 1e3,
    width_mb = width / 1e6
  )

# (optional) drop tiny blocks if you want
if (!is.null(pmd_params$min_block_kb)) {
  hypo_blocks_df <- hypo_blocks_df %>% filter(width_kb >= pmd_params$min_block_kb)
}

# --- If you want a "case/control" color anyway, use the *sample-level* case_status ---
# Blocks are global, so a standard way is to show one distribution per sample group
# by expanding blocks per sample and joining pmd_load_df/pd.
# If you just want overall block size distribution, skip this expansion.

# Simple overall boxplot / violin of block sizes:
p_size_overall <- ggplot(hypo_blocks_df, aes(x = "hypo-PMD blocks", y = width_kb)) +
  geom_violin(fill="grey90", color="grey30", trim=FALSE) +
  geom_boxplot(width=0.15, outlier.size=0.5, fill="white") +
  labs(title="Hypo-PMD block size distribution",
       x=NULL, y="Block width (kb)") +
  theme_bw(base_size=14) +
  theme(panel.grid.major.x = element_blank())

ggsave(file.path(fig_dir, "PMD_block_sizes_overall.pdf"),
       p_size_overall, width=4, height=5)

# --- OPTIONAL: block sizes *shown per sample case/control* ---
# This repeats the same block widths for each sample,
# but lets you compare how "cases vs controls" distribute PMD sizes (purely a display choice).

# --- OPTIONAL: block sizes *shown per sample case/control* ---

if (!exists("hypo_blocks_df")) {
  stopifnot(exists("hypo_blocks_gr"))
  hypo_blocks_df <- as.data.frame(hypo_blocks_gr)
}

# add widths safely (recompute from start/end)
hypo_blocks_df <- hypo_blocks_df %>%
  dplyr::mutate(
    width_bp = as.numeric(.data$end - .data$start + 1),
    width_kb = width_bp / 1e3
  )

block_group_df <- tidyr::crossing(
  hypo_blocks_df,
  Sample = pmd_load_df$Sample
) %>%
  dplyr::left_join(
    pmd_load_df %>% dplyr::select(Sample, Sample_Group),
    by = "Sample"
  ) %>%
  dplyr::mutate(
    case_status = dplyr::if_else(
      Sample_Group %in% c("1_Tumor", "2_OPL"),
      "Case", "Control"
    )
  )

p_size_by_case <- ggplot(
  block_group_df,
  aes(x = case_status, y = width_kb, fill = case_status)
) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.2, outlier.size = 0.3, alpha = 0.9) +
  labs(
    title = "Hypo-PMD block sizes (displayed by case/control)",
    x = NULL, y = "Block width (kb)"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

ggsave(
  file.path(fig_dir, "PMD_block_sizes_by_case_control.pdf"),
  p_size_by_case, width = 5.5, height = 5
)

p_size_by_group <- ggplot(
  block_group_df,
  aes(x = Sample_Group, y = width_kb, fill = Sample_Group)
) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.2, outlier.size = 0.3, alpha = 0.9) +
  labs(
    title = "Hypo-PMD block sizes stratified by sample group",
    x = NULL, y = "Block width (kb)"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

ggsave(
  file.path(fig_dir, "PMD_block_sizes_by_group.pdf"),
  p_size_by_group, width = 6, height = 5
)

########################
## 14h. Save objects
########################
saveRDS(list(
  hypo_blocks_gr = hypo_blocks_gr,
  pmd_sample_stats = pmd_sample_stats,
  pmd_load_df = pmd_load_df,
  block_beta_mat = block_beta_mat,
  eroded_bool = eroded_bool
), file.path(results_dir, "PMD_erosion_objects.rds"))

message("PMD erosion module complete.")





















############################
## Publication-quality PMD figures
## Inputs:
##   hypo_PMD_blocks.csv
##   PMD_inside_outside_stats_per_sample.csv
##   PMD_load_per_sample.csv
############################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(scales)
})

# --- Paths (edit if needed) ---
hypo_blocks_csv <- file.path(results_dir, "hypo_PMD_blocks.csv")
inside_out_csv  <- file.path(results_dir, "PMD_inside_outside_stats_per_sample.csv")
load_csv        <- file.path(results_dir, "PMD_load_per_sample.csv")

# --- 0) Load with base R to avoid readr sysdata issues ---
hypo_blocks_df <- read.csv(hypo_blocks_csv, stringsAsFactors = FALSE)
pmd_sample_stats <- read.csv(inside_out_csv, stringsAsFactors = FALSE)
pmd_load_df <- read.csv(load_csv, stringsAsFactors = FALSE)

# --- 1) Standardize hypo_blocks_df column names/types ---
# expected columns from your file: seqnames, start, end, width, value, ...
if (!"chr" %in% names(hypo_blocks_df)) {
  if ("seqnames" %in% names(hypo_blocks_df)) {
    hypo_blocks_df <- hypo_blocks_df %>% rename(chr = seqnames)
  }
}

# force numeric start/end/width/value
hypo_blocks_df <- hypo_blocks_df %>%
  mutate(
    start = as.numeric(start),
    end   = as.numeric(end),
    width = if ("width" %in% names(.)) as.numeric(width) else (end - start + 1),
    value = if ("value" %in% names(.)) as.numeric(value) else NA_real_
  )

# add widths in bp/kb/Mb
hypo_blocks_df <- hypo_blocks_df %>%
  mutate(
    width_bp = width,
    width_kb = width_bp / 1e3,
    width_Mb = width_bp / 1e6
  )

# --- 2) Standardize per-sample tables + derive case_status ---
# Use Sample_Group if present; otherwise fall back to existing case_status or Control/Case labels.

derive_case_status <- function(group_vec){
  # OPL and Tumor are both Case (disease); Normal is Control
  out <- ifelse(group_vec %in% c("1_Tumor", "2_OPL"), "Case", "Control")
  factor(out, levels=c("Control","Case"))
}

if (!"case_status" %in% names(pmd_load_df)) {
  if ("Sample_Group" %in% names(pmd_load_df)) {
    pmd_load_df$case_status <- derive_case_status(pmd_load_df$Sample_Group)
  } else {
    pmd_load_df$case_status <- factor("Case") # last resort
  }
} else {
  pmd_load_df$case_status <- factor(pmd_load_df$case_status, levels=c("Control","Case"))
}

if (!"case_status" %in% names(pmd_sample_stats)) {
  if ("Sample_Group" %in% names(pmd_sample_stats)) {
    pmd_sample_stats$case_status <- derive_case_status(pmd_sample_stats$Sample_Group)
  } else if ("Sample" %in% names(pmd_sample_stats) && "Sample_Group" %in% names(pmd_load_df)) {
    # join group from load table if needed
    pmd_sample_stats <- pmd_sample_stats %>%
      left_join(pmd_load_df %>% select(Sample, Sample_Group, case_status),
                by="Sample")
  }
}

# Ensure Sample is character everywhere
pmd_load_df$Sample <- as.character(pmd_load_df$Sample)
pmd_sample_stats$Sample <- as.character(pmd_sample_stats$Sample)

# Match ordering between tables
pmd_load_df <- pmd_load_df %>% arrange(desc(PMD_load_bp))
pmd_sample_stats <- pmd_sample_stats %>%
  left_join(pmd_load_df %>% select(Sample, PMD_load_bp, Sample_Group, case_status),
            by="Sample")

############################
## FIGURE 1. Genome-wide PMD erosion map (Nature-ready)
############################

############################
## FIGURE 1. Genome-wide PMD erosion map (Nature-ready, NO Cairo)
############################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
})

# ==============================================================================
# 1. DATA CLEANING & PREPARATION (autosomals only)
# ==============================================================================
target_chromosomes <- paste0("chr", 1:22)

clean_df <- hypo_blocks_df %>%
  mutate(
    chr   = as.character(chr),
    chr   = ifelse(!grepl("^chr", chr), paste0("chr", chr), chr),
    start = as.numeric(start),
    end   = as.numeric(end),
    value = as.numeric(value)
  ) %>%
  filter(chr %in% target_chromosomes) %>%
  mutate(chr = factor(chr, levels = target_chromosomes)) %>%
  arrange(chr, start)

# ==============================================================================
# 2. CALCULATE CHR OFFSETS + CUMULATIVE COORDS (numeric safe)
# ==============================================================================
chr_stats <- clean_df %>%
  group_by(chr) %>%
  summarise(chr_len = max(end, na.rm = TRUE), .groups = "drop") %>%
  arrange(chr) %>%
  mutate(
    chr_len = as.numeric(chr_len),
    tot     = cumsum(lag(chr_len, default = 0)),
    x_center = tot + chr_len/2
  )

plot_df <- clean_df %>%
  left_join(chr_stats[, c("chr","tot")], by="chr") %>%
  mutate(
    cum_start = start + tot,
    cum_end   = end   + tot
  )

# ==============================================================================
# 3. AXIS + ZEBRA STRIPES
# ==============================================================================
axis_df <- chr_stats %>%
  mutate(
    is_odd = as.integer(chr) %% 2 == 1,
    xmin   = tot,
    xmax   = tot + chr_len
  )

x_breaks <- axis_df$x_center
x_labels <- gsub("chr", "", axis_df$chr)

# ==============================================================================
# 4. LABELS + LIMITS
# ==============================================================================
# Use math expressions but avoid unicode minus issues in themes
y_label_math <- expression(Delta*beta~"(block mean)")
legend_title <- expression(Delta*beta)

y_range  <- range(plot_df$value, na.rm = TRUE)
y_buffer <- diff(y_range) * 0.06
y_limits <- c(y_range[1]-y_buffer, y_range[2]+y_buffer)

# ==============================================================================
# 5. PLOT (vector points; no rasterization needed)
# ==============================================================================
p_manhattan <- ggplot() +

  # zebra chromosome bands
  geom_rect(
    data = axis_df,
    aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf, fill=is_odd),
    alpha = 0.10,
    show.legend = FALSE
  ) +
  scale_fill_manual(values=c("TRUE"="grey85", "FALSE"="white")) +

  # baseline at 0
  geom_hline(yintercept=0, linewidth=0.25, color="grey35") +

  # blocks
  geom_point(
    data = plot_df,
    aes(x=cum_start, y=value, color=value),
    size = 1.1,
    alpha = 0.9
  ) +

  # colorblind-safe diverging scale centered at 0
  scale_color_gradient2(
    low  = "#2166AC",
    mid  = "grey92",
    high = "#B2182B",
    midpoint = 0,
    limits = y_limits,
    oob = squish,
    name = legend_title
  ) +

  scale_x_continuous(
    breaks = x_breaks,
    labels = x_labels,
    expand = expansion(mult=c(0, 0.002))
  ) +
  coord_cartesian(ylim = y_limits, clip="off") +

  labs(
    title = "Genome-wide hypo-PMD erosion map",
    x = "Chromosome",
    y = y_label_math
  ) +

  theme_classic(base_size = 8) +
  theme(
    text = element_text(family="sans"),

    axis.text.x  = element_text(size=7, color="black"),
    axis.text.y  = element_text(size=7, color="black"),
    axis.title.x = element_text(size=8, margin=margin(t=4)),
    axis.title.y = element_text(size=8, margin=margin(r=4)),
    axis.ticks.x = element_blank(),

    panel.border = element_rect(color="black", fill=NA, linewidth=0.35),

    legend.position = "right",
    legend.title = element_text(size=8),
    legend.text  = element_text(size=7),
    legend.key.height = unit(0.35, "cm"),

    plot.title = element_text(size=9, face="plain", hjust=0)
  )

print(p_manhattan)

# ==============================================================================
# 6. EXPORT (NO Cairo)
# ==============================================================================
# Nature-ish sizing in mm:
# width = 180 mm (full single-column width), height = 70 mm
ggsave(
  filename = file.path(fig_dir, "Fig1_PMD_manhattan_nature.pdf"),
  plot = p_manhattan,
  width = 300/25.4,
  height = 70/25.4,
  device = grDevices::pdf,
  dpi = 600,
  useDingbats = FALSE  # avoids font weirdness in some journal pipelines
)

ggsave(
  filename = file.path(fig_dir, "Fig1_PMD_manhattan_nature.png"),
  plot = p_manhattan,
  width = 180/25.4,
  height = 70/25.4,
  dpi = 600,
  bg = "white"
)

############################
## FIGURE 2. PMD load by 3 biological groups
############################

library(dplyr)
library(ggplot2)
library(scales)
library(stringr)

# 0) Make sure Sample_Group is clean and ordered
pmd_load_df2 <- pmd_load_df %>%
  mutate(
    Sample_Group = as.character(Sample_Group),
    # standardize just in case
    Sample_Group = case_when(
      str_detect(Sample_Group, regex("tumor", ignore_case=TRUE))  ~ "1_Tumor",
      str_detect(Sample_Group, regex("opl", ignore_case=TRUE))    ~ "2_OPL",
      str_detect(Sample_Group, regex("normal|control", ignore_case=TRUE)) ~ "3_Normal",
      TRUE ~ Sample_Group
    ),
    Sample_Group = factor(
      Sample_Group,
      levels = c("3_Normal", "2_OPL", "1_Tumor")  # biological progression order
    )
  )

# Colorblind-safe, Nature-ish palette
pal_groups <- c(
  "3_Normal" = "#3C5488FF",
  "2_OPL"    = "#00A087FF",
  "1_Tumor"  = "#E64B35FF"
)

# 1) Main text figure: distribution by group
p_fig2 <- ggplot(pmd_load_df2,
                 aes(x = Sample_Group, y = PMD_load_Mb, fill = Sample_Group)) +
  geom_violin(trim = FALSE, alpha = 0.55, linewidth = 0.3, color="grey20") +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.95, linewidth = 0.35) +
  geom_jitter(width = 0.10, size = 1.6, alpha = 0.75) +
  scale_fill_manual(values = pal_groups) +
  scale_y_continuous(
    name = "Hypo-PMD coverage (Mb)",
    expand = expansion(mult = c(0.02, 0.06))
  ) +
  labs(
    title = "Cumulative hypo-PMD erosion by sample group",
    x = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face="bold"),
    axis.text.x = element_text(face="bold")
  )

ggsave(file.path(fig_dir, "Fig2_PMD_load_by_group.pdf"),
       p_fig2, width = 5.0, height = 4.6, useDingbats = FALSE)

print(p_fig2)

# 2) Optional companion panel (good for supplement):
#    ranked PMD load, colored by group, no sample labels
p_rank <- pmd_load_df2 %>%
  arrange(desc(PMD_load_Mb)) %>%
  mutate(rank = row_number())

p_fig2b <- ggplot(p_rank, aes(x = rank, y = PMD_load_Mb, color = Sample_Group)) +
  geom_point(size = 2.0, alpha = 0.9) +
  scale_color_manual(values = pal_groups) +
  labs(
    title = "Ranked hypo-PMD load across samples",
    x = "Sample rank (highest to lowest)",
    y = "Hypo-PMD coverage (Mb)",
    color = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face="bold"))

ggsave(file.path(fig_dir, "Fig2b_PMD_load_ranked.pdf"),
       p_fig2b, width = 6.2, height = 4.4, useDingbats = FALSE)

print(p_fig2b)

############################
## FIGURE 3. Inside vs outside PMD methylation per sample
## (delta heatmap + boxplot)  -- 3-group aware
############################

library(dplyr)
library(stringr)
library(ggplot2)

# ------------------------------------------------------------------------------
# 0) Harmonize columns: collapse x/y duplicates into one clean grouping
# ------------------------------------------------------------------------------

pmd_sample_stats2 <- pmd_sample_stats %>%
  mutate(
    Sample = as.character(Sample),

    # Prefer Sample_Group.x, otherwise Sample_Group.y
    Sample_Group = dplyr::coalesce(as.character(Sample_Group.x),
                                   as.character(Sample_Group.y)),

    # Optional: keep a single case_status if you need it elsewhere
    case_status = dplyr::coalesce(as.character(case_status.x),
                                  as.character(case_status.y))
  ) %>%
  # standardize naming to the 3 biological groups
  mutate(
    Sample_Group = case_when(
      str_detect(Sample_Group, regex("tumor", ignore_case = TRUE))  ~ "1_Tumor",
      str_detect(Sample_Group, regex("opl",   ignore_case = TRUE))  ~ "2_OPL",
      str_detect(Sample_Group, regex("normal|control", ignore_case = TRUE)) ~ "3_Normal",
      TRUE ~ Sample_Group
    ),
    Sample_Group = factor(Sample_Group,
                          levels = c("3_Normal", "2_OPL", "1_Tumor")) # progression order
  )

# sanity check
table(pmd_sample_stats2$Sample_Group, useNA="ifany")

# ------------------------------------------------------------------------------
# 1) Order samples for heatmap (by group then delta)
# ------------------------------------------------------------------------------

pmd_sample_stats2 <- pmd_sample_stats2 %>%
  arrange(Sample_Group, delta_in_minus_out)

sample_order <- pmd_sample_stats2$Sample

# ------------------------------------------------------------------------------
# 2) Heatmap of delta (inside - outside)
# ------------------------------------------------------------------------------

delta_mat <- matrix(
  pmd_sample_stats2$delta_in_minus_out,
  nrow = 1,
  dimnames = list("Delta_beta_inside_minus_outside", sample_order)
)

# annotation for columns (3-group only)
ann_col <- data.frame(
  Sample_Group = pmd_sample_stats2$Sample_Group,
  row.names = sample_order
)

# Use pheatmap if installed, otherwise ggplot fallback
if (requireNamespace("pheatmap", quietly = TRUE)) {
  library(pheatmap)

  pdf(file.path(fig_dir, "Fig3_PMD_delta_in_out_heatmap.pdf"),
      width = 10, height = 2.3, useDingbats = FALSE)

  pheatmap(
    delta_mat,
    cluster_cols = FALSE,       # keep biological ordering
    cluster_rows = FALSE,
    annotation_col = ann_col,
    show_colnames = FALSE,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
    breaks = seq(min(delta_mat, na.rm=TRUE),
                 max(delta_mat, na.rm=TRUE), length.out=101),
    main = "Inside vs outside hypo-PMD methylation"
  )

  dev.off()

} else {

  delta_df <- data.frame(
    Sample = factor(sample_order, levels = sample_order),
    delta  = pmd_sample_stats2$delta_in_minus_out,
    Sample_Group = pmd_sample_stats2$Sample_Group
  )

  p_delta_tile <- ggplot(delta_df, aes(x = Sample, y = 1, fill = delta)) +
    geom_tile() +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
      name = "Delta_beta"
    ) +
    facet_grid(~Sample_Group, scales="free_x", space="free_x") +
    theme_void(base_size = 12) +
    theme(
      strip.text = element_text(size=10, face="bold"),
      legend.position = "right"
    ) +
    labs(title = "Inside vs outside hypo-PMD methylation (Δβ)")

  ggsave(file.path(fig_dir, "Fig3_PMD_delta_in_out_heatmap_fallback.pdf"),
         p_delta_tile, width = 10, height = 2.3, dpi=300, useDingbats=FALSE)
}

# ------------------------------------------------------------------------------
# 3) Companion violin+boxplot of delta across the 3 groups
# ------------------------------------------------------------------------------

p_delta_box <- ggplot(
  pmd_sample_stats2,
  aes(x = Sample_Group, y = delta_in_minus_out, fill = Sample_Group)
) +
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
  geom_boxplot(width = 0.18, outlier.size = 0.7, alpha = 0.9) +
  scale_fill_manual(values = c(
    "3_Normal" = "#4DAF4A",
    "2_OPL"    = "#FFB000",
    "1_Tumor"  = "#E41A1C"
  )) +
  labs(
    title = "Δβ (inside − outside) hypo-PMD methylation by group",
    x = NULL,
    y = "Δβ inside − outside"
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

ggsave(file.path(fig_dir, "Fig3b_PMD_delta_in_out_boxplot.pdf"),
       p_delta_box, width = 4.8, height = 4.5, dpi=300, useDingbats=FALSE)

print(p_delta_box)

############################
## FIGURE 4. PMD load by biological group (Normal / OPL / Tumor)
############################

library(dplyr)
library(ggplot2)
library(stringr)

# ------------------------------------------------------------------------------
# 0) Clean Sample_Group (your dataset contains only "Sample_Group")
# ------------------------------------------------------------------------------

pmd_load_df2 <- pmd_load_df %>%
  mutate(
    Sample = as.character(Sample),

    # Clean Sample_Group text
    Sample_Group = as.character(Sample_Group),
    Sample_Group = case_when(
      str_detect(Sample_Group, regex("tumor",  ignore_case=TRUE)) ~ "1_Tumor",
      str_detect(Sample_Group, regex("opl",    ignore_case=TRUE)) ~ "2_OPL",
      str_detect(Sample_Group, regex("normal|control", ignore_case=TRUE)) ~ "3_Normal",
      TRUE ~ Sample_Group
    ),

    # Ordered factor (Normal → OPL → Tumor)
    Sample_Group = factor(
      Sample_Group,
      levels = c("3_Normal", "2_OPL", "1_Tumor")
    )
  )

# Check the recoding succeeded
print(table(pmd_load_df2$Sample_Group, useNA = "ifany"))


# ------------------------------------------------------------------------------
# 1) Plot: PMD load by biological group
# ------------------------------------------------------------------------------

colors_groups <- c(
  "3_Normal" = "#3C5488FF",  # blue
  "2_OPL"    = "#00A087FF",  # green
  "1_Tumor"  = "#E64B35FF"   # red
)

p_load_group <- ggplot(
  pmd_load_df2,
  aes(x = Sample_Group, y = PMD_load_Mb, fill = Sample_Group)
) +
  geom_violin(trim = FALSE, alpha = 0.55, color = NA) +
  geom_boxplot(width = 0.18, outlier.shape = NA,
               alpha = 0.95, linewidth = 0.35) +
  geom_jitter(width = 0.15, size = 1.4, alpha = 0.75) +
  scale_fill_manual(values = colors_groups) +
  labs(
    title = "Hypo-PMD load across biological groups",
    x = NULL,
    y = "Hypo-PMD load (Mb)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(fig_dir, "Fig4_PMD_load_by_group.pdf"),
  p_load_group,
  width = 6, height = 5, dpi = 300, useDingbats = FALSE
)

print(p_load_group)


# ------------------------------------------------------------------------------
# 2) Optional statistical testing
# ------------------------------------------------------------------------------

# Global differences
kw <- kruskal.test(PMD_load_Mb ~ Sample_Group, data = pmd_load_df2)
print(kw)

# Pairwise Wilcoxon tests
pairwise <- pairwise.wilcox.test(
  pmd_load_df2$PMD_load_Mb,
  pmd_load_df2$Sample_Group,
  p.adjust.method = "BH"
)
print(pairwise)

# ------------------------------------------------------------------------------
# 2) Optional: Statistical testing (recommended for paper)
# ------------------------------------------------------------------------------

# Global difference
kw <- kruskal.test(PMD_load_Mb ~ Sample_Group, data = pmd_load_df2)
print(kw)

# Pairwise comparisons
pairwise <- pairwise.wilcox.test(
  pmd_load_df2$PMD_load_Mb,
  pmd_load_df2$Sample_Group,
  p.adjust.method = "BH"
)
print(pairwise)

############################
## FIGURE 5. Distribution of hypo-PMD block sizes
############################

library(dplyr)
library(ggplot2)
library(scales)

# ------------------------------------------------------------------------------
# 0) Ensure hypo_blocks_df exists and has numeric width
# ------------------------------------------------------------------------------

if (!exists("hypo_blocks_df")) {
  stopifnot(exists("hypo_blocks_gr"))
  hypo_blocks_df <- as.data.frame(hypo_blocks_gr)
}

# Some GRanges -> data.frame exports call width something else or it can get masked.
# So we define width_bp explicitly from the GRanges object if needed.
if (!("width" %in% names(hypo_blocks_df))) {
  hypo_blocks_df$width <- IRanges::width(hypo_blocks_gr)
}

# Force numeric width (protect against width being a function / generic)
hypo_blocks_df <- hypo_blocks_df %>%
  mutate(
    width_bp = as.numeric(width),
    width_kb = width_bp / 1e3
  ) %>%
  filter(!is.na(width_kb), width_kb > 0)

# ------------------------------------------------------------------------------
# 1) Plot: log10 distribution with violin + box + points
# ------------------------------------------------------------------------------

p_block_sizes <- ggplot(hypo_blocks_df, aes(x = "", y = width_kb)) +
  geom_violin(
    fill = "grey85", color = NA, trim = FALSE, alpha = 0.9
  ) +
  geom_boxplot(
    width = 0.12, outlier.shape = NA, alpha = 0.95, linewidth = 0.35
  ) +
  geom_jitter(
    width = 0.08, size = 1.0, alpha = 0.35
  ) +
  scale_y_log10(
    labels = scales::comma_format(accuracy = 1)
  ) +
  labs(
    title = "Hypo-PMD block size distribution",
    x = NULL,
    y = "Block width (kb, log10 scale)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title   = element_text(face = "bold")
  )

ggsave(
  file.path(fig_dir, "Fig5_block_size_distribution.pdf"),
  p_block_sizes,
  width = 4.3, height = 5,
  dpi = 300, useDingbats = FALSE
)

print(p_block_sizes)

message("Figure 5 saved to: ", fig_dir)









############################
## Chromosome-level hypo-PMD "breakdown" summary + plot
## Uses hypo_blocks_df (or reads hypo_PMD_blocks.csv)
## Breakdown = how much of each chromosome is covered by hypo-PMD blocks
## (optionally severity-weighted by |Δβ|)
############################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(scales)
})

# ---------------------------
# 0) Load / standardize blocks df
# ---------------------------
if (!exists("hypo_blocks_df")) {
  hypo_blocks_df <- read_csv(hypo_blocks_csv, show_col_types = FALSE)
}

safe_rename <- function(df, new, old) {
  old <- intersect(old, names(df))
  if (length(old) == 1) {
    df <- dplyr::rename(df, !!new := !!sym(old))
  }
  df
}

# Standardize key column names
hypo_blocks_df <- hypo_blocks_df %>%
  safe_rename("chr",      c("chr","seqnames","chromosome")) %>%
  safe_rename("start",    c("start","start_position","bp_start")) %>%
  safe_rename("end",      c("end","end_position","bp_end")) %>%
  safe_rename("width_bp", c("width_bp","width","len","widths")) %>%
  safe_rename("delta_beta", c("value","mean_diff","delta_beta","mean_delta")) %>%
  mutate(
    chr        = as.character(chr),
    chr        = ifelse(!stringr::str_detect(chr, "^chr"), paste0("chr", chr), chr),
    start      = as.numeric(start),
    end        = as.numeric(end),
    width_bp   = as.numeric(width_bp),
    delta_beta = as.numeric(delta_beta)
  )

# Keep autosomes only for clean ordering (optional)
autosomes <- paste0("chr", 1:22)
hypo_blocks_df2 <- hypo_blocks_df %>%
  filter(chr %in% autosomes) %>%
  mutate(chr = factor(chr, levels = autosomes)) %>%
  # drop tiny 1bp artifacts if present
  filter(!is.na(width_bp), width_bp >= 1000)

# ---------------------------
# 1) Chromosome summary
# ---------------------------
chr_breakdown <- hypo_blocks_df2 %>%
  group_by(chr) %>%
  summarise(
    n_blocks            = n(),
    total_hypo_PMD_bp   = sum(width_bp, na.rm = TRUE),
    total_hypo_PMD_Mb   = total_hypo_PMD_bp / 1e6,
    mean_delta_beta     = mean(delta_beta, na.rm = TRUE),
    # severity-weighted coverage (more negative Δβ counts more)
    severity_weighted_bp = sum(width_bp * abs(delta_beta), na.rm = TRUE),
    severity_weighted_Mb = severity_weighted_bp / 1e6,
    .groups = "drop"
  ) %>%
  arrange(desc(total_hypo_PMD_bp))

write_csv(chr_breakdown, file.path(results_dir, "PMD_breakdown_by_chromosome.csv"))

# ---------------------------
# 2) Plot A: total hypo-PMD coverage per chromosome
# ---------------------------
p_chr_cov <- ggplot(chr_breakdown,
                    aes(x = reorder(chr, total_hypo_PMD_bp),
                        y = total_hypo_PMD_Mb)) +
  geom_col(width = 0.85, fill = "grey25") +
  coord_flip() +
  labs(
    title = "Chromosome-level hypo-PMD breakdown",
    subtitle = "Total genomic coverage of hypo-PMD blocks",
    x = NULL,
    y = "Hypo-PMD coverage (Mb)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  )

ggsave(file.path(fig_dir, "Fig_chr_hypoPMD_coverage_Mb.pdf"),
       p_chr_cov, width=6.5, height=6, useDingbats=FALSE)

# ---------------------------
# 3) Plot B: severity-weighted breakdown (optional but nice)
# ---------------------------
p_chr_sev <- ggplot(chr_breakdown,
                    aes(x = reorder(chr, severity_weighted_bp),
                        y = severity_weighted_Mb)) +
  geom_col(width = 0.85, fill = "#1f78b4") +
  coord_flip() +
  labs(
    title = "Chromosome-level hypo-PMD breakdown (severity-weighted)",
    subtitle = "Coverage weighted by |Δβ| of each block",
    x = NULL,
    y = "Severity-weighted hypo-PMD (Mb·|Δβ|)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  )

ggsave(file.path(fig_dir, "Fig_chr_hypoPMD_severity_weighted.pdf"),
       p_chr_sev, width=6.5, height=6, useDingbats=FALSE)

# ---------------------------
# 4) Plot C: two-panel lollipop (coverage + mean Δβ)
# ---------------------------
chr_breakdown_long <- chr_breakdown %>%
  mutate(chr_lab = gsub("chr","", as.character(chr))) %>%
  arrange(total_hypo_PMD_bp) %>%
  mutate(chr_lab = factor(chr_lab, levels = chr_lab))

p_lolli <- ggplot(chr_breakdown_long, aes(y = chr_lab)) +
  geom_segment(aes(x = 0, xend = total_hypo_PMD_Mb, yend = chr_lab),
               linewidth=0.8, color="grey50") +
  geom_point(aes(x = total_hypo_PMD_Mb), size=2.2, color="black") +
  geom_point(aes(x = severity_weighted_Mb), size=2.2, color="#1f78b4") +
  labs(
    title = "Chromosome breakdown (coverage vs severity)",
    x = "Mb (black = total coverage, blue = severity-weighted)",
    y = "Chromosome"
  ) +
  theme_classic(base_size=14)

ggsave(file.path(fig_dir, "Fig_chr_breakdown_lollipop.pdf"),
       p_lolli, width=7, height=6, useDingbats=FALSE)

message("Chromosome breakdown table + figures saved to: ", fig_dir)















############################
## Per-sample PMD breakdown per chromosome
## Requires objects from your pipeline:
##   - probe_gr, beta_gr
##   - hypo_blocks_gr
##   - pmd_params
##   - pmd_load_df (for ordering / groups)
############################

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

stopifnot(exists("probe_gr"),
          exists("beta_gr"),
          exists("hypo_blocks_gr"),
          exists("mean_out"),
          exists("pmd_params"))

# autosomes only (recommended for clean interpretation)
autosomes <- paste0("chr", 1:22)

# Ensure probe GRanges has seqnames
probe_chr <- as.character(seqnames(probe_gr))

# Ensure blocks have seqnames
block_chr <- as.character(seqnames(hypo_blocks_gr))

# Keep autosomal blocks
hypo_blocks_auto <- hypo_blocks_gr[block_chr %in% autosomes]

# --- Helper: compute per-sample erosion load for a set of blocks ---
compute_load_for_blocks <- function(blocks_gr, probe_gr, beta_gr, mean_out_vec, pmd_params){

  # overlap probes with blocks
  ol <- findOverlaps(probe_gr, blocks_gr)
  if (length(ol)==0) {
    return(rep(0, ncol(beta_gr)))
  }

  block_probe_list <- split(queryHits(ol), subjectHits(ol))
  block_probe_list <- lapply(seq_along(blocks_gr), function(i){
    block_probe_list[[as.character(i)]] %||% integer(0)
  })

  # block x sample mean beta
  block_beta_mat <- sapply(block_probe_list, function(idx){
    if (length(idx)==0) rep(NA_real_, ncol(beta_gr))
    else colMeans(beta_gr[idx, , drop=FALSE], na.rm=TRUE)
  })
  block_beta_mat <- t(block_beta_mat)

  # erosion logic (same as your 14f)
  mean_out_vec <- mean_out_vec[colnames(beta_gr)]
  eroded_bool <- sweep(block_beta_mat, 2, mean_out_vec, FUN = "<=") &
    (block_beta_mat < pmd_params$sample_beta_cutoff) &
    (sweep(block_beta_mat, 2, mean_out_vec, FUN = "-") <= -pmd_params$sample_delta_from_out)

  widths_bp <- width(blocks_gr)
  load_bp <- colSums(eroded_bool * widths_bp, na.rm=TRUE)

  return(load_bp)
}

# --- Main loop: chromosome by chromosome ---
mean_out_vec <- mean_out[colnames(beta_gr)]  # align names/order

chr_sample_load <- lapply(autosomes, function(chr_i){

  blocks_chr <- hypo_blocks_auto[as.character(seqnames(hypo_blocks_auto)) == chr_i]

  load_bp_chr <- compute_load_for_blocks(
    blocks_chr, probe_gr, beta_gr, mean_out_vec, pmd_params
  )

  data.frame(
    Sample = colnames(beta_gr),
    chr = chr_i,
    PMD_load_bp_chr = as.numeric(load_bp_chr),
    PMD_load_Mb_chr = as.numeric(load_bp_chr)/1e6,
    stringsAsFactors = FALSE
  )
})

chr_sample_load_df <- bind_rows(chr_sample_load)

# join metadata + total PMD load for ordering
chr_sample_load_df <- chr_sample_load_df %>%
  left_join(
    pmd_load_df %>% select(Sample, PMD_load_Mb, Sample_Group),
    by="Sample"
  ) %>%
  mutate(
    Sample_Group = case_when(
      str_detect(Sample_Group, regex("tumor", ignore_case=TRUE))  ~ "1_Tumor",
      str_detect(Sample_Group, regex("opl", ignore_case=TRUE))    ~ "2_OPL",
      str_detect(Sample_Group, regex("normal|control", ignore_case=TRUE)) ~ "3_Normal",
      TRUE ~ Sample_Group
    ),
    Sample_Group = factor(Sample_Group,
                          levels=c("3_Normal","2_OPL","1_Tumor"))
  )

write.csv(chr_sample_load_df,
          file.path(results_dir, "PMD_load_per_sample_per_chromosome.csv"),
          row.names=FALSE)

head(chr_sample_load_df)



############################
## FIGURE: Cohort-wide chromosome breakdown
############################
library(ggplot2)
library(scales)

chr_total_df <- chr_sample_load_df %>%
  group_by(chr) %>%
  summarise(
    total_Mb = sum(PMD_load_Mb_chr, na.rm=TRUE),
    mean_Mb  = mean(PMD_load_Mb_chr, na.rm=TRUE),
    n_samples_erosion = sum(PMD_load_bp_chr > 0, na.rm=TRUE),
    .groups="drop"
  ) %>%
  arrange(desc(total_Mb)) %>%
  mutate(chr = factor(chr, levels=chr))

p_chr_total <- ggplot(chr_total_df,
                      aes(x=chr, y=total_Mb)) +
  geom_col(width=0.8, fill="grey25") +
  coord_flip() +
  labs(
    title="Chromosomes with greatest hypo-PMD breakdown",
    subtitle="Total eroded hypo-PMD coverage summed across samples",
    x=NULL, y="Total hypo-PMD erosion (Mb)"
  ) +
  theme_classic(base_size=14)

ggsave(file.path(fig_dir, "Fig_chr_total_breakdown.pdf"),
       p_chr_total, width=6.5, height=6, dpi=300, useDingbats=FALSE)

p_chr_total







############################
## FIGURE: Heatmap (samples x chromosomes)
############################
library(pheatmap)

# wide matrix: samples x chr
heat_df <- chr_sample_load_df %>%
  select(Sample, chr, PMD_load_Mb_chr, PMD_load_Mb, Sample_Group) %>%
  distinct()

wide_df <- heat_df %>%
  select(Sample, chr, PMD_load_Mb_chr) %>%
  pivot_wider(names_from = chr, values_from = PMD_load_Mb_chr, values_fill = 0)

heat_mat <- as.matrix(wide_df[,-1])
rownames(heat_mat) <- wide_df$Sample

# order samples by total load (descending)
ord_samp <- heat_df %>%
  distinct(Sample, PMD_load_Mb, Sample_Group) %>%
  arrange(desc(PMD_load_Mb))

heat_mat <- heat_mat[ord_samp$Sample, , drop=FALSE]

ann_row <- data.frame(
  Group = ord_samp$Sample_Group,
  row.names = ord_samp$Sample
)

pdf(file.path(fig_dir, "Fig_chr_sample_heatmap.pdf"),
    width=9, height=7, useDingbats=FALSE)

pheatmap(
  heat_mat,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  annotation_row = ann_row,
  color = colorRampPalette(c("white","#56B1F7","#132B43"))(100),
  fontsize = 9,
  main = "Per-sample hypo-PMD erosion across chromosomes (Mb)"
)

dev.off()




############################
## Correlate chromosome erosion with sample severity rank
############################
library(dplyr)

rank_df <- pmd_load_df %>%
  arrange(desc(PMD_load_Mb)) %>%
  mutate(rank = row_number()) %>%
  select(Sample, rank)

corr_df <- chr_sample_load_df %>%
  left_join(rank_df, by="Sample") %>%
  group_by(chr) %>%
  summarise(
    spearman_rho = cor(rank, PMD_load_Mb_chr, method="spearman"),
    spearman_p   = cor.test(rank, PMD_load_Mb_chr, method="spearman")$p.value,
    .groups="drop"
  ) %>%
  mutate(p_adj = p.adjust(spearman_p, method="BH")) %>%
  arrange(p_adj)

write.csv(corr_df,
          file.path(results_dir, "PMD_chr_rank_spearman.csv"),
          row.names=FALSE)

corr_df




library(ggplot2)

corr_df <- corr_df %>%
  mutate(sig = p_adj < 0.05,
         chr = factor(chr, levels=rev(chr)))

p_corr <- ggplot(corr_df,
                 aes(y=chr, x=spearman_rho, color=sig)) +
  geom_segment(aes(x=0, xend=spearman_rho, yend=chr),
               linewidth=0.8) +
  geom_point(size=2.6) +
  scale_color_manual(values=c("TRUE"="#E64B35FF","FALSE"="grey40")) +
  labs(
    title="Do chromosome-specific erosions track with sample severity rank?",
    x="Spearman rho (rank vs chr erosion)",
    y="Chromosome",
    color=NULL
  ) +
  theme_classic(base_size=14)

ggsave(file.path(fig_dir, "Fig_chr_rank_correlation_lollipop.pdf"),
       p_corr, width=7, height=6, dpi=300, useDingbats=FALSE)

p_corr
