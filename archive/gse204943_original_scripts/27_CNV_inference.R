# ==============================================================================
# Script 27: Copy Number Inference from Methylation (conumee2) — ROBUST
# Project: OPMD Methylation Analysis (GSE204943)
# Author: Michael E. Troka
#
# IMPORTANT: Your sanity test proved CNV.fit/bin/detail/segment works.
# The error you see in the loop happens AFTER the 3x “(100%)” messages,
# i.e. most often inside CNV.write() (because the sanity test didn’t call CNV.write).
# This version wraps CNV.write() safely so a write failure can’t kill the loop.
# ==============================================================================

suppressPackageStartupMessages({
  library(minfi)
  library(conumee2)
  library(dplyr)
  library(readr)
})

# ============================================================
# 0. Paths
# ============================================================
results_dir <- "GSE204943_data/GSE204943/results"
fig_dir     <- "GSE204943_data/GSE204943/figures"
source_dir  <- file.path(results_dir, "source_data")
dir.create(source_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. Load data
# ============================================================
cat("Loading analysis objects...\n")
obj <- readRDS(file.path(results_dir, "analysis_objects.rds"))
pd  <- obj$pd

group <- factor(pd$Sample_Group, levels = c("3_Normal", "2_OPL", "1_Tumor"))
sample_ids <- basename(pd$Basename)

if (!is.null(obj$RGset.qc)) {
  cat("Using preprocessIllumina(RGset.qc) for intensities...\n")
  Mset <- preprocessIllumina(obj$RGset.qc)
  sampleNames(Mset) <- sample_ids
} else {
  cat("RGset.qc not available; using stored Mset.noob...\n")
  Mset <- obj$Mset.noob
  sampleNames(Mset) <- sample_ids
}

normal_idx <- which(group == "3_Normal")
query_idx  <- which(group != "3_Normal")

cat("\n--- Running conumee2 CNV inference ---\n")
cat("Building reference from", length(normal_idx), "Normal samples\n")
cat("Query samples:", length(query_idx), "\n")

# ============================================================
# 2. Build ONE CNV.data object, then subset
# ============================================================
cnv_all_f <- CNV.load(tot_int_f)
names(cnv_all_f) <- sample_ids

cnv_ref <- cnv_all_f[normal_idx]

data(detail_regions)
anno <- CNV.create_anno(
  array_type     = "EPIC",
  detail_regions = detail_regions,
  features       = rownames(tot_int_f)
)

# ============================================================
# 3. Process each query sample
# ============================================================
all_segments <- list()
all_details  <- list()

for (i in seq_along(query_idx)) {
  idx <- query_idx[i]
  sample_name  <- sample_ids[idx]
  sample_group <- as.character(group[idx])
  
  if (i %% 10 == 1) {
    cat("  Processing sample", i, "of", length(query_idx), ":", sample_name, "\n")
  }
  
  tryCatch({
    cnv_query <- cnv_all_f[idx]
    
    x <- CNV.fit(cnv_query, cnv_ref, anno)
    x <- CNV.bin(x)
    x <- CNV.detail(x)
    x <- CNV.segment(x)
    
    segs <- CNV.write(x, what = "segments")
    if (is.data.frame(segs) && NROW(segs) > 0) {
      segs$Sample <- sample_name
      segs$Group  <- sample_group
      all_segments[[sample_name]] <- segs
    }
    
    det <- CNV.write(x, what = "detail")
    if (is.data.frame(det) && NROW(det) > 0) {
      det$Sample <- sample_name
      det$Group  <- sample_group
      all_details[[sample_name]] <- det
    }
    
  }, error = function(e) {
    cat("  ERROR for ", sample_name, ": ", e$message, "\n", sep = "")
  })
}

# ============================================================
# 5. Export compiled results
# ============================================================
segments_df <- if (length(all_segments)) bind_rows(all_segments) else tibble()
details_df  <- if (length(all_details))  bind_rows(all_details)  else tibble()

out_segments <- file.path(results_dir, "CNV_segments_all.csv")
out_details  <- file.path(results_dir, "CNV_detail_regions_all.csv")

if (nrow(segments_df) > 0) {
  readr::write_csv(segments_df, out_segments)
  cat("Wrote:", out_segments, "\n")
} else {
  cat("No segments to write.\n")
}

if (nrow(details_df) > 0) {
  readr::write_csv(details_df, out_details)
  cat("Wrote:", out_details, "\n")
} else {
  cat("No detail regions to write.\n")
}

# Optional: per-sample burden
if (nrow(segments_df) > 0) {
  cnv_burden <- segments_df %>%
    mutate(
      seg_length = end - start,
      is_gain = seg.mean > 0.1,
      is_loss = seg.mean < -0.1
    ) %>%
    group_by(Sample, Group) %>%
    summarise(
      n_segments = n(),
      total_gain_Mb = sum(seg_length[is_gain], na.rm = TRUE) / 1e6,
      total_loss_Mb = sum(seg_length[is_loss], na.rm = TRUE) / 1e6,
      total_CNV_Mb  = total_gain_Mb + total_loss_Mb,
      mean_abs_seg_mean = mean(abs(seg.mean), na.rm = TRUE),
      .groups = "drop"
    )
  
  out_burden <- file.path(results_dir, "CNV_burden_by_sample.csv")
  readr::write_csv(cnv_burden, out_burden)
  cat("Wrote:", out_burden, "\n")
}

# ============================================================
# 4. Frequency plot: recurrent gains/losses
# ============================================================
cat("\n--- Generating CNV frequency plot ---\n")

# Bin genome into windows and compute gain/loss frequency per group
window_size <- 1e6  # 1Mb windows

# Create genomic windows
chr_lengths <- c(chr1=249250621, chr2=243199373, chr3=198022430,
                 chr4=191154276, chr5=180915260, chr6=171115067,
                 chr7=159138663, chr8=146364022, chr9=141213431,
                 chr10=135534747, chr11=135006516, chr12=133851895,
                 chr13=115169878, chr14=107349540, chr15=102531392,
                 chr16=90354753, chr17=81195210, chr18=78077248,
                 chr19=59128983, chr20=63025520, chr21=48129895,
                 chr22=51304566)

windows <- do.call(rbind, lapply(names(chr_lengths), function(chr) {
  starts <- seq(1, chr_lengths[chr], by = window_size)
  data.frame(chr = chr, start = starts,
             end = pmin(starts + window_size - 1, chr_lengths[chr]))
}))

# Map segments to windows and compute frequency
compute_freq <- function(segs, group_name, total_n) {
  segs_group <- segs %>% filter(Group == group_name)
  if (nrow(segs_group) == 0) return(NULL)

  segs_gr <- GRanges(seqnames = segs_group$chrom,
                      ranges = IRanges(start = segs_group$loc.start,
                                       end = segs_group$loc.end))
  segs_gr$seg.mean <- segs_group$seg.mean
  segs_gr$Sample <- segs_group$Sample

  win_gr <- GRanges(seqnames = windows$chr,
                     ranges = IRanges(start = windows$start, end = windows$end))

  ol <- findOverlaps(win_gr, segs_gr)
  if (length(ol) == 0) return(NULL)

  # For each window, count unique samples with gain or loss
  gain_samples <- tapply(
    mcols(segs_gr)$Sample[subjectHits(ol)][mcols(segs_gr)$seg.mean[subjectHits(ol)] > 0.1],
    queryHits(ol)[mcols(segs_gr)$seg.mean[subjectHits(ol)] > 0.1],
    function(x) length(unique(x))
  )

  loss_samples <- tapply(
    mcols(segs_gr)$Sample[subjectHits(ol)][mcols(segs_gr)$seg.mean[subjectHits(ol)] < -0.1],
    queryHits(ol)[mcols(segs_gr)$seg.mean[subjectHits(ol)] < -0.1],
    function(x) length(unique(x))
  )

  windows$gain_freq <- 0
  windows$loss_freq <- 0
  if (length(gain_samples) > 0) {
    windows$gain_freq[as.integer(names(gain_samples))] <- as.numeric(gain_samples) / total_n
  }
  if (length(loss_samples) > 0) {
    windows$loss_freq[as.integer(names(loss_samples))] <- -as.numeric(loss_samples) / total_n
  }
  windows$Group <- group_name
  windows
}

tumor_n <- sum(group == "1_Tumor")
opl_n <- sum(group == "2_OPL")

freq_tumor <- compute_freq(segments_df, "1_Tumor", tumor_n)
freq_opl <- compute_freq(segments_df, "2_OPL", opl_n)

freq_all <- bind_rows(freq_tumor, freq_opl)

if (!is.null(freq_all) && nrow(freq_all) > 0) {
  freq_all <- freq_all %>%
    mutate(
      chr_num = as.numeric(gsub("chr", "", chr)),
      mid = (start + end) / 2,
      Group_Label = ifelse(Group == "1_Tumor", "Tumor", "OPL")
    )

  p_freq <- ggplot(freq_all, aes(x = mid / 1e6)) +
    geom_bar(aes(y = gain_freq), stat = "identity", fill = "#E64B35", alpha = 0.6) +
    geom_bar(aes(y = loss_freq), stat = "identity", fill = "#3C5488", alpha = 0.6) +
    facet_grid(Group_Label ~ chr, scales = "free_x", space = "free_x") +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    labs(
      title = "Copy Number Alterations Inferred from EPIC Array",
      subtitle = "Conumee | Reference: 22 Normal oral mucosa samples",
      x = "Genomic Position (Mb)",
      y = "Frequency (gain ↑ / loss ↓)"
    ) +
    theme_classic(base_size = 8) +
    theme(
      strip.text.x = element_text(size = 6),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.spacing.x = unit(0.1, "lines")
    )

  ggsave(file.path(fig_dir, "CNV_frequency_plot.pdf"), p_freq,
         width = 16, height = 5, dpi = 300)
  ggsave(file.path(fig_dir, "CNV_frequency_plot.png"), p_freq,
         width = 16, height = 5, dpi = 300)
}

# ============================================================
# 5. Overlap with PMD blocks
# ============================================================
cat("\n--- CNV-PMD overlap analysis ---\n")

hypo_blocks <- read.csv(file.path(results_dir, "hypo_PMD_blocks.csv"),
                         stringsAsFactors = FALSE)
pmd_gr <- GRanges(seqnames = hypo_blocks$seqnames,
                   ranges = IRanges(start = hypo_blocks$start, end = hypo_blocks$end))

# Check if CNV losses overlap PMDs
if (nrow(segments_df) > 0) {
  loss_segs <- segments_df %>% filter(seg.mean < -0.1)
  if (nrow(loss_segs) > 0) {
    loss_gr <- GRanges(seqnames = loss_segs$chrom,
                        ranges = IRanges(start = loss_segs$loc.start,
                                         end = loss_segs$loc.end))

    ol <- findOverlaps(loss_gr, pmd_gr)
    cat("CNV loss segments overlapping PMDs:", length(unique(queryHits(ol))),
        "of", length(loss_gr), "\n")

    # Fisher test
    total_genome <- sum(chr_lengths)
    pmd_coverage <- sum(width(pmd_gr))
    loss_coverage <- sum(width(loss_gr))
    overlap_coverage <- sum(width(pintersect(loss_gr[queryHits(ol)], pmd_gr[subjectHits(ol)])))

    cat("PMD genome coverage:", round(pmd_coverage / total_genome * 100, 1), "%\n")
    cat("Loss segment coverage:", round(loss_coverage / total_genome * 100, 1), "%\n")
    cat("Overlap:", round(overlap_coverage / total_genome * 100, 2), "%\n")
  }
}

# ============================================================
# 6. Correlate CNV burden with PMD load
# ============================================================
cat("\n--- CNV burden vs PMD load ---\n")

pmd_load <- read.csv(file.path(results_dir, "PMD_load_per_sample.csv"),
                      stringsAsFactors = FALSE)

cnv_pmd <- cnv_burden %>%
  left_join(pmd_load %>% select(Sample, PMD_load_Mb), by = "Sample") %>%
  filter(!is.na(PMD_load_Mb))

if (nrow(cnv_pmd) > 5) {
  cor_test <- cor.test(cnv_pmd$total_CNV_Mb, cnv_pmd$PMD_load_Mb, method = "spearman")
  cat("Total CNV burden vs PMD load: rho =", round(cor_test$estimate, 3),
      ", p =", format.pval(cor_test$p.value), "\n")

  cnv_pmd$Group_Label <- factor(cnv_pmd$Group,
                                 levels = c("2_OPL", "1_Tumor"),
                                 labels = c("OPL", "Tumor"))

  p_cnv_pmd <- ggplot(cnv_pmd, aes(x = PMD_load_Mb, y = total_CNV_Mb, color = Group_Label)) +
    geom_point(size = 2, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed",
                linewidth = 0.5) +
    scale_color_manual(values = c("OPL" = "#00A087", "Tumor" = "#E64B35")) +
    annotate("text", x = Inf, y = Inf,
             label = paste0("rho = ", round(cor_test$estimate, 3),
                            "\np = ", format.pval(cor_test$p.value, digits = 2)),
             hjust = 1.1, vjust = 1.1, size = 4) +
    labs(
      title = "CNV Burden vs PMD Erosion Load",
      x = "PMD Load (Mb)",
      y = "Total CNV Burden (Mb)",
      color = "Group"
    ) +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom")

  ggsave(file.path(fig_dir, "CNV_vs_PMD_scatter.pdf"), p_cnv_pmd,
         width = 7, height = 6, dpi = 300)
  ggsave(file.path(fig_dir, "CNV_vs_PMD_scatter.png"), p_cnv_pmd,
         width = 7, height = 6, dpi = 300)
}

# ============================================================
# 7. Check specific loci
# ============================================================
cat("\n--- Checking specific loci ---\n")

loci_of_interest <- data.frame(
  Gene = c("CDKN2A", "CCND1", "ZNF471", "TFAP2A"),
  Chr = c("chr9", "chr11", "chr19", "chr6"),
  Start = c(21967751, 69455873, 57016865, 10395687),
  End = c(21995301, 69469242, 57041611, 10541768),
  stringsAsFactors = FALSE
)

for (i in 1:nrow(loci_of_interest)) {
  gene <- loci_of_interest$Gene[i]
  chr <- loci_of_interest$Chr[i]
  s <- loci_of_interest$Start[i]
  e <- loci_of_interest$End[i]

  gene_segs <- segments_df %>%
    filter(chrom == chr, loc.start <= e, loc.end >= s)

  if (nrow(gene_segs) > 0) {
    gains <- sum(gene_segs$seg.mean > 0.1)
    losses <- sum(gene_segs$seg.mean < -0.1)
    cat(gene, "(", chr, "):",
        gains, "gains,", losses, "losses in",
        length(unique(gene_segs$Sample)), "samples\n")
  } else {
    cat(gene, "(", chr, "): No CNV segments overlap\n")
  }
}

# ============================================================
# 8. Save results
# ============================================================
cat("\n--- Saving results ---\n")

write.csv(segments_df, file.path(results_dir, "CNV_segments_all.csv"),
          row.names = FALSE)
write.csv(cnv_burden, file.path(results_dir, "CNV_frequency_by_group.csv"),
          row.names = FALSE)

# ============================================================
# 9. Source data XLSX
# ============================================================
wb <- createWorkbook()

addWorksheet(wb, "README")
writeData(wb, "README", data.frame(
  Field = c("Analysis", "Method", "Reference", "Date"),
  Value = c("Copy Number Inference from EPIC Methylation Array",
            "conumee (combined intensity-based CNV calling)",
            "22 Normal oral mucosa samples as reference",
            format(Sys.time(), "%Y-%m-%d"))
))

addWorksheet(wb, "CNV_Burden_PerSample")
writeData(wb, "CNV_Burden_PerSample", cnv_burden)

addWorksheet(wb, "EDFig_CNV_Frequency")
if (exists("freq_all") && !is.null(freq_all)) {
  freq_summary <- freq_all %>%
    group_by(chr, Group_Label) %>%
    summarise(
      max_gain_freq = max(gain_freq, na.rm = TRUE),
      max_loss_freq = min(loss_freq, na.rm = TRUE),
      .groups = "drop"
    )
  writeData(wb, "EDFig_CNV_Frequency", freq_summary)
}

addWorksheet(wb, "Loci_of_Interest")
writeData(wb, "Loci_of_Interest", loci_of_interest)

saveWorkbook(wb, file.path(source_dir, "SourceData_ED_CNV.xlsx"),
             overwrite = TRUE)

cat("\n=== Script 27 complete ===\n")
