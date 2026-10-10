# ==============================================================================
# Script 01: Data Processing, QC, and Differential Methylation Analysis
# Project: OPMD Methylation Analysis (GSE204943)
# Author: Michael E. Troka
#
# Purpose:
#   End-to-end EPIC array pipeline: download IDATs from GEO, QC filtering,
#   Noob normalization, ComBat batch correction, DMP analysis (3 contrasts),
#   gene-centric summaries (TFAP2A, ZNF471), DMR calling (DMRcate), elastic-net
#   classifier for OPL tumor-likeness scoring
#
# Inputs:
#   - GSE204943 IDAT files (auto-downloaded from GEO)
#   - GSE204943_Sample_sheet_118_samples_MBM_MI_modified.csv
#
# Outputs:
#   - analysis_objects.rds
#   - DMPs_*_annot.csv
#   - DMRs_Tumor_vs_Normal.csv/.rds
#   - drivers_sorted.csv
#   - progressive_*_CpGs.csv
#   - QC/PCA/gene methylation figures
#
# Dependencies:
#   minfi, limma, sva, glmnet, DMRcate,
#   IlluminaHumanMethylationEPICanno.ilm10b4.hg19, ggplot2, dplyr
#
# Upstream: None (entry point)
# ==============================================================================

############################################################
# GSE204943 EPIC methylation analysis
# Author: Michael E. Troka
# Date: 2025-11-22
# Purpose:
#   1) Download & unpack IDATs from GEO
#   2) QC via detection P filtering
#   3) Noob normalization + ComBat batch correction
#   4) DMP analysis (Tumor vs Normal, OPL vs Normal, Tumor vs OPL)
#   5) Gene-centric summaries (TFAP2A, ZNF471)
#   6) DMR calling (DMRcate)
#   7) Elastic-net classifier Tumor vs Normal; score OPL tumor-likeness
############################################################

########################
## 0. Parameters
########################
params <- list(
  gse_id         = "GSE204943",
  base_dir       = "GSE204943_data",
  idat_subdir    = "idat",
  targets_csv    = "GSE204943_Sample_sheet_118_samples_MBM_MI_modified.csv",

  # QC cutoffs
  detP_probe_cut   = 0.01,
  detP_samp_frac   = 0.01,  # remove samples with >1% failed probes
  detP_probe_frac  = 0.01,  # remove probes with >1% failed samples

  # limma thresholds
  fdr_cut        = 0.05,
  logfc_cut      = 1,

  # tumor-likeness thresholds
  tumor_like_prob = 0.7,
  normal_like_prob = 0.3,

  seed = 123
)

########################
## 1. Packages
########################
suppressPackageStartupMessages({
  library(minfi)
  library(GEOquery)
  library(R.utils)
  library(S4Vectors)
  library(sva)
  library(limma)
  library(DMRcate)
  library(glmnet)
  library(ggplot2)
})

########################
## 2. Helpers
########################

ensure_dir <- function(path) if (!dir.exists(path)) dir.create(path, recursive = TRUE)

# Download and unpack GEO IDAT tarball
download_and_unpack_idats <- function(gse_id, base_dir, idat_subdir) {
  getGEOSuppFiles(gse_id, baseDir = base_dir, makeDirectory = TRUE)
  gse_dir <- file.path(base_dir, gse_id)
  idat_dir <- file.path(gse_dir, idat_subdir)
  ensure_dir(idat_dir)

  tar_file <- file.path(gse_dir, paste0(gse_id, "_RAW.tar"))
  untar(tar_file, exdir = idat_dir)

  # gunzip all gz
  gz_files <- list.files(idat_dir, pattern = "\\.gz$", full.names = TRUE, recursive = TRUE)
  for (f in gz_files) gunzip(f, overwrite = TRUE, remove = TRUE)

  return(idat_dir)
}

# Attach targets to RGset safely
attach_targets <- function(RGset, targets) {
  sn <- sampleNames(RGset)

  # Use full basename; only strip color suffix if present
  barcode <- sub("_(Grn|Red)$", "", sn)

  idx <- match(barcode, targets$Basename)

  if (any(is.na(idx))) {
    unmatched <- barcode[is.na(idx)]
    message("Example sampleNames(RGset):"); print(head(sn, 6))
    message("Example extracted barcodes:"); print(head(barcode, 6))
    message("Example targets$Basename:");  print(head(targets$Basename, 6))
    message("Unmatched examples:");        print(head(unmatched, 20))
    stop("Some RGset samples not found in targets$Basename.")
  }

  targets_matched <- targets[idx, , drop = FALSE]
  colData(RGset) <- S4Vectors::DataFrame(targets_matched)
  RGset
}

# QC filter based on detection P
qc_filter_detectionP <- function(RGset, detP, samp_frac_cut, probe_frac_cut, detP_cut) {
  # samples
  failed_frac_samp <- colMeans(detP > detP_cut)
  keep_samples <- failed_frac_samp < samp_frac_cut
  RGset <- RGset[, keep_samples]
  detP  <- detP[, keep_samples]

  # probes
  failed_frac_probe <- rowMeans(detP > detP_cut)
  keep_probes <- failed_frac_probe < probe_frac_cut
  RGset <- RGset[keep_probes, ]

  list(RGset = RGset, detP = detP,
       failed_frac_samp = failed_frac_samp,
       failed_frac_probe = failed_frac_probe)
}

# Run limma DMP for a coefficient or contrast
run_limma <- function(Mmat, design, coef = NULL, contrast = NULL, adjust = "BH") {
  fit <- lmFit(Mmat, design)
  if (!is.null(contrast)) {
    fit <- contrasts.fit(fit, contrast)
  }
  fit <- eBayes(fit)
  topTable(fit, coef = coef, number = Inf, adjust.method = adjust)
}

# Add EPIC annotation columns
annotate_dmps <- function(DMP_table, ann) {
  common <- intersect(rownames(DMP_table), rownames(ann))
  cbind(
    DMP_table[common, ],
    ann[common, c("chr", "pos", "UCSC_RefGene_Name",
                  "UCSC_RefGene_Group", "Relation_to_Island")]
  )
}

# Pull CpGs by gene + region
get_gene_cpgs <- function(ann, gene, regions) {
  rownames(ann)[ann$UCSC_RefGene_Name %in% gene &
                  ann$UCSC_RefGene_Group %in% regions]
}

########################
## 3. I/O setup
########################
gse_dir   <- file.path(params$base_dir, params$gse_id)
idat_dir  <- file.path(gse_dir, params$idat_subdir)
results_dir <- file.path(gse_dir, "results")
fig_dir     <- file.path(gse_dir, "figures")
ensure_dir(results_dir); ensure_dir(fig_dir)

########################
## 4. Download / load IDATs + targets
########################
if (!dir.exists(idat_dir) || length(list.files(idat_dir, "\\.idat$")) == 0) {
  idat_dir <- download_and_unpack_idats(params$gse_id, params$base_dir, params$idat_subdir)
}

targets <- read.csv(file.path(gse_dir, params$targets_csv), stringsAsFactors = FALSE)

RGset <- read.metharray.exp(base = idat_dir, force = TRUE)

RGset <- attach_targets(RGset, targets)

########################
## 5. QC: detection p-values filter
########################
detP <- detectionP(RGset)

qc <- qc_filter_detectionP(
  RGset, detP,
  samp_frac_cut = params$detP_samp_frac,
  probe_frac_cut = params$detP_probe_frac,
  detP_cut = params$detP_probe_cut
)
RGset.qc <- qc$RGset
detP.qc  <- qc$detP

# Save QC plots
pdf(file.path(fig_dir, "QC_failed_probe_fraction_per_sample.pdf"))
barplot(qc$failed_frac_samp, las = 2,
        main = "Fraction failed probes per sample")
abline(h = params$detP_samp_frac, col = "red")
dev.off()

########################
## 6. Normalization + batch correction
########################
Mset.noob <- preprocessNoob(RGset.qc)
beta <- getBeta(Mset.noob)
Mval <- getM(Mset.noob)

pd <- as.data.frame(colData(Mset.noob))
group <- factor(pd$Sample_Group)
batch <- factor(pd$Sample_Plate)

# Set reference group once
group <- relevel(group, ref = "3_Normal")
design <- model.matrix(~ group)

Mval.bc <- ComBat(dat = as.matrix(Mval), batch = batch, mod = design)

########################
## 7. DMPs (limma)
########################
ann <- getAnnotation(Mset.noob)

# 7a) Tumor vs Normal (coef 2)
DMPs_TvN <- run_limma(Mval.bc, design, coef = 2)
sig_TvN  <- subset(DMPs_TvN, adj.P.Val < params$fdr_cut & abs(logFC) > params$logfc_cut)

# 7b) OPL vs Normal (coef 3)
DMPs_OvN <- run_limma(Mval.bc, design, coef = 3)
sig_OvN  <- subset(DMPs_OvN, adj.P.Val < params$fdr_cut & abs(logFC) > params$logfc_cut)

# 7c) Tumor vs OPL (contrast)
contrast.matrix <- makeContrasts(Tumor_vs_OPL = group1_Tumor - group2_OPL,
                                 levels = design)
DMPs_TvO <- run_limma(Mval.bc, design, coef = 1, contrast = contrast.matrix)
sig_TvO  <- subset(DMPs_TvO, adj.P.Val < params$fdr_cut & abs(logFC) > params$logfc_cut)

# Annotate + save
write.csv(annotate_dmps(DMPs_TvN, ann), file.path(results_dir, "DMPs_Tumor_vs_Normal_annot.csv"))
write.csv(annotate_dmps(DMPs_OvN, ann), file.path(results_dir, "DMPs_OPL_vs_Normal_annot.csv"))
write.csv(annotate_dmps(DMPs_TvO, ann), file.path(results_dir, "DMPs_Tumor_vs_OPL_annot.csv"))

########################
## 8. TFAP2A / ZNF471 summaries
########################
tfap2a_prom <- get_gene_cpgs(ann, "TFAP2A", c("TSS200","TSS1500"))
znf471_prom <- get_gene_cpgs(ann, "ZNF471", c("TSS200","TSS1500"))
tfap2a_body <- get_gene_cpgs(ann, "TFAP2A", c("Body","Exon1","3'UTR"))

pd$TFAP2A_prom_beta <- colMeans(beta[tfap2a_prom, , drop=FALSE], na.rm=TRUE)
pd$ZNF471_prom_beta <- colMeans(beta[znf471_prom, , drop=FALSE], na.rm=TRUE)
pd$TFAP2A_body_beta <- colMeans(beta[tfap2a_body, , drop=FALSE], na.rm=TRUE)

# Significant TFAP2A promoter CpGs only (Tumor vs Normal)
sig_tfap2a_cpgs <- intersect(tfap2a_prom, rownames(sig_TvN))
pd$TFAP2A_sig_prom_beta <- colMeans(beta[sig_tfap2a_cpgs, , drop=FALSE], na.rm=TRUE)

# Plots
p1 <- ggplot(pd, aes(group, TFAP2A_sig_prom_beta, fill=group)) +
  geom_boxplot(alpha=.7) + geom_jitter(width=.15, size=2, alpha=.7) +
  labs(title="TFAP2A Promoter (Significant CpGs Only)", x="Group", y="Mean β") +
  theme_bw(base_size=14) + theme(legend.position="none")
ggsave(file.path(fig_dir, "TFAP2A_sig_promoter_beta.pdf"), p1, width=6, height=4)

p2 <- ggplot(pd, aes(group, TFAP2A_body_beta, fill=group)) +
  geom_boxplot(alpha=.7) + geom_jitter(width=.15, size=2, alpha=.7) +
  labs(title="TFAP2A Gene Body Methylation", x="Group", y="Mean β") +
  theme_bw(base_size=14) + theme(legend.position="none")
ggsave(file.path(fig_dir, "TFAP2A_body_beta.pdf"), p2, width=6, height=4)

p3 <- ggplot(pd, aes(group, ZNF471_prom_beta, fill=group)) +
  geom_boxplot(alpha=.7) + geom_jitter(width=.15, size=2, alpha=.7) +
  labs(title="ZNF471 Promoter Methylation", x="Group", y="Mean β") +
  theme_bw(base_size=14) + theme(legend.position="none")
ggsave(file.path(fig_dir, "ZNF471_promoter_beta.pdf"), p3, width=6, height=4)

p4 <- ggplot(pd, aes(ZNF471_prom_beta, TFAP2A_body_beta, color=group)) +
  geom_point(size=3, alpha=.8) +
  labs(title="Epigenetic Switch: ZNF471 Silencing vs TFAP2A Activation",
       x="ZNF471 promoter β", y="TFAP2A body β") +
  theme_bw(base_size=14)
ggsave(file.path(fig_dir, "Epigenetic_switch_scatter.pdf"), p4, width=6, height=4)

########################
## 9. Progressive CpGs across Normal→OPL→Tumor
########################
mean_by_group <- sapply(levels(group), function(g)
  rowMeans(beta[, group==g, drop=FALSE], na.rm=TRUE)
)
colnames(mean_by_group) <- levels(group)

Normal_col <- "3_Normal"; OPL_col <- "2_OPL"; Tumor_col <- "1_Tumor"

prog_inc_cpgs <- rownames(mean_by_group)[
  mean_by_group[,Normal_col] < mean_by_group[,OPL_col] &
    mean_by_group[,OPL_col]   < mean_by_group[,Tumor_col]
]
prog_dec_cpgs <- rownames(mean_by_group)[
  mean_by_group[,Normal_col] > mean_by_group[,OPL_col] &
    mean_by_group[,OPL_col]   > mean_by_group[,Tumor_col]
]

write.csv(prog_inc_cpgs, file.path(results_dir, "progressive_increasing_CpGs.csv"), row.names=FALSE)
write.csv(prog_dec_cpgs, file.path(results_dir, "progressive_decreasing_CpGs.csv"), row.names=FALSE)

########################
## 10. DMRcate (Tumor vs Normal only)
########################
myAnnot_TvN <- cpg.annotate(
  datatype="array", object=Mval.bc, what="M",
  analysis.type="differential", design=design, coef=2,
  fdr=params$fdr_cut, arraytype="EPICv1"
)
dmr_TvN <- dmrcate(myAnnot_TvN, lambda=1000, C=2)
dmr_TvN_ranges <- extractRanges(dmr_TvN)
saveRDS(dmr_TvN_ranges, file.path(results_dir, "DMRs_Tumor_vs_Normal.rds"))
write.csv(dmr_TvN_ranges, file.path(results_dir, "DMRs_Tumor_vs_Normal.csv"), row.names=FALSE)

########################
## 11. Elastic-net classifier Tumor vs Normal; score OPL
########################
sig_TvN_strict <- subset(DMPs_TvN, adj.P.Val < 0.01)
sig_cpgs_TvN <- rownames(sig_TvN_strict)

beta_sig <- beta[sig_cpgs_TvN, ]

tumor_idx  <- which(group == "1_Tumor")
normal_idx <- which(group == "3_Normal")
train_idx  <- c(tumor_idx, normal_idx)

X_train <- t(beta_sig[, train_idx, drop=FALSE])
y_train <- ifelse(group[train_idx]=="1_Tumor", 1, 0)

set.seed(params$seed)
cvfit <- cv.glmnet(X_train, y_train, family="binomial", alpha=0.5, nfolds=10)

fit_final <- glmnet(X_train, y_train, family="binomial",
                    alpha=0.5, lambda=cvfit$lambda.min)

opl_idx <- which(group=="2_OPL")
X_opl <- t(beta_sig[, opl_idx, drop=FALSE])
opl_prob_tumor <- as.numeric(predict(fit_final, newx=X_opl, type="response"))

pd$TumorProb <- NA
pd$TumorProb[opl_idx] <- opl_prob_tumor

# Category labels
tumorlikeness_category <- ifelse(
  opl_prob_tumor >= params$tumor_like_prob, "Tumor-like",
  ifelse(opl_prob_tumor <= params$normal_like_prob, "Normal-like", "Borderline")
)
tumorlikeness_category <- factor(tumorlikeness_category,
                                 levels=c("Normal-like","Borderline","Tumor-like"))

# Plot OPL tumor-likeness
p_opl <- ggplot(pd[opl_idx,], aes(x="", y=TumorProb)) +
  geom_jitter(width=0.1, size=2) +
  geom_hline(yintercept=params$tumor_like_prob, linetype="dashed") +
  geom_hline(yintercept=params$normal_like_prob, linetype="dashed") +
  labs(title="Tumor-likeness of OPL samples by methylation",
       y="Predicted P(Tumor pattern)", x="") +
  theme_bw(base_size=14)
ggsave(file.path(fig_dir, "OPL_tumorlikeness_probs.pdf"), p_opl, width=5, height=4)

# PCA on TvN CpGs
X_all <- t(beta_sig)
pca <- prcomp(X_all, scale.=TRUE)
df_pca <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=group)

p_pca <- ggplot(df_pca, aes(PC1, PC2, color=group)) +
  geom_point(size=3) + theme_bw(base_size=14) +
  labs(title="PCA on Tumor-vs-Normal DMP CpGs")
ggsave(file.path(fig_dir, "PCA_TvN_sigCpGs.pdf"), p_pca, width=5, height=4)

########################
## 12. Driver CpGs (nonzero elastic-net coefficients)
########################
coef_vec <- coef(fit_final)
coef_df <- data.frame(CpG=rownames(coef_vec),
                      Coefficient=as.numeric(coef_vec)) |>
  subset(CpG!="(Intercept)" & Coefficient!=0)

drivers <- cbind(coef_df,
                 ann[coef_df$CpG, c("chr","pos","UCSC_RefGene_Name",
                                    "UCSC_RefGene_Group","Relation_to_Island")])
drivers_sorted <- drivers[order(abs(drivers$Coefficient), decreasing=TRUE), ]

write.csv(drivers_sorted, file.path(results_dir, "drivers_sorted.csv"), row.names=FALSE)

########################
## 13. Save analysis objects + session
########################
saveRDS(list(
  RGset.qc=RGset.qc, Mset.noob=Mset.noob,
  beta=beta, Mval.bc=Mval.bc,
  DMPs=list(TvN=DMPs_TvN, OvN=DMPs_OvN, TvO=DMPs_TvO),
  sig=list(TvN=sig_TvN, OvN=sig_OvN, TvO=sig_TvO),
  drivers=drivers_sorted,
  pd=pd
), file.path(results_dir, "analysis_objects.rds"))

writeLines(capture.output(sessionInfo()),
           file.path(results_dir, "sessionInfo.txt"))
