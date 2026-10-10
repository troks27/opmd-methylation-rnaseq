## Pilot exploratory analysis (not used in the manuscript figures): paired lesion-versus-
## adjacent differential methylation in the three patients with both tissues (P3-P5).
## With n = 3 pairs no single CpG passes FDR control; DMRs are reported as candidate
## regions only (DMRcate at a permissive cpg.annotate FDR of 0.5).
##
## Input : outputs/processed.rds
## Output: outputs/extended/tables/dmp_lesion_vs_normal.csv, dmrs_lesion_vs_normal.csv
## Run from the repository root:  Rscript scripts/02_pilot_methylation/03_exploratory_dmp_dmr.R

source("R/lib/config.R")
suppressPackageStartupMessages({
  library(limma); library(DMRcate); library(IlluminaHumanMethylationEPICv2anno.20a1.hg38); library(tidyverse)
})
source(file.path(REPO_ROOT, "R/lib/utils.R"))
p <- readRDS(file.path(OUT_DIR, "processed.rds")); ss <- p$sample_sheet

keep <- ss$Patient %in% c("Patient 3", "Patient 4", "Patient 5") & ss$Group %in% c("Normal", "Lesion")
ss_p <- droplevels(ss[keep, ])
Group <- factor(as.character(ss_p$Group), levels = c("Normal", "Lesion")); Patient <- factor(ss_p$Patient)
design <- model.matrix(~ 0 + Group + Patient); colnames(design) <- make.names(colnames(design))
contrasts <- makeContrasts(LesionVsNormal = GroupLesion - GroupNormal, levels = design)

## DMPs (collapsed CpGs)
fit2 <- eBayes(contrasts.fit(lmFit(p$mvals[, ss_p$Sample_Name], design), contrasts), robust = TRUE, trend = TRUE)
dmp <- topTable(fit2, coef = "LesionVsNormal", number = Inf, sort.by = "P", adjust.method = "BH")
B_p <- p$betas[, ss_p$Sample_Name]
dmp$delta_beta <- (rowMeans(B_p[, ss_p$Group == "Lesion"]) - rowMeans(B_p[, ss_p$Group == "Normal"]))[rownames(dmp)]
anno <- minfi::getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
anno$cpg <- sub("_[^_]+$", "", anno$Name); anno <- anno[!duplicated(anno$cpg), ]
m <- match(rownames(dmp), anno$cpg)
dmp$chr <- anno$chr[m]; dmp$pos <- anno$pos[m]; dmp$relation_to_island <- anno$Relation_to_Island[m]
dmp$gene <- vapply(strsplit(anno$UCSC_RefGene_Name[m], ";"), function(x) if (!length(x)) NA_character_ else paste(unique(x), collapse = ";"), "")
write.csv(dmp, TAB("dmp_lesion_vs_normal.csv"))
cat(sprintf("DMPs: %d CpGs tested; minimum FDR %.2g\n", nrow(dmp), min(dmp$adj.P.Val)))

## candidate DMRs (uncollapsed EPIC v2 probe IDs)
myanno <- cpg.annotate(object = beta_to_m(p$betas_uncollapsed[, ss_p$Sample_Name]), datatype = "array", what = "M",
                       arraytype = "EPICv2", analysis.type = "differential", design = design, contrasts = TRUE,
                       cont.matrix = contrasts, coef = "LesionVsNormal", fdr = 0.5)
dmrs <- tryCatch(dmrcate(myanno, lambda = 1000, C = 2), error = function(e)
  tryCatch(dmrcate(myanno, lambda = 1000, C = 2, pcutoff = 1e-3), error = function(e) NULL))
if (!is.null(dmrs)) {
  dmr_df <- as.data.frame(extractRanges(dmrs, genome = "hg38"))
  write.csv(dmr_df, TAB("dmrs_lesion_vs_normal.csv"), row.names = FALSE)
  cat(sprintf("Candidate regions: %d\n", nrow(dmr_df)))
}
