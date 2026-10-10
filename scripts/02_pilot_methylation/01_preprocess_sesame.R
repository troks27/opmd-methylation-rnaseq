## Pilot EPIC v2 chip: preprocessing, filtering and PCA.
##
## Input : pilot IDATs + GenomeStudio sample sheet (config: pilot_idat_dir, pilot_sample_sheet)
## Output: outputs/processed.rds   list(betas, mvals, betas_uncollapsed, sample_sheet, pca, pca_nl, qc,
##                                       prep_passes, session)
##         outputs/extended/tables/qc_summary.csv
##
## Steps
##   1. Read IDAT pairs; sesame preprocessing "QCDPB" (quality masking, infer Infinium-I
##      channel, non-linear dye-bias correction, pOOBAH detection masking, noob).
##   2. Drop probes masked in > 25% of samples, probes off chr1-22, and incomplete rows;
##      collapse EPIC v2 replicate probes to one beta per CpG (betasCollapseToPfx).
##   3. PCA on the 5,000 most variable CpGs (all 8 samples), and on Normal + Lesion only.
##
## Preprocessing passes. The intended method applies QCDPB once (default). The first
## version of this analysis (exploratory report, April-October 2026) called
## openSesame() on already-preprocessed SigDFs, which in sesame 1.28 applies QCDPB a
## second time. Set SESAME_PREP_PASSES=2 to reproduce those betas exactly
## (docs/CHANGELOG.md documents the effect on reported numbers).
##
## Run from the repository root:  Rscript scripts/02_pilot_methylation/01_preprocess_sesame.R

source("R/lib/config.R")
suppressPackageStartupMessages({
  library(sesame); library(sesameData); library(GenomicRanges); library(tidyverse)
})
source(file.path(REPO_ROOT, "R/lib/sample_sheet.R"))
source(file.path(REPO_ROOT, "R/lib/utils.R"))

PREP_PASSES <- as.integer(Sys.getenv("SESAME_PREP_PASSES", "1"))
stopifnot(PREP_PASSES %in% 1:2)
N_TOP <- 5000                     # most variable CpGs for PCA
MAX_MASKED <- 0.25                # drop probes masked in more than this fraction of samples

require_files(PILOT_SHEET)
ss <- build_sample_sheet(PILOT_SHEET, PILOT_IDAT_DIR)
cat(sprintf("Pilot samples: %d (%s)\n", nrow(ss), paste(table(ss$Group), names(table(ss$Group)), collapse = ", ")))

## ---------------------------------------------------------------- 1. preprocessing
sesameDataCacheAll()   # one-time download of sesame annotation; idempotent
raw <- lapply(setNames(ss$Basename, ss$Sample_Name), readIDATpair, platform = "EPICv2")
sdfs <- lapply(raw, prepSesame, prep = "QCDPB")
if (PREP_PASSES == 2) sdfs <- lapply(sdfs, prepSesame, prep = "QCDPB")   # legacy behaviour
betas <- sapply(sdfs, getBetas)
colnames(betas) <- ss$Sample_Name
cat(sprintf("QCDPB passes: %d; probes x samples: %d x %d\n", PREP_PASSES, nrow(betas), ncol(betas)))

safe_sex <- function(b) tryCatch(as.character(inferSex(b, platform = "EPICv2")), error = function(e) NA_character_)
qc <- tibble(Sample = ss$Sample_Name, Group = as.character(ss$Group), Patient = ss$Patient,
             Unmasked_fraction = round(colMeans(!is.na(betas)), 4),
             pOOBAH_pass_rate = round(vapply(raw, function(s) mean(pOOBAH(s, return.pval = TRUE) < 0.05, na.rm = TRUE), 0), 4),
             Mean_intensity = round(vapply(raw, meanIntensity, 0), 1),
             Inferred_sex = vapply(seq_len(ncol(betas)), function(i) safe_sex(betas[, i]), ""))
print(qc)
write_csv(qc, TAB("qc_summary.csv"))

## ---------------------------------------------------------------- 2. filtering
betas_keep <- betas[rowMeans(is.na(betas)) <= MAX_MASKED, , drop = FALSE]
manifest <- sesameData_getManifestGRanges("EPICv2")
auto <- names(manifest)[as.character(seqnames(manifest)) %in% paste0("chr", 1:22)]
betas_keep <- betas_keep[rownames(betas_keep) %in% auto, , drop = FALSE]
betas_u <- betas_keep[complete.cases(betas_keep), , drop = FALSE]   # replicate probes kept (DMRcate)
betas_f <- betasCollapseToPfx(betas_keep)
betas_f <- betas_f[complete.cases(betas_f), , drop = FALSE]
mvals <- beta_to_m(betas_f)
cat(sprintf("CpGs after filtering and collapsing: %d (of %d probes)\n", nrow(betas_f), nrow(betas)))

## ---------------------------------------------------------------- 3. PCA
v <- matrixStats::rowVars(betas_f)
top <- order(v, decreasing = TRUE)[seq_len(min(N_TOP, length(v)))]
pca <- prcomp(t(betas_f[top, ]), scale. = FALSE)
nl <- ss$Group %in% c("Normal", "Lesion")
top_nl <- order(matrixStats::rowVars(betas_f[, nl]), decreasing = TRUE)[seq_len(N_TOP)]
pca_nl <- prcomp(t(betas_f[top_nl, nl]), scale. = FALSE)
pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
cat(sprintf("PC1 %.1f%%, PC2 %.1f%% of variance\n", pct[1], pct[2]))

saveRDS(list(betas = betas_f, mvals = mvals, betas_uncollapsed = betas_u, sample_sheet = ss, pca = pca, pca_nl = pca_nl, qc = qc,
             prep_passes = PREP_PASSES, session = sessionInfo()),
        file.path(OUT_DIR, "processed.rds"))
cat("Saved", file.path(OUT_DIR, "processed.rds"), "\n")
