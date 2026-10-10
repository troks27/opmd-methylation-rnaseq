## Stage 01 / step 1. GSE204943 (Infinium MethylationEPIC v1; Tata Memorial, Mumbai):
## raw IDATs -> detection-P filtering -> noob-normalised beta values.
##
## Inputs : IDATs (config: gse204943_idat; downloaded by scripts/00_download) and the GEO
##          sample sheet GSE204943_Sample_sheet_118_samples_MBM_MI_modified.csv(.gz),
##          expected in the parent folder of the IDAT directory.
## Outputs: <outputs>/gse204943/crucial_objects/beta.rds     CpG x sample betas (colnames = Basename)
##                                             pd.rds       sample sheet, same order as beta columns
##                                             ann_auto.rds hg19 annotation of the autosomal CpGs in beta
##          <outputs>/gse204943/intermediate/RGset_qc.rds, Mset_noob.rds  (inputs to steps 3 and 4)
##
## Filtering (as described in Methods): a probe "fails" in a sample when its detection
## P > 0.01. Samples with >= 1% failed probes are removed; probes failing in >= 1% of the
## remaining samples are removed after noob normalisation.
##
## Note on the original pipeline: the first implementation applied the probe-level keep
## vector to the RGChannelSet, whose rows are bead addresses rather than probes, so the
## logical index was silently recycled. Here the filter is applied to the MethylSet, whose
## rows are probes (docs/CHANGELOG.md).
##
## Run from the repository root:  Rscript scripts/01_gse204943/01_minfi_preprocess.R  (~20 min, ~10 GB RAM)

source("R/lib/config.R")
suppressPackageStartupMessages({
  library(minfi); library(IlluminaHumanMethylationEPICmanifest)
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19); library(tidyverse)
})
DETP_CUT <- 0.01; MAX_FAIL_SAMPLE <- 0.01; MAX_FAIL_PROBE <- 0.01

G <- file.path(OUT_DIR, "gse204943")
for (d in c("crucial_objects", "results", "intermediate")) dir.create(file.path(G, d), recursive = TRUE, showWarnings = FALSE)

sheet <- list.files(dirname(GSE204943_IDAT), pattern = "Sample_sheet.*\\.csv(\\.gz)?$", full.names = TRUE)
require_files(sheet[1])
targets <- read.csv(sheet[1], stringsAsFactors = FALSE)
targets <- targets[nzchar(targets$Sample_Name), ]
cat(sprintf("Sample sheet: %d samples (%s)\n", nrow(targets), paste(names(table(targets$Sample_Group)), table(targets$Sample_Group), collapse = ", ")))

RGset <- read.metharray.exp(base = GSE204943_IDAT, force = TRUE)
barcode <- sub("^GSM[0-9]+_", "", sampleNames(RGset))          # GEO prefixes IDATs with the GSM accession
idx <- match(barcode, targets$Basename)
if (anyNA(idx)) stop("IDATs without a sample-sheet row: ", paste(head(barcode[is.na(idx)]), collapse = ", "))
colData(RGset) <- S4Vectors::DataFrame(targets[idx, ])
sampleNames(RGset) <- barcode

## sample filter
detP <- detectionP(RGset)
fail_sample <- colMeans(detP > DETP_CUT)
keep_s <- fail_sample < MAX_FAIL_SAMPLE
cat(sprintf("Samples removed for >= %.0f%% failed probes: %d\n", 100 * MAX_FAIL_SAMPLE, sum(!keep_s)))
RGset <- RGset[, keep_s]; detP <- detP[, keep_s]

## noob, then probe filter on the MethylSet (rows = probes)
Mset_noob <- preprocessNoob(RGset)
detP <- detP[featureNames(Mset_noob), ]
keep_p <- rowMeans(detP > DETP_CUT) < MAX_FAIL_PROBE
Mset_noob <- Mset_noob[keep_p, ]
cat(sprintf("Probes removed for failing in >= %.0f%% of samples: %d; retained: %d\n",
            100 * MAX_FAIL_PROBE, sum(!keep_p), sum(keep_p)))

beta <- getBeta(Mset_noob)
pd <- as.data.frame(colData(Mset_noob))
stopifnot(identical(colnames(beta), pd$Basename))
ann <- as.data.frame(getAnnotation(Mset_noob))
ann_auto <- ann[ann$chr %in% paste0("chr", 1:22), ]

## ZNF471 promoter mean: used downstream as a column-alignment check
znf <- rownames(ann)[ann$UCSC_RefGene_Name == "ZNF471" & ann$UCSC_RefGene_Group %in% c("TSS200", "TSS1500")]
pd$ZNF471_prom_beta <- colMeans(beta[intersect(znf, rownames(beta)), , drop = FALSE], na.rm = TRUE)

saveRDS(beta, file.path(G, "crucial_objects", "beta.rds"))
saveRDS(pd, file.path(G, "crucial_objects", "pd.rds"))
saveRDS(ann_auto, file.path(G, "crucial_objects", "ann_auto.rds"))
saveRDS(RGset, file.path(G, "intermediate", "RGset_qc.rds"))
saveRDS(Mset_noob, file.path(G, "intermediate", "Mset_noob.rds"))
cat(sprintf("Saved betas: %d CpGs x %d samples\n", nrow(beta), ncol(beta)))
