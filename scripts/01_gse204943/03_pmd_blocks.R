## Stage 01 / step 3. GSE204943 partially methylated domain (PMD) blocks and per-sample erosion.
##
## 1. Blocks. CpGs are collapsed into clusters (minfi::cpgCollapse) and hypomethylated blocks
##    (no minimum length is imposed) are called with minfi::blockFinder on open-sea clusters
##    of all chromosomes, comparing OPL + OSCC
##    ("Case") with healthy mucosa ("Control") on the beta scale (cutoff 0.10, smoothing on).
##    Blocks with a negative smoothed Case - Control difference are kept.
## 2. Per-sample erosion, computed on individual CpGs:
##      delta_in_minus_out = mean beta inside blocks - mean beta outside blocks
##      a block is "eroded" in a sample when its mean beta is < 0.50 and at least 0.05
##      below that sample's mean outside blocks; PMD_load_Mb = total width of eroded blocks
##      PMD_erosion_score = mean of z(PMD_load_Mb) and z(-delta_in_minus_out) across samples
##    The score is used only as a consistency check of the PMD index in later scripts.
##
## Input : <outputs>/gse204943/intermediate/Mset_noob.rds, crucial_objects/pd.rds
## Output: <outputs>/gse204943/results/hypo_PMD_blocks.csv, PMD_erosion_score_per_sample.csv
## Permutations (B) only produce block P values; block coordinates do not depend on them and
## only coordinates are used downstream. The original run used B = 100
## (PMD_PERMUTATIONS=100 Rscript ...); the default of 0 skips them.
## Run from the repository root:  Rscript scripts/01_gse204943/03_pmd_blocks.R  (~30 min)

source("R/lib/config.R")
suppressPackageStartupMessages({ library(minfi); library(GenomicRanges); library(tidyverse) })
CUTOFF <- 0.10; SAMPLE_BETA_CUT <- 0.50; SAMPLE_DELTA_OUT <- 0.05
B <- as.integer(Sys.getenv("PMD_PERMUTATIONS", "0"))

G <- file.path(OUT_DIR, "gse204943")
Mset <- readRDS(file.path(G, "intermediate", "Mset_noob.rds"))
pd <- readRDS(file.path(G, "crucial_objects", "pd.rds"))
grSet <- mapToGenome(ratioConvert(Mset, what = "beta"))
colnames(grSet) <- pd$Basename
rm(Mset); invisible(gc())

## ---- 1. blocks
case <- factor(ifelse(pd$Sample_Group %in% c("1_Tumor", "2_OPL"), "Case", "Control"), levels = c("Control", "Case"))
design <- model.matrix(~ case)
cc <- cpgCollapse(grSet, what = "Beta", returnBlockInfo = FALSE)
collapsed <- if (is.list(cc)) cc$object else cc
bf <- blockFinder(collapsed, design = design, coef = 2, what = "Beta", cluster = rowRanges(collapsed)$blockgroup,
                  cutoff = CUTOFF, pickCutoff = FALSE, smooth = TRUE, B = B, verbose = FALSE)
rm(cc, collapsed); invisible(gc())
hypo <- bf$table[bf$table$value < 0, ]
blocks <- GRanges(hypo$chr, IRanges(hypo$start, hypo$end))
mcols(blocks) <- hypo[, intersect(c("value", "area", "L", "p.value", "fwer", "cluster"), names(hypo))]
write_csv(as.data.frame(blocks), file.path(G, "results", "hypo_PMD_blocks.csv"))
cat(sprintf("Hypomethylated blocks: %d (%.0f Mb)\n", length(blocks), sum(width(blocks)) / 1e6))

## ---- 2. per-sample erosion on individual CpGs
b <- getBeta(grSet)
ov <- findOverlaps(rowRanges(grSet), blocks)
inside <- unique(queryHits(ov))
mean_in <- colMeans(b[inside, ], na.rm = TRUE)
mean_out <- colMeans(b[-inside, ], na.rm = TRUE)
block_mean <- rowsum(b[queryHits(ov), ], subjectHits(ov), na.rm = TRUE) /
  as.vector(table(factor(subjectHits(ov), levels = sort(unique(subjectHits(ov))))))
diff_out <- sweep(block_mean, 2, mean_out)
eroded <- (block_mean < SAMPLE_BETA_CUT) & (diff_out <= -SAMPLE_DELTA_OUT)
w <- width(blocks)[as.integer(rownames(block_mean))]

score <- tibble(Sample = colnames(b),
                PMD_load_bp = colSums(eroded * w, na.rm = TRUE),
                PMD_load_Mb = PMD_load_bp / 1e6,
                n_eroded_blocks = colSums(eroded, na.rm = TRUE),
                mean_beta_in_PMD = mean_in, mean_beta_out_PMD = mean_out,
                delta_in_minus_out = mean_in - mean_out) |>
  mutate(z_load = as.numeric(scale(PMD_load_Mb)), z_delta = as.numeric(scale(-delta_in_minus_out)),
         PMD_erosion_score = (z_load + z_delta) / 2,
         Sample_Group = pd$Sample_Group[match(Sample, pd$Basename)])
write_csv(score, file.path(G, "results", "PMD_erosion_score_per_sample.csv"))
print(score |> group_by(Sample_Group) |> summarise(median_load_Mb = median(PMD_load_Mb), median_score = median(PMD_erosion_score)))
