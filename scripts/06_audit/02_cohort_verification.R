## Cohort verification checks (October 2026); see
## the dataset verification note kept with the manuscript files.
##   A. Mumbai leukoplakia (GSE204943) outcome and histology by methylome group.
##      Per-lesion grade and status at last follow-up are published in Bhosale
##      et al. 2017 (Transl Oncol) Table S5-A and Inchanalkar et al. 2023
##      (Clin Epigenetics) Table S2; transcribed below. No follow-up times exist.
##   B. GSE204943 signature re-derived with sex as a covariate (normals 17F/5M,
##      OSCC 18F/56M).
##   C. GSE234379 transfer and KZFP replication in HPV-negative patients.
##   D. PMD index with blocks hypomethylated in OSCC alone (the GSE204943
##      blocks were called on OPL + OSCC vs healthy mucosa).
##   E. GSE227919 proxy AUC in disease-specific subsets.
##   F. KZFP promoter CpG accounting (ZNF382 and ZNF529 share a CpG island).
## Run from the repository root:  Rscript scripts/06_audit/02_cohort_verification.R  (or make audit; reads stage 02-04 outputs; ~5 min)

source("R/lib/config.R")   # paths and helpers (TAB, FIG, OUT_DIR, GSE_DIR, ...)
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(limma); library(GenomicRanges)
})
source(file.path(REPO_ROOT, "R/lib/utils.R"))
auc <- function(pos, neg) mean(outer(pos, neg, ">") + 0.5 * outer(pos, neg, "=="))

pd <- readRDS(file.path(GSE_DIR, "crucial_objects/pd.rds"))
pd$Group <- factor(c(`1_Tumor` = "Tumor", `2_OPL` = "OPL", `3_Normal` = "Normal")[pd$Sample_Group],
                   levels = c("Normal", "OPL", "Tumor"))
st <- readRDS(file.path(EXT_DIR, "story_analyses.rds"))

## ================================================================== A. Mumbai OPL outcome
opl <- st$opl |> mutate(id = pd$Sample_Name[match(Sample, pd$Basename)])
mild_dysplasia <- c("OC2061", "OC2067", "OC2300")
transformed    <- c("OC2061", "OC2067", "OC2300", "OC1884")
recurred       <- c("OC2204", "OC2202", "OC2203")
opl <- opl |> mutate(histology = ifelse(id %in% mild_dysplasia, "Mild dysplasia", "Hyperplasia"),
                     outcome = case_when(id %in% transformed ~ "Transformed to SCC",
                                         id %in% recurred ~ "Lesion recurrence",
                                         TRUE ~ "No evidence of disease"))
stopifnot(nrow(opl) == 22, sum(opl$outcome == "Transformed to SCC") == 4)
tr <- opl$outcome == "Transformed to SCC"; ev <- opl$outcome != "No evidence of disease"
tl <- opl$methylome_state == "tumour-like"
outA <- tibble(test = c("transformed", "transformed_or_recurred"),
               tumour_like = c(sprintf("%d/%d", sum(tr & tl), sum(tl)), sprintf("%d/%d", sum(ev & tl), sum(tl))),
               normal_like = c(sprintf("%d/%d", sum(tr & !tl), sum(!tl)), sprintf("%d/%d", sum(ev & !tl), sum(!tl))),
               fisher_p = c(fisher.test(table(tr, tl))$p.value, fisher.test(table(ev, tl))$p.value),
               auc_score = c(auc(opl$tumor_likeness[tr], opl$tumor_likeness[!tr]),
                             auc(opl$tumor_likeness[ev], opl$tumor_likeness[!ev])))
cat("A. Mumbai OPL outcome by methylome group\n"); print(outA |> mutate(across(where(is.numeric), ~ signif(.x, 3))))
print(table(histology = opl$histology, group = opl$methylome_state))
write_csv(opl |> select(id, tumor_likeness, methylome_state, histology, outcome) |> arrange(desc(tumor_likeness)) |>
            mutate(tumor_likeness = round(tumor_likeness, 3)), TAB("verification_mumbai_opl_outcome.csv"))
write_csv(outA, TAB("verification_mumbai_opl_outcome_tests.csv"))

## ================================================================== B. sex-adjusted signature
ORD <- c("P2-NAT", "P2-Tumor", "P3-Adj", "P3-Lesion", "P4-Adj", "P4-Lesion", "P5-Adj", "P5-Lesion")
p <- readRDS(file.path(OUT_DIR, "processed.rds")); ss <- p$sample_sheet
ss$Short <- paste0("P", sub("Patient ", "", ss$Patient), "-",
                   c(Normal = "Adj", Lesion = "Lesion", Tumor = "Tumor")[as.character(ss$Group)])
ss$Short[ss$Short == "P2-Adj"] <- "P2-NAT"
Bp <- p$betas; colnames(Bp) <- ss$Short[match(colnames(Bp), ss$Sample_Name)]; Bp <- Bp[, ORD]; rm(p); gc()
Bg <- readRDS(file.path(GSE_DIR, "crucial_objects/beta.rds")); colnames(Bg) <- pd$Basename
cp <- read_csv(file.path(GSE_DIR, "results/cell_type_proportions.csv"), show_col_types = FALSE)
pd <- pd |> left_join(cp |> select(Basename = Sample, Epi = Epithelial, IC = Immune), by = "Basename")
common <- intersect(rownames(Bp), rownames(Bg)); Bp <- Bp[common, ]; Bg <- Bg[common, ]
sig0 <- readRDS(file.path(EXT_DIR, "gse_signature.rds"))
pick <- function(s) list(hyper = s |> filter(fdr_TvN < 1e-6, d_TvN > 0.2) |> slice_max(t_TvN, n = 2000) |> pull(cpg),
                         hypo  = s |> filter(fdr_TvN < 1e-6, d_TvN < -0.2) |> slice_min(t_TvN, n = 2000) |> pull(cpg))
S0 <- pick(sig0)
plate <- ifelse(pd$Sample_Plate == "", "none", pd$Sample_Plate)
keep <- pd$Group %in% c("Normal", "Tumor"); g <- droplevels(pd$Group[keep])
d <- model.matrix(~ g + pd$Epi[keep] + pd$IC[keep] + plate[keep] + pd$Sex[keep])
f <- eBayes(lmFit(beta_to_m(Bg[, keep]), d))
tt <- topTable(f, coef = 2, number = Inf, sort.by = "none")
sexc <- topTable(f, coef = ncol(d), number = Inf, sort.by = "none")
sig1 <- sig0 |> select(cpg, muN, muT, d_TvN) |>
  left_join(tibble(cpg = rownames(tt), t_TvN = tt$t, fdr_TvN = tt$adj.P.Val), by = "cpg")
S1 <- pick(sig1)
score <- function(B, S, s) {
  sc <- function(cp) { m <- s[match(cp, s$cpg), ]; apply((B[cp, , drop = FALSE] - m$muN) / (m$muT - m$muN), 2, median, na.rm = TRUE) }
  (sc(S$hyper) + sc(S$hypo)) / 2 }
o <- pd$Group == "OPL"; g0 <- score(Bg, S0, sig0); g1 <- score(Bg, S1, sig1)
p0 <- score(Bp, S0, sig0); p1 <- score(Bp, S1, sig1)
shift <- function(x) c(P2 = x[["P2-Tumor"]] - x[["P2-NAT"]], P3 = x[["P3-Lesion"]] - x[["P3-Adj"]],
                       P4 = x[["P4-Lesion"]] - x[["P4-Adj"]], P5 = x[["P5-Lesion"]] - x[["P5-Adj"]])
outB <- tibble(retained = length(intersect(S0$hyper, S1$hyper)) + length(intersect(S0$hypo, S1$hypo)),
               sig_cpgs_sex_fdr05 = sum(sexc[c(S0$hyper, S0$hypo), "adj.P.Val"] < 0.05),
               opl_r = cor(g0[o], g1[o]), opl_same_group = sum((g0[o] > 0.35) == (g1[o] > 0.35)),
               max_pilot_shift_change = max(abs(shift(p1) - shift(p0))))
cat("\nB. Sex-adjusted signature\n"); print(outB)
write_csv(outB, TAB("verification_sex_adjusted_signature.csv"))

## ================================================================== C. HPV in GSE234379
s3 <- readRDS(file.path(EXT_DIR, "gse234379_kzfp.rds"))$s3
s3$kzfp <- rowMeans(s3[, c("ZNF471", "ZNF382", "ZNF529", "ZNF582", "ZNF671")])
pairs3 <- s3 |> count(patient) |> filter(n == 2) |> pull(patient)
pair_test <- function(d) { w <- tidyr::pivot_wider(d |> select(patient, group, kzfp), names_from = group, values_from = kzfp)
  dd <- w$Tumor - w$Normal; tibble(n_pairs = length(dd), frac_up = mean(dd > 0), exact_p = wilcox.test(dd)$p.value) }
outC <- bind_rows(lapply(c("all", "Negative", "Positive"), function(h) {
  d <- if (h == "all") s3 else filter(s3, hpv == h)
  tibble(hpv = h, n = nrow(d), auc = auc(d$tumor_likeness[d$group == "Tumor"], d$tumor_likeness[d$group == "Normal"])) |>
    bind_cols(pair_test(filter(d, patient %in% pairs3))) }))
cat("\nC. GSE234379 by HPV status\n"); print(outC |> mutate(across(where(is.numeric), ~ signif(.x, 3))))
write_csv(outC, TAB("verification_gse234379_hpv.csv"))

## ================================================================== D. tumour-only PMD blocks
ann <- readRDS(file.path(GSE_DIR, "crucial_objects/ann_auto.rds")); ann <- ann[intersect(rownames(ann), rownames(Bg)), ]
os <- ann$Relation_to_Island == "OpenSea"
cg <- GRanges(ann$chr[os], IRanges(ann$pos[os], width = 1)); names(cg) <- rownames(ann)[os]
pmd <- read_csv(file.path(GSE_DIR, "results/hypo_PMD_blocks.csv"), show_col_types = FALSE)
pg <- GRanges(pmd$seqnames, IRanges(pmd$start, pmd$end))
ov <- findOverlaps(cg, pg); cgs <- names(cg)[queryHits(ov)]; blk <- subjectHits(ov)
bm <- rowsum(Bg[cgs, ], blk) / as.vector(table(blk))
dT <- rowMeans(bm[, pd$Group == "Tumor"]) - rowMeans(bm[, pd$Group == "Normal"])
in_t <- unique(cgs[blk %in% as.integer(rownames(bm))[dT < -0.10]])
out_cg <- names(cg)[!overlapsAny(cg, pg)]
idx0 <- colMeans(Bg[out_cg, ]) - colMeans(Bg[unique(cgs), ]); idx1 <- colMeans(Bg[out_cg, ]) - colMeans(Bg[in_t, ])
oo <- match(opl$Sample, colnames(Bg))
outD <- tibble(blocks_with_cpgs = nrow(bm), blocks_tumour_only = sum(dT < -0.10), r_index = cor(idx0, idx1),
               opl_p_original = wilcox.test(idx0[oo] ~ opl$methylome_state)$p.value,
               opl_p_tumour_only = wilcox.test(idx1[oo] ~ opl$methylome_state)$p.value)
cat("\nD. PMD blocks defined in OSCC alone\n"); print(outD)
write_csv(outD, TAB("verification_pmd_tumour_only_blocks.csv"))

## ================================================================== E. GSE227919 subsets
a2 <- readRDS(file.path(EXT_DIR, "gse26549_progression.rds"))$gse227919_proxy |>
  filter(Class != "Cancer", !is.na(Progression_status)) |> mutate(y = Progression_status != "Stable")
sub_auc <- function(lab, d) tibble(subset = lab, n = nrow(d), events = sum(d$y), auc = auc(d$proxy[d$y], d$proxy[!d$y]))
## MK_EOG_004/005/012 share age, sex, site and outcome with MK_B_0004/0005/0012
outE <- bind_rows(sub_auc("all non-cancer with outcome", a2),
                  sub_auc("hyperkeratosis + dysplasia", filter(a2, Class != "Control")),
                  sub_auc("without possible repeat patients", filter(a2, !Sample %in% c("MK_EOG_004", "MK_EOG_005", "MK_EOG_012"))),
                  sub_auc("SCC events vs stable", filter(a2, Progression_status != "Progressed-Dys")))
cat("\nE. GSE227919 proxy AUC by subset\n"); print(outE |> mutate(auc = round(auc, 2)))
write_csv(outE, TAB("verification_gse227919_subsets.csv"))

## ================================================================== F. KZFP CpG accounting
KZ <- c("ZNF471", "ZNF382", "ZNF529", "ZNF582", "ZNF671")
ann <- readRDS(file.path(GSE_DIR, "crucial_objects/ann_auto.rds"))   # full EPIC v1 annotation (section D subset it)
promk <- grepl("TSS200|TSS1500|1stExon", ann$UCSC_RefGene_Group) & ann$Relation_to_Island == "Island"
kset <- setNames(lapply(KZ, function(g) rownames(ann)[promk & grepl(paste0("(^|;)", g, "(;|$)"), ann$UCSC_RefGene_Name)]), KZ)
ext <- readRDS(file.path(EXT_DIR, "external_methylation.rds"))
ku <- unique(unlist(kset))
outF <- tibble(per_gene = paste(sprintf("%s %d", KZ, lengths(kset)), collapse = "; "),
               sum_per_gene = sum(lengths(kset)), unique_cpgs = length(ku),
               znf382_in_znf529 = length(intersect(kset$ZNF382, kset$ZNF529)),
               unique_in_signature = sum(ku %in% c(ext$tum_hyper, ext$tum_hypo)))
cat("\nF. KZFP promoter CpGs\n"); print(as.data.frame(outF))
write_csv(outF, TAB("verification_kzfp_cpg_accounting.csv"))
