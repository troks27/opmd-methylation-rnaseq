## Pilot (EPIC v2, n = 8) placed against an independent oral cohort:
## GSE204943 (EPIC v1; 22 normal, 22 OPL, 74 OSCC; gingivobuccal, India).
##
##   1. cell-type deconvolution (EpiDISH, same reference as the GSE204943 work)
##   2. composition-adjusted tumour signature learned in GSE204943
##   3. "tumour-likeness" score -> reproducible replacement for manuscript Fig. 3B
##   4. pilot samples projected into the GSE204943 methylation PCA
##   5. do the pilot's lesion-vs-adjacent changes replicate the external program?
##   6. PMD (partially methylated domain) erosion and ZNF471 promoter methylation
##   7. promoter methylation of genes named in the manuscript
##
## Run from the repository root:  Rscript scripts/03_methylation_cohorts/01_composition_signature_scores.R  (or make cohorts)

source("R/lib/config.R")   # paths and helpers (TAB, FIG, OUT_DIR, GSE_DIR, ...)
suppressPackageStartupMessages({
  library(tidyverse); library(limma); library(EpiDISH); library(GenomicRanges)
  library(IlluminaHumanMethylationEPICv2anno.20a1.hg38); library(ggrepel); library(patchwork)
})
source(file.path(REPO_ROOT, "R/lib/utils.R"))

ORD <- c("P2-NAT", "P2-Tumor", "P3-Adj", "P3-Lesion", "P4-Adj", "P4-Lesion", "P5-Adj", "P5-Lesion")
HIST <- c(`P2-Tumor` = "T2 OSCC", `P3-Lesion` = "PVL", `P4-Lesion` = "Atypical hyperplasia",
          `P5-Lesion` = "Mod-severe dysplasia (progressed)")
grp_cols <- c(Normal = "#3b6fb6", OPL = "#e08a1e", Tumor = "#c0392b")

## ------------------------------------------------------------------ pilot
p <- readRDS(file.path(OUT_DIR, "processed.rds"))
ss <- p$sample_sheet
ss$Short <- paste0("P", sub("Patient ", "", ss$Patient), "-",
                   c(Normal = "Adj", Lesion = "Lesion", Tumor = "Tumor")[as.character(ss$Group)])
ss$Short[ss$Short == "P2-Adj"] <- "P2-NAT"
Bp <- p$betas
colnames(Bp) <- ss$Short[match(colnames(Bp), ss$Sample_Name)]
Bp <- Bp[, ORD]
ss <- ss[match(ORD, ss$Short), ]
own_adj <- c(`P2-Tumor` = "P2-NAT", `P3-Lesion` = "P3-Adj", `P4-Lesion` = "P4-Adj", `P5-Lesion` = "P5-Adj")

## ------------------------------------------------------------------ GSE204943
pd <- readRDS(file.path(GSE_DIR, "crucial_objects/pd.rds"))
Bg <- readRDS(file.path(GSE_DIR, "crucial_objects/beta.rds"))
colnames(Bg) <- pd$Basename                       # saved without colnames; same order as pd
pd$Group <- factor(c(`1_Tumor` = "Tumor", `2_OPL` = "OPL", `3_Normal` = "Normal")[pd$Sample_Group],
                   levels = c("Normal", "OPL", "Tumor"))

## alignment check: recompute the ZNF471 promoter mean stored in pd
ann_g <- readRDS(file.path(GSE_DIR, "crucial_objects/ann_auto.rds"))
znf <- rownames(ann_g)[grepl("(^|;)ZNF471(;|$)", ann_g$UCSC_RefGene_Name) &
                         grepl("TSS200|TSS1500", ann_g$UCSC_RefGene_Group)]
chk <- cor(colMeans(Bg[intersect(znf, rownames(Bg)), ]), pd$ZNF471_prom_beta)
cat(sprintf("GSE204943 column-alignment check (ZNF471 promoter r): %.4f\n", chk))
stopifnot(chk > 0.999)

cp <- read_csv(file.path(GSE_DIR, "results/cell_type_proportions.csv"), show_col_types = FALSE)
pd <- pd |> left_join(cp |> select(Basename = Sample, Epi = Epithelial, Fib = Fibroblast, IC = Immune),
                      by = "Basename")

common <- intersect(rownames(Bp), rownames(Bg))
cat(sprintf("CpGs: pilot %d, GSE204943 %d, shared %d\n", nrow(Bp), nrow(Bg), length(common)))
Bp <- Bp[common, ]; Bg <- Bg[common, ]

## platform/pipeline concordance: pilot adjacent vs GSE normal, per CpG means
conc_r <- cor(rowMeans(Bp[, c("P3-Adj", "P4-Adj", "P5-Adj", "P2-NAT")]),
              rowMeans(Bg[, pd$Group == "Normal"]))
cat(sprintf("Per-CpG mean beta, pilot adjacent vs GSE normal: r = %.3f\n", conc_r))

## ================================================================== 1. EpiDISH
data(centEpiFibIC.m); data(centBloodSub.m)
est_p <- epidish(beta.m = Bp, ref.m = centEpiFibIC.m, method = "RPC")$estF
hep <- hepidish(beta.m = Bp, ref1.m = centEpiFibIC.m, ref2.m = centBloodSub.m,
                h.CT.idx = 3, method = "RPC")
decon <- as_tibble(hep, rownames = "Sample") |>
  mutate(Epi_1st = est_p[Sample, "Epi"], Fib_1st = est_p[Sample, "Fib"], IC_1st = est_p[Sample, "IC"],
         Group = as.character(ss$Group[match(Sample, ss$Short)]))
pc1 <- p$pca$x[, 1]; names(pc1) <- ss$Short[match(rownames(p$pca$x), ss$Sample_Name)]
decon$PC1_fig3A <- pc1[decon$Sample]
write_csv(decon |> mutate(across(where(is.numeric), ~ round(.x, 3))), TAB("deconvolution_epidish.csv"))
print(decon |> select(Sample, Epi_1st, Fib_1st, IC_1st, PC1_fig3A) |> mutate(across(where(is.numeric), ~ round(.x, 2))))
r_pc1 <- cor(decon$PC1_fig3A, decon$Epi_1st)
cat(sprintf("Fig. 3A PC1 vs epithelial fraction: r = %.2f\n", r_pc1))

## ================================================================== 2. signature
## limma on M-values, Normal vs Tumor, adjusted for epithelial + immune fraction
## and array plate. CpGs whose difference is explained by tissue composition
## (epithelium vs stroma) are therefore not selected as 'tumour' CpGs.
SIG_FILE <- file.path(EXT_DIR, "gse_signature.rds")   # per-CpG statistics, reused by later scripts
{
  plate <- ifelse(pd$Sample_Plate == "", "none", pd$Sample_Plate)
  fit_contrast <- function(groups) {
    keep <- pd$Group %in% groups
    g <- droplevels(pd$Group[keep])
    d <- model.matrix(~ g + pd$Epi[keep] + pd$IC[keep] + plate[keep])
    M <- beta_to_m(Bg[, keep])
    f <- eBayes(lmFit(M, d))
    tt <- topTable(f, coef = 2, number = Inf, sort.by = "none")
    tibble(cpg = rownames(tt), t = tt$t, p = tt$P.Value, fdr = tt$adj.P.Val)
  }
  tvn <- fit_contrast(c("Normal", "Tumor"))
  ovn <- fit_contrast(c("Normal", "OPL"))
  mu <- tibble(cpg = common,
               muN = rowMeans(Bg[, pd$Group == "Normal"]),
               muO = rowMeans(Bg[, pd$Group == "OPL"]),
               muT = rowMeans(Bg[, pd$Group == "Tumor"]))
  sig <- mu |> left_join(tvn |> dplyr::rename(t_TvN = t, p_TvN = p, fdr_TvN = fdr), by = "cpg") |>
    left_join(ovn |> dplyr::rename(t_OvN = t, p_OvN = p, fdr_OvN = fdr), by = "cpg") |>
    mutate(d_TvN = muT - muN, d_OvN = muO - muN)
  saveRDS(sig, SIG_FILE)
}

N_EACH <- 2000
tum_hyper <- sig |> filter(fdr_TvN < 1e-6, d_TvN >  0.2) |> slice_max(t_TvN, n = N_EACH) |> pull(cpg)
tum_hypo  <- sig |> filter(fdr_TvN < 1e-6, d_TvN < -0.2) |> slice_min(t_TvN, n = N_EACH) |> pull(cpg)
cat(sprintf("Tumour signature: %d hyper + %d hypo CpGs (composition- and plate-adjusted)\n",
            length(tum_hyper), length(tum_hypo)))
write_csv(sig |> filter(cpg %in% c(tum_hyper, tum_hypo)) |>
            mutate(direction = ifelse(d_TvN > 0, "hyper", "hypo")),
          TAB("gse204943_tumor_signature_cpgs.csv"))

## ================================================================== 3. score
## Per CpG, position between the GSE normal mean (0) and tumour mean (1);
## sample score = median over CpGs. Hyper and hypo arms are scored separately
## so they cannot cancel (the flaw in the original 'field defect' mean-beta).
score_fun <- function(B, cpgs) {
  s <- sig[match(cpgs, sig$cpg), ]
  z <- (B[cpgs, , drop = FALSE] - s$muN) / (s$muT - s$muN)
  apply(z, 2, median, na.rm = TRUE)
}
score_tab <- function(B, ids) tibble(Sample = ids,
  hyper = score_fun(B, tum_hyper)[ids], hypo = score_fun(B, tum_hypo)[ids]) |>
  mutate(tumor_likeness = (hyper + hypo) / 2)
sc_g <- score_tab(Bg, colnames(Bg)) |> mutate(Group = pd$Group, Epi = pd$Epi, cohort = "GSE204943")
sc_p <- score_tab(Bp, ORD) |> mutate(Group = ss$Group, Epi = decon$Epi_1st[match(ORD, decon$Sample)],
                                     cohort = "Pilot")
## percentile of each pilot sample within the GSE204943 OPL distribution
## (OPL samples were not used to choose the CpGs, so this is out-of-sample)
opl_scores <- sc_g$tumor_likeness[sc_g$Group == "OPL"]
sc_p <- sc_p |> mutate(pct_of_GSE_OPL = round(100 * ecdf(opl_scores)(tumor_likeness)),
                       histology = HIST[Sample])
sc_p$delta_vs_own_adj <- NA_real_
for (s in names(own_adj)) sc_p$delta_vs_own_adj[sc_p$Sample == s] <-
  sc_p$tumor_likeness[sc_p$Sample == s] - sc_p$tumor_likeness[sc_p$Sample == own_adj[s]]
## Purity-scaled score. Among GSE204943 tumours the raw score rises with
## epithelial fraction (more tumour cells, more signal). Dividing by the
## tumour trend at the sample's own epithelial fraction gives the fraction of
## the full tumour program carried by the epithelium that is present.
trend_T <- loess(tumor_likeness ~ Epi, data = sc_g |> filter(Group == "Tumor"), span = 0.9,
                 control = loess.control(surface = "direct"))
scale_fun <- function(d) d$tumor_likeness / pmax(predict(trend_T, newdata = data.frame(Epi = d$Epi)), 0.15)
sc_g$purity_scaled <- scale_fun(sc_g); sc_p$purity_scaled <- scale_fun(sc_p)
sc_p$purity_scaled_pct_GSE_OPL <- round(100 * ecdf(sc_g$purity_scaled[sc_g$Group == "OPL"])(sc_p$purity_scaled))
r_epi_norm <- with(sc_g |> filter(Group == "Normal"), cor(tumor_likeness, Epi, method = "spearman"))
cat(sprintf("Spearman(score, Epi) within GSE normals = %.2f\n", r_epi_norm))
write_csv(sc_p |> mutate(across(where(is.numeric), ~ round(.x, 3))), TAB("tumor_likeness_pilot.csv"))
write_csv(sc_g |> mutate(across(where(is.numeric), ~ round(.x, 3))), TAB("tumor_likeness_gse204943.csv"))
cat("\nGSE204943 group medians:\n")
print(sc_g |> group_by(Group) |> summarise(n = n(), across(c(hyper, hypo, tumor_likeness), ~ round(median(.x), 3))))
cat("\nPilot:\n"); print(sc_p |> select(Sample, hyper, hypo, tumor_likeness, delta_vs_own_adj, pct_of_GSE_OPL, Epi) |>
                           mutate(across(where(is.numeric), ~ round(.x, 3))))

## is the score just epithelial content? check in GSE normals+OPL and overall
r_epi_g <- with(sc_g, cor(tumor_likeness, Epi, method = "spearman"))
cat(sprintf("GSE204943: Spearman(tumour-likeness, Epi) = %.2f\n", r_epi_g))

g_cont <- ggplot() +
  geom_boxplot(data = sc_g, aes(tumor_likeness, Group, fill = Group), alpha = 0.35, outlier.shape = NA, width = 0.6) +
  geom_jitter(data = sc_g, aes(tumor_likeness, Group, colour = Group), height = 0.15, size = 1, alpha = 0.6) +
  geom_point(data = sc_p |> mutate(y = "Pilot"), aes(tumor_likeness, y), size = 3,
             shape = 21, fill = ifelse(sc_p$Group == "Normal", "white", "#2d2d2d")) +
  geom_text_repel(data = sc_p |> mutate(y = "Pilot"), aes(tumor_likeness, y, label = Sample),
                  size = 2.8, direction = "y", nudge_y = 0.35, segment.size = 0.2, max.overlaps = 30) +
  scale_fill_manual(values = grp_cols, guide = "none") + scale_colour_manual(values = grp_cols, guide = "none") +
  scale_y_discrete(limits = c("Pilot", "Normal", "OPL", "Tumor")) +
  labs(x = "Tumour-likeness (0 = GSE204943 normal mean, 1 = tumour mean)", y = NULL,
       title = "Normal-to-tumour methylation continuum, anchored on an external cohort",
       subtitle = sprintf("%d hyper + %d hypo CpGs, composition- and plate-adjusted (GSE204943); open = adjacent tissue",
                          length(tum_hyper), length(tum_hypo))) +
  theme_bw(base_size = 10)
g_epi <- ggplot() +
  geom_point(data = sc_g, aes(Epi, tumor_likeness, colour = Group), size = 1.3, alpha = 0.6) +
  geom_point(data = sc_p, aes(Epi, tumor_likeness), size = 2.6, shape = 21,
             fill = ifelse(sc_p$Group == "Normal", "white", "#2d2d2d")) +
  geom_text_repel(data = sc_p, aes(Epi, tumor_likeness, label = Sample), size = 2.6, max.overlaps = 30) +
  scale_colour_manual(values = grp_cols) +
  labs(x = "Epithelial fraction (EpiDISH)", y = "Tumour-likeness", colour = NULL,
       title = "Score vs tissue composition") + theme_bw(base_size = 10)
ggsave(FIG("continuum_external.png"), g_cont / g_epi + plot_layout(heights = c(1, 1.1)),
       width = 8.5, height = 9, dpi = 200)

## ================================================================== 4. PCA
vg <- matrixStats::rowVars(Bg)
top <- order(vg, decreasing = TRUE)[1:5000]
pca_g <- prcomp(t(Bg[top, ]), center = TRUE)
proj_p <- predict(pca_g, t(Bp[top, ]))
pct <- round(100 * pca_g$sdev^2 / sum(pca_g$sdev^2), 1)
pcg <- tibble(PC1 = pca_g$x[, 1], PC2 = pca_g$x[, 2], Group = pd$Group)
pcp <- tibble(Sample = ORD, PC1 = proj_p[, 1], PC2 = proj_p[, 2], Group = ss$Group)
arrows <- tibble(from = unname(own_adj), to = names(own_adj)) |>
  left_join(pcp |> select(from = Sample, x = PC1, y = PC2), by = "from") |>
  left_join(pcp |> select(to = Sample, xend = PC1, yend = PC2), by = "to")
g_pca <- ggplot() +
  geom_point(data = pcg, aes(PC1, PC2, colour = Group), size = 1.4, alpha = 0.55) +
  stat_ellipse(data = pcg, aes(PC1, PC2, colour = Group), level = 0.8, linewidth = 0.4) +
  geom_segment(data = arrows, aes(x, y, xend = xend, yend = yend), arrow = arrow(length = unit(2, "mm")),
               colour = "grey25", linewidth = 0.4) +
  geom_point(data = pcp, aes(PC1, PC2), size = 2.8, shape = 21,
             fill = ifelse(pcp$Group == "Normal", "white", "#2d2d2d")) +
  geom_text_repel(data = pcp, aes(PC1, PC2, label = Sample), size = 2.8, max.overlaps = 30) +
  scale_colour_manual(values = grp_cols) +
  labs(x = sprintf("PC1 (%s%%)", pct[1]), y = sprintf("PC2 (%s%%)", pct[2]), colour = "GSE204943",
       title = "Pilot samples projected into the GSE204943 methylation space",
       subtitle = "PCA fit on GSE204943 only (top 5,000 variable shared CpGs); arrows: adjacent -> lesion") +
  theme_bw(base_size = 10)
ggsave(FIG("pca_projection_gse204943.png"), g_pca, width = 8, height = 6, dpi = 200)
write_csv(pcp |> mutate(across(where(is.numeric), ~ round(.x, 2))), TAB("pca_projection_pilot.csv"))

## ================================================================== 5. replication
## (a) Per patient: correlation of the pilot lesion-minus-own-adjacent change with
##     the external tumour-minus-normal change, across the signature CpGs,
##     and the fraction moving in the external direction by > 0.05.
sigset <- c(tum_hyper, tum_hypo)
ext_d <- sig$d_TvN[match(sigset, sig$cpg)]
rep_tab <- map_dfr(names(own_adj), function(s) {
  d <- Bp[sigset, s] - Bp[sigset, own_adj[s]]
  ## background: same number of random shared CpGs
  set.seed(1)
  bg <- sample(setdiff(common, sigset), length(sigset))
  dbg <- Bp[bg, s] - Bp[bg, own_adj[s]]
  ext_bg <- sig$d_TvN[match(bg, sig$cpg)]
  tibble(Sample = s, n_cpg = length(sigset),
         r_with_external_TvN = cor(d, ext_d),
         frac_concordant_gt0.05 = mean(sign(d) == sign(ext_d) & abs(d) > 0.05),
         frac_discordant_gt0.05 = mean(sign(d) != sign(ext_d) & abs(d) > 0.05),
         r_background = cor(dbg, ext_bg),
         median_shift_hyper = median(d[seq_along(tum_hyper)]),
         median_shift_hypo = median(d[-seq_along(tum_hyper)]))
})
write_csv(rep_tab |> mutate(across(where(is.numeric), ~ round(.x, 3))), TAB("replication_vs_gse204943.csv"))
cat("\nPilot lesion-vs-adjacent change vs external tumour program:\n")
print(rep_tab |> mutate(across(where(is.numeric), ~ round(.x, 3))), width = 200)

## (b) Pilot 'consistent' CpGs: same direction, |delta| > 0.1 in all 3 OPL pairs.
##     Replaces DMRcate-at-FDR-0.5 as the pilot's honest discovery list.
dl <- sapply(c("P3-Lesion", "P4-Lesion", "P5-Lesion"), function(s) Bp[, s] - Bp[, own_adj[s]])
cons_up <- rownames(dl)[apply(dl, 1, function(x) all(x > 0.1))]
cons_dn <- rownames(dl)[apply(dl, 1, function(x) all(x < -0.1))]
## null: how many 'consistent' CpGs arise if each patient's lesion/adjacent labels are randomly flipped
flip_counts <- map_int(1:7, function(k) {
  sgn <- as.integer(intToBits(k)[1:3]) * 2 - 1
  if (all(sgn == 1)) return(NA_integer_)
  x <- sweep(dl, 2, sgn, `*`)
  sum(apply(x, 1, function(v) all(v > 0.1)) | apply(x, 1, function(v) all(v < -0.1)))
})
cons <- tibble(cpg = c(cons_up, cons_dn), direction = rep(c("hyper", "hypo"), c(length(cons_up), length(cons_dn)))) |>
  left_join(sig |> select(cpg, d_OvN, p_OvN, fdr_OvN, d_TvN, p_TvN, fdr_TvN), by = "cpg") |>
  mutate(rep_OPL = sign(d_OvN) == ifelse(direction == "hyper", 1, -1) & p_OvN < 0.05,
         rep_Tumor = sign(d_TvN) == ifelse(direction == "hyper", 1, -1) & p_TvN < 0.05)
bg_rep <- sig |> filter(!cpg %in% cons$cpg) |>
  summarise(OPL_hyper = mean(d_OvN > 0 & p_OvN < 0.05), OPL_hypo = mean(d_OvN < 0 & p_OvN < 0.05),
            T_hyper = mean(d_TvN > 0 & p_TvN < 0.05), T_hypo = mean(d_TvN < 0 & p_TvN < 0.05))
cons_sum <- cons |> group_by(direction) |>
  summarise(n = n(), replicated_in_GSE_OPL = mean(rep_OPL), replicated_in_GSE_Tumor = mean(rep_Tumor)) |>
  mutate(background_OPL = ifelse(direction == "hyper", bg_rep$OPL_hyper, bg_rep$OPL_hypo),
         background_Tumor = ifelse(direction == "hyper", bg_rep$T_hyper, bg_rep$T_hypo))
anno_v2 <- minfi::getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
anno_v2$cpg <- sub("_[^_]+$", "", anno_v2$Name)
anno_v2 <- anno_v2[!duplicated(anno_v2$cpg), ]
cons <- cons |> mutate(gene = anno_v2$UCSC_RefGene_Name[match(cpg, anno_v2$cpg)],
                       gene = map_chr(strsplit(coalesce(gene, ""), ";"), ~ paste(unique(.x), collapse = ";")),
                       island = anno_v2$Relation_to_Island[match(cpg, anno_v2$cpg)])
write_csv(cons, TAB("pilot_consistent_cpgs_replication.csv"))
write_csv(cons_sum |> mutate(across(where(is.numeric), ~ round(.x, 3))), TAB("pilot_consistent_cpgs_summary.csv"))
cat(sprintf("\nConsistent CpGs (all 3 OPL pairs, |d|>0.1): %d hyper, %d hypo; label-flip null: %s\n",
            length(cons_up), length(cons_dn), paste(na.omit(flip_counts), collapse = ", ")))
print(cons_sum |> mutate(across(where(is.numeric), ~ round(.x, 3))))

## ================================================================== 6. PMD + ZNF471
pmd <- read_csv(file.path(GSE_DIR, "results/hypo_PMD_blocks.csv"), show_col_types = FALSE)
pmd_gr <- GRanges(pmd$seqnames, IRanges(pmd$start, pmd$end))
ag <- ann_g[intersect(common, rownames(ann_g)), ]
os <- ag$Relation_to_Island == "OpenSea"
cg_gr <- GRanges(ag$chr, IRanges(ag$pos, width = 1))     # hg19, same build as blocks
in_pmd <- rownames(ag)[overlapsAny(cg_gr, pmd_gr) & os]
out_pmd <- rownames(ag)[!overlapsAny(cg_gr, pmd_gr) & os]
pmd_score <- function(B) colMeans(B[in_pmd, ]) - colMeans(B[out_pmd, ])   # within-sample contrast
pmd_g <- tibble(Sample = colnames(Bg), Group = pd$Group, pmd_contrast = pmd_score(Bg),
                pmd_mean = colMeans(Bg[in_pmd, ]))
gse_pmd <- read_csv(file.path(GSE_DIR, "results/PMD_erosion_score_per_sample.csv"), show_col_types = FALSE)
r_pmd <- cor(pmd_g$pmd_contrast, gse_pmd$PMD_erosion_score[match(pmd_g$Sample, gse_pmd$Sample)],
             use = "complete.obs", method = "spearman")
cat(sprintf("\nPMD contrast vs published GSE204943 PMD_erosion_score: Spearman %.2f (n CpGs in PMD %d)\n",
            r_pmd, length(in_pmd)))
pmd_p <- tibble(Sample = ORD, Group = ss$Group, pmd_contrast = pmd_score(Bp)[ORD], pmd_mean = colMeans(Bp[in_pmd, ORD]))

znf_cpgs <- read_csv(file.path(GSE_DIR, "results/ZNF471_CpG_details.csv"), show_col_types = FALSE) |>
  filter(Region_Class == "Promoter") |> pull(Name) |> intersect(common)
znf_g <- colMeans(Bg[znf_cpgs, ]); znf_p <- colMeans(Bp[znf_cpgs, ORD])
pz_p <- pmd_p |> mutate(ZNF471_promoter = znf_p,
                        pmd_pct_GSE_OPL = round(100 * ecdf(pmd_g$pmd_contrast[pmd_g$Group == "OPL"])(pmd_contrast)),
                        znf_pct_GSE_OPL = round(100 * ecdf(znf_g[pd$Group == "OPL"])(ZNF471_promoter)))
pz_g <- pmd_g |> mutate(ZNF471_promoter = znf_g)
write_csv(pz_p |> mutate(across(where(is.numeric), ~ round(.x, 4))), TAB("pmd_znf471_pilot.csv"))
cat("\nGSE204943 group medians (PMD contrast, ZNF471):\n")
print(pz_g |> group_by(Group) |> summarise(pmd = round(median(pmd_contrast), 4), znf = round(median(ZNF471_promoter), 3)))
print(pz_p |> mutate(across(where(is.numeric), ~ round(.x, 4))))

g_pz <- ggplot() +
  geom_point(data = pz_g, aes(pmd_contrast, ZNF471_promoter, colour = Group), size = 1.3, alpha = 0.6) +
  geom_point(data = pz_p, aes(pmd_contrast, ZNF471_promoter), size = 2.6, shape = 21,
             fill = ifelse(pz_p$Group == "Normal", "white", "#2d2d2d")) +
  geom_text_repel(data = pz_p, aes(pmd_contrast, ZNF471_promoter, label = Sample), size = 2.6) +
  scale_colour_manual(values = grp_cols) +
  labs(x = "PMD erosion (mean beta in tumour-hypomethylated PMD blocks minus outside; open sea)",
       y = "ZNF471 promoter mean beta", colour = "GSE204943",
       title = "The two GSE204943 programs, with the pilot overlaid") + theme_bw(base_size = 10)
ggsave(FIG("pmd_znf471.png"), g_pz, width = 7.5, height = 5.5, dpi = 200)

## ================================================================== 7. genes
gene_lists <- strsplit(coalesce(anno_v2$UCSC_RefGene_Name, ""), ";")
## CGI/shore TSS-proximal CpGs; genes without a CGI promoter (LAMC2, COL17A1,
## CASP8, FAP, CD163, ...) fall back to all TSS200/TSS1500 CpGs
prom_cpgs <- function(gene) {
  hit <- vapply(gene_lists, function(x) gene %in% x, logical(1))
  grp <- grepl("TSS200|TSS1500|1stExon|5'UTR", anno_v2$UCSC_RefGene_Group)
  cgi <- intersect(anno_v2$cpg[hit & grp & anno_v2$Relation_to_Island %in% c("Island", "N_Shore", "S_Shore")], common)
  if (length(cgi) > 0) return(structure(cgi, type = "CGI"))
  tss <- intersect(anno_v2$cpg[hit & grepl("TSS200|TSS1500", anno_v2$UCSC_RefGene_Group)], common)
  structure(tss, type = "TSS non-CGI")
}
genes <- c(
  ## CNA / tumour suppressors
  "CDKN2A", "CDKN2B", "MTAP", "FHIT", "MGMT", "RASSF1", "TP53",
  ## genes reported mutated
  "CASP8", "FAT1", "NOTCH1", "PIK3CA", "KMT2D",
  ## manuscript transcriptional programs
  "LAMC2", "COL17A1", "ITGA6", "ITGB4", "PTHLH", "HMGA2", "ITGB6", "MAGEA11",
  "LRRC15", "FAP", "FOLR2", "CD163", "MRC1", "PI16",
  ## top pilot DMR genes / GSE204943 program
  "ZNF382", "ZNF529", "RELN", "NELL1", "HOXA5", "HOXC9", "ZNF471", "TFAP2A")
gene_tab <- map_dfr(genes, function(g) {
  cps <- prom_cpgs(g)
  if (length(cps) == 0) return(tibble(gene = g, n_cpg = 0L))
  bp <- colMeans(Bp[cps, , drop = FALSE])
  s <- sig[match(cps, sig$cpg), ]
  tibble(gene = g, n_cpg = length(cps), cpg_type = attr(cps, "type"),
         GSE_Normal = mean(s$muN), GSE_OPL = mean(s$muO), GSE_Tumor = mean(s$muT),
         GSE_TvN_minP = min(s$p_TvN), GSE_OvN_minP = min(s$p_OvN),
         !!!as.list(bp[ORD]),
         d_P2 = bp["P2-Tumor"] - bp["P2-NAT"], d_P3 = bp["P3-Lesion"] - bp["P3-Adj"],
         d_P4 = bp["P4-Lesion"] - bp["P4-Adj"], d_P5 = bp["P5-Lesion"] - bp["P5-Adj"])
})
write_csv(gene_tab |> mutate(across(where(is.numeric), ~ signif(.x, 3))), TAB("candidate_gene_promoter_methylation.csv"))
cat("\nPromoter (CGI/shore, TSS-proximal) methylation, candidate genes:\n")
print(gene_tab |> select(gene, n_cpg, cpg_type, GSE_Normal, GSE_Tumor, d_P2, d_P3, d_P4, d_P5) |>
        mutate(across(where(is.double), ~ round(.x, 3))), n = 60)

saveRDS(list(decon = decon, sc_p = sc_p, sc_g = sc_g, rep_tab = rep_tab, cons = cons, cons_sum = cons_sum,
             pz_p = pz_p, pz_g = pz_g, gene_tab = gene_tab, pcp = pcp, pct = pct,
             tum_hyper = tum_hyper, tum_hypo = tum_hypo, flip_counts = flip_counts,
             conc_r = conc_r, n_common = length(common), r_pc1 = r_pc1, r_epi_g = r_epi_g, r_epi_norm = r_epi_norm, r_pmd = r_pmd),
        file.path(EXT_DIR, "external_methylation.rds"))
cat("done\n")
