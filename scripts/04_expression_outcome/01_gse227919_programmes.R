## External transcriptomic test of the manuscript's programs and of the
## methylation findings, in GSE227919 (oral biopsies, n = 66: Control 18 [benign
## mucosa or non-dysplastic lesions such as fibroma, varix or non-dysplastic
## leukoplakia; not healthy volunteers], HkNR 17, Dysplasia 22 [any grade],
## Cancer 9; progression follow-up for 25 non-cancer biopsies: 15 stable,
## 10 progressed, 3 of them controls).
##
## The supplied "raw" matrix holds non-integer normalised counts (0 to ~650k),
## so it is log2(x + 1) transformed (as the earlier OPMD_methyl pipeline did)
## and analysed with limma-trend.
##
## Run from the repository root:  Rscript scripts/04_expression_outcome/01_gse227919_programmes.R  (or make outcome)

source("R/lib/config.R")   # paths and helpers (TAB, FIG, OUT_DIR, GSE_DIR, ...)
suppressPackageStartupMessages({
  library(tidyverse); library(limma); library(patchwork); library(ggrepel)
  library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
})
cls_cols <- c(Control = "#7f8c8d", HkNR = "#3b6fb6", Dysplasia = "#e08a1e", Cancer = "#c0392b")

## ------------------------------------------------------------------ data
## whitespace-delimited; header has one field fewer than rows (gene IDs = row names)
ex <- read.table(gzfile(file.path(RNA_DIR, "GSE227919_2023_20_03_expr_mat_raw.txt.gz")),
                 header = TRUE, check.names = FALSE)
ann <- read.csv(gzfile(file.path(RNA_DIR, "GSE227919_2023_21_03_pml_annotations.csv.gz")),
                check.names = FALSE, na.strings = c("#N/A", "NA", ""))
names(ann)[1] <- "Sample"
ann <- ann[match(colnames(ex), ann$Sample), ]
stopifnot(identical(ann$Sample, colnames(ex)))
ann$Class <- factor(ann$Class, levels = c("Control", "HkNR", "Dysplasia", "Cancer"))
ann$Progressed <- ifelse(is.na(ann$Progression_status), NA, ann$Progression_status != "Stable")
cat(sprintf("Expression: %d genes x %d samples; value range %.1f-%.1f\n",
            nrow(ex), ncol(ex), min(ex), max(ex)))
ex <- log2(as.matrix(ex) + 1)
ex <- ex[rowMeans(ex > 1) >= 0.2, ]                  # expressed in >= 20% of samples
cat(sprintf("After log2(x+1) and expression filter: %d genes\n", nrow(ex)))
z <- t(scale(t(ex)))

## ------------------------------------------------------------------ gene sets
ext <- readRDS(file.path(EXT_DIR, "external_methylation.rds"))
anno_v2 <- minfi::getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
anno_v2$cpg <- sub("_[^_]+$", "", anno_v2$Name)
anno_v2 <- anno_v2[!duplicated(anno_v2$cpg), ]
prom_genes <- function(cpgs, cgi_only = TRUE) {
  a <- anno_v2[match(cpgs, anno_v2$cpg), ]
  keep <- grepl("TSS200|TSS1500|1stExon", a$UCSC_RefGene_Group) &
    (!cgi_only | a$Relation_to_Island == "Island")
  unique(unlist(strsplit(a$UCSC_RefGene_Name[keep], ";")))
}
## Methylation-defined sets from the external methylation cohort (GSE204943),
## restricted to genes whose promoter CpGs also moved in the same direction
## in >= 2 of the 3 pilot pairs that follow the tumour program (P2, P4, P5).
sig <- readRDS(file.path(EXT_DIR, "gse_signature.rds"))
p <- readRDS(file.path(OUT_DIR, "processed.rds"))
nm <- setNames(p$sample_sheet$Sample_Name, paste0(p$sample_sheet$Patient, p$sample_sheet$Group))
dpil <- cbind(P2 = p$betas[, nm["Patient 2Tumor"]]  - p$betas[, nm["Patient 2Normal"]],
              P4 = p$betas[, nm["Patient 4Lesion"]] - p$betas[, nm["Patient 4Normal"]],
              P5 = p$betas[, nm["Patient 5Lesion"]] - p$betas[, nm["Patient 5Normal"]],
              P3 = p$betas[, nm["Patient 3Lesion"]] - p$betas[, nm["Patient 3Normal"]])
hyp_c <- intersect(ext$tum_hyper, rownames(dpil)); hypo_c <- intersect(ext$tum_hypo, rownames(dpil))
hyp_pil  <- hyp_c[rowSums(dpil[hyp_c, 1:3] > 0.1) >= 2]
hypo_pil <- hypo_c[rowSums(dpil[hypo_c, 1:3] < -0.1) >= 2]
sets <- list(
  Invasive_epithelial = c("LAMC2", "COL17A1", "PTHLH", "ITGA6", "ITGB4", "HMGA2"),
  Dysplastic_epithelial = c("MAGEA11", "ITGB6", "LAMC2", "TNFRSF12A", "HBEGF", "ULBP2"),
  Cycling_epithelial = c("HMMR", "ASPM", "BUB1", "NEIL3", "BUB1B", "PBK"),
  LRRC15_CAF = c("LRRC15", "FAP", "COL11A1"),
  PI16_fibroblast = c("PI16"),
  FOLR2_TAM = c("FOLR2", "CD163", "MRC1"),
  Treg = c("FOXP3", "IL2RA", "CTLA4", "TNFRSF18"),
  mregDC = c("LAMP3", "CCR7", "FSCN1", "CD274"),
  Meth_silenced_promoters = prom_genes(hyp_pil),       # expected DOWN with progression
  Meth_activated_promoters = prom_genes(hypo_pil, cgi_only = FALSE),  # expected UP
  KZFP_chr19_hyper = c("ZNF382", "ZNF529", "ZNF471"),
  p16_9p21 = c("CDKN2A", "CDKN2B", "MTAP")
)
sets <- lapply(sets, intersect, rownames(z))
cat("Gene-set sizes in GSE227919:\n"); print(sapply(sets, length))
writeLines(c(sprintf("Meth_silenced_promoters (%d): %s", length(sets$Meth_silenced_promoters),
                     paste(sort(sets$Meth_silenced_promoters), collapse = ", ")),
             sprintf("Meth_activated_promoters (%d): %s", length(sets$Meth_activated_promoters),
                     paste(sort(sets$Meth_activated_promoters), collapse = ", "))),
           TAB("rnaseq_methylation_gene_sets.txt"))

scores <- sapply(sets, function(g) colMeans(z[g, , drop = FALSE]))
sdat <- bind_cols(ann |> select(Sample, Class, Progression_status, Progressed, Age, Sex,
                                Smoking = imputed_smoking_label), as_tibble(scores))

## ------------------------------------------------------------------ tests
## (1) trend across Control -> HkNR -> Dysplasia -> Cancer
## (2) progressed vs stable among non-cancer lesions with follow-up, with
##     AUC and a logistic model adjusted for histologic class (does the
##     program add information beyond the histology?)
auc <- function(x, y) { r <- rank(x); n1 <- sum(y); n0 <- sum(!y)
  (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
prog <- sdat |> filter(Class != "Cancer", !is.na(Progressed))
cat(sprintf("Progression subset: %d lesions (%d progressed)\n", nrow(prog), sum(prog$Progressed)))
res <- map_dfr(names(sets), function(s) {
  trend <- cor.test(sdat[[s]], as.integer(sdat$Class), method = "spearman", exact = FALSE)
  w <- wilcox.test(prog[[s]][prog$Progressed], prog[[s]][!prog$Progressed], exact = FALSE)
  set.seed(7)
  boots <- replicate(2000, { i <- sample(nrow(prog), replace = TRUE)
    if (length(unique(prog$Progressed[i])) < 2) NA else auc(prog[[s]][i], prog$Progressed[i]) })
  glm_adj <- glm(Progressed ~ prog[[s]] + Class, data = prog, family = binomial)
  tibble(set = s, n_genes = length(sets[[s]]),
         trend_rho = unname(trend$estimate), trend_p = trend$p.value,
         Control = mean(sdat[[s]][sdat$Class == "Control"]), HkNR = mean(sdat[[s]][sdat$Class == "HkNR"]),
         Dysplasia = mean(sdat[[s]][sdat$Class == "Dysplasia"]), Cancer = mean(sdat[[s]][sdat$Class == "Cancer"]),
         prog_AUC = auc(prog[[s]], prog$Progressed),
         prog_AUC_lo = quantile(boots, 0.025, na.rm = TRUE), prog_AUC_hi = quantile(boots, 0.975, na.rm = TRUE),
         prog_wilcox_p = w$p.value,
         prog_OR_per_SD_adj_class = exp(coef(glm_adj)[2] * sd(prog[[s]])),
         prog_p_adj_class = summary(glm_adj)$coefficients[2, 4])
}) |> mutate(trend_fdr = p.adjust(trend_p, "BH"), prog_fdr = p.adjust(prog_wilcox_p, "BH"))
write_csv(res |> mutate(across(where(is.numeric), ~ signif(.x, 3))), TAB("rnaseq_program_tests.csv"))
print(res |> select(set, n_genes, trend_rho, trend_fdr, prog_AUC, prog_AUC_lo, prog_AUC_hi, prog_wilcox_p, prog_p_adj_class) |>
        mutate(across(where(is.numeric), ~ signif(.x, 2))), width = 200)

## histology alone as a predictor, for reference
hist_auc <- auc(as.integer(prog$Class), prog$Progressed)
cat(sprintf("Histologic class alone (ordinal) AUC for progression: %.2f\n", hist_auc))

## per-gene table (programs + CNA genes)
genes_all <- unique(c(unlist(sets[c("Invasive_epithelial", "Dysplastic_epithelial", "LRRC15_CAF",
                                    "PI16_fibroblast", "FOLR2_TAM", "Treg", "mregDC",
                                    "KZFP_chr19_hyper", "p16_9p21")]), "NELL1", "HOXA5", "FHIT"))
genes_all <- intersect(genes_all, rownames(ex))
d <- model.matrix(~ 0 + Class, data = ann); colnames(d) <- levels(ann$Class)
fit <- eBayes(contrasts.fit(lmFit(ex, d), makeContrasts(Dys = Dysplasia - Control, Ca = Cancer - Control, levels = d)), trend = TRUE)
gene_tab <- tibble(gene = genes_all,
                   logFC_Dysplasia_vs_Control = fit$coefficients[genes_all, "Dys"],
                   p_Dys = fit$p.value[genes_all, "Dys"],
                   logFC_Cancer_vs_Control = fit$coefficients[genes_all, "Ca"],
                   p_Ca = fit$p.value[genes_all, "Ca"],
                   prog_AUC = sapply(genes_all, function(g) auc(ex[g, prog$Sample], prog$Progressed)),
                   prog_p = sapply(genes_all, function(g) wilcox.test(ex[g, prog$Sample][prog$Progressed],
                                                                      ex[g, prog$Sample][!prog$Progressed], exact = FALSE)$p.value))
write_csv(gene_tab |> mutate(across(where(is.numeric), ~ signif(.x, 3))), TAB("rnaseq_gene_tests.csv"))
print(gene_tab |> mutate(across(where(is.numeric), ~ signif(.x, 2))), n = 60)

## ------------------------------------------------------------------ figure
show <- c("Invasive_epithelial", "Dysplastic_epithelial", "LRRC15_CAF", "FOLR2_TAM", "Treg",
          "Meth_silenced_promoters", "Meth_activated_promoters", "KZFP_chr19_hyper", "p16_9p21")
lg <- sdat |> pivot_longer(all_of(show), names_to = "set", values_to = "score") |>
  mutate(set = factor(set, show))
g1 <- ggplot(lg, aes(Class, score, fill = Class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) + geom_jitter(width = 0.15, size = 0.7, alpha = 0.6) +
  facet_wrap(~set, scales = "free_y", ncol = 3) + scale_fill_manual(values = cls_cols, guide = "none") +
  labs(x = NULL, y = "Program score (mean gene z)", title = "GSE227919: programs across histologic classes") +
  theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
lp <- prog |> pivot_longer(all_of(show), names_to = "set", values_to = "score") |>
  mutate(set = factor(set, show), Outcome = ifelse(Progressed, "Progressed", "Stable")) |>
  left_join(res |> select(set, prog_AUC, prog_wilcox_p), by = "set") |>
  mutate(lab = sprintf("%s\nAUC %.2f, p = %.2g", set, prog_AUC, prog_wilcox_p))
g2 <- ggplot(lp, aes(Outcome, score)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90") +
  geom_jitter(aes(colour = Class), width = 0.15, size = 1.2) +
  facet_wrap(~lab, scales = "free_y", ncol = 3) + scale_colour_manual(values = cls_cols) +
  labs(x = NULL, y = "Program score", colour = "Histology",
       title = sprintf("GSE227919: progressed (n = %d) vs stable (n = %d) non-cancer lesions",
                       sum(prog$Progressed), sum(!prog$Progressed))) +
  theme_bw(base_size = 9)
ggsave(FIG("rnaseq_programs_by_class.png"), g1, width = 9, height = 8, dpi = 200)
ggsave(FIG("rnaseq_programs_progression.png"), g2, width = 9, height = 8.5, dpi = 200)
saveRDS(list(res = res, gene_tab = gene_tab, sets = sets, hist_auc = hist_auc, n_prog = nrow(prog),
             n_events = sum(prog$Progressed)), file.path(EXT_DIR, "external_rnaseq.rds"))
cat("done\n")
