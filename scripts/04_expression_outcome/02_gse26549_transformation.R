## Outcome test in GSE26549: 86 oral premalignant lesions from an MD Anderson
## retinoid chemoprevention trial, Affymetrix HuGene 1.0 ST, oral cancer-free
## survival (OCFS), 35 cancers, median follow-up ~6 years.
##
## Features were fixed before this cohort was opened (all defined elsewhere):
##   - the manuscript's transcriptional programs (pilot Xenium / Fig. 4-6)
##   - arm-level expression dosage at the oral LOH arms (3p, 9p, 13q, 4q, ...)
##   - methylation-linked gene sets (GSE204943 + pilot, from 04_expression_outcome/01_gse227919_programmes.R)
##   - an expression proxy of the OPL methylome state, trained ONLY in the
##     Mumbai OPLs (GSE85195 x GSE204943, 03_methylation_cohorts/03_mumbai_crosscheck.R) and applied here unchanged
## Cox models per SD, unadjusted and adjusted for baseline histology
## (dysplasia vs hyperplasia); C-index and likelihood-ratio test of the added
## value over histology. BH-FDR across all features.
##
## Run from the repository root:  Rscript scripts/04_expression_outcome/02_gse26549_transformation.R  (or make outcome)

source("R/lib/config.R")   # paths and helpers (TAB, FIG, OUT_DIR, GSE_DIR, ...)
suppressPackageStartupMessages({
  library(tidyverse); library(survival); library(sesameData); library(GenomicRanges)
  library(hugene10sttranscriptcluster.db)
})
## AnnotationDbi / S4Vectors mask these dplyr verbs
select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename

## ------------------------------------------------------------------ data
meta <- readLines(gzfile(file.path(XDATA, "GSE26549_series_matrix.txt.gz")))
ch <- grep("^!Sample_characteristics_ch1", meta, value = TRUE)
gsm <- gsub('"', "", strsplit(grep("^!Sample_geo_accession", meta, value = TRUE), "\t")[[1]][-1])
clin <- map(ch, ~ gsub('"', "", strsplit(.x, "\t")[[1]][-1]))
key <- map_chr(clin, ~ sub(":.*", "", .x[1]))
clin <- setNames(map(clin, ~ sub("^[^:]*: ", "", .x)), make.names(key))
clin <- as_tibble(clin) |> mutate(gsm = gsm)
names(clin) <- c("sex", "race", "alcohol", "smoking", "age", "arm", "histology", "histology_detail",
                 "p63", "podoplanin", "ocfs_years", "outcome", "biopsy_time", "gsm")
clin <- clin |> mutate(ocfs_years = as.numeric(ocfs_years), event = as.integer(outcome == "oral cancer development"),
                       dysplasia = as.integer(histology == "dysplasia"), age = as.numeric(age))
b <- grep("series_matrix_table_begin", meta); e <- grep("series_matrix_table_end", meta)
X <- as.matrix(read.delim(text = meta[(b + 1):(e - 1)], check.names = FALSE, row.names = 1))
X <- X[, clin$gsm]
cat(sprintf("GSE26549: %d probesets x %d lesions; %d cancers; median follow-up %.1f y; value range %.1f-%.1f\n",
            nrow(X), ncol(X), sum(clin$event), median(clin$ocfs_years), min(X, na.rm = TRUE), max(X, na.rm = TRUE)))
print(table(histology = clin$histology_detail, cancer = clin$event))
sym <- mapIds(hugene10sttranscriptcluster.db, rownames(X), "SYMBOL", "PROBEID", multiVals = "first")
ok <- !is.na(sym)
X <- X[ok, ]; sym <- sym[ok]
keep <- tapply(seq_along(sym), sym, function(i) i[which.max(rowMeans(X[i, , drop = FALSE]))])
X <- X[unlist(keep), ]; rownames(X) <- names(keep)
Z <- t(scale(t(X)))
cat(sprintf("Genes after symbol mapping: %d\n", nrow(X)))

## ------------------------------------------------------------------ features
rna <- readRDS(file.path(EXT_DIR, "external_rnaseq.rds"))
sets <- rna$sets[c("Invasive_epithelial", "Dysplastic_epithelial", "Cycling_epithelial", "LRRC15_CAF",
                   "FOLR2_TAM", "Treg", "mregDC", "Meth_silenced_promoters", "Meth_activated_promoters",
                   "KZFP_chr19_hyper", "p16_9p21")]
sets <- lapply(sets, intersect, rownames(Z))
feat <- as_tibble(sapply(sets, function(g) colMeans(Z[g, , drop = FALSE])))

## arm dosage (relative to the cohort median; no normal tissue in this series)
arms38 <- sesameData_getGenomeInfo("hg38")$cytoBand |> mutate(across(c(chrom, name), as.character)) |>
  filter(chrom %in% paste0("chr", 1:22)) |> mutate(arm = substr(name, 1, 1)) |>
  group_by(chrom, arm) |> summarise(start = min(chromStart), end = max(chromEnd), .groups = "drop") |>
  mutate(arm_label = paste0(sub("chr", "", chrom), arm))
tx <- sesameData_getTxnGRanges("hg38")
gpos <- tibble(gene = mcols(tx)$gene_name, chr = as.character(seqnames(tx)), pos = start(tx)) |>
  filter(chr %in% paste0("chr", 1:22), gene %in% rownames(X)) |>
  group_by(gene) |> summarise(chr = dplyr::first(chr), pos = min(pos), .groups = "drop") |>
  filter(!(chr == "chr6" & pos > 28e6 & pos < 34e6))
gpos$arm_label <- map2_chr(gpos$chr, gpos$pos, function(ch, p) {
  a <- arms38$arm_label[arms38$chrom == ch & arms38$start <= p & arms38$end > p]; if (length(a)) a[1] else NA })
expr_ok <- rowMeans(X[gpos$gene, ]) > quantile(rowMeans(X), 0.3)
gpos <- gpos[expr_ok & !is.na(gpos$arm_label), ]
zz <- pmax(pmin(Z[gpos$gene, ], 4), -4)
dos <- sapply(c("3p", "9p", "13q", "4q", "8p", "11q", "17p"), function(a) colMeans(zz[gpos$gene[gpos$arm_label == a], , drop = FALSE]))
all_arm <- sapply(unique(gpos$arm_label), function(a) colMeans(zz[gpos$gene[gpos$arm_label == a], , drop = FALSE]))
dos <- dos - rowMeans(all_arm)
feat <- feat |> mutate(loss_3p = -dos[, "3p"], loss_9p = -dos[, "9p"], loss_13q = -dos[, "13q"], loss_4q = -dos[, "4q"],
                       loss_core_3p9p = -(dos[, "3p"] + dos[, "9p"]) / 2,
                       loss_methylome_arms = -rowMeans(dos[, c("3p", "9p", "13q", "4q")]))

## methylome-state expression proxy, trained in Mumbai OPLs only (GSE85195 x GSE204943)
mum <- readRDS(file.path(EXT_DIR, "mumbai_crosscheck.rds"))
XDATA_E <- file.path(XDATA, "GSE85195_series_matrix.txt.gz")
read_series <- function(f) {
  meta <- readLines(gzfile(f))
  ttl <- strsplit(grep("^!Sample_title", meta, value = TRUE), "\t")[[1]][-1]
  gsm <- strsplit(grep("^!Sample_geo_accession", meta, value = TRUE), "\t")[[1]][-1]
  b <- grep("series_matrix_table_begin", meta); e <- grep("series_matrix_table_end", meta)
  m <- read.delim(text = meta[(b + 1):(e - 1)], check.names = FALSE, row.names = 1)
  list(m = as.matrix(m), map = tibble(gsm = gsub('"', "", gsm), id = gsub(" ", "", sub(".*_", "", gsub('"', "", ttl)))))
}
ex <- read_series(XDATA_E)
gl <- readLines(file.path(XDATA, "GPL6480.annot.txt")); h <- grep("^ID\t", gl)[1]
g6480 <- read.delim(text = gl[h:length(gl)], quote = "", check.names = FALSE)
names(g6480) <- make.unique(names(g6480))
g6480 <- g6480[g6480$GENE_SYMBOL != "" & g6480$ID != "!platform_table_end", c("ID", "GENE_SYMBOL")]
E <- ex$m[intersect(rownames(ex$m), g6480$ID), ]
es <- g6480$GENE_SYMBOL[match(rownames(E), g6480$ID)]
kk <- tapply(seq_along(es), es, function(i) i[which.max(rowMeans(E[i, , drop = FALSE], na.rm = TRUE))])
E <- E[unlist(kk), ]; rownames(E) <- names(kk); colnames(E) <- ex$map$id[match(colnames(E), ex$map$gsm)]
tr <- mum$opl_e                                    # 15 Mumbai OPLs with methylome tumour-likeness
Etr <- E[intersect(rownames(E), rownames(X)), tr$id]
Etr <- Etr[apply(Etr, 1, sd, na.rm = TRUE) > 0.3, ]
rr <- apply(Etr, 1, function(v) suppressWarnings(cor(v, tr$tumor_likeness, method = "spearman")))
up <- names(sort(rr, decreasing = TRUE))[1:50]; dn <- names(sort(rr))[1:50]
writeLines(c(paste("UP:", paste(up, collapse = ", ")), paste("DOWN:", paste(dn, collapse = ", "))),
           TAB("methylome_state_expression_proxy_genes.txt"))
feat$Methylome_state_proxy <- colMeans(Z[up, ]) - colMeans(Z[dn, ])
cat(sprintf("Methylome-state proxy: 50 up / 50 down genes from %d Mumbai OPLs (|rho| range %.2f-%.2f)\n",
            nrow(tr), min(abs(rr[c(up, dn)])), max(abs(rr[c(up, dn)]))))

## ------------------------------------------------------------------ Cox
d <- bind_cols(clin, feat)
fit_h <- coxph(Surv(ocfs_years, event) ~ dysplasia, data = d)
c_h <- concordance(fit_h)$concordance
cox_tab <- map_dfr(names(feat), function(f) {
  x <- as.numeric(scale(d[[f]]))
  m1 <- coxph(Surv(d$ocfs_years, d$event) ~ x)
  m2 <- coxph(Surv(d$ocfs_years, d$event) ~ x + d$dysplasia)
  lr <- anova(fit_h, coxph(Surv(ocfs_years, event) ~ dysplasia + x, data = mutate(d, x = x)))
  s1 <- summary(m1)$coefficients; s2 <- summary(m2)$coefficients
  tibble(feature = f, n_genes = if (f %in% names(sets)) length(sets[[f]]) else NA_integer_,
         HR_per_SD = s1[1, 2], p = s1[1, 5], C = concordance(m1)$concordance,
         HR_adj_histology = s2[1, 2], p_adj_histology = s2[1, 5],
         C_with_histology = concordance(coxph(Surv(ocfs_years, event) ~ dysplasia + x, data = mutate(d, x = x)))$concordance,
         LR_p_added_to_histology = lr$`Pr(>|Chi|)`[2])
}) |> mutate(fdr = p.adjust(p, "BH"), fdr_adj = p.adjust(p_adj_histology, "BH"))
write_csv(cox_tab |> mutate(across(where(is.numeric), ~ signif(.x, 3))), TAB("gse26549_cox.csv"))
cat(sprintf("\nHistology alone (dysplasia vs hyperplasia): HR %.2f, p = %.2g, C = %.2f\n",
            exp(coef(fit_h)), summary(fit_h)$coefficients[1, 5], c_h))
print(cox_tab |> arrange(p) |> mutate(across(where(is.numeric), ~ signif(.x, 2))), width = 220, n = 30)

## Kaplan-Meier for the leading features (median split), saved for the figure
km_feats <- cox_tab |> arrange(p) |> slice_head(n = 3) |> pull(feature)
km <- map_dfr(km_feats, function(f) {
  grp <- ifelse(d[[f]] > median(d[[f]]), "High", "Low")
  sf <- survfit(Surv(d$ocfs_years, d$event) ~ grp)
  lr <- survdiff(Surv(d$ocfs_years, d$event) ~ grp)
  tibble(feature = f, time = sf$time, surv = sf$surv, strata = rep(sub("grp=", "", names(sf$strata)), sf$strata),
         logrank_p = 1 - pchisq(lr$chisq, 1))
})
## ================================================================== robustness
cox_p <- function(score) { x <- as.numeric(scale(score))
  s <- summary(coxph(Surv(d$ocfs_years, d$event) ~ x))$coefficients; c(HR = s[1, 2], p = s[1, 5]) }
proxy_from <- function(rho, k) { u <- names(sort(rho, decreasing = TRUE))[1:k]; v <- names(sort(rho))[1:k]
  colMeans(Z[u, , drop = FALSE]) - colMeans(Z[v, , drop = FALSE]) }
obs <- cox_p(feat$Methylome_state_proxy)

## (1) multivariable model
d2 <- d |> mutate(proxy = as.numeric(scale(Methylome_state_proxy)), ever_smoker = as.integer(smoking != "never"),
                  baseline_biopsy = as.integer(grepl("baseline", biopsy_time)), arm = factor(arm))
mv <- coxph(Surv(ocfs_years, event) ~ proxy + dysplasia + ever_smoker + age + arm + baseline_biopsy, data = d2)
mv_tab <- as_tibble(summary(mv)$coefficients, rownames = "term") |> select(term, HR = `exp(coef)`, p = `Pr(>|z|)`)
cat("\nMultivariable Cox (proxy per SD):\n"); print(mv_tab |> mutate(across(where(is.numeric), ~ signif(.x, 2))))

## (2) is the proxy just dysplasia/proliferation/keratinocyte content?
assoc_proxy <- map_dfr(c(names(sets)[1:7], "loss_3p", "loss_core_3p9p"), function(f)
  tibble(feature = f, rho_with_proxy = cor(d$Methylome_state_proxy, d[[f]], method = "spearman")))
w_dys <- wilcox.test(d$Methylome_state_proxy ~ d$dysplasia)$p.value
cat(sprintf("\nProxy, dysplasia vs hyperplasia: Wilcoxon p = %.2g\n", w_dys)); print(assoc_proxy |> mutate(rho_with_proxy = round(rho_with_proxy, 2)))
adj_prog <- coxph(Surv(ocfs_years, event) ~ proxy + scale(Cycling_epithelial) + scale(Dysplastic_epithelial) + dysplasia, data = d2)
cat("Proxy adjusted for cycling + dysplastic programs + histology:\n"); print(signif(summary(adj_prog)$coefficients[, c(2, 5)], 2))

## (2b) is the proxy a biopsy-composition axis (epithelium vs stroma / immune)?
comp_sets <- list(stromal = c("DCN", "LUM", "COL1A1", "COL1A2", "COL3A1", "PDGFRB", "FBLN1"),
                  immune = c("PTPRC", "CD3E", "CD2", "CD19", "MS4A1", "CD68", "CD14"),
                  keratinocyte = c("KRT5", "KRT14", "KRT13", "KRT4", "TP63", "SFN"))
comp_sc <- sapply(comp_sets, function(g) colMeans(Z[intersect(g, rownames(Z)), , drop = FALSE]))
d2 <- bind_cols(d2, as_tibble(scale(comp_sc)))
comp_rho <- sapply(colnames(comp_sc), function(k) cor(d2$proxy, d2[[k]], method = "spearman"))
adj_comp <- coxph(Surv(ocfs_years, event) ~ proxy + stromal + immune + keratinocyte + dysplasia, data = d2)
comp_alone <- coxph(Surv(ocfs_years, event) ~ stromal + immune + keratinocyte + dysplasia, data = d2)
cat("\nProxy vs composition scores (Spearman):\n"); print(round(comp_rho, 2))
cat("Proxy adjusted for stromal, immune, keratinocyte scores + histology:\n"); print(signif(summary(adj_comp)$coefficients[, c(2, 5)], 2))
cat(sprintf("Composition scores alone: C = %.2f; with proxy: C = %.2f; LR p for adding proxy = %.2g\n",
            concordance(comp_alone)$concordance, concordance(adj_comp)$concordance, anova(comp_alone, adj_comp)$`Pr(>|Chi|)`[2]))
comp_check <- list(rho = comp_rho, adj = summary(adj_comp)$coefficients, lr_p = anova(comp_alone, adj_comp)$`Pr(>|Chi|)`[2])

## (3) size sensitivity, and up-only / down-only
pool_rho <- rr
size_tab <- map_dfr(c(10, 25, 50, 100, 200), function(k) { r <- cox_p(proxy_from(pool_rho, k)); tibble(k_each = k, HR = r["HR"], p = r["p"]) })
size_tab <- bind_rows(size_tab,
  tibble(k_each = NA, HR = cox_p(colMeans(Z[up, ]))["HR"], p = cox_p(colMeans(Z[up, ]))["p"], note = "up genes only"),
  tibble(k_each = NA, HR = cox_p(-colMeans(Z[dn, ]))["HR"], p = cox_p(-colMeans(Z[dn, ]))["p"], note = "down genes only (sign-flipped)"))
cat("\nSize sensitivity:\n"); print(size_tab |> mutate(across(where(is.numeric), ~ signif(.x, 2))))

## (4) random-signature null: 1000 random 50/50 sets from the same training gene pool
set.seed(11)
pool <- names(pool_rho)
rand_p <- replicate(1000, { u <- sample(pool, 50); v <- sample(setdiff(pool, u), 50)
  cox_p(colMeans(Z[u, ]) - colMeans(Z[v, ]))["p"] })
emp_rand <- mean(rand_p <= obs["p"])
## (5) full-pipeline label-permutation null: shuffle the 15 training labels, retrain, test
perm_p <- replicate(1000, { rp <- apply(Etr, 1, function(v) suppressWarnings(cor(v, sample(tr$tumor_likeness), method = "spearman")))
  rp <- rp[!is.na(rp)]; cox_p(proxy_from(rp, 50))["p"] })
emp_perm <- mean(perm_p <= obs["p"])
cat(sprintf("\nObserved proxy p = %.2g. Random-signature null: %.1f%% of 1000 random 50/50 sets reach it (empirical p = %.3f).\nLabel-permutation null (retrain on shuffled methylome labels): empirical p = %.3f; median null p = %.2f\n",
            obs["p"], 100 * emp_rand, emp_rand, emp_perm, median(perm_p)))

## (6) individual KZFP genes (+ ZNF582, an oral-rinse methylation marker)
kz <- map_dfr(intersect(c("ZNF471", "ZNF382", "ZNF529", "ZNF582", "ZNF154", "ZNF671"), rownames(Z)), function(g) {
  r <- cox_p(Z[g, ]); tibble(gene = g, HR_per_SD_expression = r["HR"], p = r["p"]) })
cat("\nKZFP genes, expression vs OCFS:\n"); print(kz |> mutate(across(where(is.numeric), ~ signif(.x, 2))))

## (7) replication of the proxy in GSE227919 (25 lesions with follow-up)
ex2 <- read.table(gzfile(file.path(RNA_DIR, "GSE227919_2023_20_03_expr_mat_raw.txt.gz")), header = TRUE, check.names = FALSE)
an2 <- read.csv(gzfile(file.path(RNA_DIR, "GSE227919_2023_21_03_pml_annotations.csv.gz")),
                check.names = FALSE, na.strings = c("#N/A", "NA", "")); names(an2)[1] <- "Sample"
an2 <- an2[match(colnames(ex2), an2$Sample), ]
L2 <- log2(as.matrix(ex2) + 1); L2 <- L2[rowMeans(L2 > 1) >= 0.2, ]; Z2 <- t(scale(t(L2)))
u2 <- intersect(up, rownames(Z2)); v2 <- intersect(dn, rownames(Z2))
an2$proxy <- colMeans(Z2[u2, ]) - colMeans(Z2[v2, ])
an2$Class <- factor(an2$Class, levels = c("Control", "HkNR", "Dysplasia", "Cancer"))
pr2 <- an2 |> filter(Class != "Cancer", !is.na(Progression_status)) |> mutate(prog = Progression_status != "Stable")
auc <- function(x, y) { r <- rank(x); n1 <- sum(y); n0 <- sum(!y); (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
rep2 <- tibble(genes_used = paste0(length(u2), " up / ", length(v2), " down"),
               trend_rho = cor(an2$proxy, as.integer(an2$Class), method = "spearman"),
               trend_p = cor.test(an2$proxy, as.integer(an2$Class), method = "spearman", exact = FALSE)$p.value,
               prog_AUC = auc(pr2$proxy, pr2$prog),
               prog_p = wilcox.test(pr2$proxy[pr2$prog], pr2$proxy[!pr2$prog], exact = FALSE)$p.value,
               prog_p_adj_class = summary(glm(prog ~ proxy + Class, data = pr2, family = binomial))$coefficients["proxy", 4])
cat("\nProxy in GSE227919:\n"); print(rep2 |> mutate(across(where(is.numeric), ~ signif(.x, 2))))

saveRDS(list(clin = clin, feat = feat, cox_tab = cox_tab, c_h = c_h, fit_h = summary(fit_h)$coefficients,
             km = km, proxy_genes = list(up = up, dn = dn), mv_tab = mv_tab, assoc_proxy = assoc_proxy,
             w_dys = w_dys, size_tab = size_tab, rand_p = rand_p, perm_p = perm_p, obs = obs,
             emp_rand = emp_rand, emp_perm = emp_perm, kz = kz, rep2 = rep2, comp_check = comp_check,
             gse227919_proxy = an2 |> select(Sample, Class, Progression_status, proxy)),
        file.path(EXT_DIR, "gse26549_progression.rds"))
cat("done\n")
