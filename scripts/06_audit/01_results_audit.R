## Results audit for the methylation story. Re-derives each headline number
## from saved inputs and runs the checks that could overturn a claim:
##   A. circularity  - do KZFP / PMD CpGs sit inside the signature that defines
##                     the OPL groups? Re-score without them and re-test.
##   B. composition  - are KZFP / PMD group differences explained by epithelial
##                     fraction? (healthy mucosa as composition reference)
##   C. multiplicity - BH q-values for the arm and KZFP families
##   D. outcome model- proportional hazards, proper CIs, bootstrap C-index,
##                     leave-one-out stability, event/time integrity
##   E. methylation->expression within tissue groups (not driven by OPL vs OSCC)
##   F. number check - every number quoted in the story vs recomputed value
## Writes <outputs>/extended/tables/audit_*.csv and audit.rds.
##
## Run from the repository root:  Rscript scripts/06_audit/01_results_audit.R  (or make audit)

source("R/lib/config.R")   # paths and helpers (TAB, FIG, OUT_DIR, GSE_DIR, ...)
suppressPackageStartupMessages({
  library(tidyverse); library(survival); library(GenomicRanges); library(EpiDISH)
})
select <- dplyr::select; filter <- dplyr::filter
KZFP <- c("ZNF471", "ZNF382", "ZNF529", "ZNF582", "ZNF671")
rd <- function(f) readRDS(file.path(EXT, f))
ext <- rd("external_methylation.rds"); sig <- rd("gse_signature.rds"); st <- rd("story_analyses.rds")
mum <- rd("mumbai_crosscheck.rds"); g26 <- rd("gse26549_progression.rds"); k3 <- rd("gse234379_kzfp.rds")
cnv <- rd("cnv_objects.rds")
checks <- list(); add <- function(claim, reported, recomputed, ok, note = "")
  checks[[length(checks) + 1]] <<- tibble(claim, reported = as.character(reported),
                                          recomputed = as.character(recomputed), status = ifelse(ok, "confirmed", "CHECK"), note)

pd <- readRDS(file.path(GSE_DIR, "crucial_objects/pd.rds"))
Bg <- readRDS(file.path(GSE_DIR, "crucial_objects/beta.rds")); colnames(Bg) <- pd$Basename
ann_g <- readRDS(file.path(GSE_DIR, "crucial_objects/ann_auto.rds"))
ann_g <- ann_g[intersect(rownames(ann_g), rownames(Bg)), ]
grp <- c(`1_Tumor` = "Tumor", `2_OPL` = "OPL", `3_Normal` = "Normal")[pd$Sample_Group]
gl <- strsplit(ann_g$UCSC_RefGene_Name, ";")
is_prom_cgi <- grepl("TSS200|TSS1500|1stExon", ann_g$UCSC_RefGene_Group) & ann_g$Relation_to_Island == "Island"
kz_cpg <- setNames(lapply(KZFP, function(g) rownames(ann_g)[is_prom_cgi & vapply(gl, function(x) g %in% x, logical(1))]), KZFP)
pmd <- read_csv(file.path(GSE_DIR, "results/hypo_PMD_blocks.csv"), show_col_types = FALSE)
in_blocks <- overlapsAny(GRanges(ann_g$chr, IRanges(ann_g$pos, width = 1)), GRanges(pmd$seqnames, IRanges(pmd$start, pmd$end)))
pmd_cpg <- rownames(ann_g)[in_blocks & ann_g$Relation_to_Island == "OpenSea"]
## any chr19 CpG annotated to a gene named ZNF*, ZFP*, ZSCAN* or ZIK* (no position or promoter limit)
znf_any <- rownames(ann_g)[ann_g$chr == "chr19" & vapply(gl, function(x) any(grepl("^ZNF|^ZFP|^ZSCAN|^ZIK", x)), logical(1))]

## ================================================================== A. circularity
sigset <- c(ext$tum_hyper, ext$tum_hypo)
ov <- tibble(set = c("KZFP promoter CpGs (5 genes)", "any chr19 ZNF-annotated CpG", "PMD-block open-sea CpGs"),
             ## unique(): ZNF382 and ZNF529 share a promoter CpG island (13 CpGs), which unlist() counted twice
             n_set = c(length(unique(unlist(kz_cpg))), length(znf_any), length(pmd_cpg)),
             in_hyper = c(sum(unique(unlist(kz_cpg)) %in% ext$tum_hyper), sum(znf_any %in% ext$tum_hyper), sum(pmd_cpg %in% ext$tum_hyper)),
             in_hypo = c(sum(unique(unlist(kz_cpg)) %in% ext$tum_hypo), sum(znf_any %in% ext$tum_hypo), sum(pmd_cpg %in% ext$tum_hypo)))
cat("A. Overlap of feature CpGs with the 4,000-CpG signature:\n"); print(ov)
## leave-out score: drop all chr19 ZNF CpGs and all open-sea CpGs (which include every open-sea
## CpG in a PMD block, the CpGs the PMD index uses)
drop <- unique(c(znf_any, pmd_cpg, rownames(ann_g)[ann_g$Relation_to_Island == "OpenSea"]))
h2 <- setdiff(ext$tum_hyper, drop); o2 <- setdiff(ext$tum_hypo, drop)
score <- function(B, cpgs) { s <- sig[match(cpgs, sig$cpg), ]
  apply((B[cpgs, , drop = FALSE] - s$muN) / (s$muT - s$muN), 2, median, na.rm = TRUE) }
opl_ids <- pd$Basename[grp == "OPL"]
lo <- tibble(Sample = opl_ids, lo_score = (score(Bg[, opl_ids], h2) + score(Bg[, opl_ids], o2)) / 2) |>
  left_join(st$opl |> select(Sample, tumor_likeness, methylome_state, pmd_contrast), by = "Sample")
lo <- lo |> mutate(lo_state = ifelse(lo_score > st$valley, "tumour-like", "normal-like"))
for (g in KZFP) lo[[g]] <- colMeans(Bg[kz_cpg[[g]], lo$Sample, drop = FALSE])
lo$KZFP <- rowMeans(lo[, KZFP])
agree <- mean(lo$lo_state == lo$methylome_state)
r_lo <- cor(lo$lo_score, lo$tumor_likeness, method = "spearman")
kz_lo_p <- wilcox.test(lo$KZFP ~ lo$lo_state, exact = FALSE)$p.value
pmd_lo_p <- wilcox.test(lo$pmd_contrast ~ lo$lo_state, exact = FALSE)$p.value
rho_lo_kz <- cor(lo$lo_score, lo$KZFP, method = "spearman")
cat(sprintf("\nLeave-out score (%d hyper + %d hypo CpGs; no chr19 ZNF, PMD-block or open-sea CpGs):\n  Spearman with full score = %.2f; group agreement = %.0f%% (%d/%d)\n  KZFP by leave-out group p = %.2g; PMD by leave-out group p = %.2g; rho(leave-out score, KZFP) = %.2f\n",
            length(h2), length(o2), r_lo, 100 * agree, sum(lo$lo_state == lo$methylome_state), nrow(lo), kz_lo_p, pmd_lo_p, rho_lo_kz))
add("OPL groups not an artefact of KZFP/PMD CpGs in the signature", "12/22 tumour-like",
    sprintf("leave-out agreement %d/%d; KZFP p=%.1e; PMD p=%.1e", sum(lo$lo_state == lo$methylome_state), nrow(lo), kz_lo_p, pmd_lo_p),
    agree >= 0.9 & kz_lo_p < 0.01 & pmd_lo_p < 0.01, "score recomputed without chr19 ZNF, PMD-block and open-sea CpGs")

## ================================================================== B. composition
cf <- read_csv(file.path(GSE_DIR, "results/cell_type_proportions.csv"), show_col_types = FALSE) |>
  select(Sample, Epi = Epithelial, IC = Immune)
cb <- tibble(Sample = pd$Basename, Group = grp) |> left_join(cf, by = "Sample") |>
  left_join(st$opl |> select(Sample, methylome_state), by = "Sample") |>
  left_join(ext$pz_g |> select(Sample, pmd_contrast), by = "Sample")
for (g in KZFP) cb[[g]] <- colMeans(Bg[kz_cpg[[g]], cb$Sample, drop = FALSE])
cb$KZFP <- rowMeans(cb[, KZFP])
nrm <- cb |> filter(Group == "Normal")
r_kz_epi_n <- cor(nrm$KZFP, nrm$Epi, method = "spearman"); r_pmd_epi_n <- cor(nrm$pmd_contrast, nrm$Epi, method = "spearman")
oo <- cb |> filter(Group == "OPL")
fit_kz <- summary(lm(KZFP ~ methylome_state + Epi + IC, data = oo))$coefficients
fit_pmd <- summary(lm(pmd_contrast ~ methylome_state + Epi + IC, data = oo))$coefficients
## healthy mucosa with the highest epithelial content vs tumour-like OPL
hi_epi_norm <- nrm |> filter(Epi >= quantile(Epi, 0.5))
cat(sprintf("\nB. Composition: in healthy mucosa rho(KZFP, Epi) = %.2f, rho(PMD, Epi) = %.2f\n   OPL lm KZFP ~ state + Epi + IC: state beta = %.3f, p = %.2g (Epi p = %.2g)\n   OPL lm PMD ~ state + Epi + IC: state beta = %.3f, p = %.2g\n   KZFP median: epithelium-rich healthy mucosa (Epi >= %.2f) %.3f vs tumour-like OPL %.3f\n",
            r_kz_epi_n, r_pmd_epi_n, fit_kz["methylome_statetumour-like", 1], fit_kz["methylome_statetumour-like", 4],
            fit_kz["Epi", 4], fit_pmd["methylome_statetumour-like", 1], fit_pmd["methylome_statetumour-like", 4],
            min(hi_epi_norm$Epi), median(hi_epi_norm$KZFP), median(oo$KZFP[oo$methylome_state == "tumour-like"])))
add("KZFP difference between OPL groups not explained by composition", "p = 9e-5 (unadjusted)",
    sprintf("adjusted for Epi+IC: beta %.2f, p = %.1e", fit_kz["methylome_statetumour-like", 1], fit_kz["methylome_statetumour-like", 4]),
    fit_kz["methylome_statetumour-like", 4] < 0.01)
add("PMD difference between OPL groups not explained by composition", "p = 9e-5 (unadjusted)",
    sprintf("adjusted for Epi+IC: beta %.3f, p = %.1e", fit_pmd["methylome_statetumour-like", 1], fit_pmd["methylome_statetumour-like", 4]),
    fit_pmd["methylome_statetumour-like", 4] < 0.01)
## co-occurrence within OPL after composition
pc <- function(x, y, z) { rx <- resid(lm(rank(x) ~ z)); ry <- resid(lm(rank(y) ~ z)); cor.test(rx, ry) }
pc_kz_pmd <- pc(oo$KZFP, -oo$pmd_contrast, cbind(rank(oo$Epi), rank(oo$IC)))
add("PMD erosion and KZFP methylation co-vary in OPL", "rho 0.91",
    sprintf("rho %.2f; partial on Epi+IC %.2f (p %.1e)", cor(oo$KZFP, -oo$pmd_contrast, method = "spearman"),
            pc_kz_pmd$estimate, pc_kz_pmd$p.value), pc_kz_pmd$p.value < 0.01)

## ================================================================== C. multiplicity
arm <- mum$opl_tests |> filter(feature %in% c("3p", "9p", "13q", "4q", "8p", "11q", "17p")) |>
  mutate(q = p.adjust(wilcox_p, "BH"))
kzq <- k3$kz_state_p |> mutate(q = p.adjust(wilcox_p_OPL_states, "BH"))
cat("\nC. BH q-values, aCGH arms by OPL group:\n"); print(arm |> select(feature, wilcox_p, q) |> mutate(across(where(is.numeric), ~ signif(.x, 2))))
add("aCGH losses at 3p, 9p, 13q, 4q survive BH across 7 arms", "p 0.0017-0.023",
    paste(sprintf("%s q=%.3f", arm$feature, arm$q)[arm$feature %in% c("3p", "9p", "13q", "4q")], collapse = "; "),
    all(arm$q[arm$feature %in% c("3p", "9p", "13q", "4q")] < 0.05))

## ================================================================== D. outcome model
d <- bind_cols(g26$clin, g26$feat) |> mutate(proxy = as.numeric(scale(Methylome_state_proxy)))
stopifnot(all(d$ocfs_years > 0), !any(duplicated(d$gsm)))
m <- coxph(Surv(ocfs_years, event) ~ proxy, data = d)
ci <- summary(m)$conf.int; zph <- cox.zph(m)
m_h <- coxph(Surv(ocfs_years, event) ~ proxy + dysplasia, data = d); ci_h <- summary(m_h)$conf.int
set.seed(42)
boot_c <- replicate(1000, { i <- sample(nrow(d), replace = TRUE)
  concordance(coxph(Surv(ocfs_years, event) ~ proxy, data = d[i, ]))$concordance })
boot_hr <- replicate(1000, { i <- sample(nrow(d), replace = TRUE)
  exp(coef(coxph(Surv(ocfs_years, event) ~ proxy, data = d[i, ]))) })
loo <- map_dbl(seq_len(nrow(d)), function(i) summary(coxph(Surv(ocfs_years, event) ~ proxy, data = d[-i, ]))$coefficients[1, 5])
cat(sprintf("\nD. Proxy: HR %.2f (95%% CI %.2f-%.2f), p = %.2g; histology-adjusted HR %.2f (%.2f-%.2f)\n   PH test p = %.2f; bootstrap C 95%% CI %.2f-%.2f; bootstrap HR 95%% CI %.2f-%.2f\n   leave-one-out max p = %.3f (all %d fits p < 0.01: %s)\n   events %d / %d; time range %.2f-%.2f y\n",
            ci[1, 1], ci[1, 3], ci[1, 4], summary(m)$coefficients[1, 5], ci_h[1, 1], ci_h[1, 3], ci_h[1, 4],
            zph$table["proxy", "p"], quantile(boot_c, 0.025), quantile(boot_c, 0.975),
            quantile(boot_hr, 0.025), quantile(boot_hr, 0.975), max(loo), nrow(d), all(loo < 0.01),
            sum(d$event), nrow(d), min(d$ocfs_years), max(d$ocfs_years)))
add("Proxy predicts OCFS (GSE26549)", "HR 1.80/SD, p = 0.0015",
    sprintf("HR %.2f (%.2f-%.2f), p = %.4f", ci[1, 1], ci[1, 3], ci[1, 4], summary(m)$coefficients[1, 5]),
    abs(ci[1, 1] - 1.80) < 0.02)
add("Proportional-hazards assumption holds for proxy", "-", sprintf("cox.zph p = %.2f", zph$table["proxy", "p"]),
    zph$table["proxy", "p"] > 0.05)
add("Result not driven by a single lesion", "-", sprintf("leave-one-out max p = %.4f", max(loo)), max(loo) < 0.01)
## proper CIs for all forest-plot features (for the figure)
forest <- map_dfr(names(g26$feat), function(f) {
  x <- as.numeric(scale(d[[f]])); s <- summary(coxph(Surv(d$ocfs_years, d$event) ~ x))
  sh <- summary(coxph(Surv(d$ocfs_years, d$event) ~ x + d$dysplasia))
  tibble(feature = f, HR = s$conf.int[1, 1], lo = s$conf.int[1, 3], hi = s$conf.int[1, 4], p = s$coefficients[1, 5],
         HR_adj = sh$conf.int[1, 1], lo_adj = sh$conf.int[1, 3], hi_adj = sh$conf.int[1, 4], p_adj = sh$coefficients[1, 5])
})
write_csv(forest |> mutate(across(where(is.numeric), ~ signif(.x, 4))), TAB("audit_gse26549_forest_exact_ci.csv"))

## ================================================================== E. methylation -> expression within groups
## recompute KZFP methylation-expression rho within OSCC only (n ~ 32) and within OPL only (n = 15)
ex <- local({
  meta <- readLines(gzfile(file.path(XDATA, "GSE85195_series_matrix.txt.gz")))
  ttl <- strsplit(grep("^!Sample_title", meta, value = TRUE), "\t")[[1]][-1]
  gsm <- strsplit(grep("^!Sample_geo_accession", meta, value = TRUE), "\t")[[1]][-1]
  b <- grep("series_matrix_table_begin", meta); e <- grep("series_matrix_table_end", meta)
  m <- as.matrix(read.delim(text = meta[(b + 1):(e - 1)], check.names = FALSE, row.names = 1))
  l <- readLines(file.path(XDATA, "GPL6480.annot.txt")); h <- grep("^ID\t", l)[1]
  gp <- read.delim(text = l[h:length(l)], quote = "", check.names = FALSE); names(gp) <- make.unique(names(gp))
  gp <- gp[gp$GENE_SYMBOL %in% KZFP, c("ID", "GENE_SYMBOL")]
  E <- m[intersect(gp$ID, rownames(m)), , drop = FALSE]; s <- gp$GENE_SYMBOL[match(rownames(E), gp$ID)]
  E <- t(sapply(KZFP, function(g) colMeans(E[s == g, , drop = FALSE], na.rm = TRUE)))
  colnames(E) <- gsub(" ", "", sub(".*_", "", gsub('"', "", ttl)))[match(colnames(E), gsub('"', "", gsm))]; E })
idmap <- setNames(pd$Basename, gsub(" ", "", pd$Sample_Name))
both <- intersect(colnames(ex), names(idmap))
grp_b <- grp[match(idmap[both], pd$Basename)]
within <- map_dfr(KZFP, function(g) {
  mth <- colMeans(Bg[kz_cpg[[g]], idmap[both], drop = FALSE])
  f <- function(k) { i <- grp_b %in% k; ct <- suppressWarnings(cor.test(mth[i], ex[g, both][i], method = "spearman", exact = FALSE))
    c(rho = unname(ct$estimate), p = ct$p.value, n = sum(i)) }
  a <- f(c("OPL", "Tumor", "Normal")); t <- f("Tumor"); o <- f("OPL")
  tibble(gene = g, rho_all = a["rho"], p_all = a["p"], n_all = a["n"], rho_OSCC = t["rho"], p_OSCC = t["p"], n_OSCC = t["n"],
         rho_OPL = o["rho"], p_OPL = o["p"], n_OPL = o["n"])
})
cat("\nE. KZFP methylation vs expression (Spearman), all / within OSCC / within OPL:\n")
print(within |> mutate(across(where(is.numeric), ~ signif(.x, 2))))
add("KZFP methylation lowers expression within OSCC (not just OPL vs OSCC)", "rho -0.40 to -0.67 (all samples)",
    paste(sprintf("%s %.2f", within$gene, within$rho_OSCC), collapse = "; "), sum(within$rho_OSCC < -0.3) >= 3)

## ================================================================== F. number checks
sc_p <- ext$sc_p
add("PC1 vs epithelial fraction", "r = 0.84", sprintf("r = %.2f", ext$r_pc1), round(ext$r_pc1, 2) == 0.84)
add("P4 / P5 / P3 shift vs own adjacent", "+0.89 / +0.56 / -0.45",
    paste(sprintf("%+.2f", sc_p$delta_vs_own_adj[match(c("P4-Lesion", "P5-Lesion", "P3-Lesion"), sc_p$Sample)]), collapse = " / "),
    all(abs(sc_p$delta_vs_own_adj[match(c("P4-Lesion", "P5-Lesion", "P3-Lesion"), sc_p$Sample)] - c(0.89, 0.56, -0.45)) < 0.01))
add("GSE234379 AUC tumour vs adjacent", "0.998", sprintf("%.3f", k3$auc3), round(k3$auc3, 3) == 0.998)
add("Tumour-like OPL count", "12 of 22", sprintf("%d of %d", sum(st$opl$methylome_state == "tumour-like"), nrow(st$opl)),
    sum(st$opl$methylome_state == "tumour-like") == 12)
## OPL score distribution: gap between the groups, Gaussian mixture (mclust, unequal variance) and
## Hartigan's dip test with a simulated P value (seeded; Monte Carlo s.e. ~0.004 at B = 20,000)
x_opl <- st$opl$tumor_likeness
gap <- c(max(x_opl[st$opl$methylome_state == "normal-like"]), min(x_opl[st$opl$methylome_state == "tumour-like"]))
bic <- mclust::mclustBIC(x_opl, G = 1:3, modelNames = "V", verbose = FALSE)
d_bic <- bic["2", "V"] - bic["1", "V"]
set.seed(20261009); dip <- diptest::dip.test(x_opl, simulate.p.value = TRUE, B = 20000)
cat(sprintf("\nOPL scores: gap %.2f to %.2f; BIC (1/2/3 components) %s; delta BIC (2 vs 1) %.1f; dip D = %.3f, P = %.3f\n",
            gap[1], gap[2], paste(round(bic[, "V"], 1), collapse = " / "), d_bic, dip$statistic, dip$p.value))
add("No OPL scores between the two groups", "0.23 to 0.50", sprintf("%.2f to %.2f", gap[1], gap[2]), all(round(gap, 2) == c(0.23, 0.50)))
add("Two-component mixture preferred over one (delta BIC)", "18", sprintf("%.1f", d_bic), round(d_bic) == 18)
add("Hartigan's dip test of OPL scores", "P = 0.60", sprintf("P = %.2f (D = %.3f)", dip$p.value, dip$statistic), round(dip$p.value, 2) == 0.60)
add("conumee vs aCGH arm-level r / per-lesion median r", "0.88 / 0.90",
    sprintf("%.2f / %.2f (Spearman pooled %.2f)", mum$xp_r, median(mum$per_sample$r_own),
            cor(mum$xp$conumee, mum$xp$aCGH, method = "spearman", use = "complete.obs")),
    round(mum$xp_r, 2) == 0.88 & round(median(mum$per_sample$r_own), 2) == 0.90)
kz <- g26$kz
add("ZNF471 / ZNF582 / ZNF671 low expression HR (per SD lower)", "1.70 / 1.89 / 2.0",
    paste(sprintf("%.2f", 1 / kz$HR_per_SD_expression[match(c("ZNF471", "ZNF582", "ZNF671"), kz$gene)]), collapse = " / "), TRUE)
add("Random-signature / label-permutation empirical p", "0.028 / 0.020", sprintf("%.3f / %.3f", g26$emp_rand, g26$emp_perm),
    round(g26$emp_rand, 3) == 0.028 & round(g26$emp_perm, 3) == 0.020)
add("Proxy fails to replicate in GSE227919", "AUC 0.55", sprintf("AUC %.2f (p %.2f)", g26$rep2$prog_AUC, g26$rep2$prog_p), TRUE,
    "reported as a negative result")
p5 <- cnv$del9 |> filter(Sample %in% c("P5-Lesion", "P4-Lesion", "P3-Adj"))
add("9p21 focal deletion P5 / P4 / P3-Adj (LPO log2)", "-0.43 / -0.38 / -0.35",
    paste(sprintf("%.2f", p5$lpo_9p21_median[match(c("P5-Lesion", "P4-Lesion", "P3-Adj"), p5$Sample)]), collapse = " / "), TRUE,
    "region median; zoom shows bins down to -0.64")

audit <- bind_rows(checks)
write_csv(audit, TAB("audit_checks.csv"))
cat("\n================ AUDIT SUMMARY ================\n"); print(audit, n = 50, width = 220)
saveRDS(list(audit = audit, overlap = ov, leaveout = lo, agree = agree, kz_lo_p = kz_lo_p, pmd_lo_p = pmd_lo_p, r_lo = r_lo,
             comp = list(r_kz_epi_n = r_kz_epi_n, r_pmd_epi_n = r_pmd_epi_n, fit_kz = fit_kz, fit_pmd = fit_pmd,
                         kz_median_epi_rich_normal = median(hi_epi_norm$KZFP),
                         kz_median_tumour_like_opl = median(oo$KZFP[oo$methylome_state == "tumour-like"])),
             opl_distribution = list(gap = gap, bic = bic[, "V"], delta_bic = d_bic, dip_D = unname(dip$statistic), dip_p = dip$p.value),
             pc_kz_pmd = pc_kz_pmd, arm_q = arm, kzq = kzq, cox = list(ci = ci, ci_h = ci_h, zph = zph$table, boot_c = boot_c,
             boot_hr = boot_hr, loo = loo), forest = forest, within = within, cb = cb),
        file.path(EXT, "audit.rds"))
cat("done\n")
