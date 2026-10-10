## Every number quoted in the manuscript text, recomputed from the saved outputs and formatted
## as in the text, in reading order. After any change to the analysis, rerun (make audit) and
## compare the new table with the previous one; changed rows are the sentences to update:
##   Rscript tools/compare_tables.R <old>/extended/tables <new>/extended/tables
## It reads values produced by stages 01-06. A few summary statistics are computed here from saved
## tables (e.g. correlations of saved scores, PCA distances, group medians), as noted per line.
## Output: <outputs>/extended/tables/manuscript_numbers.csv (section, statement, value)
##
## Run from the repository root:  Rscript scripts/06_audit/03_manuscript_numbers.R  (or make audit)

source("R/lib/config.R")
suppressPackageStartupMessages(library(tidyverse))
rd <- function(f) readRDS(file.path(EXT_DIR, f))
tb <- function(f) read_csv(TAB(f), show_col_types = FALSE)
p <- readRDS(file.path(OUT_DIR, "processed.rds"))
ext <- rd("external_methylation.rds"); st <- rd("story_analyses.rds"); mum <- rd("mumbai_crosscheck.rds")
k3 <- rd("gse234379_kzfp.rds"); au <- rd("audit.rds"); rs <- rd("rnaseq_story.rds"); cnv <- rd("cnv_objects.rds")
g26 <- rd("gse26549_progression.rds"); rna <- rd("external_rnaseq.rds")

N <- list()
add <- function(section, statement, value) N[[length(N) + 1]] <<- tibble(section, statement, value = paste(value, collapse = "; "))
f1 <- function(x) sprintf("%.1f", x); f2 <- function(x) sprintf("%.2f", x); f3 <- function(x) sprintf("%.3f", x)
pc0 <- function(x) sprintf("%.0f%%", 100 * x); sg <- function(x, d = 2) formatC(signif(x, d), format = "g")
rng <- function(x, f = f2) paste(f(min(x)), "to", f(max(x)))
sx <- function(x) { x <- signif(x, 2); e <- floor(log10(x)); sprintf("%.1f x 10^%d", x / 10^e, e) }
big <- function(x) format(x, big.mark = ",")
sample_of <- function(s, x) x[match(s, names(x))]

## ------------------------------------------------------------------ Results 1: composition
sc_p <- ext$sc_p; dec <- ext$decon
add("R1", "pilot autosomal CpGs after filtering and collapsing (Methods)", big(nrow(p$betas)))
add("R1", "PC1 share of variance", pc0(p$pca$sdev[1]^2 / sum(p$pca$sdev^2)))
add("R1", "epithelial fraction: P3 (PVL) lesion / P4 adjacent", paste(pc0(dec$Epi_1st[dec$Sample == "P3-Lesion"]), pc0(dec$Epi_1st[dec$Sample == "P4-Adj"]), sep = " / "))
add("R1", "immune fraction of the P3 (PVL) lesion", pc0(dec$IC_1st[dec$Sample == "P3-Lesion"]))
add("R1", "fold range of epithelial fraction", sprintf("%.1f-fold", max(dec$Epi_1st) / min(dec$Epi_1st)))
add("R1", "PC1 vs epithelial fraction (Pearson r)", f2(ext$r_pc1))

## ------------------------------------------------------------------ Results 2: score
hpv <- tb("verification_gse234379_hpv.csv")
add("R2", "CpGs shared by pilot and GSE204943; r adjacent vs healthy", paste(big(ext$n_common), f2(ext$conc_r), sep = "; r = "))
add("R2", "GSE234379 AUC tumour vs adjacent / HPV-negative", paste(f3(k3$auc3), f2(hpv$auc[hpv$hpv == "Negative"]), sep = " / "))
add("R2", "score vs epithelial fraction: healthy mucosa / tumour-adjacent (Spearman)", paste(f2(ext$r_epi_norm), f2(k3$r_epi_n), sep = " / "))
od <- au$opl_distribution
add("R2", "no OPL between", paste(f2(od$gap[1]), "and", f2(od$gap[2])))
add("R2", "OPL groups normal-like / tumour-like", paste(sum(st$opl$methylome_state == "normal-like"), sum(st$opl$methylome_state == "tumour-like"), sep = " / "))
add("R2", "delta BIC (2 vs 1 component); dip test P", paste(sprintf("%.0f", od$delta_bic), f2(od$dip_p), sep = "; P = "))
add("R2", "group boundary (KDE minimum)", f2(st$valley))
sg_n <- ext$sc_g |> filter(Group == "Normal")
add("R2", "GSE204943 healthy-mucosa score range", rng(sg_n$tumor_likeness))
adj <- sc_p |> filter(Sample %in% c("P2-NAT", "P4-Adj", "P5-Adj"))
add("R2", "pilot adjacent scores P2 / P4 / P5", paste(f2(sample_of(c("P2-NAT", "P4-Adj", "P5-Adj"), setNames(sc_p$tumor_likeness, sc_p$Sample))), collapse = " / "))
sh <- setNames(sc_p$delta_vs_own_adj, sc_p$Sample)
add("R2", "shift vs own adjacent: P2 tumour / P4 AH / P5 dysplasia / P3 PVL", paste(sprintf("%+.2f", sample_of(c("P2-Tumor", "P4-Lesion", "P5-Lesion", "P3-Lesion"), sh)), collapse = " / "))
tl <- setNames(sc_p$tumor_likeness, sc_p$Sample)
add("R2", "score P5 lesion vs P2 tumour; P3 adjacent", paste(f2(tl["P5-Lesion"]), "vs", f2(tl["P2-Tumor"]), "; P3-Adj", f2(tl["P3-Adj"])))
short <- function(sn) { pt <- sub("Patient ", "P", sub("_.*", "", sn))   # as in 07_figures/01_figures_1_to_4.R
  paste0(pt, "-", ifelse(grepl("tumor", tolower(sn)), "Tumor", ifelse(grepl("lesion", tolower(sn)), "Lesion", ifelse(pt == "P2", "NAT", "Adj")))) }
pcx <- p$pca$x; lab <- short(rownames(pcx))
d2t <- sqrt((pcx[, 1] - pcx[lab == "P2-Tumor", 1])^2 + (pcx[, 2] - pcx[lab == "P2-Tumor", 2])^2); names(d2t) <- lab
les <- c("P3-Lesion", "P4-Lesion", "P5-Lesion")
add("R2", "Fig. 1f unadjusted PCA distance to P2 tumour: P3 / P4 / P5", paste(sprintf("%.0f", d2t[les]), collapse = " / "))
add("R2", "Fig. 1f rank by distance (1 = closest) / by adjusted shift: P3, P4, P5",
    paste(paste(rank(d2t[les]), collapse = ","), paste(rank(-sh[les]), collapse = ","), sep = " / "))

## ------------------------------------------------------------------ Results 3: PMD and KZFP
aud <- setNames(au$audit$recomputed, au$audit$claim)
fk <- au$comp$fit_kz; fp <- au$comp$fit_pmd
add("R3", "rho PMD vs KZFP in OPL", aud["PMD erosion and KZFP methylation co-vary in OPL"])
add("R3", "KZFP group difference adjusted for Epi+IC (beta, P)", paste(f2(fk["methylome_statetumour-like", 1]), sx(fk["methylome_statetumour-like", 4]), sep = ", P = "))
add("R3", "PMD group difference adjusted for Epi+IC (beta, P)", paste(f2(-fp["methylome_statetumour-like", 1]), sx(fp["methylome_statetumour-like", 4]), sep = ", P = "))
add("R3", "KZFP median: epithelium-rich healthy mucosa vs tumour-like OPL",
    paste(f2(au$comp$kz_median_epi_rich_normal), "vs", f2(au$comp$kz_median_tumour_like_opl)))
kza <- tb("verification_kzfp_cpg_accounting.csv")
add("R3", "KZFP promoter CpGs in signature / unique; ZNF382 island CpGs also ZNF529 (Methods)", paste0(kza$unique_in_signature, " of ", kza$unique_cpgs, "; ", kza$znf382_in_znf529))
add("R3", "leave-out score vs full score (Spearman); same group", paste0(f2(au$r_lo), "; ", sum(au$leaveout$lo_state == au$leaveout$methylome_state), " of ", nrow(au$leaveout)))
add("R3", "GSE234379 KZFP pairs up / P; HPV-negative pairs / P",
    with(hpv, paste0(n_pairs[hpv == "all"], " pairs, P = ", sx(exact_p[hpv == "all"]), "; ", n_pairs[hpv == "Negative"], " pairs, P = ", sx(exact_p[hpv == "Negative"]))))
kd <- tb("kzfp_pilot_deltas.csv")
add("R3", "P5 lesion KZFP delta beta range", paste(sprintf("%+.2f", min(kd$d_P5_progressor)), "to", sprintf("%+.2f", max(kd$d_P5_progressor))))
add("R3", "KZFP methylation-expression rho, all 48 samples", rng(au$within$rho_all))
w3 <- au$within |> filter(gene %in% c("ZNF471", "ZNF529", "ZNF671"))
add("R3", "rho within OSCC for ZNF471 / ZNF529 / ZNF671", paste(f2(w3$rho_OSCC), collapse = " / "))

## ------------------------------------------------------------------ Results 4: copy number
ps <- mum$per_sample
add("R4", "conumee vs aCGH median per-lesion r: OSCC / OPL; shuffled null", paste0(f2(median(ps$r_own[ps$Group == "Tumor"])), " / ", f2(median(ps$r_own[ps$Group == "OPL"])), "; null ", f2(median(mum$null_r))))
aq <- au$arm_q
add("R4", "aCGH by OPL group q: 3p, 9p, 13q, 4q (range); 8p, 11q, 17p (min)",
    paste0(rng(aq$q[aq$feature %in% c("3p", "9p", "13q", "4q")], f3), "; >= ", f2(min(aq$q[aq$feature %in% c("8p", "11q", "17p")]))))
d9 <- cnv$del9; s9 <- c("P5-Lesion", "P4-Lesion", "P3-Adj")
add("R4", "9p21 log2 ratio P5 / P4 / P3-Adj", paste(f2(d9$lpo_9p21_median[match(s9, d9$Sample)]), collapse = " / "))
add("R4", "9p21 probe-failure ratio P5 / P4 / P3-Adj", paste(f1(d9$fail_ratio[match(s9, d9$Sample)]), collapse = " / "))
cg <- tb("candidate_gene_promoter_methylation.csv") |> filter(gene == "CDKN2A")
add("R4", "CDKN2A promoter beta, max over P5 / P4 lesion and P3 adjacent", f3(max(unlist(cg[, s9]))))
ar <- cnv$arms |> filter(call != "")
for (s in c("P5-Lesion", "P4-Adj", "P4-Lesion", "P2-Tumor", "P3-Adj"))
  add("R4", paste("array arm calls", s), paste(sprintf("%s %s", ar$arm_label[ar$Sample == s], ar$call[ar$Sample == s]), collapse = ", "))

## ------------------------------------------------------------------ Results 5: outcome
ot <- tb("verification_mumbai_opl_outcome_tests.csv")
cl <- g26$clin; ev26 <- cl$outcome == "oral cancer development"
add("R5", "GSE26549 lesions; oral cancers; median follow-up of the cancer-free (years)", sprintf("%d; %d; %.1f", nrow(cl), sum(ev26), median(cl$ocfs_years[!ev26])))
add("R5", "transformed 4/12 vs 0/10 (Fisher P); transformed or recurred (P)", paste(f3(ot$fisher_p[ot$test == "transformed"]), f3(ot$fisher_p[ot$test == "transformed_or_recurred"]), sep = "; "))
fo <- tb("audit_gse26549_forest_exact_ci.csv") |> filter(feature == "Methylome_state_proxy")
add("R5", "proxy HR per s.d. (95% CI), P", sprintf("%.2f (%.2f-%.2f), P = %s", fo$HR, fo$lo, fo$hi, sx(fo$p)))
add("R5", "proxy HR adjusted for histology (95% CI)", sprintf("%.2f (%.2f-%.2f)", fo$HR_adj, fo$lo_adj, fo$hi_adj))
mv <- g26$mv_tab |> filter(grepl("proxy", term, ignore.case = TRUE))
add("R5", "proxy HR in multivariable model, P", sprintf("%.2f, P = %s", mv$HR[1], sx(mv$p[1])))
ct <- g26$cox_tab |> filter(feature == "Methylome_state_proxy")
add("R5", "C-index proxy + histology vs histology", paste(f2(ct$C_with_histology), "vs", f2(g26$c_h)))
add("R5", "random-set / label-permutation empirical P", paste(f3(g26$emp_rand), f3(g26$emp_perm), sep = " / "))
cc <- g26$comp_check$adj; rcc <- grep("proxy", rownames(cc), ignore.case = TRUE)[1]
add("R5", "proxy HR adjusted for composition scores, P", sprintf("%.1f, P = %s", exp(cc[rcc, 1]), sx(cc[rcc, ncol(cc)])))
kz <- g26$kz
add("R5", "low-expression HR ZNF471 / ZNF582 / ZNF529 (P)", paste(sprintf("%.2f (P %s)", 1 / kz$HR_per_SD_expression[match(c("ZNF471", "ZNF582", "ZNF529"), kz$gene)],
                                                                         sg(kz$p[match(c("ZNF471", "ZNF582", "ZNF529"), kz$gene)])), collapse = " / "))
l3 <- g26$cox_tab |> filter(feature == "loss_3p")
add("R5", "expression-inferred 3p loss HR, P", sprintf("%.2f, P = %s", l3$HR_per_SD, f3(l3$p)))
spat <- g26$cox_tab |> filter(feature %in% c("Invasive_epithelial", "Dysplastic_epithelial", "Cycling_epithelial", "LRRC15_CAF", "FOLR2_TAM", "Treg", "mregDC"))
add("R5", "spatial-profiling programmes, minimum Cox P", f2(min(spat$p)))
add("R5", "proxy in GSE227919: rho with grade (P); progression AUC", with(g26$rep2, sprintf("%.2f (P = %s); AUC %.2f", trend_rho, f3(trend_p), prog_AUC)))
sp7 <- rna$res |> filter(set %in% c("Invasive_epithelial", "Dysplastic_epithelial", "Cycling_epithelial", "LRRC15_CAF", "Treg", "mregDC", "FOLR2_TAM"))
inc <- sp7 |> filter(trend_fdr < 0.05)
add("R5", "spatial programmes rising with grade in GSE227919: n of 7, rho range, max FDR; exception",
    sprintf("%d of 7; %s; FDR <= %.3f; %s", nrow(inc), rng(inc$trend_rho), max(inc$trend_fdr), paste(setdiff(sp7$set, inc$set), collapse = ", ")))

## ------------------------------------------------------------------ Results 6-7: RNA vs methylation
ss <- rs$ss
add("R6", "read pairs per library (millions); pseudoaligned (%)", sprintf("%.1f-%.1f; %.1f-%.1f", min(ss$n_pairs) / 1e6, max(ss$n_pairs) / 1e6, min(ss$pct_pseudoaligned), max(ss$pct_pseudoaligned)))
cm <- tb("rnaseq_composition_vs_epidish.csv") |> filter(!is.na(Epi_meth))
add("R6", "RNA immune score vs EpiDISH immune (Spearman); non-epithelial vs epithelial",
    paste(f2(cor(cm$IC_meth, cm$Immune, method = "spearman")), f2(cor(cm$Epi_meth, -(cm$Immune + cm$Fibroblast), method = "spearman")), sep = "; "))
lpo <- tb("rnaseq_arm_dosage_vs_cn_lpo.csv") |> filter(!is.na(cn_lpo))
rr <- lpo |> group_by(sample) |> summarise(r = cor(rna_lpo, cn_lpo), .groups = "drop")
add("R6", "arm dosage vs array copy number r: P2 tumour / P4 adj / P4 lesion / P5 lesion",
    paste(f2(rr$r[match(c("P2-Tumor", "P4-Adj", "P4-Lesion", "P5-Lesion"), rr$sample)]), collapse = " / "))
n9 <- tb("rnaseq_9p21.csv"); g9 <- c("MTAP", "CDKN2A", "CDKN2B")
for (s9 in c("P4-Lesion", "P5-Lesion")) { x9 <- n9 |> filter(sample == s9)
  add("R6", paste("9p21 log2FC vs own adjacent, MTAP / CDKN2A / CDKN2B:", s9), paste(f2(x9$log2FC_vs_own_adj[match(g9, x9$gene)]), collapse = " / ")) }
add("R6", "all three 9p21 genes higher than adjacent: P2 / P6 / P7",
    paste(sapply(c("P2-Tumor", "P6-Lesion", "P7-Lesion"), function(s9) all(n9$log2FC_vs_own_adj[n9$sample == s9 & n9$gene %in% g9] > 0)), collapse = " / "))
ev <- read_csv(file.path(REPO_ROOT, "resources", "cnv_events_curated.csv"), show_col_types = FALSE)
dl <- function(s, a) mean(tb("rnaseq_arm_dosage_vs_cn_lpo.csv") |> filter(!is.na(cn_lpo), sample == s, arm_label %in% a) |> pull(rna_lpo))   # as in Fig. 5f
rna_call <- function(s, arms, dir) { if (is.na(arms)) return(NA_character_); v <- dl(s, strsplit(arms, ";")[[1]])
  if ((dir == "loss" && v < -rs$floor99) || (dir == "gain" && v > rs$floor99)) "Detected" else "Not detected" }
ev <- ev |> mutate(rna = coalesce(rna, pmap_chr(list(sample, arms, direction), rna_call)))
add("R6", "array events detected by RNA / assessable", paste(sum(ev$rna == "Detected"), "of", sum(ev$rna != "Not assessed")))
pp <- tb("rnaseq_proxy_programmes.csv"); m3 <- pp |> left_join(tb("rnaseq_composition_vs_epidish.csv") |> select(sample, Immune), by = "sample") |> filter(!is.na(tumor_likeness))
fh <- summary(lm(rank(proxy) ~ rank(tumor_likeness) + rank(Immune), data = m3))$coefficients
add("R7", "proxy vs methylation tumour-likeness (Spearman, n = 7); rank coefficient given immune, P",
    sprintf("%.2f; %.2f, P = %.3f", cor(m3$proxy, m3$tumor_likeness, method = "spearman"), fh[2, 1], fh[2, 4]))
add("R7", "genes excluded as composition-tracking", big(length(rs$nonepi)))
mc <- rs$mxs_c
for (s in c("P5-Lesion", "P4-Lesion", "P2-Tumor")) with(mc[mc$sample == s, ],
  add("R7", paste("gain >= 0.2 vs stable, median log2FC and P:", s), sprintf("%+.2f vs %+.2f; P = %s", med_expr_gain, med_expr_stable, sx(p_gain_vs_stable))))
ks <- rs$kz_shared
add("R7", "CpG-island promoters measured in P4 and P5; shared gain >= 0.2; of which KZFP", paste(big(rs$univ_n), nrow(rs$shared), sum(ks$kzfp), sep = "; "))
add("R7", "KZFP enrichment among shared gains: OR, P", sprintf("%.0f, P = %s", rs$f_gain$estimate, sx(rs$f_gain$p.value)))
add("R7", "KZFPs lower in both / KZFPs; other genes lower in both / other genes",
    sprintf("%d of %d; %d of %d", sum(ks$kzfp & ks$down_both), sum(ks$kzfp), sum(!ks$kzfp & ks$down_both), sum(!ks$kzfp)))
add("R7", "hypermethylated and less expressed in both: KZFP / all; OR, P",
    sprintf("%d of %d; %.0f, P = %s", sum(ks$kzfp & ks$down_both), sum(ks$down_both), rs$f_dn$estimate, sx(rs$f_dn$p.value)))
add("R7", "composition-tracking among these KZFPs", sum(ks$kzfp & ks$down_both & ks$composition_tracking))
kme <- tb("rnaseq_kzfp_meth_vs_expr.csv")
z471 <- kme |> filter(gene == "ZNF471", delta_beta > 0, !is.na(delta_log2_tpm))
add("R7", "ZNF471: lesions with promoter gain; log2FC range", paste(nrow(z471), "lesions;", rng(z471$delta_log2_tpm)))
ok <- !is.na(kme$delta_log2_tpm)
add("R7", "KZFP methylation vs expression change across gene-lesion pairs (Spearman, n)", sprintf("%.2f (n = %d)", cor(kme$delta_beta[ok], kme$delta_log2_tpm[ok], method = "spearman"), sum(ok)))

## ------------------------------------------------------------------ Results 8: generic lesion programme
de <- rs$de
add("R8", "genes tested; higher / lower at FDR < 0.05", sprintf("%s; %d / %d", big(nrow(de)), sum(de$FDR < 0.05 & de$logFC > 0), sum(de$FDR < 0.05 & de$logFC < 0)))
pdl <- setNames(rs$progA_delta$Proliferation, sub("\\..*$", "", rs$progA_delta$sample))
add("R8", "proliferation rise: P5; other lesions/tumours (range); library with highest score",
    sprintf("%+.2f; %+.2f to %+.2f; %s", pdl["P5-Lesion"], min(pdl[names(pdl) != "P5-Lesion"]), max(pdl[names(pdl) != "P5-Lesion"]),
            rownames(rs$score)[which.max(rs$score[, "Proliferation"])]))
add("R8", "lesion programme, GSE26549: HR per s.d. (95% CI), P; empirical P vs random sets",
    sprintf("%.2f (%.2f-%.2f), P = %.2f; %.2f", rs$obs["HR"], rs$obs["lo"], rs$obs["hi"], rs$obs["p"], rs$emp))
pr7 <- rs$a2$Class != "Cancer" & !is.na(rs$a2$Progression_status); yv7 <- rs$a2$Progression_status[pr7] != "Stable"
add("R8", "lesion programme, GSE227919: rho with grade (P); progression AUC (P)",
    sprintf("%.2f (P = %.4f); AUC %.2f (P = %.2f)", rs$gr$estimate, rs$gr$p.value, rs$auc27,
            wilcox.test(rs$prog27[pr7][yv7], rs$prog27[pr7][!yv7], exact = FALSE)$p.value))

## ------------------------------------------------------------------ Methods
pmdb <- read_csv(file.path(GSE_DIR, "results", "hypo_PMD_blocks.csv"), show_col_types = FALSE)
add("Methods", "GSE204943 hypomethylated PMD blocks: number; total Mb", sprintf("%d; %.1f", nrow(pmdb), sum(pmdb$width) / 1e6))
add("Methods", "RNA arm-dosage noise floor (99th percentile of absolute dosage in adjacent tissues)", f2(rs$floor99))
gb <- readRDS(file.path(GSE_DIR, "crucial_objects", "beta.rds"))
add("Methods", "GSE204943 CpGs and samples after detection-P filtering", paste(big(nrow(gb)), "CpGs;", ncol(gb), "samples")); rm(gb)
sx_sig <- tb("verification_sex_adjusted_signature.csv"); pm <- tb("verification_pmd_tumour_only_blocks.csv")
add("Methods", "sex-adjusted signature: CpGs retained; max pilot shift change", paste0(big(sx_sig$retained), " of 4,000; ", f3(sx_sig$max_pilot_shift_change)))
add("Methods", "tumour-only PMD blocks; Wilcoxon P (both block sets)", paste0(pm$blocks_tumour_only, " of ", pm$blocks_with_cpgs, "; P = ", sx(pm$opl_p_original), " / ", sx(pm$opl_p_tumour_only)))
add("Methods", "PMD index vs published erosion score (|rho|)", f2(abs(ext$r_pmd)))

## ------------------------------------------------------------------ Supplementary Table 2
o2 <- tb("verification_mumbai_opl_outcome.csv") |> arrange(desc(tumor_likeness))
add("Supp. Table 2", "OPL tumour-likeness (descending)", paste(sprintf("%s %.2f", o2$id, o2$tumor_likeness), collapse = ", "))

out <- bind_rows(N)
write_csv(out, TAB("manuscript_numbers.csv"))
options(width = 200); print(as.data.frame(out), right = FALSE)
