## Same-patient cross-checks in the Mumbai gingivobuccal cohort, where three
## platforms profiled overlapping lesions:
##   GSE204943  EPIC methylation (22 OPL, 74 OSCC, 22 normal)
##   GSE85514   array CGH vs pooled blood (24 OPL, 38 OSCC)
##   GSE85195   Agilent expression (15 OPL, 34 OSCC, 1 control)
##
##   A. Does copy number inferred from methylation arrays (conumee) agree with
##      array CGH in the same lesions?  (validates the method used in the pilot)
##   B. In OPL, is the tumour-like methylome state accompanied by copy-number
##      loss (3p/9p/13q...) measured on an independent platform?
##   C. Does promoter hypermethylation lower expression in matched lesions
##      (ZNF471, ZNF382, ... and the signature genes genome-wide)?
##   D. Do the manuscript's transcriptional programs differ by methylome state?
##
## Note: GPL22311 coordinates appear to be hg18; all comparisons are arm-level
## or >= 300 kb, where hg18/hg19 differences are immaterial.
##
## Run from the repository root:  Rscript scripts/03_methylation_cohorts/03_mumbai_crosscheck.R  (or make cohorts)

source("R/lib/config.R")   # paths and helpers (TAB, FIG, OUT_DIR, GSE_DIR, ...)
suppressPackageStartupMessages({
  library(tidyverse); library(sesameData); library(GenomicRanges)
})
LOH_ARMS <- c("3p", "9p", "13q", "4q", "8p", "11q", "17p")

read_series <- function(f) {
  meta <- readLines(gzfile(f))
  ttl <- strsplit(grep("^!Sample_title", meta, value = TRUE), "\t")[[1]][-1]
  gsm <- strsplit(grep("^!Sample_geo_accession", meta, value = TRUE), "\t")[[1]][-1]
  b <- grep("series_matrix_table_begin", meta); e <- grep("series_matrix_table_end", meta)
  m <- read.delim(text = meta[(b + 1):(e - 1)], check.names = FALSE, row.names = 1)
  oc <- gsub('"', "", ttl)
  list(m = as.matrix(m), map = tibble(gsm = gsub('"', "", gsm), title = oc,
       id = gsub(" ", "", sub(".*_", "", oc)),
       group = ifelse(grepl("OPL", oc), "OPL", ifelse(grepl("OSCC", oc), "Tumor", "Normal"))))
}
read_gpl <- function(f) {
  l <- readLines(f); h <- grep("^ID\t", l)[1]
  d <- read.delim(text = l[h:length(l)], quote = "", check.names = FALSE)
  names(d) <- make.unique(names(d))
  d[d$ID != "!platform_table_end", ]
}
arm_table <- function(build) {
  sesameData_getGenomeInfo(build)$cytoBand |>
    mutate(across(c(chrom, name), as.character)) |>
    filter(chrom %in% paste0("chr", 1:22)) |>
    mutate(arm = substr(name, 1, 1)) |>
    group_by(chrom, arm) |> summarise(start = min(chromStart), end = max(chromEnd), .groups = "drop") |>
    mutate(arm_label = paste0(sub("chr", "", chrom), arm))
}

## ------------------------------------------------------------------ GSE204943 side
st <- readRDS(file.path(EXT_DIR, "story_analyses.rds"))
pd <- readRDS(file.path(GSE_DIR, "crucial_objects/pd.rds"))
id_map <- tibble(Sample = pd$Basename, id = gsub(" ", "", pd$Sample_Name))
gse <- st$gse |> left_join(id_map, by = "Sample") |>
  left_join(st$opl |> select(Sample, methylome_state), by = "Sample")

## ================================================================== A/B: array CGH
cgh <- read_series(file.path(XDATA, "GSE85514_series_matrix.txt.gz"))
gpl <- read_gpl(file.path(XDATA, "GPL22311.annot.txt")) |>
  filter(grepl("^chr[0-9]+:", CHROMOSOMAL_LOCATION)) |>
  transmute(ID, chr = sub(":.*", "", CHROMOSOMAL_LOCATION),
            pos = as.numeric(sub(".*:([0-9]+)-.*", "\\1", CHROMOSOMAL_LOCATION)))
arms <- arm_table("hg19")
gpl$arm_label <- {
  gr <- GRanges(gpl$chr, IRanges(gpl$pos, width = 1))
  ar <- GRanges(arms$chrom, IRanges(arms$start, arms$end))
  h <- findOverlaps(gr, ar, select = "first"); arms$arm_label[h]
}
gpl <- gpl |> filter(!is.na(arm_label))
M <- cgh$m[intersect(gpl$ID, rownames(cgh$m)), ]
g2 <- gpl[match(rownames(M), gpl$ID), ]
colnames(M) <- cgh$map$id[match(colnames(M), cgh$map$gsm)]
cgh_grp <- setNames(cgh$map$group, cgh$map$id)
M <- sweep(M, 2, apply(M, 2, median, na.rm = TRUE))
cgh_arm <- sapply(split(seq_len(nrow(M)), g2$arm_label), function(i) apply(M[i, , drop = FALSE], 2, median, na.rm = TRUE))
## focal 9p21 (MTAP-CDKN2A-CDKN2B, ~21.8-22.1 Mb in hg18/hg19)
f9 <- g2$chr == "chr9" & g2$pos > 21.75e6 & g2$pos < 22.15e6
cgh_9p21 <- apply(M[f9, , drop = FALSE], 2, median, na.rm = TRUE)
cgh_df <- as_tibble(cgh_arm, rownames = "id") |> mutate(group = cgh_grp[id], focal_9p21 = cgh_9p21[id])
cat(sprintf("Array CGH: %d probes on autosomes, %d samples (%s)\n", nrow(M), ncol(M),
            paste(names(table(cgh_df$group)), table(cgh_df$group), collapse = ", ")))

## positive control: OSCC arm medians (expected 3q, 8q, 20 gain; 3p, 9p, 4q, 18q loss)
pc <- cgh_df |> filter(group == "Tumor") |> select(-id, -group) |>
  summarise(across(everything(), median)) |> pivot_longer(everything(), names_to = "arm", values_to = "median_OSCC") |>
  arrange(median_OSCC)
cat("OSCC arm medians (array CGH), most lost and most gained:\n")
print(bind_rows(head(pc, 6), tail(pc, 6)) |> mutate(median_OSCC = round(median_OSCC, 3)))

## A. cross-platform agreement, samples profiled by both
both <- gse |> filter(id %in% cgh_df$id)
arm_cols <- intersect(colnames(cgh_arm), names(st$gse))
long_m <- both |> select(id, Group, all_of(arm_cols)) |> pivot_longer(all_of(arm_cols), names_to = "arm", values_to = "conumee")
long_c <- cgh_df |> select(id, all_of(arm_cols)) |> pivot_longer(all_of(arm_cols), names_to = "arm", values_to = "aCGH")
xp <- inner_join(long_m, long_c, by = c("id", "arm"))
xp_r <- cor(xp$conumee, xp$aCGH, use = "complete.obs")
## sample-wise: does each sample's arm profile match its own aCGH profile better than other samples'?
per_sample <- xp |> group_by(id, Group) |> summarise(r_own = cor(conumee, aCGH), .groups = "drop")
null_r <- map_dbl(1:500, function(k) { set.seed(k)
  sh <- xp |> group_by(arm) |> mutate(aCGH = sample(aCGH)) |> ungroup()
  sh |> group_by(id) |> summarise(r = cor(conumee, aCGH), .groups = "drop") |> pull(r) |> median() })
## restrict to clearly altered arms for an amplitude-independent check
big <- xp |> filter(abs(aCGH) > 0.1)
conc_dir <- mean(sign(big$aCGH) == sign(big$conumee))
cat(sprintf("\nA. methylation-array CN (conumee) vs array CGH, same lesions: n = %d samples (%s)\n",
            n_distinct(xp$id), paste(names(table(both$Group)), table(both$Group), collapse = ", ")))
cat(sprintf("   pooled arm-level r = %.2f; median per-sample r = %.2f (label-shuffled null median %.2f);\n   direction agreement on arms with |aCGH| > 0.1: %.0f%% (n = %d)\n",
            xp_r, median(per_sample$r_own), median(null_r), 100 * conc_dir, nrow(big)))
write_csv(per_sample |> mutate(r_own = round(r_own, 3)), TAB("mumbai_conumee_vs_acgh_per_sample.csv"))

## B. OPL: copy number (aCGH) by methylome state
opl_c <- cgh_df |> filter(group == "OPL") |>
  left_join(gse |> filter(Group == "OPL") |> select(id, methylome_state, tumor_likeness, Epi), by = "id") |>
  mutate(instability = apply(across(all_of(arm_cols)), 1, function(v) mean(abs(v) > 0.1)),
         core_loss = -(`3p` + `9p`) / 2)
cat(sprintf("\nB. OPL with array CGH: %d (%d with methylome state)\n", nrow(opl_c), sum(!is.na(opl_c$methylome_state))))
opl_tests <- map_dfr(c(LOH_ARMS, "focal_9p21", "core_loss", "instability"), function(v) {
  d <- opl_c |> filter(!is.na(methylome_state))
  tibble(feature = v,
         normal_like = median(d[[v]][d$methylome_state == "normal-like"]),
         tumour_like = median(d[[v]][d$methylome_state == "tumour-like"]),
         wilcox_p = wilcox.test(d[[v]] ~ d$methylome_state, exact = FALSE)$p.value,
         rho_tumour_likeness = cor(d[[v]], d$tumor_likeness, method = "spearman"),
         n_lost_lt_minus0.1 = if (v %in% c(LOH_ARMS, "focal_9p21")) sum(d[[v]] < -0.1) else NA_integer_)
})
write_csv(opl_tests |> mutate(across(where(is.numeric), ~ signif(.x, 3))), TAB("mumbai_opl_acgh_by_methylome_state.csv"))
print(opl_tests |> mutate(across(where(is.numeric), ~ signif(.x, 2))))
write_csv(opl_c |> mutate(across(where(is.numeric), ~ round(.x, 3))), TAB("mumbai_opl_acgh_arms.csv"))

## ================================================================== C/D: expression
ex <- read_series(file.path(XDATA, "GSE85195_series_matrix.txt.gz"))
g6480 <- read_gpl(file.path(XDATA, "GPL6480.annot.txt")) |> select(ID, GENE_SYMBOL) |> filter(GENE_SYMBOL != "")
E <- ex$m[intersect(rownames(ex$m), g6480$ID), ]
sym <- g6480$GENE_SYMBOL[match(rownames(E), g6480$ID)]
## one probe per gene: highest mean
keep <- tapply(seq_along(sym), sym, function(i) i[which.max(rowMeans(E[i, , drop = FALSE], na.rm = TRUE))])
E <- E[unlist(keep), ]; rownames(E) <- names(keep)
colnames(E) <- ex$map$id[match(colnames(E), ex$map$gsm)]
cat(sprintf("\nGSE85195 expression: %d genes x %d samples; value range %.1f to %.1f\n", nrow(E), ncol(E),
            min(E, na.rm = TRUE), max(E, na.rm = TRUE)))
common_ids <- intersect(colnames(E), gse$id)
cat(sprintf("Samples with both methylation (GSE204943) and expression (GSE85195): %d (%s)\n", length(common_ids),
            paste(names(table(gse$Group[match(common_ids, gse$id)])), table(gse$Group[match(common_ids, gse$id)]), collapse = ", ")))

Bg <- readRDS(file.path(GSE_DIR, "crucial_objects/beta.rds")); colnames(Bg) <- pd$Basename
ann_g <- readRDS(file.path(GSE_DIR, "crucial_objects/ann_auto.rds"))
ann_g <- ann_g[intersect(rownames(ann_g), rownames(Bg)), ]
prom <- grepl("TSS200|TSS1500|1stExon", ann_g$UCSC_RefGene_Group) & ann_g$Relation_to_Island == "Island"
gl <- strsplit(ann_g$UCSC_RefGene_Name[prom], ";")
prom_map <- tibble(cpg = rep(rownames(ann_g)[prom], lengths(gl)), gene = unlist(gl)) |> distinct()
basename_of <- setNames(gse$Sample, gse$id)
Bm <- Bg[, basename_of[common_ids]]; colnames(Bm) <- common_ids
Em <- E[, common_ids]
gene_meth <- prom_map |> filter(gene %in% rownames(Em)) |>
  group_by(gene) |> summarise(cpgs = list(cpg), .groups = "drop")
gmeth <- t(sapply(gene_meth$cpgs, function(c) colMeans(Bm[c, , drop = FALSE])))
rownames(gmeth) <- gene_meth$gene
## methylation-expression correlation across matched samples, per gene
rho_me <- sapply(rownames(gmeth), function(g) suppressWarnings(cor(gmeth[g, ], Em[g, ], method = "spearman")))
sig <- readRDS(file.path(EXT_DIR, "gse_signature.rds"))
ext <- readRDS(file.path(EXT_DIR, "external_methylation.rds"))
hyper_genes <- prom_map |> filter(cpg %in% ext$tum_hyper) |> pull(gene) |> unique() |> intersect(names(rho_me))
me_tab <- tibble(gene = names(rho_me), rho = rho_me, sig_hyper = gene %in% hyper_genes,
                 meth_sd = apply(gmeth, 1, sd), expr_mean = rowMeans(Em[names(rho_me), ])) |>
  filter(meth_sd > 0.05)
cat(sprintf("\nC. Promoter-CGI methylation vs expression across %d matched samples (genes with methylation SD > 0.05: %d)\n",
            length(common_ids), nrow(me_tab)))
print(me_tab |> group_by(sig_hyper) |> summarise(n = n(), median_rho = round(median(rho, na.rm = TRUE), 3),
                                                  frac_rho_lt_minus0.3 = round(mean(rho < -0.3, na.rm = TRUE), 3)))
w_me <- wilcox.test(rho ~ sig_hyper, data = me_tab)
cat(sprintf("   signature-hyper genes vs other variable genes: Wilcoxon p = %.2g\n", w_me$p.value))
focus <- c("ZNF471", "ZNF382", "ZNF529", "NELL1", "HOXA5", "TFAP2A", "CDKN2A", "MGMT", "FAT1", "HMX3", "TRIM71")
focus_tab <- me_tab |> filter(gene %in% focus) |> mutate(across(where(is.numeric), ~ round(.x, 3)))
print(focus_tab)
write_csv(me_tab |> mutate(across(where(is.numeric), ~ round(.x, 3))), TAB("mumbai_methylation_expression_rho.csv"))

## D. programs by methylome state among OPLs with expression
progs <- list(Invasive_epithelial = c("LAMC2", "COL17A1", "PTHLH", "ITGA6", "ITGB4", "HMGA2"),
              Dysplastic_epithelial = c("MAGEA11", "ITGB6", "LAMC2", "TNFRSF12A", "HBEGF", "ULBP2"),
              Cycling_epithelial = c("HMMR", "ASPM", "BUB1", "NEIL3", "BUB1B", "PBK"),
              LRRC15_CAF = c("LRRC15", "FAP", "COL11A1"), FOLR2_TAM = c("FOLR2", "CD163", "MRC1"),
              Treg = c("FOXP3", "IL2RA", "CTLA4", "TNFRSF18"), mregDC = c("LAMP3", "CCR7", "FSCN1", "CD274"))
Z <- t(scale(t(E)))
psc <- sapply(progs, function(g) colMeans(Z[intersect(g, rownames(Z)), , drop = FALSE], na.rm = TRUE))
opl_e <- tibble(id = rownames(psc)) |> bind_cols(as_tibble(psc)) |>
  inner_join(gse |> filter(Group == "OPL") |> select(id, methylome_state, tumor_likeness), by = "id")
cat(sprintf("\nD. OPL with expression and methylome state: %d (%s)\n", nrow(opl_e),
            paste(names(table(opl_e$methylome_state)), table(opl_e$methylome_state), collapse = ", ")))
prog_state <- map_dfr(names(progs), function(p) tibble(program = p,
  normal_like = median(opl_e[[p]][opl_e$methylome_state == "normal-like"]),
  tumour_like = median(opl_e[[p]][opl_e$methylome_state == "tumour-like"]),
  wilcox_p = wilcox.test(opl_e[[p]] ~ opl_e$methylome_state, exact = FALSE)$p.value,
  rho_tumour_likeness = cor(opl_e[[p]], opl_e$tumor_likeness, method = "spearman")))
print(prog_state |> mutate(across(where(is.numeric), ~ signif(.x, 2))))
write_csv(prog_state |> mutate(across(where(is.numeric), ~ signif(.x, 3))), TAB("mumbai_opl_programs_by_methylome_state.csv"))

saveRDS(list(xp = xp, xp_r = xp_r, per_sample = per_sample, null_r = null_r, conc_dir = conc_dir,
             pc = pc, opl_c = opl_c, opl_tests = opl_tests, me_tab = me_tab, w_me_p = w_me$p.value,
             focus_tab = focus_tab, prog_state = prog_state, opl_e = opl_e, n_common = length(common_ids)),
        file.path(EXT_DIR, "mumbai_crosscheck.rds"))
cat("done\n")
