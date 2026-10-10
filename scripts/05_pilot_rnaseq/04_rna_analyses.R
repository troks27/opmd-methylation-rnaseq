## Pilot bulk RNA-seq: analyses for the RNA section of the manuscript.
## Inputs: gene matrices (02_gene_matrices.R), pilot and cohort methylation outputs (stages 02-03), GSE26549 and GSE227919 (stage 04).
## Descriptive throughout: 11 libraries, 5 lesion/adjacent pairs (P2, P4-P7), P3 lesion only,
## 4 patients with methylation; P6 and P7 have RNA only.
##   A. programme / cell-type scores across the 11 libraries
##   B. paired lesion-versus-adjacent differential expression (edgeR QL) + hallmark sets (camera)
##   C. genome-wide promoter methylation change vs expression change within patients (P2, P4, P5)
##   D. RNA-inferred arm-level copy number for all libraries, incl. P6 and P7
##   E. external test of the pilot lesion programme: GSE26549 (oral cancer-free survival) and
##      GSE227919 (grade, progression), against 1,000 random gene sets
##   F. robustness: coupling after removing composition-tracking genes; genes hypermethylated and
##      silenced in both premalignant lesions; skeletal-muscle content of biopsies
##   G. enrichment of chromosome 19 KZFP genes among promoters hypermethylated (and silenced) in both
##      premalignant lesions (P4 AH, P5 MSD)
## Run from the repository root:  Rscript scripts/05_pilot_rnaseq/04_rna_analyses.R  (or make rnaseq; ~5 min)
source("R/lib/config.R")   # paths and helpers (TAB, FIG, OUT_DIR, GSE_DIR, ...)
suppressPackageStartupMessages({
  library(tidyverse); library(edgeR); library(limma); library(survival); library(msigdbr)
  library(sesameData); library(GenomicRanges); library(hugene10sttranscriptcluster.db)
})
select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename; first <- dplyr::first
set.seed(20261009)

## ------------------------------------------------------------------ data
ss <- read_tsv(file.path(MAT, "samples.tsv"), show_col_types = FALSE)
ss <- ss |> mutate(sample = paste0(patient, "-", c(Tumor = "Tumor", Lesion = "Lesion", Adjacent = "Adj")[tissue]),
                   sample = ifelse(sample == "P2-Adj", "P2-NAT", sample),
                   lesion = as.integer(tissue != "Adjacent"),
                   diagnosis = case_when(tissue == "Adjacent" ~ "Adjacent tissue",
                                         patient == "P2" ~ "OSCC (T2)", patient == "P3" ~ "PVL",
                                         patient == "P4" ~ "Atypical hyperplasia",
                                         patient == "P5" ~ "Moderate-severe dysplasia (progressed)",
                                         patient == "P6" ~ "OSCC (per multimodal draft)",
                                         patient == "P7" ~ "Not available"))
ann <- read_tsv(file.path(MAT, "gene_annotation.tsv"), show_col_types = FALSE)
rd_mat <- function(f) read_tsv(file.path(MAT, f), show_col_types = FALSE) |> column_to_rownames("gene_id") |> as.matrix()
tpm <- rd_mat("gene_tpm.tsv.gz"); cnt <- rd_mat("gene_counts_lengthScaledTPM.tsv.gz")
## symbol level: keep the most expressed gene_id per symbol
sym <- ann$gene_name[match(rownames(tpm), ann$gene_id)]
keep_id <- unlist(tapply(seq_along(sym), sym, function(i) i[which.max(rowMeans(tpm[i, , drop = FALSE]))]))
gtype <- setNames(ann$gene_type[match(rownames(tpm)[keep_id], ann$gene_id)], sym[keep_id])
tpm <- tpm[keep_id, ]; cnt <- cnt[keep_id, ]; rownames(tpm) <- rownames(cnt) <- sym[keep_id]
colnames(tpm) <- colnames(cnt) <- ss$sample[match(colnames(tpm), ss$library)]
lt <- log2(tpm + 1); zs <- t(scale(t(lt[rowSums(tpm >= 1) >= 2, ])))
adj_of <- setNames(ifelse(ss$patient == "P2", "P2-NAT", paste0(ss$patient, "-Adj")), ss$sample)
pairs <- ss |> filter(lesion == 1, adj_of[sample] %in% ss$sample) |> pull(sample)   # P2, P4, P5, P6, P7
cat(sprintf("Genes (symbols): %d; libraries %d; lesion/adjacent pairs: %s\n", nrow(tpm), ncol(tpm), paste(pairs, collapse = ", ")))

## ------------------------------------------------------------------ A. programme scores
rna_ext <- readRDS(file.path(EXT, "external_rnaseq.rds"))
spatial <- rna_ext$sets[c("Invasive_epithelial", "Dysplastic_epithelial", "Cycling_epithelial", "LRRC15_CAF",
                          "FOLR2_TAM", "Treg", "mregDC")]
cells <- list(T_cell = c("CD3D", "CD3E", "CD2", "CD5", "CD6", "CD247"),
              Cytotoxic = c("CD8A", "CD8B", "GZMB", "GZMA", "PRF1", "NKG7"),
              B_cell = c("MS4A1", "CD79A", "CD19", "CD22", "BANK1", "FCRL5"),
              Plasma_cell = c("MZB1", "JCHAIN", "IGKC", "IGHG1", "XBP1", "TNFRSF17"),
              Myeloid = c("CD68", "CD14", "LYZ", "CD163", "C1QA", "C1QB"),
              IFNG_response = c("CXCL9", "CXCL10", "CXCL11", "IDO1", "STAT1", "GBP1"),
              Checkpoint = c("PDCD1", "CD274", "LAG3", "HAVCR2", "TIGIT", "CTLA4"),
              Proliferation = c("MKI67", "TOP2A", "CCNB1", "CDK1", "BIRC5", "UBE2C", "CENPF", "TPX2", "AURKA", "MCM2"),
              Differentiation = c("KRT4", "KRT13", "CRNN", "SPRR3", "TGM3", "CRCT1", "MAL"),
              Basal = c("KRT5", "KRT14", "TP63", "KRT15"),
              pEMT = c("LAMC2", "LAMB3", "LAMA3", "TGFBI", "PTHLH", "ITGA5", "MMP10", "SERPINE1", "PDPN", "TNC"),
              CAF = c("FAP", "COL1A1", "COL3A1", "POSTN", "COL11A1", "LRRC15"),
              Endothelial = c("PECAM1", "VWF", "CDH5", "CLDN5"))
sets <- c(cells, spatial)
score <- sapply(sets, function(g) colMeans(zs[intersect(g, rownames(zs)), , drop = FALSE]))
progA <- as_tibble(score, rownames = "sample") |> left_join(ss |> select(sample, patient, tissue, diagnosis), by = "sample")
delta <- function(v) sapply(pairs, function(s) unname(v[s] - v[adj_of[s]]))
progA_delta <- as_tibble(apply(score, 2, delta), rownames = "sample")
cat("\nA. Programme scores (z), sorted by proliferation\n")
print(progA |> select(sample, Proliferation, pEMT, Invasive_epithelial, T_cell, B_cell, Plasma_cell, Cytotoxic, IFNG_response, Differentiation) |>
        arrange(desc(Proliferation)) |> mutate(across(where(is.double), ~ round(.x, 2))), width = 200)
write_csv(progA, TAB("rnaseq_programme_scores.csv")); write_csv(progA_delta, TAB("rnaseq_programme_deltas.csv"))

## ------------------------------------------------------------------ B. paired DE
y <- DGEList(round(cnt[, c(pairs, adj_of[pairs])]))
pt <- factor(sub("-.*", "", colnames(y))); les <- as.integer(!grepl("Adj|NAT", colnames(y)))
design <- model.matrix(~ pt + les)
y <- y[filterByExpr(y, design = design), , keep.lib.sizes = FALSE]; y <- normLibSizes(y)
y <- estimateDisp(y, design); fit <- glmQLFit(y, design, robust = TRUE)
qlf <- glmQLFTest(fit, coef = "les")
de <- as_tibble(topTags(qlf, n = Inf)$table, rownames = "gene") |> mutate(gene_type = gtype[gene])
cat(sprintf("\nB. Paired lesion vs adjacent (n = %d pairs): %d genes tested; FDR < 0.05: %d up, %d down; FDR < 0.10: %d\n",
            length(pairs), nrow(de), sum(de$FDR < 0.05 & de$logFC > 0), sum(de$FDR < 0.05 & de$logFC < 0), sum(de$FDR < 0.1)))
print(de |> filter(gene_type == "protein_coding") |> slice_head(n = 25) |> select(gene, logFC, logCPM, PValue, FDR) |>
        mutate(across(where(is.double), ~ signif(.x, 3))), n = 25)
write_csv(de, TAB("rnaseq_paired_de.csv"))
hm <- tryCatch(msigdbr(species = "Homo sapiens", collection = "H"), error = function(e) msigdbr(species = "Homo sapiens", category = "H"))
hidx <- ids2indices(split(hm$gene_symbol, hm$gs_name), rownames(y))
cam <- camera(y, hidx, design, contrast = ncol(design)) |> as_tibble(rownames = "set") |> arrange(PValue)
cat("\nHallmark gene sets (camera, paired design), FDR < 0.1:\n")
print(cam |> filter(FDR < 0.1) |> mutate(across(where(is.double), ~ signif(.x, 2))), n = 30)
write_csv(cam, TAB("rnaseq_paired_hallmarks.csv"))

## ------------------------------------------------------------------ C. methylation vs expression
p <- readRDS(file.path(OUT_DIR, "processed.rds"))
short_m <- function(sn) { pt <- sub("Patient ", "P", sub("_.*", "", sn))
  paste0(pt, "-", ifelse(grepl("tumor", tolower(sn)), "Tumor", ifelse(grepl("lesion", tolower(sn)), "Lesion", ifelse(pt == "P2", "NAT", "Adj")))) }
B <- p$betas; colnames(B) <- short_m(colnames(B)); rm(p); invisible(gc())
ann_g <- readRDS(file.path(GSE_DIR, "crucial_objects/ann_auto.rds"))
ann_p <- ann_g[intersect(rownames(ann_g), rownames(B)), ]
is_prom <- grepl("TSS200|TSS1500|1stExon", ann_p$UCSC_RefGene_Group) & ann_p$Relation_to_Island == "Island"
gl <- strsplit(ann_p$UCSC_RefGene_Name[is_prom], ";")
pm <- tibble(cpg = rep(rownames(ann_p)[is_prom], lengths(gl)), gene = unlist(gl)) |> distinct() |> filter(gene %in% rownames(lt))
pmB <- rowsum(B[pm$cpg, ], pm$gene) / as.vector(table(pm$gene)[sort(unique(pm$gene))])
ext <- readRDS(file.path(EXT, "external_methylation.rds"))
sig_hyper_genes <- pm |> filter(cpg %in% ext$tum_hyper) |> pull(gene) |> unique()
mx <- map_dfr(c("P2-Tumor", "P4-Lesion", "P5-Lesion"), function(s) { a <- adj_of[s]
  tibble(sample = s, gene = rownames(pmB), d_beta = pmB[, s] - pmB[, a], beta_adj = pmB[, a],
         d_expr = lt[rownames(pmB), s] - lt[rownames(pmB), a], tpm_adj = tpm[rownames(pmB), a], tpm_les = tpm[rownames(pmB), s]) }) |>
  filter(pmax(tpm_adj, tpm_les) >= 1) |>
  mutate(cls = cut(d_beta, c(-Inf, -0.2, -0.05, 0.05, 0.2, Inf), labels = c("Loss >= 0.2", "Loss 0.05-0.2", "Stable", "Gain 0.05-0.2", "Gain >= 0.2")),
         sig_hyper = gene %in% sig_hyper_genes)
mxs <- mx |> group_by(sample) |> summarise(
  n_genes = n(), n_gain = sum(d_beta >= 0.2), n_loss = sum(d_beta <= -0.2),
  rho = cor(d_beta, d_expr, method = "spearman"),
  med_expr_gain = median(d_expr[d_beta >= 0.2]), med_expr_stable = median(d_expr[abs(d_beta) < 0.05]),
  med_expr_loss = median(d_expr[d_beta <= -0.2]),
  p_gain_vs_stable = wilcox.test(d_expr[d_beta >= 0.2], d_expr[abs(d_beta) < 0.05])$p.value,
  p_loss_vs_stable = wilcox.test(d_expr[d_beta <= -0.2], d_expr[abs(d_beta) < 0.05])$p.value,
  ## silencing needs baseline expression: genes expressed in adjacent tissue (TPM >= 5)
  med_expr_gain_expressed = median(d_expr[d_beta >= 0.2 & tpm_adj >= 5]), n_gain_expressed = sum(d_beta >= 0.2 & tpm_adj >= 5),
  med_expr_sigHyper_gain = median(d_expr[sig_hyper & d_beta >= 0.1]), n_sigHyper_gain = sum(sig_hyper & d_beta >= 0.1),
  .groups = "drop")
cat("\nC. Promoter-CGI methylation change vs expression change (lesion minus own adjacent; genes with TPM >= 1)\n")
print(mxs |> mutate(across(where(is.double), ~ signif(.x, 3))), width = 220)
write_csv(mx, TAB("rnaseq_meth_expr_pairs.csv")); write_csv(mxs, TAB("rnaseq_meth_expr_summary.csv"))

## ------------------------------------------------------------------ D. RNA-inferred arm dosage
arms <- sesameData_getGenomeInfo("hg38")$cytoBand |> mutate(across(c(chrom, name), as.character)) |>
  filter(chrom %in% paste0("chr", 1:22)) |> mutate(arm = substr(name, 1, 1)) |>
  group_by(chrom, arm) |> summarise(start = min(chromStart), end = max(chromEnd), .groups = "drop") |>
  mutate(arm_label = paste0(sub("chr", "", chrom), arm))
tx <- sesameData_getTxnGRanges("hg38")
gp <- tibble(gene = mcols(tx)$gene_name, chr = as.character(seqnames(tx)), pos = start(tx)) |>
  filter(chr %in% paste0("chr", 1:22), gene %in% rownames(lt)) |> group_by(gene) |>
  summarise(chr = first(chr), pos = min(pos), .groups = "drop") |> filter(!(chr == "chr6" & pos > 28e6 & pos < 34e6))
gp$arm <- map2_chr(gp$chr, gp$pos, function(ch, q) { a <- arms$arm_label[arms$chrom == ch & arms$start <= q & arms$end > q]; if (length(a)) a[1] else NA })
gp <- gp[rowMeans(lt[gp$gene, ] > 2) == 1 & !is.na(gp$arm), ]
ua <- names(which(table(gp$arm) >= 30))
adjs <- grep("Adj|NAT", colnames(lt), value = TRUE)
dos <- map_dfr(colnames(lt), function(s) {
  ref <- setdiff(adjs, grep(paste0("^", sub("-.*", "", s), "-"), adjs, value = TRUE))
  d <- lt[gp$gene, s] - apply(lt[gp$gene, ref], 1, median)
  a <- sapply(ua, function(x) mean(d[gp$arm == x])); tibble(sample = s, arm_label = ua, dosage = a - mean(a)) })
## noise floor from adjacent tissues of other patients: 99th percentile of |dosage| across adjacent arms
floor99 <- quantile(abs(dos$dosage[dos$sample %in% adjs]), 0.99)
loh_arms <- c("3p", "9p", "13q", "4q", "8p", "11q", "17p")
cat(sprintf("\nD. RNA arm dosage (vs other patients' adjacent); noise floor (99th pct of adjacent |dosage|) = %.3f\n", floor99))
print(dos |> filter(arm_label %in% loh_arms) |> mutate(dosage = round(dosage, 2)) |>
        pivot_wider(names_from = arm_label, values_from = dosage), width = 200)
write_csv(dos |> mutate(beyond_noise = abs(dosage) > floor99), TAB("rnaseq_arm_dosage_all.csv"))

## ------------------------------------------------------------------ E. external test of the pilot lesion programme
up <- de |> filter(gene_type == "protein_coding", logFC > 0) |> slice_min(PValue, n = 50) |> pull(gene)
dn <- de |> filter(gene_type == "protein_coding", logFC < 0) |> slice_min(PValue, n = 50) |> pull(gene)
## GSE26549
meta <- readLines(gzfile(file.path(XDATA, "GSE26549_series_matrix.txt.gz")))
ch <- grep("^!Sample_characteristics_ch1", meta, value = TRUE)
gsm <- gsub('"', "", strsplit(grep("^!Sample_geo_accession", meta, value = TRUE), "\t")[[1]][-1])
cl <- map(ch, ~ gsub('"', "", strsplit(.x, "\t")[[1]][-1]))
key <- map_chr(cl, ~ sub(":.*", "", .x[1])); cl <- setNames(map(cl, ~ sub("^[^:]*: ", "", .x)), make.names(key))
clin <- as_tibble(cl) |> mutate(gsm = gsm)
names(clin) <- c("sex", "race", "alcohol", "smoking", "age", "arm", "histology", "histology_detail", "p63", "podoplanin",
                 "ocfs_years", "outcome", "biopsy_time", "gsm")
clin <- clin |> mutate(ocfs_years = as.numeric(ocfs_years), event = as.integer(outcome == "oral cancer development"),
                       dysplasia = as.integer(histology == "dysplasia"))
b <- grep("series_matrix_table_begin", meta); e <- grep("series_matrix_table_end", meta)
X <- as.matrix(read.delim(text = meta[(b + 1):(e - 1)], check.names = FALSE, row.names = 1))[, clin$gsm]
sx <- mapIds(hugene10sttranscriptcluster.db, rownames(X), "SYMBOL", "PROBEID", multiVals = "first")
X <- X[!is.na(sx), ]; sx <- sx[!is.na(sx)]
kx <- tapply(seq_along(sx), sx, function(i) i[which.max(rowMeans(X[i, , drop = FALSE]))])
X <- X[unlist(kx), ]; rownames(X) <- names(kx); Zx <- t(scale(t(X)))
cox_hr <- function(sc, adj = FALSE) { x <- as.numeric(scale(sc))
  f <- if (adj) coxph(Surv(clin$ocfs_years, clin$event) ~ x + clin$dysplasia) else coxph(Surv(clin$ocfs_years, clin$event) ~ x)
  s <- summary(f); c(HR = s$coefficients[1, 2], lo = s$conf.int[1, 3], hi = s$conf.int[1, 4], p = s$coefficients[1, 5]) }
u26 <- intersect(up, rownames(Zx)); d26 <- intersect(dn, rownames(Zx))
prog26 <- colMeans(Zx[u26, ]) - colMeans(Zx[d26, ])
obs <- cox_hr(prog26); obs_adj <- cox_hr(prog26, adj = TRUE)
pool <- intersect(de$gene[de$gene_type == "protein_coding"], rownames(Zx))
null_p <- replicate(1000, { r <- sample(pool, length(u26) + length(d26))
  cox_hr(colMeans(Zx[r[seq_along(u26)], ]) - colMeans(Zx[r[-seq_along(u26)], ]))["p"] })
emp <- (sum(null_p <= obs["p"]) + 1) / 1001
cat(sprintf("\nE. Pilot lesion programme (%d up / %d down genes found) in GSE26549: HR per s.d. %.2f (%.2f-%.2f), P = %.3g; histology-adjusted HR %.2f, P = %.3g; random-set empirical P = %.3f\n",
            length(u26), length(d26), obs["HR"], obs["lo"], obs["hi"], obs["p"], obs_adj["HR"], obs_adj["p"], emp))
## proxy correlation in GSE26549 (is the pilot programme the same axis as the methylome proxy?)
pg <- readLines(TAB("methylome_state_expression_proxy_genes.txt"))
pup <- trimws(strsplit(sub("^UP: ", "", grep("^UP", pg, value = TRUE)), ",")[[1]])
pdn <- trimws(strsplit(sub("^DOWN: ", "", grep("^DOWN", pg, value = TRUE)), ",")[[1]])
proxy26 <- colMeans(Zx[intersect(pup, rownames(Zx)), ]) - colMeans(Zx[intersect(pdn, rownames(Zx)), ])
rho_prog_proxy <- cor(prog26, proxy26, method = "spearman")
cat(sprintf("   Spearman with the methylome-state proxy in GSE26549: %.2f; overlap of gene lists: %d up, %d down\n",
            rho_prog_proxy, length(intersect(up, pup)), length(intersect(dn, pdn))))
## GSE227919
ex <- read.table(gzfile(file.path(RNA_DIR, "GSE227919_2023_20_03_expr_mat_raw.txt.gz")), header = TRUE, check.names = FALSE)
a2 <- read.csv(gzfile(file.path(RNA_DIR, "GSE227919_2023_21_03_pml_annotations.csv.gz")), check.names = FALSE, na.strings = c("#N/A", "NA", ""))
names(a2)[1] <- "Sample"; a2 <- a2[match(colnames(ex), a2$Sample), ]
a2$Class <- factor(a2$Class, levels = c("Control", "HkNR", "Dysplasia", "Cancer"))
ex <- log2(as.matrix(ex) + 1); ex <- ex[rowMeans(ex > 1) >= 0.2, ]; z2 <- t(scale(t(ex)))
prog27 <- colMeans(z2[intersect(up, rownames(z2)), ]) - colMeans(z2[intersect(dn, rownames(z2)), ])
gr <- cor.test(prog27, as.integer(a2$Class), method = "spearman", exact = FALSE)
pr <- a2$Class != "Cancer" & !is.na(a2$Progression_status); yv <- a2$Progression_status[pr] != "Stable"
auc27 <- mean(outer(prog27[pr][yv], prog27[pr][!yv], ">") + 0.5 * outer(prog27[pr][yv], prog27[pr][!yv], "=="))
cat(sprintf("   GSE227919: Spearman with histological class %.2f (P = %.2g); progression AUC %.2f (%d of %d progressed)\n",
            gr$estimate, gr$p.value, auc27, sum(yv), length(yv)))
outE <- tibble(cohort = c("GSE26549", "GSE26549 (adj. histology)", "GSE227919 grade", "GSE227919 progression"),
               estimate = c(obs["HR"], obs_adj["HR"], gr$estimate, auc27), p = c(obs["p"], obs_adj["p"], gr$p.value, NA),
               empirical_p_random_sets = c(emp, NA, NA, NA))
write_csv(outE, TAB("rnaseq_lesion_programme_external.csv"))
writeLines(c(paste("UP:", paste(up, collapse = ", ")), paste("DOWN:", paste(dn, collapse = ", "))), TAB("rnaseq_lesion_programme_genes.txt"))

## ------------------------------------------------------------------ F. robustness
cmp <- cbind(Immune = rowMeans(score[, c("T_cell", "B_cell", "Myeloid", "Plasma_cell")]), CAF = score[, "CAF"], Endo = score[, "Endothelial"])
gg <- unique(mx$gene)
rr <- sapply(colnames(cmp), function(k) apply(lt[gg, rownames(cmp)], 1, function(v) suppressWarnings(cor(v, cmp[, k], method = "spearman"))))
nonepi <- gg[apply(rr, 1, function(v) any(v > 0.5, na.rm = TRUE))]
mxs_c <- mx |> filter(!gene %in% nonepi) |> group_by(sample) |> summarise(n_genes = n(), n_gain = sum(d_beta >= 0.2),
  med_expr_gain = median(d_expr[d_beta >= 0.2]), med_expr_stable = median(d_expr[abs(d_beta) < 0.05]),
  p_gain_vs_stable = wilcox.test(d_expr[d_beta >= 0.2], d_expr[abs(d_beta) < 0.05])$p.value, .groups = "drop")
cat(sprintf("\nF. Coupling after removing %d genes that track immune/CAF/endothelial content (rho > 0.5):\n", length(nonepi)))
print(mxs_c |> mutate(across(where(is.double), ~ signif(.x, 3))))
shared <- mx |> filter(sample %in% c("P4-Lesion", "P5-Lesion")) |> select(sample, gene, d_beta, d_expr) |>
  pivot_wider(names_from = sample, values_from = c(d_beta, d_expr)) |>
  filter(`d_beta_P4-Lesion` >= 0.2, `d_beta_P5-Lesion` >= 0.2) |>
  mutate(down_both = `d_expr_P4-Lesion` < 0 & `d_expr_P5-Lesion` < 0)
cat(sprintf("Promoters gaining >= 0.2 in both P4 and P5 lesions: %d; expression lower in both lesions: %d\n", nrow(shared), sum(shared$down_both)))
muscle <- intersect(c("ACTA1", "MYH2", "MYH7", "CKM", "DES"), rownames(lt))
muscle_score <- colMeans(lt[muscle, ])
cat("Skeletal-muscle expression (mean log2 TPM+1 of ACTA1, MYH2, MYH7, CKM, DES):\n"); print(round(sort(muscle_score), 1))
write_csv(mxs_c, TAB("rnaseq_meth_expr_summary_composition_controlled.csv")); write_csv(shared, TAB("rnaseq_shared_hypermethylated_P4_P5.csv"))

## ------------------------------------------------------------------ G. chr19 KZFP enrichment
gpos19 <- tibble(gene = mcols(tx)$gene_name, chr = as.character(seqnames(tx)), pos = start(tx)) |>
  group_by(gene) |> summarise(chr = first(chr), pos = min(pos), .groups = "drop")
is_kz <- function(g) { q <- gpos19[match(g, gpos19$gene), ]; !is.na(q$chr) & q$chr == "chr19" & grepl("^(ZNF|ZFP|ZSCAN|ZIK|ZKSCAN|ZIM)", g) }
univ <- mx |> filter(sample %in% c("P4-Lesion", "P5-Lesion")) |> count(gene) |> filter(n == 2) |> pull(gene)
kz_univ <- is_kz(univ)
f_gain <- fisher.test(table(univ %in% shared$gene, kz_univ)); f_dn <- fisher.test(table(univ %in% shared$gene[shared$down_both], kz_univ))
kz_shared <- shared |> mutate(kzfp = is_kz(gene), chr19_mb = gpos19$pos[match(gene, gpos19$gene)] / 1e6, composition_tracking = gene %in% nonepi)
cat(sprintf("\nG. Promoters tested in both premalignant lesions: %d (%d chr19 KZFP). Shared gain >= 0.2: %d (%d KZFP; OR %.1f, P = %.2g).\n   Gain and lower expression in both: %d (%d KZFP; OR %.1f, P = %.2g). KZFPs lower in both: %d of %d; other genes: %d of %d.\n   Composition-tracking genes among the gained-and-silenced KZFPs: %d\n",
            length(univ), sum(kz_univ), nrow(shared), sum(kz_shared$kzfp), f_gain$estimate, f_gain$p.value,
            sum(shared$down_both), sum(kz_shared$kzfp & kz_shared$down_both), f_dn$estimate, f_dn$p.value,
            sum(kz_shared$kzfp & kz_shared$down_both), sum(kz_shared$kzfp), sum(!kz_shared$kzfp & kz_shared$down_both), sum(!kz_shared$kzfp),
            sum(kz_shared$kzfp & kz_shared$down_both & kz_shared$composition_tracking)))
print(kz_shared |> filter(down_both) |> arrange(desc(kzfp), chr19_mb) |> select(gene, kzfp, chr19_mb, starts_with("d_")) |> mutate(across(where(is.double), ~ round(.x, 2))), n = 30)
write_csv(kz_shared, TAB("rnaseq_shared_hypermethylated_P4_P5_kzfp.csv"))
pair_mean <- mx |> filter(sample %in% c("P4-Lesion", "P5-Lesion"), gene %in% univ) |> group_by(gene) |>
  summarise(d_beta = mean(d_beta), d_expr = mean(d_expr), .groups = "drop") |> mutate(kzfp = is_kz(gene), shared_gain = gene %in% shared$gene)

saveRDS(list(ss = ss, lt = lt, tpm = tpm, score = score, progA = progA, progA_delta = progA_delta, de = de, cam = cam,
             mx = mx, mxs = mxs, dos = dos, floor99 = floor99, up = up, dn = dn, prog26 = prog26, proxy26 = proxy26,
             clin = clin, obs = obs, obs_adj = obs_adj, null_p = null_p, emp = emp, prog27 = prog27, a2 = a2[, c("Sample", "Class", "Progression_status")],
             gr = gr, auc27 = auc27, rho_prog_proxy = rho_prog_proxy,
             nonepi = nonepi, mxs_c = mxs_c, shared = shared, muscle_score = muscle_score,
             kz_shared = kz_shared, f_gain = f_gain, f_dn = f_dn, univ_n = length(univ), univ_kz = sum(kz_univ), pair_mean = pair_mean), file.path(EXT, "rnaseq_story.rds"))
cat("\nSaved", file.path(EXT_DIR, "rnaseq_story.rds"), "\n")
