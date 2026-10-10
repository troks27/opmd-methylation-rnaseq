## Pilot RNA-seq (kallisto -> gene matrices from 02_gene_matrices.R) integrated with the pilot methylation results.
## Descriptive only (4 methylation patients, 1 progressor). Run from the repository root:  Rscript scripts/05_pilot_rnaseq/03_rna_methylation_validation.R  (or make rnaseq)
##   A. QC and sex (all pilot samples female by array)
##   B. composition: RNA marker scores vs EpiDISH fractions; the P3 label question
##   C. 9p21: MTAP / CDKN2A / CDKN2B expression where the arrays call deletion
##   D. arm-level expression dosage vs methylation-array copy number (paired to own adjacent)
##   E. KZFP promoter methylation change vs expression change, lesion vs own adjacent
##   F. the GSE26549 methylome-state expression proxy and spatial programmes on the pilot
##   G. arm dosage vs array CN against OTHER patients' adjacent tissue (shared field events, e.g. P4 chr12)
##   H. composition checks: proxy vs immune content; EpiDISH Epi vs non-epithelial RNA scores
source("R/lib/config.R")   # paths and helpers (TAB, FIG, OUT_DIR, GSE_DIR, ...)
suppressPackageStartupMessages({
  library(tidyverse); library(sesameData); library(GenomicRanges)
})
ss <- read_tsv(file.path(MAT, "samples.tsv"), show_col_types = FALSE)
ann <- read_tsv(file.path(MAT, "gene_annotation.tsv"), show_col_types = FALSE)
tpm <- read_tsv(file.path(MAT, "gene_tpm.tsv.gz"), show_col_types = FALSE) |> column_to_rownames("gene_id") |> as.matrix()
## one row per symbol: the most expressed gene_id
sym <- ann$gene_name[match(rownames(tpm), ann$gene_id)]
keep <- tapply(seq_along(sym), sym, function(i) i[which.max(rowMeans(tpm[i, , drop = FALSE]))])
tpm <- tpm[unlist(keep), ]; rownames(tpm) <- names(keep)
lt <- log2(tpm + 1)
short <- setNames(paste0(ss$patient, "-", c(Tumor = "Tumor", Lesion = "Lesion", Adjacent = "Adj")[ss$tissue]), ss$library)
short[short == "P2-Adj"] <- "P2-NAT"
colnames(lt) <- short[colnames(lt)]
adj_of <- function(s) { p <- sub("-.*", "", s); a <- if (p == "P2") "P2-NAT" else paste0(p, "-Adj"); if (a %in% colnames(lt)) a else NA }
les <- colnames(lt)[!grepl("Adj|NAT", colnames(lt))]

## ---------------------------------------------------------------- A. QC / sex
sexg <- c("XIST", "RPS4Y1", "DDX3Y", "UTY", "KDM5D", "EIF1AY")
qcA <- ss |> transmute(sample = short[library], library, n_pairs, pct_pseudoaligned) |>
  left_join(as_tibble(t(lt[intersect(sexg, rownames(lt)), ]), rownames = "sample"), by = "sample")
cat("A. QC and sex genes (log2 TPM+1)\n"); print(qcA |> mutate(across(where(is.double), ~ round(.x, 2))), width = 200)
write_csv(qcA, TAB("rnaseq_qc_sex.csv"))

## ---------------------------------------------------------------- B. composition
markers <- list(Epithelial = c("KRT5", "KRT14", "KRT13", "KRT4", "TP63", "SPRR3", "KRT6A"),
                Immune = c("PTPRC", "CD3E", "CD2", "MS4A1", "CD79A", "LYZ", "CD68", "CD14"),
                Lymphoid = c("CD3E", "CD3D", "CD2", "MS4A1", "CD79A", "CD19", "IGHM"),
                Fibroblast = c("COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA", "COL3A1"),
                Endothelial = c("PECAM1", "VWF", "CDH5", "CLDN5"))
zs <- t(scale(t(lt)))
comp <- sapply(markers, function(g) colMeans(zs[intersect(g, rownames(zs)), , drop = FALSE]))
epi <- read_csv(TAB("deconvolution_epidish.csv"), show_col_types = FALSE)
compB <- as_tibble(comp, rownames = "sample") |> left_join(epi |> select(sample = Sample, Epi_meth = Epi_1st, Fib_meth = Fib_1st, IC_meth = IC_1st), by = "sample")
cat("\nB. RNA composition scores (z) vs methylation EpiDISH fractions\n"); print(compB |> mutate(across(where(is.double), ~ round(.x, 2))))
mm <- compB |> filter(!is.na(Epi_meth))
cat(sprintf("Spearman RNA epithelial vs EpiDISH Epi: %.2f; RNA immune vs EpiDISH IC: %.2f (n = %d)\n",
            cor(mm$Epithelial, mm$Epi_meth, method = "spearman"), cor(mm$Immune, mm$IC_meth, method = "spearman"), nrow(mm)))
cat(sprintf("P3-Lesion immune score rank: %d of %d (methylation: %.0f%% immune); lymphoid rank %d of %d\n",
            rank(-compB$Immune)[compB$sample == "P3-Lesion"], nrow(compB), 100 * compB$IC_meth[compB$sample == "P3-Lesion"],
            rank(-compB$Lymphoid)[compB$sample == "P3-Lesion"], nrow(compB)))
write_csv(compB, TAB("rnaseq_composition_vs_epidish.csv"))

## ---------------------------------------------------------------- C. 9p21
cn9 <- read_csv(TAB("cnv_9p21_detection_check.csv"), show_col_types = FALSE)
g9 <- c("MTAP", "CDKN2A", "CDKN2B", "DMRTA1", "IFNB1")
c9 <- map_dfr(les, function(s) { a <- adj_of(s); if (is.na(a)) return(NULL)
  tibble(sample = s, gene = g9, log2FC_vs_own_adj = lt[g9, s] - lt[g9, a], tpm_lesion = tpm[g9, names(short)[short == s]]) }) |>
  left_join(cn9 |> select(sample = Sample, cn_9p21_log2 = lpo_9p21_median, probe_fail_ratio = fail_ratio), by = "sample")
cat("\nC. 9p21 genes, lesion vs own adjacent (log2 fold change) with array copy number\n")
print(c9 |> select(sample, gene, log2FC_vs_own_adj, cn_9p21_log2) |> pivot_wider(names_from = gene, values_from = log2FC_vs_own_adj) |>
        mutate(across(where(is.double), ~ round(.x, 2))))
write_csv(c9, TAB("rnaseq_9p21.csv"))

## ---------------------------------------------------------------- D. arm dosage vs array CN
arms <- sesameData_getGenomeInfo("hg38")$cytoBand |> mutate(across(c(chrom, name), as.character)) |>
  filter(chrom %in% paste0("chr", 1:22)) |> mutate(arm = substr(name, 1, 1)) |>
  group_by(chrom, arm) |> summarise(start = min(chromStart), end = max(chromEnd), .groups = "drop") |>
  mutate(arm_label = paste0(sub("chr", "", chrom), arm))
tx <- sesameData_getTxnGRanges("hg38")
gpos <- tibble(gene = mcols(tx)$gene_name, chr = as.character(seqnames(tx)), pos = start(tx)) |>
  filter(chr %in% paste0("chr", 1:22), gene %in% rownames(lt)) |>
  group_by(gene) |> summarise(chr = dplyr::first(chr), pos = min(pos), .groups = "drop") |>
  filter(!(chr == "chr6" & pos > 28e6 & pos < 34e6))
gpos$arm_label <- map2_chr(gpos$chr, gpos$pos, function(ch, p) {
  a <- arms$arm_label[arms$chrom == ch & arms$start <= p & arms$end > p]; if (length(a)) a[1] else NA })
expr_ok <- rowMeans(lt[gpos$gene, ] > 2) == 1               # expressed in every library
gpos <- gpos[expr_ok & !is.na(gpos$arm_label), ]
use_arms <- names(which(table(gpos$arm_label) >= 30))
dosD <- map_dfr(les, function(s) { a <- adj_of(s); if (is.na(a)) return(NULL)
  d <- lt[gpos$gene, s] - lt[gpos$gene, a]
  arm <- sapply(use_arms, function(x) mean(d[gpos$arm_label == x]))
  tibble(sample = s, arm_label = use_arms, rna_dosage = arm - mean(arm)) })
cna <- read_csv(TAB("cnv_arm_level.csv"), show_col_types = FALSE) |> select(sample = Sample, arm_label, cn_paired = paired)
dosD <- dosD |> left_join(cna, by = c("sample", "arm_label"))
cat("\nD. Expression dosage vs methylation-array paired copy number, per lesion (Pearson across arms)\n")
print(dosD |> filter(!is.na(cn_paired)) |> group_by(sample) |>
        summarise(n_arms = n(), r = cor(rna_dosage, cn_paired), .groups = "drop") |> mutate(r = round(r, 2)))
focus <- tibble(sample = c("P5-Lesion", "P5-Lesion", "P4-Lesion", "P4-Lesion", "P2-Tumor", "P2-Tumor", "P2-Tumor"),
                arm_label = c("3p", "13q", "12p", "12q", "20p", "20q", "4q"))
print(focus |> left_join(dosD, by = c("sample", "arm_label")) |> mutate(across(where(is.double), ~ round(.x, 3))))
write_csv(dosD, TAB("rnaseq_arm_dosage_vs_cn.csv"))

## ---------------------------------------------------------------- E. KZFP methylation vs expression
kz <- read_csv(TAB("kzfp_pilot_deltas.csv"), show_col_types = FALSE) |>
  pivot_longer(starts_with("d_"), names_to = "pair", values_to = "delta_beta") |>
  mutate(sample = c(d_P2_tumor = "P2-Tumor", d_P3 = "P3-Lesion", d_P4 = "P4-Lesion", d_P5_progressor = "P5-Lesion")[pair])
kz$delta_log2_tpm <- map2_dbl(kz$gene, kz$sample, function(g, s) { a <- adj_of(s)
  if (is.na(a) || !g %in% rownames(lt)) NA else lt[g, s] - lt[g, a] })
cat("\nE. KZFP promoter delta-beta vs expression log2 fold change (lesion vs own adjacent)\n")
print(kz |> select(gene, sample, delta_beta, delta_log2_tpm) |> mutate(across(where(is.double), ~ round(.x, 2))), n = 30)
ok <- !is.na(kz$delta_log2_tpm)
cat(sprintf("Spearman across gene x lesion pairs: %.2f (n = %d)\n", cor(kz$delta_beta[ok], kz$delta_log2_tpm[ok], method = "spearman"), sum(ok)))
write_csv(kz, TAB("rnaseq_kzfp_meth_vs_expr.csv"))

## ---------------------------------------------------------------- F. proxy and programmes
pg <- readLines(TAB("methylome_state_expression_proxy_genes.txt"))
up <- trimws(strsplit(sub("^UP: ", "", grep("^UP", pg, value = TRUE)), ",")[[1]])
dn <- trimws(strsplit(sub("^DOWN: ", "", grep("^DOWN", pg, value = TRUE)), ",")[[1]])
progs <- list(Invasive_epithelial = c("LAMC2", "COL17A1", "PTHLH", "ITGA6", "ITGB4", "HMGA2"),
              Cycling_epithelial = c("HMMR", "ASPM", "BUB1", "NEIL3", "BUB1B", "PBK"),
              LRRC15_CAF = c("LRRC15", "FAP", "COL11A1"), FOLR2_TAM = c("FOLR2", "CD163", "MRC1"))
sc <- tibble(sample = colnames(zs),
             proxy = colMeans(zs[intersect(up, rownames(zs)), ]) - colMeans(zs[intersect(dn, rownames(zs)), ])) |>
  bind_cols(as_tibble(sapply(progs, function(g) colMeans(zs[intersect(g, rownames(zs)), , drop = FALSE]))))
tl <- read_csv(TAB("tumor_likeness_pilot.csv"), show_col_types = FALSE) |> select(sample = Sample, tumor_likeness)
scF <- sc |> left_join(tl, by = "sample")
cat(sprintf("\nF. Proxy genes found: %d/%d up, %d/%d down\n", sum(up %in% rownames(zs)), length(up), sum(dn %in% rownames(zs)), length(dn)))
print(scF |> arrange(desc(proxy)) |> mutate(across(where(is.double), ~ round(.x, 2))))
m2 <- scF |> filter(!is.na(tumor_likeness))
cat(sprintf("Spearman proxy vs methylation tumour-likeness (n = %d): %.2f\n", nrow(m2), cor(m2$proxy, m2$tumor_likeness, method = "spearman")))
write_csv(scF, TAB("rnaseq_proxy_programmes.csv"))

## ---------------------------------------------------------------- G. leave-patient-out dosage
adjs <- grep("Adj|NAT", colnames(lt), value = TRUE)
lpo <- map_dfr(colnames(lt), function(s) {
  ref <- setdiff(adjs, grep(paste0("^", sub("-.*", "", s), "-"), adjs, value = TRUE))
  d <- lt[gpos$gene, s] - apply(lt[gpos$gene, ref], 1, median)
  a <- sapply(use_arms, function(x) mean(d[gpos$arm_label == x])); tibble(sample = s, arm_label = use_arms, rna_lpo = a - mean(a)) })
lpo <- lpo |> left_join(read_csv(TAB("cnv_arm_level.csv"), show_col_types = FALSE) |> select(sample = Sample, arm_label, cn_lpo = LPO),
                        by = c("sample", "arm_label"))
cat("\nG. RNA vs array copy number against other patients' adjacent tissue (Pearson across arms)\n")
print(lpo |> filter(!is.na(cn_lpo)) |> group_by(sample) |> summarise(r = round(cor(rna_lpo, cn_lpo), 2), .groups = "drop"))
print(lpo |> filter(arm_label %in% c("12p", "12q"), !is.na(cn_lpo)) |> mutate(across(where(is.double), ~ round(.x, 3))))
write_csv(lpo, TAB("rnaseq_arm_dosage_vs_cn_lpo.csv"))

## ---------------------------------------------------------------- H. composition checks
x <- scF |> left_join(compB, by = "sample")
m3 <- x |> filter(!is.na(tumor_likeness))
fitH <- lm(rank(proxy) ~ rank(tumor_likeness) + rank(Immune), data = m3)
cat(sprintf("\nH. proxy vs RNA immune score rho %.2f (n = %d); proxy ~ tumour-likeness given immune: rank coef %.2f, p %.3f (n = %d)\n",
            cor(x$proxy, x$Immune, method = "spearman"), nrow(x), coef(fitH)[2], summary(fitH)$coefficients[2, 4], nrow(m3)))
cat(sprintf("   EpiDISH Epi vs -(immune + fibroblast) RNA score rho %.2f; vs keratinocyte-marker score rho %.2f (n = %d)\n",
            cor(mm$Epi_meth, -(mm$Immune + mm$Fibroblast), method = "spearman"), cor(mm$Epi_meth, mm$Epithelial, method = "spearman"), nrow(mm)))
