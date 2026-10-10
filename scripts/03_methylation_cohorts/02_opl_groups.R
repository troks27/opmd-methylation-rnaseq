## Follow-up tests behind the central story:
##   A. GSE204943 OPLs: do copy-number losses at the classic oral-LOH arms
##      (3p, 9p, 13q, 4q, 8p, 11q, 17p) travel with the tumour-like methylome?
##      (conumee segments from the GSE204943 project, hg19)
##   B. GSE227919: arm-level expression dosage as a copy-number proxy.
##      Positive control: cancer vs control should show the known HNSCC arm
##      pattern (3p, 9p, 4q loss; 3q, 8q, 20 gain). Then: does 3p/9p/13q
##      dosage loss mark the lesions that progressed?
##
## Run from the repository root:  Rscript scripts/03_methylation_cohorts/02_opl_groups.R  (or make cohorts)

source("R/lib/config.R")   # paths and helpers (TAB, FIG, OUT_DIR, GSE_DIR, ...)
suppressPackageStartupMessages({
  library(tidyverse); library(sesameData); library(GenomicRanges)
})
LOH_ARMS <- c("3p", "9p", "13q", "4q", "8p", "11q", "17p")

arm_table <- function(build) {
  sesameData_getGenomeInfo(build)$cytoBand |>
    mutate(across(c(chrom, name), as.character)) |>
    filter(chrom %in% paste0("chr", 1:22)) |>
    mutate(arm = substr(name, 1, 1)) |>
    group_by(chrom, arm) |>
    summarise(start = min(chromStart), end = max(chromEnd), .groups = "drop") |>
    mutate(arm_label = paste0(sub("chr", "", chrom), arm))
}

## ================================================================== A
arms19 <- arm_table("hg19")
segs <- read_csv(file.path(GSE_DIR, "results/CNV_segments_all.csv"), show_col_types = FALSE)
seg_gr <- GRanges(segs$chrom, IRanges(segs$loc.start, segs$loc.end))
arm_gr <- GRanges(arms19$chrom, IRanges(arms19$start, arms19$end))
ov <- findOverlaps(seg_gr, arm_gr)
w <- width(pintersect(seg_gr[queryHits(ov)], arm_gr[subjectHits(ov)]))
gse_arm <- tibble(Sample = segs$ID[queryHits(ov)], Group = segs$Group[queryHits(ov)],
                  arm_label = arms19$arm_label[subjectHits(ov)],
                  seg.mean = segs$seg.mean[queryHits(ov)], w = w) |>
  group_by(Sample, Group, arm_label) |>
  summarise(arm_lr = weighted.mean(seg.mean, w), .groups = "drop")
## conumee centres each profile; re-centre on the sample's median arm
gse_arm <- gse_arm |> group_by(Sample) |> mutate(arm_lr = arm_lr - median(arm_lr)) |> ungroup()
gse_wide <- gse_arm |> select(-Group) |> pivot_wider(names_from = arm_label, values_from = arm_lr)

sc_g <- read_csv(TAB("tumor_likeness_gse204943.csv"), show_col_types = FALSE)
ext <- readRDS(file.path(EXT_DIR, "external_methylation.rds"))
pz_g <- ext$pz_g |> select(Sample, pmd_contrast, ZNF471_promoter)
gse <- sc_g |> left_join(pz_g, by = "Sample") |> left_join(gse_wide, by = "Sample") |>
  mutate(core_loss = -(`3p` + `9p`) / 2,                       # higher = more 3p/9p loss
         loh_arms_lost = rowSums(across(all_of(LOH_ARMS), ~ .x <= -0.1)),
         any_3p9p_loss = `3p` <= -0.1 | `9p` <= -0.1)
write_csv(gse |> mutate(across(where(is.numeric), ~ round(.x, 3))), TAB("gse204943_cn_vs_methylome.csv"))

opl <- gse |> filter(Group == "OPL")
cat(sprintf("GSE204943 OPL with CN segments: %d of 22\n", sum(!is.na(opl$`3p`))))
## OPL methylome is bimodal: call 'tumour-like' at the valley between modes
dens <- density(opl$tumor_likeness, bw = 0.15)
valley <- with(dens, x[x > 0.1 & x < 0.8][which.min(y[x > 0.1 & x < 0.8])])
opl <- opl |> mutate(methylome_state = ifelse(tumor_likeness > valley, "tumour-like", "normal-like"))
cat(sprintf("OPL tumour-likeness valley = %.2f; tumour-like %d / %d\n",
            valley, sum(opl$methylome_state == "tumour-like"), nrow(opl)))
cat("Epithelial fraction by OPL methylome state:\n")
print(opl |> group_by(methylome_state) |> summarise(n = n(), Epi = round(median(Epi), 2),
        IC = NA, loh_arms = round(mean(loh_arms_lost, na.rm = TRUE), 2),
        any_3p9p = sum(any_3p9p_loss, na.rm = TRUE), pmd = round(median(pmd_contrast), 3),
        znf471 = round(median(ZNF471_promoter), 3)))
## association tests inside OPL (n = 22) - Spearman, and partial on Epi
pcor <- function(x, y, z) {
  rx <- resid(lm(rank(x) ~ rank(z))); ry <- resid(lm(rank(y) ~ rank(z)))
  ct <- cor.test(rx, ry); c(r = unname(ct$estimate), p = ct$p.value)
}
assoc <- map_dfr(c(LOH_ARMS, "core_loss", "loh_arms_lost"), function(a) {
  d <- opl |> filter(!is.na(.data[[a]]))
  ct <- cor.test(d$tumor_likeness, d[[a]], method = "spearman", exact = FALSE)
  pc <- pcor(d$tumor_likeness, d[[a]], d$Epi)
  tibble(feature = a, n = nrow(d), rho = unname(ct$estimate), p = ct$p.value,
         partial_rho_given_Epi = pc["r"], partial_p = pc["p"])
})
write_csv(assoc |> mutate(across(where(is.numeric), ~ signif(.x, 3))), TAB("gse204943_opl_cn_methylome_assoc.csv"))
cat("\nOPL: methylome tumour-likeness vs arm copy number (negative rho for an arm = more loss in tumour-like OPL):\n")
print(assoc |> mutate(across(where(is.numeric), ~ signif(.x, 2))))
## do the two GSE204943 programs (PMD erosion, ZNF471) co-occur within OPL?
pz_opl <- cor.test(opl$pmd_contrast, opl$ZNF471_promoter, method = "spearman", exact = FALSE)
pz_opl_adj <- pcor(opl$pmd_contrast, opl$ZNF471_promoter, opl$tumor_likeness)
cat(sprintf("OPL: PMD contrast vs ZNF471 rho = %.2f (p = %.2g); partial on tumour-likeness rho = %.2f (p = %.2g)\n",
            pz_opl$estimate, pz_opl$p.value, pz_opl_adj["r"], pz_opl_adj["p"]))
st_tests <- map_dfr(c("pmd_contrast", "ZNF471_promoter", "Epi", "loh_arms_lost", "core_loss"), function(v)
  tibble(variable = v, normal_like = median(opl[[v]][opl$methylome_state == "normal-like"], na.rm = TRUE),
         tumour_like = median(opl[[v]][opl$methylome_state == "tumour-like"], na.rm = TRUE),
         wilcox_p = wilcox.test(opl[[v]] ~ opl$methylome_state, exact = FALSE)$p.value))
write_csv(st_tests |> mutate(across(where(is.numeric), ~ signif(.x, 3))), TAB("gse204943_opl_state_tests.csv"))
print(st_tests |> mutate(across(where(is.numeric), ~ signif(.x, 2))))
wt <- wilcox.test(loh_arms_lost ~ methylome_state, data = opl, exact = FALSE)
cat(sprintf("LOH-arm losses, tumour-like vs normal-like OPL: Wilcoxon p = %.3g\n", wt$p.value))

## ================================================================== B
arms38 <- arm_table("hg38")
ex <- read.table(gzfile(file.path(RNA_DIR, "GSE227919_2023_20_03_expr_mat_raw.txt.gz")),
                 header = TRUE, check.names = FALSE)
ann <- read.csv(gzfile(file.path(RNA_DIR, "GSE227919_2023_21_03_pml_annotations.csv.gz")),
                check.names = FALSE, na.strings = c("#N/A", "NA", ""))
names(ann)[1] <- "Sample"
ann <- ann[match(colnames(ex), ann$Sample), ]
ann$Class <- factor(ann$Class, levels = c("Control", "HkNR", "Dysplasia", "Cancer"))
ann$Progressed <- ifelse(is.na(ann$Progression_status), NA, ann$Progression_status != "Stable")
lx <- log2(as.matrix(ex) + 1)
## dosage needs well-measured genes: mean log2 > 3 and expressed in 90% of samples
lx <- lx[rowMeans(lx) > 3 & rowMeans(lx > 1) >= 0.9, ]

tx <- sesameData_getTxnGRanges("hg38")
gpos <- tibble(gene = mcols(tx)$gene_name, chr = as.character(seqnames(tx)), pos = start(tx)) |>
  filter(chr %in% paste0("chr", 1:22)) |>
  group_by(gene) |> summarise(chr = dplyr::first(chr), pos = min(pos), .groups = "drop")
gpos <- gpos |> filter(gene %in% rownames(lx))
gpos$arm_label <- map2_chr(gpos$chr, gpos$pos, function(ch, p) {
  a <- arms38[arms38$chrom == ch & arms38$start <= p & arms38$end > p, "arm_label"]
  if (nrow(a)) a$arm_label[1] else NA_character_ })
gpos <- gpos |> filter(!is.na(arm_label))
## exclude the HLA/immune-dense 6p21 block, whose expression tracks infiltration
gpos <- gpos |> filter(!(chr == "chr6" & pos > 28e6 & pos < 34e6))
## z relative to the Control (benign mucosa / non-dysplastic lesion) median per gene,
## then arm score = mean z on the arm minus the sample's genome-wide mean z
ctrl <- ann$Class == "Control"
zz <- (lx[gpos$gene, ] - matrixStats::rowMedians(lx[gpos$gene, ctrl])) /
  pmax(matrixStats::rowMads(lx[gpos$gene, ctrl]), 0.1)
zz <- pmax(pmin(zz, 4), -4)
arm_n <- table(gpos$arm_label)
use_arms <- names(arm_n)[arm_n >= 40]
dos <- sapply(use_arms, function(a) colMeans(zz[gpos$gene[gpos$arm_label == a], , drop = FALSE]))
dos <- dos - rowMeans(dos)                                  # remove sample-wide shift
dos_df <- bind_cols(ann |> select(Sample, Class, Progression_status, Progressed), as_tibble(dos))
write_csv(dos_df |> mutate(across(where(is.numeric), ~ round(.x, 3))), TAB("gse227919_arm_dosage.csv"))

## positive control: cancer vs control
pc_tab <- map_dfr(use_arms, function(a) {
  w <- wilcox.test(dos_df[[a]][dos_df$Class == "Cancer"], dos_df[[a]][dos_df$Class == "Control"], exact = FALSE)
  tibble(arm = a, n_genes = as.integer(arm_n[a]),
         cancer_minus_control = median(dos_df[[a]][dos_df$Class == "Cancer"]) -
           median(dos_df[[a]][dos_df$Class == "Control"]), p = w$p.value)
}) |> mutate(fdr = p.adjust(p, "BH")) |> arrange(p)
write_csv(pc_tab |> mutate(across(where(is.numeric), ~ signif(.x, 3))), TAB("gse227919_arm_dosage_cancer_vs_control.csv"))
cat("\nGSE227919 arm dosage, cancer vs control (positive control; top 12):\n")
print(pc_tab |> head(12) |> mutate(across(where(is.numeric), ~ signif(.x, 2))))

## progression
auc <- function(x, y) { r <- rank(x); n1 <- sum(y); n0 <- sum(!y); (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
prog <- dos_df |> filter(Class != "Cancer", !is.na(Progressed)) |>
  mutate(core_loss = -(`3p` + `9p`) / 2,
         loh_loss = -rowMeans(across(all_of(intersect(LOH_ARMS, use_arms)))))
prog_tab <- map_dfr(c(intersect(LOH_ARMS, use_arms), "core_loss", "loh_loss"), function(a) {
  x <- if (a %in% c("core_loss", "loh_loss")) prog[[a]] else -prog[[a]]   # higher = more loss
  g <- glm(prog$Progressed ~ x + prog$Class, family = binomial)
  tibble(feature = a, AUC_loss_predicts_progression = auc(x, prog$Progressed),
         p = wilcox.test(x[prog$Progressed], x[!prog$Progressed], exact = FALSE)$p.value,
         p_adj_histology = summary(g)$coefficients[2, 4])
})
write_csv(prog_tab |> mutate(across(where(is.numeric), ~ signif(.x, 3))), TAB("gse227919_arm_dosage_progression.csv"))
cat(sprintf("\nGSE227919 progression (n = %d, %d events): arm-dosage LOSS as predictor\n", nrow(prog), sum(prog$Progressed)))
print(prog_tab |> mutate(across(where(is.numeric), ~ signif(.x, 2))))

saveRDS(list(gse = gse, opl = opl, valley = valley, assoc = assoc, wt_p = wt$p.value, st_tests = st_tests,
             pz_opl = c(rho = unname(pz_opl$estimate), p = pz_opl$p.value), pz_opl_adj = pz_opl_adj,
             dos_df = dos_df, pc_tab = pc_tab, prog = prog, prog_tab = prog_tab, use_arms = use_arms),
        file.path(EXT_DIR, "story_analyses.rds"))
cat("done\n")
