## A third methylation cohort and the chr19 KRAB-zinc-finger (KZFP) thread.
##
##   A. GSE234379 (450K; 46 matched OSCC / normal-adjacent pairs, Brazil):
##      does the GSE204943-derived, composition-adjusted tumour signature
##      transfer to an independent population and platform? How "normal" is
##      tumour-adjacent tissue (field)? KZFP promoter methylation and PMD erosion.
##   B. KZFP promoter methylation (ZNF471, ZNF382, ZNF529, ZNF582, ZNF671) by
##      OPL methylome state (GSE204943) and in the pilot pairs.
##
## Run from the repository root:  Rscript scripts/03_methylation_cohorts/04_gse234379_kzfp.R  (or make cohorts)

source("R/lib/config.R")   # paths and helpers (TAB, FIG, OUT_DIR, GSE_DIR, ...)
suppressPackageStartupMessages({
  library(tidyverse); library(EpiDISH); library(GenomicRanges); library(data.table)
})
select <- dplyr::select; filter <- dplyr::filter
KZFP <- c("ZNF471", "ZNF382", "ZNF529", "ZNF582", "ZNF671")

sig <- readRDS(file.path(EXT_DIR, "gse_signature.rds"))
ext <- readRDS(file.path(EXT_DIR, "external_methylation.rds"))
ann_g <- readRDS(file.path(GSE_DIR, "crucial_objects/ann_auto.rds"))
prom_cgi <- function(genes, ids) {
  a <- ann_g[intersect(ids, rownames(ann_g)), ]
  keep <- grepl("TSS200|TSS1500|1stExon", a$UCSC_RefGene_Group) & a$Relation_to_Island == "Island"
  gl <- strsplit(a$UCSC_RefGene_Name, ";")
  setNames(lapply(genes, function(g) rownames(a)[keep & vapply(gl, function(x) g %in% x, logical(1))]), genes)
}
score_fun <- function(B, cpgs) {
  cpgs <- intersect(cpgs, rownames(B))
  s <- sig[match(cpgs, sig$cpg), ]
  z <- (B[cpgs, , drop = FALSE] - s$muN) / (s$muT - s$muN)
  apply(z, 2, median, na.rm = TRUE)
}
pmd <- read_csv(file.path(GSE_DIR, "results/hypo_PMD_blocks.csv"), show_col_types = FALSE)
pmd_contrast <- function(B) {
  a <- ann_g[intersect(rownames(B), rownames(ann_g)), ]
  os <- a$Relation_to_Island == "OpenSea"
  inb <- overlapsAny(GRanges(a$chr, IRanges(a$pos, width = 1)), GRanges(pmd$seqnames, IRanges(pmd$start, pmd$end)))
  colMeans(B[rownames(a)[inb & os], ], na.rm = TRUE) - colMeans(B[rownames(a)[!inb & os], ], na.rm = TRUE)
}

## ================================================================== A
meta <- readLines(gzfile(file.path(XDATA, "GSE234379_series_matrix.txt.gz")))
grab <- function(tag, k = 1) gsub('"', "", strsplit(grep(tag, meta, value = TRUE)[k], "\t")[[1]][-1])
sm <- tibble(title = grab("^!Sample_title"), tissue = grab("^!Sample_description", 1),
             col = grab("^!Sample_description", 2), hpv = sub(".*: ", "", grab("hpv status", 1))) |>
  mutate(col = sub("^X", "", gsub("\\.", "-", col)), group = ifelse(grepl("Tumor", tissue), "Tumor", "Normal"),
         patient = sub("[NT]_.*", "", title))
dt <- fread(file.path(XDATA, "GSE234379_normalized_data.txt.gz"))
B3 <- as.matrix(dt[, -1]); rownames(B3) <- dt[[1]]; rm(dt)
B3 <- B3[, sm$col]; colnames(B3) <- paste0(sm$patient, "_", sm$group)
B3 <- B3[rowMeans(is.na(B3)) < 0.1, ]
cat(sprintf("GSE234379: %d CpGs x %d samples (%s); HPV+ %d\n", nrow(B3), ncol(B3),
            paste(names(table(sm$group)), table(sm$group), collapse = ", "), sum(sm$hpv == "Positive")))
pairs3 <- sm |> count(patient) |> filter(n == 2) |> pull(patient)
cat(sprintf("Complete tumour/adjacent pairs: %d\n", length(pairs3)))

data(centEpiFibIC.m)
cf3 <- epidish(beta.m = B3[complete.cases(B3), ], ref.m = centEpiFibIC.m, method = "RPC")$estF
th <- intersect(ext$tum_hyper, rownames(B3)); to <- intersect(ext$tum_hypo, rownames(B3))
cat(sprintf("Signature CpGs on 450K: %d of %d hyper, %d of %d hypo\n", length(th), length(ext$tum_hyper), length(to), length(ext$tum_hypo)))
s3 <- tibble(sample = colnames(B3), patient = sm$patient, group = sm$group, hpv = sm$hpv,
             hyper = score_fun(B3, th), hypo = score_fun(B3, to)) |>
  mutate(tumor_likeness = (hyper + hypo) / 2, Epi = cf3[sample, "Epi"], IC = cf3[sample, "IC"],
         pmd_contrast = pmd_contrast(B3)[sample])
kz3 <- prom_cgi(KZFP, rownames(B3))
for (g in KZFP) s3[[g]] <- if (length(kz3[[g]])) colMeans(B3[kz3[[g]], s3$sample, drop = FALSE], na.rm = TRUE) else NA
auc <- function(x, y) { r <- rank(x); n1 <- sum(y); n0 <- sum(!y); (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
auc3 <- auc(s3$tumor_likeness, s3$group == "Tumor")
pw <- s3 |> filter(patient %in% pairs3) |> select(patient, group, tumor_likeness, pmd_contrast, all_of(KZFP), Epi) |>
  pivot_wider(names_from = group, values_from = c(tumor_likeness, pmd_contrast, all_of(KZFP), Epi))
paired <- map_dfr(c("tumor_likeness", "pmd_contrast", KZFP), function(v) {
  dT <- pw[[paste0(v, "_Tumor")]] - pw[[paste0(v, "_Normal")]]
  tibble(variable = v, n_pairs = sum(!is.na(dT)), median_normal = median(pw[[paste0(v, "_Normal")]], na.rm = TRUE),
         median_tumor = median(pw[[paste0(v, "_Tumor")]], na.rm = TRUE), median_paired_delta = median(dT, na.rm = TRUE),
         frac_pairs_in_expected_direction = mean(if (v == "pmd_contrast") dT < 0 else dT > 0, na.rm = TRUE),
         wilcox_paired_p = wilcox.test(dT, exact = FALSE)$p.value)
})
gse_norm <- ext$sc_g |> filter(Group == "Normal")
field <- s3 |> filter(group == "Normal") |>
  summarise(n = n(), median_adjacent = median(tumor_likeness),
            frac_above_GSE204943_normal_max = mean(tumor_likeness > max(gse_norm$tumor_likeness)),
            frac_above_0.35 = mean(tumor_likeness > 0.35))
cat(sprintf("\nTumour signature transfer: AUC tumour vs adjacent = %.3f\n", auc3))
print(paired |> mutate(across(where(is.numeric), ~ signif(.x, 3))))
cat("Tumour-adjacent 'normal' tissue vs GSE204943 healthy normals:\n"); print(field |> mutate(across(where(is.numeric), ~ round(.x, 3))))
r_epi_n <- with(s3 |> filter(group == "Normal"), cor(tumor_likeness, Epi, method = "spearman"))
cat(sprintf("Within adjacent tissue: Spearman(score, Epi) = %.2f\n", r_epi_n))
write_csv(s3 |> mutate(across(where(is.numeric), ~ round(.x, 4))), TAB("gse234379_scores.csv"))
write_csv(paired |> mutate(across(where(is.numeric), ~ signif(.x, 3))), TAB("gse234379_paired_tests.csv"))

## ================================================================== B
pd <- readRDS(file.path(GSE_DIR, "crucial_objects/pd.rds"))
Bg <- readRDS(file.path(GSE_DIR, "crucial_objects/beta.rds")); colnames(Bg) <- pd$Basename
st <- readRDS(file.path(EXT_DIR, "story_analyses.rds"))
kzg <- prom_cgi(KZFP, rownames(Bg))
gk <- tibble(Sample = colnames(Bg), Group = as.character(c(`1_Tumor` = "Tumor", `2_OPL` = "OPL", `3_Normal` = "Normal")[pd$Sample_Group])) |>
  left_join(st$opl |> select(Sample, methylome_state), by = "Sample") |>
  mutate(state = case_when(Group == "OPL" ~ paste("OPL", methylome_state), TRUE ~ Group))
for (g in KZFP) gk[[g]] <- colMeans(Bg[kzg[[g]], gk$Sample, drop = FALSE])
kz_state <- gk |> group_by(state) |> summarise(n = n(), across(all_of(KZFP), ~ round(median(.x), 3)), .groups = "drop")
kz_state_p <- map_dfr(KZFP, function(g) { d <- gk |> filter(Group == "OPL")
  tibble(gene = g, n_cpg = length(kzg[[g]]), wilcox_p_OPL_states = wilcox.test(d[[g]] ~ d$methylome_state, exact = FALSE)$p.value) })
cat("\nKZFP promoter methylation by group / OPL methylome state (GSE204943):\n"); print(kz_state); print(kz_state_p |> mutate(wilcox_p_OPL_states = signif(wilcox_p_OPL_states, 2)))

## pilot pairs (EPIC v2 collapsed IDs share EPIC v1 cg names)
p <- readRDS(file.path(OUT_DIR, "processed.rds"))
nm <- setNames(p$sample_sheet$Sample_Name, paste0(p$sample_sheet$Patient, p$sample_sheet$Group))
kzp <- prom_cgi(KZFP, rownames(p$betas))
pil <- map_dfr(KZFP, function(g) {
  b <- colMeans(p$betas[kzp[[g]], , drop = FALSE])
  tibble(gene = g, n_cpg = length(kzp[[g]]),
         d_P2_tumor = b[nm["Patient 2Tumor"]] - b[nm["Patient 2Normal"]],
         d_P3 = b[nm["Patient 3Lesion"]] - b[nm["Patient 3Normal"]],
         d_P4 = b[nm["Patient 4Lesion"]] - b[nm["Patient 4Normal"]],
         d_P5_progressor = b[nm["Patient 5Lesion"]] - b[nm["Patient 5Normal"]])
})
cat("\nPilot KZFP promoter delta-beta (lesion minus own adjacent):\n"); print(pil |> mutate(across(where(is.double), ~ round(.x, 3))))
write_csv(kz_state, TAB("kzfp_by_opl_state_gse204943.csv"))
write_csv(pil |> mutate(across(where(is.double), ~ round(.x, 3))), TAB("kzfp_pilot_deltas.csv"))

saveRDS(list(s3 = s3, auc3 = auc3, paired = paired, field = field, r_epi_n = r_epi_n, n_pairs = length(pairs3),
             gk = gk, kz_state = kz_state, kz_state_p = kz_state_p, pil = pil),
        file.path(EXT_DIR, "gse234379_kzfp.rds"))
cat("done\n")
