## Integrity of the GSE204943 objects built in stage 01: are the inputs the deposited files,
## are the sample labels right, and do any doubtful samples affect the results?
##   A. files      local IDATs = the IDATs in GEO's file list (and the RAW archive size, if present)
##   B. labels     for every sample, name, group and sex from the sample sheet = GEO metadata,
##                 joined through the IDAT barcode (GEO series matrix)
##   C. quality    detection-P failures and median intensity per array
##   D. sex        recorded sex vs the data: Y-chromosome intensity (two-cluster split of
##                 median log2 Y minus X intensity) and X CpG-island methylation (X inactivation
##                 gives intermediate island methylation, mean beta > 0.3, in females). minfi::getSex is not
##                 used: its fixed cutoff mis-calls many females in this noob-normalised set.
##   E. identity   genotype concordance on the 59 SNP probes; pairs > 0.85 would indicate the same
##                 individual twice
##   F. profile    each sample's correlation with the median profile of each group (labels that
##                 look like another group; tumour-like OPLs are expected to resemble OSCC)
##   G. robustness the tumour signature re-derived without sex-discordant OSCCs and without OSCCs
##                 scoring in the healthy range: OPL scores, OPL groups and pilot shifts
##                 (same limma model as 03_methylation_cohorts/01, with group means re-estimated
##                 without the dropped samples)
##
## Inputs : <gse>/crucial_objects, <gse>/intermediate (stage 01), GEO series matrix and file list
##          next to the GSE204943 sample sheet (scripts/00_download), signature and pilot betas.
## Outputs: <ext>/tables/integrity_gse204943_samples.csv (one row per sample)
##          <ext>/tables/integrity_gse204943_summary.csv  (check, result, status: pass / note / FAIL)
## Run from the repository root:  Rscript scripts/06_audit/04_gse204943_integrity.R  (or make audit; ~10 min)

source("R/lib/config.R"); source("R/lib/utils.R")
suppressPackageStartupMessages({ library(minfi); library(limma); library(tidyverse) })
G <- GSE_DIR; GEOM <- dirname(GSE204943_IDAT)
SM <- file.path(GEOM, "GSE204943_series_matrix.txt.gz"); FL <- file.path(GEOM, "GSE204943_filelist.txt")
require_files(SM, FL, file.path(G, "crucial_objects", c("pd.rds", "beta.rds")),
              file.path(G, "intermediate", c("RGset_qc.rds", "Mset_noob.rds")), file.path(EXT_DIR, "gse_signature.rds"))
pd <- readRDS(file.path(G, "crucial_objects", "pd.rds"))
## status: pass / FAIL, "note" for findings that need a look but not a fix, "info" for descriptive rows (ok = NA)
chk <- list(); check <- function(name, result, ok, note_only = FALSE)
  chk[[length(chk) + 1]] <<- tibble(check = name, result = result,
                                    status = if (is.na(ok)) "info" else if (isTRUE(ok)) "pass" else if (note_only) "note" else "FAIL")

## ------------------------------------------------------------------ A. files
fl <- read_tsv(FL, show_col_types = FALSE, col_types = cols(.default = "c"))
geo_idat <- fl$Name[grepl("\\.idat\\.gz$", fl$Name)] |> sub(pattern = "^GSM[0-9]+_", replacement = "") |> sub(pattern = "\\.gz$", replacement = "")
loc_idat <- list.files(GSE204943_IDAT, pattern = "_(Grn|Red)\\.idat(\\.gz)?$") |> sub(pattern = "^GSM[0-9]+_", replacement = "") |> sub(pattern = "\\.gz$", replacement = "")
check("Local IDATs are the IDATs deposited on GEO",
      sprintf("%d local, %d on GEO; %d missing, %d not on GEO", length(loc_idat), length(geo_idat), length(setdiff(geo_idat, loc_idat)), length(setdiff(loc_idat, geo_idat))),
      setequal(geo_idat, loc_idat))
tar <- file.path(GEOM, "GSE204943_RAW.tar")
if (file.exists(tar)) check("RAW archive size equals GEO's", sprintf("%s bytes (GEO %s)", format(file.size(tar), big.mark = ","), format(as.numeric(fl$Size[fl$Name == "GSE204943_RAW.tar"]), big.mark = ",")),
                            file.size(tar) == as.numeric(fl$Size[fl$Name == "GSE204943_RAW.tar"]))

## ------------------------------------------------------------------ B. labels vs GEO metadata
sm <- readLines(gzfile(SM))
field <- function(key) lapply(sm[startsWith(sm, key)], function(z) gsub('"', "", strsplit(z, "\t")[[1]][-1]))
ch <- field("!Sample_characteristics_ch1")
char <- function(prefix) { for (v in ch) if (any(startsWith(v, prefix))) return(sub(paste0("^", prefix, ":\\s*"), "", v)); NA }
geo <- tibble(gsm = field("!Sample_geo_accession")[[1]], geo_name = field("!Sample_title")[[1]],
              geo_sex = char("Sex"), geo_group = char("sample_group"),
              Basename = sub("^GSM[0-9]+_(.*)_Grn\\.idat\\.gz$", "\\1", basename(field("!Sample_supplementary_file")[[1]])))
S <- pd |> transmute(Basename, name = gsub(" ", "", Sample_Name), group = Sample_Group, sex = Sex) |> left_join(geo, by = "Basename")
check("Every sample is a GEO sample (matched by IDAT barcode)", sprintf("%d of %d (GEO lists %d)", sum(!is.na(S$gsm)), nrow(S), nrow(geo)), all(!is.na(S$gsm)) && nrow(geo) == nrow(S))
check("Sample name, group and sex agree with GEO",
      sprintf("name %d, group %d, sex %d of %d", sum(S$name == S$geo_name, na.rm = TRUE), sum(S$group == S$geo_group, na.rm = TRUE), sum(S$sex == S$geo_sex, na.rm = TRUE), nrow(S)),
      isTRUE(all(S$name == S$geo_name & S$group == S$geo_group & S$sex == S$geo_sex)))

## ------------------------------------------------------------------ C-D. quality and sex (noob MethylSet)
M <- readRDS(file.path(G, "intermediate", "Mset_noob.rds")); stopifnot(identical(sampleNames(M), pd$Basename))
q <- getQC(M); S$median_log2_intensity <- (q$mMed + q$uMed) / 2
gm <- mapToGenome(M); rm(M); invisible(gc())
an <- getAnnotation(gm); cn <- getCN(gm); bx <- getBeta(gm)
S$y_minus_x <- matrixStats::colMedians(cn[an$chr == "chrY", ], na.rm = TRUE) - matrixStats::colMedians(cn[an$chr == "chrX", ], na.rm = TRUE)
S$x_island_beta <- colMeans(bx[an$chr == "chrX" & an$Relation_to_Island == "Island", ], na.rm = TRUE)
rm(gm, cn, bx); invisible(gc())
km <- kmeans(S$y_minus_x, centers = range(S$y_minus_x))
S <- S |> mutate(sex_from_Y = ifelse(km$cluster == which.max(km$centers), "Male", "Female"),
                 sex_from_X = ifelse(x_island_beta > 0.3, "Female", "Male"),
                 sex_discordant = sex_from_Y != sex & sex_from_X != sex)
dis <- S |> filter(sex_discordant)
check("Recorded sex agrees with X/Y data (both signals)",
      sprintf("%d of %d agree%s", sum(!S$sex_discordant), nrow(S), if (nrow(dis)) paste0("; discordant: ", paste(sprintf("%s (%s, recorded %s)", dis$name, dis$group, dis$sex), collapse = ", ")) else ""),
      nrow(dis) == 0, note_only = TRUE)

## ------------------------------------------------------------------ C, E. detection failures and genotypes (RGChannelSet)
RG <- readRDS(file.path(G, "intermediate", "RGset_qc.rds")); stopifnot(identical(sampleNames(RG), pd$Basename))
S$failed_probe_fraction <- colMeans(detectionP(RG) > 0.01)
check("Detection failures (P > 0.01) below 1% in every array", sprintf("max %.2f%% (%s)", 100 * max(S$failed_probe_fraction), S$name[which.max(S$failed_probe_fraction)]),
      max(S$failed_probe_fraction) < 0.01)
check("Median log2 intensity", sprintf("%.2f-%.2f; %d arrays below 10.5 (%s)", min(S$median_log2_intensity), max(S$median_log2_intensity),
      sum(S$median_log2_intensity < 10.5), paste(S$name[S$median_log2_intensity < 10.5], collapse = ", ")), all(S$median_log2_intensity >= 10.5), note_only = TRUE)
sb <- getSnpBeta(RG); rm(RG); invisible(gc())
gt <- ifelse(sb < 0.2, 0, ifelse(sb > 0.8, 2, 1))
conc <- outer(seq_len(ncol(gt)), seq_len(ncol(gt)), Vectorize(function(i, j) mean(gt[, i] == gt[, j]))); diag(conc) <- NA
S$max_genotype_concordance <- apply(conc, 2, max, na.rm = TRUE); S$most_similar_sample <- S$name[apply(conc, 2, which.max)]
dup <- which(conc > 0.85 & upper.tri(conc), arr.ind = TRUE)
check("No individual appears twice (SNP-probe concordance > 0.85)",
      sprintf("%d SNP probes; max pairwise concordance %.2f; %d pairs above 0.85", nrow(sb), max(conc, na.rm = TRUE), nrow(dup)), nrow(dup) == 0)

## ------------------------------------------------------------------ F. profile vs group
Bg <- readRDS(file.path(G, "crucial_objects", "beta.rds")); colnames(Bg) <- pd$Basename
set.seed(1); sub_cpg <- sample(nrow(Bg), 50000)
gmed <- sapply(split(seq_len(ncol(Bg)), pd$Sample_Group), function(i) matrixStats::rowMedians(Bg[sub_cpg, i]))
rr <- cor(Bg[sub_cpg, ], gmed)
S$r_own_group <- rr[cbind(seq_len(nrow(rr)), match(S$group, colnames(gmed)))]; S$nearest_group_profile <- colnames(gmed)[max.col(rr)]
S <- S |> left_join(read_csv(TAB("tumor_likeness_gse204943.csv"), show_col_types = FALSE) |> select(Basename = Sample, tumor_likeness, epithelial_fraction = Epi), by = "Basename")
cnv_f <- file.path(G, "results", "CNV_segments_all.csv")
if (file.exists(cnv_f)) S <- S |> left_join(read_csv(cnv_f, show_col_types = FALSE) |> mutate(len = loc.end - loc.start) |> group_by(Basename = Sample) |>
                                              summarise(fraction_genome_altered = sum(len[abs(seg.mean) >= 0.1]) / sum(len), .groups = "drop"), by = "Basename")
off <- S |> filter(nearest_group_profile != group)
check("Samples whose profile is nearest another group's median", sprintf("%d: %s", nrow(off), paste(sprintf("%s %s->%s", off$name, sub("^[0-9]_", "", off$group), sub("^[0-9]_", "", off$nearest_group_profile)), collapse = ", ")),
      NA)   # informational: tumour-like OPLs resemble OSCC by design; normal-like OSCCs are listed for review

## ------------------------------------------------------------------ G. robustness of the signature to doubtful OSCCs
p <- readRDS(file.path(OUT_DIR, "processed.rds")); ss <- p$sample_sheet
ss$Short <- paste0("P", sub("Patient ", "", ss$Patient), "-", c(Normal = "Adj", Lesion = "Lesion", Tumor = "Tumor")[as.character(ss$Group)])
ss$Short[ss$Short == "P2-Adj"] <- "P2-NAT"
ORD <- c("P2-NAT", "P2-Tumor", "P3-Adj", "P3-Lesion", "P4-Adj", "P4-Lesion", "P5-Adj", "P5-Lesion")
Bp <- p$betas; colnames(Bp) <- ss$Short[match(colnames(Bp), ss$Sample_Name)]; Bp <- Bp[, ORD]; rm(p); invisible(gc())
common <- intersect(rownames(Bp), rownames(Bg)); Bp <- Bp[common, ]; Bg <- Bg[common, ]
cp <- read_csv(file.path(G, "results", "cell_type_proportions.csv"), show_col_types = FALSE)
grp <- factor(c(`1_Tumor` = "Tumor", `2_OPL` = "OPL", `3_Normal` = "Normal")[pd$Sample_Group], levels = c("Normal", "OPL", "Tumor"))
epi <- cp$Epithelial[match(pd$Basename, cp$Sample)]; ic <- cp$Immune[match(pd$Basename, cp$Sample)]
plate <- ifelse(pd$Sample_Plate == "", "none", pd$Sample_Plate)
sig0 <- readRDS(file.path(EXT_DIR, "gse_signature.rds"))
pick <- function(s) list(hyper = s |> filter(fdr_TvN < 1e-6, d_TvN > 0.2) |> slice_max(t_TvN, n = 2000) |> pull(cpg),
                         hypo  = s |> filter(fdr_TvN < 1e-6, d_TvN < -0.2) |> slice_min(t_TvN, n = 2000) |> pull(cpg))
score <- function(B, Sg, s) { sc <- function(cpg) { m <- s[match(cpg, s$cpg), ]; apply((B[cpg, , drop = FALSE] - m$muN) / (m$muT - m$muN), 2, median, na.rm = TRUE) }
  (sc(Sg$hyper) + sc(Sg$hypo)) / 2 }
shift <- function(x) c(P2 = x[["P2-Tumor"]] - x[["P2-NAT"]], P3 = x[["P3-Lesion"]] - x[["P3-Adj"]], P4 = x[["P4-Lesion"]] - x[["P4-Adj"]], P5 = x[["P5-Lesion"]] - x[["P5-Adj"]])
refit <- function(drop) {
  keep <- grp %in% c("Normal", "Tumor") & !pd$Basename %in% drop; g <- droplevels(grp[keep])
  tt <- topTable(eBayes(lmFit(beta_to_m(Bg[, keep]), model.matrix(~ g + epi[keep] + ic[keep] + plate[keep]))), coef = 2, number = Inf, sort.by = "none")
  muN <- rowMeans(Bg[, keep & grp == "Normal"]); muT <- rowMeans(Bg[, keep & grp == "Tumor"])
  tibble(cpg = rownames(tt), t_TvN = tt$t, fdr_TvN = tt$adj.P.Val, muN = muN[rownames(tt)], muT = muT[rownames(tt)], d_TvN = muT - muN)
}
S0 <- pick(sig0); o <- grp == "OPL"; g0 <- score(Bg, S0, sig0); p0 <- score(Bp, S0, sig0)
drops <- list(`sex-discordant OSCCs` = S$Basename[S$sex_discordant & S$group == "1_Tumor"],
              `OSCCs scoring in the healthy range (< 0.35)` = pd$Basename[grp == "Tumor" & g0 < 0.35])
drops <- drops[lengths(drops) > 0]
for (nm in names(drops)) {
  s1 <- refit(drops[[nm]]); S1 <- pick(s1); g1 <- score(Bg, S1, s1); d <- max(abs(shift(score(Bp, S1, s1)) - shift(p0)))
  same <- sum((g0[o] > 0.35) == (g1[o] > 0.35))
  check(paste("Signature without", nm),
        sprintf("%d removed (%s); OPL score r = %.4f; OPL groups unchanged %d of %d; max pilot shift change %.3f",
                length(drops[[nm]]), paste(S$name[match(drops[[nm]], S$Basename)], collapse = ", "), cor(g0[o], g1[o]), same, sum(o), d),
        same == sum(o) && d < 0.05)
}

## ------------------------------------------------------------------ output
summary_tab <- bind_rows(chk)
write_csv(S |> mutate(across(where(is.double), ~ signif(.x, 4))), TAB("integrity_gse204943_samples.csv"))
write_csv(summary_tab, TAB("integrity_gse204943_summary.csv"))
options(width = 220); print(as.data.frame(summary_tab), right = FALSE)
if (any(summary_tab$status == "FAIL")) message("\nGSE204943 integrity: at least one check FAILED (see table above).")
