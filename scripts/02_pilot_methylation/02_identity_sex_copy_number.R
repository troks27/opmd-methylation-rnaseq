## Raw-signal analyses that need the SigDFs (not just betas):
##   1. sex inference (fixes the NA calls in the original report: inferSex()
##      in sesame >= 1.2x takes betas, not a SigDF)
##   2. genotype concordance on the 65 rs-probes -> confirms sample pairing
##   3. copy-number inference from total intensity (conumee-style), used as an
##      orthogonal check of the WES/XHMM CNA calls in manuscript Fig. 2
##
## Run from the repository root:  Rscript scripts/02_pilot_methylation/02_identity_sex_copy_number.R  (or make methylation)

source("R/lib/config.R")   # paths and helpers (TAB, FIG, OUT_DIR, GSE_DIR, ...)
suppressPackageStartupMessages({
  library(sesame); library(sesameData); library(DNAcopy)
  library(GenomicRanges); library(tidyverse)
  library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
})
source(file.path(REPO_ROOT, "R/lib/sample_sheet.R"))

require_files(PILOT_SHEET)
dir.create(file.path(EXT_DIR, "tables"),  recursive = TRUE, showWarnings = FALSE)

ss <- build_sample_sheet(PILOT_SHEET, PILOT_IDAT_DIR)
ss$Short <- paste0("P", sub("Patient ", "", ss$Patient), "-",
                   c(Normal = "Adj", Lesion = "Lesion", Tumor = "Tumor")[as.character(ss$Group)])
ss$Short[ss$Short == "P2-Adj"] <- "P2-NAT"
ORD <- c("P2-NAT", "P2-Tumor", "P3-Adj", "P3-Lesion", "P4-Adj", "P4-Lesion", "P5-Adj", "P5-Lesion")

sdfs <- lapply(setNames(ss$Basename, ss$Short), readIDATpair, platform = "EPICv2")
## same preprocessing as the main report (openSesame default for EPICv2)
sdfs_p <- lapply(sdfs, prepSesame, prep = "QCDPB")

## ---------------------------------------------------------------- 1. sex
betas_all <- sapply(sdfs_p, getBetas)
sex <- sapply(colnames(betas_all), function(s)
  tryCatch(as.character(inferSex(betas_all[, s], platform = "EPICv2")),
           error = function(e) NA_character_))

## independent check: chrY total intensity relative to autosomes
mft <- sesameData_getManifestGRanges("EPICv2")
chr_of <- setNames(as.character(seqnames(mft)), names(mft))
ti_raw <- sapply(sdfs, function(s) totalIntensities(noob(dyeBiasNL(s))))
auto_p <- names(chr_of)[chr_of %in% paste0("chr", 1:22)]
y_p    <- names(chr_of)[chr_of == "chrY"]
x_p    <- names(chr_of)[chr_of == "chrX"]
y_ratio <- log2(colMeans(ti_raw[intersect(y_p, rownames(ti_raw)), ], na.rm = TRUE) /
                colMeans(ti_raw[intersect(auto_p, rownames(ti_raw)), ], na.rm = TRUE))
## X-linked methylation: ~0.5+ in females (Xi), low in males
x_beta <- colMeans(betas_all[intersect(x_p, rownames(betas_all)), ], na.rm = TRUE)

## pOOBAH-only detection rate (the report's 'pass rate' also counted design masks)
pvals <- sapply(sdfs, function(s) pOOBAH(s, return.pval = TRUE))
pooba <- colMeans(pvals < 0.05, na.rm = TRUE)

qc2 <- tibble(Sample = ss$Short, Patient = ss$Patient, Group = as.character(ss$Group),
              Sex_sesame = sex[ss$Short],
              chrY_log2_vs_auto = round(y_ratio[ss$Short], 2),
              chrX_mean_beta = round(x_beta[ss$Short], 3),
              pOOBAH_detection_rate = round(pooba[ss$Short], 4))
write_csv(qc2, file.path(EXT_DIR, "tables", "qc_sex_detection.csv"))
print(qc2)

## ---------------------------------------------------------------- 2. genotype
rs <- grep("^rs", rownames(betas_all), value = TRUE)
rs_b <- sapply(sdfs_p, function(s) getBetas(s, mask = FALSE))[rs, ]
gt <- apply(rs_b, 2, function(b) ifelse(b < 0.3, 0L, ifelse(b > 0.7, 2L, 1L)))
conc <- outer(colnames(gt), colnames(gt), Vectorize(function(a, b)
  mean(gt[, a] == gt[, b], na.rm = TRUE)))
dimnames(conc) <- list(colnames(gt), colnames(gt))
conc <- conc[ORD, ORD]
write.csv(round(conc, 3), file.path(EXT_DIR, "tables", "genotype_concordance_rs_probes.csv"))
print(round(conc, 2))

## ---------------------------------------------------------------- 3. CNV
## log2 total intensity of autosomal CpG probes (masked probes dropped)
ti <- sapply(sdfs, function(s) totalIntensities(noob(dyeBiasNL(qualityMask(s))), mask = TRUE))
ti <- ti[intersect(rownames(ti), auto_p[grepl("^cg", auto_p)]), ]
ti <- ti[complete.cases(ti) & apply(ti, 1, min) > 500, ]
lti <- log2(ti)
lti <- sweep(lti, 2, apply(lti, 2, median))          # per-sample centering

pos <- start(mft[rownames(lti)])
chr <- chr_of[rownames(lti)]
ord <- order(match(chr, paste0("chr", 1:22)), pos)
lti <- lti[ord, ]; pos <- pos[ord]; chr <- chr[ord]

cyto <- sesameData_getGenomeInfo("hg38")$cytoBand |>
  mutate(across(c(chrom, name), as.character))
arm_bounds <- cyto |> filter(chrom %in% paste0("chr", 1:22)) |>
  group_by(chrom) |>
  summarise(cen = max(chromEnd[startsWith(name, "p")]), len = max(chromEnd))
## probes in methylation-dense clusters whose total intensity tracks
## methylation, not copy number (HOX A-D, ZIC1/4) -> excluded from CN
meth_clusters <- GRanges(c("chr7", "chr17", "chr12", "chr2", "chr3"),
                         IRanges(c(27.09e6, 48.52e6, 53.93e6, 176.08e6, 147.38e6),
                                 c(27.25e6, 48.63e6, 54.06e6, 176.23e6, 147.42e6)))
in_cluster <- overlapsAny(GRanges(chr, IRanges(pos, width = 1)), meth_clusters)
cl <- arm_bounds$len[match(chr, arm_bounds$chrom)]
terminal <- pos < 5e6 | pos > cl - 5e6               # GC-wave-prone ends
keep <- !in_cluster
lti <- lti[keep, ]; pos <- pos[keep]; chr <- chr[keep]; terminal <- terminal[keep]

adj_samples <- ss$Short[ss$Group == "Normal"]

## Mode A ("LPO"): leave-patient-out reference = median of the *other* patients'
##   adjacent tissue. Detects field events shared by a patient's adjacent and
##   lesion (e.g. P4 chr12). Median, not mean, so one aberrant reference
##   (P3-Adj 9p21, P4-Adj chr12) cannot leak into every other profile.
## Mode B ("paired"): lesion/tumor vs its own adjacent. Cancels germline CNVs
##   and patient-specific signal artefacts, but is blind to shared field events.
fit_lr <- function(q, ref) resid(lm(q ~ ref))        # conumee-style affinity fit
lr_lpo <- sapply(ss$Short, function(s) {
  refs <- setdiff(adj_samples, ss$Short[ss$Patient == ss$Patient[ss$Short == s]])
  fit_lr(lti[, s], matrixStats::rowMedians(lti[, refs, drop = FALSE]))
})
les <- ss$Short[ss$Group != "Normal"]
own_adj <- setNames(ss$Short[match(paste(ss$Patient[match(les, ss$Short)], "Normal"),
                                   paste(ss$Patient, ss$Group))], les)
lr_pair <- sapply(les, function(s) fit_lr(lti[, s], lti[, own_adj[s]]))

## 50-probe bins, median per bin
bin_id <- ave(seq_along(chr), chr, FUN = function(i) ceiling(seq_along(i) / 50))
bin_key <- paste(chr, bin_id)
make_bins <- function(lr) {
  tibble(chr = chr, pos = pos, bin = bin_key, term = terminal) |>
    bind_cols(as_tibble(lr)) |>
    group_by(chr, bin) |>
    summarise(start = min(pos), end = max(pos), n = n(), terminal = mean(term) > 0.5,
              across(all_of(colnames(lr)), median), .groups = "drop") |>
    mutate(mid = (start + end) / 2,
           arm = ifelse(mid < arm_bounds$cen[match(chr, arm_bounds$chrom)], "p", "q")) |>
    arrange(match(chr, paste0("chr", 1:22)), mid) |>
    mutate(across(all_of(colnames(lr)), ~ .x - median(.x[!terminal])))
}
segment_bins <- function(b, ids) {
  set.seed(1)   # segment() assesses change points by permutation; fix the seed so calls are reproducible
  cna <- CNA(as.matrix(b[, ids]), chrom = b$chr, maploc = b$mid,
             data.type = "logratio", sampleid = make.names(ids))
  segment(smooth.CNA(cna), alpha = 0.001, min.width = 3,
          undo.splits = "sdundo", undo.SD = 2, verbose = 0)$output |>
    as_tibble() |> mutate(ID = ids[match(ID, make.names(ids))])
}
bins_lpo  <- make_bins(lr_lpo);  segs_lpo  <- segment_bins(bins_lpo, ss$Short)
bins_pair <- make_bins(lr_pair); segs_pair <- segment_bins(bins_pair, les)
write_csv(segs_lpo,  file.path(EXT_DIR, "tables", "cnv_segments_LPO.csv"))
write_csv(segs_pair, file.path(EXT_DIR, "tables", "cnv_segments_paired.csv"))

## ---- arm-level calls (terminal 5 Mb excluded, arms with >= 10 bins).
## Rule: |arm median| >= 0.07 vs other patients' adjacent tissue (LPO).
## (A noise-scaled rule such as 3 x MAD-of-arm-medians is circular for an
## aneuploid tumour: the more arms altered, the higher its own threshold.)
## Each LPO call is then annotated:
##   paired_confirmed - same-direction >= 0.07 vs the patient's own adjacent
##   field_shared     - own adjacent tissue carries the same event (LPO)
ARM_T <- 0.07
arm_medians <- function(b, ids, mode) {
  b |> filter(!terminal) |>
    pivot_longer(all_of(ids), names_to = "Sample", values_to = "lr") |>
    group_by(Sample, chr, arm) |>
    summarise(n_bins = n(), arm_median = median(lr), .groups = "drop") |>
    filter(n_bins >= 10) |>
    mutate(mode = mode, arm_label = paste0(sub("chr", "", chr), arm))
}
arms <- bind_rows(arm_medians(bins_lpo, ss$Short, "LPO"), arm_medians(bins_pair, les, "paired"))
arm_tab <- arms |> select(Sample, arm_label, mode, arm_median) |>
  pivot_wider(names_from = mode, values_from = arm_median) |>
  mutate(Patient = ss$Patient[match(Sample, ss$Short)],
         own_adj = own_adj[Sample]) |>
  left_join(arms |> filter(mode == "LPO") |> select(own_adj = Sample, arm_label, adj_LPO = arm_median),
            by = c("own_adj", "arm_label")) |>
  mutate(call = case_when(LPO >=  ARM_T ~ "gain", LPO <= -ARM_T ~ "loss", TRUE ~ ""),
         paired_confirmed = call != "" & !is.na(paired) & sign(paired) == sign(LPO) & abs(paired) >= ARM_T,
         field_shared = call != "" & !is.na(adj_LPO) & sign(adj_LPO) == sign(LPO) & abs(adj_LPO) >= ARM_T)
write_csv(arm_tab, file.path(EXT_DIR, "tables", "cnv_arm_level.csv"))
cat("\nArm-level calls (|LPO| >= 0.07):\n")
print(arm_tab |> filter(call != "") |>
        select(Sample, arm_label, LPO, paired, call, paired_confirmed, field_shared) |>
        mutate(across(c(LPO, paired), ~ round(.x, 3))) |>
        arrange(factor(Sample, ORD)), n = 100)

## ---- regions / genes named in manuscript Fig. 2 (hg38)
band_range <- function(chrom, bands) {
  r <- cyto[cyto$chrom == chrom & cyto$name %in% bands, ]
  c(min(r$chromStart), max(r$chromEnd))
}
txns <- sesameData_getTxnGRanges("hg38")
gene_range <- function(g) {
  r <- range(txns[mcols(txns)$gene_name == g])
  c(as.character(seqnames(r))[1], start(r)[1], end(r)[1])
}
genes <- c("CDKN2A", "CDK4", "MDM2", "ERBB2", "FHIT", "CCND1", "EGFR", "MYC", "PIK3CA", "TP63", "SOX2")
gene_tab <- map_dfr(genes, function(g) {
  r <- tryCatch(gene_range(g), error = function(e) c(NA, NA, NA))
  band <- if (!is.na(r[1])) paste0(sub("chr", "", r[1]),
    cyto$name[cyto$chrom == r[1] & cyto$chromStart <= as.numeric(r[2]) & cyto$chromEnd >= as.numeric(r[2])]) else NA
  tibble(region = g, chr = r[1], start = as.numeric(r[2]), end = as.numeric(r[3]), cytoband = band)
})
write_csv(gene_tab, file.path(EXT_DIR, "tables", "gene_cytobands_hg38.csv"))
print(gene_tab)

fig2_bands <- tribble(~region, ~chr, ~bands,
  "2q31.2", "chr2", "q31.2", "2q32.1", "chr2", "q32.1", "19p13.2", "chr19", "p13.2",
  "12q13.13", "chr12", "q13.13", "12p13.2", "chr12", "p13.2", "9p21.3", "chr9", "p21.3",
  "12q13.2", "chr12", "q13.2", "12q24.31", "chr12", "q24.31", "3p21.31", "chr3", "p21.31",
  "13q12.13", "chr13", "q12.13", "19p13.3", "chr19", "p13.3", "17q12", "chr17", "q12",
  "17q21.31", "chr17", "q21.31", "1p22.2", "chr1", "p22.2", "3p14.3", "chr3", "p14.3",
  "13q33.1", "chr13", "q33.1") |>
  rowwise() |> mutate(start = band_range(chr, bands)[1], end = band_range(chr, bands)[2],
                      cytoband = region) |>
  ungroup() |> select(region, chr, start, end, cytoband)
regions <- bind_rows(fig2_bands, gene_tab) |>
  mutate(s = ifelse(end - start < 1e6, (start + end) / 2 - 5e5, start),
         e = ifelse(end - start < 1e6, (start + end) / 2 + 5e5, end))
region_summary <- function(b, ids) {
  m <- sapply(ids, function(smp) sapply(seq_len(nrow(regions)), function(k) {
    sel <- b$chr == regions$chr[k] & b$end >= regions$s[k] & b$start <= regions$e[k]
    if (!isTRUE(any(sel))) NA else round(median(b[[smp]][sel %in% TRUE]), 3)
  }))
  bind_cols(regions[, c("region", "chr", "start", "end")], as_tibble(m))
}
reg_lpo  <- region_summary(bins_lpo, ORD)
reg_pair <- region_summary(bins_pair, intersect(ORD, les))
write_csv(reg_lpo,  file.path(EXT_DIR, "tables", "cnv_fig2_regions_LPO.csv"))
write_csv(reg_pair, file.path(EXT_DIR, "tables", "cnv_fig2_regions_paired.csv"))
cat("\nFig. 2 regions, LPO bin medians:\n"); print(reg_lpo |> select(-start, -end), n = 40, width = 200)
cat("\nFig. 2 regions, paired bin medians:\n"); print(reg_pair |> select(-start, -end), n = 40, width = 200)

## ---- 9p21 deletion check: a homozygous deletion should also drop probes
## below detection (pOOBAH) and push intensity far down at probe level.
p9 <- names(chr_of)[chr_of == "chr9"]
p9 <- p9[start(mft[p9]) > 21.6e6 & start(mft[p9]) < 22.25e6]
del9 <- tibble(Sample = colnames(pvals),
               n_probes = length(intersect(p9, rownames(pvals))),
               frac_fail_9p21 = colMeans(pvals[intersect(p9, rownames(pvals)), ] >= 0.05, na.rm = TRUE),
               frac_fail_genome = colMeans(pvals >= 0.05, na.rm = TRUE),
               lpo_9p21_median = sapply(colnames(pvals), function(s)
                 median(lr_lpo[chr == "chr9" & pos > 21.6e6 & pos < 22.25e6, s]))) |>
  mutate(fail_ratio = frac_fail_9p21 / frac_fail_genome) |>
  arrange(factor(Sample, ORD))
write_csv(del9, file.path(EXT_DIR, "tables", "cnv_9p21_detection_check.csv"))
print(del9)

## burden: fraction of non-terminal genome in segments >= 3 Mb with |seg| > 0.12,
## plus focal (< 3 Mb) events with |seg| > 0.3
burden_tab <- function(segs, mode) {
  segs |> mutate(len = loc.end - loc.start,
                 cl = arm_bounds$len[match(chrom, arm_bounds$chrom)],
                 term = loc.end < 5e6 | loc.start > cl - 5e6,
                 broad = len >= 3e6 & abs(seg.mean) > 0.12 & !term,
                 focal = len < 3e6 & abs(seg.mean) > 0.3 & num.mark >= 4 & !term) |>
    group_by(Sample = ID) |>
    summarise(frac_genome_altered = sum(len[broad]) / 2.88e9,
              n_broad_events = sum(broad), n_focal_events = sum(focal), .groups = "drop") |>
    mutate(mode = mode)
}
burden <- bind_rows(burden_tab(segs_lpo, "LPO"), burden_tab(segs_pair, "paired")) |>
  arrange(mode, factor(Sample, ORD))
write_csv(burden, file.path(EXT_DIR, "tables", "cnv_burden_array.csv"))
print(burden)

## ---- genome-wide plot (terminal bins greyed)
if (!requireNamespace("ggnewscale", quietly = TRUE))
  install.packages("ggnewscale", repos = "https://cloud.r-project.org")
plot_cn <- function(b, segs, ids, file, title) {
  chr_len <- arm_bounds |> arrange(match(chrom, paste0("chr", 1:22))) |>
    mutate(offset = lag(cumsum(as.numeric(len)), default = 0))
  gb <- b |> left_join(chr_len, by = c("chr" = "chrom")) |> mutate(g = mid + offset) |>
    pivot_longer(all_of(ids), names_to = "Sample", values_to = "lr") |>
    mutate(Sample = factor(Sample, ORD))
  gs <- segs |> left_join(chr_len, by = c("chrom" = "chrom")) |>
    mutate(gs = loc.start + offset, ge = loc.end + offset, Sample = factor(ID, ORD),
           state = ifelse(seg.mean > 0.12, "gain", ifelse(seg.mean < -0.12, "loss", "neutral")))
  p <- ggplot() +
    geom_point(data = gb, aes(g, pmax(pmin(lr, 0.8), -0.8), colour = terminal),
               size = 0.15, alpha = 0.5, show.legend = FALSE) +
    scale_colour_manual(values = c(`FALSE` = "grey55", `TRUE` = "grey85")) +
    ggnewscale::new_scale_colour() +
    geom_segment(data = gs, aes(x = gs, xend = ge, y = seg.mean, yend = seg.mean, colour = state),
                 linewidth = 0.9) +
    geom_vline(xintercept = chr_len$offset, colour = "grey85", linewidth = 0.3) +
    scale_colour_manual(values = c(gain = "#c0392b", loss = "#2e6db4", neutral = "grey20"), name = NULL) +
    scale_x_continuous(breaks = chr_len$offset + chr_len$len / 2, labels = sub("chr", "", chr_len$chrom),
                       expand = c(0.005, 0)) +
    coord_cartesian(ylim = c(-0.8, 0.8)) +
    facet_wrap(~Sample, ncol = 1, strip.position = "right") +
    labs(x = "Chromosome", y = "log2 ratio", title = title) +
    theme_bw(base_size = 9) + theme(panel.grid = element_blank(), legend.position = "top",
                                    axis.text.x = element_text(size = 6))
  ggsave(file, p, width = 11, height = 1.25 * length(ids) + 1, dpi = 200)
}
plot_cn(bins_lpo, segs_lpo, ss$Short, file.path(EXT_DIR, "figures", "cnv_genome_LPO.png"),
        "Array copy number vs median of other patients' adjacent tissue (terminal 5 Mb in light grey)")
plot_cn(bins_pair, segs_pair, les, file.path(EXT_DIR, "figures", "cnv_genome_paired.png"),
        "Array copy number, lesion/tumor vs own adjacent tissue")

## 9p21 zoom
z <- bins_lpo |> filter(chr == "chr9", mid > 19e6, mid < 25e6) |>
  pivot_longer(all_of(ss$Short), names_to = "Sample", values_to = "lr") |>
  mutate(Sample = factor(Sample, ORD))
gz <- ggplot(z, aes(mid / 1e6, lr)) +
  annotate("rect", xmin = 21.967, xmax = 21.995, ymin = -Inf, ymax = Inf, fill = "#f4d03f", alpha = 0.5) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_point(size = 1) + geom_line(linewidth = 0.3) +
  facet_wrap(~Sample, ncol = 2) +
  labs(x = "chr9 position (Mb, hg38); yellow = CDKN2A", y = "log2 ratio (50-probe bins)",
       title = "9p21.3 focal copy number") + theme_bw(base_size = 9)
ggsave(file.path(EXT_DIR, "figures", "cnv_9p21_zoom.png"), gz, width = 7, height = 6, dpi = 200)

saveRDS(list(bins_lpo = bins_lpo, segs_lpo = segs_lpo, bins_pair = bins_pair,
             segs_pair = segs_pair, arms = arm_tab, ss = ss, qc2 = qc2, conc = conc,
             del9 = del9), file.path(EXT_DIR, "cnv_objects.rds"))
cat("done\n")
