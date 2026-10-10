## Publication figures for the RNA-seq sections (Nature Communications format), continuing
## the methylation figures (scripts/07_figures/01_figures_1_to_4.R):
##   Fig. 5  RNA of the same specimens validates the methylation-array readouts
##           (composition, copy number, 9p21) - cross-assay concordance
##   Fig. 6  The methylation-defined state in the transcriptome: proxy, programme shifts,
##           promoter methylation -> expression, and outcome specificity
##   Supplementary Fig. 3  RNA-seq QC, strandedness, sex, paired DE, muscle content, RNA-only dosage
## Same format as Figs 1-4: 183 mm wide, <= 170 mm tall, Helvetica 5-7 pt via Quartz PDF,
## 8 pt bold lowercase panel letters, no in-image titles, a key for every encoding.
## Run from the repository root:  Rscript scripts/07_figures/02_figures_5_6.R  (or make figures)
source("R/lib/config.R")   # paths and helpers (TAB, FIG, OUT_DIR, GSE_DIR, ...)
suppressPackageStartupMessages({ library(tidyverse); library(patchwork); library(ggrepel); library(scales); library(ggtext) })
select <- dplyr::select; filter <- dplyr::filter
ML <- scales::label_number(style_negative = "minus")
MLp <- scales::label_percent(accuracy = 1, style_negative = "minus")
OUT <- NATURE_DIR
r <- readRDS(file.path(EXT, "rnaseq_story.rds")); g26 <- readRDS(file.path(EXT, "gse26549_progression.rds"))

## ------------------------------------------------------------------ style (as Figs 1-4)
PT <- function(x) x / ggplot2::.pt
INK <- "#1a1a1a"; INK2 <- "#4d4c48"; MUTED <- "#8f8d86"; HAIR <- "#d9d7cf"
C_N <- "#2a78d6"; C_O <- "#eb6834"; C_T <- "#1baf7a"; BLUE <- "#2a78d6"; RED <- "#e34948"; MID <- "#f0efec"
theme_nat <- function() {
  theme_classic(base_size = 6, base_family = "Helvetica") +
    theme(line = element_line(linewidth = PT(0.5), colour = INK), axis.line = element_line(linewidth = PT(0.5), colour = INK),
          axis.ticks = element_line(linewidth = PT(0.5), colour = INK), axis.ticks.length = unit(0.9, "mm"),
          axis.text = element_text(size = 5.5, colour = INK), axis.title = element_text(size = 6, colour = INK),
          legend.text = element_text(size = 5.2, colour = INK), legend.title = element_text(size = 5.5, colour = INK),
          legend.key.size = unit(2.4, "mm"), legend.key.spacing.y = unit(0.2, "mm"), legend.spacing = unit(0.5, "mm"),
          legend.margin = margin(0, 0, 0, 0), legend.background = element_blank(), legend.box.spacing = unit(1, "mm"),
          strip.background = element_blank(), strip.text = element_text(size = 5.8, colour = INK),
          plot.margin = margin(1.2, 1.5, 1.2, 1.5, "mm"), panel.background = element_blank(), plot.background = element_blank())
}
tag_theme <- theme(plot.tag = element_text(size = 8, face = "bold", family = "Helvetica"))
## Quartz PDF embeds Helvetica (macOS); elsewhere fall back to cairo_pdf.
qpdf <- function(filename, width, height, ...) {
  if (capabilities("aqua")) grDevices::quartz(type = "pdf", file = filename, width = width, height = height)
  else grDevices::cairo_pdf(filename = filename, width = width, height = height)
}
save_nat <- function(p, name, h_mm, w_mm = 183) {
  stopifnot(w_mm <= 183, h_mm <= 170)
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w_mm / 25.4, height = h_mm / 25.4, device = qpdf)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w_mm / 25.4, height = h_mm / 25.4, dpi = 600, device = ragg::agg_png, bg = "white")
}
sup <- c("-" = "⁻", "0" = "⁰", "1" = "¹", "2" = "²", "3" = "³", "4" = "⁴", "5" = "⁵",
         "6" = "⁶", "7" = "⁷", "8" = "⁸", "9" = "⁹")
ptxt <- function(pv) vapply(pv, function(x) { if (is.na(x)) return("")
  if (x < 0.001) { e <- floor(log10(x)); m <- x / 10^e; sprintf("%.1f × 10%s", m, paste(sup[strsplit(as.character(e), "")[[1]]], collapse = "")) }
  else if (x >= 0.1) sprintf("%.2f", floor(x * 100) / 100) else formatC(signif(x, 2), format = "fg", digits = 2) }, "")
SH_T <- c(Adjacent = 21, Lesion = 24, Tumour = 22)
FILL_T <- c(Adjacent = "white", Lesion = C_O, Tumour = C_T)
tissue_of <- function(s) ifelse(grepl("Adj|NAT", s), "Adjacent", ifelse(grepl("Tumor", s), "Tumour", "Lesion"))
## column labels: patient + diagnosis abbreviation (keys in legends)
LAB <- c(`P2-NAT` = "P2 Adj", `P2-Tumor` = "P2 OSCC", `P3-Lesion` = "P3 PVL", `P4-Adj` = "P4 Adj", `P4-Lesion` = "P4 AH",
         `P5-Adj` = "P5 Adj", `P5-Lesion` = "P5 MSD†", `P6-Adj` = "P6 Adj", `P6-Lesion` = "P6 OSCC", `P7-Adj` = "P7 Adj", `P7-Lesion` = "P7 Les")
ORD11 <- names(LAB)
PCOL <- c(P2 = "#1baf7a", P3 = "#8f8d86", P4 = "#eda100", P5 = "#e34948", P6 = "#2a78d6", P7 = "#7d3cb5")
DIV <- scale_fill_gradient2(low = BLUE, mid = MID, high = RED, midpoint = 0, limits = c(-2.2, 2.2), oob = squish, name = "Score (z)", labels = ML)
cor_lab <- function(x, y, m = "spearman", sym = "ρ") sub("-", "−", sprintf("%s = %.2f", sym, cor(x, y, method = m)), fixed = TRUE)

## ================================================================== Fig. 5
cb <- read_csv(TAB("rnaseq_composition_vs_epidish.csv"), show_col_types = FALSE) |> filter(!is.na(Epi_meth)) |>
  mutate(Tissue = tissue_of(sample), lab = LAB[sample], nonepi = -(Immune + Fibroblast))
## labels are repelled from the points and from the fitted line (unlabelled obstacles along it)
line_obst <- function(d, x, y, n = 60) {
  f <- lm(reformulate(x, y), data = d); g <- tibble(!!x := seq(min(d[[x]]), max(d[[x]]), length.out = n))
  g[[y]] <- predict(f, g); g$lab <- ""; g }
f5a <- ggplot(cb, aes(IC_meth, Immune)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = MUTED, linewidth = PT(0.5)) +
  geom_point(aes(shape = Tissue, fill = Tissue), size = 1.6, stroke = PT(0.6)) +
  geom_text_repel(data = bind_rows(cb, line_obst(cb, "IC_meth", "Immune")), aes(label = lab), size = PT(5), family = "Helvetica",
                  min.segment.length = 0, segment.size = PT(0.3), seed = 1, box.padding = 0.45, force = 5, max.iter = 2e4, max.overlaps = Inf) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3, label = cor_lab(cb$IC_meth, cb$Immune), size = PT(5.5), family = "Helvetica") +
  scale_shape_manual(values = SH_T) + scale_fill_manual(values = FILL_T) +
  scale_x_continuous(labels = MLp, expand = expansion(mult = c(0.05, 0.08))) +
  scale_y_continuous(labels = ML, expand = expansion(mult = c(0.2, 0.06))) +
  labs(x = "Immune fraction (EpiDISH)", y = "Immune marker score (RNA, z)", shape = NULL, fill = NULL) +
  theme_nat() + theme(legend.position = "bottom")
f5b <- ggplot(cb, aes(Epi_meth, nonepi)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = MUTED, linewidth = PT(0.5)) +
  geom_point(aes(shape = Tissue, fill = Tissue), size = 1.6, stroke = PT(0.6), show.legend = FALSE) +
  geom_text_repel(data = bind_rows(cb, line_obst(cb, "Epi_meth", "nonepi")), aes(label = lab), size = PT(5), family = "Helvetica",
                  min.segment.length = 0, segment.size = PT(0.3), seed = 2, box.padding = 0.35, force = 4, max.iter = 2e4, max.overlaps = Inf) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3, label = cor_lab(cb$Epi_meth, cb$nonepi), size = PT(5.5), family = "Helvetica") +
  scale_shape_manual(values = SH_T) + scale_fill_manual(values = FILL_T) + scale_x_continuous(labels = MLp) + scale_y_continuous(labels = ML) +
  labs(x = "Epithelial fraction (EpiDISH)", y = "Non-epithelial score (RNA, −z)") + theme_nat()
## c. cell-type programmes across all 11 libraries
musc_z <- as.numeric(scale(r$muscle_score[rownames(r$score)]))
hm <- as_tibble(r$score[, c("T_cell", "Cytotoxic", "B_cell", "Plasma_cell", "Myeloid", "IFNG_response", "Checkpoint", "CAF", "Endothelial")], rownames = "sample") |>
  mutate(Skeletal_muscle = musc_z) |> pivot_longer(-sample, names_to = "prog", values_to = "z") |>
  mutate(sample = factor(LAB[sample], levels = LAB[ORD11]),
         prog = factor(prog, levels = rev(c("T_cell", "Cytotoxic", "B_cell", "Plasma_cell", "Myeloid", "IFNG_response", "Checkpoint", "CAF", "Endothelial", "Skeletal_muscle")),
                       labels = rev(c("T cell", "Cytotoxic", "B cell", "Plasma cell", "Myeloid", "IFN-γ response", "Immune checkpoint", "Fibroblast (CAF)", "Endothelial", "Skeletal muscle"))))
f5c <- ggplot(hm, aes(sample, prog, fill = z)) + geom_tile(colour = "white", linewidth = PT(0.4)) + DIV +
  annotate("rect", xmin = 2.5, xmax = 3.5, ymin = 0.5, ymax = 10.5, fill = NA, colour = INK, linewidth = PT(0.7)) +
  scale_x_discrete(position = "top") + labs(x = NULL, y = NULL) + theme_nat() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0),
        legend.key.height = unit(5, "mm"), legend.key.width = unit(1.6, "mm"))
## d. 9p21 expression
c9 <- read_csv(TAB("rnaseq_9p21.csv"), show_col_types = FALSE) |> filter(gene %in% c("MTAP", "CDKN2A", "CDKN2B")) |>
  mutate(evidence = case_when(sample == "P5-Lesion" ~ "Array and exome deletion", sample == "P4-Lesion" ~ "Array deletion only",
                              sample == "P2-Tumor" ~ "No array deletion", TRUE ~ "No methylation array"),
         evidence = factor(evidence, levels = c("Array and exome deletion", "Array deletion only", "No array deletion", "No methylation array")),
         gene = factor(gene, levels = c("MTAP", "CDKN2A", "CDKN2B")), lab = LAB[sample])
EV <- c("Array and exome deletion" = RED, "Array deletion only" = C_O, "No array deletion" = C_T, "No methylation array" = MUTED)
f5d <- ggplot(c9, aes(gene, log2FC_vs_own_adj, group = sample, colour = evidence)) +
  geom_hline(yintercept = 0, colour = HAIR, linewidth = PT(0.5)) +
  geom_line(linewidth = PT(0.6)) + geom_point(size = 1.4) +
  geom_text_repel(data = c9 |> filter(gene == "CDKN2B"), aes(label = lab), nudge_x = 0.45, direction = "y", hjust = 0, size = PT(5),
                  family = "Helvetica", segment.size = PT(0.3), min.segment.length = 0, show.legend = FALSE, seed = 6) +
  scale_colour_manual(values = EV) + scale_x_discrete(expand = expansion(add = c(0.3, 1.1))) + scale_y_continuous(labels = ML) +
  labs(x = NULL, y = "Expression change vs own adjacent\n(log₂ fold change)", colour = "9p21 call") +
  theme_nat() + theme(axis.text.x = element_text(face = "italic"), legend.position = "bottom") + guides(colour = guide_legend(ncol = 2, title.position = "top"))
## e. arm dosage: RNA vs array, both against other patients' adjacent tissue
lpo <- read_csv(TAB("rnaseq_arm_dosage_vs_cn_lpo.csv"), show_col_types = FALSE) |> filter(!is.na(cn_lpo))
show_s <- c("P2-Tumor", "P4-Adj", "P4-Lesion", "P5-Lesion")
key_arms <- list(`P2-Tumor` = c("20p", "20q", "4q", "3p"), `P4-Adj` = c("12p", "12q"), `P4-Lesion` = c("12p", "12q"), `P5-Lesion` = c("3p", "13q"))
le <- lpo |> filter(sample %in% show_s) |> mutate(key = map2_lgl(sample, arm_label, ~ .y %in% key_arms[[.x]]),
                                                    panel = factor(LAB[sample], levels = LAB[show_s]))
rlab <- le |> group_by(panel) |> summarise(r = sprintf("r = %.2f", cor(cn_lpo, rna_lpo)), .groups = "drop")
f5e <- ggplot(le, aes(cn_lpo, rna_lpo)) +
  geom_hline(yintercept = 0, colour = HAIR, linewidth = PT(0.4)) + geom_vline(xintercept = 0, colour = HAIR, linewidth = PT(0.4)) +
  geom_point(aes(colour = key), size = 0.9) +
  geom_text_repel(data = le |> filter(key), aes(label = arm_label), size = PT(5), family = "Helvetica", min.segment.length = 0,
                  segment.size = PT(0.3), seed = 3, colour = RED) +
  geom_text(data = rlab, aes(x = -Inf, y = Inf, label = r), hjust = -0.1, vjust = 1.3, size = PT(5.5), family = "Helvetica") +
  facet_wrap(~panel, nrow = 1) + scale_colour_manual(values = c(`FALSE` = MUTED, `TRUE` = RED), labels = c("Other arms", "Arms called by array"), name = NULL) +
  scale_x_continuous(breaks = c(-0.2, 0, 0.2), labels = ML) + scale_y_continuous(labels = ML) +
  labs(x = "Arm copy number, methylation array (log₂ ratio)", y = "Arm expression dosage (RNA)") +
  theme_nat() + theme(legend.position = "bottom", panel.spacing.x = unit(3, "mm"))
## f. cross-assay concordance of copy-number events (exome calls from the multimodal draft, Fig. 2)
dl <- function(s, a) { v <- lpo$rna_lpo[lpo$sample == s & lpo$arm_label %in% a]; mean(v) }
rna_fl <- r$floor99
## array and exome calls are curated inputs (resources/cnv_events_curated.csv); an empty rna
## cell is called here from arm dosage against the RNA noise floor
rna_call <- function(s, arms, dir) {
  if (is.na(arms)) return(NA_character_)
  v <- dl(s, strsplit(arms, ";")[[1]])
  if ((dir == "loss" && v < -rna_fl) || (dir == "gain" && v > rna_fl)) "Detected" else "Not detected"
}
ev <- read_csv(file.path(REPO_ROOT, "resources", "cnv_events_curated.csv"), show_col_types = FALSE) |>
  mutate(rna = coalesce(rna, pmap_chr(list(sample, arms, direction), rna_call))) |> select(event, array, exome, rna)
evl <- ev |> pivot_longer(-event, names_to = "assay", values_to = "call") |>
  mutate(assay = factor(assay, levels = c("array", "exome", "rna"), labels = c("Methylation array", "Exome sequencing", "RNA sequencing")),
         event = factor(event, levels = rev(ev$event)), call = factor(call, levels = c("Detected", "Not detected", "Not assessed")))
f5f <- ggplot(evl, aes(assay, event, fill = call)) + geom_tile(colour = "white", linewidth = PT(0.8)) +
  scale_fill_manual(values = c(Detected = "#1c5cab", `Not detected` = "#e6effb", `Not assessed` = "#d9d7cf"), name = NULL, drop = FALSE) +
  scale_x_discrete(position = "top") + labs(x = NULL, y = NULL) + theme_nat() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), legend.position = "right")
fig5 <- (f5a + f5b + f5c + plot_layout(widths = c(1, 1, 1.75))) / f5e /
  (f5d + f5f + plot_layout(widths = c(0.85, 1.6))) +
  plot_layout(heights = c(1.05, 0.85, 0.95)) + plot_annotation(tag_levels = "a") & tag_theme
save_nat(fig5, "Fig5", 168)

## ================================================================== Fig. 6
pf <- read_csv(TAB("rnaseq_proxy_programmes.csv"), show_col_types = FALSE) |> filter(!is.na(tumor_likeness)) |> mutate(Tissue = tissue_of(sample), lab = LAB[sample])
lab6a <- bind_rows(pf, line_obst(pf, "tumor_likeness", "proxy")) |>
  mutate(ny = case_when(sample == "P5-Lesion" ~ 0.3, sample == "P2-Tumor" ~ -0.3, TRUE ~ 0))
f6a <- ggplot(pf, aes(tumor_likeness, proxy)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = MUTED, linewidth = PT(0.5)) +
  geom_point(aes(shape = Tissue, fill = Tissue), size = 1.6, stroke = PT(0.6)) +
  ## P5 and the P2 tumour nearly coincide: send their labels to opposite sides
  geom_text_repel(data = lab6a, aes(label = lab), size = PT(5), family = "Helvetica", nudge_y = lab6a$ny,
                  min.segment.length = 0, segment.size = PT(0.3), seed = 4, box.padding = 0.35, force = 4, max.iter = 2e4, max.overlaps = Inf) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3, label = cor_lab(pf$tumor_likeness, pf$proxy), size = PT(5.5), family = "Helvetica") +
  scale_shape_manual(values = SH_T) + scale_fill_manual(values = FILL_T) +
  scale_x_continuous(labels = ML) + scale_y_continuous(labels = ML) +
  labs(x = "Methylation tumour-likeness", y = "Methylome-state expression proxy (RNA)", shape = NULL, fill = NULL) +
  theme_nat() + theme(legend.position = "bottom")
## b. programme shifts lesion minus own adjacent
pdl <- read_csv(TAB("rnaseq_proxy_programmes.csv"), show_col_types = FALSE) |> select(sample, proxy) |> deframe()
dd <- r$progA_delta |> mutate(sample = sub("\\..*", "", sample))
dd$proxy <- map_dbl(dd$sample, function(s) { a <- if (grepl("^P2", s)) "P2-NAT" else paste0(sub("-.*", "", s), "-Adj"); pdl[s] - pdl[a] })
ddl <- dd |> select(sample, proxy, Proliferation, IFNG_response, pEMT, Differentiation) |> pivot_longer(-sample, names_to = "prog", values_to = "d") |>
  mutate(prog = factor(prog, levels = c("proxy", "Proliferation", "IFNG_response", "pEMT", "Differentiation"),
                       labels = c("Methylome-state\nproxy", "Proliferation", "IFN-γ\nresponse", "Partial EMT", "Squamous\ndifferentiation")),
         pt = sub("-.*", "", sample), lab = LAB[sample])
ddl$prog <- factor(ddl$prog, levels = rev(levels(ddl$prog)))
f6b <- ggplot(ddl, aes(d, prog, colour = pt)) + geom_vline(xintercept = 0, colour = HAIR, linewidth = PT(0.5)) +
  geom_point(size = 1.4, position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = PCOL[c("P2", "P4", "P5", "P6", "P7")], breaks = c("P2", "P4", "P5", "P6", "P7"),
                      labels = unname(LAB[c("P2-Tumor", "P4-Lesion", "P5-Lesion", "P6-Lesion", "P7-Lesion")]), name = NULL) +
  scale_x_continuous(labels = ML) +
  labs(y = NULL, x = "Lesion minus own adjacent (score)") + theme_nat() + theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 2))
## c. hallmark sets
cam <- r$cam |> filter(FDR < 0.05) |> mutate(sgn = ifelse(Direction == "Up", 1, -1) * -log10(FDR))
cam <- bind_rows(cam |> filter(Direction == "Up") |> slice_max(sgn, n = 9), cam |> filter(Direction == "Down") |> slice_min(sgn, n = 4)) |>
  mutate(set = str_to_sentence(gsub("_", " ", sub("HALLMARK_", "", set))),
         set = recode(set, "Interferon alpha response" = "Interferon-α response", "Interferon gamma response" = "Interferon-γ response",
                      "Il6 jak stat3 signaling" = "IL-6/JAK/STAT3 signalling", "Tnfa signaling via nfkb" = "TNFα signalling via NF-κB",
                      "Il2 stat5 signaling" = "IL-2/STAT5 signalling", "Epithelial mesenchymal transition" = "Epithelial–mesenchymal transition",
                      "E2f targets" = "E2F targets", "Kras signaling up" = "KRAS signalling up", "G2m checkpoint" = "G2M checkpoint"),
         set = fct_reorder(set, sgn))
f6c <- ggplot(cam, aes(sgn, set, fill = Direction)) + geom_col(width = 0.7) + geom_vline(xintercept = 0, colour = INK, linewidth = PT(0.5)) +
  scale_fill_manual(values = c(Up = RED, Down = BLUE), labels = c(Up = "Higher in lesion", Down = "Lower in lesion"), name = NULL) +
  scale_x_continuous(labels = function(x) abs(x)) +
  labs(x = "−log₁₀ FDR (camera, paired)", y = NULL) + theme_nat() + theme(legend.position = "bottom")
## d. promoter methylation change -> expression change (composition-tracking genes removed)
mxc <- r$mx |> filter(!gene %in% r$nonepi, d_beta >= 0.2 | abs(d_beta) < 0.05) |>
  mutate(cls = factor(ifelse(d_beta >= 0.2, "Gain ≥ 0.2", "Stable"), levels = c("Stable", "Gain ≥ 0.2")),
         panel = factor(LAB[sample], levels = LAB[c("P2-Tumor", "P4-Lesion", "P5-Lesion")]))
pann <- r$mxs_c |> mutate(panel = factor(LAB[sample], levels = levels(mxc$panel)), lab = paste0("Gain vs stable\n", "P = ", ptxt(p_gain_vs_stable)))
nann <- mxc |> count(panel, cls)
f6d <- ggplot(mxc, aes(cls, d_expr, fill = cls)) + geom_hline(yintercept = 0, colour = HAIR, linewidth = PT(0.5)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = PT(0.4)) +
  geom_text(data = nann, aes(cls, -3.6, label = paste0("n=", n)), size = PT(4.8), family = "Helvetica", inherit.aes = FALSE) +
  geom_text(data = pann, aes(x = 1.5, y = 3.6, label = lab), size = PT(5), family = "Helvetica", inherit.aes = FALSE, lineheight = 0.9) +
  facet_wrap(~panel, nrow = 1) + coord_cartesian(ylim = c(-3.9, 4.1)) +
  scale_fill_manual(values = c("Stable" = MID, "Gain ≥ 0.2" = "#f4a6a6"), guide = "none") + scale_y_continuous(labels = ML) +
  labs(x = "Promoter methylation change vs own adjacent (Δβ; stable, |Δβ| < 0.05)", y = "Expression change\n(log₂ fold change)") +
  theme_nat() + theme(axis.text.x = element_text(size = 5))
## e. promoters hypermethylated in both premalignant lesions (P4 AH, P5 MSD): chr19 KZFP cluster
pm <- r$pair_mean |> mutate(cat = case_when(shared_gain & kzfp ~ "19q13 KZFP, gain in both lesions", shared_gain ~ "Other gene, gain in both lesions", TRUE ~ "All other promoters"),
                            cat = factor(cat, levels = c("All other promoters", "Other gene, gain in both lesions", "19q13 KZFP, gain in both lesions")))
kzl <- r$kz_shared |> filter(kzfp, down_both) |> pull(gene)
ann_e <- sprintf("Gain ≥ 0.2 in both lesions: %d promoters\n19q13 KZFP: %d (OR %.0f, P = %s)\nLower in both: %d of %d KZFP vs %d of %d other",
                 nrow(r$kz_shared), sum(r$kz_shared$kzfp), r$f_gain$estimate, ptxt(r$f_gain$p.value),
                 sum(r$kz_shared$kzfp & r$kz_shared$down_both), sum(r$kz_shared$kzfp), sum(!r$kz_shared$kzfp & r$kz_shared$down_both), sum(!r$kz_shared$kzfp))
f6e <- ggplot(pm |> arrange(cat), aes(d_beta, d_expr)) + geom_hline(yintercept = 0, colour = HAIR, linewidth = PT(0.5)) +
  geom_point(aes(colour = cat, size = cat, alpha = cat)) +
  ## gene labels right-aligned in one column at x = 0.17 and separated vertically only
  geom_text_repel(data = pm |> filter(gene %in% kzl), aes(label = gene), size = PT(4.6), family = "Helvetica", fontface = "italic",
                  colour = RED, min.segment.length = 0, segment.size = PT(0.25), max.overlaps = Inf, seed = 8,
                  nudge_x = 0.17 - (pm |> filter(gene %in% kzl))$d_beta, box.padding = 0.3, point.padding = 0.2, force = 3,
                  direction = "y", hjust = 1, xlim = c(NA, 0.2), ylim = c(-3.7, NA), max.iter = 2e4) +
  annotate("text", x = -Inf, y = -Inf, hjust = -0.03, vjust = -0.25, label = ann_e, size = PT(4.8), family = "Helvetica", lineheight = 0.95) +
  scale_colour_manual(values = c("All other promoters" = HAIR, "Other gene, gain in both lesions" = C_O, "19q13 KZFP, gain in both lesions" = RED), name = NULL) +
  scale_size_manual(values = c(0.25, 0.9, 1.2), guide = "none") + scale_alpha_manual(values = c(0.5, 1, 1), guide = "none") +
  scale_x_continuous(labels = ML) + scale_y_continuous(labels = ML) + coord_cartesian(ylim = c(-5.6, 4.2)) +
  labs(x = "Promoter Δβ vs own adjacent (mean of P4 AH and P5 MSD†)", y = "Expression change (log₂ FC, mean)") +
  theme_nat() + theme(legend.position = "bottom", legend.direction = "vertical") + guides(colour = guide_legend(override.aes = list(size = 1.4)))
## f. outcome specificity: methylation-anchored proxy vs pilot lesion programme
pr26 <- g26$cox_tab |> filter(feature == "Methylome_state_proxy")
ci_proxy <- read_csv(TAB("audit_gse26549_forest_exact_ci.csv"), show_col_types = FALSE) |>   # exact Cox CI,
  filter(feature == "Methylome_state_proxy") |> select(lo, hi) |> unlist()                      # scripts/06_audit/01_results_audit.R
rep2 <- g26$rep2
pr7 <- r$a2$Class != "Cancer" & !is.na(r$a2$Progression_status); yv7 <- r$a2$Progression_status[pr7] != "Stable"
oc <- tribble(~metric, ~score, ~est, ~lo, ~hi, ~p,
  "GSE26549\nHR per s.d. (oral cancer)", "Methylome-state proxy", pr26$HR_per_SD, ci_proxy[1], ci_proxy[2], pr26$p,
  "GSE26549\nHR per s.d. (oral cancer)", "Pilot lesion programme", r$obs["HR"], r$obs["lo"], r$obs["hi"], r$obs["p"],
  "GSE227919\nρ with histological grade", "Methylome-state proxy", rep2$trend_rho, NA, NA, rep2$trend_p,
  "GSE227919\nρ with histological grade", "Pilot lesion programme", unname(r$gr$estimate), NA, NA, r$gr$p.value,
  "GSE227919\nAUC for progression", "Methylome-state proxy", rep2$prog_AUC, NA, NA, rep2$prog_p,
  "GSE227919\nAUC for progression", "Pilot lesion programme", r$auc27, NA, NA, wilcox.test(r$prog27[pr7][yv7], r$prog27[pr7][!yv7], exact = FALSE)$p.value) |>
  mutate(score = factor(score, levels = c("Pilot lesion programme", "Methylome-state proxy")),
         metric = factor(metric, levels = unique(metric)), ref = c(1, 1, 0, 0, 0.5, 0.5),
         plab = ifelse(is.na(p), "", paste0("P = ", ptxt(p))))
f6f <- ggplot(oc, aes(est, score, colour = score)) + geom_vline(aes(xintercept = ref), colour = HAIR, linewidth = PT(0.5)) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2, linewidth = PT(0.6), na.rm = TRUE) + geom_point(size = 1.6) +
  geom_text(aes(label = plab), vjust = -0.9, size = PT(4.8), family = "Helvetica", show.legend = FALSE) +
  facet_wrap(~metric, nrow = 1, scales = "free_x") + scale_x_continuous(expand = expansion(mult = 0.18)) +
  scale_colour_manual(values = c(`Methylome-state proxy` = RED, `Pilot lesion programme` = MUTED), guide = "none") +
  scale_y_discrete(expand = expansion(add = 0.7)) + labs(x = NULL, y = NULL) + theme_nat()
fig6 <- (f6a + f6b + f6c + plot_layout(widths = c(0.95, 1.1, 1.2))) /
  (f6d + f6e + plot_layout(widths = c(1.35, 1))) / free(f6f) +
  plot_layout(heights = c(1.1, 1, 0.5)) + plot_annotation(tag_levels = "a") & tag_theme
save_nat(fig6, "Fig6", 165)

## ================================================================== Supplementary Fig. 3
qc <- read_csv(TAB("rnaseq_qc_sex.csv"), show_col_types = FALSE) |> mutate(lab = factor(LAB[sample], levels = LAB[ORD11]), Tissue = tissue_of(sample))
s3a <- ggplot(qc, aes(lab, n_pairs / 1e6, fill = Tissue)) + geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", pct_pseudoaligned)), vjust = -0.4, size = PT(4.6), family = "Helvetica") +
  scale_fill_manual(values = c(Adjacent = "#b7d3f6", Lesion = C_O, Tumour = C_T), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) + labs(x = NULL, y = "Read pairs (millions)\n(label: % pseudoaligned)") +
  theme_nat() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
## strandedness check: scripts/05_pilot_rnaseq/00_strandedness_check.sh
strand <- read_tsv(file.path(RNASEQ_DIR, "strandedness_check.tsv"), show_col_types = FALSE) |>
  transmute(mode = factor(c(unstranded = "Unstranded", `--fr-stranded` = "Forward (--fr-stranded)", `--rf-stranded` = "Reverse (--rf-stranded)")[mode],
                          levels = c("Unstranded", "Forward (--fr-stranded)", "Reverse (--rf-stranded)")),
            pct = p_pseudoaligned) |> arrange(mode)
s3b <- ggplot(strand, aes(mode, pct)) + geom_col(width = 0.6, fill = c(MUTED, BLUE, RED)) +
  geom_text(aes(label = sprintf("%.1f%%", pct)), vjust = -0.4, size = PT(5), family = "Helvetica") +
  scale_y_continuous(limits = c(0, 90), expand = expansion(mult = c(0, 0))) + labs(x = NULL, y = "Pairs pseudoaligned (%)\n(300,000-pair subsample)") +
  theme_nat() + theme(axis.text.x = element_text(angle = 20, hjust = 1))
qc$Y <- rowMeans(qc[, c("RPS4Y1", "DDX3Y", "KDM5D", "EIF1AY")])
s3c <- ggplot(qc, aes(Y, XIST)) + geom_point(aes(shape = Tissue, fill = Tissue), size = 1.5, stroke = PT(0.6)) +
  scale_shape_manual(values = SH_T) + scale_fill_manual(values = FILL_T) +
  coord_cartesian(xlim = c(0, 8), ylim = c(0, 8)) +
  annotate("text", x = 4, y = 4, label = "All 11 libraries:\nfemale pattern", size = PT(5), family = "Helvetica", colour = INK2) +
  labs(x = "Y-linked genes (log<sub>2</sub> TPM + 1)", y = "*XIST* (log<sub>2</sub> TPM + 1)", shape = NULL, fill = NULL) +
  theme_nat() + theme(legend.position = "bottom", axis.title.x = element_markdown(), axis.title.y = element_markdown())
de <- r$de |> mutate(sig = FDR < 0.05, cls = case_when(sig & logFC > 0 ~ "Higher in lesion", sig & logFC < 0 ~ "Lower in lesion", TRUE ~ "FDR ≥ 0.05"))
labg <- c("UBD", "IDO1", "CXCL9", "IFIT1", "BATF2", "KLHDC7B", "CASQ1", "CKMT2", "MMP10", "MKI67")
s3d <- ggplot(de, aes(logFC, -log10(PValue), colour = cls)) + geom_point(size = 0.25, alpha = 0.6) +
  geom_text_repel(data = de |> filter(gene %in% labg), aes(label = gene), size = PT(4.8), family = "Helvetica", fontface = "italic",
                  colour = INK, min.segment.length = 0, segment.size = PT(0.3), seed = 5) +
  scale_colour_manual(values = c("Higher in lesion" = RED, "Lower in lesion" = BLUE, "FDR ≥ 0.05" = HAIR), name = NULL) +
  scale_x_continuous(labels = ML) +
  labs(x = "log₂ fold change, lesion vs own adjacent (5 pairs)", y = "−log₁₀ P") + theme_nat() +
  theme(legend.position = "bottom") + guides(colour = guide_legend(override.aes = list(size = 1.2, alpha = 1)))
mu <- tibble(sample = names(r$muscle_score), m = r$muscle_score) |> mutate(lab = factor(LAB[sample], levels = LAB[ORD11]), Tissue = tissue_of(sample))
s3e <- ggplot(mu, aes(lab, m, fill = Tissue)) + geom_col(width = 0.7, show.legend = FALSE) +
  scale_fill_manual(values = c(Adjacent = "#b7d3f6", Lesion = C_O, Tumour = C_T)) +
  labs(x = NULL, y = "Skeletal-muscle genes\n(mean log₂ TPM + 1)") + theme_nat() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dos <- read_csv(TAB("rnaseq_arm_dosage_all.csv"), show_col_types = FALSE) |> filter(arm_label %in% c("3p", "9p", "13q", "4q", "8p", "11q", "17p", "12p", "12q", "20q")) |>
  mutate(lab = factor(LAB[sample], levels = LAB[ORD11]), arm_label = factor(arm_label, levels = rev(c("3p", "9p", "13q", "4q", "8p", "11q", "17p", "12p", "12q", "20q"))))
s3f <- ggplot(dos, aes(lab, arm_label, fill = dosage)) + geom_tile(colour = "white", linewidth = PT(0.4)) +
  geom_point(data = dos |> filter(beyond_noise) |> mutate(k = "Beyond adjacent-tissue\nnoise floor"), aes(shape = k), size = 0.6, colour = INK) +
  scale_shape_manual(values = 16, name = NULL) +
  scale_fill_gradient2(low = BLUE, mid = MID, high = RED, midpoint = 0, limits = c(-0.6, 0.6), oob = squish, name = "Expression\ndosage", labels = ML) +
  labs(x = NULL, y = NULL) + theme_nat() + theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
sup3 <- (s3a + s3b + s3c + plot_layout(widths = c(1.6, 0.9, 0.9))) / (s3d + s3e + plot_layout(widths = c(1.3, 1))) / s3f +
  plot_layout(heights = c(1, 1, 0.75)) + plot_annotation(tag_levels = "a") & tag_theme
save_nat(sup3, "Supplementary_Fig3", 170)
cat("Saved Fig5, Fig6, Supplementary_Fig3 (PDF + PNG)\n")
