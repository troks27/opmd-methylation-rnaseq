## Publication figures for the methylation story, Nature / Nature Communications
## format. Four main figures + two Supplementary Figures; the panel order
## follows the Results narrative of the manuscript.
##
##   Fig. 1  composition problem -> composition-aware external axis
##   Fig. 2  the tumour-like state: PMD + 19q13 KZFP programme, silencing
##   Fig. 3  copy number from the same arrays; LOH-arm loss; pilot 9p21
##   Fig. 4  transformation prediction; integrated pilot profile; model
##   Supplementary Fig. 1  robustness / audit;  Supplementary Fig. 2  pilot CN + identity
##
## Format: 183 mm wide, <= 170 mm tall; Helvetica 5-7 pt (embedded via the
## Quartz PDF device); 8 pt bold lowercase panel letters; no in-image titles;
## italic gene symbols and P; every colour/shape encoding has a key in its panel.
## Colour roles (validated reference palette): healthy/adjacent #2a78d6,
## OPL/lesion #eb6834, OSCC #1baf7a; composition #eda100/#e87ba4/#008300;
## magnitude = single-hue blue ramp; polarity = blue-grey-red.
##
## Run from the repository root:  Rscript scripts/07_figures/01_figures_1_to_4.R  (or make figures)

source("R/lib/config.R")   # paths and helpers (TAB, FIG, OUT_DIR, GSE_DIR, ...)
suppressPackageStartupMessages({
  library(tidyverse); library(patchwork); library(ggrepel); library(survival); library(scales); library(ggtext)
  library(hugene10sttranscriptcluster.db)
})
select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename
OUT <- NATURE_DIR
rd <- function(f) readRDS(file.path(EXT, f))
ext <- rd("external_methylation.rds"); cnv <- rd("cnv_objects.rds"); st <- rd("story_analyses.rds")
mum <- rd("mumbai_crosscheck.rds"); g26 <- rd("gse26549_progression.rds"); k3 <- rd("gse234379_kzfp.rds")
rna <- rd("external_rnaseq.rds"); au <- rd("audit.rds"); p <- readRDS(file.path(OUT_DIR, "processed.rds"))

## ------------------------------------------------------------------ style
PT <- function(x) x / ggplot2::.pt
INK <- "#1a1a1a"; INK2 <- "#4d4c48"; MUTED <- "#8f8d86"; HAIR <- "#d9d7cf"
C_N <- "#2a78d6"; C_O <- "#eb6834"; C_T <- "#1baf7a"
C_EPI <- "#eda100"; C_FIB <- "#e87ba4"; C_IC <- "#008300"
BLUE <- "#2a78d6"; RED <- "#e34948"; MID <- "#f0efec"
SEQ <- c("#e6effb", "#b7d3f6", "#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281")
CLS <- c("Healthy / adjacent" = C_N, "OPL" = C_O, "OSCC" = C_T)
cls_of <- function(g) c(Normal = "Healthy / adjacent", OPL = "OPL", Tumor = "OSCC")[as.character(g)]
SHP <- c("Healthy / adjacent" = 16, "OPL normal-like" = 1, "OPL tumour-like" = 17, "OSCC" = 16)
theme_nat <- function() {
  theme_classic(base_size = 6, base_family = "Helvetica") +
    theme(line = element_line(linewidth = PT(0.5), colour = INK),
          axis.line = element_line(linewidth = PT(0.5), colour = INK),
          axis.ticks = element_line(linewidth = PT(0.5), colour = INK), axis.ticks.length = unit(0.9, "mm"),
          axis.text = element_text(size = 5.5, colour = INK), axis.title = element_text(size = 6, colour = INK),
          legend.text = element_text(size = 5.2, colour = INK), legend.title = element_text(size = 5.5, colour = INK),
          legend.key.size = unit(2.4, "mm"), legend.key.spacing.y = unit(0.2, "mm"), legend.spacing = unit(0.5, "mm"),
          legend.margin = margin(0, 0, 0, 0), legend.background = element_blank(),
          legend.box.spacing = unit(1, "mm"), legend.box.background = element_blank(),
          strip.background = element_blank(), strip.text = element_text(size = 5.8, colour = INK),
          plot.margin = margin(1.2, 1.5, 1.2, 1.5, "mm"), panel.background = element_blank(), plot.background = element_blank())
}
cbar <- function(h = 9) guide_colourbar(barheight = unit(h, "mm"), barwidth = unit(1.6, "mm"), title.position = "top")
tag_theme <- theme(plot.tag = element_text(size = 8, face = "bold", family = "Helvetica"))
## signed value with two decimals, a true minus sign and no negative zero
sgn2 <- function(x) { s <- sprintf("%+.2f", x); s[s %in% c("+0.00", "-0.00")] <- "0.00"; sub("-", "\u2212", s, fixed = TRUE) }
## Quartz PDF embeds Helvetica (macOS); elsewhere fall back to cairo_pdf.
qpdf <- function(filename, width, height, ...) {
  if (capabilities("aqua")) grDevices::quartz(type = "pdf", file = filename, width = width, height = height)
  else grDevices::cairo_pdf(filename = filename, width = width, height = height)
}
save_nat <- function(p, name, h_mm, w_mm = 183) {
  stopifnot(w_mm <= 183, h_mm <= 170)
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w_mm / 25.4, height = h_mm / 25.4, device = qpdf)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w_mm / 25.4, height = h_mm / 25.4, dpi = 600,
         device = ragg::agg_png, bg = "white")
}
## P-value text: 2 significant figures; < 0.001 as a x 10^-b (Unicode superscripts);
## >= 0.1 to two decimals (so 0.995 prints 0.99, not 1)
sup <- c("-" = "\u207b", "0" = "\u2070", "1" = "\u00b9", "2" = "\u00b2", "3" = "\u00b3", "4" = "\u2074",
         "5" = "\u2075", "6" = "\u2076", "7" = "\u2077", "8" = "\u2078", "9" = "\u2079")
ptxt <- function(pv) vapply(pv, function(x) {
  if (x < 0.001) { xs <- signif(x, 2); e <- floor(log10(xs)); m <- xs / 10^e
    sprintf("%.1f \u00d7 10%s", m, paste(sup[strsplit(as.character(e), "")[[1]]], collapse = "")) }
  else if (x >= 0.1) sprintf("%.2f", round(x, 2)) else formatC(signif(x, 2), format = "fg", digits = 2, flag = "#") }, "")
ORD <- c("P2-NAT", "P2-Tumor", "P3-Adj", "P3-Lesion", "P4-Adj", "P4-Lesion", "P5-Adj", "P5-Lesion")
own_adj <- c(`P2-Tumor` = "P2-NAT", `P3-Lesion` = "P3-Adj", `P4-Lesion` = "P4-Adj", `P5-Lesion` = "P5-Adj")
tissue_of <- function(s) ifelse(grepl("Adj|NAT", s), "Adjacent", ifelse(grepl("Tumor", s), "Tumour", "Lesion"))
short <- function(sn) { pt <- sub("Patient ", "P", sub("_.*", "", sn))
  paste0(pt, "-", ifelse(grepl("tumor", tolower(sn)), "Tumor", ifelse(grepl("lesion", tolower(sn)), "Lesion",
                                                                        ifelse(pt == "P2", "NAT", "Adj")))) }
KZFP <- c("ZNF471", "ZNF382", "ZNF529", "ZNF582", "ZNF671")
SH_T <- c(Adjacent = 21, Lesion = 24, Tumour = 22)
PILOT_FILL <- c("Pilot adjacent" = "white", "Pilot lesion / tumour" = INK)

## pilot KZFP promoter means (same CGI-promoter CpGs as the cohorts)
ann_g <- readRDS(file.path(GSE_DIR, "crucial_objects/ann_auto.rds"))
ann_p <- ann_g[intersect(rownames(ann_g), rownames(p$betas)), ]
keep <- grepl("TSS200|TSS1500|1stExon", ann_p$UCSC_RefGene_Group) & ann_p$Relation_to_Island == "Island"
glp <- strsplit(ann_p$UCSC_RefGene_Name, ";")
kz_p <- sapply(KZFP, function(g) colMeans(p$betas[rownames(ann_p)[keep & vapply(glp, function(x) g %in% x, logical(1))], , drop = FALSE]))
rownames(kz_p) <- short(rownames(kz_p))

## shared cohort tables
sg <- ext$sc_g |> mutate(Group = as.character(Group)) |>
  left_join(st$opl |> select(Sample, methylome_state), by = "Sample") |>
  mutate(Class = cls_of(Group), State = ifelse(Group == "OPL", paste("OPL", methylome_state), Class))
pil <- ext$sc_p |> mutate(Tissue = tissue_of(Sample), Patient = sub("-.*", "", Sample),
                          PilotFill = ifelse(Tissue == "Adjacent", "Pilot adjacent", "Pilot lesion / tumour"))

## ================================================================== FIGURE 1
pc <- tibble(Sample = short(rownames(p$pca$x)), PC1 = p$pca$x[, 1], PC2 = p$pca$x[, 2]) |>
  left_join(ext$decon |> select(Sample, Epi = Epi_1st, Fib = Fib_1st, IC = IC_1st), by = "Sample") |>
  mutate(Tissue = tissue_of(Sample))
pct <- round(100 * p$pca$sdev^2 / sum(p$pca$sdev^2))
f1a <- ggplot(pc, aes(PC1, PC2)) +
  geom_point(aes(fill = Epi, shape = Tissue), size = 1.9, colour = INK, stroke = PT(0.4)) +
  geom_text_repel(aes(label = Sample), size = PT(5), colour = INK2, seed = 1, box.padding = 0.3, point.padding = 0.2,
                  min.segment.length = 0.25, segment.size = PT(0.3)) +
  scale_fill_gradientn(colours = SEQ, limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0%", "50%", "100%"),
                       name = "Epithelial\nfraction", guide = cbar(8)) +
  scale_shape_manual(values = SH_T, name = "Tissue") +
  guides(shape = guide_legend(override.aes = list(fill = "white", size = 1.5), order = 1)) +
  scale_y_continuous(expand = expansion(mult = 0.12)) +
  labs(x = sprintf("PC1 (%d%%)", pct[1]), y = sprintf("PC2 (%d%%)", pct[2])) + theme_nat()
## labels are repelled from the points and from the fitted line (unlabelled obstacles along it)
fit_b <- lm(PC1 ~ Epi, data = pc)
obst_b <- tibble(Epi = seq(min(pc$Epi), max(pc$Epi), length.out = 60)) |> mutate(PC1 = predict(fit_b, pick(Epi)), Sample = "")
f1b <- ggplot(pc, aes(Epi, PC1)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = MUTED, linewidth = PT(0.6)) +
  geom_point(aes(shape = Tissue), size = 1.6, fill = INK2, colour = "white", stroke = PT(0.3)) +
  geom_text_repel(data = bind_rows(pc, obst_b), aes(label = Sample), size = PT(5), colour = INK2, seed = 2, segment.size = PT(0.3),
                  box.padding = 0.4, point.padding = 0.15, force = 4, max.iter = 2e4, min.segment.length = 0.1, max.overlaps = Inf) +
  scale_shape_manual(values = SH_T, name = "Tissue") +
  guides(shape = guide_legend(override.aes = list(fill = "white", colour = INK, size = 1.5))) +
  annotate("text", x = 0.02, y = 29, hjust = 0, vjust = 1, size = PT(5.5), parse = TRUE,
           label = sprintf("italic(r)=='%.2f'", ext$r_pc1)) +
  scale_x_continuous(limits = c(0, 1), labels = percent) +
  labs(x = "Epithelial fraction (EpiDISH)", y = "PC1") + theme_nat() + theme(legend.position = "right")
comp <- pc |> select(Sample, Epithelial = Epi, Fibroblast = Fib, Immune = IC) |>
  pivot_longer(-Sample, names_to = "Cell", values_to = "Fraction") |>
  mutate(Sample = factor(Sample, ORD), Cell = factor(Cell, c("Immune", "Fibroblast", "Epithelial")))
## horizontal bars: sample names read left-to-right and add no height to the row
comp <- comp |> mutate(Sample = factor(Sample, rev(ORD)))
f1c <- ggplot(comp, aes(Fraction, Sample, fill = Cell)) +
  geom_col(width = 0.78, colour = "white", linewidth = PT(0.5)) +
  scale_fill_manual(values = c(Epithelial = C_EPI, Fibroblast = C_FIB, Immune = C_IC), name = NULL,
                    breaks = c("Epithelial", "Fibroblast", "Immune")) +
  scale_x_continuous(labels = percent, expand = c(0, 0)) +
  labs(x = "Estimated fraction (EpiDISH)", y = NULL) + theme_nat() +
  theme(legend.position = "top", legend.justification = "left", legend.margin = margin(0, 0, 0.5, 0, "mm"))
## d. external continuum
s3 <- k3$s3 |> mutate(Group = group, Class = cls_of(group), State = Class)
PHIST <- c(P2 = "OSCC", P3 = "PVL", P4 = "atypical hyperplasia", P5 = "dysplasia, progressed")
rows <- c(paste0("Pilot ", c("P5", "P4", "P3", "P2"), " (", PHIST[c("P5", "P4", "P3", "P2")], ")"),
          "GSE234379 OSCC", "GSE234379 adjacent", "GSE204943 OSCC†", "GSE204943 OPL", "GSE204943 healthy†")
pts <- bind_rows(sg |> mutate(row = recode(Group, Normal = rows[9], OPL = rows[8], Tumor = rows[7])) |>
                   select(row, tumor_likeness, Class, State),
                 s3 |> mutate(row = ifelse(group == "Tumor", rows[5], rows[6])) |> select(row, tumor_likeness, Class, State)) |>
  mutate(y = match(row, rows))
pil_d <- pil |> mutate(y = match(paste0("Pilot ", Patient, " (", PHIST[Patient], ")"), rows))
arr <- tibble(lesion = names(own_adj), adj = unname(own_adj)) |>
  left_join(pil_d |> select(lesion = Sample, xend = tumor_likeness, y), by = "lesion") |>
  left_join(pil_d |> select(adj = Sample, x = tumor_likeness), by = "adj")
set.seed(3)
f1d <- ggplot() +
  annotate("rect", xmin = st$valley, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#fdf3ec") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = 4.5, fill = "#f4f3ef") +
  geom_vline(xintercept = c(0, 1), colour = HAIR, linewidth = PT(0.4)) +
  geom_vline(xintercept = st$valley, linetype = "22", colour = MUTED, linewidth = PT(0.5)) +
  geom_jitter(data = pts, aes(tumor_likeness, y, colour = Class, shape = State), height = 0.17, width = 0,
              size = 0.7, stroke = PT(0.5), alpha = 0.85) +
  geom_segment(data = arr, aes(x = x, xend = xend, y = y, yend = y), colour = INK2, linewidth = PT(0.6),
               arrow = arrow(length = unit(0.9, "mm"), type = "closed")) +
  geom_point(data = pil_d, aes(tumor_likeness, y, fill = PilotFill), shape = 21, size = 1.4, colour = INK, stroke = PT(0.4)) +
  annotate("text", x = st$valley + 0.03, y = 9.62, label = "Tumour-like range", hjust = 0, size = PT(5), colour = C_O) +
  scale_colour_manual(values = CLS, name = NULL) +
  scale_shape_manual(values = SHP, name = NULL, breaks = c("OPL normal-like", "OPL tumour-like")) +
  scale_fill_manual(values = PILOT_FILL, name = NULL) +
  guides(colour = guide_legend(order = 1, override.aes = list(size = 1.3, alpha = 1)),
         shape = guide_legend(order = 2, override.aes = list(size = 1.3, colour = C_O)),
         fill = guide_legend(order = 3, override.aes = list(shape = 21, size = 1.4))) +
  scale_y_continuous(breaks = seq_along(rows), labels = rows, expand = expansion(add = c(0.5, 0.7))) +
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2)) + coord_cartesian(xlim = c(-0.5, 2.15)) +
  labs(x = "Tumour-likeness (0 = healthy-mucosa mean, 1 = OSCC mean)", y = NULL) + theme_nat() +
  theme(axis.ticks.y = element_blank(), axis.line.y = element_blank(), legend.position = "right")
## e. composition independence
f1e <- ggplot() +
  geom_point(data = sg, aes(Epi, tumor_likeness, colour = Class, shape = State), size = 0.75, stroke = PT(0.5), alpha = 0.85) +
  geom_point(data = pil, aes(Epi, tumor_likeness, fill = PilotFill), shape = 21, size = 1.4, colour = INK, stroke = PT(0.4)) +
  geom_text_repel(data = bind_rows(pil |> select(Epi, tumor_likeness, Sample), sg |> transmute(Epi, tumor_likeness, Sample = "")),
                  aes(Epi, tumor_likeness, label = Sample), size = PT(4.8), colour = INK2, seed = 4, segment.size = PT(0.3),
                  min.segment.length = 0.1, box.padding = 0.35, point.padding = 0.45, force = 3, max.iter = 2e4, max.overlaps = Inf) +
  scale_colour_manual(values = CLS, name = NULL) +
  scale_shape_manual(values = SHP, name = NULL, breaks = c("OPL normal-like", "OPL tumour-like")) +
  scale_fill_manual(values = PILOT_FILL, name = NULL) +
  guides(colour = guide_legend(order = 1, override.aes = list(size = 1.3, alpha = 1)),
         shape = guide_legend(order = 2, override.aes = list(size = 1.3, colour = C_O)),
         fill = guide_legend(order = 3, override.aes = list(shape = 21, size = 1.4))) +
  scale_x_continuous(labels = percent, limits = c(0, 0.92)) +
  labs(x = "Epithelial fraction", y = "Tumour-likeness") + theme_nat() + theme(legend.position = "right")
## f. pilot ranking change
tum <- pc |> filter(Sample == "P2-Tumor")
rk <- pc |> filter(Sample %in% c("P3-Lesion", "P4-Lesion", "P5-Lesion")) |>
  mutate(dist = sqrt((PC1 - tum$PC1)^2 + (PC2 - tum$PC2)^2)) |>
  left_join(ext$sc_p |> select(Sample, delta = delta_vs_own_adj), by = "Sample") |>
  mutate(naive = rank(dist), adj = rank(-delta),
         lab = c(`P3-Lesion` = "P3 (PVL)", `P4-Lesion` = "P4 (AH)", `P5-Lesion` = "P5 (progressed)")[Sample])
sl <- rk |> pivot_longer(c(naive, adj), names_to = "view", values_to = "rank") |>
  mutate(x = ifelse(view == "naive", 1, 2), val = ifelse(view == "naive", sprintf("%.1f", dist), sgn2(delta)))
f1f <- ggplot(sl, aes(x, rank, group = Sample)) +
  geom_line(aes(colour = Sample == "P5-Lesion"), linewidth = PT(1.2)) +
  geom_point(aes(colour = Sample == "P5-Lesion"), size = 1.6) +
  geom_text(data = filter(sl, x == 1), aes(label = val), nudge_x = -0.1, hjust = 1, size = PT(5), colour = INK2) +
  geom_text(data = filter(sl, x == 2), aes(label = val), nudge_x = 0.1, hjust = 0, size = PT(5), colour = INK2) +
  geom_text(data = filter(sl, x == 2), aes(x = 2.62, label = lab), hjust = 0, size = PT(5.2), colour = INK) +
  scale_colour_manual(values = c(`TRUE` = C_O, `FALSE` = MUTED), guide = "none") +
  scale_y_reverse(breaks = 1:3, expand = expansion(add = 0.35)) +
  scale_x_continuous(breaks = 1:2, labels = c("PCA distance\nto tumour", "Shift vs own\nadjacent"), limits = c(0.55, 3.75)) +
  labs(x = NULL, y = "Rank (1 = closest to tumour)") + theme_nat() +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank())
## free(): f1d's long row labels must not push the other panels' axes inward
fig1 <- (f1a + f1b + f1c + plot_layout(widths = c(1.05, 1.05, 1))) /
  free(f1d) / (f1e + f1f + plot_layout(widths = c(1.1, 0.95))) +
  plot_layout(heights = c(1, 1.15, 1)) + plot_annotation(tag_levels = "a") & tag_theme
save_nat(fig1, "Fig1", 170)

## ================================================================== FIGURE 2
gk <- au$cb |> mutate(Class = cls_of(Group), State = ifelse(Group == "OPL", paste("OPL", methylome_state), Class))
pil_k <- ext$pz_p |> select(Sample, pmd_contrast) |>
  mutate(KZFP = rowMeans(kz_p)[Sample], Tissue = tissue_of(Sample),
         PilotFill = ifelse(Tissue == "Adjacent", "Pilot adjacent", "Pilot lesion / tumour"))
arr_k <- tibble(lesion = names(own_adj), adj = unname(own_adj)) |>
  left_join(pil_k |> select(lesion = Sample, xend = pmd_contrast, yend = KZFP), by = "lesion") |>
  left_join(pil_k |> select(adj = Sample, x = pmd_contrast, y = KZFP), by = "adj")
rho_o <- with(filter(gk, Group == "OPL"), cor(-pmd_contrast, KZFP, method = "spearman"))
f2a <- ggplot() +
  geom_point(data = gk, aes(-pmd_contrast, KZFP, colour = Class, shape = State), size = 0.75, stroke = PT(0.5), alpha = 0.85) +
  geom_segment(data = arr_k, aes(-x, y, xend = -xend, yend = yend), colour = INK2, linewidth = PT(0.5),
               arrow = arrow(length = unit(0.9, "mm"), type = "closed")) +
  geom_point(data = pil_k, aes(-pmd_contrast, KZFP, fill = PilotFill), shape = 21, size = 1.4, colour = INK, stroke = PT(0.4)) +
  geom_text(data = filter(pil_k, Tissue != "Adjacent") |>
              mutate(dx = c(`P2-Tumor` = 0.012, `P3-Lesion` = 0.012, `P4-Lesion` = 0.012, `P5-Lesion` = 0.012)[Sample],
                     dy = c(`P2-Tumor` = -0.035, `P3-Lesion` = -0.03, `P4-Lesion` = -0.03, `P5-Lesion` = 0.03)[Sample]),
            aes(-pmd_contrast + dx, KZFP + dy, label = sub("-.*", "", Sample)), hjust = 0, size = PT(5), colour = INK) +
  scale_colour_manual(values = CLS, name = NULL) +
  scale_shape_manual(values = SHP, name = NULL, breaks = c("OPL normal-like", "OPL tumour-like")) +
  scale_fill_manual(values = PILOT_FILL, name = NULL) +
  guides(colour = guide_legend(order = 1, override.aes = list(size = 1.3, alpha = 1)),
         shape = guide_legend(order = 2, override.aes = list(size = 1.3, colour = C_O)),
         fill = guide_legend(order = 3, override.aes = list(shape = 21, size = 1.4))) +
  annotate("text", x = 0.31, y = 0.04, hjust = 1, size = PT(5.5), parse = TRUE, label = sprintf("'OPL '*rho=='%.2f'", rho_o)) +
  labs(x = "PMD hypomethylation (outside − inside PMDs)", y = "KZFP promoter methylation (β)") + theme_nat() +
  theme(legend.position = "right")
f2b <- ggplot(gk, aes(Epi, KZFP, colour = Class, shape = State)) +
  geom_point(size = 0.75, stroke = PT(0.5), alpha = 0.85) +   # identical keys to f2a so the legends merge
  scale_colour_manual(values = CLS, name = NULL) +
  scale_shape_manual(values = SHP, name = NULL, breaks = c("OPL normal-like", "OPL tumour-like")) +
  guides(colour = guide_legend(order = 1, override.aes = list(size = 1.3, alpha = 1)),
         shape = guide_legend(order = 2, override.aes = list(size = 1.3, colour = C_O))) +
  scale_x_continuous(labels = percent) +
  labs(x = "Epithelial fraction", y = "KZFP promoter methylation (\u03b2)") + theme_nat() +
  theme(legend.position = "right")
hm <- k3$gk |> mutate(col = recode(state, Normal = "Healthy\nmucosa", `OPL normal-like` = "OPL\nnormal-like",
                                   `OPL tumour-like` = "OPL\ntumour-like", Tumor = "OSCC")) |>
  select(col, all_of(KZFP)) |> pivot_longer(all_of(KZFP), names_to = "gene", values_to = "beta") |>
  group_by(col, gene) |> summarise(beta = median(beta), .groups = "drop") |>
  mutate(col = factor(col, c("Healthy\nmucosa", "OPL\nnormal-like", "OPL\ntumour-like", "OSCC")),
         gene = factor(paste0("*", gene, "*"), rev(paste0("*", KZFP, "*"))))
f2c <- ggplot(hm, aes(col, gene, fill = beta)) + geom_tile(colour = "white", linewidth = PT(0.8)) +
  geom_text(aes(label = sprintf("%.2f", beta), colour = beta > 0.38), size = PT(5)) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = INK), guide = "none") +
  scale_fill_gradientn(colours = SEQ, limits = c(0, 0.6), breaks = c(0, 0.3, 0.6), name = "Median β", guide = cbar(9)) +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  labs(x = "GSE204943", y = NULL) + theme_nat() +
  theme(axis.text.y = element_markdown(), axis.line = element_blank(), axis.ticks = element_blank())
pk <- k3$pil |> select(gene, P2 = d_P2_tumor, P3 = d_P3, P4 = d_P4, P5 = d_P5_progressor) |>
  pivot_longer(-gene, names_to = "Patient", values_to = "delta") |>
  mutate(gene = factor(paste0("*", gene, "*"), rev(paste0("*", KZFP, "*"))),
         Patient = factor(Patient, c("P2", "P3", "P4", "P5"), c("P2\nOSCC", "P3\nPVL", "P4\nAH", "P5\nprogr.")))
f2d <- ggplot(pk, aes(Patient, gene, fill = delta)) + geom_tile(colour = "white", linewidth = PT(0.8)) +
  geom_text(aes(label = sgn2(delta)), size = PT(5), colour = INK) +
  scale_fill_gradient2(low = BLUE, mid = MID, high = RED, midpoint = 0, limits = c(-0.6, 0.6), breaks = c(-0.5, 0, 0.5),
                       name = "Δβ vs own\nadjacent", guide = cbar(9)) +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  labs(x = "Pilot", y = NULL) + theme_nat() +
  theme(axis.text.y = element_markdown(), axis.line = element_blank(), axis.ticks = element_blank())
pr <- k3$s3 |> mutate(KZFP = rowMeans(across(all_of(KZFP)))) |> group_by(patient) |> filter(n() == 2) |> ungroup() |>
  mutate(group = factor(ifelse(group == "Tumor", "OSCC", "Adjacent"), c("Adjacent", "OSCC")))
pr_up <- pr |> select(patient, group, KZFP) |> pivot_wider(names_from = group, values_from = KZFP)
pr_p <- wilcox.test(pr_up$OSCC, pr_up$Adjacent, paired = TRUE)$p.value
f2e <- ggplot(pr, aes(group, KZFP)) +
  geom_line(aes(group = patient), colour = "#c9c7bf", linewidth = PT(0.5)) +
  geom_point(aes(colour = group), size = 0.9) +
  stat_summary(aes(group = 1), fun = median, geom = "line", colour = INK, linewidth = PT(1.2)) +
  scale_colour_manual(values = c(Adjacent = C_N, OSCC = C_T), name = NULL) +
  annotate("text", x = 1.5, y = max(pr$KZFP) * 1.04, size = PT(5.5), parse = TRUE,
           label = sprintf("'%d pairs, '*italic(P)*' = %s'", nrow(pr_up), ptxt(pr_p))) +
  labs(x = "GSE234379 (São Paulo)", y = "KZFP promoter methylation (β)") + theme_nat() +
  theme(legend.position = "none")
wi <- au$within |> select(gene, all = rho_all, OSCC = rho_OSCC, OPL = rho_OPL) |>
  pivot_longer(-gene, names_to = "set", values_to = "rho") |>
  mutate(set = factor(set, c("all", "OPL", "OSCC"), with(au$within[1, ], sprintf(c("All (n = %d)", "OPL (n = %d)", "OSCC (n = %d)"),
                                                                                   c(n_all, n_OPL, n_OSCC)))),
         gene = factor(paste0("*", gene, "*"), rev(paste0("*", KZFP, "*"))))
f2f <- ggplot(wi, aes(rho, gene, colour = set, shape = set)) +
  geom_vline(xintercept = 0, colour = MUTED, linewidth = PT(0.5)) +
  geom_vline(xintercept = median(mum$me_tab$rho, na.rm = TRUE), colour = MUTED, linewidth = PT(0.5), linetype = "22") +
  geom_point(size = 1.2, stroke = PT(0.6), position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = c(INK, C_O, C_T), name = "GSE85195 samples") +
  scale_shape_manual(values = c(16, 17, 15), name = "GSE85195 samples") +
  annotate("text", x = median(mum$me_tab$rho, na.rm = TRUE) + 0.015, y = 5.45, hjust = 0, size = PT(4.8), colour = MUTED,
           label = "median, all\nvariable\npromoters", lineheight = 0.9) +
  scale_x_continuous(limits = c(-0.9, 0.12), breaks = c(-0.8, -0.4, 0)) +
  coord_cartesian(clip = "off") +
  labs(x = "Spearman ρ, methylation vs expression", y = NULL) + theme_nat() +
  theme(axis.text.y = element_markdown(), legend.position = "right")
fig2 <- (f2a + f2b + plot_layout(widths = c(1.15, 1), guides = "collect")) /
  (f2c + f2d + plot_layout(widths = c(1.15, 1))) /
  (f2e + f2f + plot_layout(widths = c(0.75, 1.25))) +
  plot_layout(heights = c(1.05, 0.85, 0.95)) + plot_annotation(tag_levels = "a") & tag_theme
save_nat(fig2, "Fig2", 160)

## ================================================================== FIGURE 3
xp <- mum$xp |> mutate(Class = cls_of(Group))
r_by <- mum$per_sample |> group_by(Group) |> summarise(r = median(r_own), n = n(), .groups = "drop")
f3a <- ggplot(xp, aes(conumee, aCGH)) +
  geom_abline(slope = 1, intercept = 0, colour = MUTED, linetype = "22", linewidth = PT(0.5)) +
  geom_hline(yintercept = 0, colour = HAIR, linewidth = PT(0.4)) + geom_vline(xintercept = 0, colour = HAIR, linewidth = PT(0.4)) +
  geom_point(aes(colour = Class), size = 0.45, alpha = 0.6, stroke = 0) +
  scale_colour_manual(values = CLS[c("OPL", "OSCC")], name = NULL,
                      labels = c(OPL = sprintf("OPL (%d)", r_by$n[r_by$Group == "OPL"]), OSCC = sprintf("OSCC (%d)", r_by$n[r_by$Group == "Tumor"]))) +
  annotate("text", x = 0.62, y = -0.62, hjust = 1, vjust = 0, size = PT(5), lineheight = 0.95,
           label = sprintf("Median per-lesion r\nOSCC %.2f, OPL %.2f\nshuffled null %.2f",
                           r_by$r[r_by$Group == "Tumor"], r_by$r[r_by$Group == "OPL"], median(mum$null_r))) +
  coord_equal(xlim = c(-0.62, 0.62), ylim = c(-0.62, 0.62)) +
  guides(colour = guide_legend(override.aes = list(size = 1.3, alpha = 1))) +
  labs(x = "Methylation array (log2)", y = "Array CGH (log2)") + theme_nat() +
  theme(legend.position = "inside", legend.position.inside = c(0.22, 0.9))
risk_arms <- c("3p", "9p", "13q", "4q", "8p", "11q", "17p")
qtab <- au$arm_q |> select(arm = feature, q)
lev <- sprintf("%s\nq = %s", risk_arms, formatC(signif(qtab$q[match(risk_arms, qtab$arm)], 2), format = "fg", digits = 2))
ob <- mum$opl_c |> filter(!is.na(methylome_state)) |> select(id, methylome_state, all_of(risk_arms)) |>
  pivot_longer(all_of(risk_arms), names_to = "arm", values_to = "lr") |> left_join(qtab, by = "arm") |>
  mutate(lab = factor(sprintf("%s\nq = %s", arm, formatC(signif(q, 2), format = "fg", digits = 2)), lev),
         State = paste("OPL", methylome_state))
f3b <- ggplot(ob, aes(State, lr)) +
  geom_hline(yintercept = 0, colour = HAIR, linewidth = PT(0.4)) +
  geom_jitter(aes(shape = State), width = 0.18, size = 0.85, colour = C_O, stroke = PT(0.5)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.55, linewidth = PT(0.5), colour = INK) +
  scale_shape_manual(values = SHP, name = NULL) +
  coord_cartesian(ylim = c(-0.12, 0.08)) + facet_wrap(~lab, nrow = 1) +
  labs(x = "OPL methylome group (GSE204943); arm copy number by array CGH", y = "Arm log2 ratio (GSE85514)") + theme_nat() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "top",
        legend.justification = "right", legend.margin = margin(0, 0, 0, 0, "mm"), panel.spacing = unit(1, "mm"))
arm_order <- paste0(rep(1:22, each = 2), c("p", "q"))
cn <- cnv$arms |> mutate(Sample = factor(Sample, rev(ORD)), arm_label = factor(arm_label, arm_order)) |> filter(!is.na(arm_label))
f3c <- ggplot(cn, aes(arm_label, Sample, fill = pmax(pmin(LPO, 0.25), -0.25))) +
  geom_tile(colour = "white", linewidth = PT(0.5)) +
  geom_point(data = filter(cn, call != ""), aes(shape = "Called arm"), size = 0.4, colour = INK) +
  scale_shape_manual(values = c(`Called arm` = 16), name = NULL) +
  scale_fill_gradient2(low = BLUE, mid = MID, high = RED, midpoint = 0, limits = c(-0.25, 0.25), breaks = c(-0.2, 0, 0.2),
                       name = "log2 ratio", guide = cbar(8)) +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  labs(x = "Chromosome arm", y = NULL) + theme_nat() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, size = 4.8))
z9 <- cnv$bins_lpo |> filter(chr == "chr9", mid > 20e6, mid < 24e6) |>
  pivot_longer(all_of(ORD), names_to = "Sample", values_to = "lr") |>
  mutate(Line = ifelse(Sample %in% c("P3-Adj", "P4-Lesion", "P5-Lesion"), Sample, "Other pilot samples"))
f3d <- ggplot(z9, aes(mid / 1e6, lr, group = Sample, colour = Line, linewidth = Line)) +
  annotate("rect", xmin = 21.80, xmax = 22.01, ymin = -Inf, ymax = Inf, fill = "#fdf1d6") +
  annotate("text", x = 21.905, y = 0.44, label = "italic(MTAP)*'-'*italic(CDKN2A)*'-'*italic(CDKN2B)", parse = TRUE,
           size = PT(5), colour = INK2) +
  geom_hline(yintercept = 0, colour = HAIR, linewidth = PT(0.4)) +
  geom_line() +
  scale_colour_manual(values = c(`P5-Lesion` = C_O, `P4-Lesion` = C_T, `P3-Adj` = C_N, `Other pilot samples` = "#c9c7bf"),
                      breaks = c("P5-Lesion", "P4-Lesion", "P3-Adj", "Other pilot samples"), name = NULL) +
  scale_linewidth_manual(values = c(`P5-Lesion` = PT(1), `P4-Lesion` = PT(1), `P3-Adj` = PT(1), `Other pilot samples` = PT(0.6)),
                         guide = "none") +
  coord_cartesian(ylim = c(-0.7, 0.5)) +
  labs(x = "Chromosome 9 (Mb, hg38)", y = "log2 ratio") + theme_nat() +
  theme(legend.position = "inside", legend.position.inside = c(0.22, 0.2))
cd <- ext$gene_tab |> filter(gene == "CDKN2A") |> select(all_of(ORD)) |>
  pivot_longer(everything(), names_to = "Sample", values_to = "prom") |>
  left_join(cnv$del9 |> select(Sample, lr = lpo_9p21_median, fail = fail_ratio), by = "Sample") |> mutate(Tissue = tissue_of(Sample))
f3e <- ggplot(cd, aes(lr, prom)) +
  geom_point(aes(shape = Tissue, fill = fail), size = 1.7, colour = INK, stroke = PT(0.4)) +
  geom_text_repel(aes(label = Sample), size = PT(4.8), colour = INK2, seed = 6, segment.size = PT(0.3),
                  box.padding = 0.35, max.overlaps = 30, min.segment.length = 0.1) +
  scale_shape_manual(values = SH_T, name = "Tissue") +
  scale_fill_gradientn(colours = SEQ[c(1, 4, 7)], name = "9p21 probe\nfailure (×)", guide = cbar(8)) +
  guides(shape = guide_legend(order = 1, override.aes = list(fill = "white", size = 1.4))) +
  scale_y_continuous(limits = c(0, 0.3)) + scale_x_continuous(limits = c(-0.5, 0.25)) +
  labs(x = "9p21 log2 ratio", y = "*CDKN2A* promoter methylation (β)") + theme_nat() +
  theme(axis.title.y = element_markdown(), legend.position = "right")
fig3 <- (f3a + f3b + plot_layout(widths = c(0.78, 2))) / f3c / (f3d + f3e + plot_layout(widths = c(1, 1.05))) +
  plot_layout(heights = c(1.05, 0.62, 1)) + plot_annotation(tag_levels = "a") & tag_theme
save_nat(fig3, "Fig3", 160)

## ================================================================== FIGURE 4
bx <- tibble(x = c(0.6, 2.2, 3.8), xmin = x - 0.72, xmax = x + 0.72,
  head = c("GSE204943 (Mumbai)", "GSE85195 (same OPLs)", "GSE26549 (Houston)"),
  body = c("EPIC methylation, 22 OPL:\nmethylome tumour-likeness", "Expression, 15 OPL:\n50 up / 50 down gene proxy",
           "86 OPL, 35 oral cancers:\nproxy applied unchanged"),
  fill = c("#fdf3ec", "#fdf3ec", "#eaf2fc"))
f4a <- ggplot(bx) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 0.2, ymax = 1.8, fill = fill), colour = INK, linewidth = PT(0.5)) +
  geom_text(aes(x, 1.42, label = head), fontface = "bold", size = PT(5.8), colour = INK) +
  geom_text(aes(x, 0.78, label = body), size = PT(5.3), colour = INK, lineheight = 0.95) +
  geom_segment(data = tibble(x = c(1.34, 2.94), xend = c(1.46, 3.06)), aes(x = x, xend = xend, y = 1, yend = 1),
               arrow = arrow(length = unit(1.1, "mm"), type = "closed"), colour = INK2, linewidth = PT(0.6)) +
  annotate("text", x = 1.4, y = 1.98, label = "train", size = PT(5), colour = INK2, fontface = "italic") +
  annotate("text", x = 3.0, y = 1.98, label = "test", size = PT(5), colour = INK2, fontface = "italic") +
  scale_fill_identity() + coord_cartesian(xlim = c(-0.15, 4.55), ylim = c(0.1, 2.1), expand = FALSE) + theme_void()
d26 <- bind_cols(g26$clin, g26$feat) |> mutate(grp = ifelse(Methylome_state_proxy > median(Methylome_state_proxy), "Above median", "Below median"))
sf <- survfit(Surv(ocfs_years, event) ~ grp, data = d26)
kmd <- tibble(time = sf$time, surv = sf$surv, n.censor = sf$n.censor, strata = rep(sub("grp=", "", names(sf$strata)), sf$strata))
kmd <- bind_rows(kmd |> group_by(strata) |> summarise(time = 0, surv = 1, n.censor = 0), kmd) |> arrange(strata, time)
cx <- au$cox$ci; cxh <- au$cox$ci_h
p_cont <- summary(coxph(Surv(ocfs_years, event) ~ scale(Methylome_state_proxy), data = d26))$coefficients[1, 5]
KMC <- c("Above median" = C_O, "Below median" = C_N)
f4b_km <- ggplot(kmd, aes(time, surv, colour = strata)) +
  geom_step(linewidth = PT(1)) +
  geom_point(data = filter(kmd, n.censor > 0), shape = 3, size = 0.7, stroke = PT(0.5)) +
  scale_colour_manual(values = KMC, name = "Methylome-state proxy") +
  scale_y_continuous(labels = percent, expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(0, 12, 2), expand = expansion(add = c(0.8, 0.3))) +
  coord_cartesian(xlim = c(0, 12.2), ylim = c(0, 1)) +
  annotate("text", x = 0.4, y = 0.03, hjust = 0, vjust = 0, size = PT(5.2), lineheight = 1,
           label = sprintf("HR %.2f (95%% CI %.2f–%.2f) per s.d., P = %.4f\nHistology-adjusted HR %.2f (%.2f–%.2f)",
                           cx[1, 1], cx[1, 3], cx[1, 4], p_cont, cxh[1, 1], cxh[1, 3], cxh[1, 4])) +
  labs(x = NULL, y = "Oral cancer-free") + theme_nat() +
  theme(legend.position = "inside", legend.position.inside = c(0.72, 0.88), axis.text.x = element_blank())
times <- seq(0, 12, 2)
risk <- map_dfr(names(KMC), function(g) tibble(strata = g, time = times, n = sapply(times, function(t) sum(d26$ocfs_years[d26$grp == g] >= t))))
f4b_risk <- ggplot(risk, aes(time, factor(strata, rev(names(KMC))), label = n, colour = strata)) +
  geom_text(size = PT(5)) +
  scale_colour_manual(values = KMC, guide = "none") +
  scale_x_continuous(breaks = seq(0, 12, 2), expand = expansion(add = c(0.8, 0.3))) +
  coord_cartesian(xlim = c(0, 12.2)) +
  scale_y_discrete(labels = c("Below median" = "Below", "Above median" = "Above")) +
  labs(x = "Years", y = "No. at risk") + theme_nat() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())
f4b <- f4b_km / f4b_risk + plot_layout(heights = c(1, 0.26))
kzc <- local({
  meta <- readLines(gzfile(file.path(XDATA, "GSE26549_series_matrix.txt.gz")))
  b <- grep("series_matrix_table_begin", meta); e <- grep("series_matrix_table_end", meta)
  X <- as.matrix(read.delim(text = meta[(b + 1):(e - 1)], check.names = FALSE, row.names = 1))[, g26$clin$gsm]
  sy <- suppressMessages(mapIds(hugene10sttranscriptcluster.db, rownames(X), "SYMBOL", "PROBEID", multiVals = "first"))
  map_dfr(KZFP, function(g) { x <- -as.numeric(scale(colMeans(X[which(sy == g), , drop = FALSE])))
    s <- summary(coxph(Surv(g26$clin$ocfs_years, g26$clin$event) ~ x))
    tibble(feature = g, HR = s$conf.int[1, 1], lo = s$conf.int[1, 3], hi = s$conf.int[1, 4], p = s$coefficients[1, 5]) })
})
lab_map <- c(Methylome_state_proxy = "Methylome-state proxy", KZFP_chr19_hyper = "KZFP set (low)",
             loss_3p = "3p dosage loss", loss_9p = "9p dosage loss", loss_core_3p9p = "3p + 9p dosage loss",
             Dysplastic_epithelial = "Dysplastic epithelial", Invasive_epithelial = "Invasive epithelial",
             Cycling_epithelial = "Cycling epithelial", LRRC15_CAF = "*LRRC15*<sup>+</sup> CAF", FOLR2_TAM = "*FOLR2*<sup>+</sup> TAM",
             Treg = "Regulatory T cell", mregDC = "mregDC")
fo <- au$forest |> filter(feature %in% names(lab_map)) |>
  mutate(inv = feature == "KZFP_chr19_hyper", HR0 = HR,
         HR = ifelse(inv, 1 / HR0, HR0), lo2 = ifelse(inv, 1 / hi, lo), hi2 = ifelse(inv, 1 / lo, hi), lo = lo2, hi = hi2,
         group = case_when(feature %in% c("Methylome_state_proxy", "KZFP_chr19_hyper") ~ "Methylation",
                           grepl("loss", feature) ~ "Copy number", TRUE ~ "Xenium\nprogrammes"),
         label = lab_map[feature]) |> select(label, group, HR, lo, hi, p)
fo <- bind_rows(fo,
  kzc |> mutate(label = paste0("*", feature, "* (low)", ifelse(feature == "ZNF671", "<sup>‡</sup>", "")), group = "Methylation") |>
    select(label, group, HR, lo, hi, p),
  tibble(label = "Dysplasia vs hyperplasia", group = "Histology", HR = exp(g26$fit_h[1, 1]),
         lo = exp(g26$fit_h[1, 1] - 1.96 * g26$fit_h[1, 3]), hi = exp(g26$fit_h[1, 1] + 1.96 * g26$fit_h[1, 3]), p = g26$fit_h[1, 5])) |>
  mutate(group = factor(group, c("Methylation", "Copy number", "Xenium\nprogrammes", "Histology")),
         label = fct_reorder(label, HR), Significance = ifelse(p < 0.05, "P < 0.05", "P ≥ 0.05"), pt = ptxt(p))
f4c <- ggplot(fo, aes(HR, label)) +
  geom_vline(xintercept = 1, colour = MUTED, linewidth = PT(0.5), linetype = "22") +
  geom_linerange(aes(xmin = lo, xmax = hi), colour = INK, linewidth = PT(0.6)) +
  geom_point(aes(fill = Significance), shape = 21, size = 1.4, colour = INK, stroke = PT(0.4)) +
  geom_text(aes(x = 6.2, label = pt), hjust = 1, size = PT(5), colour = INK2) +
  geom_text(data = tibble(group = factor("Methylation", levels(fo$group)), HR = 6.2, label = "italic(P)"),
            aes(x = HR, y = Inf, label = label), parse = TRUE, hjust = 1, vjust = -0.4, size = PT(5.5), inherit.aes = FALSE) +
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c("P < 0.05" = INK, "P ≥ 0.05" = "white"), name = NULL) +
  scale_x_log10(limits = c(0.4, 6.3), breaks = c(0.5, 1, 2, 4)) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Hazard ratio per s.d. (95% CI)", y = NULL) + theme_nat() +
  theme(axis.text.y = element_markdown(), strip.text.y = element_text(angle = 0, hjust = 0, size = 5.3),
        panel.spacing = unit(0.8, "mm"), legend.position = "top", legend.justification = "left",
        legend.margin = margin(0, 0, 1.5, 0, "mm"))
nul <- bind_rows(tibble(null = "Random 50/50 gene sets", lp = -log10(g26$rand_p)),
                 tibble(null = "Retrained on shuffled labels", lp = -log10(g26$perm_p)))
emp <- tibble(null = c("Random 50/50 gene sets", "Retrained on shuffled labels"), e = c(g26$emp_rand, g26$emp_perm))
f4d <- ggplot(nul, aes(lp)) +
  geom_histogram(aes(fill = "Null (1,000 each)"), bins = 34, colour = "white", linewidth = PT(0.3)) +
  geom_vline(aes(xintercept = -log10(g26$obs["p"]), colour = "Observed proxy"), linewidth = PT(1)) +
  geom_text(data = emp, aes(x = -log10(g26$obs["p"]) + 0.12, y = Inf, label = sprintf("Empirical P = %.3f", e)),
            hjust = 0, vjust = 1.5, size = PT(5)) +
  scale_x_continuous(breaks = 0:6) + coord_cartesian(xlim = c(0, 6.3)) +
  scale_fill_manual(values = c("Null (1,000 each)" = "#b7d3f6"), name = NULL) +
  scale_colour_manual(values = c("Observed proxy" = C_O), name = NULL) +
  facet_wrap(~null, ncol = 1, scales = "free_y") + scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "−log10 Cox P in GSE26549", y = "Count") + theme_nat() +
  theme(legend.position = "top", legend.justification = "left", legend.margin = margin(0, 0, 0.5, 0, "mm"))
## e. integrated pilot profile
pz <- ext$pz_p |> select(Sample, pmd = pmd_contrast)
cn_pick <- cnv$arms |> filter(arm_label %in% c("3p", "13q", "12p", "12q")) |>
  select(Sample, arm_label, LPO) |> pivot_wider(names_from = arm_label, values_from = LPO) |> mutate(chr12 = (`12p` + `12q`) / 2)
HIST <- c(P2 = "OSCC", P3 = "PVL", P4 = "AH", P5 = "MSD, progressed")
mx <- tibble(Sample = ORD) |> left_join(ext$sc_p |> select(Sample, Epi, tumor_likeness), by = "Sample") |>
  left_join(pz, by = "Sample") |> mutate(KZFP = rowMeans(kz_p)[Sample]) |>
  left_join(cnv$del9 |> select(Sample, del9 = lpo_9p21_median), by = "Sample") |>
  left_join(cn_pick |> select(Sample, `3p`, `13q`, chr12), by = "Sample") |> mutate(Tissue = tissue_of(Sample))
fd <- tribble(~col, ~label, ~block, ~dir, ~fmt,
  "Epi", "Epith.\nfraction", "Cells", NA, "%.2f",
  "tumor_likeness", "Tumour-\nlikeness", "Methylome", 1, "%.2f",
  "KZFP", "KZFP β", "Methylome", 1, "%.2f",
  "pmd", "PMD", "Methylome", -1, "%.3f",
  "3p", "3p", "Copy number", -1, "%.2f", "del9", "9p21", "Copy number", -1, "%.2f",
  "13q", "13q", "Copy number", -1, "%.2f", "chr12", "Chr 12", "Copy number", 1, "%.2f")
long <- map_dfr(seq_len(nrow(fd)), function(i) {
  f <- fd[i, ]; v <- mx[[f$col]]; ref <- median(v[mx$Tissue == "Adjacent"], na.rm = TRUE)
  s <- if (is.na(f$dir)) rep(NA_real_, length(v)) else f$dir * (v - ref)
  if (!is.na(f$dir)) s <- s / max(abs(s), na.rm = TRUE)
  txt <- sprintf(f$fmt, v); txt[grepl("^-0\\.0+$", txt)] <- sub("-", "", txt[grepl("^-0\\.0+$", txt)])
  tibble(Sample = mx$Sample, label = f$label, block = f$block, fill = s, txt = sub("-", "\u2212", txt, fixed = TRUE))
}) |> mutate(Sample = factor(Sample, rev(ORD)), label = factor(label, fd$label),
             block = factor(block, c("Cells", "Methylome", "Copy number")))
f4e <- ggplot(long, aes(label, Sample)) +
  geom_tile(aes(fill = fill), colour = "white", linewidth = PT(0.8)) +
  geom_text(aes(label = txt), size = PT(4.8), colour = INK) +
  scale_fill_gradient2(low = BLUE, mid = MID, high = RED, midpoint = 0, limits = c(-1, 1), na.value = "#f4f3ef",
                       name = "Toward\ntumour", breaks = c(-1, 0, 1), labels = c("Away", "Adjacent", "Max"), guide = cbar(9)) +
  facet_grid(~block, scales = "free_x", space = "free_x") +
  scale_y_discrete(labels = function(s) ifelse(grepl("Adj|NAT", s), s, paste0(s, " (", HIST[sub("-.*", "", s)], ")")), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) + labs(x = NULL, y = NULL) + theme_nat() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), panel.spacing = unit(0.8, "mm"))
## f. model
st_df <- tibble(x = c(0.55, 1.95, 3.35, 4.75), lab = c("Healthy\nmucosa", "OPL, normal-like\nmethylome", "OPL, tumour-like\nmethylome", "OSCC"),
                fill = c("#eaf2fc", "#fbf4ee", "#f6dccc", "#d8f2e7"))
notes <- tibble(y = c(-0.05, -0.55, -1.05), lab = c(
  "Epigenetic: PMD hypomethylation; 19q13 KZFP promoter hypermethylation and silencing (ZNF471, ZNF529, ZNF671)",
  "Genomic: early low-level loss at 3p, 9p, 13q and 4q; focal MTAP–CDKN2A–CDKN2B deletion (pilot P4, P5)",
  "Clinical: expression proxy predicts transformation beyond histology; microenvironment programmes track grade only"))
f4f <- ggplot() +
  geom_segment(data = tibble(x = st_df$x[-4] + 0.52, xend = st_df$x[-1] - 0.54), aes(x = x, xend = xend, y = 0.75, yend = 0.75),
               colour = INK2, linewidth = PT(0.8), arrow = arrow(length = unit(1.2, "mm"), type = "closed")) +
  geom_rect(data = st_df, aes(xmin = x - 0.5, xmax = x + 0.5, ymin = 0.42, ymax = 1.08, fill = fill), colour = INK, linewidth = PT(0.5)) +
  geom_text(data = st_df, aes(x, 0.75, label = lab), size = PT(5.3), lineheight = 0.9, colour = INK) +
  geom_segment(aes(x = 3.35, xend = 3.35, y = 0.42, yend = 0.15), colour = MUTED, linewidth = PT(0.5)) +
  geom_text(data = notes, aes(x = 2.65, y = y, label = lab), size = PT(5.2), colour = INK) +
  scale_fill_identity() + coord_cartesian(xlim = c(-0.05, 5.35), ylim = c(-1.3, 1.2), expand = FALSE) + theme_void()
fig4 <- f4a / (wrap_elements(full = f4b) + f4c + plot_layout(widths = c(1, 1.05))) / (f4d + f4e + plot_layout(widths = c(0.62, 1.38))) / f4f +
  plot_layout(heights = c(0.27, 1.2, 0.95, 0.42)) + plot_annotation(tag_levels = "a") & tag_theme
save_nat(fig4, "Fig4", 170)

## ================================================================== SUPPLEMENTARY FIGURE 1
lo <- au$leaveout |> mutate(State = paste("OPL", methylome_state))
s1a <- ggplot(lo, aes(tumor_likeness, lo_score, shape = State)) +
  geom_abline(colour = MUTED, linetype = "22", linewidth = PT(0.5)) +
  geom_point(size = 1.2, colour = C_O, stroke = PT(0.6)) +
  scale_shape_manual(values = SHP, name = NULL) +
  annotate("text", x = -0.25, y = 1.75, hjust = 0, vjust = 1, size = PT(5.2), parse = TRUE,
           label = sprintf("rho=='%.2f'*'; same group in 22/22'", au$r_lo)) +
  labs(x = "Tumour-likeness, full signature", y = "Tumour-likeness, without chr19 ZNF,\nPMD-block and open-sea CpGs") + theme_nat() +
  theme(legend.position = "inside", legend.position.inside = c(0.72, 0.18))
sz <- g26$size_tab |> mutate(lab = ifelse(is.na(k_each), note, paste0(k_each, " + ", k_each, " genes")), lab = factor(lab, rev(unique(lab))))
s1b <- ggplot(sz, aes(-log10(p), lab)) +
  geom_vline(xintercept = -log10(0.05), colour = MUTED, linetype = "22", linewidth = PT(0.5)) +
  geom_segment(aes(x = 0, xend = -log10(p), yend = lab), colour = C_N, linewidth = PT(1.2)) +
  geom_point(size = 1.3, colour = C_N) +
  geom_text(aes(label = sprintf("HR %.2f", HR)), nudge_x = 0.2, hjust = 0, size = PT(5)) +
  scale_x_continuous(limits = c(0, 4.8), expand = c(0, 0)) +
  labs(x = "−log10 Cox P (GSE26549)", y = NULL) + theme_nat()
s1c <- ggplot(tibble(hr = au$cox$boot_hr), aes(hr)) +
  geom_histogram(aes(fill = "Bootstrap (1,000)"), bins = 40, colour = "white", linewidth = PT(0.3)) +
  geom_vline(aes(xintercept = au$cox$ci[1, 1], colour = "Observed"), linewidth = PT(1)) +
  geom_vline(xintercept = 1, colour = MUTED, linetype = "22", linewidth = PT(0.5)) +
  scale_fill_manual(values = c("Bootstrap (1,000)" = "#b7d3f6"), name = NULL) +
  scale_colour_manual(values = c(Observed = C_O), name = NULL) +
  annotate("text", x = 2.05, y = Inf, hjust = 0, vjust = 1.3, size = PT(5),
           label = sprintf("95%% interval %.2f\u2013%.2f\nLeave-one-out max P = %.3f\nPH test P = %.2f",
                           quantile(au$cox$boot_hr, 0.025), quantile(au$cox$boot_hr, 0.975), max(au$cox$loo), au$cox$zph["proxy", "p"])) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.6))) + coord_cartesian(xlim = c(0.7, 4.6), clip = "off") +
  labs(x = "Proxy hazard ratio per s.d.", y = "Count") + theme_nat() +
  theme(legend.position = "top", legend.justification = "left", legend.margin = margin(0, 0, 0.5, 0, "mm"))
cc <- as_tibble(g26$comp_check$adj, rownames = "term") |>
  transmute(term = recode(term, proxy = "Methylome-state proxy", stromal = "Stromal score", immune = "Immune score",
                          keratinocyte = "Keratinocyte score", dysplasia = "Dysplasia"),
            HR = `exp(coef)`, lo = exp(coef - 1.96 * `se(coef)`), hi = exp(coef + 1.96 * `se(coef)`), p = `Pr(>|z|)`) |>
  mutate(term = factor(term, rev(term)), Significance = ifelse(p < 0.05, "P < 0.05", "P ≥ 0.05"))
s1d <- ggplot(cc, aes(HR, term)) +
  geom_vline(xintercept = 1, colour = MUTED, linetype = "22", linewidth = PT(0.5)) +
  geom_linerange(aes(xmin = lo, xmax = hi), linewidth = PT(0.6)) +
  geom_point(aes(fill = Significance), shape = 21, size = 1.4, stroke = PT(0.4)) +
  scale_fill_manual(values = c("P < 0.05" = INK, "P ≥ 0.05" = "white"), name = NULL) +
  scale_x_log10() + labs(x = "Hazard ratio (95% CI), joint Cox model", y = NULL) + theme_nat() +
  theme(legend.position = "top", legend.justification = "left", legend.margin = margin(0, 0, 0.5, 0, "mm"))
rr <- bind_rows(rna$res |> transmute(feature = set, AUC = prog_AUC, lo = prog_AUC_lo, hi = prog_AUC_hi, p = prog_wilcox_p),
                tibble(feature = "Methylome_state_proxy", AUC = g26$rep2$prog_AUC, lo = NA, hi = NA, p = g26$rep2$prog_p)) |>
  mutate(feature = recode(feature, Methylome_state_proxy = "Methylome-state proxy", Dysplastic_epithelial = "Dysplastic epithelial",
                          Invasive_epithelial = "Invasive epithelial", Cycling_epithelial = "Cycling epithelial",
                          LRRC15_CAF = "LRRC15+ CAF", FOLR2_TAM = "FOLR2+ TAM", PI16_fibroblast = "PI16+ fibroblast",
                          Meth_silenced_promoters = "Methylation-silenced set", Meth_activated_promoters = "Methylation-activated set",
                          KZFP_chr19_hyper = "KZFP set", p16_9p21 = "CDKN2A/B, MTAP", Treg = "Regulatory T cell"),
         feature = fct_reorder(feature, AUC), Significance = ifelse(p < 0.05, "Nominal P < 0.05", "P ≥ 0.05"))
s1e <- ggplot(rr, aes(AUC, feature)) +
  geom_vline(xintercept = 0.5, colour = MUTED, linewidth = PT(0.5)) +
  geom_vline(aes(xintercept = rna$hist_auc, linetype = "Histology alone"), colour = MUTED, linewidth = PT(0.5)) +
  geom_linerange(aes(xmin = lo, xmax = hi), linewidth = PT(0.6), na.rm = TRUE) +
  geom_point(aes(fill = Significance), shape = 21, size = 1.4, stroke = PT(0.4)) +
  scale_fill_manual(values = c("Nominal P < 0.05" = INK, "P ≥ 0.05" = "white"), name = NULL) +
  scale_linetype_manual(values = c("Histology alone" = "22"), name = NULL) +
  labs(x = "AUC for progression, GSE227919 (10 of 25 progressed)", y = NULL) + theme_nat() +
  theme(legend.position = "top", legend.justification = "left", legend.margin = margin(0, 0, 0.5, 0, "mm"))
pcx <- st$pc_tab |> mutate(arm = factor(arm, arm_order), Direction = ifelse(cancer_minus_control > 0, "Higher in cancer", "Lower in cancer")) |>
  filter(!is.na(arm))
s1f <- ggplot(pcx, aes(arm, cancer_minus_control, fill = Direction)) +
  geom_col(width = 0.75) + geom_hline(yintercept = 0, linewidth = PT(0.5)) +
  geom_text(data = filter(pcx, fdr < 0.05), aes(label = "*", y = cancer_minus_control + sign(cancer_minus_control) * 0.07), size = PT(7)) +
  scale_fill_manual(values = c("Higher in cancer" = RED, "Lower in cancer" = BLUE), name = NULL) +
  labs(x = "Chromosome arm (* FDR < 0.05)", y = "Expression dosage,\ncancer − control (z)") + theme_nat() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 4.8), legend.position = "inside", legend.position.inside = c(0.88, 0.25))
sup1 <- (s1a + s1b + s1c) / (free(s1d) + s1e) / s1f + plot_layout(heights = c(1, 1.05, 0.65)) +
  plot_annotation(tag_levels = "a") & tag_theme
save_nat(sup1, "Supplementary_Fig1", 165)

## ================================================================== SUPPLEMENTARY FIGURE 2
arm_b <- cnv$bins_lpo
chr_len <- arm_b |> group_by(chr) |> summarise(len = max(end)) |> arrange(match(chr, paste0("chr", 1:22))) |>
  mutate(offset = lag(cumsum(as.numeric(len)), default = 0))
gb <- arm_b |> left_join(chr_len, by = "chr") |> mutate(g = mid + offset) |>
  pivot_longer(all_of(ORD), names_to = "Sample", values_to = "lr") |> mutate(Sample = factor(Sample, ORD),
                                                                            Bin = ifelse(terminal, "Terminal 5 Mb (not called)", "50-probe bin"))
gs <- cnv$segs_lpo |> left_join(chr_len, by = c("chrom" = "chr")) |>
  mutate(gs = loc.start + offset, ge = loc.end + offset, Sample = factor(ID, ORD),
         Segment = ifelse(seg.mean > 0.12, "Gain segment", ifelse(seg.mean < -0.12, "Loss segment", "Neutral segment")))
s2a <- ggplot() +
  ggrastr::rasterise(geom_point(data = gb, aes(g, pmax(pmin(lr, 0.8), -0.8), colour = Bin), size = 0.04, alpha = 0.5), dpi = 600, dev = "ragg") +
  scale_colour_manual(values = c("50-probe bin" = "#9a978e", "Terminal 5 Mb (not called)" = "#dcdad3"), name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 1.2, alpha = 1), order = 1)) +
  ggnewscale::new_scale_colour() +
  geom_segment(data = gs, aes(x = gs, xend = ge, y = seg.mean, yend = seg.mean, colour = Segment), linewidth = PT(1)) +
  geom_vline(xintercept = chr_len$offset, colour = HAIR, linewidth = PT(0.3)) +
  scale_colour_manual(values = c("Gain segment" = RED, "Loss segment" = BLUE, "Neutral segment" = INK), name = NULL) +
  scale_x_continuous(breaks = chr_len$offset + chr_len$len / 2, labels = sub("chr", "", chr_len$chr), expand = c(0.003, 0)) +
  coord_cartesian(ylim = c(-0.8, 0.8)) + facet_wrap(~Sample, ncol = 1, strip.position = "right") +
  labs(x = "Chromosome", y = "log2 ratio (vs other patients' adjacent tissue)") + theme_nat() +
  theme(axis.text.x = element_text(size = 4.5), strip.text.y = element_text(angle = 0, hjust = 0), legend.position = "top",
        panel.spacing = unit(0.5, "mm"))
conc <- as_tibble(cnv$conc, rownames = "a") |> pivot_longer(-a, names_to = "b", values_to = "c") |>
  mutate(a = factor(a, rev(ORD)), b = factor(b, ORD))
s2b <- ggplot(conc, aes(b, a, fill = c)) + geom_tile(colour = "white", linewidth = PT(0.6)) +
  geom_text(aes(label = sprintf("%.2f", c), colour = c > 0.75), size = PT(4.5)) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = INK), guide = "none") +
  scale_fill_gradientn(colours = SEQ, limits = c(0.3, 1), name = "Genotype\nconcordance", guide = cbar(9)) +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + labs(x = NULL, y = NULL) + theme_nat() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.line = element_blank(), axis.ticks = element_blank())
dl <- cnv$del9 |> mutate(Sample = factor(Sample, ORD))
s2c <- ggplot(dl, aes(Sample, fail_ratio, fill = lpo_9p21_median)) +
  geom_col(width = 0.7, colour = INK, linewidth = PT(0.3)) + geom_hline(yintercept = 1, linetype = "22", colour = MUTED, linewidth = PT(0.5)) +
  scale_fill_gradient2(low = BLUE, mid = MID, high = RED, midpoint = 0, name = "9p21\nlog2 ratio", guide = cbar(9)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "9p21 probe-failure rate\n(× genome-wide)") + theme_nat() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
sup2 <- s2a / (s2b + s2c + plot_layout(widths = c(1.1, 1))) + plot_layout(heights = c(2.3, 1)) +
  plot_annotation(tag_levels = "a") & tag_theme
save_nat(sup2, "Supplementary_Fig2", 170)
cat("done\n")
