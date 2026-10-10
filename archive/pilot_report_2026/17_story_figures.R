## Multi-panel figures for the central story. Reads only saved results from
## R/10-R/16; writes PDF (vector) + PNG (300 dpi) to output/story_figures/.
##
## Colour roles (validated reference palette, light surface):
##   tissue identity  Normal #2a78d6 | OPL/lesion #eb6834 | Tumour #1baf7a  (slots 1-3)
##   composition      Epithelial #eda100 | Fibroblast #e87ba4 | Immune #008300 (slots 4-6)
##   magnitude        single-hue blue ramp; polarity  blue - grey - red
##
## Run from analysis/:  Rscript R/17_story_figures.R

suppressPackageStartupMessages({
  library(tidyverse); library(patchwork); library(ggrepel); library(survival); library(scales)
})
select <- dplyr::select; filter <- dplyr::filter
EXT <- "output/extended"
OUT <- "output/story_figures"; dir.create(OUT, showWarnings = FALSE)
rd <- function(f) readRDS(file.path(EXT, f))

## ------------------------------------------------------------------ style
INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"; GRID <- "#e1e0d9"; AXIS <- "#c3c2b7"
C_N <- "#2a78d6"; C_O <- "#eb6834"; C_T <- "#1baf7a"
C_EPI <- "#eda100"; C_FIB <- "#e87ba4"; C_IC <- "#008300"
BLUE <- "#2a78d6"; RED <- "#e34948"; MID <- "#f0efec"
SEQ <- c("#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b")
grp_cols <- c(Normal = C_N, OPL = C_O, Lesion = C_O, Tumor = C_T)
theme_fig <- function(base = 7.5) {
  theme_minimal(base_size = base, base_family = "Helvetica") +
    theme(panel.grid.major = element_line(colour = GRID, linewidth = 0.25),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = AXIS, linewidth = 0.3),
          axis.ticks = element_line(colour = AXIS, linewidth = 0.3),
          axis.text = element_text(colour = INK2), axis.title = element_text(colour = INK),
          plot.title = element_text(face = "bold", size = base + 0.5, colour = INK),
          plot.subtitle = element_text(colour = INK2, size = base - 0.5),
          legend.text = element_text(colour = INK2), legend.title = element_text(colour = INK),
          strip.text = element_text(colour = INK, face = "bold"),
          plot.tag = element_text(face = "bold", size = base + 3))
}
save_fig <- function(p, name, w, h) {
  ## plain pdf(): the Cairo PDF device is unavailable here, so figure text is kept ASCII/Latin-1
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, device = pdf, useDingbats = FALSE)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300, device = ragg::agg_png, bg = "white")
}
ORD <- c("P2-NAT", "P2-Tumor", "P3-Adj", "P3-Lesion", "P4-Adj", "P4-Lesion", "P5-Adj", "P5-Lesion")
own_adj <- c(`P2-Tumor` = "P2-NAT", `P3-Lesion` = "P3-Adj", `P4-Lesion` = "P4-Adj", `P5-Lesion` = "P5-Adj")
tissue_of <- function(s) ifelse(grepl("Adj|NAT", s), "Adjacent", ifelse(grepl("Tumor", s), "Tumor", "Lesion"))
HIST <- c(P2 = "T2 OSCC", P3 = "PVL (high grade)", P4 = "Atypical hyperplasia (low grade)",
          P5 = "Mod-severe dysplasia (progressed)")

ext <- rd("external_methylation.rds"); cnv <- rd("cnv_objects.rds"); st <- rd("story_analyses.rds")
mum <- rd("mumbai_crosscheck.rds"); g26 <- rd("gse26549_progression.rds"); k3 <- rd("gse234379_kzfp.rds")
rna <- rd("external_rnaseq.rds"); p <- readRDS("output/processed.rds")

## ================================================================== FIGURE 1
ss <- p$sample_sheet
short <- function(sn) { pt <- sub("Patient ", "P", sub("_.*", "", sn))
  paste0(pt, "-", ifelse(grepl("tumor", tolower(sn)), "Tumor", ifelse(grepl("lesion", tolower(sn)), "Lesion",
                                                                        ifelse(pt == "P2", "NAT", "Adj")))) }
pc <- tibble(Sample = short(rownames(p$pca$x)), PC1 = p$pca$x[, 1], PC2 = p$pca$x[, 2]) |>
  left_join(ext$decon |> select(Sample, Epi = Epi_1st, Fib = Fib_1st, IC = IC_1st), by = "Sample") |>
  mutate(Tissue = tissue_of(Sample))
pct <- round(100 * p$pca$sdev^2 / sum(p$pca$sdev^2), 1)
f1a <- ggplot(pc, aes(PC1, PC2)) +
  geom_point(aes(fill = Epi, shape = Tissue), size = 3, colour = "white", stroke = 0.6) +
  geom_text_repel(aes(label = Sample), size = 2.3, colour = INK2, seed = 1, box.padding = 0.35) +
  scale_fill_gradientn(colours = SEQ, limits = c(0, 1), name = "Epithelial fraction",
                       guide = guide_colourbar(barheight = unit(2.5, "mm"), barwidth = unit(22, "mm"))) +
  scale_shape_manual(values = c(Adjacent = 21, Lesion = 24, Tumor = 22), name = NULL) +
  guides(shape = guide_legend(override.aes = list(fill = INK2, size = 2.2))) +
  labs(x = sprintf("PC1 (%s%%)", pct[1]), y = sprintf("PC2 (%s%%)", pct[2]),
       title = "Methylation PCA (manuscript Fig. 3A)") + theme_fig() +
  theme(legend.position = "bottom", legend.box = "horizontal")
r1 <- cor(pc$PC1, pc$Epi)
f1b <- ggplot(pc, aes(Epi, PC1)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = MUTED, linewidth = 0.4) +
  geom_point(aes(shape = Tissue), size = 2.6, fill = INK, colour = "white", stroke = 0.5) +
  geom_text_repel(aes(label = Sample), size = 2.3, colour = INK2, seed = 2) +
  scale_shape_manual(values = c(Adjacent = 21, Lesion = 24, Tumor = 22), guide = "none") +
  annotate("text", x = 0.05, y = max(pc$PC1), hjust = 0, vjust = 1, size = 2.8, colour = INK,
           label = sprintf("r = %.2f", r1)) +
  labs(x = "Epithelial fraction (EpiDISH)", y = "PC1", title = "PC1 is tissue composition") + theme_fig()
comp <- pc |> select(Sample, Epithelial = Epi, Fibroblast = Fib, Immune = IC) |>
  pivot_longer(-Sample, names_to = "Cell", values_to = "Fraction") |>
  mutate(Sample = factor(Sample, ORD), Cell = factor(Cell, c("Immune", "Fibroblast", "Epithelial")))
f1c <- ggplot(comp, aes(Sample, Fraction, fill = Cell)) +
  geom_col(width = 0.75, colour = "white", linewidth = 0.4) +
  scale_fill_manual(values = c(Epithelial = C_EPI, Fibroblast = C_FIB, Immune = C_IC), name = NULL) +
  scale_y_continuous(labels = percent, expand = c(0, 0)) +
  labs(x = NULL, y = "Estimated fraction", title = "Cell composition differs >8-fold in epithelium") +
  theme_fig() + theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "top")
## naive vs adjusted ranking of the three premalignant lesions
tum <- pc |> filter(Sample == "P2-Tumor")
rk <- pc |> filter(Sample %in% c("P3-Lesion", "P4-Lesion", "P5-Lesion")) |>
  mutate(dist = sqrt((PC1 - tum$PC1)^2 + (PC2 - tum$PC2)^2)) |>
  left_join(ext$sc_p |> select(Sample, delta = delta_vs_own_adj), by = "Sample") |>
  mutate(naive_rank = rank(dist), adj_rank = rank(-delta),
         lab = c(`P3-Lesion` = "P3  PVL", `P4-Lesion` = "P4  atypical hyperplasia",
                 `P5-Lesion` = "P5  dysplasia, progressed")[Sample])
sl <- rk |> select(Sample, lab, naive_rank, adj_rank, dist, delta) |>
  pivot_longer(c(naive_rank, adj_rank), names_to = "view", values_to = "rank") |>
  mutate(x = ifelse(view == "naive_rank", 1, 2),
         value = ifelse(view == "naive_rank", sprintf("%.0f", dist), sprintf("%+.2f", delta)))
f1d <- ggplot(sl, aes(x, rank, group = Sample)) +
  geom_line(aes(colour = Sample == "P5-Lesion"), linewidth = 0.8) +
  geom_point(aes(colour = Sample == "P5-Lesion"), size = 2.6) +
  geom_text(data = filter(sl, x == 1), aes(label = value), hjust = 1.6, size = 2.4, colour = INK2) +
  geom_text(data = filter(sl, x == 2), aes(label = value), hjust = -0.35, size = 2.4, colour = INK2) +
  geom_text(data = filter(sl, x == 2), aes(x = 2.62, label = lab), hjust = 0, size = 2.4, colour = INK) +
  scale_colour_manual(values = c(`TRUE` = C_O, `FALSE` = MUTED), guide = "none") +
  scale_y_reverse(breaks = 1:3, labels = c("1 (closest)", "2", "3 (farthest)")) +
  scale_x_continuous(breaks = 1:2, labels = c("Naive: PCA distance\nto tumour", "Adjusted: shift\ntoward tumour"),
                     limits = c(0.45, 4.0)) +
  labs(x = NULL, y = "Rank among lesions", title = "Composition adjustment re-orders the lesions") +
  theme_fig() + theme(panel.grid.major.x = element_blank())
fig1 <- (f1a + f1b) / (f1c + f1d) + plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A", title = "Figure 1. Bulk methylation of oral biopsies is dominated by tissue composition",
                  theme = theme(plot.title = element_text(family = "Helvetica", face = "bold", size = 10)))
save_fig(fig1, "Figure1_composition", 7.2, 6.2)

## ================================================================== FIGURE 2
sg <- ext$sc_g |> left_join(st$opl |> select(Sample, methylome_state), by = "Sample") |>
  mutate(Group = as.character(Group),
         row = recode(Group, Normal = "GSE204943 healthy mucosa", OPL = "GSE204943 OPL", Tumor = "GSE204943 OSCC"),
         State = ifelse(Group == "OPL", methylome_state, Group))
s3 <- k3$s3 |> mutate(row = ifelse(group == "Tumor", "GSE234379 OSCC", "GSE234379 tumour-adjacent"),
                      Group = group, State = group)
pil <- ext$sc_p |> mutate(row = "Pilot (this study)", Tissue = tissue_of(Sample),
                          Patient = sub("-.*", "", Sample))
rows <- c("Pilot (this study)", "GSE234379 OSCC", "GSE234379 tumour-adjacent", "GSE204943 OSCC",
          "GSE204943 OPL", "GSE204943 healthy mucosa")
pts <- bind_rows(sg |> select(row, tumor_likeness, Group, State), s3 |> select(row, tumor_likeness, Group, State)) |>
  mutate(y = match(row, rows))
pil_y <- tibble(Patient = c("P2", "P3", "P4", "P5"), dy = c(0.27, 0.09, -0.09, -0.27))
pil2 <- pil |> left_join(pil_y, by = "Patient") |> mutate(y = match(row, rows) + dy)
arr <- tibble(lesion = names(own_adj), adj = unname(own_adj)) |>
  left_join(pil2 |> select(lesion = Sample, xend = tumor_likeness, y), by = "lesion") |>
  left_join(pil2 |> select(adj = Sample, x = tumor_likeness), by = "adj")
set.seed(3)
f2a <- ggplot() +
  annotate("rect", xmin = st$valley, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#fbf3ee", alpha = 0.8) +
  geom_vline(xintercept = st$valley, linetype = "22", colour = MUTED, linewidth = 0.35) +
  geom_jitter(data = pts, aes(tumor_likeness, y, colour = Group, shape = State), height = 0.18, width = 0, size = 1.1, alpha = 0.75) +
  geom_segment(data = arr, aes(x = x, xend = xend, y = y, yend = y), colour = INK2, linewidth = 0.35,
               arrow = arrow(length = unit(1.4, "mm"), type = "closed")) +
  geom_point(data = pil2, aes(tumor_likeness, y, fill = Tissue == "Adjacent"), shape = 21, size = 2, colour = INK, stroke = 0.4) +
  geom_text(data = filter(pil2, Tissue != "Adjacent"), aes(tumor_likeness, y, label = Patient), hjust = -0.5, size = 2.3, colour = INK) +
  scale_fill_manual(values = c(`TRUE` = "white", `FALSE` = INK), guide = "none") +
  scale_colour_manual(values = grp_cols, guide = "none") +
  scale_shape_manual(values = c(Normal = 16, Tumor = 16, `normal-like` = 1, `tumour-like` = 17), guide = "none") +
  annotate("text", x = st$valley + 0.03, y = 6.45, label = "tumour-like state", hjust = 0, size = 2.4, colour = C_O) +
  scale_y_continuous(breaks = seq_along(rows), labels = rows) +
  coord_cartesian(xlim = c(-0.5, 2.15), clip = "off") +
  labs(x = "Tumour-likeness (0 = healthy mucosa, 1 = OSCC; composition-adjusted)", y = NULL,
       title = sprintf("OPLs separate into normal-like and tumour-like groups; signature transfers to a third cohort (AUC %.3f)", k3$auc3)) +
  theme_fig()
f2b <- ggplot() +
  geom_point(data = sg, aes(Epi, tumor_likeness, colour = Group, shape = State), size = 1.3, alpha = 0.8) +
  geom_point(data = pil, aes(Epi, tumor_likeness, fill = Tissue == "Adjacent"), shape = 21, size = 2.2, colour = INK) +
  geom_text_repel(data = pil, aes(Epi, tumor_likeness, label = Sample), size = 2.1, colour = INK2, seed = 4) +
  scale_fill_manual(values = c(`TRUE` = "white", `FALSE` = INK), guide = "none") +
  scale_colour_manual(values = grp_cols, name = NULL) +
  scale_shape_manual(values = c(Normal = 16, Tumor = 16, `normal-like` = 1, `tumour-like` = 17), name = "OPL state",
                     breaks = c("normal-like", "tumour-like")) +
  labs(x = "Epithelial fraction", y = "Tumour-likeness", title = "Not explained by epithelial content") +
  theme_fig() + theme(legend.position = "right")
gk <- k3$gk |> left_join(ext$pz_g |> select(Sample, pmd_contrast), by = "Sample") |>
  mutate(KZFP = rowMeans(across(c(ZNF471, ZNF382, ZNF529, ZNF582, ZNF671))),
         State = ifelse(Group == "OPL", sub("OPL ", "", state), Group))
pil_k <- ext$pz_p |> select(Sample, pmd_contrast) |>
  left_join(tibble(Sample = ORD, KZFP = NA_real_), by = "Sample")
## pilot KZFP means per sample, recomputed from the same promoter CpGs
kz_cpg <- local({
  a <- readRDS(path.expand("~/Desktop/OPMD_methyl/GSE204943_data/GSE204943/crucial_objects/ann_auto.rds"))
  a <- a[intersect(rownames(a), rownames(p$betas)), ]
  keep <- grepl("TSS200|TSS1500|1stExon", a$UCSC_RefGene_Group) & a$Relation_to_Island == "Island"
  gl <- strsplit(a$UCSC_RefGene_Name, ";")
  unlist(lapply(c("ZNF471", "ZNF382", "ZNF529", "ZNF582", "ZNF671"), function(g)
    list(rownames(a)[keep & vapply(gl, function(x) g %in% x, logical(1))])), recursive = FALSE)
})
kz_p <- sapply(kz_cpg, function(c) colMeans(p$betas[c, , drop = FALSE]))
pil_k$KZFP <- rowMeans(kz_p)[match(pil_k$Sample, short(rownames(kz_p)))]
pil_k$Tissue <- tissue_of(pil_k$Sample)
arr_k <- tibble(lesion = names(own_adj), adj = unname(own_adj)) |>
  left_join(pil_k |> select(lesion = Sample, xend = pmd_contrast, yend = KZFP), by = "lesion") |>
  left_join(pil_k |> select(adj = Sample, x = pmd_contrast, y = KZFP), by = "adj")
rho_o <- with(filter(gk, Group == "OPL"), cor(pmd_contrast, KZFP, method = "spearman"))
f2c <- ggplot() +
  geom_point(data = gk, aes(-pmd_contrast, KZFP, colour = Group, shape = State), size = 1.3, alpha = 0.8) +
  geom_segment(data = arr_k, aes(-x, y, xend = -xend, yend = yend), colour = INK2, linewidth = 0.35,
               arrow = arrow(length = unit(1.4, "mm"), type = "closed")) +
  geom_point(data = pil_k, aes(-pmd_contrast, KZFP, fill = Tissue == "Adjacent"), shape = 21, size = 2.2, colour = INK) +
  geom_text_repel(data = filter(pil_k, Tissue != "Adjacent"), aes(-pmd_contrast, KZFP, label = Sample),
                  size = 2.1, colour = INK2, seed = 5) +
  scale_fill_manual(values = c(`TRUE` = "white", `FALSE` = INK), guide = "none") +
  scale_colour_manual(values = grp_cols, guide = "none") +
  scale_shape_manual(values = c(Normal = 16, Tumor = 16, `normal-like` = 1, `tumour-like` = 17), guide = "none") +
  annotate("text", x = -0.06, y = 0.62, hjust = 0, size = 2.5, colour = INK, label = sprintf("OPL: rho = %.2f", -rho_o)) +
  labs(x = "PMD erosion (outside minus inside PMD blocks)", y = "chr19q13 KZFP promoter methylation\n(mean of 5 genes)",
       title = "PMD erosion and KZFP methylation co-occur") + theme_fig()
## KZFP heat: medians by group (GSE204943, GSE234379) + pilot deltas
kz_genes <- c("ZNF471", "ZNF382", "ZNF529", "ZNF582", "ZNF671")
hm <- bind_rows(
  k3$gk |> mutate(col = recode(state, Normal = "Healthy", `OPL normal-like` = "OPL\nnormal-like",
                               `OPL tumour-like` = "OPL\ntumour-like", Tumor = "OSCC"), cohort = "GSE204943 (India)"),
  k3$s3 |> mutate(col = ifelse(group == "Tumor", "OSCC ", "Adjacent"), cohort = "GSE234379 (Brazil)")) |>
  select(cohort, col, all_of(kz_genes)) |> pivot_longer(all_of(kz_genes), names_to = "gene", values_to = "beta") |>
  group_by(cohort, col, gene) |> summarise(beta = median(beta, na.rm = TRUE), .groups = "drop") |>
  mutate(col = factor(col, c("Healthy", "OPL\nnormal-like", "OPL\ntumour-like", "OSCC", "Adjacent", "OSCC ")),
         gene = factor(gene, rev(kz_genes)))
f2d1 <- ggplot(hm, aes(col, gene, fill = beta)) + geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", beta), colour = beta > 0.35), size = 2.2) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = INK), guide = "none") +
  scale_fill_gradientn(colours = SEQ, limits = c(0, 0.6), name = "Median beta") +
  facet_grid(~cohort, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = NULL, title = "KZFP promoter hypermethylation marks the tumour-like OPLs") + theme_fig() +
  theme(panel.grid = element_blank())
pk <- k3$pil |> select(gene, P2 = d_P2_tumor, P3 = d_P3, P4 = d_P4, P5 = d_P5_progressor) |>
  pivot_longer(-gene, names_to = "Patient", values_to = "delta") |>
  mutate(gene = factor(gene, rev(kz_genes)),
         Patient = factor(Patient, c("P2", "P3", "P4", "P5"), c("P2\nOSCC", "P3\nPVL", "P4\nAH", "P5\nprogressor")))
f2d2 <- ggplot(pk, aes(Patient, gene, fill = delta)) + geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%+.2f", delta)), size = 2.2, colour = INK) +
  scale_fill_gradient2(low = BLUE, mid = MID, high = RED, midpoint = 0, limits = c(-0.6, 0.6), name = "delta-beta vs\nown adjacent") +
  labs(x = NULL, y = NULL, title = "Pilot pairs") + theme_fig() +
  theme(panel.grid = element_blank(), axis.text.y = element_blank())
fig2 <- f2a / (f2b + f2c) / (f2d1 + f2d2 + plot_layout(widths = c(2.4, 1))) +
  plot_layout(heights = c(1.15, 1, 0.8)) +
  plot_annotation(tag_levels = list(c("A", "B", "C", "D", "")),
                  title = "Figure 2. A composition-adjusted external continuum separates normal-like from tumour-like OPLs",
                  theme = theme(plot.title = element_text(family = "Helvetica", face = "bold", size = 10)))
save_fig(fig2, "Figure2_methylome_groups", 7.4, 9)

## ================================================================== FIGURE 3
xp <- mum$xp |> mutate(Group = ifelse(Group == "Tumor", "Tumor", "OPL"))
f3a <- ggplot(xp, aes(conumee, aCGH)) +
  geom_abline(slope = 1, intercept = 0, colour = MUTED, linetype = "22", linewidth = 0.3) +
  geom_point(aes(colour = Group), size = 0.7, alpha = 0.55) +
  scale_colour_manual(values = grp_cols, name = NULL) +
  annotate("text", x = 0.6, y = -0.6, hjust = 1, vjust = 0, size = 2.4, colour = INK,
           label = sprintf("arm-level r = %.2f\nper-lesion median r = %.2f\n(shuffled null %.2f)\n%d lesions, same patients",
                           mum$xp_r, median(mum$per_sample$r_own), median(mum$null_r), n_distinct(xp$id))) +
  coord_equal(xlim = c(-0.6, 0.6), ylim = c(-0.6, 0.6)) +
  labs(x = "Copy number from methylation array (log2)", y = "Copy number from array CGH (log2)",
       title = "Methylation arrays measure copy number") + theme_fig() +
  theme(legend.position = "inside", legend.position.inside = c(0.15, 0.88))
risk_arms <- c("3p", "9p", "13q", "4q", "8p", "11q", "17p")
ob <- mum$opl_c |> filter(!is.na(methylome_state)) |> select(id, methylome_state, all_of(risk_arms)) |>
  pivot_longer(all_of(risk_arms), names_to = "arm", values_to = "lr") |>
  left_join(mum$opl_tests |> select(arm = feature, wilcox_p), by = "arm") |>
  mutate(arm = factor(arm, risk_arms), lab = sprintf("%s\np = %s", arm, formatC(signif(wilcox_p, 2), format = "fg", digits = 2)),
         lab = factor(lab, unique(lab[order(arm)])))
f3b <- ggplot(ob, aes(methylome_state, lr)) +
  geom_hline(yintercept = 0, colour = AXIS, linewidth = 0.3) +
  geom_jitter(aes(shape = methylome_state), width = 0.15, size = 1.3, colour = C_O) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, linewidth = 0.25, colour = INK) +
  scale_shape_manual(values = c(`normal-like` = 1, `tumour-like` = 17), name = "OPL methylome state (GSE204943)") +
  facet_wrap(~lab, nrow = 1) +
  labs(x = NULL, y = "Array CGH arm log2 ratio (GSE85514)",
       title = "Tumour-like OPLs carry early losses at the oral LOH arms 3p, 9p, 13q, 4q") +
  theme_fig() + theme(panel.grid.major.x = element_blank(), axis.text.x = element_blank(),
                      legend.position = "bottom")
arm_order <- paste0(rep(1:22, each = 2), c("p", "q"))
cn <- cnv$arms |> mutate(Sample = factor(Sample, rev(ORD)), arm_label = factor(arm_label, arm_order)) |>
  filter(!is.na(arm_label))
f3c <- ggplot(cn, aes(arm_label, Sample, fill = pmax(pmin(LPO, 0.25), -0.25))) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_point(data = filter(cn, call != ""), shape = 16, size = 0.7, colour = INK) +
  scale_fill_gradient2(low = BLUE, mid = MID, high = RED, midpoint = 0, limits = c(-0.25, 0.25), name = "log2 ratio") +
  labs(x = NULL, y = NULL, title = "Pilot arm-level copy number from the same arrays (dot = called)") +
  theme_fig() + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, size = 5.5))
z9 <- cnv$bins_lpo |> filter(chr == "chr9", mid > 20e6, mid < 24e6) |>
  pivot_longer(all_of(ORD), names_to = "Sample", values_to = "lr") |>
  mutate(hl = Sample %in% c("P3-Adj", "P4-Lesion", "P5-Lesion"))
genes9 <- tibble(g = c("MTAP", "CDKN2A", "CDKN2B"), s = c(21.80, 21.967, 22.002), e = c(21.94, 21.995, 22.009))
f3d <- ggplot(z9, aes(mid / 1e6, lr, group = Sample)) +
  annotate("rect", xmin = 21.80, xmax = 22.01, ymin = -Inf, ymax = Inf, fill = "#fdf1d6") +
  geom_hline(yintercept = 0, colour = AXIS, linewidth = 0.3) +
  geom_line(data = filter(z9, !hl), colour = GRID, linewidth = 0.5) +
  geom_line(data = filter(z9, hl), aes(colour = Sample), linewidth = 0.8) +
  annotate("text", x = 21.905, y = 0.5, label = "MTAP-CDKN2A-CDKN2B", size = 2.1, colour = INK2) +
  scale_colour_manual(values = c(`P5-Lesion` = C_O, `P4-Lesion` = C_T, `P3-Adj` = C_N), name = NULL) +
  coord_cartesian(ylim = c(-0.7, 0.55)) +
  labs(x = "chr9 (Mb, hg38)", y = "log2 ratio", title = "Focal 9p21 deletion (P3-Adj, P4-Lesion, P5-Lesion)") +
  theme_fig() + theme(legend.position = "inside", legend.position.inside = c(0.2, 0.2))
cd <- ext$gene_tab |> filter(gene == "CDKN2A") |> select(all_of(ORD)) |> pivot_longer(everything(), names_to = "Sample", values_to = "prom") |>
  left_join(cnv$del9 |> select(Sample, lr = lpo_9p21_median), by = "Sample") |> mutate(Tissue = tissue_of(Sample))
f3e <- ggplot(cd, aes(lr, prom)) +
  geom_point(aes(shape = Tissue), size = 2.2, fill = INK, colour = "white") +
  geom_text_repel(aes(label = Sample), size = 2, colour = INK2, seed = 6) +
  scale_shape_manual(values = c(Adjacent = 21, Lesion = 24, Tumor = 22), guide = "none") +
  coord_cartesian(ylim = c(0, 0.5)) +
  labs(x = "9p21 log2 ratio", y = "CDKN2A promoter beta", title = "Deleted, not methylated") + theme_fig()
fig3 <- (f3a + f3b + plot_layout(widths = c(1, 2))) / f3c / (f3d + f3e + plot_layout(widths = c(1.6, 1))) +
  plot_layout(heights = c(1.1, 0.75, 1)) +
  plot_annotation(tag_levels = "A",
                  title = "Figure 3. The same methylation array reads copy number and LOH-arm loss",
                  theme = theme(plot.title = element_text(family = "Helvetica", face = "bold", size = 10)))
save_fig(fig3, "Figure3_copy_number", 7.4, 8.4)

## ================================================================== FIGURE 4
## A: design schematic
boxes <- tibble(x = c(1, 2.15, 3.3), y = 1,
  lab = c("GSE204943 (Mumbai)\nEPIC methylation\n22 OPL: methylome state",
          "GSE85195 (same OPLs)\nexpression, n = 15\n50 up / 50 down gene proxy",
          "GSE26549 (MD Anderson)\n86 OPL, 35 cancers, ~6 y\nproxy applied unchanged"))
f4a <- ggplot(boxes) +
  geom_segment(data = tibble(x = c(1.45, 2.6), xend = c(1.7, 2.85)), aes(x = x, xend = xend, y = 1, yend = 1),
               arrow = arrow(length = unit(1.6, "mm"), type = "closed"), colour = INK2) +
  geom_label(aes(x, y, label = lab), size = 2.3, label.size = 0.25, label.padding = unit(1.6, "mm"),
             fill = c("#fbf3ee", "#fbf3ee", "#eaf2fc"), colour = INK, family = "Helvetica") +
  xlim(0.45, 3.85) + ylim(0.8, 1.2) + theme_void() +
  labs(title = "Train on methylation-defined state in one cohort; test on transformation in another") +
  theme(plot.title = element_text(face = "bold", size = 8, family = "Helvetica", colour = INK))
km <- g26$km |> filter(feature == "Methylome_state_proxy")
km0 <- km |> group_by(strata) |> summarise(time = 0, surv = 1, .groups = "drop")
km <- bind_rows(km0, km |> select(strata, time, surv)) |> arrange(strata, time)
ct <- g26$cox_tab |> filter(feature == "Methylome_state_proxy")
f4b <- ggplot(km, aes(time, surv, colour = strata)) + geom_step(linewidth = 0.7) +
  scale_colour_manual(values = c(High = C_O, Low = C_N), name = "Proxy (median split)") +
  scale_y_continuous(labels = percent, limits = c(0, 1)) +
  annotate("text", x = 0.3, y = 0.12, hjust = 0, size = 2.5, colour = INK,
           label = sprintf("HR %.2f per SD, p = %.3f\nadjusted for histology p = %.3f\nlog-rank (median split) p = %.3f\nC-index %.2f vs histology %.2f",
                           ct$HR_per_SD, ct$p, ct$p_adj_histology, g26$km$logrank_p[g26$km$feature == "Methylome_state_proxy"][1],
                           ct$C_with_histology, g26$c_h)) +
  labs(x = "Years", y = "Oral cancer-free", title = "Methylome-state proxy predicts transformation") +
  theme_fig() + theme(legend.position = "inside", legend.position.inside = c(0.75, 0.85))
## forest: HR per SD with Wald 95% CI from HR and p; KZFP genes shown per SD LOWER expression
wald_ci <- function(hr, p) { se <- abs(log(hr)) / qnorm(1 - p / 2); exp(log(hr) + c(-1.96, 1.96) * se) }
lab_map <- c(Methylome_state_proxy = "Methylome-state proxy", KZFP_chr19_hyper = "KZFP set (ZNF471/382/529), low",
             loss_3p = "3p dosage loss", loss_9p = "9p dosage loss", loss_core_3p9p = "3p + 9p dosage loss",
             Dysplastic_epithelial = "Dysplastic epithelial", Invasive_epithelial = "Invasive epithelial",
             Cycling_epithelial = "Cycling epithelial", LRRC15_CAF = "LRRC15+ CAF", FOLR2_TAM = "FOLR2+ TAM",
             Treg = "Treg", mregDC = "mregDC")
fo <- g26$cox_tab |> filter(feature %in% names(lab_map)) |>
  mutate(HR = ifelse(feature == "KZFP_chr19_hyper", 1 / HR_per_SD, HR_per_SD),
         group = case_when(feature %in% c("Methylome_state_proxy", "KZFP_chr19_hyper") ~ "Methylation-derived",
                           grepl("loss", feature) ~ "Copy number (dosage)", TRUE ~ "Manuscript programs"),
         label = lab_map[feature]) |> select(label, group, HR, p, p_adj = p_adj_histology)
fo <- bind_rows(fo, g26$kz |> filter(gene %in% c("ZNF471", "ZNF582", "ZNF529", "ZNF671")) |>
  transmute(label = paste0(gene, " low", ifelse(gene == "ZNF671", " *", "")), group = "Methylation-derived",
            HR = 1 / HR_per_SD_expression, p, p_adj = NA_real_),
  tibble(label = "Dysplasia vs hyperplasia", group = "Histology", HR = exp(g26$fit_h[1, 1]), p = g26$fit_h[1, 5], p_adj = NA))
fo <- fo |> rowwise() |> mutate(lo = wald_ci(HR, p)[1], hi = wald_ci(HR, p)[2]) |> ungroup() |>
  mutate(group = factor(group, c("Methylation-derived", "Copy number (dosage)", "Manuscript programs", "Histology")),
         label = fct_reorder(label, HR), sig = p < 0.05)
f4c <- ggplot(fo, aes(HR, label)) +
  geom_vline(xintercept = 1, colour = AXIS, linewidth = 0.35) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, colour = INK2, linewidth = 0.4) +
  geom_point(aes(fill = sig), shape = 21, size = 2, colour = INK) +
  geom_text(aes(x = 4.4, label = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))), size = 2.1, colour = INK2, hjust = 1) +
  scale_fill_manual(values = c(`TRUE` = INK, `FALSE` = "white"), guide = "none") +
  scale_x_log10(limits = c(0.4, 4.5), breaks = c(0.5, 1, 2, 4)) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Hazard ratio per SD (95% CI); KZFP genes per SD lower expression", y = NULL,
       title = "GSE26549, oral cancer-free survival (p at right; * post hoc)") +
  theme_fig() + theme(strip.text.y = element_text(angle = 0, hjust = 0, size = 6))
nul <- bind_rows(tibble(null = "1,000 random 50/50 gene sets", lp = -log10(g26$rand_p)),
                 tibble(null = "1,000 retrainings on shuffled\nmethylome labels", lp = -log10(g26$perm_p)))
emp <- tibble(null = unique(nul$null), e = c(g26$emp_rand, g26$emp_perm))
f4d <- ggplot(nul, aes(lp)) + geom_histogram(bins = 40, fill = "#b7d3f6", colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = -log10(g26$obs["p"]), colour = C_O, linewidth = 0.7) +
  geom_text(data = emp, aes(x = -log10(g26$obs["p"]), y = Inf, label = sprintf("observed\nempirical p = %.3f", e)),
            hjust = -0.05, vjust = 1.2, size = 2.2, colour = INK) +
  facet_wrap(~null, ncol = 1, scales = "free_y") +
  labs(x = "-log10 Cox p in GSE26549", y = "count", title = "Not a random-signature artefact") + theme_fig()
me <- mum$me_tab |> filter(gene %in% c("ZNF471", "ZNF382", "ZNF529", "ZNF582", "ZNF671", "NELL1", "MGMT")) |>
  mutate(gene = fct_reorder(gene, rho))
bg_rho <- median(mum$me_tab$rho, na.rm = TRUE)
f4e <- ggplot(me, aes(rho, gene)) +
  geom_vline(xintercept = 0, colour = AXIS, linewidth = 0.3) +
  geom_vline(xintercept = bg_rho, colour = MUTED, linetype = "22", linewidth = 0.3) +
  geom_segment(aes(x = 0, xend = rho, yend = gene), colour = C_N, linewidth = 0.8) +
  geom_point(size = 2, colour = C_N) +
  annotate("text", x = bg_rho, y = nrow(me) + 0.45, label = "median, all variable promoters", size = 2, colour = MUTED, hjust = 1.03) +
  labs(x = "Spearman rho, promoter methylation vs expression", y = NULL,
       title = sprintf("Methylation silences KZFPs (Mumbai, n = %d matched)", mum$n_common)) + theme_fig()
fig4 <- f4a / (f4b + f4c + plot_layout(widths = c(1, 1.25))) / (f4d + f4e) +
  plot_layout(heights = c(0.35, 1.25, 0.9)) +
  plot_annotation(tag_levels = "A",
                  title = "Figure 4. An expression proxy of the tumour-like methylome predicts transformation in an independent cohort",
                  theme = theme(plot.title = element_text(family = "Helvetica", face = "bold", size = 10)))
save_fig(fig4, "Figure4_transformation", 7.4, 8.6)

## ================================================================== FIGURE 5
pz <- ext$pz_p |> select(Sample, pmd = pmd_contrast)
cn_pick <- cnv$arms |> filter(arm_label %in% c("3p", "13q", "12p", "12q")) |>
  select(Sample, arm_label, LPO) |> pivot_wider(names_from = arm_label, values_from = LPO) |>
  mutate(chr12 = (`12p` + `12q`) / 2)
mx <- tibble(Sample = ORD) |>
  left_join(ext$sc_p |> select(Sample, Epi, tumor_likeness), by = "Sample") |>
  left_join(pz, by = "Sample") |>
  left_join(tibble(Sample = short(rownames(kz_p)), KZFP = rowMeans(kz_p)), by = "Sample") |>
  left_join(cnv$del9 |> select(Sample, del9 = lpo_9p21_median), by = "Sample") |>
  left_join(cn_pick |> select(Sample, `3p`, `13q`, chr12), by = "Sample") |>
  mutate(Patient = sub("-.*", "", Sample), Tissue = tissue_of(Sample))
feat_def <- tribble(~col, ~label, ~block, ~risk_sign, ~fmt,
  "Epi", "Epithelial\nfraction", "Composition", NA, "%.2f",
  "tumor_likeness", "Tumour-\nlikeness", "Methylome", 1, "%.2f",
  "KZFP", "KZFP\npromoter beta", "Methylome", 1, "%.2f",
  "pmd", "PMD\n(in - out)", "Methylome", -1, "%.3f",
  "3p", "3p", "Copy number", -1, "%.2f",
  "del9", "9p21\nfocal", "Copy number", -1, "%.2f",
  "13q", "13q", "Copy number", -1, "%.2f",
  "chr12", "chr12", "Copy number", 1, "%.2f")
## colour = position between the adjacent-tissue median (0) and the most extreme
## risk-direction value (1) per column; text shows the raw value
long <- map_dfr(seq_len(nrow(feat_def)), function(i) {
  f <- feat_def[i, ]; v <- mx[[f$col]]
  ref <- median(v[mx$Tissue == "Adjacent"], na.rm = TRUE)
  ## composition is context, not risk: drawn uncoloured
  s <- if (is.na(f$risk_sign)) rep(NA_real_, length(v)) else f$risk_sign * (v - ref)
  s <- if (is.na(f$risk_sign)) s else s / max(abs(s), na.rm = TRUE)
  tibble(Sample = mx$Sample, label = f$label, block = f$block, value = v, fill = s, txt = sprintf(f$fmt, v))
}) |> mutate(Sample = factor(Sample, rev(ORD)), label = factor(label, feat_def$label),
             block = factor(block, c("Composition", "Methylome", "Copy number")))
f5a <- ggplot(long, aes(label, Sample)) +
  geom_tile(aes(fill = fill), colour = "white", linewidth = 0.6) +
  geom_text(aes(label = txt), size = 2.1, colour = INK) +
  scale_fill_gradient2(low = BLUE, mid = MID, high = RED, midpoint = 0, limits = c(-1, 1), na.value = "#f9f9f7",
                       name = "Toward tumour\n(vs adjacent)", breaks = c(-1, 0, 1), labels = c("away", "adjacent", "max")) +
  facet_grid(~block, scales = "free_x", space = "free_x") +
  scale_y_discrete(labels = function(s) paste0(s, ifelse(grepl("Adj|NAT", s), "", paste0("  (", HIST[sub("-.*", "", s)], ")")))) +
  labs(x = NULL, y = NULL, title = "Pilot: integrated molecular profile vs histology") +
  theme_fig() + theme(panel.grid = element_blank())
## B: model schematic
stage <- tibble(x = 1:4, y = 3,
  lab = c("Mucosa", "OPL\nnormal-like\nmethylome (10/22)", "OPL\ntumour-like\nmethylome (12/22)", "OSCC"),
  fill = c("#eaf2fc", "#fbf3ee", "#f6d9c9", "#d6f2e7"))
notes <- tibble(x = c(2.5, 2.5, 2.5), y = c(2.15, 1.45, 0.75),
  lab = c("Tumour-like methylome: PMD erosion + chr19q13 KZFP / Polycomb-target promoter hypermethylation\n(ZNF471, ZNF382, ZNF529, ZNF582, ZNF671, HOXA, NELL1); silences KZFP expression",
          "Travels with early low-level losses at 3p, 9p, 13q, 4q; focal MTAP-CDKN2A-CDKN2B deletion (P4, P5)\nPilot progressor P5: tumour-like methylome + 3p / 9p21 / 13q loss",
          "Microenvironment (LRRC15+ CAF, Treg, mregDC) rises with grade\nbut did not predict transformation (GSE26549)"))
meas <- "Measure in one FFPE assay: EPIC methylation array gives state score + copy number;  MTAP IHC as 9p21 proxy;\nH&E deep learning trained with these as auxiliary labels; endpoint = transformation, not grade"
f5b <- ggplot() +
  geom_segment(data = tibble(x = 1.33, xend = 3.67), aes(x = x, xend = xend, y = 3, yend = 3), colour = AXIS, linewidth = 1.2,
               arrow = arrow(length = unit(2, "mm"), type = "closed")) +
  geom_label(data = stage, aes(x, y, label = lab), fill = stage$fill, size = 2.4, label.size = 0.25,
             label.padding = unit(1.5, "mm"), colour = INK, family = "Helvetica", lineheight = 0.95) +
  geom_text(data = notes, aes(x, y, label = lab), size = 2.15, colour = INK, family = "Helvetica", lineheight = 1) +
  geom_segment(data = tibble(x = c(3, 3, 3), y = c(2.6, 1.85, 1.15), yend = c(2.42, 1.67, 0.97)),
               aes(x = x, xend = x, y = y, yend = yend), colour = MUTED, linewidth = 0.3) +
  annotate("label", x = 2.5, y = 0, label = meas, size = 2.15, fill = "#f9f9f7",
           colour = INK, family = "Helvetica", lineheight = 1) +
  xlim(0.4, 4.6) + ylim(-0.4, 3.45) + theme_void() +
  labs(title = "Proposed model and measurement strategy") +
  theme(plot.title = element_text(face = "bold", size = 8, family = "Helvetica", colour = INK))
fig5 <- f5a / f5b + plot_layout(heights = c(1, 1.05)) +
  plot_annotation(tag_levels = "A", title = "Figure 5. Integrated profiles and a two-axis model of transformation risk",
                  theme = theme(plot.title = element_text(family = "Helvetica", face = "bold", size = 10)))
save_fig(fig5, "Figure5_model", 7.4, 7.2)

## ================================================================== SUPPLEMENTARY S1
pcx <- st$pc_tab |> mutate(arm = factor(arm, arm_order)) |> filter(!is.na(arm))
s1a <- ggplot(pcx, aes(arm, cancer_minus_control, fill = cancer_minus_control > 0)) +
  geom_col(width = 0.75) + geom_hline(yintercept = 0, colour = AXIS) +
  geom_text(data = filter(pcx, fdr < 0.05), aes(label = "*", y = cancer_minus_control + sign(cancer_minus_control) * 0.05),
            size = 3, colour = INK) +
  scale_fill_manual(values = c(`TRUE` = RED, `FALSE` = BLUE), guide = "none") +
  labs(x = NULL, y = "Arm dosage, cancer - control (z)",
       title = "GSE227919: expression dosage recovers the OSCC copy-number pattern (* FDR < 0.05)") +
  theme_fig() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 5.5))
s1b_d <- bind_rows(
  rna$res |> transmute(feature = set, AUC = prog_AUC, lo = prog_AUC_lo, hi = prog_AUC_hi, p = prog_wilcox_p),
  st$prog_tab |> filter(feature %in% c("3p", "9p", "core_loss")) |>
    transmute(feature = paste0("dosage loss ", feature), AUC = AUC_loss_predicts_progression, lo = NA, hi = NA, p),
  tibble(feature = "Methylome-state proxy", AUC = g26$rep2$prog_AUC, lo = NA, hi = NA, p = g26$rep2$prog_p)) |>
  mutate(feature = recode(feature, Dysplastic_epithelial = "Dysplastic epithelial", Invasive_epithelial = "Invasive epithelial",
                          Cycling_epithelial = "Cycling epithelial", LRRC15_CAF = "LRRC15+ CAF", FOLR2_TAM = "FOLR2+ TAM",
                          PI16_fibroblast = "PI16+ fibroblast", Meth_silenced_promoters = "Methylation-silenced promoters",
                          Meth_activated_promoters = "Methylation-activated promoters", KZFP_chr19_hyper = "KZFP set (expression)",
                          p16_9p21 = "CDKN2A/B, MTAP expression", `dosage loss core_loss` = "dosage loss 3p + 9p"),
         feature = fct_reorder(feature, AUC))
s1b <- ggplot(s1b_d, aes(AUC, feature)) +
  geom_vline(xintercept = 0.5, colour = AXIS) + geom_vline(xintercept = rna$hist_auc, colour = MUTED, linetype = "22") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, colour = INK2, na.rm = TRUE) +
  geom_point(aes(fill = p < 0.05), shape = 21, size = 2, colour = INK) +
  scale_fill_manual(values = c(`TRUE` = INK, `FALSE` = "white"), guide = "none") +
  annotate("text", x = rna$hist_auc, y = 0.6, label = "histology", size = 2, colour = MUTED, hjust = -0.1) +
  labs(x = "AUC for progression (25 lesions, 10 events)", y = NULL,
       title = "GSE227919 progression: small cohort, no feature survives correction") + theme_fig()
figs1 <- s1a / s1b + plot_layout(heights = c(0.8, 1.2)) +
  plot_annotation(tag_levels = "A", title = "Supplementary Figure S1. GSE227919 checks",
                  theme = theme(plot.title = element_text(family = "Helvetica", face = "bold", size = 10)))
save_fig(figs1, "FigureS1_GSE227919", 7.2, 6.4)
cat("done\n")
