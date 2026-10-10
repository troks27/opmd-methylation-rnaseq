## Small helpers shared by the report.

beta_to_m <- function(b, eps = 1e-6) {
  b <- pmin(pmax(b, eps), 1 - eps)
  log2(b / (1 - b))
}

## Per-CpG line plot across Normal -> Lesion -> Tumor, with patient-level pairing.
## `cpg`: row name to plot. `betas`: CpG x sample matrix. `ss`: sample sheet.
plot_cpg_trajectory <- function(cpg, betas, ss, gene = NULL) {
  stopifnot(cpg %in% rownames(betas))
  d <- data.frame(
    Sample  = colnames(betas),
    Beta    = as.numeric(betas[cpg, ]),
    Patient = ss$Patient[match(colnames(betas), ss$Sample_Name)],
    Group   = ss$Group  [match(colnames(betas), ss$Sample_Name)]
  )
  ttl <- if (is.null(gene) || !nzchar(gene)) cpg else paste0(cpg, "  (", gene, ")")
  ggplot2::ggplot(d, ggplot2::aes(Group, Beta, group = Patient)) +
    ggplot2::geom_line(colour = "grey60") +
    ggplot2::geom_point(ggplot2::aes(colour = Group), size = 3) +
    ggplot2::geom_text(ggplot2::aes(label = Patient),
                       hjust = -0.15, size = 3, colour = "grey30") +
    ggplot2::scale_colour_manual(values = c(Normal = "#1f78b4",
                                            Lesion = "#ff7f00",
                                            Tumor  = "#e31a1c"),
                                 guide = "none") +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(title = ttl, x = NULL, y = expression(beta)) +
    ggplot2::theme_bw(base_size = 11)
}
