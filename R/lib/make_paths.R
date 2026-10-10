## Writes the configured output locations as Make variable assignments to the file
## given as the first argument (the Makefile includes it). Paths inside the repository
## are written relative to it, so Make never sees the absolute path of the checkout.
## The file is rewritten only when its content changes.
source("R/lib/config.R")
rel <- function(p) {
  root <- paste0(normalizePath(REPO_ROOT), "/")
  p <- normalizePath(p, mustWork = FALSE)
  if (startsWith(p, root)) substring(p, nchar(root) + 1) else p
}
out <- rel(OUT_DIR)
gse <- rel(GSE_DIR)
## GSE204943 objects are Make targets only when stage 01 builds them inside the outputs
## folder; a frozen copy elsewhere is treated as a fixed input.
lines <- c(sprintf("OUT := %s", out),
           sprintf("EXT := %s/extended", out),
           sprintf("FIGDIR := %s/figures_nature", out),
           sprintf("MAT := %s", rel(MAT)),
           sprintf("KAL := %s", rel(KAL)),
           sprintf("RNASEQ := %s", rel(RNASEQ_DIR)),
           sprintf("GSE_BUILT := %s", if (gse == paste0(out, "/gse204943")) "yes" else ""))
f <- commandArgs(TRUE)[1]
if (is.na(f)) writeLines(lines) else if (!file.exists(f) || !identical(readLines(f), lines)) writeLines(lines, f)
