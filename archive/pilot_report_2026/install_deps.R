## Install the Bioconductor packages this report needs.
## Run once:  Rscript analysis/install_deps.R

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

bioc_pkgs <- c(
  "sesame",
  "sesameData",
  "DMRcate",
  "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
  "missMethyl",
  "org.Hs.eg.db"
)
cran_pkgs <- c("matrixStats", "patchwork", "RColorBrewer", "msigdbr")

to_install_bioc <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace,
                                     logical(1), quietly = TRUE)]
to_install_cran <- cran_pkgs[!vapply(cran_pkgs, requireNamespace,
                                     logical(1), quietly = TRUE)]

if (length(to_install_cran)) {
  install.packages(to_install_cran, repos = "https://cloud.r-project.org")
}
if (length(to_install_bioc)) {
  BiocManager::install(to_install_bioc, ask = FALSE, update = FALSE)
}

## Pre-cache sesame's EPIC v2 manifest/mask resources so the report doesn't
## download them mid-render.
suppressPackageStartupMessages(library(sesameData))
sesameDataCacheAll()

cat("Done.\n")
