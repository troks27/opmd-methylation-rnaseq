## Central configuration, sourced at the top of every analysis script.
##
## Resolves the repository root, reads config/paths.yml (merged with the optional,
## git-ignored config/paths.local.yml) and defines the path objects and helpers
## that the scripts use. No analysis code lives here.
## Any entry can also be overridden for one run with an environment variable named
## OPMD_<ENTRY>, e.g. OPMD_OUTPUTS=validation/run1 make all.
##
##   GEO_DIR / XDATA   public GEO downloads
##   GSE_DIR           GSE204943 processed objects (crucial_objects/, results/)
##   GSE204943_IDAT    GSE204943 IDAT files (stage 01 only)
##   RNA_DIR           GSE227919 supplementary files
##   PILOT_IDAT_DIR, PILOT_SHEET   pilot EPIC v2 IDATs and sample sheet
##   RNASEQ_DIR, KAL, MAT          pilot RNA-seq: kallisto outputs and gene matrices
##   OUT_DIR           all outputs; EXT_DIR = OUT_DIR/extended (objects, tables, figures)
##   NATURE_DIR        publication figures
##   TAB(f), FIG(f)    file paths inside EXT_DIR/tables and EXT_DIR/figures

local({
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, "config", "paths.yml"))) {
    parent <- dirname(d)
    if (parent == d) stop("Run scripts from inside the repository (config/paths.yml not found).")
    d <- parent
  }
  assign("REPO_ROOT", d, envir = globalenv())
})

local({
  cfg <- yaml::read_yaml(file.path(REPO_ROOT, "config", "paths.yml"))
  local_cfg <- file.path(REPO_ROOT, "config", "paths.local.yml")
  if (file.exists(local_cfg)) cfg <- utils::modifyList(cfg, yaml::read_yaml(local_cfg))
  for (k in names(cfg)) {
    v <- Sys.getenv(paste0("OPMD_", toupper(k)))
    if (nzchar(v)) cfg[[k]] <- v
  }
  resolve <- function(x) {
    x <- path.expand(x)
    if (!grepl("^/", x)) x <- file.path(REPO_ROOT, x)
    normalizePath(x, mustWork = FALSE)
  }
  assign("PATHS", lapply(cfg, resolve), envir = globalenv())
})

GEO_DIR        <- PATHS$geo
XDATA          <- GEO_DIR                       # name used by the analysis scripts
GSE_DIR        <- PATHS$gse204943_objects
GSE204943_IDAT <- PATHS$gse204943_idat
RNA_DIR        <- PATHS$gse227919
PILOT_IDAT_DIR <- PATHS$pilot_idat_dir
PILOT_SHEET    <- PATHS$pilot_sample_sheet
RNASEQ_DIR     <- PATHS$rnaseq
KAL            <- file.path(RNASEQ_DIR, "kallisto")
MAT            <- file.path(RNASEQ_DIR, "matrices")
KALLISTO_INDEX <- PATHS$kallisto_index
OUT_DIR        <- PATHS$outputs
EXT_DIR        <- file.path(OUT_DIR, "extended")
EXT            <- EXT_DIR                       # alias used by some scripts
NATURE_DIR     <- file.path(OUT_DIR, "figures_nature")

TAB <- function(f) file.path(EXT_DIR, "tables", f)
FIG <- function(f) file.path(EXT_DIR, "figures", f)

for (d in c(OUT_DIR, EXT_DIR, file.path(EXT_DIR, "tables"), file.path(EXT_DIR, "figures"), NATURE_DIR))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

## Fail early, with a readable message, when a required input is missing.
require_files <- function(...) {
  f <- c(...)
  miss <- f[!file.exists(f)]
  if (length(miss)) stop("Missing input(s):\n  ", paste(miss, collapse = "\n  "),
                         "\nRun the upstream stage (see README) or check config/paths*.yml.", call. = FALSE)
  invisible(TRUE)
}
