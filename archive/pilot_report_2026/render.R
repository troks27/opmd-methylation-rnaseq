## Render the methylation report to HTML.
## Run from the analysis/ directory:
##   cd analysis && Rscript render.R
## or from the project root:
##   Rscript -e 'rmarkdown::render("analysis/methylation_report.Rmd", output_file = "output/methylation_report.html")'

if (!file.exists("methylation_report.Rmd") && file.exists("analysis/methylation_report.Rmd")) {
  setwd("analysis")
}
stopifnot(file.exists("methylation_report.Rmd"))

## use the project-local pandoc when no system pandoc is on the PATH
if (!rmarkdown::pandoc_available() && file.exists(".bin/pandoc"))
  Sys.setenv(RSTUDIO_PANDOC = normalizePath(".bin"))

## Extended analyses (array copy number, GSE204943 methylation, GSE227919
## RNA-seq) write to output/extended/ and are read by report sections 9-14.
## Skip with:  Rscript render.R --skip-extended
if (!"--skip-extended" %in% commandArgs(trailingOnly = TRUE)) {
  for (f in c("R/10_sigdf_qc_cnv.R", "R/11_external_methylation.R", "R/12_external_rnaseq.R")) {
    message("Running ", f)
    stopifnot(system2("Rscript", f) == 0)
  }
}

rmarkdown::render(
  input       = "methylation_report.Rmd",
  output_file = "output/methylation_report.html",
  envir       = new.env()
)
