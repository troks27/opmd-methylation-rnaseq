# Software environment

| Component | Version | Used by |
|---|---|---|
| R | 4.5.2 | all analysis scripts |
| Bioconductor | 3.22 | sesame 1.28.1, minfi 1.56.0, EpiDISH 2.26.0, limma 3.66.0, edgeR 4.8.2, tximport 1.38.2, DMRcate 3.6.0 |
| conumee2 | 2.1.2 (GitHub hovestadtlab/conumee2, commit recorded in `renv.lock`) | stage 01 copy number |
| msigdbr | 26.1.0 | MSigDB hallmark sets (stage 05) |
| kallisto | 0.52.0 | stage 05 quantification |
| pandoc / python-docx | 3.9 / 1.2.0 | `tools/build_manuscript_docx.py` |

## R packages

`renv.lock` records every R package (and its dependencies) used by `scripts/` and `R/lib/`,
with its source (CRAN, Bioconductor or GitHub). To recreate the library in a fresh R 4.5.2:

```r
install.packages("renv")
renv::restore(lockfile = "env/renv.lock")
```

After editing scripts to use a new package, update the lockfile from the repository root:

```r
renv::snapshot(project = ".", lockfile = "env/renv.lock", type = "implicit")
```

sesame downloads its EPIC v2 manifests and masks on first use; cache them once with
`sesameData::sesameDataCacheAll()`.

## Command-line tools

```bash
conda env create -f env/environment.yml
conda activate opmd-tools
```

Figures are written with `quartz(type = "pdf")` (embedded Helvetica), which needs macOS.
On other systems the figure scripts fall back to `cairo_pdf`; spacing can differ slightly.
