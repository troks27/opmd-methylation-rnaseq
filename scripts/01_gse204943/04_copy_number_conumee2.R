## Stage 01 / step 4. GSE204943 copy number from methylation-array intensities (conumee 2.0).
##
## Total intensities (methylated + unmethylated) after preprocessIllumina with its defaults
## (Illumina background correction and normalisation to the internal controls), as in the
## original analysis. Each OPL and OSCC is fitted against the healthy mucosae as the reference,
## then binned and segmented with conumee2 defaults for EPIC arrays (hg19). Segmentation uses
## DNAcopy permutations, so the seed is fixed for every sample.
## Input : <outputs>/gse204943/intermediate/RGset_qc.rds, crucial_objects/pd.rds
## Output: <outputs>/gse204943/results/CNV_segments_all.csv  (conumee2 segment table + Sample, Group)
## conumee2 is installed from GitHub: remotes::install_github("hovestadtlab/conumee2", subdir = "conumee2")
## Run from the repository root:  Rscript scripts/01_gse204943/04_copy_number_conumee2.R  (~30 min)

source("R/lib/config.R")
suppressPackageStartupMessages({ library(minfi); library(conumee2); library(tidyverse) })
G <- file.path(OUT_DIR, "gse204943")
RGset <- readRDS(file.path(G, "intermediate", "RGset_qc.rds"))
pd <- readRDS(file.path(G, "crucial_objects", "pd.rds"))

Mset <- preprocessIllumina(RGset)
tot <- getMeth(Mset) + getUnmeth(Mset)
colnames(tot) <- pd$Basename
rm(RGset, Mset); invisible(gc())

cn <- CNV.load(tot)
data(detail_regions)
anno <- CNV.create_anno(array_type = "EPIC", detail_regions = detail_regions, features = rownames(tot))
ref <- pd$Sample_Group == "3_Normal"
cat(sprintf("Reference: %d healthy mucosae; query: %d lesions\n", sum(ref), sum(!ref)))

segs <- map_dfr(which(!ref), function(i) {
  set.seed(1)   # DNAcopy assesses change points by permutation
  x <- CNV.segment(CNV.detail(CNV.bin(CNV.fit(cn[i], cn[ref], anno))))
  CNV.write(x, what = "segments") |> bind_rows() |>          # conumee2 returns one table per query sample
    mutate(Sample = pd$Basename[i], Group = pd$Sample_Group[i])
})
write_csv(segs, file.path(G, "results", "CNV_segments_all.csv"))
cat(sprintf("Segments: %d across %d samples\n", nrow(segs), n_distinct(segs$Sample)))
