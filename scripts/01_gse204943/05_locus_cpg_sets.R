## Stage 01 / step 5. ZNF471 CpG table.
##
## All CpGs annotated to ZNF471 (UCSC RefGene), each with the RefGene group of its first
## ZNF471 entry ("Region"), a coarse class (Promoter = TSS1500, TSS200, 5'UTR or 1stExon;
## Gene Body = Body or ExonBnd; 3'UTR; Other) and its mean beta per GSE204943 group.
## Later scripts use the promoter CpGs as the ZNF471 promoter set.
## Input : <outputs>/gse204943/crucial_objects/{beta.rds, pd.rds, ann_auto.rds}
## Output: <outputs>/gse204943/results/ZNF471_CpG_details.csv
## Run from the repository root:  Rscript scripts/01_gse204943/05_locus_cpg_sets.R

source("R/lib/config.R")
suppressPackageStartupMessages(library(tidyverse))
G <- file.path(OUT_DIR, "gse204943")
beta <- readRDS(file.path(G, "crucial_objects", "beta.rds"))
pd <- readRDS(file.path(G, "crucial_objects", "pd.rds"))
ann <- as.data.frame(readRDS(file.path(G, "crucial_objects", "ann_auto.rds")))
colnames(beta) <- pd$Basename
GENE <- "ZNF471"

region_of <- function(names, groups) {      # RefGene group of the first entry for GENE
  n <- strsplit(names, ";")[[1]]; g <- strsplit(groups, ";")[[1]]; i <- which(n == GENE)[1]
  if (!is.na(i) && i <= length(g)) g[i] else g[1]
}
class_of <- function(r) case_when(r %in% c("TSS200", "TSS1500", "5'UTR", "1stExon") ~ "Promoter",
                                  r %in% c("Body", "ExonBnd") ~ "Gene Body", r == "3'UTR" ~ "3'UTR", TRUE ~ "Other")

cpg <- ann |> filter(grepl(paste0("(^|;)", GENE, "(;|$)"), UCSC_RefGene_Name), Name %in% rownames(beta)) |>
  mutate(Gene = GENE, Region = map2_chr(UCSC_RefGene_Name, UCSC_RefGene_Group, region_of)) |>
  select(Name, Gene, Region, chr, pos, strand, Relation_to_Island, UCSC_RefGene_Name, UCSC_RefGene_Group) |>
  mutate(Region_Class = class_of(Region))
grp_mean <- function(g) rowMeans(beta[cpg$Name, pd$Sample_Group == g, drop = FALSE], na.rm = TRUE)
out <- cpg |> mutate(Normal = grp_mean("3_Normal"), OPL = grp_mean("2_OPL"), Tumor = grp_mean("1_Tumor"),
                     Delta_TvN = Tumor - Normal, Progressive = sign(OPL - Normal) == sign(Tumor - Normal)) |>
  arrange(Region, pos)
write_csv(out, file.path(G, "results", "ZNF471_CpG_details.csv"))
cat(sprintf("%s CpGs: %d (%d promoter)\n", GENE, nrow(out), sum(out$Region_Class == "Promoter")))
