## Stage 01 / step 2. GSE204943 cell composition with EpiDISH.
##
## Robust partial correlation (RPC) against the centEpiFibIC reference (epithelial,
## fibroblast and immune cells; Zheng et al. 2018). EpiDISH uses the reference CpGs (716) that
## are present in the beta matrix.
## Input : <outputs>/gse204943/crucial_objects/{beta.rds, pd.rds}
## Output: <outputs>/gse204943/results/cell_type_proportions.csv
##         columns Epithelial, Fibroblast, Immune, Sample (= Basename), Group (Sample_Group), Group_Label
## Run from the repository root:  Rscript scripts/01_gse204943/02_cell_composition.R

source("R/lib/config.R")
suppressPackageStartupMessages({ library(EpiDISH); library(tidyverse) })
G <- file.path(OUT_DIR, "gse204943")
beta <- readRDS(file.path(G, "crucial_objects", "beta.rds"))
pd <- readRDS(file.path(G, "crucial_objects", "pd.rds"))
colnames(beta) <- pd$Basename

data(centEpiFibIC.m)
est <- epidish(beta.m = beta, ref.m = centEpiFibIC.m, method = "RPC")$estF
props <- tibble(Epithelial = est[, "Epi"], Fibroblast = est[, "Fib"], Immune = est[, "IC"], Sample = rownames(est)) |>
  mutate(Group = pd$Sample_Group[match(Sample, pd$Basename)],
         Group_Label = factor(Group, levels = c("3_Normal", "2_OPL", "1_Tumor"), labels = c("Normal", "OPL", "Tumor")))
write_csv(props, file.path(G, "results", "cell_type_proportions.csv"))
print(props |> group_by(Group_Label) |> summarise(across(Epithelial:Immune, ~ round(median(.x), 2))))
