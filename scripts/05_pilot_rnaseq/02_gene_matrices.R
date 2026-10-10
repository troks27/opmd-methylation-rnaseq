## Gene-level matrices from the sequential kallisto run (scripts/05_pilot_rnaseq/01_kallisto_sequential.sh).
## Index: GENCODE v50 transcriptome (kallisto 0.52.0). Transcript -> gene map is taken
## from the GENCODE FASTA headers carried in the kallisto target names.
## Outputs (<rnaseq>/matrices/):
##   gene_counts.tsv.gz              tximport estimated counts (countsFromAbundance = "no"; use with
##                                   the length offset in gene_tx_length.tsv.gz for DESeq2/edgeR)
##   gene_counts_lengthScaledTPM.tsv.gz  counts scaled to library size and average transcript length
##   gene_tpm.tsv.gz                 gene TPM (sum of transcript TPM)
##   gene_tx_length.tsv.gz           tximport average transcript length per gene and sample
##   gene_annotation.tsv             gene_id, gene_name, gene_type
##   samples.tsv                     library, patient, tissue, methylation-array match, QC
source("R/lib/config.R")   # paths and helpers (TAB, FIG, OUT_DIR, GSE_DIR, ...)
suppressPackageStartupMessages({ library(tximport); library(readr); library(dplyr) })

dir.create(MAT, showWarnings = FALSE)
libs <- list.dirs(KAL, recursive = FALSE, full.names = FALSE)
libs <- libs[file.exists(file.path(KAL, libs, "abundance.tsv.gz"))]
libs <- libs[order(as.integer(sub(".*_S", "", libs)))]
stopifnot(length(libs) > 0)

files <- setNames(file.path(KAL, libs, "abundance.tsv.gz"), libs)
tx <- read_tsv(files[1], col_types = "c____", col_names = "target_id", skip = 1)$target_id
parts <- strsplit(tx, "|", fixed = TRUE)
tx2gene <- data.frame(tx = tx, gene = vapply(parts, `[`, "", 2))
ann <- tibble(gene_id = vapply(parts, `[`, "", 2), gene_name = vapply(parts, `[`, "", 6),
              gene_type = vapply(parts, function(p) p[length(p)], "")) |> distinct(gene_id, .keep_all = TRUE)

txi  <- tximport(files, type = "kallisto", tx2gene = tx2gene, countsFromAbundance = "no")
txiL <- tximport(files, type = "kallisto", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")

wr <- function(m, f) write_tsv(as_tibble(round(m, 3), rownames = "gene_id"), file.path(MAT, f))
wr(txi$counts, "gene_counts.tsv.gz"); wr(txiL$counts, "gene_counts_lengthScaledTPM.tsv.gz")
wr(txi$abundance, "gene_tpm.tsv.gz"); wr(txi$length, "gene_tx_length.tsv.gz")
write_tsv(ann, file.path(MAT, "gene_annotation.tsv"))

## sample sheet: library -> patient / tissue / matching methylation sample (same specimen IDs)
qc <- lapply(libs, function(l) jsonlite::fromJSON(file.path(KAL, l, "run_info.json")))
ss <- tibble(library = libs,
             patient = paste0("P", sub("^Patient-(\\d+)-.*", "\\1", libs)),
             tissue = ifelse(grepl("Tumor", libs), "Tumor", ifelse(grepl("Lesion", libs), "Lesion", "Adjacent")),
             specimen_id = sub(".*-(\\d+T?\\d+)_S\\d+$", "\\1", libs),
             n_pairs = vapply(qc, `[[`, 0, "n_processed"),
             n_pseudoaligned = vapply(qc, `[[`, 0, "n_pseudoaligned"),
             pct_pseudoaligned = vapply(qc, `[[`, 0, "p_pseudoaligned")) |>
  mutate(methylation_sample = c(P2_Tumor = "P2-Tumor", P2_Adjacent = "P2-NAT", P3_Lesion = "P3-Lesion",
                                P3_Adjacent = "P3-Adj", P4_Lesion = "P4-Lesion", P4_Adjacent = "P4-Adj",
                                P5_Lesion = "P5-Lesion", P5_Adjacent = "P5-Adj")[paste(patient, tissue, sep = "_")])
write_tsv(ss, file.path(MAT, "samples.tsv"))
print(ss, width = 200)
cat(sprintf("\n%d genes x %d libraries; median %.1f%% pseudoaligned\n", nrow(txi$counts), ncol(txi$counts),
            median(ss$pct_pseudoaligned)))
