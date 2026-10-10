## Build a tidy sample sheet from the GenomeStudio-style CSV shipped with the project.
## Returns a data.frame with: Sample_Name, Patient, Group (factor Normal<Lesion<Tumor),
## Sentrix_ID, Sentrix_Position, Basename (full IDAT prefix path).

build_sample_sheet <- function(csv_path, idat_dir) {
  raw <- readLines(csv_path)
  data_start <- grep("^\\[Data\\]", raw)
  if (length(data_start) != 1L) stop("Could not locate [Data] section in ", csv_path)

  ss <- read.csv(csv_path, skip = data_start, stringsAsFactors = FALSE,
                 strip.white = TRUE, blank.lines.skip = TRUE)
  ss <- ss[nzchar(ss$Sample_Name), , drop = FALSE]

  ## Derive Patient ("Patient 2", "Patient 3", ...)
  ss$Patient <- sub("^(Patient\\s*\\d+).*", "\\1", ss$Sample_Name)

  ## Derive Group: anything mentioning "tumor" or "lesion" -> those; otherwise Normal.
  name_lc <- tolower(ss$Sample_Name)
  ss$Group <- ifelse(grepl("tumor", name_lc), "Tumor",
              ifelse(grepl("lesion", name_lc), "Lesion", "Normal"))
  ss$Group <- factor(ss$Group, levels = c("Normal", "Lesion", "Tumor"), ordered = TRUE)

  ## Path prefix used by sesame / minfi to find paired Grn+Red IDATs.
  ss$Basename <- file.path(idat_dir, ss$Sentrix_ID,
                           paste0(ss$Sentrix_ID, "_", ss$Sentrix_Position))

  ## Sanity: confirm both IDAT files exist for every sample.
  missing_idats <- vapply(ss$Basename, function(b) {
    !(file.exists(paste0(b, "_Grn.idat")) && file.exists(paste0(b, "_Red.idat")))
  }, logical(1))
  if (any(missing_idats)) {
    stop("Missing IDATs for: ",
         paste(ss$Sample_Name[missing_idats], collapse = ", "))
  }

  ss[, c("Sample_Name", "Patient", "Group", "Sentrix_ID",
         "Sentrix_Position", "Sample_Plate", "Sample_Well", "Basename")]
}
