## Compare the CSV tables of two output folders, e.g. before and after an edit:
##   Rscript tools/compare_tables.R <old>/extended/tables <new>/extended/tables [tolerance]
## For each table present in both folders it reports whether it is identical, numerically
## equal within the tolerance (default 1e-6, relative), or different, with the largest
## absolute difference per differing column. Writes nothing; prints a summary.

suppressPackageStartupMessages(library(readr))
args <- commandArgs(TRUE)
if (length(args) < 2) stop("usage: Rscript tools/compare_tables.R <old_dir> <new_dir> [tolerance]")
old_dir <- args[1]; new_dir <- args[2]; tol <- if (length(args) > 2) as.numeric(args[3]) else 1e-6

f_old <- list.files(old_dir, pattern = "\\.csv$"); f_new <- list.files(new_dir, pattern = "\\.csv$")
only_old <- setdiff(f_old, f_new); only_new <- setdiff(f_new, f_old)
rd <- function(f) suppressWarnings(read_csv(f, show_col_types = FALSE, guess_max = 1e5, progress = FALSE))

compare <- function(f) {
  a <- rd(file.path(old_dir, f)); b <- rd(file.path(new_dir, f))
  if (identical(a, b)) return(data.frame(table = f, status = "identical", detail = ""))
  if (!identical(dim(a), dim(b)) || !identical(names(a), names(b)))
    return(data.frame(table = f, status = "DIFFERENT",
                      detail = sprintf("shape %s -> %s; columns %s", paste(dim(a), collapse = "x"),
                                       paste(dim(b), collapse = "x"),
                                       if (identical(names(a), names(b))) "same" else "differ")))
  msg <- character()
  for (k in names(a)) {
    x <- a[[k]]; y <- b[[k]]
    if (is.numeric(x) && is.numeric(y)) {
      d <- abs(x - y); s <- pmax(abs(x), abs(y), 1e-12)
      bad <- (is.na(x) != is.na(y)) | (!is.na(d) & d / s > tol)
      if (any(bad)) msg <- c(msg, sprintf("%s: %d values, max |diff| %.3g", k, sum(bad), max(d[bad], na.rm = TRUE)))
    } else if (!identical(as.character(x), as.character(y))) {
      msg <- c(msg, sprintf("%s: %d values", k, sum(as.character(x) != as.character(y), na.rm = TRUE)))
    }
  }
  if (!length(msg)) data.frame(table = f, status = "equal (tolerance)", detail = "")
  else data.frame(table = f, status = "DIFFERENT", detail = paste(msg, collapse = "; "))
}

res <- do.call(rbind, lapply(intersect(f_old, f_new), compare))
cat(sprintf("%d tables compared: %s\n", nrow(res), paste(names(table(res$status)), table(res$status), collapse = ", ")))
if (length(only_old)) cat("Only in", old_dir, ":", paste(only_old, collapse = ", "), "\n")
if (length(only_new)) cat("Only in", new_dir, ":", paste(only_new, collapse = ", "), "\n")
d <- res[res$status == "DIFFERENT", ]
if (nrow(d)) { cat("\nDifferent tables:\n"); for (i in seq_len(nrow(d))) cat(sprintf("  %s\n      %s\n", d$table[i], d$detail[i])) }
