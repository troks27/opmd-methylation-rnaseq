#!/bin/bash
## Download the public GEO data used by the analysis into data/geo (config key: geo).
##
##   GSE204943  EPIC v1 methylation, Tata Memorial (22 healthy, 22 leukoplakia, 74 OSCC):
##              raw IDAT archive (1.9 GB) + sample sheet
##   GSE85514   array CGH, same Mumbai patients        (series matrix + GPL22311 table)
##   GSE85195   Agilent expression, same Mumbai patients (series matrix + GPL6480 table)
##   GSE234379  450K methylation, 46 OSCC/adjacent pairs, Sao Paulo (series matrix + betas)
##   GSE26549   HuGene 1.0 ST expression, 86 leukoplakias with oral cancer outcome
##   GSE227919  total RNA-seq, 66 oral biopsies, Boston/Buffalo (matrix + annotation)
##
## Files already present with the expected MD5 are skipped. A checksum mismatch is
## reported (GEO occasionally regenerates platform tables) but does not stop the run.
## Usage: bash scripts/00_download/download_geo.sh        (from the repository root)
set -euo pipefail
cd "$(dirname "$0")/../.."
GEO=$(Rscript -e 'source("R/lib/config.R"); cat(GEO_DIR)')
G227=$(Rscript -e 'source("R/lib/config.R"); cat(RNA_DIR)')
GIDAT=$(Rscript -e 'source("R/lib/config.R"); cat(GSE204943_IDAT)')
mkdir -p "$GEO" "$G227" "$GIDAT"
FTP=https://ftp.ncbi.nlm.nih.gov/geo
ACC="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?targ=self&form=text&view=full&acc="

## file | url | md5 of the copy used for the manuscript ("-" = not recorded)
get() {
  local dest="$1" url="$2" md5="$3"
  if [ -s "$dest" ] && { [ "$md5" = "-" ] || [ "$(md5sum_ "$dest")" = "$md5" ]; }; then
    echo "ok       $(basename "$dest")"; return; fi
  echo "download $(basename "$dest")"
  curl -fsSL --retry 3 -o "$dest.part" "$url" && mv "$dest.part" "$dest"
  if [ "$md5" != "-" ] && [ "$(md5sum_ "$dest")" != "$md5" ]; then
    echo "  WARNING: MD5 differs from the copy used for the manuscript ($md5)"; fi
}
md5sum_() { if command -v md5 >/dev/null; then md5 -q "$1"; else md5sum "$1" | cut -d' ' -f1; fi; }

get "$GEO/GSE85514_series_matrix.txt.gz"  "$FTP/series/GSE85nnn/GSE85514/matrix/GSE85514_series_matrix.txt.gz"   05d6f69b01f2850db441818501a28138
get "$GEO/GSE85195_series_matrix.txt.gz"  "$FTP/series/GSE85nnn/GSE85195/matrix/GSE85195_series_matrix.txt.gz"   52947f1845defd44dc9b6cf182658e62
get "$GEO/GSE234379_series_matrix.txt.gz" "$FTP/series/GSE234nnn/GSE234379/matrix/GSE234379_series_matrix.txt.gz" 82b1f63f5759f05826877a951e561cf7
get "$GEO/GSE26549_series_matrix.txt.gz"  "$FTP/series/GSE26nnn/GSE26549/matrix/GSE26549_series_matrix.txt.gz"   c03a23d6f3bca62c21a96fbad0237425
get "$GEO/GSE234379_normalized_data.txt.gz" "$FTP/series/GSE234nnn/GSE234379/suppl/GSE234379_normalized_data.txt.gz" 928554f783a59ccf06912ab76bf8115c
get "$GEO/GPL22311.annot.txt" "${ACC}GPL22311" aa26900924bcc440b669bbd51c304941
get "$GEO/GPL6480.annot.txt"  "${ACC}GPL6480"  7aa299582ebf8b653e4af097696ba447
get "$G227/GSE227919_2023_20_03_expr_mat_raw.txt.gz"     "$FTP/series/GSE227nnn/GSE227919/suppl/GSE227919_2023_20_03_expr_mat_raw.txt.gz"     887f7933c94116fa27c3b65441edc58f
get "$G227/GSE227919_2023_21_03_pml_annotations.csv.gz"  "$FTP/series/GSE227nnn/GSE227919/suppl/GSE227919_2023_21_03_pml_annotations.csv.gz"  f93be2362cac3f14cc5c531e34550f91

## GSE204943: sample sheet + IDAT archive (stage 01); GEO sample metadata and file list
## (06_audit/04_gse204943_integrity.R)
G2=$(dirname "$GIDAT")
get "$G2/GSE204943_Sample_sheet_118_samples_MBM_MI_modified.csv.gz" \
    "$FTP/series/GSE204nnn/GSE204943/suppl/GSE204943_Sample_sheet_118_samples_MBM_MI_modified.csv.gz" -
get "$G2/GSE204943_series_matrix.txt.gz" "$FTP/series/GSE204nnn/GSE204943/matrix/GSE204943_series_matrix.txt.gz" 675df351df3db19d1698ef6a0f56b061
get "$G2/GSE204943_filelist.txt"         "$FTP/series/GSE204nnn/GSE204943/suppl/filelist.txt"                   4dc59ddfece2029e61caffd76b6ac167
if [ -z "$(ls "$GIDAT"/*_Grn.idat* 2>/dev/null | head -1)" ]; then
  get "$G2/GSE204943_RAW.tar" "$FTP/series/GSE204nnn/GSE204943/suppl/GSE204943_RAW.tar" -
  echo "extract  GSE204943_RAW.tar -> $GIDAT"
  tar -xf "$G2/GSE204943_RAW.tar" -C "$GIDAT"
  echo "  (the 1.9 GB archive can be deleted once extraction is verified: $G2/GSE204943_RAW.tar)"
else
  echo "ok       GSE204943 IDATs ($(ls "$GIDAT" | wc -l | tr -d ' ') files)"
fi
echo "Done. Public data in $GEO"
