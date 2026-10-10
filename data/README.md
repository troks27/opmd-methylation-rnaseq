# Data

Nothing in this folder is tracked by git except this file. Data locations are set in
`config/paths.yml`; machine-specific locations go in `config/paths.local.yml`.

## Public data (`make download`)

`scripts/00_download/download_geo.sh` downloads into `data/geo/` and checks MD5 sums:

| Accession | Files |
|---|---|
| GSE204943 | GEO sample sheet and `GSE204943_RAW.tar` (IDATs, extracted to `data/geo/GSE204943/idat/`) |
| GSE85514, GSE85195, GSE26549 | series matrices; platform tables GPL22311 and GPL6480 |
| GSE234379 | series matrix and normalised beta values |
| GSE227919 | expression matrix and sample annotation (`data/geo/GSE227919/`) |

## Pilot data (not public yet)

These will be deposited in GEO on publication (accession to be added). Until then they are
available from the corresponding author.

| Data | Config entry | Expected content |
|---|---|---|
| EPIC v2 IDATs | `pilot_idat_dir` | one sub-folder per Sentrix ID, `<SentrixID>_<Position>_{Grn,Red}.idat` |
| Sample sheet | `pilot_sample_sheet` | GenomeStudio CSV (parsed by `R/lib/sample_sheet.R`) |
| RNA-seq FASTQs | `FASTQ_DIR` (environment variable) | `<library>_L00{1,2}_R{1,2}_001.fastq.gz`, two lanes per library |
| FASTQ transfer manifest | `data/pilot/rnaseq/transfer_manifest.json` | optional; file names, sizes, MD5 sums and download IDs used by `01_kallisto_sequential.sh` |

Specimen identifiers in file names and manifests must stay out of the repository. Only the
`data/README.md` file is tracked.

## Reference

`kallisto_index`: GENCODE v50 transcriptome index built with kallisto 0.52.0. The build
command is in the header of `scripts/05_pilot_rnaseq/01_kallisto_sequential.sh`.
