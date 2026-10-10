#!/bin/bash
## Stage 05 / step 0. Library strandedness, checked before quantification.
##
## The first 300,000 read pairs of the first library (S1, lane 1) are pseudoaligned with
## kallisto in unstranded, --fr-stranded and --rf-stranded modes. A reverse-stranded (dUTP)
## library pseudoaligns almost as well with --rf-stranded as unstranded, and hardly at all
## with --fr-stranded. 01_kallisto_sequential.sh therefore uses --rf-stranded.
## Output: <rnaseq>/strandedness_check.tsv  (mode, n_processed, n_pseudoaligned, p_pseudoaligned),
##         shown in Supplementary Fig. 3b.
##
## FASTQ source as in 01_kallisto_sequential.sh: FASTQ_DIR=/path/to/fastqs, or the manifest
## (only the first 40 MiB of each file is downloaded).
## Run from the repository root:  bash scripts/05_pilot_rnaseq/00_strandedness_check.sh
set -euo pipefail

cd "$(dirname "$0")/../.."
eval "$(Rscript -e 'source("R/lib/config.R"); cat(sprintf("RNA=%s\nIDX=%s\n", shQuote(RNASEQ_DIR), shQuote(KALLISTO_INDEX)))')"
KALLISTO="${KALLISTO:-kallisto}"
MAN="${MANIFEST:-data/pilot/rnaseq/transfer_manifest.json}"
N_PAIRS=300000
URL="https://drive.usercontent.google.com/download?export=download&confirm=t&id="
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT

if [ -n "${FASTQ_DIR:-}" ]; then
  r1=$(ls "$FASTQ_DIR"/*_S1_L001_R1_001.fastq.gz); r2=${r1/_R1_/_R2_}
  gzip -dc "$r1" | head -n $((4 * N_PAIRS)) > "$W/R1.fq" || true
  gzip -dc "$r2" | head -n $((4 * N_PAIRS)) > "$W/R2.fq" || true
else
  for r in R1 R2; do
    id=$(python3 -c "import json,sys; print([x['drive_id'] for x in json.load(open(sys.argv[1])) if '_S1_L001_'+sys.argv[2]+'_' in x['name']][0])" "$MAN" "$r")
    curl -sSL -r 0-41943039 -o "$W/$r.part" "${URL}${id}"
    gzip -dc "$W/$r.part" 2>/dev/null | head -n $((4 * N_PAIRS)) > "$W/$r.fq" || true
  done
fi
[ "$(wc -l < "$W/R1.fq")" -eq $((4 * N_PAIRS)) ] || { echo "fewer than $N_PAIRS reads extracted"; exit 1; }

printf "mode\tn_processed\tn_pseudoaligned\tp_pseudoaligned\n" > "$RNA/strandedness_check.tsv"
for m in unstranded --fr-stranded --rf-stranded; do
  flag=$([ "$m" = unstranded ] && echo "" || echo "$m")
  "$KALLISTO" quant -i "$IDX" -o "$W/$m" -t 8 $flag "$W/R1.fq" "$W/R2.fq" > "$W/$m.log" 2>&1
  python3 -c "import json,sys; r=json.load(open(sys.argv[1])); print('\t'.join(map(str, [sys.argv[2], r['n_processed'], r['n_pseudoaligned'], r['p_pseudoaligned']])))" "$W/$m/run_info.json" "$m" >> "$RNA/strandedness_check.tsv"
done
cat "$RNA/strandedness_check.tsv"
