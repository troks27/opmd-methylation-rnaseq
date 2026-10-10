#!/bin/bash
## Stage 05 / step 1. Pilot bulk RNA-seq quantification with kallisto, one library at a time.
##
## Library: paired-end, two lanes per library, reverse-stranded (dUTP). Strandedness was
## checked on 300,000 read pairs of S1: --rf-stranded 77.0%, --fr-stranded 3.6%,
## unstranded 78.2% pseudoaligned.
## Index:   GENCODE v50 transcripts, kallisto 0.52.0, k = 31 (build once, see below).
## Output:  <rnaseq>/kallisto/<library>/{abundance.tsv.gz, run_info.json, kallisto.log}
##          (<rnaseq> from config/paths.yml; abundance.h5 is not kept.)
## Resumable: libraries that already have abundance.tsv.gz are skipped.
##
## FASTQ source (the FASTQs are not public; see data/README.md), either
##   FASTQ_DIR=/path/to/fastqs bash scripts/05_pilot_rnaseq/01_kallisto_sequential.sh
##       local files named <library>_L00{1,2}_R{1,2}_001.fastq.gz; nothing is deleted, or
##   bash scripts/05_pilot_rnaseq/01_kallisto_sequential.sh
##       downloads each library's four files listed in data/pilot/rnaseq/transfer_manifest.json
##       (fields name, bytes, md5, drive_id; git-ignored), checks size and MD5, quantifies,
##       then deletes the FASTQs, so that only one library (<= 2.7 GB) is on disk at a time.
##
## Build the index once (about 15 min, 1.3 GB):
##   curl -O https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.transcripts.fa.gz
##   kallisto index -i <kallisto_index> gencode.v50.transcripts.fa.gz
##
## Run from the repository root. Needs kallisto 0.52.0 on PATH (env/environment.yml).
set -euo pipefail

cd "$(dirname "$0")/../.."
eval "$(Rscript -e 'source("R/lib/config.R"); cat(sprintf("KAL=%s\nIDX=%s\n", shQuote(KAL), shQuote(KALLISTO_INDEX)))')"
KALLISTO="${KALLISTO:-kallisto}"
MAN="${MANIFEST:-data/pilot/rnaseq/transfer_manifest.json}"
TMP="${TMPDIR_FASTQ:-$HOME/rnaseq_work/tmp}"      # keep FASTQs outside cloud-synced folders
LOG="logs/kallisto_run.log"
THREADS="${THREADS:-8}"
STRAND="--rf-stranded"
MIN_FREE_GB=5
URL="https://drive.usercontent.google.com/download?export=download&confirm=t&id="
mkdir -p "$KAL" logs

[ -s "$IDX" ] || { echo "kallisto index not found: $IDX (see the header of this script)"; exit 1; }
"$KALLISTO" version | grep -q "0.52.0" || echo "warning: the published results used kallisto 0.52.0"

log() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }
free_gb() { df -k "$HOME" | awk 'NR==2 {printf "%.1f", $4/1048576}'; }
field() { python3 -c "import json,sys; d={x['name']:x for x in json.load(open(sys.argv[3]))}; print(d[sys.argv[1]][sys.argv[2]])" "$1" "$2" "$MAN"; }

fetch() {  # download one file from the manifest and verify size + MD5
  local f="$1" id bytes md5 got
  id=$(field "$f" drive_id); bytes=$(field "$f" bytes); md5=$(field "$f" md5)
  for attempt in 1 2 3; do
    curl -sSL --retry 3 --retry-delay 5 -o "$TMP/$f" "${URL}${id}" || true
    got=$(wc -c < "$TMP/$f" 2>/dev/null | tr -d ' ' || echo 0)
    if [ "$got" = "$bytes" ] && [ "$(md5sum_of "$TMP/$f")" = "$md5" ]; then return 0; fi
    log "   retry $attempt for $f (size $got of $bytes)"; rm -f "$TMP/$f"; sleep 10
  done
  log "FAILED to download a verified copy of $f"; return 1
}
md5sum_of() { if command -v md5 >/dev/null; then md5 -q "$1"; else md5sum "$1" | cut -d' ' -f1; fi; }

if [ -n "${FASTQ_DIR:-}" ]; then
  SRC="$FASTQ_DIR"
  libs=$(ls "$FASTQ_DIR" | sed -n 's/_L001_R1_001\.fastq\.gz$//p' | python3 -c "import sys; print('\n'.join(sorted(sys.stdin.read().split(), key=lambda s: int(s.rsplit('_S', 1)[-1]))))")
else
  [ -s "$MAN" ] || { echo "Set FASTQ_DIR or provide $MAN"; exit 1; }
  SRC="$TMP"; mkdir -p "$TMP"
  libs=$(python3 -c "import json,sys; m=json.load(open(sys.argv[1])); print('\n'.join(sorted({x['name'].split('_L00')[0] for x in m}, key=lambda s: int(s.rsplit('_S', 1)[-1]))))" "$MAN")
fi

for s in $libs; do
  if [ -s "$KAL/$s/abundance.tsv.gz" ]; then log "skip $s (done)"; continue; fi
  if [ -z "${FASTQ_DIR:-}" ]; then
    need=$(python3 -c "import json,sys; m=json.load(open(sys.argv[1])); print(round(sum(x['bytes'] for x in m if x['name'].startswith(sys.argv[2]+'_L00'))/1073741824,2))" "$MAN" "$s")
    fr=$(free_gb)
    python3 -c "import sys; sys.exit(0 if $fr - $need >= $MIN_FREE_GB else 1)" || {
      log "STOP: $s needs ${need} GB, only ${fr} GB free (floor ${MIN_FREE_GB} GB)"; exit 2; }
    log "== $s: downloading ${need} GB; ${fr} GB free"
    for lane in L001 L002; do for r in R1 R2; do fetch "${s}_${lane}_${r}_001.fastq.gz"; done; done
  else
    log "== $s"
  fi

  tmpo="$KAL/.$s.tmp"; rm -rf "$tmpo"
  "$KALLISTO" quant -i "$IDX" -o "$tmpo" -t "$THREADS" $STRAND \
    "$SRC/${s}_L001_R1_001.fastq.gz" "$SRC/${s}_L001_R2_001.fastq.gz" \
    "$SRC/${s}_L002_R1_001.fastq.gz" "$SRC/${s}_L002_R2_001.fastq.gz" > "$tmpo.log" 2>&1 \
    || { log "kallisto FAILED for $s"; cat "$tmpo.log" >> "$LOG"; exit 4; }
  gzip -f "$tmpo/abundance.tsv"; rm -f "$tmpo/abundance.h5"; mv "$tmpo.log" "$tmpo/kallisto.log"
  rm -rf "$KAL/$s"; mv "$tmpo" "$KAL/$s"
  [ -z "${FASTQ_DIR:-}" ] && rm -f "$TMP/${s}"_L00*_R*_001.fastq.gz
  qc=$(python3 -c "import json,sys; r=json.load(open(sys.argv[1])); print('{:,} pairs, {}% pseudoaligned'.format(r['n_processed'], r['p_pseudoaligned']))" "$KAL/$s/run_info.json")
  log "   $qc"
done
log "all libraries done"
