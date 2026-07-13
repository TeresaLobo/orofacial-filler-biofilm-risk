#!/usr/bin/env bash
set -euo pipefail

mkdir -p data/raw_genomes
mkdir -p data/unzipped_genomes
mkdir -p logs

SPECIES_FILE="species_list.txt"

if [ ! -f "$SPECIES_FILE" ]; then
  echo "[ERROR] species_list.txt not found."
  exit 1
fi

while IFS= read -r species || [ -n "$species" ]; do

  # skip empty lines
  [ -z "$species" ] && continue

  safe_name=$(echo "$species" | tr '[:upper:]' '[:lower:]' | sed 's/ /_/g')
  zipfile="data/raw_genomes/${safe_name}.zip"
  outdir="data/unzipped_genomes/${safe_name}"
  logfile="logs/${safe_name}_datasets.log"

  echo "[INFO] Downloading up to 5 RefSeq complete genomes for: $species"

  rm -f "$zipfile"
  rm -rf "$outdir"

  if datasets download genome taxon "$species" \
      --assembly-source RefSeq \
      --assembly-level complete \
      --include genome,gff3,protein,cds,seq-report \
      --limit 5 \
      --filename "$zipfile" \
      > "$logfile" 2>&1; then

      if unzip -t "$zipfile" > /dev/null 2>&1; then
          echo "[INFO] Valid ZIP downloaded for: $species"
          mkdir -p "$outdir"
          unzip -q "$zipfile" -d "$outdir"
      else
          echo "[WARN] Invalid ZIP for $species. Removing file. Check $logfile"
          rm -f "$zipfile"
          continue
      fi

  else
      echo "[WARN] Download failed for $species. Check $logfile"
      rm -f "$zipfile"
      continue
  fi

done < "$SPECIES_FILE"

echo "[INFO] Download step finished."
echo "[INFO] ZIP files:"
ls -lh data/raw_genomes/*.zip 2>/dev/null || true
