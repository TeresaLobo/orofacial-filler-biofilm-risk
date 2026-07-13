#!/usr/bin/env bash
set -euo pipefail

mkdir -p results/amr
mkdir -p logs

INDEX="results/tables/downloaded_genomes_index.tsv"

if [ ! -f "$INDEX" ]; then
  echo "[ERROR] Genome index not found: $INDEX"
  echo "Run: python scripts/00_index_downloaded_genomes.py"
  exit 1
fi

tail -n +2 "$INDEX" | while IFS=$'\t' read -r species accession faa gff fna; do

  if [ -z "$fna" ] || [ ! -f "$fna" ]; then
    echo "[WARN] Missing genomic FNA for $species $accession"
    continue
  fi

  if [ -z "$faa" ] || [ ! -f "$faa" ]; then
    echo "[WARN] Missing protein FAA for $species $accession"
    continue
  fi

  if [ -z "$gff" ] || [ ! -f "$gff" ]; then
    echo "[WARN] Missing GFF for $species $accession"
    continue
  fi

  out="results/amr/${species}_${accession}_amrfinder.tsv"
  log="logs/${species}_${accession}_amrfinder.log"

  echo "[INFO] AMRFinderPlus: $species $accession"

  amrfinder \
    -n "$fna" \
    -p "$faa" \
    -g "$gff" \
    -o "$out" \
    > "$log" 2>&1 || {
      echo "[WARN] AMRFinder failed for $species $accession. Check $log"
      continue
    }

done

echo "[INFO] AMRFinderPlus finished."
