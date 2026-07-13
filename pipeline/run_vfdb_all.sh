#!/usr/bin/env bash
set -euo pipefail

mkdir -p results/virulence
mkdir -p logs

INDEX="results/tables/downloaded_genomes_index.tsv"

if [ ! -f "$INDEX" ]; then
  echo "[ERROR] Genome index not found: $INDEX"
  exit 1
fi

tail -n +2 "$INDEX" | while IFS=$'\t' read -r species accession faa gff fna; do

  if [ -z "$fna" ] || [ ! -f "$fna" ]; then
    echo "[WARN] Missing genomic FNA for $species $accession"
    continue
  fi

  out="results/virulence/${species}_${accession}_vfdb.tsv"
  log="logs/${species}_${accession}_vfdb.log"

  echo "[INFO] VFDB/ABRicate: $species $accession"

  abricate \
    --db vfdb \
    --minid 80 \
    --mincov 70 \
    "$fna" \
    > "$out" 2> "$log" || {
      echo "[WARN] ABRicate VFDB failed for $species $accession. Check $log"
      continue
    }

done

echo "[INFO] VFDB screening finished."
