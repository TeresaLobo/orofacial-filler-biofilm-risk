#!/usr/bin/env bash
set -euo pipefail

mkdir -p results/abricate_extra
mkdir -p logs

INDEX="results/tables/downloaded_genomes_index.tsv"

DBS=("card" "resfinder" "bacmet2" "victors" "plasmidfinder")

for db in "${DBS[@]}"; do
  mkdir -p "results/abricate_extra/${db}"

  tail -n +2 "$INDEX" | while IFS=$'\t' read -r species accession faa gff fna; do

    if [ -z "$fna" ] || [ ! -f "$fna" ]; then
      echo "[WARN] Missing genomic FNA for $species $accession"
      continue
    fi

    out="results/abricate_extra/${db}/${species}_${accession}_${db}.tsv"
    log="logs/${species}_${accession}_${db}.log"

    echo "[INFO] ABRicate $db: $species $accession"

    abricate \
      --db "$db" \
      --minid 80 \
      --mincov 70 \
      "$fna" \
      > "$out" 2> "$log" || {
        echo "[WARN] ABRicate $db failed for $species $accession. Check $log"
        continue
      }

  done
done

echo "[INFO] Extra ABRicate databases finished."
