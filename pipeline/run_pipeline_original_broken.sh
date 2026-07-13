#!/usr/bin/env bash
set -euo pipefail

# Orofacial filler-associated biofilm risk: comparative genomics pipeline
# Usage:
#   bash run_pipeline.sh
# Requirements:
#   conda/mamba env created from environment.yml
#   internet access for NCBI downloads

THREADS=${THREADS:-8}
MAX_GENOMES=${MAX_GENOMES:-30}
PROJECT_DIR=${PROJECT_DIR:-$PWD}
SPECIES_FILE=${SPECIES_FILE:-species_list.txt}

mkdir -p data/raw_genomes data/genomes_fna data/proteins data/metadata \
  data/databases results/qc results/amr results/virulence results/biofilm \
  results/matrices results/figures logs

sanitize_name() {
  echo "$1" | tr '[:upper:]' '[:lower:]' | sed 's/ /_/g; s/[^a-z0-9_]/_/g'
}

while IFS= read -r SPECIES || [[ -n "$SPECIES" ]]; do
  [[ -z "$SPECIES" ]] && continue
  SAFE=$(sanitize_name "$SPECIES")
  echo "[INFO] Downloading genomes for: $SPECIES"
 datasets download genome taxon "$SPECIES" \
    --assembly-source RefSeq \
    --include genome,protein,gff3,seq-report \
    --filename "data/raw_genomes/${SAFE}.zip" \
    > "logs/${SAFE}_datasets.log" 2>&1 || {
      echo "[WARN] Download failed for $SPECIES. Check logs/${SAFE}_datasets.log" >&2
      continue
    }
  rm -rf "data/raw_genomes/${SAFE}"
  unzip -q "data/raw_genomes/${SAFE}.zip" -d "data/raw_genomes/${SAFE}"
  if [[ -f "data/raw_genomes/${SAFE}/ncbi_dataset/data/assembly_data_report.jsonl" ]]; then
    dataformat tsv genome \
      --package "data/raw_genomes/${SAFE}/ncbi_dataset/data/assembly_data_report.jsonl" \
      > "data/metadata/${SAFE}_metadata.tsv" || true
 fi
  find "data/raw_genomes/${SAFE}/ncbi_dataset/data" -name "*.fna" | head -n "$MAX_GENOMES" | while read -r FNA; do
    ID=$(basename "$(dirname "$FNA")")
    cp "$FNA" "data/genomes_fna/${SAFE}_${ID}.fna"
  done
  find "data/raw_genomes/${SAFE}/ncbi_dataset/data" -name "*.faa" | head -n "$MAX_GENOMES" | while read -r FAA; do
    ID=$(basename "$(dirname "$FAA")")
    cp "$FAA" "data/proteins/${SAFE}_${ID}.faa"
  done
done < "$SPECIES_FILE"

# QC
if compgen -G "data/genomes_fna/*.fna" > /dev/null; then
  quast.py data/genomes_fna/*.fna -o results/qc/quast_all -t "$THREADS" || true
fi

# AMRFinderPlus
for FNA in data/genomes_fna/*.fna; do
  [[ -e "$FNA" ]] || continue
  BASE=$(basename "$FNA" .fna)
  echo "[INFO] AMRFinderPlus: $BASE"
  amrfinder -n "$FNA" -o "results/amr/${BASE}_amrfinder.tsv" --threads "$THREADS" || true
done

# VFDB with ABRicate; requires `abricate --setupdb` and VFDB database installed.
for FNA in data/genomes_fna/*.fna; do
  [[ -e "$FNA" ]] || continue
  BASE=$(basename "$FNA" .fna)
  echo "[INFO] VFDB/ABRicate: $BASE"
  abricate --db vfdb "$FNA" > "results/virulence/${BASE}_vfdb.tsv" || true
done

# Biofilm gene panel placeholder: user should place curated FASTA at data/databases/biofilm_genes.faa
if [[ -f data/databases/biofilm_genes.faa ]]; then
  diamond makedb --in data/databases/biofilm_genes.faa -d data/databases/biofilm_genes
  for FAA in data/proteins/*.faa; do
    [[ -e "$FAA" ]] || continue
    BASE=$(basename "$FAA" .faa)
    diamond blastp -q "$FAA" -d data/databases/biofilm_genes.dmnd \
      -o "results/biofilm/${BASE}_biofilm.tsv" \
      --threads "$THREADS" --id 30 --query-cover 60 --subject-cover 60 \
      --evalue 1e-5 --max-target-seqs 1 || true
  done
else
  echo "[WARN] data/databases/biofilm_genes.faa not found; skipping DIAMOND biofilm screen."
fi

python scripts/01_merge_results.py
python scripts/02_make_figures.py

echo "[DONE] Results are in results/."
