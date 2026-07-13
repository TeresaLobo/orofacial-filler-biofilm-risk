from pathlib import Path
import pandas as pd

base = Path("data/unzipped_genomes")
rows = []

for species_dir in sorted(base.iterdir()):
    if not species_dir.is_dir():
        continue

    species = species_dir.name
    data_dir = species_dir / "ncbi_dataset" / "data"

    if not data_dir.exists():
        continue

    for acc_dir in sorted(data_dir.glob("GCF_*")):
        faa = acc_dir / "protein.faa"
        gff = acc_dir / "genomic.gff"

        # Important: exclude cds_from_genomic.fna and keep the true genomic assembly FASTA
        fna_files = [
            f for f in acc_dir.glob("*_genomic.fna")
            if f.name != "cds_from_genomic.fna"
        ]

        rows.append({
            "species_folder": species,
            "accession": acc_dir.name,
            "protein_faa": str(faa) if faa.exists() else "",
            "genomic_gff": str(gff) if gff.exists() else "",
            "genomic_fna": str(fna_files[0]) if fna_files else ""
        })

df = pd.DataFrame(rows)
Path("results/tables").mkdir(parents=True, exist_ok=True)
df.to_csv("results/tables/downloaded_genomes_index.tsv", sep="\t", index=False)

print(df.head())
print(f"\nTotal genomes indexed: {len(df)}")
print("\nMissing genomic_fna:", (df["genomic_fna"] == "").sum())
print("Missing protein_faa:", (df["protein_faa"] == "").sum())
print("Missing genomic_gff:", (df["genomic_gff"] == "").sum())
