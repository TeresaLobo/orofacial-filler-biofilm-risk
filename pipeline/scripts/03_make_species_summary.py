from pathlib import Path
import pandas as pd

Path("results/tables").mkdir(parents=True, exist_ok=True)

# Genome counts
idx = pd.read_csv("results/tables/downloaded_genomes_index.tsv", sep="\t")
genome_counts = (
    idx.groupby("species_folder")
    .size()
    .reset_index(name="n_genomes")
)

# AMR counts
amr_path = Path("results/tables/amrfinder_species_summary.tsv")
if amr_path.exists():
    amr = pd.read_csv(amr_path, sep="\t")
else:
    amr = pd.DataFrame(columns=["species_folder", "amr_hit_count"])

# VFDB counts
vfdb_path = Path("results/tables/vfdb_species_summary.tsv")
if vfdb_path.exists():
    vfdb = pd.read_csv(vfdb_path, sep="\t")
else:
    vfdb = pd.DataFrame(columns=["species_folder", "vfdb_hit_count"])

summary = genome_counts.merge(amr, on="species_folder", how="left")
summary = summary.merge(vfdb, on="species_folder", how="left")

summary["amr_hit_count"] = summary["amr_hit_count"].fillna(0).astype(int)
summary["vfdb_hit_count"] = summary["vfdb_hit_count"].fillna(0).astype(int)

summary["amr_hits_per_genome"] = summary["amr_hit_count"] / summary["n_genomes"]
summary["vfdb_hits_per_genome"] = summary["vfdb_hit_count"] / summary["n_genomes"]

summary = summary.sort_values(
    ["amr_hit_count", "vfdb_hit_count"],
    ascending=False
)

summary.to_csv("results/tables/final_species_amr_vfdb_summary.tsv", sep="\t", index=False)

print(summary)
print("\n[INFO] Saved: results/tables/final_species_amr_vfdb_summary.tsv")
