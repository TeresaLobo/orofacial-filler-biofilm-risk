from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

figdir = Path("results/figures")
figdir.mkdir(parents=True, exist_ok=True)

# Summary integrated
summary_path = Path("results/tables/final_species_integrated_summary.tsv")

# Build integrated summary
idx = pd.read_csv("results/tables/downloaded_genomes_index.tsv", sep="\t")
genome_counts = idx.groupby("species_folder").size().reset_index(name="n_genomes")

def read_summary(path, count_col):
    p = Path(path)
    if p.exists():
        df = pd.read_csv(p, sep="\t")
        return df
    return pd.DataFrame(columns=["species_folder", count_col])

amr = read_summary("results/tables/amrfinder_species_summary.tsv", "amr_hit_count")
vfdb = read_summary("results/tables/vfdb_species_summary.tsv", "vfdb_hit_count")
bio = read_summary("results/tables/biofilm_adhesion_species_summary.tsv", "biofilm_adhesion_hit_count")

summary = genome_counts.merge(amr, on="species_folder", how="left")
summary = summary.merge(vfdb, on="species_folder", how="left")
summary = summary.merge(bio, on="species_folder", how="left")

for col in ["amr_hit_count", "vfdb_hit_count", "biofilm_adhesion_hit_count"]:
    if col not in summary.columns:
        summary[col] = 0
    summary[col] = summary[col].fillna(0).astype(int)

summary["amr_hits_per_genome"] = summary["amr_hit_count"] / summary["n_genomes"]
summary["vfdb_hits_per_genome"] = summary["vfdb_hit_count"] / summary["n_genomes"]
summary["biofilm_adhesion_hits_per_genome"] = summary["biofilm_adhesion_hit_count"] / summary["n_genomes"]

summary["species_label"] = summary["species_folder"].str.replace("_", " ", regex=False).str.capitalize()

summary.to_csv("results/tables/final_species_integrated_summary.tsv", sep="\t", index=False)

# Barplots
for metric, title, xlabel, fname in [
    ("amr_hits_per_genome", "Normalized antimicrobial resistance burden", "AMR hits per genome", "Figure_AMR_normalized.png"),
    ("vfdb_hits_per_genome", "Normalized virulence-associated burden", "VFDB hits per genome", "Figure_VFDB_normalized.png"),
    ("biofilm_adhesion_hits_per_genome", "Normalized biofilm/adhesion burden", "Biofilm/adhesion hits per genome", "Figure_Biofilm_Adhesion_normalized.png"),
]:
    df = summary.sort_values(metric, ascending=True)
    plt.figure(figsize=(8, 6))
    plt.barh(df["species_label"], df[metric])
    plt.xlabel(xlabel)
    plt.ylabel("Species")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(figdir / fname, dpi=300)
    plt.savefig(figdir / fname.replace(".png", ".pdf"))
    plt.close()

# Heatmap integrated by species
heat = summary.set_index("species_label")[
    ["amr_hits_per_genome", "vfdb_hits_per_genome", "biofilm_adhesion_hits_per_genome"]
]

# z-score by column for visualization
heat_z = heat.copy()
for col in heat_z.columns:
    std = heat_z[col].std()
    if std and std > 0:
        heat_z[col] = (heat_z[col] - heat_z[col].mean()) / std
    else:
        heat_z[col] = 0

plt.figure(figsize=(7, 7))
plt.imshow(heat_z.values, aspect="auto")
plt.xticks(range(len(heat_z.columns)), ["AMR", "Virulence", "Biofilm/adhesion"], rotation=45, ha="right")
plt.yticks(range(len(heat_z.index)), heat_z.index)
plt.colorbar(label="Z-score")
plt.title("Integrated genomic risk profile by species")
plt.tight_layout()
plt.savefig(figdir / "Figure_Integrated_species_heatmap.png", dpi=300)
plt.savefig(figdir / "Figure_Integrated_species_heatmap.pdf")
plt.close()

print("[INFO] Integrated summary and figures saved.")
print(summary.sort_values(["amr_hits_per_genome", "vfdb_hits_per_genome", "biofilm_adhesion_hits_per_genome"], ascending=False))
