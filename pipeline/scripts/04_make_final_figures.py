from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

Path("results/figures").mkdir(parents=True, exist_ok=True)

summary = pd.read_csv("results/tables/final_species_amr_vfdb_summary.tsv", sep="\t")

# Melhorar nomes das espécies
summary["species_label"] = (
    summary["species_folder"]
    .str.replace("_", " ", regex=False)
    .str.capitalize()
)

# Figura 1: AMR hits por espécie
df = summary.sort_values("amr_hit_count", ascending=True)

plt.figure(figsize=(8, 6))
plt.barh(df["species_label"], df["amr_hit_count"])
plt.xlabel("AMR hits")
plt.ylabel("Species")
plt.title("Antimicrobial resistance gene hits by species")
plt.tight_layout()
plt.savefig("results/figures/Figure_AMR_hits_by_species.png", dpi=300)
plt.savefig("results/figures/Figure_AMR_hits_by_species.pdf")
plt.close()

# Figura 2: VFDB hits por espécie
df = summary.sort_values("vfdb_hit_count", ascending=True)

plt.figure(figsize=(8, 6))
plt.barh(df["species_label"], df["vfdb_hit_count"])
plt.xlabel("VFDB hits")
plt.ylabel("Species")
plt.title("Virulence-associated gene hits by species")
plt.tight_layout()
plt.savefig("results/figures/Figure_VFDB_hits_by_species.png", dpi=300)
plt.savefig("results/figures/Figure_VFDB_hits_by_species.pdf")
plt.close()

# Figura 3: AMR hits por genoma
df = summary.sort_values("amr_hits_per_genome", ascending=True)

plt.figure(figsize=(8, 6))
plt.barh(df["species_label"], df["amr_hits_per_genome"])
plt.xlabel("AMR hits per genome")
plt.ylabel("Species")
plt.title("Normalized antimicrobial resistance burden")
plt.tight_layout()
plt.savefig("results/figures/Figure_AMR_hits_per_genome.png", dpi=300)
plt.savefig("results/figures/Figure_AMR_hits_per_genome.pdf")
plt.close()

# Figura 4: VFDB hits por genoma
df = summary.sort_values("vfdb_hits_per_genome", ascending=True)

plt.figure(figsize=(8, 6))
plt.barh(df["species_label"], df["vfdb_hits_per_genome"])
plt.xlabel("VFDB hits per genome")
plt.ylabel("Species")
plt.title("Normalized virulence-associated gene burden")
plt.tight_layout()
plt.savefig("results/figures/Figure_VFDB_hits_per_genome.png", dpi=300)
plt.savefig("results/figures/Figure_VFDB_hits_per_genome.pdf")
plt.close()

# Tabela limpa para manuscrito
clean = summary.copy()
clean["Species"] = clean["species_label"]
clean = clean[
    [
        "Species",
        "n_genomes",
        "amr_hit_count",
        "vfdb_hit_count",
        "amr_hits_per_genome",
        "vfdb_hits_per_genome",
    ]
]
clean.columns = [
    "Species",
    "Number of genomes",
    "AMR hits",
    "VFDB hits",
    "AMR hits/genome",
    "VFDB hits/genome",
]
clean.to_csv("results/tables/Table_species_AMR_VFDB_for_manuscript.tsv", sep="\t", index=False)

print("[INFO] Figures saved in results/figures/")
print("[INFO] Manuscript table saved: results/tables/Table_species_AMR_VFDB_for_manuscript.tsv")
print(clean)
