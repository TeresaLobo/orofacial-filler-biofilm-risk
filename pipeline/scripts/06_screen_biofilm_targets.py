from pathlib import Path
import pandas as pd
import re

idx = pd.read_csv("results/tables/downloaded_genomes_index.tsv", sep="\t")
targets = pd.read_csv("data/curated_targets/biofilm_adhesion_targets.tsv", sep="\t")

rows = []

def read_text(path):
    try:
        return Path(path).read_text(errors="ignore")
    except Exception:
        return ""

for _, genome in idx.iterrows():
    species = genome["species_folder"]
    accession = genome["accession"]
    faa = genome["protein_faa"]
    gff = genome["genomic_gff"]

    faa_text = read_text(faa)
    gff_text = read_text(gff)
    combined = faa_text + "\n" + gff_text

    for _, target in targets.iterrows():
        pattern = str(target["target_regex"])
        category = target["category"]
        description = target["description"]

        matches = re.findall(pattern, combined, flags=re.IGNORECASE)

        if matches:
            rows.append({
                "species_folder": species,
                "accession": accession,
                "category": category,
                "target_regex": pattern,
                "description": description,
                "match_count": len(matches)
            })

outdir = Path("results/tables")
outdir.mkdir(parents=True, exist_ok=True)

bio = pd.DataFrame(rows)

if bio.empty:
    bio = pd.DataFrame(columns=[
        "species_folder", "accession", "category", "target_regex",
        "description", "match_count"
    ])

bio.to_csv(outdir / "biofilm_adhesion_screen_merged.tsv", sep="\t", index=False)

species_summary = (
    bio.groupby(["species_folder", "category"])
    ["match_count"]
    .sum()
    .reset_index()
    .sort_values(["species_folder", "match_count"], ascending=[True, False])
)

species_summary.to_csv(outdir / "biofilm_adhesion_species_category_summary.tsv", sep="\t", index=False)

species_total = (
    bio.groupby("species_folder")["match_count"]
    .sum()
    .reset_index(name="biofilm_adhesion_hit_count")
    .sort_values("biofilm_adhesion_hit_count", ascending=False)
)

species_total.to_csv(outdir / "biofilm_adhesion_species_summary.tsv", sep="\t", index=False)

print("[INFO] Biofilm/adhesion screening finished.")
print(f"[INFO] Total hits: {bio['match_count'].sum() if not bio.empty else 0}")
print(species_total.head(20))
