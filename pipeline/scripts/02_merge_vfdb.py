from pathlib import Path
import pandas as pd

vfdb_dir = Path("results/virulence")
files = sorted(vfdb_dir.glob("*_vfdb.tsv"))
rows = []

for f in files:
    name = f.name.replace("_vfdb.tsv", "")

    if "_GCF_" in name:
        species, acc_part = name.split("_GCF_", 1)
        accession = "GCF_" + acc_part
    else:
        species = name
        accession = ""

    try:
        df = pd.read_csv(f, sep="\t", comment="#")
    except Exception as e:
        print(f"[WARN] Could not read {f}: {e}")
        continue

    if df.empty:
        continue

    df.insert(0, "species_folder", species)
    df.insert(1, "accession", accession)
    df.insert(2, "source_file", str(f))
    rows.append(df)

Path("results/tables").mkdir(parents=True, exist_ok=True)

if rows:
    merged = pd.concat(rows, ignore_index=True)
    merged.to_csv("results/tables/vfdb_merged.tsv", sep="\t", index=False)

    print(f"[INFO] VFDB files found: {len(files)}")
    print(f"[INFO] VFDB files with hits: {len(rows)}")
    print(f"[INFO] Total VFDB hits: {len(merged)}")

    if "GENE" in merged.columns:
        gene_summary = (
            merged.groupby(["species_folder", "GENE"])
            .size()
            .reset_index(name="count")
            .sort_values(["species_folder", "count"], ascending=[True, False])
        )
        gene_summary.to_csv("results/tables/vfdb_gene_summary.tsv", sep="\t", index=False)

    species_summary = (
        merged.groupby("species_folder")
        .size()
        .reset_index(name="vfdb_hit_count")
        .sort_values("vfdb_hit_count", ascending=False)
    )
    species_summary.to_csv("results/tables/vfdb_species_summary.tsv", sep="\t", index=False)

else:
    pd.DataFrame(columns=["species_folder", "accession", "source_file"]).to_csv(
        "results/tables/vfdb_merged.tsv", sep="\t", index=False
    )
    print("[WARN] No VFDB hits found.")
