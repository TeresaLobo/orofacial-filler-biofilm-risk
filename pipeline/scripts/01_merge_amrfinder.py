from pathlib import Path
import pandas as pd

amr_dir = Path("results/amr")
files = sorted(amr_dir.glob("*_amrfinder.tsv"))

rows = []

for f in files:
    name = f.name.replace("_amrfinder.tsv", "")

    if "_GCF_" in name:
        species, acc_part = name.split("_GCF_", 1)
        accession = "GCF_" + acc_part
    else:
        species = name
        accession = ""

    try:
        df = pd.read_csv(f, sep="\t")
    except Exception as e:
        print(f"[WARN] Could not read {f}: {e}")
        continue

    if df.empty:
        continue

    # Alguns arquivos podem conter apenas cabeçalho, sem hits
    if len(df.columns) == 0:
        continue

    df.insert(0, "species_folder", species)
    df.insert(1, "accession", accession)
    df.insert(2, "source_file", str(f))
    rows.append(df)

Path("results/tables").mkdir(parents=True, exist_ok=True)

if rows:
    merged = pd.concat(rows, ignore_index=True)
    merged.to_csv("results/tables/amrfinder_merged.tsv", sep="\t", index=False)

    print(f"[INFO] AMRFinder files found: {len(files)}")
    print(f"[INFO] AMRFinder files with hits: {len(rows)}")
    print(f"[INFO] Total AMR hits: {len(merged)}")
    print("[INFO] Saved: results/tables/amrfinder_merged.tsv")

    # Detecta automaticamente a coluna de gene
    possible_gene_cols = ["Gene symbol", "Element symbol", "Gene", "Element name"]
    gene_col = next((c for c in possible_gene_cols if c in merged.columns), None)

    if gene_col:
        summary = (
            merged.groupby(["species_folder", gene_col])
            .size()
            .reset_index(name="count")
            .sort_values(["species_folder", "count"], ascending=[True, False])
        )
        summary.to_csv("results/tables/amrfinder_gene_summary.tsv", sep="\t", index=False)
        print(f"[INFO] Saved: results/tables/amrfinder_gene_summary.tsv using column: {gene_col}")

        species_summary = (
            merged.groupby("species_folder")
            .size()
            .reset_index(name="amr_hit_count")
            .sort_values("amr_hit_count", ascending=False)
        )
        species_summary.to_csv("results/tables/amrfinder_species_summary.tsv", sep="\t", index=False)
        print("[INFO] Saved: results/tables/amrfinder_species_summary.tsv")
    else:
        print("[WARN] No recognized gene column found. Merged table was still saved.")
        print("[INFO] Available columns:")
        print(list(merged.columns))

else:
    empty = pd.DataFrame(columns=["species_folder", "accession", "source_file"])
    empty.to_csv("results/tables/amrfinder_merged.tsv", sep="\t", index=False)

    print(f"[INFO] AMRFinder files found: {len(files)}")
    print("[WARN] No AMR hits found in the AMRFinder outputs.")
    print("[INFO] Empty merged table saved: results/tables/amrfinder_merged.tsv")
