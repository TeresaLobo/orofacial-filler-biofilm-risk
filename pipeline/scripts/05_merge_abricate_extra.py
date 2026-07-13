from pathlib import Path
import pandas as pd

base = Path("results/abricate_extra")
outdir = Path("results/tables")
outdir.mkdir(parents=True, exist_ok=True)

all_rows = []

for db_dir in sorted(base.iterdir()):
    if not db_dir.is_dir():
        continue

    db = db_dir.name
    files = sorted(db_dir.glob(f"*_{db}.tsv"))
    rows = []

    for f in files:
        name = f.name.replace(f"_{db}.tsv", "")

        if "_GCF_" in name:
            species, acc_part = name.split("_GCF_", 1)
            accession = "GCF_" + acc_part
        else:
            species = name
            accession = ""

        try:
            df = pd.read_csv(f, sep="\t", comment="#")
        except Exception:
            continue

        if df.empty:
            continue

        df.insert(0, "database", db)
        df.insert(1, "species_folder", species)
        df.insert(2, "accession", accession)
        df.insert(3, "source_file", str(f))
        rows.append(df)

    if rows:
        merged = pd.concat(rows, ignore_index=True)
        merged.to_csv(outdir / f"abricate_{db}_merged.tsv", sep="\t", index=False)

        species_summary = (
            merged.groupby(["database", "species_folder"])
            .size()
            .reset_index(name=f"{db}_hit_count")
            .sort_values(f"{db}_hit_count", ascending=False)
        )
        species_summary.to_csv(outdir / f"abricate_{db}_species_summary.tsv", sep="\t", index=False)

        all_rows.append(merged)
        print(f"[INFO] {db}: {len(files)} files; {len(merged)} hits")
    else:
        print(f"[WARN] {db}: no hits")

if all_rows:
    all_merged = pd.concat(all_rows, ignore_index=True)
    all_merged.to_csv(outdir / "abricate_all_extra_merged.tsv", sep="\t", index=False)

    all_summary = (
        all_merged.groupby(["database", "species_folder"])
        .size()
        .reset_index(name="hit_count")
        .sort_values(["database", "hit_count"], ascending=[True, False])
    )
    all_summary.to_csv(outdir / "abricate_all_extra_species_summary.tsv", sep="\t", index=False)

    print("[INFO] Saved combined ABRicate extra results.")
