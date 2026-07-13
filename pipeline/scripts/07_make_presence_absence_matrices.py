from pathlib import Path
import pandas as pd

outdir = Path("results/matrices")
outdir.mkdir(parents=True, exist_ok=True)

idx = pd.read_csv("results/tables/downloaded_genomes_index.tsv", sep="\t")
idx["genome_id"] = idx["species_folder"] + "|" + idx["accession"]

# AMR matrix
amr_path = Path("results/tables/amrfinder_merged.tsv")
if amr_path.exists():
    amr = pd.read_csv(amr_path, sep="\t")
    if not amr.empty and "Element symbol" in amr.columns:
        amr["genome_id"] = amr["species_folder"] + "|" + amr["accession"]
        amr_mat = pd.crosstab(amr["genome_id"], amr["Element symbol"])
        amr_mat = (amr_mat > 0).astype(int)
        amr_mat.to_csv(outdir / "matrix_genome_by_amr_gene.tsv", sep="\t")

# VFDB matrix
vfdb_path = Path("results/tables/vfdb_merged.tsv")
if vfdb_path.exists():
    vfdb = pd.read_csv(vfdb_path, sep="\t")
    if not vfdb.empty and "GENE" in vfdb.columns:
        vfdb["genome_id"] = vfdb["species_folder"] + "|" + vfdb["accession"]
        vfdb_mat = pd.crosstab(vfdb["genome_id"], vfdb["GENE"])
        vfdb_mat = (vfdb_mat > 0).astype(int)
        vfdb_mat.to_csv(outdir / "matrix_genome_by_vfdb_gene.tsv", sep="\t")

# Biofilm matrix
bio_path = Path("results/tables/biofilm_adhesion_screen_merged.tsv")
if bio_path.exists():
    bio = pd.read_csv(bio_path, sep="\t")
    if not bio.empty:
        bio["genome_id"] = bio["species_folder"] + "|" + bio["accession"]
        bio_mat = pd.crosstab(bio["genome_id"], bio["category"])
        bio_mat = (bio_mat > 0).astype(int)
        bio_mat.to_csv(outdir / "matrix_genome_by_biofilm_category.tsv", sep="\t")

print("[INFO] Presence/absence matrices saved in results/matrices/")
