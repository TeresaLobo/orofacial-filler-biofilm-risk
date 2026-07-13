from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

figdir = Path("results/figures")
figdir.mkdir(parents=True, exist_ok=True)

matrices = {
    "AMR": "results/matrices/matrix_genome_by_amr_gene.tsv",
    "VFDB": "results/matrices/matrix_genome_by_vfdb_gene.tsv",
    "Biofilm": "results/matrices/matrix_genome_by_biofilm_category.tsv",
}

for label, path in matrices.items():
    p = Path(path)
    if not p.exists():
        print(f"[WARN] Missing matrix: {path}")
        continue

    mat = pd.read_csv(p, sep="\t", index_col=0)

    if mat.shape[0] < 3 or mat.shape[1] < 2:
        print(f"[WARN] Matrix too small for PCA: {label} {mat.shape}")
        continue

    X = mat.values
    Xs = StandardScaler(with_mean=True, with_std=True).fit_transform(X)

    pca = PCA(n_components=2)
    coords = pca.fit_transform(Xs)

    plot_df = pd.DataFrame({
        "genome_id": mat.index,
        "PC1": coords[:, 0],
        "PC2": coords[:, 1],
    })
    plot_df["species_folder"] = plot_df["genome_id"].str.split("|").str[0]
    plot_df.to_csv(f"results/tables/pca_{label.lower()}_coordinates.tsv", sep="\t", index=False)

    plt.figure(figsize=(7, 6))
    for species, sub in plot_df.groupby("species_folder"):
        plt.scatter(sub["PC1"], sub["PC2"], label=species.replace("_", " "), s=35)

    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    plt.title(f"PCA based on {label} presence/absence profile")
    plt.legend(fontsize=6, bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(f"results/figures/Figure_PCA_{label}.png", dpi=300)
    plt.savefig(f"results/figures/Figure_PCA_{label}.pdf")
    plt.close()

    print(f"[INFO] PCA {label}: {mat.shape}; explained variance = {pca.explained_variance_ratio_}")
