from pathlib import Path
import pandas as pd
import networkx as nx

outdir = Path("results/networks")
outdir.mkdir(parents=True, exist_ok=True)

edges = []

# AMR edges
amr_path = Path("results/tables/amrfinder_merged.tsv")
if amr_path.exists():
    amr = pd.read_csv(amr_path, sep="\t")
    if not amr.empty and "Element symbol" in amr.columns:
        for _, r in amr.iterrows():
            edges.append({
                "source": r["species_folder"],
                "target": r["Element symbol"],
                "interaction": "AMR",
                "category": r.get("Class", "AMR")
            })

# VFDB edges
vfdb_path = Path("results/tables/vfdb_merged.tsv")
if vfdb_path.exists():
    vfdb = pd.read_csv(vfdb_path, sep="\t")
    if not vfdb.empty and "GENE" in vfdb.columns:
        for _, r in vfdb.iterrows():
            edges.append({
                "source": r["species_folder"],
                "target": r["GENE"],
                "interaction": "Virulence",
                "category": "VFDB"
            })

# Biofilm edges
bio_path = Path("results/tables/biofilm_adhesion_screen_merged.tsv")
if bio_path.exists():
    bio = pd.read_csv(bio_path, sep="\t")
    if not bio.empty:
        for _, r in bio.iterrows():
            edges.append({
                "source": r["species_folder"],
                "target": r["category"],
                "interaction": "Biofilm_adhesion",
                "category": r["description"]
            })

edges_df = pd.DataFrame(edges)
edges_df = (
    edges_df.groupby(["source", "target", "interaction", "category"])
    .size()
    .reset_index(name="weight")
)

edges_df.to_csv(outdir / "species_gene_function_edges.tsv", sep="\t", index=False)

G = nx.from_pandas_edgelist(
    edges_df,
    source="source",
    target="target",
    edge_attr=["interaction", "category", "weight"]
)

nodes = []
for node in G.nodes():
    if node in set(edges_df["source"]):
        ntype = "species"
    else:
        ntype = "gene_or_function"
    nodes.append({
        "id": node,
        "label": node.replace("_", " "),
        "type": ntype,
        "degree": G.degree(node)
    })

nodes_df = pd.DataFrame(nodes)
nodes_df.to_csv(outdir / "species_gene_function_nodes.tsv", sep="\t", index=False)

nx.write_graphml(G, outdir / "species_gene_function_network.graphml")

print("[INFO] Network files saved in results/networks/")
print(f"[INFO] Nodes: {G.number_of_nodes()} Edges: {G.number_of_edges()}")
