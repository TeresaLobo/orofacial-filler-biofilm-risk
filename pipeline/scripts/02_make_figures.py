from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

root=Path.cwd()
figdir=root/'results'/'figures'
figdir.mkdir(parents=True, exist_ok=True)
for name, title in [('amr_presence_absence.tsv','AMR gene richness'),('virulence_presence_absence.tsv','Virulence gene richness'),('biofilm_presence_absence.tsv','Biofilm gene richness')]:
    path=root/'results'/'matrices'/name
    if not path.exists():
        continue
    df=pd.read_csv(path, sep='\t', index_col=0)
    if df.empty or df.shape[1] == 0:
        continue
    richness=df.sum(axis=1).sort_values(ascending=False).head(30)
    plt.figure(figsize=(10, max(4, len(richness)*0.22)))
    richness.sort_values().plot(kind='barh')
    plt.title(title)
    plt.xlabel('Number of detected genes')
    plt.tight_layout()
    plt.savefig(figdir/(name.replace('.tsv','_richness.png')), dpi=300)
    plt.savefig(figdir/(name.replace('.tsv','_richness.pdf')))
    plt.close()
