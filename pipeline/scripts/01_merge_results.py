from pathlib import Path
import pandas as pd

root = Path.cwd()
out = root / 'results' / 'matrices'
out.mkdir(parents=True, exist_ok=True)

# Merge AMRFinder outputs
amr_rows=[]
for path in (root/'results'/'amr').glob('*_amrfinder.tsv'):
    sample = path.name.replace('_amrfinder.tsv','')
    try:
        df = pd.read_csv(path, sep='\t')
        if df.empty:
            continue
        df.insert(0, 'sample', sample)
        amr_rows.append(df)
    except Exception as e:
        print(f'[WARN] Could not parse {path}: {e}')
if amr_rows:
    amr=pd.concat(amr_rows, ignore_index=True)
    amr.to_csv(out/'amrfinder_merged.tsv', sep='\t', index=False)
    gene_col='Gene symbol' if 'Gene symbol' in amr.columns else ('Element symbol' if 'Element symbol' in amr.columns else None)
    if gene_col:
        mat=(amr.assign(present=1).pivot_table(index='sample', columns=gene_col, values='present', aggfunc='max', fill_value=0))
        mat.to_csv(out/'amr_presence_absence.tsv', sep='\t')
else:
    pd.DataFrame(columns=['sample','gene','present']).to_csv(out/'amr_presence_absence.tsv', sep='\t', index=False)

# Merge VFDB outputs from ABRicate
vf_rows=[]
for path in (root/'results'/'virulence').glob('*_vfdb.tsv'):
    sample=path.name.replace('_vfdb.tsv','')
    try:
        df=pd.read_csv(path, sep='\t', comment='#')
        if df.empty: continue
        df.insert(0,'sample',sample)
        vf_rows.append(df)
    except Exception as e:
        print(f'[WARN] Could not parse {path}: {e}')
if vf_rows:
    vf=pd.concat(vf_rows, ignore_index=True)
    vf.to_csv(out/'vfdb_merged.tsv', sep='\t', index=False)
    gene_col='GENE' if 'GENE' in vf.columns else None
    if gene_col:
        mat=(vf.assign(present=1).pivot_table(index='sample', columns=gene_col, values='present', aggfunc='max', fill_value=0))
        mat.to_csv(out/'virulence_presence_absence.tsv', sep='\t')
else:
    pd.DataFrame(columns=['sample','gene','present']).to_csv(out/'virulence_presence_absence.tsv', sep='\t', index=False)

# Merge DIAMOND biofilm hits
bio_rows=[]
for path in (root/'results'/'biofilm').glob('*_biofilm.tsv'):
    sample=path.name.replace('_biofilm.tsv','')
    try:
        df=pd.read_csv(path, sep='\t', header=None)
        if df.empty: continue
        df.columns=['query','subject','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore'][:df.shape[1]]
        df.insert(0,'sample',sample)
        bio_rows.append(df)
    except Exception as e:
        print(f'[WARN] Could not parse {path}: {e}')
if bio_rows:
    bio=pd.concat(bio_rows, ignore_index=True)
    bio.to_csv(out/'biofilm_merged.tsv', sep='\t', index=False)
    mat=(bio.assign(present=1).pivot_table(index='sample', columns='subject', values='present', aggfunc='max', fill_value=0))
    mat.to_csv(out/'biofilm_presence_absence.tsv', sep='\t')
else:
    pd.DataFrame(columns=['sample','gene','present']).to_csv(out/'biofilm_presence_absence.tsv', sep='\t', index=False)
