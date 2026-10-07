from pathlib import Path
import zipfile,io,re,json,hashlib
import pandas as pd
import argparse
p=argparse.ArgumentParser();p.add_argument('--galaxy',type=Path,required=True);p.add_argument('--reference',type=Path,required=True);p.add_argument('--genome-index',type=Path,required=True);p.add_argument('--out',type=Path,required=True);args=p.parse_args()
ROOT=args.galaxy.parent; OLD=args.reference.parent
OUT=args.out; OUT.mkdir(parents=True,exist_ok=True)
idx=pd.read_csv(args.genome_index,sep='\t').set_index('accession')
report=[]; counts=[]
for path in sorted(args.galaxy.glob('*.zip')):
 db={'CARD':'card','ResFinder':'resfinder','VFDB':'vfdb','PlasmidFinder':'plasmidfinder','BacMet':'bacmet2'}[path.stem]
 frames=[]; entries=[]
 with zipfile.ZipFile(path) as z:
  assert z.testzip() is None
  for name in z.namelist():
   if name.endswith('/'):continue
   acc=re.search(r'GCF_\d+\.\d+',name).group(); assert acc in idx.index
   raw=z.read(name); d=pd.read_csv(io.BytesIO(raw),sep='\t'); n=len(d)
   assert {'%IDENTITY','%COVERAGE','GENE','DATABASE'}.issubset(d.columns)
   d=d[d['%IDENTITY'].ge(80)&d['%COVERAGE'].ge(70)].copy()
   frames.append(d.assign(species_folder=idx.loc[acc,'species_folder'],accession=acc))
   counts.append(dict(database=db,accession=acc,species_folder=idx.loc[acc,'species_folder'],raw_records=n,filtered_records=len(d)))
   entries.append(dict(accession=acc,member=name,sha256=hashlib.sha256(raw).hexdigest(),raw_records=n,filtered_records=len(d)))
  assert len(entries)==82 and len({e['accession'] for e in entries})==82
 fresh=pd.concat(frames,ignore_index=True); fresh.to_csv(OUT/f'S6_{db}_new_Galaxy_records.tsv',sep='\t',index=False)
 old=pd.read_csv(args.reference/f'S6_{db}_corrected_records.tsv',sep='\t')
 cols=[c for c in old.columns if c not in ['#FILE','FILE'] and c in fresh.columns]
 def signatures(df):return df[cols].fillna('').astype(str).apply(tuple,axis=1).value_counts().to_dict()
 a,b=signatures(old),signatures(fresh)
 added=sum(max(n-a.get(k,0),0) for k,n in b.items()); removed=sum(max(n-b.get(k,0),0) for k,n in a.items())
 report.append(dict(database=db,genomes=82,raw_records=sum(e['raw_records'] for e in entries),filtered_records=len(fresh),archived_filtered_records=len(old),added_records=added,removed_records=removed,identical_records=(a==b),zip_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),members=entries))
pd.DataFrame(counts).to_csv(OUT/'Galaxy_records_per_genome.csv',index=False)
(OUT/'Galaxy_comparison.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps([{k:v for k,v in r.items() if k!='members'} for r in report],indent=2))
