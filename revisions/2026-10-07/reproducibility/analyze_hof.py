from pathlib import Path
from urllib.parse import unquote
import json,re,csv,hashlib,collections,datetime
import pandas as pd
import numpy as np

import argparse
parser=argparse.ArgumentParser(description='Reconstruct corrected HOF results from archived inputs; does not repeat BLAST searches.')
parser.add_argument('--source',type=Path,default=Path('inputs'))
parser.add_argument('--out',type=Path,default=Path('results_reproduced'))
args=parser.parse_args()
SRC=args.source
OUT=args.out;TABLE=OUT/'tables';TABLE.mkdir(parents=True,exist_ok=True)
EVID=OUT/'provenance';EVID.mkdir(exist_ok=True)
idx=pd.read_csv(SRC/'results/tables/downloaded_genomes_index.tsv',sep='\t')
oldtargets=pd.read_csv(SRC/'data/curated_targets/biofilm_adhesion_targets.tsv',sep='\t')

def fasta(path):
 records=[];header=None;seq=[]
 for line in Path(path).read_text(encoding='utf-8',errors='strict').splitlines():
  if line.startswith('>'):
   if header is not None:records.append((header,''.join(seq).upper()))
   header=line[1:];seq=[]
  else:seq.append(line.strip())
 if header is not None:records.append((header,''.join(seq).upper()))
 return records

def attrs(text):return {k:unquote(v) for part in text.split(';') if '=' in part for k,v in [part.split('=',1)]}
def parse_gff(path):
 genes={};cds={};comments=[]
 for line in Path(path).read_text(encoding='utf-8',errors='strict').splitlines():
  if line.startswith('#'):comments.append(line);continue
  p=line.split('\t')
  if len(p)!=9:continue
  a=attrs(p[8]);typ=p[2]
  if typ in ['gene','pseudogene']:genes[a.get('ID','')]={'type':typ,**a}
  if typ!='CDS':continue
  key=(p[0],a.get('locus_tag',a.get('Parent',a.get('ID'))),a.get('protein_id',''))
  if key not in cds:cds[key]={'seqid':p[0],'start':int(p[3]),'end':int(p[4]),'strand':p[6],**a}
  else:cds[key]['start']=min(cds[key]['start'],int(p[3]));cds[key]['end']=max(cds[key]['end'],int(p[4]))
 for c in cds.values():
  parent=genes.get(c.get('Parent',''),{})
  c['gene_name']=c.get('gene',parent.get('gene',''))
  c['gene_synonyms']=c.get('gene_synonym',parent.get('gene_synonym',''))
  c['is_pseudo']=c.get('pseudo')=='true' or parent.get('pseudo')=='true' or parent.get('type')=='pseudogene'
 return list(cds.values()),comments

def restriction(cat):
 if cat.startswith('pseudomonas_') and cat!='pseudomonas_pili':return 'Pseudomonas'
 if cat.startswith('streptococcus_') or cat in ['biofilm_brpA','regulation_vicRK','competence_com']:return 'Streptococcus'
 if cat.startswith('enterococcus_'):return 'Enterococcus'
 if cat.startswith('pg_'):return 'Porphyromonas'
 if cat=='fn_adhesin':return 'Fusobacterium'
 if cat=='tf_bspA':return 'Tannerella'
 if cat=='td_virulence':return 'Treponema'
 if cat in ['biofilm_ica','biofilm_accumulation','biofilm_bap','biofilm_sasG','quorum_agr','regulation_sarA','stress_sigB','adhesion_sdr','adhesion_clumping']:return 'Staphylococcus'
 return ''

panel=[]
for t in oldtargets.to_dict('records'):
 cat=t['category'];alternatives=t['target_regex'].split('|');names=[x for x in alternatives if re.fullmatch(r'[A-Za-z0-9]+',x)];phrases=[x for x in alternatives if x not in names]
 # Broad functional labels must not be treated as short gene symbols.
 for word in ['autolysin','glucosyltransferase','sortase','enolase','gelatinase']:
  if word in names:names.remove(word);phrases.append(word)
 if cat=='adhesion_fibronectin' and 'fibronectin' in names:names.remove('fibronectin');phrases.append('fibronectin-binding')
 if cat=='adhesion_collagen':phrases=[x for x in phrases if x!='collagen-binding']+['collagen-binding adhesin']
 if cat=='pg_gingipain':phrases+=['gingipain']
 if cat=='pg_fimbriae':phrases+=['fimbrial subunit protein FimA','minor fimbrial subunit protein Mfa1']
 if cat=='fn_adhesin':phrases+=['FadA family adhesin','Fap2','RadD','FomA']
 if cat=='tf_bspA':phrases+=['BspA']
 if cat=='pseudomonas_pili':cat='type_iv_pili';phrases=['type IV pilus'];names=[]
 panel.append({'category':cat,'gene_names':';'.join(names),'product_phrases':';'.join(sorted(set(phrases))),'restricted_genus':restriction(t['category']),'description':t['description'],'origin':'Original curated target list, corrected field matching and genus restrictions'})
paneldf=pd.DataFrame(panel);paneldf.to_csv(TABLE/'S2_annotation_screen_panel.csv',index=False)

metadata={}
for p in (SRC/'data/unzipped_genomes').glob('*/ncbi_dataset/data/assembly_data_report.jsonl'):
 for line in p.read_text(encoding='utf-8').splitlines():
  if line.strip():
   m=json.loads(line);metadata[m['accession']]=m
qc=[];hits=[];comparison=[];input_manifest=[];selection=[]
for i,row in enumerate(idx.to_dict('records'),1):
 sp=row['species_folder'];acc=row['accession'];genus=sp.split('_')[0].capitalize();m=metadata.get(acc,{})
 records=fasta(SRC/row['genomic_fna']);proteins=fasta(SRC/row['protein_faa']);cds,comments=parse_gff(SRC/row['genomic_gff'])
 a=m.get('assemblyInfo',{});s=m.get('assemblyStats',{});check=m.get('checkmInfo',{});bio=a.get('biosample',{});bsattrs={x['name']:x['value'] for x in bio.get('attributes',[])}
 total=sum(len(seq) for _,seq in records);N=sum(seq.count('N') for _,seq in records);ACGT=sum(sum(seq.count(b) for b in 'ACGT') for _,seq in records);GC=sum(seq.count('G')+seq.count('C') for _,seq in records)
 lens=sorted([len(seq) for _,seq in records],reverse=True);cum=np.cumsum(lens);n50=lens[int(np.searchsorted(cum,total/2))]
 protein_ids={h.split()[0] for h,_ in proteins};valid=[c for c in cds if not c['is_pseudo'] and c.get('protein_id')]
 digest=hashlib.sha256('\n'.join(sorted(seq for _,seq in records)).encode()).hexdigest()
 qc.append({'species_folder':sp,'species':sp.replace('_',' ').capitalize(),'accession':acc,'taxid':m.get('organism',{}).get('taxId'),'ncbi_organism':m.get('organism',{}).get('organismName'),'assembly_level':a.get('assemblyLevel'),'assembly_status':a.get('assemblyStatus'),'release_date':a.get('releaseDate'),'annotation_release_date':m.get('annotationInfo',{}).get('releaseDate'),'annotation_software':m.get('annotationInfo',{}).get('softwareVersion'),'biosample':bio.get('accession'),'strain':bio.get('strain',bsattrs.get('strain','')),'isolate':bio.get('isolate',bsattrs.get('isolate','')),'isolation_source':bio.get('isolationSource',bsattrs.get('isolation_source','')),'host':bio.get('host',bsattrs.get('host','')),'country':bio.get('geoLocName',bsattrs.get('geo_loc_name','')),'collection_date':bio.get('collectionDate',bsattrs.get('collection_date','')),'sequence_length_bp':total,'metadata_length_bp':int(s.get('totalSequenceLength',0)),'length_match':total==int(s.get('totalSequenceLength',0)),'replicons':len(records),'N_bases':N,'ACGT_bases':ACGT,'ambiguous_bases':total-ACGT,'GC_percent':GC/ACGT*100 if ACGT else None,'N50_bp':n50,'protein_sequences':len(proteins),'GFF_nonpseudo_CDS':len(valid),'CDS_missing_protein_sequence':sum(c['protein_id'] not in protein_ids for c in valid),'checkm_completeness':check.get('completeness'),'checkm_contamination':check.get('contamination'),'checkm_version':check.get('checkmVersion'),'taxonomy_check':m.get('averageNucleotideIdentity',{}).get('taxonomyCheckStatus'),'canonical_sequence_sha256':digest,'NCBI_URL':f'https://www.ncbi.nlm.nih.gov/datasets/genome/{acc}/'})
 for c in valid:
  gene_tokens=set(re.split(r'[,;\s]+',c['gene_name'].lower()+' '+c['gene_synonyms'].lower()))
  product=c.get('product','')
  for t in panel:
   if t['restricted_genus'] and t['restricted_genus']!=genus:continue
   names=t['gene_names'].split(';') if t['gene_names'] else [];phrases=t['product_phrases'].split(';') if t['product_phrases'] else []
   matched_gene=[name for name in names if name.lower() in gene_tokens]
   matched_symbol=[name for name in names if re.search(r'(?<![A-Za-z0-9])'+re.escape(name)+r'(?![A-Za-z0-9])',product,re.I)]
   matched_phrase=[phrase for phrase in phrases if re.search(r'(?<![A-Za-z0-9])'+re.escape(phrase)+r'(?![A-Za-z0-9])',product,re.I)]
   evidence=matched_gene or matched_symbol or matched_phrase
   if not evidence:continue
   hits.append({'species_folder':sp,'accession':acc,'category':t['category'],'contig':c['seqid'],'locus_tag':c.get('locus_tag',''),'protein_id':c['protein_id'],'start':c['start'],'end':c['end'],'strand':c['strand'],'gene_name':c['gene_name'],'product':product,'evidence_tier':'named_symbol' if matched_gene or matched_symbol else 'product_description','matched_terms':';'.join(evidence),'field':'gene_or_synonym' if matched_gene else 'product','restricted_genus':t['restricted_genus']})
 faa=(SRC/row['protein_faa']).read_text(encoding='utf-8');gff=(SRC/row['genomic_gff']).read_text(encoding='utf-8');headers='\n'.join(h for h,_ in proteins);sequence_text='\n'.join(seq for _,seq in proteins)
 for t in oldtargets.to_dict('records'):
  comparison.append({'species_folder':sp,'accession':acc,'original_category':t['category'],'legacy_total_matches':len(re.findall(t['target_regex'],faa+'\n'+gff,re.I)),'sequence_only_matches':len(re.findall(t['target_regex'],sequence_text,re.I)),'text_annotation_matches_without_dedup':len(re.findall(t['target_regex'],headers+'\n'+gff,re.I))})
 for key in ['genomic_fna','protein_faa','genomic_gff']:
  p=SRC/row[key];input_manifest.append({'accession':acc,'kind':key,'relative_path':row[key],'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size})
 if i%10==0:print(f'QC and corrected annotation screen: {i}/82',flush=True)
q=pd.DataFrame(qc);q['posthoc_strict_QC']=(q.length_match & q.CDS_missing_protein_sequence.eq(0) & q.checkm_completeness.ge(95)&q.checkm_contamination.le(5))
q.to_csv(TABLE/'S1_genome_metadata_QC.csv',index=False)
pd.DataFrame(input_manifest).to_csv(EVID/'input_sha256.csv',index=False)
h=pd.DataFrame(hits);h.to_csv(TABLE/'S3_corrected_annotation_candidates.csv',index=False)
pd.DataFrame(comparison).to_csv(TABLE/'S4_legacy_regex_audit_all_genomes.csv',index=False)
dup=q[q.duplicated('canonical_sequence_sha256',keep=False)];dup.to_csv(TABLE/'S5_exact_sequence_duplicates.csv',index=False)
strain=q[q.duplicated(['species_folder','strain'],keep=False)&q.strain.notna()&q.strain.ne('')];strain.to_csv(TABLE/'S5_repeated_strain_labels.csv',index=False)
for sp,g in idx.groupby('species_folder',sort=False):
 accessionfile=SRC/'data/accessions'/f'{sp}_accessions.txt';summ=SRC/'logs'/f'{sp}_summary.jsonl'
 listed=accessionfile.read_text(encoding='utf-8').splitlines() if accessionfile.exists() else []
 returned=[json.loads(line).get('accession') for line in summ.read_text(encoding='utf-8').splitlines() if line.startswith('{')] if summ.exists() else []
 selection.append({'species_folder':sp,'selected_count':len(g),'selected_accessions':'|'.join(listed),'matches_first_five_returned':listed==returned[:5],'summary_returned':len(returned),'matches_index_set':set(listed)==set(g.accession)})
pd.DataFrame(selection).to_csv(EVID/'accession_selection_audit.csv',index=False)

# Verify every archived per-genome analysis, not merely merged files.
runrows=[];counts=[];merged={}
for r in idx.to_dict('records'):
 sp=r['species_folder'];acc=r['accession']
 for db in ['amrfinder','vfdb','card','resfinder','bacmet2','victors','plasmidfinder']:
  log=SRC/'logs'/f'{sp}_{acc}_{db}.log';txt=log.read_text(encoding='utf-8') if log.exists() else ''
  if db=='amrfinder':p=SRC/'results/amr'/f'{sp}_{acc}_amrfinder.tsv'
  elif db=='vfdb':p=SRC/'results/virulence'/f'{sp}_{acc}_vfdb.tsv'
  else:p=SRC/'results/abricate_extra'/db/f'{sp}_{acc}_{db}.tsv'
  d=pd.read_csv(p,sep='\t') if p.exists() else pd.DataFrame()
  archived_count=len(d)
  if db!='amrfinder' and len(d):
   d=d[d['%IDENTITY'].ge(80)&d['%COVERAGE'].ge(70)].copy()
  merged.setdefault(db,[]).append(d.assign(species_folder=sp,accession=acc))
  info=re.search(r'Using (\w+) database (\S+):\s+(\d+) sequences\s+-\s+(.+)',txt)
  ver=re.search(r'Software version:\s+(\S+)',txt);dbver=re.search(r'Database version:\s+(\S+)',txt)
  complete='amrfinder took' in txt if db=='amrfinder' else 'Done.' in txt
  runrows.append({'species_folder':sp,'accession':acc,'database':db,'output_present':p.exists(),'rows':len(d),'log_complete':complete,'db_type':info.group(1) if info else None,'db_sequence_count':int(info.group(3)) if info else None,'db_timestamp':info.group(4).strip() if info else dbver.group(1) if dbver else None,'software_version':ver.group(1) if ver else '1.4.0 (environment manifest)','min_identity':None if db=='amrfinder' else 80,'min_coverage':None if db=='amrfinder' else 70,'minimum_identity_observed':float(d['%IDENTITY'].min()) if len(d) and '%IDENTITY' in d else None,'minimum_coverage_observed':float(d['%COVERAGE'].min()) if len(d) and '%COVERAGE' in d else None,'output_sha256':hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else None,'log_sha256':hashlib.sha256(log.read_bytes()).hexdigest() if log.exists() else None})
  counts.append({'species_folder':sp,'accession':acc,'database':db,'archived_records':archived_count,'total_records':len(d),'removed_below_identity_or_coverage':archived_count-len(d),'AMR_records':int(d.Type.eq('AMR').sum()) if 'Type' in d else None,'STRESS_records':int(d.Type.eq('STRESS').sum()) if 'Type' in d else None})
for db,frames in merged.items():pd.concat(frames,ignore_index=True).to_csv(TABLE/f'S6_{db}_corrected_records.tsv',sep='\t',index=False)
rr=pd.DataFrame(runrows);rr['rows_definition']='after explicit filtering; hashes identify archived files';rr.to_csv(EVID/'analysis_run_audit.csv',index=False);pd.DataFrame(counts).to_csv(TABLE/'S6_database_records_per_genome.csv',index=False)
cat_order=list(paneldf.category)
matrix=pd.DataFrame(0,index=idx.accession,columns=cat_order,dtype=int)
named=matrix.copy()
for (acc,cat),g in h.groupby(['accession','category']):
 matrix.loc[acc,cat]=int(g.protein_id.nunique()>0);named.loc[acc,cat]=int((g.evidence_tier=='named_symbol').any())
matrix.to_csv(TABLE/'S7_annotation_category_presence_matrix.csv',index_label='accession');named.to_csv(TABLE/'S7_named_symbol_presence_matrix.csv',index_label='accession')
summary=q[['species_folder','accession','sequence_length_bp','checkm_completeness','checkm_contamination','posthoc_strict_QC']].set_index('accession').join(pd.DataFrame(counts).pivot(index='accession',columns='database',values='total_records')).join(pd.DataFrame(counts).query("database=='amrfinder'").set_index('accession')[['AMR_records','STRESS_records']])
annotation_total=h.groupby('accession').apply(lambda x:len(x.drop_duplicates(['contig','locus_tag','protein_id'])),include_groups=False)
summary['annotation_loci']=annotation_total.reindex(summary.index).fillna(0).astype(int)
summary['named_symbol_loci']=h[h.evidence_tier=='named_symbol'].groupby('accession').apply(lambda x:len(x.drop_duplicates(['contig','locus_tag','protein_id'])),include_groups=False).reindex(summary.index).fillna(0).astype(int)
summary.to_csv(TABLE/'S8_integrated_per_genome.csv',index_label='accession')
species=summary.groupby('species_folder').agg(n_genomes=('AMR_records','size'),AMR_records=('AMR_records','sum'),STRESS_records=('STRESS_records','sum'),AMR_genomes=('AMR_records',lambda x:x.gt(0).sum()),annotation_loci=('annotation_loci','sum'),named_symbol_loci=('named_symbol_loci','sum'),sequence_length_mean=('sequence_length_bp','mean'),strict_QC_genomes=('posthoc_strict_QC','sum'))
for col in ['AMR_records','STRESS_records','annotation_loci','named_symbol_loci','vfdb','card','resfinder','bacmet2','victors','plasmidfinder']:
 if col not in species:species[col]=summary.groupby('species_folder')[col].sum()
 species[col+'_mean']=summary.groupby('species_folder')[col].mean();species[col+'_min']=summary.groupby('species_folder')[col].min();species[col+'_max']=summary.groupby('species_folder')[col].max()
species.to_csv(TABLE/'Table2_species_database_and_annotation_summary.csv',index_label='species_folder')

# Scores are descriptive; independent biological domains are not assumed.
features=['AMR_records','card','resfinder','bacmet2','vfdb','victors','plasmidfinder']
def standardize(frame):
 std=frame.std(ddof=1);return (frame-frame.mean()).div(std.replace(0,np.nan),axis=1).fillna(0)
burden=summary.groupby('species_folder')[features].mean();zs=standardize(burden)
scores=pd.DataFrame(index=burden.index);scores['seven_databases_AMR_only']=zs.mean(axis=1)
orig=pd.DataFrame(counts).pivot(index='accession',columns='database',values='archived_records').join(q.set_index('accession').species_folder).groupby('species_folder').mean().rename(columns={'amrfinder':'AMR_records'})[features];scores['archived_raw_seven_databases']=standardize(orig).mean(axis=1)
scores['AMRFinder_VFDB_only']=zs[['AMR_records','vfdb']].mean(axis=1)
domains=pd.DataFrame({'AMR':zs[['AMR_records','card','resfinder']].mean(axis=1),'virulence':zs[['vfdb','victors']].mean(axis=1),'biocide_metal':zs.bacmet2,'plasmid':zs.plasmidfinder})
scores['four_domains_equal_weight']=domains.mean(axis=1)
scores['nucleotide_databases_only']=zs[[x for x in features if x!='bacmet2']].mean(axis=1)
for col in features:scores['leave_out_'+col]=zs.drop(columns=col).mean(axis=1)
strict=summary[summary.posthoc_strict_QC].groupby('species_folder')[features].mean().reindex(burden.index)
scores['strict_QC_seven_databases']=standardize(strict).mean(axis=1).mask(strict.isna().all(axis=1))
sizeadjusted=summary[features].div(summary.sequence_length_bp/1e6,axis=0).groupby(summary.species_folder).mean();scores['records_per_Mbp']=standardize(sizeadjusted).mean(axis=1)
# Remove repeated strain labels deterministically, preserving missing labels separately.
tmp=summary.reset_index().merge(q[['accession','strain']],on='accession');tmp['strain_key']=tmp.strain.fillna('').mask(tmp.strain.fillna('').isin(['','missing','not applicable']),tmp.accession)
unique=tmp.sort_values('accession').drop_duplicates(['species_folder','strain_key']).groupby('species_folder')[features].mean().reindex(burden.index)
scores['one_accession_per_strain_label']=standardize(unique).mean(axis=1)
scores.to_csv(TABLE/'S9_score_sensitivity.csv',index_label='species_folder')
ranks=scores.rank(ascending=False,method='min');ranks.to_csv(TABLE/'S9_rank_sensitivity.csv',index_label='species_folder')
zs.to_csv(TABLE/'S9_standardized_burdens.csv',index_label='species_folder')
burden.rank().corr().to_csv(TABLE/'S9_database_spearman_correlations.csv',index_label='database')

# Stratified resampling describes stability in this convenience sample only.
rng=np.random.default_rng(20261007);B=1000;bootranks=[]
groups=[summary[summary.species_folder==sp][features].to_numpy(float) for sp in burden.index]
for _ in range(B):
 means=np.array([g[rng.integers(0,len(g),len(g))].mean(axis=0) for g in groups]);z=standardize(pd.DataFrame(means,index=burden.index,columns=features));sc=z.mean(axis=1);bootranks.append(sc.rank(ascending=False,method='min').to_numpy())
br=np.array(bootranks);pd.DataFrame({'species_folder':burden.index,'rank_median':np.median(br,axis=0),'rank_p025':np.quantile(br,.025,axis=0),'rank_p975':np.quantile(br,.975,axis=0),'fraction_top5':(br<=5).mean(axis=0)}).to_csv(TABLE/'S9_within_sample_bootstrap_rank_stability.csv',index=False)

# PCA from corrected binary matrix, no independent clinical interpretation.
X=matrix.to_numpy(float);nonconstant=X.std(axis=0)>0;X=X[:,nonconstant];X-=X.mean(axis=0);U,S,Vt=np.linalg.svd(X,full_matrices=False);ev=S*S/(len(X)-1);explained=ev/ev.sum() if ev.sum() else ev
pca=pd.DataFrame({'accession':idx.accession,'species_folder':idx.species_folder,'PC1':U[:,0]*S[0],'PC2':U[:,1]*S[1]});pca.to_csv(TABLE/'S10_corrected_annotation_PCA.csv',index=False)
pd.DataFrame({'component':np.arange(1,len(explained)+1),'explained_variance_ratio':explained}).to_csv(TABLE/'S10_PCA_explained_variance.csv',index=False)
pd.DataFrame(Vt[:2].T,index=np.array(cat_order)[nonconstant],columns=['PC1_loading','PC2_loading']).to_csv(TABLE/'S10_PCA_loadings.csv',index_label='category')
network=[]
for sp,g in idx.groupby('species_folder'):
 for cat in matrix:
  n=int(matrix.loc[g.accession,cat].sum())
  if n:network.append({'source':sp,'target':cat,'n_genomes_with_annotation':n,'n_genomes':len(g),'within_species_fraction':n/len(g),'interpretation':'annotation-category association; not microbial interaction'})
pd.DataFrame(network).to_csv(TABLE/'S11_corrected_annotation_network_edges.csv',index=False)

comp=pd.DataFrame(comparison)
report={'genomes':len(q),'species':q.species_folder.nunique(),'assembly_levels':q.assembly_level.value_counts().to_dict(),'sequence_lengths_all_match':bool(q.length_match.all()),'ambiguous_bases':int(q.ambiguous_bases.sum()),'missing_CDS_protein_sequences':int(q.CDS_missing_protein_sequence.sum()),'checkm_available':int(q.checkm_completeness.notna().sum()),'checkm_completeness_range':[q.checkm_completeness.min(),q.checkm_completeness.max()],'checkm_contamination_range':[q.checkm_contamination.min(),q.checkm_contamination.max()],'strict_QC_pass':int(q.posthoc_strict_QC.sum()),'strict_QC_failed_accessions':q.loc[~q.posthoc_strict_QC,'accession'].tolist(),'exact_duplicate_genomes':len(dup),'repeated_strain_labels_genomes':len(strain),'legacy_regex_total':int(comp.legacy_total_matches.sum()),'legacy_sequence_only_matches':int(comp.sequence_only_matches.sum()),'corrected_candidate_category_rows':len(h),'corrected_unique_loci':int(summary.annotation_loci.sum()),'named_symbol_unique_loci':int(summary.named_symbol_loci.sum()),'corrected_positive_genomes':int(summary.annotation_loci.gt(0).sum()),'corrected_positive_categories':int(matrix.sum().gt(0).sum()),'AMR_records':int(summary.AMR_records.sum()),'STRESS_records':int(summary.STRESS_records.sum()),'AMR_positive_genomes':int(summary.AMR_records.gt(0).sum()),'analysis_outputs':len(rr),'analysis_logs_complete':int(rr.log_complete.sum()),'database_details':rr[['database','db_type','db_sequence_count','db_timestamp','software_version']].drop_duplicates().fillna('').to_dict('records'),'PC_variance':explained[:2].tolist(),'score_top5':scores.seven_databases_AMR_only.nlargest(5).to_dict(),'original_top5':scores.archived_raw_seven_databases.nlargest(5).to_dict(),'rank_ranges':{sp:[float(ranks.loc[sp].min()),float(ranks.loc[sp].max())] for sp in ranks.index},'candidate_evidence_tiers':h.evidence_tier.value_counts().to_dict(),'selected_first5_confirmed_species':sum(x['matches_first_five_returned'] for x in selection)}
(OUT/'analysis_summary.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(report,ensure_ascii=False,indent=2),flush=True)
