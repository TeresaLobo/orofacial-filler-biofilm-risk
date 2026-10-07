from pathlib import Path
import os,sys,json
os.environ.setdefault('MPLCONFIGDIR',str(Path('.mplconfig').resolve()))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np,pandas as pd
import argparse
parser=argparse.ArgumentParser();parser.add_argument('--results',type=Path,default=Path('results_reproduced'));args=parser.parse_args()
OUT=args.results;T=OUT/'tables';F=OUT/'figures';F.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
def label(x):return x.replace('_',' ').capitalize()
def save(fig,name):
 title='Resistance and adhesion annotations in 82 bacterial genomes relevant to orofacial filler research'
 fig.savefig(F/(name+'.tiff'),dpi=400,format='tiff',pil_kwargs={'compression':'tiff_lzw','tiffinfo':{270:title}})
 fig.savefig(F/(name+'.png'),dpi=250)
 fig.savefig(F/(name+'.pdf'),metadata={'Title':title,'Subject':name})
 plt.close(fig)
s=pd.read_csv(T/'Table2_species_database_and_annotation_summary.csv',index_col=0)
fig,axs=plt.subplots(1,2,figsize=(6.3,4.2),layout='constrained')
a=pd.read_csv(T/'S4_legacy_regex_audit_all_genomes.csv');v=[a.legacy_total_matches.sum(),a.sequence_only_matches.sum()]
axs[0].bar(['All legacy\ntext matches','Sequence-only\nlegacy matches'],v,color=['#727f8e','#b25040']);axs[0].set_ylabel('Regex occurrences (not loci)');axs[0].set_title('A  Legacy screen audit',loc='left')
for i,x in enumerate(v):axs[0].text(i,x+1800,f'{x:,}',ha='center',fontsize=8)
axs[0].set_ylim(0,120000);axs[0].ticklabel_format(axis='y',style='plain')
hits=pd.read_csv(T/'S3_corrected_annotation_candidates.csv');keys=['accession','contig','locus_tag','protein_id'];v=[len(hits.drop_duplicates(keys)),len(hits[hits.evidence_tier.eq('named_symbol')].drop_duplicates(keys))];axs[1].bar(['All candidate\nloci','Named-symbol\ncandidate loci'],v,color=['#267c84','#164951']);axs[1].set_ylabel('Distinct annotated loci');axs[1].set_title('B  Corrected annotation screen',loc='left')
for i,x in enumerate(v):axs[1].text(i,x+18,str(x),ha='center');axs[1].set_ylim(0,1150)
save(fig,'Figure1_screen_audit')

mat=pd.read_csv(T/'S7_annotation_category_presence_matrix.csv',index_col=0);idx=pd.read_csv(T/'S8_integrated_per_genome.csv').set_index('accession')
freq=mat.join(idx.species_folder).groupby('species_folder').mean();freq=freq.loc[:,freq.sum().gt(0)]
fig,axs=plt.subplots(2,1,figsize=(6.3,7),layout='constrained')
for i,ax in enumerate(axs):
 block=freq.iloc[:,i*15:(i+1)*15];im=ax.imshow(block,vmin=0,vmax=1,cmap='YlGnBu',aspect='auto')
 ax.set_yticks(range(len(block)),[label(x) for x in block.index],fontsize=6.5);ax.set_xticks(range(len(block.columns)),[x.replace('_',' ') for x in block.columns],rotation=90,fontsize=6)
 ax.set_title(f'{chr(65+i)}  Candidate annotation categories {i*15+1}–{i*15+len(block.columns)}',loc='left',fontsize=8)
 fig.colorbar(im,ax=ax,shrink=.7,label='Fraction of sampled genomes')
save(fig,'Figure2_annotation_categories')

fig,axs=plt.subplots(1,2,figsize=(6.3,4.9),layout='constrained')
order=s.sort_values('AMR_records',ascending=True).index
axs[0].barh(range(17),s.loc[order,'AMR_records'],color='#267c84');axs[0].set_yticks(range(17),[label(x) for x in order],fontsize=6.5);axs[0].set_xlabel('AMRFinderPlus AMR records');axs[0].set_title('A  Antibiotic resistance annotations',loc='left',fontsize=8)
c=pd.read_csv(T/'S6_database_records_per_genome.csv');d=c.groupby('database')[['archived_records','total_records']].sum().loc['bacmet2']
axs[1].bar(['Archived\nprotein matches','After ≥80% identity\nand ≥70% coverage'],d.values,color=['#b25040','#267c84']);axs[1].set_title('B  BacMet2 threshold correction',loc='left',fontsize=8);axs[1].set_ylabel('Alignment records')
for i,x in enumerate(d):axs[1].text(i,x+120,f'{x:,}',ha='center',fontsize=7)
axs[1].set_ylim(0,8200);save(fig,'Figure3_resistance_records')

r=pd.read_csv(T/'S9_rank_sensitivity.csv',index_col=0);cols=['seven_databases_AMR_only','four_domains_equal_weight','AMRFinder_VFDB_only','nucleotide_databases_only','strict_QC_seven_databases','records_per_Mbp']
order=r.seven_databases_AMR_only.sort_values().index;r=r.loc[order,cols]
fig,axs=plt.subplots(2,1,figsize=(6.3,7),layout='constrained');ax=axs[0];cmap=plt.get_cmap('viridis_r').copy();cmap.set_bad('#ededed');im=ax.imshow(r.to_numpy(),vmin=1,vmax=17,cmap=cmap,aspect='auto')
ax.set_yticks(range(17),[label(x) for x in r.index],fontsize=7);ax.set_xticks(range(6),['Seven databases','Four domains','AMR + VFDB','Exclude BacMet2','Strict QC','Per Mbp'],rotation=45,ha='right',fontsize=7)
for i in range(17):
 for j in range(6):
  val=r.iloc[i,j];ax.text(j,i,'NA' if pd.isna(val) else str(int(val)),ha='center',va='center',fontsize=7,color='black' if pd.isna(val) or val<10 else 'white')
ax.set_title('A  Species ranks depend on the descriptive model',loc='left',fontsize=8);fig.colorbar(im,ax=ax,shrink=.65,label='Rank (1 = largest burden score)')
b=pd.read_csv(T/'S9_within_sample_bootstrap_rank_stability.csv').set_index('species_folder').loc[order];ax=axs[1];ax.errorbar(b.rank_median,range(17),xerr=[b.rank_median-b.rank_p025,b.rank_p975-b.rank_median],fmt='o',color='#267c84',markersize=3,capsize=2,linewidth=1)
ax.set_yticks(range(17),[label(x) for x in b.index],fontsize=6.5);ax.invert_yaxis();ax.set_xlim(.5,17.5);ax.set_xticks([1,5,10,15,17]);ax.set_xlabel('Within-sample rank median and 2.5–97.5% quantiles',fontsize=7);ax.set_title('B  Resampling stability in the convenience sample',loc='left',fontsize=8);save(fig,'Figure4_rank_sensitivity')

p=pd.read_csv(T/'S10_corrected_annotation_PCA.csv');ev=pd.read_csv(T/'S10_PCA_explained_variance.csv');fig,axs=plt.subplots(1,2,figsize=(6.3,5.1));ax=axs[0]
for i,(sp,g) in enumerate(p.groupby('species_folder')):ax.scatter(g.PC1,g.PC2,s=20,color=plt.get_cmap('tab20')(i),label=label(sp),edgecolors='white',linewidth=.3)
ax.set_xlabel(f'PC1 ({ev.iloc[0].explained_variance_ratio:.1%})');ax.set_ylabel(f'PC2 ({ev.iloc[1].explained_variance_ratio:.1%})');ax.set_title('A  Corrected annotation PCA',loc='left',fontsize=8)
axs[1].plot(ev.component,np.cumsum(ev.explained_variance_ratio),marker='o',markersize=3,color='#267c84');axs[1].set_xlabel('Principal component');axs[1].set_ylabel('Cumulative explained variance');axs[1].set_ylim(0,1.05);axs[1].set_title('B  Variance retained',loc='left',fontsize=8)
fig.legend(*ax.get_legend_handles_labels(),loc='lower center',ncol=3,fontsize=6,frameon=False);fig.subplots_adjust(left=.12,right=.98,top=.94,bottom=.34,wspace=.55);save(fig,'FigureS1_annotation_PCA')
print('Five figures written as TIFF, PNG and PDF.')
