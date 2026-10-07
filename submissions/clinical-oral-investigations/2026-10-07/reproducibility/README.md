# Reproduzir sem WSL
Python 3.12; instalar requirements.txt em ambiente separado. Descompactar data/Primary_inputs.zip em inputs/ e supplementary/ESM_3.zip em galaxy_export/. Não são necessários novos alinhamentos para reproduzir esta análise.

```
python reproducibility/analyze_hof.py --source inputs --out results_reproduced
python reproducibility/make_figures.py --results results_reproduced
python reproducibility/audit_galaxy.py --galaxy galaxy_export/galaxy --reference tables --genome-index inputs/results/tables/downloaded_genomes_index.tsv --out galaxy_audit_reproduced
```

O primeiro comando reconstrói resultados a partir dos outputs arquivados; não executa BLAST ou AMRFinderPlus. Os cinco novos conjuntos do Galaxy foram comparados por multiconjuntos, preservando duplicatas e excluindo caminhos locais. Não foram recuperados FASTA originais dos bancos. Datas/quantidades correspondentes não demonstram identidade dos arquivos de referência. AMRFinderPlus e Victors são exclusivamente arquivados. O inventário completo e hashes acompanham o pacote.
