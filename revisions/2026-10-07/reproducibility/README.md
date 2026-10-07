# Reproduzir a análise corrigida de HOF

O pacote completo contém `inputs/` com o snapshot recuperado. O pacote leve contém resultados, código e documentos e depende dos inputs para executar novamente.

Em Python 3.12, crie um ambiente separado e instale as versões de `requirements.txt`. Não é necessário WSL. A partir da raiz extraída do pacote:

```
python reproducibility/analyze_hof.py --source inputs --out results_reproduced
python reproducibility/make_figures.py --results results_reproduced
```

O primeiro comando audita os genomas, corrige a busca em anotações e reconstrói as contagens a partir dos resultados individuais. Ele não executa novamente AMRFinderPlus ou BLAST. As sequências originais dos bancos não foram recuperadas. As versões antigas são documentadas por logs e manifests; hashes permitem verificar a integridade dos arquivos recuperados.

As figuras têm TIFF com compressão sem perda a 400 dpi, PDF vetorial e PNG. Os arquivos CSV/TSV preservam cabeçalhos, resultados por genoma e as relações entre contagens, loci e categorias. As anotações candidatas não constituem validação de biofilme ou risco clínico.

Critérios corrigidos: busca somente em campos estruturados; símbolos delimitados; categorias específicas por gênero; deduplicação de CDS; separação AMR/STRESS; filtragem explícita de identidade/cobertura; sensibilidade a QC, tamanho de genoma e ponderação. O índice de score é exploratório.

Ver `../LEIA_PRIMEIRO.md` e `../GALAXY.md` para limites e próximos passos. Os documentos marcados REVISAR/CONFIRMAR/PREENCHER requerem revisão autoral antes de qualquer submissão. Nenhum experimento foi criado ou pressuposto.
