# HOF no Galaxy sem WSL

O Galaxy europeu é compatível com os formatos recuperados. Em 7 outubro 2026, sua API listou AMRFinderPlus 4.2.7+galaxy0, ABRicate 1.4.0 e 1.4.0+galaxy1, NCBI BLASTp e CheckM lineage_wf. Essa verificação confirma a presença das ferramentas, não a execução do projeto nem a disponibilidade dos mesmos bancos antigos.

## Preparar os dados

Extraia o ZIP completo. Em `inputs/results/tables/downloaded_genomes_index.tsv`, cada linha liga um acesso aos seus arquivos FNA/FAA/GFF. O manifesto `provenance/archived_input_inventory.csv` permite conferir hashes. Há 82 genomas; não há leituras FASTQ nem ensaios clínicos/experimentais neste conjunto.

Crie um histórico identificado e uma coleção de dados por genoma. Mantenha o acesso GCF no nome dos arquivos e no identificador da coleção. Use FASTA nucleotídico para FNA, FASTA proteico para FAA e GFF3 para GFF. Não concatene genomas como se fossem uma única amostra.

## Nova rodada recomendada

1. **AMRFinderPlus real.** Selecione AMRFinderPlus 4.2.7 quando disponível, confira os modos e use entradas FNA/FAA/GFF emparelhadas conforme a interface instalada. Registre banco, data/versão, opções e logs. A opção de banco NCBI no ABRicate não substitui o algoritmo AMRFinderPlus.
2. **Bancos complementares.** Execute ABRicate sobre cada FNA para CARD, ResFinder, VFDB, Victors e PlasmidFinder quando o servidor os oferecer. Registre a lista de bancos e seus metadados antes de executar. Configure 80% de identidade e 70% de cobertura e confira esses campos nos outputs.
3. **BacMet2 proteico.** Confirme que o servidor oferece BacMet2 e o tipo correto. O código ABRicate 1.4.0 usa BLASTx para banco proteico, e os outputs antigos demonstram identidades inferiores ao limiar declarado. Aplique filtragem explícita aos campos de saída. Se o banco não estiver disponível, use busca proteica com referência BacMet apropriada e arquivada; registre essa mudança de método e refaça a interpretação.
4. **QC independente, se necessário.** CheckM lineage_wf permite recalcular estimativas. Essas estimativas e versões devem ser distinguidas das recuperadas do NCBI. O critério 95%/5% é uma sensibilidade post hoc e não altera retrospectivamente a seleção original.
5. **Adesão/biofilme por homologia.** Use uma coleção de proteínas de referência cuidadosamente curada com acessos, sequências, domínios e evidência experimental documentados. BLASTp pode localizar homólogos, mas a identificação deve avaliar identidade, extensão alinhada nas duas sequências, domínios, arquitetura e contexto. Não há um limiar universal que valide todos os determinantes deste painel. Não classifique automaticamente enolase/sortase como prova de biofilme.
6. Exporte histórico, parâmetros, workflow Galaxy, datasets, bancos de referência e logs. Recalcule as tabelas, os modelos de sensibilidade e todas as figuras com os outputs da nova rodada. Resultados obtidos com bancos novos constituem uma nova análise, não uma reprodução exata da antiga.

## Limite atual

A correção local já funciona em Windows/Python sem WSL e reproduz a interpretação dos arquivos arquivados. Nenhuma tarefa foi enviada ao Galaxy neste trabalho. Para a nova rodada, falta escolher o histórico/conta e confirmar as referências disponíveis. A identidade e o arquivamento dos bancos são a principal pendência técnica para repetir os alinhamentos com rastreabilidade completa.

Fontes primárias: [Galaxy europeu](https://usegalaxy.eu/), [wrapper AMRFinderPlus](https://github.com/galaxyproject/tools-iuc/tree/main/tools/amrfinderplus), [wrapper ABRicate](https://github.com/galaxyproject/tools-iuc/tree/main/tools/abricate), [código ABRicate](https://github.com/tseemann/abricate). Os identificadores verificados estão em `provenance/galaxy_EU_tools_verified.json`.
