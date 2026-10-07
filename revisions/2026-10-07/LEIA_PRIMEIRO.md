# HOF revisado para revisão autoral

Foi preparado um manuscrito de análise genômica exploratória para o Brazilian Journal of Oral Sciences, com resultados corrigidos, tabelas, figuras, suplemento, página de título, carta e formulário oficial de ciência aberta. O material está em revisão e não foi submetido nem publicado.

## O que mudou

| Etapa | Resultado |
|---|---|
| Biofilme e adesão | A busca deixou de pesquisar sequências de aminoácidos. Foram encontrados 1.009 loci candidatos em campos de anotação, com 543 loci de evidência por símbolo nomeado. Não são genes de biofilme experimentalmente validados. |
| Erro antigo | A busca antiga produziu 106.727 ocorrências, incluindo 67.614 em sequências proteicas. Ocorrências antigas e loci novos são unidades diferentes. |
| Seleção | Confirmadas as primeiras entradas retornadas pelas consultas arquivadas para as 17 espécies. Amostra de conveniência, sem seleção aleatória ou restrição a isolados humanos. |
| QC | 82 genomas completos com tamanho concordante; 66 passam o filtro mais estrito post hoc. Prevotella intermedia não tem representante nesse subconjunto. |
| AMRFinderPlus | 168 registros separados em 158 AMR e 10 STRESS. Os 158 registros AMR ocorrem em 41 genomas. |
| ABRicate | Tabelas reconstruídas dos arquivos individuais, preservando o cabeçalho #FILE. |
| BacMet2 | 7.412 registros arquivados, dos quais 192 atendem explicitamente a identidade ≥80% e cobertura ≥70%. Foram excluídos 7.220 registros. |
| Sensibilidade | Resultados por banco, ponderação por domínios, retirada de bancos, QC, normalização por Mbp, reamostragem, PCA e associações regenerados. O score permanece descritivo. |
| Manuscrito | 97 caracteres no título, 210 palavras no resumo, 20 referências, duas tabelas, quatro figuras principais e uma suplementar. DOCX em A4, Arial 12, espaçamento 1,5 e margens de 3 cm. |
| Referências | Verificação de metadados e retirada de referência com DOI incorreto. O DOI antigo da referência 7 apontava para outro artigo. |

## Arquivos para ler

- `manuscript/HOF_BJOS_Manuscrito_REVISAR.docx` e PDF: artigo revisado.
- `HOF_Tabelas_Corrigidas.xlsx`: consulta das contagens, candidatos, QC e sensibilidade em oito abas.
- `tables/`: dados completos em CSV/TSV, inclusive registros individuais corrigidos dos bancos.
- `figures/`: TIFF sem perda a 400 dpi, PDF e PNG das cinco figuras.
- `manuscript/HOF_BJOS_Suplemento.docx`: métodos suplementares e guia dos dados.
- `manuscript/HOF_BJOS_Pagina_Titulo_CONFIRMAR.docx`: dados autorais provisórios e declaração de assistência de IA.
- `manuscript/HOF_BJOS_Carta_REVISAR.docx`: carta preparada para revisão.
- `manuscript/BJOS_Formulario_Ciencia_Aberta_PREENCHER.docx`: formulário oficial com ausência de preprint e ambas as opções de revisão aberta confirmadas pela autora; disponibilidade integral pendente.
- `reproducibility/`: scripts portáveis e versões necessárias.
- `provenance/`: hashes, logs auditados, seleção, referências, versões e verificações públicas.
- `GALAXY.md`: caminho para novas análises sem WSL.

O ZIP completo inclui os arquivos FNA/FAA/GFF dos 82 genomas, metadados, consultas, acessos, 574 resultados e 574 logs originais utilizados na análise. O ZIP leve contém a revisão e os documentos, sem esses inputs grandes. O backup integral anterior continua preservado à parte.

## O que falta antes de submeter

1. Revisão científica pela autora: aceitar o enquadramento como estudo genômico exploratório, conferir os candidatos relevantes e as limitações de seleção/QC. A origem dos genomas não é limitada a casos com preenchedores.
2. Resolver a reprodutibilidade dos bancos: versões de ferramentas e timestamps dos bancos foram recuperados, mas os FASTA originais dos bancos não. Para repetir exatamente os alinhamentos, recuperar os snapshots. Se isso não for possível, uma nova rodada no Galaxy deve usar bancos identificados e arquivados, com revisão dos resultados. A revisão atual reconstrói e filtra outputs existentes.
3. Para sustentar alegações específicas de determinantes de biofilme, acrescentar validação por homologia/domínios. O manuscrito atual restringe suas conclusões a anotações candidatas; ensaios experimentais seriam necessários para alegações fenotípicas ou clínicas.
4. Confirmar autora, afiliação, contribuições, financiamento e conflitos; completar endereço postal e telefone. O “sim” recebido autoriza manter os dados existentes provisoriamente, não fornece os campos ausentes nem substitui aprovação final.
5. Confirmar revisão profissional do inglês, originalidade e ausência de submissão simultânea, conforme exigências do periódico.
6. Preprint: não depositado. Revisão aberta: publicação dos pareceres e interação direta aceitas; ambas registradas no formulário oficial. Finalizar a seção de disponibilidade de dados após o depósito corrigido.
7. Publicar uma nova versão do repositório/Zenodo com esta revisão, conferir licença e DOI, e atualizar a declaração de disponibilidade. A versão pública v1.0.4/DOI 10.5281/zenodo.21339092 é anterior às correções e não deve ser citada como arquivo dos resultados novos.
8. Aprovação autoral final e envio ao BJOS, com conferência das exigências no dia da submissão.

Não foram encontrados dados de ensaios experimentais neste projeto. Não foram fabricados dados, medições, aprovação ética, revisão profissional, declarações autorais ou resultados executados no Galaxy.

## Fontes verificadas em 7 outubro 2026

- [Diretrizes e templates do BJOS](https://periodicos.sbu.unicamp.br/ojs/index.php/bjos/about/submissions).
- [Repositório original](https://github.com/TeresaLobo/orofacial-filler-biofilm-risk).
- [Arquivo público original no Zenodo](https://zenodo.org/records/21339092).
- [Ferramentas Galaxy verificadas por API](https://usegalaxy.eu/api/tools?in_panel=false).

As páginas do artigo, página de título, carta e suplemento foram renderizadas e inspecionadas; as figuras e planilhas também foram conferidas. Os documentos que dependem de declaração humana continuam claramente identificados para revisão.

Atualização autoral: e-mail de correspondência teresa.lobo@esenfar.ufal.br. Endereço postal e telefone ainda não foram fornecidos.

Repositório atualizado localmente em uma branch de revisão. Execução completa do script portátil: 27 tabelas idênticas às revisadas; 1.471 inputs verificados por SHA-256. Publicação no GitHub ainda depende de autenticação. Galaxy: histórico HOF exclusivo criado, ZIP com 82 genomas transferido; extração em execução e inventário dos bancos enfileirado. Nenhum novo alinhamento foi integrado.
