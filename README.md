# HOF: resistance and adhesion annotations in bacterial reference genomes

## Current submission preparation — Clinical Oral Investigations

The latest manuscript, editable tables, vector/EPS and 600-dpi TIFF figures, supplementary data, original Galaxy exports and portable code are in [submissions/clinical-oral-investigations/2026-10-07](submissions/clinical-oral-investigations/2026-10-07/LEIA_PRIMEIRO.md).

Five new Galaxy screens completed 410 outputs (82 genomes per database). All filtered alignment multisets exactly match the corrected archive: CARD 214, ResFinder 105, VFDB 172, PlasmidFinder 40 and BacMet2 192. AMRFinderPlus and Victors remain archived-only. Reference FASTA provenance remains incomplete and is disclosed. All original scientific inputs remain under pipeline/; primary-input hashes are included in the submission materials. Unpack ESM_2.zip to obtain individual processed CSV/TSV tables and ESM_3.zip to obtain original Galaxy export ZIPs.

Clinical Oral Investigations offers a subscription route without APC; optional open access is paid. The article is prepared for author review and has not been submitted. The author requested omission of telephone. The old Zenodo DOI describes the earlier deposit and does not represent these updated files.

## Corrections

- The legacy adhesion screen searched amino acid strings and counted repeated annotations. The revised screen evaluates structured GFF annotation fields, uses bounded names and genus restrictions, excludes pseudogenes and deduplicates loci: 1,009 candidate loci, including 543 named-symbol candidates. These are not experimentally validated biofilm determinants.
- The 82 assemblies cover 17 species selected as the first available entries of archived RefSeq complete-genome queries. Selection is by convenience. Sixty-six assemblies pass the post hoc strict CheckM and sequence-integrity filter; Prevotella intermedia has no retained assembly in that subset.
- AMRFinderPlus outputs contain 158 AMR and 10 STRESS records, reported separately.
- Individual ABRicate files were reconstructed with their actual #FILE header. Explicit BacMet2 protein identity >=80% and coverage >=70% filtering retains 192 of 7,412 archived alignments.
- Species summaries, rank sensitivity, bootstrap summaries, PCA and association edges were regenerated. Rankings are descriptive and cannot predict clinical outcomes.

## Read the corrected materials

[Manuscript](revisions/2026-10-07/manuscript/HOF_BJOS_Manuscrito_REVISAR.pdf), [tables](revisions/2026-10-07/tables/), [figures](revisions/2026-10-07/figures/), [workbook](revisions/2026-10-07/HOF_Tabelas_Corrigidas.xlsx), [provenance](revisions/2026-10-07/provenance/) and [Galaxy guidance](revisions/2026-10-07/GALAXY.md).

## Reconstruct the corrected outputs without WSL

Python 3.12 and the versions in the requirements file are sufficient for the archived-output reconstruction:

```sh
python -m pip install -r revisions/2026-10-07/reproducibility/requirements.txt
python revisions/2026-10-07/reproducibility/analyze_hof.py --source pipeline --out results_reproduced
python revisions/2026-10-07/reproducibility/make_figures.py --results results_reproduced
```

The archived inputs and 574 individual result/log pairs are preserved under pipeline/. Input hashes are in the corrected provenance directory. This command reconstructs and filters existing alignments; it does not execute fresh BLAST searches. Original reference database FASTA snapshots were not recovered, so exact alignment reruns remain unresolved. Tool versions and local database timestamps are documented; a timestamp is not a verified database release identifier.

The historical pipeline scripts and historical derived files are retained for audit. Do not use their biofilm counts or clinical-risk interpretation in a new submission. Use the corrected scripts and outputs linked above.

## Publication status

This revision is a draft for author review, not a submitted or accepted article. Correspondence email, no preprint deposit and acceptance of both open-review options were confirmed by the author. Postal address, telephone, professional English review, final scientific approval and the complete research-data availability declaration remain pending.

The earlier archive [Zenodo record 21339092](https://zenodo.org/records/21339092), DOI 10.5281/zenodo.21339092, predates these corrections. It must not be presented as an archive of the corrected results. This corrected review revision is publicly available on GitHub. A new Zenodo archival deposit and DOI still require verification. See [publication status](PUBLICATION_STATUS.md).

## Licensing

Code: MIT. Author-generated text, figures and tables: CC BY 4.0 under DATA_LICENSE.md. Source genomic data and external reference databases retain their original terms. No reference database redistribution rights or experimental findings are implied.
