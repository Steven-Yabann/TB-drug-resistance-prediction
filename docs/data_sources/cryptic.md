# CRyPTIC release_june2022

Source of genotypes and phenotypes for the per-drug matrices. A copy of this file is
kept as `data/raw/cryptic/MANIFEST.md`; `data/raw/cryptic/download_log.json` and
`SHA256SUMS` are written by `src/tb_resistance/download.py`.

| | |
|---|---|
| Provider | CRyPTIC Consortium, hosted by EMBL-EBI |
| Source URL | <ftp://ftp.ebi.ac.uk/pub/databases/cryptic/release_june2022/> |
| Tables used | `reproducibility/data_tables/cryptic-analysis-group/` (CRyPTIC "Release 1.1", tables dated 2020–2022) |
| Downloaded | 2026-10-29 (UTC timestamps per file below) |
| Reference genome | *M. tuberculosis* H37Rv, **NC_000962.3** (CRyPTIC Clockwork v0.8.3 pipeline; `GENOMES.CLOCKWORK_VERSION`) |
| Phenotype method | UKMYC5 / UKMYC6 broth-microdilution plates; MIC → binary call by ECOFF/ECV |
| Isolate key | `UNIQUEID` (`site.SITEID.subj.SUBJID.lab.LABID.iso.ISOLATENO`) |
| Access | Open, anonymous FTP; no registration or data-access agreement |

## Licence and citation

The release does not include a licence file and its READMEs state no licence terms.
They ask users who publish with the data to cite the compendium paper:

> The CRyPTIC Consortium. "A data compendium associating the genomes of 12,289
> *Mycobacterium tuberculosis* isolates with quantitative resistance phenotypes to
> 13 antibiotics." *PLoS Biology* 20(8): e3001721 (2022).
> Preprint: "A data compendium of *Mycobacterium tuberculosis* antibiotic
> resistance", bioRxiv, <https://doi.org/10.1101/2021.09.14.460274>.

Phenotype binarisation (ECOFFs): The CRyPTIC Consortium, "Epidemiological cut-off
values for a 96-well broth microdilution plate for high-throughput research
antibiotic susceptibility testing of *M. tuberculosis*", medRxiv,
<https://doi.org/10.1101/2021.02.24.21252386>.

**To verify before publication:** the licence terms with the consortium or EBI, and
the final journal citation details above.

## Files

| File | FTP path (under the release) | Bytes | SHA-256 | Downloaded (UTC) | Used for |
|---|---|---:|---|---|---|
| `UKMYC_PHENOTYPES.csv.gz` | `reproducibility/data_tables/cryptic-analysis-group/` | 4,029,662 | `c03b0d6ce7f9b97b9942ca341729ece4fbf3a5fad5b3596806d7f7a1f3ac0bc0` | 2026-10-29T06:22:14 | Labels |
| `MUTATIONS.csv.gz` | `reproducibility/data_tables/cryptic-analysis-group/` | 1,486,015,433 | `9bb79071c016f9fd3af40fb8f86852fdceb8ec185235fe1321031f27aada5cc5` | 2026-10-29T06:47:52 | Features |
| `GENOMES.csv.gz` | `reproducibility/data_tables/cryptic-analysis-group/` | 5,215,282 | `0ea44f01c835260f85f1cde7930d5d0608915433b6faeabe92132a6a71980be4` | 2026-10-29T06:48:06 | Genotyped-isolate list, lineage |
| `SAMPLES.csv.gz` | `reproducibility/data_tables/cryptic-analysis-group/` | 250,429 | `9958f419edeb334381726ba183dffbb3bf9aa07202a3df21780c343e20ce90d9` | 2026-10-29T06:48:09 | Not joined (no `UNIQUEID`); country metadata for later |
| `CRyPTIC_reuse_table_20240917.csv` | `reuse/` | 5,393,811 | `f785ae87f2bfb7f39f08d2ca2618a1925a290c39bfdfd5a535995b1827541c1b` | 2026-10-29T06:48:16 | Cross-check of labels |
| `README.md` | `reproducibility/data_tables/cryptic-analysis-group/` | 6,760 | `31a4862f49fad1a2bdda4c8ce05874d9dd8e9e40721b4fa1c8b76e12ae4e29dc` | 2026-10-29T06:48:17 | Documentation |
| `CHANGELOG.md` | `reproducibility/data_tables/cryptic-analysis-group/` | 4,115 | `cc7cc99f676f38945d00359efd223c0826f4201eed521a5c3eb149754c7f1b61` | 2026-10-29T06:48:19 | Documentation |
| `DRUG_CODES.csv` | `reproducibility/data_tables/cryptic-analysis-group/` | 687 | `3093ff576b4cb81fe2ab168a92f6cdfcc68999554c045d172f5dbc3c547878f9` | 2026-10-29T06:48:20 | Documentation |
| `DATA_SCHEMA.pdf` | `reproducibility/data_tables/cryptic-analysis-group/` | 103,678 | `9f49bfc12dd24e3d6947d6de637596d08411eb43daf85e94495238d757e4963f` | 2026-10-29T06:48:23 | Documentation |
| `README` | `reuse/` | 5,441 | `f111eb7607509b0a7f201aea2f2f827fa441163303d97fe0722b4fc7215ba746` | 2026-10-29T06:48:24 | Documentation (ECOFF table) |
| `CHANGELOG` | `reuse/` | 864 | `2db878ea09c9f48d123d8ec1806c5945a53fb1ee31e13c454fbb9eabe5f8a93b` | 2026-10-29T06:48:26 | Documentation |

Not downloaded: `VARIANTS.csv.gz` (3.0 GB; per-variant depth and read support,
not needed because `MUTATIONS` carries the het/null/filter flags), the VCFs, FASTQs
and the Clockwork container image.

## Notes on the release

- `UKMYC_PHENOTYPES` has 228,650 rows for 17,050 isolates and 14 drugs (AMI, BDQ,
  CFZ, DLM, EMB, ETH, INH, KAN, LEV, LZD, MXF, PAS, RFB, RIF). **There is no PZA**:
  the UKMYC plates do not test pyrazinamide.
- Each isolate has exactly one reading per drug. Reading day is 14 for most calls,
  10 or 21 for the rest.
- Late site-04 and all site-13 plates were withheld by the consortium for plate
  quality (`CHANGELOG.md` §1).
- Never modify files in `data/raw/cryptic/`; derived data goes to
  `data/interim/cryptic/` and `data/processed/`.
