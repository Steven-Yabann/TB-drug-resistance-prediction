# Data dictionary: CRyPTIC per-drug matrices

Produced by `scripts/build_matrices.py` from `configs/data.yaml`. Source provenance
is in [`data_sources/cryptic.md`](data_sources/cryptic.md).

Coordinates refer to *M. tuberculosis* H37Rv, **NC_000962.3**, as annotated by the
CRyPTIC pipeline (Clockwork v0.8.3 / gumpy-style mutation names).

## Isolate identifier

All outputs are indexed by **`UNIQUEID`**, e.g. `site.02.subj.0958.lab.22A197.iso.1`
(`site.SITEID.subj.SUBJID.lab.LABID.iso.ISOLATENO`). It is the key shared by
`UKMYC_PHENOTYPES`, `MUTATIONS` and `GENOMES`, so no ID translation is needed.
`SAMPLES` has no `UNIQUEID` and is not used in the join.

Rows in every file of a drug directory are in the same order (sorted by `UNIQUEID`).

## Outputs

```text
data/interim/cryptic/
  panel_mutations.parquet   long table: one row per panel-gene MUTATIONS record
  panel_mutations.json      cache key (panel config + source file size)
data/processed/
  build_report.json         counts for every drop/join step
  <drug>/                   rifampicin, isoniazid, levofloxacin, moxifloxacin
    features.parquet        isolate x feature presence (uint8)
    features_missing.parquet isolate x feature missing flag (uint8)
    features_het.parquet    isolate x feature heterozygous-call flag (uint8)
    labels.csv              binary phenotype
    lineage.csv             Mykrobe lineage
    feature_metadata.csv    one row per feature column
```

### `features.parquet`

| Column | Type | Meaning |
|---|---|---|
| `UNIQUEID` (index) | str | Isolate ID |
| `<GENE>@<MUTATION>` (one per feature) | uint8 | 1 = mutation called present; 0 = not called present |

Feature IDs are the source `GENE` and `MUTATION` fields joined with `@`, otherwise
verbatim, e.g. `rpoB@S450L` (amino-acid change), `fabG1@c-15t` (promoter nucleotide
change, negative = upstream of the start codon), `katG@1000_indel` (indel at
nucleotide 1000 of the gene). The `MUTATION` field alone does not name the gene,
so the prefix is required for uniqueness.

The column set is **the same for every drug**: every panel mutation called present
in at least one genotyped isolate in the release. No frequency filtering is applied
here; columns that are all-zero within a drug's isolates are retained (counted in
`build_report.json`). Any frequency threshold must be fitted on training data only.

### `features_missing.parquet`

Same shape, index and columns as `features.parquet`. 1 = the isolate has a **null**
(`IS_NULL`) or **filter-failed** (`IS_FILTER_PASS == False`) call at the feature's
site, so the 0 in `features.parquet` means "unknown", not "absent".

### `features_het.parquet`

Same shape. 1 = the isolate has a **heterozygous / minor-allele** call
(`IS_HET == True`) at the feature's site. The feature is 0 in `features.parquet`.

### `labels.csv`

| Column | Type | Meaning |
|---|---|---|
| `UNIQUEID` | str | Isolate ID |
| `label` | uint8 | 1 = resistant (R), 0 = susceptible (S) |
| `phenotype_quality` | str | `HIGH` / `MEDIUM` / `LOW` (source `PHENOTYPE_QUALITY`); empty if not recorded. **Not used as a filter.** |
| `reading_day` | str | Plate reading day of the call (`10`, `14` or `21`) |
| `n_calls` | int | Number of agreeing R/S calls for this isolate and drug (1 in this release) |

### `lineage.csv`

| Column | Type | Meaning |
|---|---|---|
| `UNIQUEID` | str | Isolate ID |
| `lineage` | str | `GENOMES.MYKROBE_LINEAGE_NAME_1`, e.g. `Lineage 4`, `Mixed`; empty if not called |
| `sublineage` | str | `GENOMES.MYKROBE_LINEAGE_NAME_2`, e.g. `lineage4.1.2`, `lineage2.2/lineage4` |

### `feature_metadata.csv`

| Column | Meaning |
|---|---|
| `feature_id` | Column name in the feature matrices |
| `gene` | Source `GENE` (`rpoB`, `katG`, `inhA`, `fabG1`, `gyrA`) |
| `mutation` | Source `MUTATION`, verbatim |
| `variant_type` | `snp` (single nucleotide change, or an amino-acid change caused by one), `mnp` (amino-acid change from ≥2 nucleotide changes in one codon), `indel` |
| `site_key` | `<gene>:aa:<codon>` for coding changes, `<gene>:nt:<position>` otherwise. Used to match null/het calls to features |
| `in_promoter` | True if upstream of the gene's start codon |
| `indel_detail` | Source `INDEL_2` for indels (e.g. `1000_del_1`); several indels at one position share one feature column |

## Phenotype definition

Labels come from `UKMYC_PHENOTYPES.BINARY_PHENOTYPE`: the CRyPTIC call made by
comparing the UKMYC5/UKMYC6 broth-microdilution MIC with the epidemiological
cut-off (ECOFF/ECV):

| Drug | Code | S (mg/L) | I (mg/L) | R (mg/L) |
|---|---|---|---|---|
| Rifampicin | RIF | ≤0.5 | – | ≥1.0 |
| Isoniazid | INH | ≤0.1 | 0.2, 0.4 | ≥0.8 |
| Levofloxacin | LEV | ≤1 | – | ≥2 |
| Moxifloxacin | MXF | ≤1 | – | ≥2 |

(ECOFF table from `reuse/README`, copying Table 2 of the CRyPTIC ECOFF paper.)

R → 1, S → 0. These are ECOFF-based (wild-type vs non-wild-type) calls, not WHO
clinical breakpoints.

INH `I` calls (MIC 0.2–0.4 mg/L, 943 rows) are dropped. Note that the curated
`reuse/CRyPTIC_reuse_table_20240917.csv` recodes these as `R` (above the 0.1 mg/L
ECOFF). Where both sources have an R/S call, the labels agree for 100% of isolates
for all four drugs (checked 2026-10-29).

## Handling rules

### Missing data

- **Phenotype missing or non-R/S** (blank, `I`): the row is dropped; counts per value
  are in `build_report.json` → `drugs.<drug>.labels.rows_dropped_values`. `I` occurs
  only for INH (MICs 0.2–0.4 mg/L).
- **Isolate with no panel records** in `MUTATIONS` but present in `GENOMES`: all
  features 0 and no missing flags (wild type at every panel site). `MUTATIONS` lists
  only differences from the reference plus null calls.
- **Null genotype call** (`IS_NULL`, alt base `x`, amino acid `X`, e.g. `katG@E195X`):
  never a feature. Every feature at the same site (codon or nucleotide) is flagged in
  `features_missing.parquet` for that isolate.
- **Lineage not called**: blank in `lineage.csv`.

### Low coverage / failed calls

Positions where the CRyPTIC variant caller could not make a confident call
(low depth or failed filters) appear as `IS_NULL` or `IS_FILTER_PASS == False`
records. Both are treated as missing (above), never as 0 = absent. Positions with no
record at all are assumed reference. Read depth (`DP`, `FRS`) lives only in
`VARIANTS.csv.gz`, which this loader does not download.

### Heteroresistance / mixed calls

Calls with `IS_HET == True` (minor-allele / mixed populations) are **not counted as
present**. The site is flagged in `features_het.parquet` and is not marked missing.
Downstream code can treat het sites as present, absent or missing as an explicit
modelling choice.

**In release_june2022, `MUTATIONS` contains no `IS_HET == True` records in the panel
genes**, so `features_het.parquet` is all zeros, even though `GENOMES.N_HET` is
non-zero genome-wide. Mixed calls at panel sites are therefore not identifiable from
`MUTATIONS`. Detecting heteroresistance would need per-variant read support
(`FRS`) from `VARIANTS.csv.gz` or the VCFs.

### Conflicting phenotypes

If an isolate has both an R and an S call for the same drug, it is dropped for that
drug and counted (`isolates_conflicting`). Agreeing duplicate calls are collapsed;
`phenotype_quality` takes the best quality among them. In release_june2022 each
isolate has exactly one reading per drug, so no conflicts arise in practice.

### Reading day

`labels.reading_day` in the config restricts calls to one plate-reading day. It is
`null` (keep all): most calls are day 14, and isolates read on day 10 or day 21 have
no other reading, so filtering would only drop isolates.

### Duplicates and join losses

- Duplicate `UNIQUEID`s in `GENOMES` are counted and collapsed to the first row.
- Exact duplicate (isolate, feature, call status) records in `MUTATIONS` are counted
  and collapsed.
- Labelled isolates absent from `GENOMES` are dropped and counted
  (`labelled_isolates_lost_in_join`).
- Panel `MUTATIONS` isolates absent from `GENOMES` are counted
  (`panel_mutation_ids_not_in_genomes`) and excluded.

## Gene panel

| Gene | Drug(s) | Included records |
|---|---|---|
| rpoB | RIF | coding + promoter |
| katG | INH | coding + promoter |
| inhA | INH | coding + promoter |
| fabG1 | INH | **promoter only**: the fabG1-inhA operon promoter, e.g. `fabG1@c-15t`, `fabG1@t-8c` |
| gyrA | LEV, MXF | coding + promoter |

Set `panel.fabG1.include: all` in `configs/data.yaml` to also include fabG1 coding
variants such as `fabG1@L203L` (g609a), which is associated with INH resistance.

Every drug's matrix uses the full panel, not only its own genes.
