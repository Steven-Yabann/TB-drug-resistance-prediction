# Genomic Prediction of TB Drug Resistance

This project investigates whether single-nucleotide polymorphisms (SNPs) in known resistance-associated genes can predict phenotypic drug resistance in *Mycobacterium tuberculosis*.

The initial scope focuses on a small, biologically interpretable gene panel:

- `rpoB` for rifampicin resistance
- `katG` and `inhA` for isoniazid resistance
- `gyrA` for fluoroquinolone resistance
- `pncA` for pyrazinamide resistance (optional extension)

The intended outcome is a reproducible pipeline that links genomic variants to laboratory-measured resistance labels, trains one binary classifier per drug, evaluates predictive performance, and identifies the SNPs that contribute most strongly to each prediction.

## Project status

- **Done:** CRyPTIC data loader (`release_june2022`): download with checksums, binary labels for rifampicin, isoniazid, levofloxacin and moxifloxacin, and isolate-by-mutation matrices for the rpoB / katG / inhA / fabG1-inhA promoter / gyrA panel. See [`docs/data_sources/cryptic.md`](docs/data_sources/cryptic.md) and [`docs/data_dictionary.md`](docs/data_dictionary.md).
- **Not possible from this release:** pyrazinamide. UKMYC plates do not test PZA, so the release has no PZA phenotypes and the `pncA` extension needs another data source.
- **Next:** train/test splitting (lineage-aware), training-set-only feature filtering, baseline models and evaluation.

## Research question

Can SNPs in established resistance-associated genes accurately predict phenotypic drug resistance in *M. tuberculosis* isolates, and which mutations contribute most to prediction accuracy across different drug classes?

## Objectives

1. Curate and preprocess a public *M. tuberculosis* whole-genome sequencing dataset, extracting SNPs in the target gene panel and matching them to phenotypic resistance labels.
2. Build and evaluate classification models, initially logistic regression and random forest, that predict resistance separately for each drug.
3. Rank predictive SNPs and compare model performance across drug classes.
4. Optionally compare model predictions with a rule-based lookup using a recognized mutation catalogue, such as the WHO catalogue of mutations associated with TB drug resistance.

## Planned analysis flow

```text
Public WGS variants + phenotypic labels
                  |
                  v
       Sample and label validation
                  |
                  v
       Target-gene SNP extraction
                  |
                  v
   Isolate-by-SNP presence/absence matrix
                  |
                  v
     Train/validation/test separation
                  |
                  v
      Per-drug classification models
                  |
                  v
 Performance metrics + important mutations
```

Each row of the modelling matrix will represent one isolate. Each SNP column will indicate whether a mutation is present (`1`) or absent (`0`), and the target column will describe whether the isolate is phenotypically resistant (`1`) or susceptible (`0`) to a particular drug.

## Repository structure

```text
.
├── configs/                 # Version-controlled pipeline and model settings
├── data/
│   ├── raw/                 # Original, immutable downloaded data (not committed)
│   ├── external/            # Reference catalogues and gene annotations
│   ├── interim/             # Intermediate validated or filtered data
│   └── processed/           # Model-ready feature matrices and labels
├── docs/                    # Data dictionary, decisions, and project documentation
├── models/                  # Serialized fitted models (not committed)
├── notebooks/               # Numbered exploratory and analysis notebooks
├── reports/
│   ├── figures/             # Final plots
│   └── tables/              # Final result tables
├── scripts/                 # Command-line pipeline entry points
├── src/
│   └── tb_resistance/       # Reusable Python package
├── tests/                   # Automated tests
├── .gitignore
└── README.md
```

### Directory conventions

- Treat `data/raw/` as immutable: scripts may read raw data but should never overwrite it.
- Put reusable logic in `src/tb_resistance/`; notebooks should call that logic instead of duplicating it.
- Use numbered notebook names such as `01_data_exploration.ipynb` and `02_feature_analysis.ipynb` so the intended order is clear.
- Keep generated models and large datasets out of Git. Small final tables and figures may be committed when they support the report.

## Candidate data sources

The dataset must contain both genomic variants and matched phenotypic drug-susceptibility results. Candidate sources include:

- CRyPTIC Consortium public releases
- NCBI or ENA collections used in published TB resistance studies
- BV-BRC curated *M. tuberculosis* antimicrobial-resistance data

Before choosing a source, verify its licence, sample identifiers, phenotype definitions, genome reference version, missing-label rate, and representation of resistant and susceptible isolates for every drug.

No dataset is bundled with this repository.

## Getting started

Requires Python 3.10 or newer.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
python -m pytest
```

Build the CRyPTIC matrices (downloads about 1.5 GB into `data/raw/cryptic/` on the first run; interrupted downloads resume):

```bash
python scripts/build_matrices.py                    # download + build
python scripts/build_matrices.py --skip-download    # rebuild from existing raw files
python scripts/build_matrices.py --rebuild-interim  # re-filter raw MUTATIONS after changing the panel
```

Settings (source files, gene panel, drugs, paths) are in `configs/data.yaml`. Outputs go to `data/processed/<drug>/`. Every column is described in `docs/data_dictionary.md`.

Place any other downloaded source files under `data/raw/` and record their origin, version, download date, licence, checksums, and phenotype definitions in `docs/`. Do not commit controlled-access data or sensitive sample metadata.

## Modelling plan

Models will be trained separately for each drug because each drug has a different resistance mechanism and may have different sample availability.

The initial baselines are:

- **Logistic regression:** a transparent model whose coefficients can show how individual SNPs affect predicted resistance.
- **Random forest:** a nonlinear model that can capture interactions among mutations, although its importance scores require careful interpretation.

Data splitting should account for related bacterial isolates. A purely random split can put genetically similar isolates in both training and test sets, producing an overly optimistic estimate of performance. Where metadata permits, the project should evaluate lineage-aware or cluster-aware splitting.

## Evaluation

Accuracy alone can be misleading when resistant isolates are uncommon. For each drug, report:

- Sensitivity (recall for resistant isolates)
- Specificity (recall for susceptible isolates)
- Precision
- F1 score
- ROC-AUC
- Precision-recall AUC
- Confusion matrix
- Calibration, if predicted probabilities will be interpreted clinically

Threshold selection must be documented. Clinical interpretation should emphasize false-susceptible errors because failing to detect resistance can lead to ineffective treatment.

## Reproducibility and leakage controls

- Fit preprocessing and feature-selection steps using training data only.
- Keep isolate identifiers consistent across genomic and phenotype tables.
- Record the reference genome and variant-normalization procedure.
- Define how missing, heteroresistant, low-coverage, and conflicting phenotypes are handled.
- Set and record random seeds.
- Save train/test sample lists and configuration files with each experiment.
- Check for duplicate or near-duplicate isolates before splitting.

## Important limitations

- A targeted gene panel may miss novel resistance mechanisms and mutations outside the selected genes.
- Phenotypic drug-susceptibility testing contains measurement uncertainty and may differ among laboratories.
- Mutation frequencies and lineage structure vary geographically, so strong performance on a global dataset does not guarantee equivalent performance in East Africa.
- Feature importance indicates predictive association, not necessarily a causal biological mechanism.
- This is a research project, not a clinically validated diagnostic tool.

## Expected outputs

- A cleaned isolate-by-SNP feature matrix with per-drug phenotype labels
- Reproducible training and evaluation scripts
- Per-drug model performance tables and plots
- Ranked resistance-associated mutations
- A discussion of biological interpretation, generalizability, and relevance to TB programs in Kenya and East Africa

## Contributing workflow

1. Create a focused branch for the change.
2. Keep raw and generated large files out of Git.
3. Add tests for reusable preprocessing and modelling functions.
4. Run the test suite before opening a pull request.
5. Document any new dataset, feature definition, or modelling assumption.

## Licence

A project licence has not yet been selected. Dataset licences and citation requirements remain applicable independently of the eventual code licence.
