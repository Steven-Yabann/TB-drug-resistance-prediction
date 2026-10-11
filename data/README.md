# Data directories

- `raw/`: immutable source data exactly as downloaded
- `external/`: reference annotations and mutation catalogues
- `interim/`: validated, joined, or partially transformed data
- `processed/`: model-ready features, labels, and saved data splits

Large data files are ignored by Git. Preserve provenance, checksums, licences, reference-genome versions, and phenotype definitions in `docs/`.
