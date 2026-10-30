"""Exploratory data analysis helpers for TB drug resistance feature matrices and labels."""

from __future__ import annotations

from pathlib import Path
from typing import NamedTuple

import numpy as np
import pandas as pd


class DrugDataset(NamedTuple):
    drug: str
    features: pd.DataFrame
    labels: pd.DataFrame
    lineage: pd.DataFrame
    metadata: pd.DataFrame
    missing: pd.DataFrame | None = None
    het: pd.DataFrame | None = None


def load_drug_dataset(
    processed_dir: str | Path,
    drug: str,
    load_quality_matrices: bool = False,
) -> DrugDataset:
    """Load the processed datasets for a single drug.

    Parameters
    ----------
    processed_dir : str or Path
        Directory containing per-drug subdirectories (e.g. data/processed).
    drug : str
        Drug name (e.g. 'rifampicin', 'isoniazid', 'levofloxacin', 'moxifloxacin').
    load_quality_matrices : bool, default False
        Whether to also load features_missing.parquet and features_het.parquet.

    Returns
    -------
    DrugDataset
        Named tuple containing features, labels, lineage, metadata, and optional quality matrices.
    """
    drug_dir = Path(processed_dir) / drug
    if not drug_dir.exists():
        raise FileNotFoundError(f"Directory not found for drug '{drug}': {drug_dir}")

    features = pd.read_parquet(drug_dir / "features.parquet")
    labels = pd.read_csv(drug_dir / "labels.csv", index_col="UNIQUEID")
    lineage = pd.read_csv(drug_dir / "lineage.csv", index_col="UNIQUEID")
    metadata = pd.read_csv(drug_dir / "feature_metadata.csv", index_col="feature_id")

    missing = None
    het = None
    if load_quality_matrices:
        missing_path = drug_dir / "features_missing.parquet"
        if missing_path.exists():
            missing = pd.read_parquet(missing_path)
        het_path = drug_dir / "features_het.parquet"
        if het_path.exists():
            het = pd.read_parquet(het_path)

    return DrugDataset(
        drug=drug,
        features=features,
        labels=labels,
        lineage=lineage,
        metadata=metadata,
        missing=missing,
        het=het,
    )


def summarize_drug_cohorts(
    processed_dir: str | Path,
    drugs: list[str] | None = None,
) -> pd.DataFrame:
    """Generate a summary table of cohort size, label balance, and quality across drugs."""
    root = Path(processed_dir)
    if drugs is None:
        drugs = [d.name for d in root.iterdir() if d.is_dir() and (d / "labels.csv").exists()]

    rows = []
    for drug in sorted(drugs):
        ds = load_drug_dataset(root, drug)
        y = ds.labels["label"]
        n_total = len(y)
        n_res = int((y == 1).sum())
        n_sus = int((y == 0).sum())
        res_pct = (n_res / n_total * 100.0) if n_total > 0 else 0.0

        # Phenotype quality counts
        qual_counts = ds.labels["phenotype_quality"].value_counts().to_dict()
        high_qual = qual_counts.get("HIGH", 0)

        # Reading days
        reading_days = ds.labels["reading_day"].value_counts().to_dict()

        rows.append({
            "drug": drug,
            "total_isolates": n_total,
            "resistant": n_res,
            "susceptible": n_sus,
            "prevalence_pct": round(res_pct, 2),
            "high_quality_n": high_qual,
            "high_quality_pct": round(high_qual / n_total * 100.0, 1) if n_total > 0 else 0.0,
            "day_14_n": reading_days.get("14", reading_days.get(14, 0)),
            "n_features": ds.features.shape[1],
        })

    return pd.DataFrame(rows).set_index("drug")


def top_mutations_by_resistance(
    features: pd.DataFrame,
    labels: pd.Series,
    metadata: pd.DataFrame,
    gene: str | None = None,
    top_n: int = 20,
    min_resistant_count: int = 5,
) -> pd.DataFrame:
    """Identify mutations most strongly enriched in resistant vs. susceptible isolates.

    Parameters
    ----------
    features : pd.DataFrame
        Binary feature matrix (isolates x mutations).
    labels : pd.Series
        Binary resistance labels (1 = resistant, 0 = susceptible).
    metadata : pd.DataFrame
        Metadata table indexed by feature_id with columns: gene, mutation, variant_type.
    gene : str, optional
        Filter results to mutations in this gene (e.g. 'rpoB').
    top_n : int, default 20
        Number of top mutations to return.
    min_resistant_count : int, default 5
        Minimum count of resistant isolates carrying the mutation.

    Returns
    -------
    pd.DataFrame
        Ranked table of mutations with counts, frequencies, risk ratio, and delta %.
    """
    if gene is not None:
        gene_feats = metadata[metadata["gene"] == gene].index.intersection(features.columns)
        features = features[gene_feats]

    is_res = (labels == 1)
    is_sus = (labels == 0)
    n_res = is_res.sum()
    n_sus = is_sus.sum()

    res_counts = features.loc[is_res].sum(axis=0)
    sus_counts = features.loc[is_sus].sum(axis=0)

    # Filter to mutations meeting min_resistant_count
    mask = res_counts >= min_resistant_count
    filtered_feats = res_counts[mask].index

    res_freq = res_counts[filtered_feats] / n_res
    sus_freq = sus_counts[filtered_feats] / n_sus
    delta = res_freq - sus_freq

    # Odds ratio with Haldane-Anscombe correction (+0.5) to avoid division by zero
    a = res_counts[filtered_feats] + 0.5
    b = (n_res - res_counts[filtered_feats]) + 0.5
    c = sus_counts[filtered_feats] + 0.5
    d = (n_sus - sus_counts[filtered_feats]) + 0.5
    odds_ratio = (a * d) / (b * c)

    df = pd.DataFrame({
        "resistant_count": res_counts[filtered_feats].astype(int),
        "resistant_freq_pct": (res_freq * 100.0).round(2),
        "susceptible_count": sus_counts[filtered_feats].astype(int),
        "susceptible_freq_pct": (sus_freq * 100.0).round(2),
        "delta_pct": (delta * 100.0).round(2),
        "odds_ratio": odds_ratio.round(2),
    }, index=filtered_feats)

    # Join metadata
    df = df.join(metadata[["gene", "mutation", "variant_type", "in_promoter"]], how="left")
    df = df.sort_values(by="delta_pct", ascending=False).head(top_n)
    return df


def lineage_distribution(
    labels: pd.DataFrame,
    lineage: pd.DataFrame,
) -> pd.DataFrame:
    """Compute isolate counts and resistance prevalence stratified by M. tuberculosis lineage."""
    merged = labels[["label"]].join(lineage[["lineage", "sublineage"]], how="inner")
    
    # Fill uncalled lineages
    merged["lineage"] = merged["lineage"].fillna("Unclassified")
    
    grouped = merged.groupby("lineage").agg(
        total_isolates=("label", "count"),
        resistant=("label", lambda s: (s == 1).sum()),
        susceptible=("label", lambda s: (s == 0).sum()),
    )
    grouped["prevalence_pct"] = (grouped["resistant"] / grouped["total_isolates"] * 100.0).round(2)
    grouped["share_of_cohort_pct"] = (grouped["total_isolates"] / len(merged) * 100.0).round(2)
    return grouped.sort_values(by="total_isolates", ascending=False)


def gene_burden_analysis(
    features: pd.DataFrame,
    labels: pd.Series,
    metadata: pd.DataFrame,
) -> pd.DataFrame:
    """Analyze the proportion of isolates carrying >=1 mutation in each panel gene."""
    genes = metadata["gene"].unique()
    is_res = (labels == 1)
    is_sus = (labels == 0)
    n_res = is_res.sum()
    n_sus = is_sus.sum()

    rows = []
    for g in sorted(genes):
        gene_feats = metadata[metadata["gene"] == g].index.intersection(features.columns)
        if len(gene_feats) == 0:
            continue
        # Any mutation in gene
        any_mut = (features[gene_feats] > 0).any(axis=1)
        res_with_mut = any_mut[is_res].sum()
        sus_with_mut = any_mut[is_sus].sum()

        rows.append({
            "gene": g,
            "n_features": len(gene_feats),
            "res_any_mut_n": int(res_with_mut),
            "res_any_mut_pct": round(res_with_mut / n_res * 100.0, 2),
            "sus_any_mut_n": int(sus_with_mut),
            "sus_any_mut_pct": round(sus_with_mut / n_sus * 100.0, 2),
            "delta_pct": round((res_with_mut / n_res - sus_with_mut / n_sus) * 100.0, 2),
        })

    return pd.DataFrame(rows).set_index("gene")


def multidrug_overlap(
    processed_dir: str | Path,
    drugs: list[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute shared isolates and multi-drug resistance cross-tabulations.

    Returns
    -------
    tuple of (shared_counts_matrix, resistance_concordance_matrix)
    """
    root = Path(processed_dir)
    if drugs is None:
        drugs = ["rifampicin", "isoniazid", "levofloxacin", "moxifloxacin"]

    label_dfs = {}
    for d in drugs:
        ds = load_drug_dataset(root, d)
        label_dfs[d] = ds.labels["label"]

    combined = pd.DataFrame(label_dfs)

    # Pairwise overlap matrix
    shared = pd.DataFrame(index=drugs, columns=drugs, dtype=int)
    both_res = pd.DataFrame(index=drugs, columns=drugs, dtype=int)

    for d1 in drugs:
        for d2 in drugs:
            s1 = combined[d1].dropna()
            s2 = combined[d2].dropna()
            common_idx = s1.index.intersection(s2.index)
            shared.loc[d1, d2] = len(common_idx)
            both_res.loc[d1, d2] = int(((s1.loc[common_idx] == 1) & (s2.loc[common_idx] == 1)).sum())

    return shared, both_res
