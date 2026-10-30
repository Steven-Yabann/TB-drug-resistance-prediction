import numpy as np
import pandas as pd
import pytest

from tb_resistance.eda import (
    DrugDataset,
    gene_burden_analysis,
    lineage_distribution,
    load_drug_dataset,
    multidrug_overlap,
    summarize_drug_cohorts,
    top_mutations_by_resistance,
)


@pytest.fixture
def dummy_processed_data(tmp_path):
    # Setup dummy processed structure for rifampicin and isoniazid
    drugs = ["rifampicin", "isoniazid"]
    isolates = ["iso1", "iso2", "iso3", "iso4"]
    features_cols = ["rpoB@S450L", "katG@S315T", "inhA@c-15t", "gyrA@D94G"]

    for d in drugs:
        drug_dir = tmp_path / d
        drug_dir.mkdir(parents=True)

        # features.parquet
        feat_df = pd.DataFrame(
            [
                [1, 0, 0, 0],
                [1, 1, 0, 0],
                [0, 0, 0, 0],
                [0, 0, 0, 1],
            ],
            index=pd.Index(isolates, name="UNIQUEID"),
            columns=features_cols,
            dtype="uint8",
        )
        feat_df.to_parquet(drug_dir / "features.parquet")

        # labels.csv
        labels_df = pd.DataFrame(
            {
                "UNIQUEID": isolates,
                "label": [1, 1, 0, 0] if d == "rifampicin" else [0, 1, 0, 0],
                "phenotype_quality": ["HIGH", "HIGH", "MEDIUM", "HIGH"],
                "reading_day": ["14", "14", "14", "21"],
                "n_calls": [1, 1, 1, 1],
            }
        ).set_index("UNIQUEID")
        labels_df.to_csv(drug_dir / "labels.csv")

        # lineage.csv
        lineage_df = pd.DataFrame(
            {
                "UNIQUEID": isolates,
                "lineage": ["Lineage 4", "Lineage 2", "Lineage 4", "Lineage 1"],
                "sublineage": ["lineage4.1", "lineage2.2", "lineage4.2", "lineage1.1"],
            }
        ).set_index("UNIQUEID")
        lineage_df.to_csv(drug_dir / "lineage.csv")

        # feature_metadata.csv
        meta_df = pd.DataFrame(
            {
                "feature_id": features_cols,
                "gene": ["rpoB", "katG", "inhA", "gyrA"],
                "mutation": ["S450L", "S315T", "c-15t", "D94G"],
                "variant_type": ["snp", "snp", "snp", "snp"],
                "site_key": ["rpoB:aa:450", "katG:aa:315", "inhA:nt:-15", "gyrA:aa:94"],
                "in_promoter": [False, False, True, False],
                "indel_detail": [None, None, None, None],
            }
        ).set_index("feature_id")
        meta_df.to_csv(drug_dir / "feature_metadata.csv")

    return tmp_path


def test_load_drug_dataset(dummy_processed_data):
    ds = load_drug_dataset(dummy_processed_data, "rifampicin")
    assert isinstance(ds, DrugDataset)
    assert ds.features.shape == (4, 4)
    assert len(ds.labels) == 4
    assert len(ds.lineage) == 4
    assert len(ds.metadata) == 4


def test_summarize_drug_cohorts(dummy_processed_data):
    summary = summarize_drug_cohorts(dummy_processed_data)
    assert "rifampicin" in summary.index
    assert "isoniazid" in summary.index
    assert summary.loc["rifampicin", "total_isolates"] == 4
    assert summary.loc["rifampicin", "resistant"] == 2
    assert summary.loc["rifampicin", "prevalence_pct"] == 50.0
    assert summary.loc["rifampicin", "high_quality_pct"] == 75.0


def test_top_mutations_by_resistance(dummy_processed_data):
    ds = load_drug_dataset(dummy_processed_data, "rifampicin")
    top = top_mutations_by_resistance(
        ds.features,
        ds.labels["label"],
        ds.metadata,
        min_resistant_count=1,
    )
    assert "rpoB@S450L" in top.index
    assert top.loc["rpoB@S450L", "resistant_count"] == 2
    assert top.loc["rpoB@S450L", "susceptible_count"] == 0
    assert top.loc["rpoB@S450L", "delta_pct"] == 100.0


def test_lineage_distribution(dummy_processed_data):
    ds = load_drug_dataset(dummy_processed_data, "rifampicin")
    dist = lineage_distribution(ds.labels, ds.lineage)
    assert "Lineage 4" in dist.index
    assert dist.loc["Lineage 4", "total_isolates"] == 2
    assert dist.loc["Lineage 4", "resistant"] == 1


def test_gene_burden_analysis(dummy_processed_data):
    ds = load_drug_dataset(dummy_processed_data, "rifampicin")
    burden = gene_burden_analysis(ds.features, ds.labels["label"], ds.metadata)
    assert "rpoB" in burden.index
    assert burden.loc["rpoB", "res_any_mut_pct"] == 100.0
    assert burden.loc["rpoB", "sus_any_mut_pct"] == 0.0


def test_multidrug_overlap(dummy_processed_data):
    shared, both_res = multidrug_overlap(dummy_processed_data, ["rifampicin", "isoniazid"])
    assert shared.loc["rifampicin", "isoniazid"] == 4
    assert both_res.loc["rifampicin", "isoniazid"] == 1

