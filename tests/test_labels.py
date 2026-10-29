from tb_resistance.labels import build_labels, read_phenotypes


def test_rs_mapping_and_conflicts(phenotypes):
    y, rep = build_labels(phenotypes, "RIF", reading_day=14)
    assert y["label"].to_dict() == {"iso1": 1, "iso2": 0, "iso4": 0, "iso9": 0}
    assert y["label"].dtype == "uint8"
    assert rep.isolates_conflicting == 1          # iso3
    assert rep.rows_dropped_values == {"<missing>": 1}
    assert (rep.resistant, rep.susceptible) == (1, 3)


def test_non_rs_values_reported_and_quality_kept_unfiltered(phenotypes):
    y, rep = build_labels(phenotypes, "INH", reading_day=14)
    assert rep.rows_dropped_values == {"I": 1}
    assert rep.rows_dropped_reading_day == {"7": 1}
    assert set(y.index) == {"iso1", "iso4"}
    # best quality among agreeing calls; missing quality does not drop the isolate
    assert y.loc["iso1", "phenotype_quality"] == "HIGH"
    assert y.loc["iso1", "n_calls"] == 2
    assert y.loc["iso1", "reading_day"] == "14"
    assert rep.isolates_with_multiple_calls == 1
    assert y.loc["iso4", "label"] == 0


def test_all_reading_days_when_unset(phenotypes):
    y, rep = build_labels(phenotypes, "INH", reading_day=None)
    # iso4 now has R (day 7) and S (day 14): a conflict
    assert "iso4" not in y.index
    assert rep.isolates_conflicting == 1


def test_read_phenotypes_keeps_strings(tmp_path, phenotypes):
    path = tmp_path / "UKMYC_PHENOTYPES.csv.gz"
    phenotypes.assign(EXTRA="x").to_csv(path, index=False)
    df = read_phenotypes(path)
    assert "EXTRA" not in df.columns
    assert df["READINGDAY"].iloc[0] == "14"
