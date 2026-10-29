"""End-to-end matrix building on synthetic raw files."""

import json
import shutil

import pandas as pd

from tb_resistance import build, config


def test_build_all(tmp_path, monkeypatch, mutations_path, phenotypes, panel):
    raw = tmp_path / "raw"
    raw.mkdir()
    shutil.copy(mutations_path, raw / "MUTATIONS.csv.gz")
    phenotypes.to_csv(raw / "UKMYC_PHENOTYPES.csv.gz", index=False)
    pd.DataFrame({
        "UNIQUEID": ["iso1", "iso2", "iso3", "iso4", "iso4"],   # iso4 duplicated
        "MYKROBE_LINEAGE_NAME_1": ["Lineage 2", "Lineage 4", "Lineage 4", None, None],
        "MYKROBE_LINEAGE_NAME_2": ["lineage2.2", "lineage4.1", "lineage4.3", None, None],
    }).to_csv(raw / "GENOMES.csv.gz", index=False)

    monkeypatch.setattr(config, "REPO_ROOT", tmp_path)
    cfg = {
        "source": {"files": {"mutations": "x/MUTATIONS.csv.gz",
                             "phenotypes": "x/UKMYC_PHENOTYPES.csv.gz",
                             "genomes": "x/GENOMES.csv.gz"}},
        "paths": {"raw_dir": "raw", "interim_dir": "interim", "processed_dir": "processed"},
        "panel": panel,
        "drugs": {"rifampicin": {"code": "RIF"}, "isoniazid": {"code": "INH"}},
        "labels": {"reading_day": 14},
    }
    report = build.build_all(cfg)

    assert report["duplicate_ids_in_genomes"] == 1
    rif = report["drugs"]["rifampicin"]
    assert rif["labelled_isolates_lost_in_join"] == 1          # iso9
    assert (rif["isolates"], rif["resistant"], rif["susceptible"]) == (3, 1, 2)

    out = tmp_path / "processed" / "rifampicin"
    X = pd.read_parquet(out / "features.parquet")
    y = pd.read_csv(out / "labels.csv", index_col="UNIQUEID")
    lineage = pd.read_csv(out / "lineage.csv", index_col="UNIQUEID")
    assert list(X.index) == list(y.index) == list(lineage.index) == ["iso1", "iso2", "iso4"]
    assert X.dtypes.eq("uint8").all()
    # the full panel is used for every drug
    X_inh = pd.read_parquet(tmp_path / "processed" / "isoniazid" / "features.parquet")
    assert list(X.columns) == list(X_inh.columns)
    assert (tmp_path / "interim" / "panel_mutations.parquet").exists()
    assert json.loads((tmp_path / "processed" / "build_report.json").read_text())

    # second run reuses the interim cache rather than the raw MUTATIONS file
    monkeypatch.setattr(build.features, "filter_panel_mutations",
                        lambda *a, **k: (_ for _ in ()).throw(AssertionError("re-read raw")))
    build.build_all(cfg)
