import pytest

from tb_resistance.features import (annotate, build_matrices, feature_metadata,
                                    filter_panel_mutations)


@pytest.fixture
def annotated(mutations_path, panel):
    return annotate(filter_panel_mutations(mutations_path, panel))


def test_gene_and_promoter_filtering(mutations_path, panel):
    df = filter_panel_mutations(mutations_path, panel)
    ids = set(df["GENE"] + "@" + df["MUTATION"])
    assert "fabG1@c-15t" in ids          # fabG1-inhA promoter kept
    assert "fabG1@L203L" not in ids      # fabG1 coding excluded in promoter mode
    assert "embB@M306V" not in ids       # outside the panel


def test_fabg1_all_mode_keeps_coding(mutations_path, panel):
    df = filter_panel_mutations(mutations_path, {**panel, "fabG1": {"include": "all"}})
    assert "L203L" in set(df.loc[df["GENE"] == "fabG1", "MUTATION"])


def test_unknown_include_mode_rejected(mutations_path):
    with pytest.raises(ValueError):
        filter_panel_mutations(mutations_path, {"rpoB": {"include": "cds"}})


def test_variant_type_and_verbatim_ids(annotated):
    types = annotated.set_index("feature_id")["variant_type"]
    assert types["katG@1000_indel"] == "indel"
    assert types["rpoB@S450L"] == "snp"
    assert types["rpoB@H445Y"] == "mnp"
    assert types["fabG1@c-15t"] == "snp"


def test_features_only_from_passing_non_het_calls(annotated):
    meta = feature_metadata(annotated)
    assert "katG@S315X" not in meta.index   # null
    assert "rpoB@S450Z" not in meta.index   # het
    assert "gyrA@D94G" in meta.index        # passing in iso3, failed in iso2


def test_missing_vs_absent(annotated):
    feats = feature_metadata(annotated).index
    X, missing, het = build_matrices(annotated, ["iso1", "iso2", "iso3", "iso4"], feats)
    assert X.dtypes.eq("uint8").all() and missing.dtypes.eq("uint8").all()

    # iso2: null call at katG codon 315 -> missing, not a confident 0
    assert X.loc["iso2", "katG@S315T"] == 0 and missing.loc["iso2", "katG@S315T"] == 1
    # the katG indel shares no site with codon 315, so it stays observed-absent
    assert missing.loc["iso2", "katG@1000_indel"] == 0
    # filter-failed D94G -> missing for iso2, present for iso3
    assert missing.loc["iso2", "gyrA@D94G"] == 1 and X.loc["iso3", "gyrA@D94G"] == 1
    # het at rpoB 450 -> not present, flagged het, not missing
    assert X.loc["iso2", "rpoB@S450L"] == 0
    assert het.loc["iso2", "rpoB@S450L"] == 1 and missing.loc["iso2", "rpoB@S450L"] == 0
    # isolate without any panel row is all observed-absent
    assert X.loc["iso4"].sum() == 0 and missing.loc["iso4"].sum() == 0


def test_matrix_values(annotated):
    feats = feature_metadata(annotated).index
    X, _, _ = build_matrices(annotated, ["iso1", "iso2", "iso3", "iso4"], feats)
    assert list(X.index) == ["iso1", "iso2", "iso3", "iso4"]
    assert X.index.name == "UNIQUEID"
    assert set(X.columns[X.loc["iso1"] == 1]) == {
        "rpoB@S450L", "katG@S315T", "fabG1@c-15t", "katG@1000_indel"}
    assert set(X.columns[X.loc["iso2"] == 1]) == {"gyrA@A90V"}
