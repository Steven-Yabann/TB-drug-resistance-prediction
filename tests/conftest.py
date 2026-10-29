"""Small synthetic CRyPTIC-shaped fixtures (no network)."""

from __future__ import annotations

import gzip

import pandas as pd
import pytest

from tb_resistance.features import MUTATION_COLUMNS

PANEL = {
    "rpoB": {"include": "all"},
    "katG": {"include": "all"},
    "inhA": {"include": "all"},
    "fabG1": {"include": "promoter"},
    "gyrA": {"include": "all"},
}


def mutation_row(uid, gene, mutation, *, aa=None, nt=None, snp=True, indel=False,
                 promoter=False, het=False, null=False, passed=True, mtype="AAM",
                 changes=1, indel2=None):
    row = dict.fromkeys(MUTATION_COLUMNS, "")
    row.update({
        "UNIQUEID": uid, "GENE": gene, "MUTATION": mutation,
        "AMINO_ACID_NUMBER": "" if aa is None else f"{aa}.0",
        "NUCLEOTIDE_NUMBER": "" if nt is None else f"{nt}.0",
        "POSITION": f"{aa if aa is not None else nt}.0",
        "IS_SNP": str(snp), "IS_INDEL": str(indel), "IN_CDS": str(not promoter),
        "IN_PROMOTER": str(promoter), "IS_HET": str(het), "IS_NULL": str(null),
        "IS_FILTER_PASS": str(passed), "MUTATION_TYPE": mtype,
        "NUMBER_NUCLEOTIDE_CHANGES": str(changes), "INDEL_2": indel2 or "",
    })
    return row


MUTATION_ROWS = [
    # iso1: classic MDR + promoter + an indel
    mutation_row("iso1", "rpoB", "S450L", aa=450),
    mutation_row("iso1", "katG", "S315T", aa=315, changes=1),
    mutation_row("iso1", "fabG1", "c-15t", nt=-15, promoter=True, mtype="SNP", changes=0),
    mutation_row("iso1", "katG", "1000_indel", aa=334, nt=1000, snp=False, indel=True,
                 mtype="INDEL", changes=0, indel2="1000_del_1"),
    # iso2: null at katG codon 315, het at rpoB 450, filter-failed gyrA D94G
    mutation_row("iso2", "katG", "S315X", aa=315, null=True),
    mutation_row("iso2", "rpoB", "S450Z", aa=450, het=True),
    mutation_row("iso2", "gyrA", "D94G", aa=94, passed=False),
    mutation_row("iso2", "gyrA", "A90V", aa=90),
    # iso3: fabG1 coding mutation (excluded), non-panel gene (excluded), gyrA D94G
    mutation_row("iso3", "fabG1", "L203L", aa=203),
    mutation_row("iso3", "embB", "M306V", aa=306),
    mutation_row("iso3", "gyrA", "D94G", aa=94),
    mutation_row("iso3", "rpoB", "H445Y", aa=445, changes=2),
    # iso4 has no panel mutations at all (only in the isolate universe)
]


@pytest.fixture
def panel():
    return PANEL


@pytest.fixture
def mutations_path(tmp_path):
    path = tmp_path / "MUTATIONS.csv.gz"
    with gzip.open(path, "wt") as fh:
        pd.DataFrame(MUTATION_ROWS, columns=MUTATION_COLUMNS).to_csv(fh, index=False)
    return path


@pytest.fixture
def phenotypes():
    rows = [
        # uid, drug, day, call, quality
        ("iso1", "RIF", "14", "R", "HIGH"),
        ("iso1", "INH", "14", "R", "MEDIUM"),
        ("iso1", "INH", "14", "R", "HIGH"),     # agreeing duplicate
        ("iso2", "RIF", "14", "S", "LOW"),
        ("iso2", "INH", "14", "I", "HIGH"),     # intermediate: dropped
        ("iso3", "RIF", "14", "R", "HIGH"),
        ("iso3", "RIF", "14", "S", "HIGH"),     # conflict: dropped
        ("iso4", "RIF", "14", "S", "HIGH"),
        ("iso4", "INH", "7", "R", "HIGH"),      # other reading day
        ("iso4", "INH", "14", "S", None),       # missing quality kept
        ("iso5", "RIF", "14", None, None),      # missing call: dropped
        ("iso9", "RIF", "14", "S", "HIGH"),     # phenotyped but not genotyped
    ]
    return pd.DataFrame(rows, columns=["UNIQUEID", "DRUG", "READINGDAY",
                                       "BINARY_PHENOTYPE", "PHENOTYPE_QUALITY"])
