"""Binary resistance labels from CRyPTIC UKMYC_PHENOTYPES binary phenotype calls.

R -> 1, S -> 0. Every other BINARY_PHENOTYPE value (e.g. I, missing) is dropped and
counted. Isolates with both R and S calls for the same drug are dropped as
conflicting. PHENOTYPE_QUALITY is carried through but never used as a filter.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

ID = "UNIQUEID"
LABEL_MAP = {"R": 1, "S": 0}
QUALITY_RANK = {"HIGH": 3, "MEDIUM": 2, "LOW": 1}
PHENOTYPE_COLUMNS = [ID, "DRUG", "READINGDAY", "BINARY_PHENOTYPE", "PHENOTYPE_QUALITY"]


@dataclass
class LabelReport:
    drug: str
    rows_for_drug: int = 0
    rows_dropped_reading_day: dict[str, int] = field(default_factory=dict)
    rows_dropped_values: dict[str, int] = field(default_factory=dict)
    isolates_with_multiple_calls: int = 0
    isolates_conflicting: int = 0
    isolates_kept: int = 0
    resistant: int = 0
    susceptible: int = 0


def read_phenotypes(path: str | Path) -> pd.DataFrame:
    """Read the columns needed for labels; values are kept as strings."""
    return pd.read_csv(path, usecols=PHENOTYPE_COLUMNS, dtype=str, keep_default_na=False,
                       na_values=[""])


def _value_counts(series: pd.Series) -> dict[str, int]:
    counts = series.fillna("<missing>").value_counts()
    return {str(k): int(v) for k, v in counts.items()}


def build_labels(pheno: pd.DataFrame, drug_code: str,
                 reading_day: int | None = None) -> tuple[pd.DataFrame, LabelReport]:
    """Return one row per isolate (indexed by UNIQUEID) with ``label`` in {0, 1}.

    Output columns: label (uint8), phenotype_quality (best quality among the
    isolate's kept calls), reading_day (of that call), n_calls (number of R/S
    calls that agreed).
    """
    report = LabelReport(drug=drug_code)
    df = pheno.loc[pheno["DRUG"] == drug_code]
    report.rows_for_drug = len(df)

    if reading_day is not None:
        off_day = df["READINGDAY"].astype(float) != float(reading_day)
        report.rows_dropped_reading_day = _value_counts(df.loc[off_day, "READINGDAY"])
        df = df.loc[~off_day]

    keep = df["BINARY_PHENOTYPE"].isin(list(LABEL_MAP))
    report.rows_dropped_values = _value_counts(df.loc[~keep, "BINARY_PHENOTYPE"])
    df = df.loc[keep].assign(label=lambda d: d["BINARY_PHENOTYPE"].map(LABEL_MAP))

    grouped = df.groupby(ID)
    n_labels = grouped["label"].nunique()
    n_calls = grouped.size()
    report.isolates_with_multiple_calls = int((n_calls > 1).sum())
    conflicting = n_labels.index[n_labels > 1]
    report.isolates_conflicting = len(conflicting)
    df = df.loc[~df[ID].isin(conflicting)]

    rank = df["PHENOTYPE_QUALITY"].map(QUALITY_RANK).fillna(0)
    best = df.assign(_rank=rank).sort_values([ID, "_rank"], ascending=[True, False])
    labels = best.drop_duplicates(ID).set_index(ID)[["label", "PHENOTYPE_QUALITY", "READINGDAY"]]
    labels = labels.rename(columns={"PHENOTYPE_QUALITY": "phenotype_quality",
                                    "READINGDAY": "reading_day"})
    labels["label"] = labels["label"].astype("uint8")
    labels["n_calls"] = n_calls.reindex(labels.index).astype(int)
    labels.index.name = ID
    labels = labels.sort_index()

    report.isolates_kept = len(labels)
    report.resistant = int(labels["label"].sum())
    report.susceptible = report.isolates_kept - report.resistant
    return labels, report
