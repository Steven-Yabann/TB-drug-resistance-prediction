"""Panel-gene mutation features from the CRyPTIC MUTATIONS table.

Pipeline:
1. ``filter_panel_mutations`` streams the raw MUTATIONS table in record batches and
   keeps rows in the panel genes (promoter-only genes keep ``IN_PROMOTER`` rows).
2. ``annotate`` adds the feature ID (``GENE@MUTATION``, source ID verbatim), the
   variant type, a call status and a site key used to match null calls to features.
3. ``build_matrices`` turns the long table into dense uint8 isolate x feature
   matrices: presence (1 = filter-passing, non-het, non-null call), a missing
   flag (null or filter-failed call at the feature's site) and a het flag.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.csv as pacsv

ID = "UNIQUEID"
MUTATION_COLUMNS = [
    ID, "GENE", "MUTATION", "POSITION", "AMINO_ACID_NUMBER", "NUCLEOTIDE_NUMBER",
    "GENOME_INDEX", "REF", "ALT", "IS_SNP", "IS_INDEL", "IN_CDS", "IN_PROMOTER",
    "IS_HET", "IS_NULL", "IS_FILTER_PASS", "MUTATION_TYPE", "INDEL_LENGTH", "INDEL_2",
    "NUMBER_NUCLEOTIDE_CHANGES",
]
BOOL_COLUMNS = ["IS_SNP", "IS_INDEL", "IN_CDS", "IN_PROMOTER", "IS_HET", "IS_NULL",
                "IS_FILTER_PASS"]

PRESENT, MISSING, HET = "present", "missing", "het"


def filter_panel_mutations(path: str | Path, panel: dict[str, dict],
                           block_size: int = 64 << 20) -> pd.DataFrame:
    """Stream ``path`` and return the rows for the panel genes (all values as strings)."""
    genes_all = [g for g, spec in panel.items() if spec.get("include", "all") == "all"]
    genes_promoter = [g for g, spec in panel.items() if spec.get("include") == "promoter"]
    unknown = {spec.get("include", "all") for spec in panel.values()} - {"all", "promoter"}
    if unknown:
        raise ValueError(f"unknown panel include mode(s): {sorted(unknown)}")

    reader = pacsv.open_csv(
        path,
        read_options=pacsv.ReadOptions(block_size=block_size),
        convert_options=pacsv.ConvertOptions(
            include_columns=MUTATION_COLUMNS,
            column_types={c: pa.string() for c in MUTATION_COLUMNS},
            strings_can_be_null=True,
        ),
    )
    kept = []
    for batch in reader:
        gene = batch.column("GENE")
        mask = pc.is_in(gene, value_set=pa.array(genes_all, pa.string()))
        if genes_promoter:
            promoter = pc.and_(
                pc.is_in(gene, value_set=pa.array(genes_promoter, pa.string())),
                pc.equal(batch.column("IN_PROMOTER"), "True"),
            )
            mask = pc.or_(mask, promoter)
        mask = pc.fill_null(mask, False)
        if pc.any(mask).as_py():
            kept.append(batch.filter(mask))
    schema = pa.schema([(c, pa.string()) for c in MUTATION_COLUMNS])
    table = pa.Table.from_batches(kept, schema=schema) if kept else schema.empty_table()
    return table.to_pandas()


def _to_bool(series: pd.Series) -> pd.Series:
    return series.map({"True": True, "False": False, True: True, False: False})


def annotate(df: pd.DataFrame) -> pd.DataFrame:
    """Add feature_id, variant_type, site_key and call_status columns."""
    df = df.copy()
    for col in BOOL_COLUMNS:
        df[col] = _to_bool(df[col]).fillna(False).astype(bool)

    df["feature_id"] = df["GENE"] + "@" + df["MUTATION"]
    n_changes = pd.to_numeric(df["NUMBER_NUCLEOTIDE_CHANGES"], errors="coerce")
    df["variant_type"] = np.select(
        [df["IS_INDEL"], n_changes > 1, df["IS_SNP"]], ["indel", "mnp", "snp"], "other"
    )

    # One site per row: the codon for coding changes, otherwise the nucleotide.
    # Null calls inside a codon are reported per codon (e.g. katG@E195X), so a
    # codon-level key lets them mask every feature in that codon.
    aa = pd.to_numeric(df["AMINO_ACID_NUMBER"], errors="coerce")
    nt = pd.to_numeric(df["NUCLEOTIDE_NUMBER"], errors="coerce")
    df["site_key"] = np.where(
        aa.notna(),
        df["GENE"] + ":aa:" + aa.astype("Int64").astype(str),
        df["GENE"] + ":nt:" + nt.astype("Int64").astype(str),
    )

    failed = df["IS_NULL"] | ~df["IS_FILTER_PASS"]
    df["call_status"] = np.select([failed, df["IS_HET"]], [MISSING, HET], PRESENT)
    return df


def feature_metadata(annotated: pd.DataFrame) -> pd.DataFrame:
    """One row per feature (present calls only): gene, mutation, type and location."""
    present = annotated.loc[annotated["call_status"] == PRESENT]
    meta = (present.groupby("feature_id")
            .agg(gene=("GENE", "first"), mutation=("MUTATION", "first"),
                 variant_type=("variant_type", "first"), site_key=("site_key", "first"),
                 in_promoter=("IN_PROMOTER", "first"), indel_detail=("INDEL_2", "first"))
            .sort_index())
    return meta


def _dense(pairs: pd.DataFrame, index: pd.Index, columns: pd.Index) -> pd.DataFrame:
    """Dense uint8 matrix with 1 at each (isolate, feature) pair in ``pairs``."""
    out = np.zeros((len(index), len(columns)), dtype=np.uint8)
    rows = index.get_indexer(pairs[ID])
    cols = columns.get_indexer(pairs["feature_id"])
    ok = (rows >= 0) & (cols >= 0)
    out[rows[ok], cols[ok]] = 1
    return pd.DataFrame(out, index=index, columns=columns)


def build_matrices(annotated: pd.DataFrame, isolate_ids, feature_ids
                   ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return (features, missing, het) uint8 matrices indexed by isolate ID.

    An isolate absent from ``annotated`` has no panel mutations: all zeros.
    ``missing`` is 1 where the isolate has a null or filter-failed call at the
    feature's site; ``features`` is 0 there and should be treated as unknown.
    ``het`` is 1 where a heterozygous (minor-allele) call was made at the site.
    """
    index = pd.Index(isolate_ids, name=ID)
    columns = pd.Index(feature_ids, name="feature_id")
    site_of = annotated.loc[annotated["call_status"] == PRESENT] \
        .drop_duplicates("feature_id").set_index("feature_id")["site_key"]
    site_of = site_of.reindex(columns)
    site_features = pd.DataFrame({"feature_id": columns, "site_key": site_of.values})

    present = annotated.loc[annotated["call_status"] == PRESENT, [ID, "feature_id"]]

    def flagged(status: str) -> pd.DataFrame:
        sites = annotated.loc[annotated["call_status"] == status, [ID, "site_key"]]
        return sites.drop_duplicates().merge(site_features, on="site_key")[[ID, "feature_id"]]

    features = _dense(present, index, columns)
    missing = _dense(flagged(MISSING), index, columns)
    het = _dense(flagged(HET), index, columns)
    # A passing call for the exact feature overrides a missing/het flag at the same site.
    missing = missing.where(features == 0, 0).astype(np.uint8)
    het = het.where(features == 0, 0).astype(np.uint8)
    return features, missing, het
