"""Join labels and panel features into per-drug model-ready matrices."""

from __future__ import annotations

import json
import logging
from dataclasses import asdict
from pathlib import Path

import pandas as pd

from tb_resistance import features, labels
from tb_resistance.config import resolve

log = logging.getLogger(__name__)

ID = "UNIQUEID"
INTERIM_NAME = "panel_mutations.parquet"
LINEAGE_COLUMNS = {"MYKROBE_LINEAGE_NAME_1": "lineage", "MYKROBE_LINEAGE_NAME_2": "sublineage"}


def raw_path(cfg: dict, key: str) -> Path:
    return resolve(cfg["paths"]["raw_dir"]) / Path(cfg["source"]["files"][key]).name


def load_panel_mutations(cfg: dict, rebuild: bool = False) -> pd.DataFrame:
    """Annotated panel mutations, cached in the interim dir keyed on the panel config."""
    interim_dir = resolve(cfg["paths"]["interim_dir"])
    cache, stamp = interim_dir / INTERIM_NAME, interim_dir / "panel_mutations.json"
    source = raw_path(cfg, "mutations")
    key = {"panel": cfg["panel"], "source": source.name,
           "source_bytes": source.stat().st_size if source.exists() else None}
    if not rebuild and cache.exists() and stamp.exists() and json.loads(stamp.read_text()) == key:
        log.info("using cached %s", cache)
        return pd.read_parquet(cache)

    log.info("filtering %s to panel genes %s", source.name, list(cfg["panel"]))
    df = features.annotate(features.filter_panel_mutations(source, cfg["panel"]))
    interim_dir.mkdir(parents=True, exist_ok=True)
    df.to_parquet(cache, index=False)
    stamp.write_text(json.dumps(key, indent=2) + "\n")
    return df


def load_genomes(cfg: dict) -> pd.DataFrame:
    cols = [ID, *LINEAGE_COLUMNS]
    return pd.read_csv(raw_path(cfg, "genomes"), usecols=cols, dtype=str)


def build_all(cfg: dict, rebuild_interim: bool = False) -> dict:
    """Write data/processed/<drug>/ for every configured drug and return a report."""
    genomes = load_genomes(cfg)
    dup_genome_ids = genomes.loc[genomes[ID].duplicated(keep=False), ID]
    genomes = genomes.drop_duplicates(ID)
    genotyped = pd.Index(genomes[ID].sort_values(), name=ID)

    muts = load_panel_mutations(cfg, rebuild=rebuild_interim)
    dup_rows = int(muts.duplicated([ID, "feature_id", "call_status"]).sum())
    muts = muts.drop_duplicates([ID, "feature_id", "call_status"])
    orphan_ids = sorted(set(muts[ID]) - set(genotyped))

    meta = features.feature_metadata(muts.loc[muts[ID].isin(genotyped)])
    feature_ids = meta.index

    pheno = labels.read_phenotypes(raw_path(cfg, "phenotypes"))
    reading_day = cfg.get("labels", {}).get("reading_day")
    out_root = resolve(cfg["paths"]["processed_dir"])

    report = {
        "genotyped_isolates": len(genotyped),
        "duplicate_ids_in_genomes": int(dup_genome_ids.nunique()),
        "duplicate_rows_in_panel_mutations": dup_rows,
        "panel_mutation_ids_not_in_genomes": len(orphan_ids),
        "n_features": len(feature_ids),
        "features_by_gene": meta["gene"].value_counts().to_dict(),
        "features_by_type": meta["variant_type"].value_counts().to_dict(),
        "drugs": {},
    }

    for drug, spec in cfg["drugs"].items():
        y, lrep = labels.build_labels(pheno, spec["code"], reading_day=reading_day)
        kept = y.index.intersection(genotyped)
        lost = y.index.difference(genotyped)
        y = y.loc[kept]
        X, missing, het = features.build_matrices(muts, kept, feature_ids)

        out = out_root / drug
        out.mkdir(parents=True, exist_ok=True)
        X.to_parquet(out / "features.parquet")
        missing.to_parquet(out / "features_missing.parquet")
        het.to_parquet(out / "features_het.parquet")
        y.to_csv(out / "labels.csv")
        meta.to_csv(out / "feature_metadata.csv")
        lineage = genomes.set_index(ID).loc[kept, list(LINEAGE_COLUMNS)] \
            .rename(columns=LINEAGE_COLUMNS)
        lineage.to_csv(out / "lineage.csv")

        report["drugs"][drug] = {
            "labels": asdict(lrep),
            "labelled_isolates_lost_in_join": len(lost),
            "isolates": len(kept),
            "resistant": int(y["label"].sum()),
            "susceptible": int(len(y) - y["label"].sum()),
            "n_features": X.shape[1],
            "features_all_zero_in_drug": int((X.sum(axis=0) == 0).sum()),
            "isolates_with_any_missing_flag": int((missing.sum(axis=1) > 0).sum()),
            "isolates_with_any_het_flag": int((het.sum(axis=1) > 0).sum()),
            "isolates_without_lineage": int(lineage["lineage"].isna().sum()),
        }
        log.info("%s: %d isolates (R=%d, S=%d), %d features, %d lost in join", drug,
                 len(kept), report["drugs"][drug]["resistant"],
                 report["drugs"][drug]["susceptible"], X.shape[1], len(lost))

    (out_root / "build_report.json").write_text(json.dumps(report, indent=2, default=str) + "\n")
    return report
