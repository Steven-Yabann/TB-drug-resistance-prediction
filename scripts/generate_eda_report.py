"""Generate EDA figures and summary tables from processed CRyPTIC matrices."""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from tb_resistance.eda import (
    gene_burden_analysis,
    lineage_distribution,
    load_drug_dataset,
    multidrug_overlap,
    summarize_drug_cohorts,
    top_mutations_by_resistance,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)


def set_plotting_theme():
    sns.set_theme(style="whitegrid", font="sans-serif")
    plt.rcParams.update({
        "figure.titlesize": 14,
        "axes.titlesize": 12,
        "axes.labelsize": 11,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
        "figure.autolayout": True,
    })


def main():
    processed_dir = Path("data/processed")
    figures_dir = Path("reports/figures")
    tables_dir = Path("reports/tables")
    figures_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)

    set_plotting_theme()
    drugs = ["rifampicin", "isoniazid", "levofloxacin", "moxifloxacin"]

    # 1. Cohort Summary
    log.info("Generating cohort summary...")
    cohort_summary = summarize_drug_cohorts(processed_dir, drugs)
    cohort_summary.to_csv(tables_dir / "cohort_summary.csv")
    # summarize_drug_cohorts returns drugs alphabetically; bars are drawn by position
    cohort_summary = cohort_summary.loc[drugs]

    fig, ax1 = plt.subplots(figsize=(8, 4.5), dpi=300)
    colors = ["#2b5c8f", "#d95f02", "#7570b3", "#1b9e77"]
    x = np.arange(len(drugs))
    width = 0.35

    ax1.bar(x - width/2, cohort_summary["resistant"], width, label="Resistant", color="#c93b2b", alpha=0.9)
    ax1.bar(x + width/2, cohort_summary["susceptible"], width, label="Susceptible", color="#4a7c59", alpha=0.9)
    ax1.set_ylabel("Number of Isolates")
    ax1.set_title("CRyPTIC Dataset: Cohort Size and Label Balance per Drug")
    ax1.set_xticks(x)
    ax1.set_xticklabels([d.capitalize() for d in drugs])
    ax1.legend(loc="upper left")

    # Annotate prevalence % just above the taller bar of each pair
    tallest = cohort_summary[["resistant", "susceptible"]].max(axis=1)
    for i, d in enumerate(drugs):
        prev = cohort_summary.loc[d, "prevalence_pct"]
        total = cohort_summary.loc[d, "total_isolates"]
        ax1.text(i, tallest[d] * 1.03, f"{prev:.1f}% Res\n(N={total:,})", ha="center",
                 va="bottom", fontsize=9,
                 bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.8, edgecolor="#ccc"))
    ax1.set_ylim(0, tallest.max() * 1.25)

    plt.savefig(figures_dir / "cohort_prevalence.png")
    plt.close()

    # 2. Multidrug Co-Resistance Heatmap
    log.info("Generating multidrug resistance heatmap...")
    shared, both_res = multidrug_overlap(processed_dir, drugs)
    shared.to_csv(tables_dir / "multidrug_shared_isolates.csv")
    both_res.to_csv(tables_dir / "multidrug_co_resistance.csv")

    # Joint resistance rate (% of isolates tested for both drugs that are resistant to both)
    concordance = (both_res / shared * 100.0).astype(float)

    fig, ax = plt.subplots(figsize=(7, 5.5), dpi=300)
    sns.heatmap(
        concordance,
        annot=True,
        fmt=".1f",
        cmap="YlOrRd",
        cbar_kws={"label": "% of Shared Isolates Resistant to Both"},
        ax=ax,
        xticklabels=[d.capitalize() for d in drugs],
        yticklabels=[d.capitalize() for d in drugs],
    )
    ax.set_title("Pairwise Co-Resistance: % of Isolates Resistant to Both Drugs\n"
                 "(diagonal = resistance prevalence for that drug)", fontsize=11)
    plt.savefig(figures_dir / "co_resistance_heatmap.png")
    plt.close()

    # 3. Lineage Distribution (Rifampicin & Isoniazid)
    log.info("Generating lineage distribution...")
    ds_rif = load_drug_dataset(processed_dir, "rifampicin")
    lin_rif = lineage_distribution(ds_rif.labels, ds_rif.lineage)
    lin_rif.to_csv(tables_dir / "lineage_summary_rifampicin.csv")

    # Focus on major lineages (N >= 50)
    major_lin = lin_rif[lin_rif["total_isolates"] >= 50].copy()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.8), dpi=300)
    # Lineage share
    sns.barplot(data=major_lin.reset_index(), x="lineage", y="share_of_cohort_pct", ax=ax1, palette="crest")
    ax1.set_title("Major Lineage Representation (% of Isolates)")
    ax1.set_xlabel("")
    ax1.set_ylabel("% of Cohort")
    ax1.tick_params(axis="x", rotation=25)

    # Lineage resistance prevalence
    sns.barplot(data=major_lin.reset_index(), x="lineage", y="prevalence_pct", ax=ax2, palette="flare")
    ax2.set_title("Rifampicin Resistance Prevalence by Lineage (%)")
    ax2.set_xlabel("")
    ax2.set_ylabel("Resistance Prevalence (%)")
    ax2.tick_params(axis="x", rotation=25)

    plt.savefig(figures_dir / "lineage_prevalence.png")
    plt.close()

    # 4. Top Mutations for Each Drug
    log.info("Extracting and plotting top mutations...")
    gene_map = {
        "rifampicin": "rpoB",
        "isoniazid": "katG",
        "levofloxacin": "gyrA",
        "moxifloxacin": "gyrA",
    }

    fig, axes = plt.subplots(2, 2, figsize=(14, 10), dpi=300)
    axes = axes.flatten()

    for idx, drug in enumerate(drugs):
        ds = load_drug_dataset(processed_dir, drug)
        target_gene = gene_map[drug]
        top = top_mutations_by_resistance(
            ds.features,
            ds.labels["label"],
            ds.metadata,
            gene=target_gene,
            top_n=10,
            min_resistant_count=3,
        )
        top.to_csv(tables_dir / f"top_mutations_{drug}.csv")

        # Plot top 8
        plot_df = top.head(8).reset_index()
        ax = axes[idx]
        y_pos = np.arange(len(plot_df))
        ax.barh(y_pos - 0.18, plot_df["resistant_freq_pct"], height=0.36, color="#c93b2b", label="Resistant %")
        ax.barh(y_pos + 0.18, plot_df["susceptible_freq_pct"], height=0.36, color="#4a7c59", label="Susceptible %")
        ax.set_yticks(y_pos)
        ax.set_yticklabels(plot_df["mutation"])
        ax.invert_yaxis()
        ax.set_xlabel("Frequency (%)")
        ax.set_title(f"{drug.capitalize()} - Top {target_gene} Mutations")
        if idx == 0:
            ax.legend(loc="lower right")

    plt.savefig(figures_dir / "top_target_mutations.png")
    plt.close()

    # 5. Gene-level mutation burden
    log.info("Computing gene-level mutation burden...")
    ds_rif = load_drug_dataset(processed_dir, "rifampicin")
    burden_rif = gene_burden_analysis(ds_rif.features, ds_rif.labels["label"], ds_rif.metadata)
    burden_rif.to_csv(tables_dir / "gene_burden_rifampicin.csv")

    ds_inh = load_drug_dataset(processed_dir, "isoniazid")
    burden_inh = gene_burden_analysis(ds_inh.features, ds_inh.labels["label"], ds_inh.metadata)
    burden_inh.to_csv(tables_dir / "gene_burden_isoniazid.csv")

    log.info("EDA report generation complete! Figures saved in %s, tables in %s", figures_dir, tables_dir)


if __name__ == "__main__":
    main()

