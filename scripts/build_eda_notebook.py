"""Generate the student-friendly 01_data_exploration.ipynb notebook."""

import json
from pathlib import Path

cells = [
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "# 01. Exploratory Data Analysis: Predicting TB Drug Resistance from Genomic SNPs\n",
            "\n",
            "**Course / Context:** University Bioinformatics Unit for Computer Science Students  \n",
            "**Core Question:** Can single nucleotide polymorphisms (SNPs) in known resistance genes accurately predict whether a tuberculosis isolate is resistant to antibiotics?\n",
            "\n",
            "---\n",
            "\n",
            "### The Big Picture (From a CS Perspective)\n",
            "In this project, we are tackling a **binary classification problem** using real-world pathogen genomics from the **CRyPTIC compendium**:\n",
            "- **The Samples (Rows):** Each row represents an isolate (a bacterial sample taken from a TB patient).\n",
            "- **The Features ($X$):** Binary flags (`1` or `0`) indicating whether a specific genetic mutation is present or absent in known target genes.\n",
            "- **The Ground Truth Labels ($y$):** Laboratory culture tests (`1` = Resistant, `0` = Susceptible) across four antibiotics.\n",
            "\n",
            "### Target Antibiotics & Genes in our Panel\n",
            "1. **First-line drugs (Core treatment):**\n",
            "   - **Rifampicin (RIF):** Targets RNA polymerase subunit B (`rpoB`).\n",
            "   - **Isoniazid (INH):** Prodrug activated by catalase-peroxidase (`katG`), targeting fatty acid synthase enzymes (`inhA` and the `fabG1` promoter).\n",
            "2. **Second-line drugs (Used when first-line fails):**\n",
            "   - **Levofloxacin (LEV)** & **Moxifloxacin (MXF):** Fluoroquinolones targeting DNA gyrase subunit A (`gyrA`).\n",
            "\n",
            "### Goals of this EDA\n",
            "Before jumping into training machine learning models (like Logistic Regression or Random Forests), we need to inspect the data:\n",
            "1. **Class Balance:** Are resistant and susceptible cases balanced, or do we have class imbalance?\n",
            "2. **Multi-Drug Resistance (MDR):** Do resistances co-occur? Could a model cheat by using an INH mutation to predict RIF resistance?\n",
            "3. **Strongest Signals:** Which individual SNPs act as clear 'smoking guns' for resistance?\n",
            "4. **The Lineage Confounder:** Are some mutations just family markers (phylogenetic lineage markers) rather than causal resistance drivers?\n",
            "5. **Sparsity & Noise:** How many features are zeros, and are there missing values?"
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# Step 1: Import standard scientific python libraries and our project helpers\n",
            "from pathlib import Path\n",
            "import numpy as np\n",
            "import pandas as pd\n",
            "import matplotlib.pyplot as plt\n",
            "import seaborn as sns\n",
            "\n",
            "from tb_resistance.eda import (\n",
            "    load_drug_dataset,\n",
            "    summarize_drug_cohorts,\n",
            "    top_mutations_by_resistance,\n",
            "    lineage_distribution,\n",
            "    gene_burden_analysis,\n",
            "    multidrug_overlap,\n",
            ")\n",
            "\n",
            "# Configure clean plots\n",
            "%matplotlib inline\n",
            "sns.set_theme(style='whitegrid', font='sans-serif')\n",
            "plt.rcParams['figure.dpi'] = 120\n",
            "\n",
            "# Processed data directory containing the per-drug matrices\n",
            "DATA_DIR = Path('../data/processed')\n",
            "print('Ready to explore! Processed data directory:', DATA_DIR.resolve())"
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "## 1. Dataset Overview: Sample Sizes & Label Balance\n",
            "\n",
            "In machine learning, checking label balance is step one. If 95% of samples are susceptible, a dumb model predicting all zeros gets 95% accuracy while being clinically useless!\n",
            "\n",
            "Let's check the number of isolates, the number of resistant ($y=1$) vs susceptible ($y=0$) cases, and the resistance rate for each drug."
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "cohort_summary = summarize_drug_cohorts(DATA_DIR)\n",
            "cohort_summary"
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# Plot sample counts and class balance\n",
            "fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4.5))\n",
            "x = np.arange(len(cohort_summary))\n",
            "width = 0.35\n",
            "\n",
            "# Left plot: Raw counts of Resistant vs Susceptible\n",
            "ax1.bar(x - width/2, cohort_summary['resistant'], width, label='Resistant (1)', color='#c93b2b', alpha=0.9)\n",
            "ax1.bar(x + width/2, cohort_summary['susceptible'], width, label='Susceptible (0)', color='#4a7c59', alpha=0.9)\n",
            "ax1.set_ylabel('Number of Isolates')\n",
            "ax1.set_title('Cohort Size & Class Breakdown')\n",
            "ax1.set_xticks(x)\n",
            "ax1.set_xticklabels([d.capitalize() for d in cohort_summary.index])\n",
            "ax1.legend()\n",
            "\n",
            "# Right plot: Resistance prevalence (% resistant)\n",
            "sns.barplot(data=cohort_summary.reset_index(), x='drug', y='prevalence_pct', ax=ax2, hue='drug', palette='mako', legend=False)\n",
            "ax2.set_ylabel('Resistant Isolates (%)')\n",
            "ax2.set_title('Resistance Prevalence Rate')\n",
            "ax2.set_xticks(range(len(cohort_summary)))\n",
            "ax2.set_xticklabels([d.capitalize() for d in cohort_summary.index])\n",
            "for i, val in enumerate(cohort_summary['prevalence_pct']):\n",
            "    ax2.text(i, val + 1, f'{val:.1f}%', ha='center', fontweight='bold')\n",
            "\n",
            "plt.tight_layout()\n",
            "plt.show()"
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "### Student Notes on Section 1:\n",
            "- **Sample size:** Each drug has between 11,500 and 12,300 tested isolates—a substantial sample size for training ML models.\n",
            "- **First-line drugs (RIF & INH):** Moderate class balance (~39% resistant for Rifampicin, ~46% for Isoniazid). This means standard cross-entropy / log-loss will train well without extreme weighting.\n",
            "- **Second-line fluoroquinolones (LEV & MXF):** Lower prevalence (~18% and ~14%). When evaluating models for LEV and MXF, we should focus on **Precision-Recall curves (PR-AUC)** and **Recall/Sensitivity** rather than simple accuracy."
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "## 2. Multi-Drug Resistance (MDR-TB) and Label Correlation\n",
            "\n",
            "In clinical practice, tuberculosis is treated with combination therapy. Patients take Rifampicin and Isoniazid together. When bacteria survive both drugs, it is classified as **Multi-Drug Resistant TB (MDR-TB)**.\n",
            "\n",
            "**Why does this matter to CS students?**  \n",
            "Because our feature matrix contains all panel genes (`rpoB`, `katG`, `inhA`, `gyrA`), if an isolate has both RIF and INH resistance, an INH mutation like `katG@S315T` will be strongly correlated with RIF resistance!  \n",
            "In data science terms, this is **multicollinearity and label leakage across tasks**. Let's measure how many isolates share resistance."
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# Compute shared isolates and co-resistant counts\n",
            "shared, both_res = multidrug_overlap(DATA_DIR)\n",
            "print('--- Number of Shared Isolates Tested for Both Drugs ---')\n",
            "display(shared)\n",
            "print('\\n--- Number of Isolates Resistant to BOTH Drugs ---')\n",
            "display(both_res)"
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# Visualize co-resistance rate as a heatmap\n",
            "concordance = (both_res / shared * 100.0).astype(float)\n",
            "plt.figure(figsize=(7, 5))\n",
            "sns.heatmap(concordance, annot=True, fmt='.1f', cmap='YlOrRd', cbar_kws={'label': '% Co-Resistant'})\n",
            "plt.title('Pairwise Co-Resistance Concordance (%)')\n",
            "plt.show()"
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "### Student Notes on Section 2:\n",
            "- **MDR-TB is common in this dataset:** Over **4,190 isolates** are resistant to both Rifampicin and Isoniazid!\n",
            "- **Fluoroquinolones travel together:** Levofloxacin and Moxifloxacin show almost identical resistance patterns (~14-17% co-resistance), which makes biological sense because they belong to the exact same drug class and target the same protein (`gyrA`).\n",
            "- **Modeling takeaway:** If we train a model to predict Rifampicin resistance using all panel genes, we should verify whether the model is learning true biological mechanisms (`rpoB`) or just picking up `katG` mutations that happen to co-occur due to MDR-TB!"
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "## 3. Gene-Level Mutation Burden\n",
            "\n",
            "Before looking at thousands of individual SNPs, let's take a bird's-eye view:  \n",
            "**If an isolate is resistant to a drug, does it carry at least one mutation in that drug's known target gene? And what about susceptible isolates?**"
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "ds_rif = load_drug_dataset(DATA_DIR, 'rifampicin')\n",
            "print('=== RIFAMPICIN: Mutation Burden across Genes ===')\n",
            "burden_rif = gene_burden_analysis(ds_rif.features, ds_rif.labels['label'], ds_rif.metadata)\n",
            "display(burden_rif)\n",
            "\n",
            "ds_inh = load_drug_dataset(DATA_DIR, 'isoniazid')\n",
            "print('\\n=== ISONIAZID: Mutation Burden across Genes ===')\n",
            "burden_inh = gene_burden_analysis(ds_inh.features, ds_inh.labels['label'], ds_inh.metadata)\n",
            "display(burden_inh)"
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "### Student Notes on Section 3:\n",
            r"- Look at `rpoB` for Rifampicin: **96.6%** of resistant isolates have at least one mutation in `rpoB`, compared to only 14.8% of susceptible ones. That is a massive difference ($\Delta = 81.8\%$), showing that `rpoB` contains a powerful predictive signal!\n",
            "- Look at `katG` for Isoniazid: **91.2%** of resistant isolates carry a `katG` mutation, compared to 44.8% of susceptible ones.\n",
            "- Notice that even in susceptible isolates, some genes have mutation rates ~40%. Why? We will uncover this in Section 5—it is caused by **phylogenetic lineage markers**!"
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "## 4. Hunting for the 'Smoking Gun' Mutations (Top SNPs)\n",
            "\n",
            "Now let's zoom in from genes to individual mutations.  \n",
            "We rank each mutation by the difference in frequency between Resistant and Susceptible isolates ($\\Delta\\% = \\%_{Resistant} - \\%_{Susceptible}$) and compute the **Odds Ratio** (how many times more likely a resistant isolate is to have this mutation compared to a susceptible one)."
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# Extract top mutations in each drug's primary target gene\n",
            "drugs_and_genes = [\n",
            "    ('rifampicin', 'rpoB'),\n",
            "    ('isoniazid', 'katG'),\n",
            "    ('levofloxacin', 'gyrA'),\n",
            "    ('moxifloxacin', 'gyrA')\n",
            "]\n",
            "\n",
            "for drug, gene in drugs_and_genes:\n",
            "    ds = load_drug_dataset(DATA_DIR, drug)\n",
            "    top = top_mutations_by_resistance(ds.features, ds.labels['label'], ds.metadata, gene=gene, top_n=6)\n",
            "    print(f'\\n=== TOP {gene.upper()} MUTATIONS FOR {drug.upper()} ===')\n",
            "    display(top[['gene', 'mutation', 'resistant_count', 'resistant_freq_pct', 'susceptible_freq_pct', 'delta_pct', 'odds_ratio']])"
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# Bar plots comparing mutation frequencies between classes\n",
            "fig, axes = plt.subplots(2, 2, figsize=(14, 9))\n",
            "axes = axes.flatten()\n",
            "\n",
            "for idx, (drug, gene) in enumerate(drugs_and_genes):\n",
            "    ds = load_drug_dataset(DATA_DIR, drug)\n",
            "    top = top_mutations_by_resistance(ds.features, ds.labels['label'], ds.metadata, gene=gene, top_n=6)\n",
            "    df_p = top.reset_index()\n",
            "    ax = axes[idx]\n",
            "    y_idx = np.arange(len(df_p))\n",
            "    \n",
            "    ax.barh(y_idx - 0.18, df_p['resistant_freq_pct'], height=0.36, color='#c93b2b', label='Resistant %' if idx==0 else '')\n",
            "    ax.barh(y_idx + 0.18, df_p['susceptible_freq_pct'], height=0.36, color='#4a7c59', label='Susceptible %' if idx==0 else '')\n",
            "    ax.set_yticks(y_idx)\n",
            "    ax.set_yticklabels(df_p['mutation'])\n",
            "    ax.invert_yaxis()\n",
            "    ax.set_xlabel('Frequency (%)')\n",
            "    ax.set_title(f'{drug.capitalize()} ({gene})')\n",
            "    if idx == 0:\n",
            "        ax.legend(loc='lower right')\n",
            "\n",
            "plt.tight_layout()\n",
            "plt.show()"
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "### Student Notes on Section 4:\n",
            "Look at the specific SNPs identified—these are famous in tuberculosis literature:\n",
            "- **Rifampicin:** `rpoB@S450L` (serine at position 450 replaced by leucine) appears in **62.1%** of all resistant isolates, but only **1.5%** of susceptible ones. That is an Odds Ratio of **104.8**! Other key players are `rpoB@D435V` and `rpoB@H445D`.\n",
            "- **Isoniazid:** `katG@S315T` is present in **86.4%** of resistant isolates and only **1.2%** of susceptible ones. Its Odds Ratio is over **500**! It is a near-deterministic predictor.\n",
            "- **Fluoroquinolones:** `gyrA@D94G` and `gyrA@A90V` account for the majority of resistance.\n",
            "\n",
            "**Conclusion for ML:** The biological signal in the genomic data is strong and concentrated in a handful of high-impact hotspot mutations."
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "## 5. The 'Lineage Trap': Beware of Confounding Variables!\n",
            "\n",
            "Bacteria reproduce asexually by dividing. Over thousands of years, distinct branches or 'families' have evolved—known as **lineages** (e.g. Lineage 2 East Asian/Beijing, Lineage 4 Euro-American).\n",
            "\n",
            "**Why is this a trap for Machine Learning?**\n",
            "1. **Lineage Markers are NOT Resistance Mutations:**  \n",
            "   Take `katG@R463L`. It is present in ~44% of resistant isolates, but ALSO in ~44% of susceptible isolates! It does not cause resistance; it is simply a genetic marker that defines Lineage 2 and Lineage 4. If a naive model gives it a high weight, it is confusing ancestry with drug resistance.\n",
            "2. **Geographical Outbreaks Cause Resistance Clumping:**  \n",
            "   Certain lineages have experienced major drug-resistant outbreaks in specific countries. Let's inspect the resistance rate per lineage."
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "ds_rif = load_drug_dataset(DATA_DIR, 'rifampicin')\n",
            "lin_df = lineage_distribution(ds_rif.labels, ds_rif.lineage)\n",
            "lin_df"
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "major_lin = lin_df[lin_df['total_isolates'] >= 50].reset_index()\n",
            "fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))\n",
            "\n",
            "# Lineage representation\n",
            "sns.barplot(data=major_lin, x='lineage', y='share_of_cohort_pct', ax=ax1, hue='lineage', palette='crest', legend=False)\n",
            "ax1.set_title('Lineage Share in Dataset (%)')\n",
            "ax1.set_ylabel('% of Isolates')\n",
            "ax1.set_xticks(range(len(major_lin)))\n",
            "ax1.set_xticklabels(major_lin['lineage'], rotation=25)\n",
            "\n",
            "# Lineage resistance rate\n",
            "sns.barplot(data=major_lin, x='lineage', y='prevalence_pct', ax=ax2, hue='lineage', palette='flare', legend=False)\n",
            "ax2.set_title('Rifampicin Resistance Rate by Lineage (%)')\n",
            "ax2.set_ylabel('Prevalence (% Resistant)')\n",
            "ax2.set_xticks(range(len(major_lin)))\n",
            "ax2.set_xticklabels(major_lin['lineage'], rotation=25)\n",
            "\n",
            "plt.tight_layout()\n",
            "plt.show()"
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "### Student Notes on Section 5:\n",
            "- **Lineage 2 (Beijing)** has a **50.8%** resistance rate, while **Lineage 1** has only **21.7%**.\n",
            "- **The Data Leakage Risk:** If we use standard random train/test splits (`train_test_split(shuffle=True)`), closely related isolates from the same transmission cluster will be in both the training set and the test set. The model can achieve inflated test accuracy by memorizing lineage markers!\n",
            "- **Best Practice:** We should use **lineage-stratified cross-validation** to test whether our models generalize across different bacterial families."
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "## 6. Data Quality, Missingness & Heteroresistance\n",
            "\n",
            "In sequencing pipelines, sometimes a genomic position cannot be called with high confidence (e.g. low sequencing depth). In CRyPTIC, these are flagged as null or filter-failed calls.\n",
            "\n",
            "Additionally, an isolate might have **heteroresistance** (a mixed population where some bacterial cells have the mutation and some have the wild-type).\n",
            "\n",
            "Let's check how clean our data is."
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "ds_rif_q = load_drug_dataset(DATA_DIR, 'rifampicin', load_quality_matrices=True)\n",
            "n_missing_isolates = (ds_rif_q.missing.sum(axis=1) > 0).sum()\n",
            "n_het_isolates = (ds_rif_q.het.sum(axis=1) > 0).sum()\n",
            "\n",
            "print(f'Total Rifampicin Isolates: {len(ds_rif_q.labels):,}')\n",
            "print(f'Isolates with >=1 missing/null feature call: {n_missing_isolates:,} ({n_missing_isolates / len(ds_rif_q.labels) * 100:.2f}%)')\n",
            "print(f'Isolates with >=1 heterozygous call: {n_het_isolates:,} ({n_het_isolates / len(ds_rif_q.labels) * 100:.2f}%)')"
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "### Student Notes on Section 6:\n",
            "- Only **~5.7%** of isolates have any missing or low-confidence calls in the panel, and heteroresistance is essentially zero in the finalized matrix.\n",
            "- This means the binary feature matrix is remarkably clean and ready for machine learning algorithms without needing complex imputation techniques."
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "## 7. Summary & Next Steps: Plan for Building Our Models\n",
            "\n",
            "Here is our project roadmap based on what we learned in this EDA:\n",
            "\n",
            "| Observation in EDA | Impact on Machine Learning | Planned Implementation |\n",
            "| :--- | :--- | :--- |\n",
            "| **High Dimensionality & Sparsity** | 6,526 mutation features, but ~5,000 are all zeros in each drug cohort. | Apply rare-variant filtering (e.g., mutation count $\\ge 3$) **fitted strictly on the training fold** to prevent data leakage. |\n",
            "| **Dominant Hotspot SNPs** | Mutations like `rpoB@S450L` and `katG@S315T` have massive odds ratios (>100). | A linear model (**Logistic Regression**) with L1/L2 penalty will provide a strong, interpretable baseline. |\n",
            "| **Phylogenetic Lineage Confounding** | Resistance rates vary drastically by lineage (50.8% in L2 vs 21.7% in L1), and benign lineage markers exist. | Implement **Lineage-stratified cross-validation** to prevent memorizing lineages and test true generalization. |\n",
            "| **MDR-TB Co-occurrence** | 4,190 isolates are resistant to both RIF and INH. | Compare models trained on **target-gene-only features** (e.g., `rpoB` for RIF) vs **full-panel features** to detect shortcut learning. |\n",
            "| **Class Imbalance in Fluoroquinolones** | Only ~14-18% of isolates are resistant to LEV/MXF. | Prioritize **ROC-AUC, PR-AUC, and Sensitivity/Recall** over accuracy. |"
        ]
    }
]

nb = {
    "cells": cells,
    "metadata": {
        "kernelspec": {
            "display_name": "Python 3 (.venv)",
            "language": "python",
            "name": "python3"
        },
        "language_info": {
            "name": "python",
            "version": "3.14.7"
        }
    },
    "nbformat": 4,
    "nbformat_minor": 5
}

Path("notebooks/01_data_exploration.ipynb").write_text(json.dumps(nb, indent=2))
print("Wrote student-friendly notebooks/01_data_exploration.ipynb successfully.")
