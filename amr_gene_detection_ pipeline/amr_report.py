#!/usr/bin/env -S micromamba run -n amrpipeline python

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import glob
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from tabulate import tabulate

os.makedirs("visualization", exist_ok=True)

# -------------------------------
# 1. Parse .res files
# -------------------------------
def parse_res_files(pattern="results/*_megares.res"):
    files = glob.glob(pattern)
    if not files:
        print("❌ No .res files found in results/! Run amr_pipeline.sh first.")
        exit(1)
    dfs = []
    for f in files:
        sample = os.path.basename(f).replace("_megares.res", "")
        df = pd.read_csv(f, comment='#', sep='\t', header=None,
                         names=['Template','Score','Expected','Template_length',
                                'Template_Identity','Template_Coverage','Query_Identity',
                                'Query_Coverage','Depth','q_value','p_value'])
        df['Sample'] = sample
        dfs.append(df)
    return pd.concat(dfs, ignore_index=True)

# -------------------------------
# 2. Data aggregation functions
# -------------------------------
def drug_class_abundance_summary(df):
    df['Drug_Class'] = df['Template'].apply(lambda x: x.split('|')[2])
    class_depth = df.groupby('Drug_Class')['Depth'].sum().reset_index()
    class_depth = class_depth.sort_values(by='Depth', ascending=False)
    total_depth = class_depth['Depth'].sum()
    class_depth['Percentage(%)'] = (class_depth['Depth'] / total_depth * 100).round(2)
    class_depth.rename(columns={'Depth': 'Total_Depth'}, inplace=True)
    return class_depth, total_depth

def per_sample_amr_abundance(df):
    sample_depth = df.groupby('Sample')['Depth'].sum().reset_index()
    sample_depth.rename(columns={'Depth':'Total_Depth'}, inplace=True)
    return sample_depth

def detailed_gene_abundance(df):
    gene_depth = df.groupby('Template')['Depth'].sum().reset_index()
    gene_depth.rename(columns={'Depth':'Total_Depth', 'Template':'Gene'}, inplace=True)
    return gene_depth.sort_values(by='Total_Depth', ascending=False).head(20)

# -------------------------------
# 3. Plotting
# -------------------------------
def plot_drug_class_abundance(df):
    plt.figure(figsize=(10, max(4, 0.3*len(df))))
    sns.barplot(data=df, y='Drug_Class', x='Total_Depth', palette='viridis', orient='h')
    plt.title("Drug Class Abundance Summary")
    plt.tight_layout()
    plt.savefig("visualization/drug_class_abundance_summary.png")
    plt.close()

# -------------------------------
# 4. Detailed text summary
# -------------------------------
def generate_summary(drug_class_df, sample_abundance_df, gene_df):
    sample_name = sample_abundance_df['Sample'].iloc[0]
    total_abundance = sample_abundance_df['Total_Depth'].iloc[0]

    top_classes = drug_class_df.head(3)
    top_genes = gene_df.head(3)

    print(f"\nSummary for sample '{sample_name}':")
    print(f"Total AMR abundance (depth): {total_abundance:.2f}")
    print("Dominant drug classes detected:")
    for _, row in top_classes.iterrows():
        print(f"  - {row['Drug_Class']} ({row['Percentage(%)']}%)")
    print("Top detected resistance genes:")
    for gene in top_genes['Gene']:
        print(f"  - {gene}")

# -------------------------------
# 5. Main
# -------------------------------
def main():
    print("Parsing .res files...")
    df = parse_res_files()

    samples = df['Sample'].unique()
    for i, sample in enumerate(samples):
        print("\n" + "="*60)
        print(f"🔬 Processing sample: {sample}")
        print("="*60)

        sample_df = df[df['Sample'] == sample]

        drug_class_df, _ = drug_class_abundance_summary(sample_df)
        sample_abundance_df = per_sample_amr_abundance(sample_df)
        gene_df = detailed_gene_abundance(sample_df)

        # Save CSVs
        drug_class_df.to_csv(f"visualization/{sample}_drug_class_abundance.csv", index=False)
        sample_abundance_df.to_csv(f"visualization/{sample}_sample_abundance.csv", index=False)
        gene_df.to_csv(f"visualization/{sample}_gene_abundance_top20.csv", index=False)

        # Plot only for first sample (drug class summary)
        if i == 0:
            plot_drug_class_abundance(drug_class_df)

        # Print beautiful tables
        print("\n📊 Drug Class Abundance Summary:")
        print(tabulate(drug_class_df, headers='keys', tablefmt='psql', showindex=False))

        print("\n📊 Per-sample AMR Abundance:")
        print(tabulate(sample_abundance_df, headers='keys', tablefmt='psql', showindex=False))

        print("\n📊 Top 20 Gene/Mechanism Abundance:")
        print(tabulate(gene_df, headers='keys', tablefmt='psql', showindex=False))

        # Print detailed summary
        generate_summary(drug_class_df, sample_abundance_df, gene_df)

    print("\n" + "="*60)
    print("🎉 Reporting complete!")
    print("👉 Check 'visualization/' for CSVs and plots.")
    print("="*60)

if __name__ == "__main__":
    main()
