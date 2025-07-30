AMR Gene Detection Metagenomics Pipeline 
==========================================

This pipeline processes paired-end metagenomic reads, detects AMR genes using the 
MegaRes database, and generates detailed summaries and plots.

-----------------------------------
1) Requirements
-----------------------------------
put your raw reads in data/raw_reads folder.

-----------------------------------
2) Running the Pipeline
-----------------------------------
Step 1: Run the main pipeline 
    ./amr_pipeline.sh

   - Prompts to download MegaRes DB (if not already present)
   - Prompts to run fastp trimming (or uses existing trimmed files)
   - Runs KMA and generates .res files in results/

Step 2: Run the reporting script
    ./amr_report.py

   - Reads all .res files
   - Prints 3 tables:
        1. Drug Class Abundance Summary
        2. Per-sample AMR Abundance
        3. Top 20 Gene/Mechanism Abundance
   - Prints a short summary for each sample
   - Saves CSV tables and plots in visualization/

-----------------------------------
3) Output Directories
-----------------------------------
data/raw_reads/       -> input FASTQs
data/trimmed_reads/   -> trimmed FASTQs (if trimming enabled)
databases/megares/    -> MegaRes DB and index
results/              -> KMA outputs (.res files)
visualization/        -> CSVs and plots

-----------------------------------
Example summary output
-----------------------------------
Summary for sample 'SRR13172349':
Total AMR abundance (depth): 2059.78
Dominant drug classes detected:
  - MLS (33.9%)
  - Aminoglycosides (27.8%)
  - Tetracyclines (14.7%)
Top detected resistance genes:
  - MEG_4223|Drugs|Oxazolidinone|...
  - MEG_3988|Drugs|MLS|...
  - MEG_3992|Drugs|MLS|...
