#!/usr/bin/env bash
set -euo pipefail

echo "==============================================="
echo "🚀 AMR Metagenomics Pipeline "
echo "==============================================="

# -------------------------------
# 0. Directory Setup
# -------------------------------
mkdir -p data/raw_reads data/trimmed_reads databases/megares results visualization

# -------------------------------
# 1. Check micromamba & env setup
# -------------------------------
if ! command -v micromamba &>/dev/null; then
    echo "❌ micromamba not found! Please install micromamba first."
    exit 1
fi

eval "$(micromamba shell hook --shell bash)"

if ! micromamba list -n amrpipeline &>/dev/null; then
    echo "🛠 Creating micromamba environment 'amrpipeline'..."
    micromamba create -y -n amrpipeline \
        -c conda-forge unzip pandas matplotlib seaborn tabulate \
        -c bioconda fastp kma
else
    echo "✅ Environment 'amrpipeline' already exists!"
fi

micromamba activate amrpipeline

# -------------------------------
# 2. Download & setup MegaRes DB
# -------------------------------
read -p "📥 Download & setup MegaRes DB? [Y/n] " mega_ans
mega_ans=${mega_ans:-Y}
if [[ "$mega_ans" =~ ^[Yy]$ ]]; then
    echo "⬇️ Downloading MegaRes database..."
    wget -qO databases/megares/megares_v3.00.zip https://www.meglab.org/downloads/megares_v3.00.zip
    unzip -o databases/megares/megares_v3.00.zip -d databases/megares
    echo "📂 Indexing MegaRes DB..."
    kma index -i databases/megares/megares_v3.00/megares_database_v3.00.fasta \
              -o databases/megares/megares_db
    echo "✅ MegaRes DB ready!"
fi

# -------------------------------
# 3. Prompt for trimming
# -------------------------------
read -p "✂️ Run adapter trimming with fastp? [Y/n] " trim_ans
trim_ans=${trim_ans:-Y}

# -------------------------------
# 4. Process Samples
# -------------------------------
for R1 in data/raw_reads/*_1.fastq.gz; do
    SAMPLE=$(basename "$R1" _1.fastq.gz)
    R2="data/raw_reads/${SAMPLE}_2.fastq.gz"

    if [[ ! -f "$R2" ]]; then
        echo "❌ Skipping $SAMPLE: missing $R2"
        continue
    fi

    echo "==============================================="
    echo "▶️  Processing sample: $SAMPLE"
    echo "==============================================="

    if [[ "$trim_ans" =~ ^[Yy]$ ]]; then
        echo "✂️  Running fastp trimming..."
        fastp -i "$R1" -I "$R2" \
              -o "data/trimmed_reads/${SAMPLE}_1.trimmed.fastq.gz" \
              -O "data/trimmed_reads/${SAMPLE}_2.trimmed.fastq.gz" \
              --detect_adapter_for_pe --thread 4 \
              --html "results/${SAMPLE}_fastp.html" \
              --json "results/${SAMPLE}_fastp.json"
    else
        TRIM1="data/trimmed_reads/${SAMPLE}_1.trimmed.fastq.gz"
        TRIM2="data/trimmed_reads/${SAMPLE}_2.trimmed.fastq.gz"
        if [[ ! -f "$TRIM1" || ! -f "$TRIM2" ]]; then
            echo "❌ Trimmed files for $SAMPLE not found in data/trimmed_reads/"
            echo "Please run trimming or provide trimmed files."
            exit 1
        else
            echo "✅ Found existing trimmed files for $SAMPLE, skipping trimming."
        fi
    fi

    # -------------------------------
    # 5. Run KMA
    # -------------------------------
    echo "🚀 Running KMA (MegaRes)..."
    kma -ipe "data/trimmed_reads/${SAMPLE}_1.trimmed.fastq.gz" \
             "data/trimmed_reads/${SAMPLE}_2.trimmed.fastq.gz" \
        -o "results/${SAMPLE}_megares" \
        -t_db databases/megares/megares_db \
        -mem_mode
    echo "✅ MegaRes .res file created: results/${SAMPLE}_megares.res"
done

# -------------------------------
# 6. Finish
# -------------------------------
echo "==============================================="
echo "🎉 All samples processed!"
echo "👉 Now run: ./amr_report.py for summary & plots"
echo "==============================================="
