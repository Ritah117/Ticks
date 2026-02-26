# Tick Genomics Project: Rhipicephalus appendiculatus

# 1.Reference Genome Acquisition

The first step was to download the reference genome for Rhipicephalus appendiculatus (the brown ear tick) using the NCBI Datasets tool. As of 2026, we are using the high-quality assembly GCA_030522465.2.

#Workspace Setup & Download (Windows PowerShell)
```bash
# Create a workspace on the D: drive
mkdir D:\TickGenome
cd D:\TickGenome

# Download the DNA, GFF3, GTF, proteins, and coding sequences
.\datasets.exe download genome accession GCA_030522465.2 --include genome,gff3,gtf,protein,cds

# Unzip the data
Expand-Archive -Path ncbi_dataset.zip -DestinationPath .\GenomeData

# Verification of files
cd .\GenomeData\ncbi_dataset\data\GCA_030522465.2\
dir
```
# 2. Pre-Assembly Quality Control
Before assembly, we moved the raw sequencing data to the HPC and verified the files.

# Data Integrity Check

Because chemosensory genes in ticks are often low-abundance, a corrupted file could lead to missing data. If the file is missing data, the assembler won't have enough "overlap" to connect the pieces of a Gustatory Receptor (GR) or Odorant Receptor (OR).

```bash
# Check for corrupted compressed files
gzip -t *.fastq.gz
```
# Sequence Cleaning with fastp

Objective: Removal of Illumina adapters and low-quality bases ($Q < 20$).

Why this is critical for Receptors:

Structural Complexity: Receptors are long, complex proteins that weave in and out of the cell membrane seven times (7-Transmembrane domains).

Avoiding Chimeras: Trimming adapters prevents the assembler from creating "chimera" sequences—artificial hybrids of tick DNA and sequencing adapters.

Read Smoothing: Trimming "blurry" bases at the ends of reads allows the assembler to stitch sequences together more accurately, giving us full-length receptors instead of fragments.

#The Cleaning Command

We used a bash loop to process Rapp1, Rapp2, and Rapp3 in one go. This ensured identical parameters were applied to all replicates.

```bash
# Navigate to the data folder
cd ~/TickProject/raw_data

# Run fastp loop
for i in 1 2 3; do
    fastp \
    -i Rapp${i}_1.fastq.gz \
    -I Rapp${i}_2.fastq.gz \
    -o Rapp${i}_clean_1.fastq.gz \
    -O Rapp${i}_clean_2.fastq.gz \
    -h Rapp${i}_report.html \
    -j fastp.json \
    --qualified_quality_phred 20 \
    --thread 4
done
```
#Data Archiving & Backup

After cleaning, we secured the data

```bash
# Backup cleaned FASTQ files
tar -cvzf ~/TickProject_Clean_Backup.tar.gz ~/TickProject/raw_data/*_clean.fastq.gz

# Backup HTML reports and the JSON summary
tar -cvzf ~/Rapp_Reports_Backup.tar.gz ~/TickProject/raw_data/*.html ~/fastp.json
```







