# Tick Genomics Project: Rhipicephalus appendiculatus

# 1. Pre-Assembly Quality Control
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
# 2. Transcriptome Assembly with Trinity
After cleaning the data, the next step was to assemble the short reads into full-length transcripts. For ticks, which have complex transcriptomes, we used the Trinity assembler.

#Data Preparation

All cleaned *.fastq.gz files and HTML reports were moved to a dedicated assembly directory to simplify file paths for the assembler.
```bash
# Navigate to working directory
cd ~/TickProject/assembly

# Move cleaned data from raw_data folder to current assembly folder
mv ~/TickProject/raw_data/*_clean_*.fastq.gz .
mv ~/TickProject/raw_data/*.html .
```
#Assembly Strategy & Script

To capture the full repertoire of chemosensory receptors (ORs, IRs, GRs, and binding proteins), we combined the reads from all three replicates (Rapp1, 2, 3). This increases coverage for low-abundance transcripts.

#Bash script 
```bash
#!/bin/bash
#SBATCH --job-name=Rapp_Trinity        # Name of your job
#SBATCH --output=trinity_%j.log        # Standard output and error log
#SBATCH --nodes=1                      # Run on a single node
#SBATCH --ntasks=1                     # Run a single task
#SBATCH --cpus-per-task=16             # Use 16 cores for faster assembly
#SBATCH --mem=64G                      # Request 64GB of RAM
#SBATCH --time=48:00:00                # Give it 2 days to finish

# 1. Define the full path to the Trinity executable
TRINITY_BIN=/opt/apps/trinity/v2.12.0/Trinity

# 2. Define the working directory
WORKDIR=/home/amukami/TickProject/assembly
cd $WORKDIR

# 3. Run Trinity using the full path
# We combine all 3 replicates (Rapp1, 2, 3) for maximum sensitivity
$TRINITY_BIN --seqType fq \
        --max_memory 60G \
        --left Rapp1_clean_1.fastq.gz,Rapp2_clean_1.fastq.gz,Rapp3_clean_1.fastq.gz \
        --right Rapp1_clean_2.fastq.gz,Rapp2_clean_2.fastq.gz,Rapp3_clean_2.fastq.gz \
        --CPU 16 \
        --output $WORKDIR/trinity_out

echo "Trinity Assembly Completed at $(date)"
```
#Execution & Monitoring
```bash
# Submit the job
sbatch submit_assembly.sh

# Monitor status (should be 'R' for Running)
squeue -u amukami

# Watch the log file for progress (Jellyfish, Inchworm, Chrysalis, Butterfly)
tail -f trinity_[JOBID].log
```

# 3.Reference Genome Acquisition

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
. Genome Database Creation
After downloading the reference genome, we need to create a searchable database. This allows us to rapidly query the genome for specific receptors and compare it against our assembled transcripts.

# 4. Genome Database Creation
After downloading the reference genome, we need to create a searchable database. This allows us to rapidly query the genome for specific receptors and compare it against our assembled transcripts.

#Why We Are Creating These Databases

Nucleotide Database (Rhip_app_genome_db): Used to map assembled transcripts back to the genome to determine exon-intron structure of ORs, IRs, and GRs.

Protein Database (Rhip_app_protein_db): Used to directly search for homologs of known binding proteins using faster protein-alignment tools (like BLASTP).

Tool: BLAST+ (Basic Local Alignment Search Tool)
We will use the BLAST+ suite to create BLAST databases from the downloaded sequences.
```bash
# Navigate to the data directory containing the FASTA files
cd ~/TickProject/GenomeData/ncbi_dataset/data/GCA_030522465.2/

# 1. Create a BLAST database from the GENOME assembly (DNA)
makeblastdb -in GCF_030522465.2_Rhipicephalus_appendiculatus_genomic.fna -dbtype nucl -out Rhip_app_genome_db

# 2. Create a BLAST database from the PROTEIN sequences
makeblastdb -in protein.faa -dbtype prot -out Rhip_app_protein_db
```
#Database Verification

We verify the database creation by checking the generated files in the directory.
```bash
# List files to check for .nhr, .nin, .nsq files (nucleotides) 
# and .phr, .pin, .psq files (proteins)
ls -lh Rhip_app_genome_db.*
ls -lh Rhip_app_protein_db.*
```
# 5. Homology Search (Protein)
Used blastx to identify potential receptors in the assembled transcripts against the custom protein database created.

#Running the search on Slurm

Results

receptor_hits.txt: BLAST output file containing potential receptor matches.

The search was run on the HPC cluster using the following script:
```bash
#!/bin/bash
#SBATCH --job-name=blastpt_job
#SBATCH --output=blast_results.log
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=04:00:00

# 1. Source environment modules (Required for Slurm)
source /etc/profile.d/modules.sh

# 2. Load the BLAST module
module load blast

# 3. Define paths based on your directory structure
PROJECT_DIR="/home/amukami/TickProject"
QUERY="$PROJECT_DIR/trinity_out/Trinity.fasta"
DB="$PROJECT_DIR/custom_database_prot"
OUTPUT="$PROJECT_DIR/receptor_hits.txt"

# 4. Run the BLASTX command
# Using 16 threads to match cpus-per-task
blastx -query $QUERY \
       -db $DB \
       -out $OUTPUT \
       -evalue 1e-10 \
       -num_threads 16 \
       -outfmt "6 qseqid sseqid pident length evalue bitscore stitle"
```






