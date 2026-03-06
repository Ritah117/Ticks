# Integrated Gene Discovery and RNA-seq Quantification Pipeline
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
# 2.Reference Genome Acquisition

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

# 3. Genome Database Creation
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

# 4. Transcriptome Assembly with Trinity
After cleaning the data, the next step was to assemble the short reads into full-length transcripts. For ticks, which have complex transcriptomes, we used the Trinity assembler.

we were able to piece together 205,016 complete  transcripts

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


# 5. Genome Mapping (HISAT2)
Raw reads were mapped against the Rhipicephalus NCBI genome to ensure sequence authenticity and confirm the tick-origin of our data.

95.39% Overall Alignment Rate

We compared the mapping efficiency across our three biological replicates (Rapp1, Rapp2, Rapp3) to ensure reproducibility.

Rapp1: 71.08%
Rapp2: 63.53%
Rapp3: 70.42%

We successfully assigned approximately 3,000,000 reads per replicate, proving high data quality and low contamination.
```bash
#!/bin/bash
#SBATCH --job-name=Tick_Genome_Map
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=LOG_genome_map_%j.log

# This line ensures the 'module' command works inside the SLURM environment
. /etc/profile.d/modules.sh || source /etc/profile.d/modules.sh

# 1. Load Tools with exact version names
module load hisat2/2.2.1
module load samtools

# 2. Setup Directories
mkdir -p genome_mapping_results                                                  
# Define paths for your OR/BP discovery
GENOME="references/GCA_030522465.2_ASM3052246v2_genomic.fna"
INDEX="references/ncbi_genome_index"

# 3. Build Index (The searchable map of the genome)
if [ ! -f "${INDEX}.1.ht2" ]; then
    echo "Building NCBI Genome Index..."
    hisat2-build -p 32 $GENOME $INDEX
fi
# 4. Map Clean Reads to Genome
# This step finds where your Trinity transcripts 'live' on the chromosomes
for i in 1 2 3; do
    echo "Mapping Replicate Rapp$i to NCBI Genome..."
    hisat2 -p 32 --dta -x $INDEX \
        -1 Rapp${i}_clean_1.fastq.gz \
        -2 Rapp${i}_clean_2.fastq.gz | \
        samtools view -@ 32 -bS - | \
        samtools sort -@ 32 -o genome_mapping_results/Rapp${i}_NCBI_Aligned.sort$
    # Create the index (.bai) so we can see the receptors in IGV
    samtools index genome_mapping_results/Rapp${i}_NCBI_Aligned.sorted.bam
done

echo "Process Complete. Your BAM files are in the 'genome_mapping_results' folde$
```
# 6. Structural Annotation (Minimap2)
We linked the Trinity transcripts to physical genomic coordinates using a custom-generated .gtf map to determine the physical location and scaffold of each sensory candidate within the tick genome

We got a 72%, 67%, and 72% match across our three samples.

```bash
#!/bin/bash
#SBATCH --job-name=Minimap_Bridge
#SBATCH --partition=debug
#SBATCH --ntasks=30
#SBATCH --mem=40G
#SBATCH --time=02:00:00
#SBATCH --output=log_minimap_%j.log

# 1. Load the specific module we found
source /etc/profile.d/modules.sh
module load minimap/2.2.22

# Define Paths
GENOME="references/GCA_030522465.2_ASM3052246v2_genomic.fna"
TRINITY_FASTA="trinity_out/Trinity.fasta"

# 2. Run Minimap2 Alignment
# -ax splice: optimized for long transcripts vs genome
# -t 30: uses the 16 threads we requested
minimap2 -ax splice -t 30 $GENOME $TRINITY_FASTA > trinity_to_genome.sam

# 3. Quick conversion to BAM (Optional but helpful for size)
# If samtools is available, it makes the file much smaller
module load samtools 2>/dev/null
samtools view -Sb trinity_to_genome.sam > trinity_to_genome.bam

echo "Bridge Step Complete. Trinity transcripts are now mapped to the genome."
```
 # 7. FeatureCounts: Expression Quantification
We calculated raw counts for all genes across the three biological replicates to measure gene activity levels.

Using featureCounts, we linked our 205,016 Trinity transcripts to specific genomic locations to quantify sensory gene expression.

Replicate Rapp3 consistently showed higher raw counts: 40,000 for the top hit vs. ~28,000 in others, reflecting higher sequencing depth or metabolic activity while maintaining a consistent expression profile across all replicates
```bash
#!/bin/bash
#SBATCH --job-name=Sensory_Counts
#SBATCH --partition=debug
#SBATCH --ntasks=16
#SBATCH --mem=30G
#SBATCH --time=01:00:00
#SBATCH --output=log_quantification.log

# 1. Load the software
source /etc/profile.d/modules.sh
module load subread

# 2. Run featureCounts
# -p: paired-end reads
# -T 8: use 8 threads
# -a: your new GTF map
# -o: the output file
featureCounts -p -T 8 \
  -t gene \
  -g gene_id \
  -a trinity_map.gtf \
  -o sensory_gene_counts.txt \
  genome_mapping_results/Rapp1_NCBI_Aligned.sorted.bam \
  genome_mapping_results/Rapp2_NCBI_Aligned.sorted.bam \
  genome_mapping_results/Rapp3_NCBI_Aligned.sorted.bam

echo "Quantification complete. Check sensory_gene_counts.txt"

```







