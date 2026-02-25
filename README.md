#Ticks

The first step is to  download the reference genome for Rhipicephalus appendiculatus (the brown ear tick), using  the NCBI Datasets tool, which is the modern interface for retrieving genomic data.

As of 2026, there is a high-quality assembly available under the accession GCA_030522465.2.
The download is through powershell 

Follow these commands to create a workspace and download the Rhipicephalus appendiculatus data:
--------------------------------------------------------------------------------

1. Move to a safe workspace

```bash
# Create a folder on your C: drive for your genomic data
mkdir C:\TickProject
cd C:\TickProject

# Move the datasets.exe tool you just downloaded into this new folder
move C:\WINDOWS\system32\datasets.exe C:\TickProject
```
#Run the download command 
This command pulls the DNA sequence, both annotation formats (GFF3 and GTF), the proteins, and the coding sequences.
```bash
.\datasets.exe download genome accession GCA_030522465.2 --include genome,gff3,gtf,protein,cds
```
#Unzip the downloaded file into a folder named 'GenomeData'
```bash
Expand-Archive -Path ncbi_dataset.zip -DestinationPath .\GenomeData
```
#Verifying Your Files
Once unzipped, you can verify everything is there by navigating to the data folder:

```bash
cd .\GenomeData\ncbi_dataset\data\GCA_030522465.2\
dir
```
# Genome assembly
copy the files to the linux working directory and
check for the presence of any corrupted file, 
this is important because chemosensory genes in ticks are often low-abundance.

If its missing half the data, the assembler won't have enough "overlap" to connect the pieces of a Gustatory Receptor (GR).

use this 
```bash
gzip -t *.fastq.gz
```
# Cleaning (QC) and Normalization
#why this is important

1.The sequencer adds artificial DNA adapters to the tick RNA. If you don't remove them, the assembler (Trinity) will mistake them for real biological sequences.

The risk associated with this is that you might "discover" a new protein that is actually just a piece of the Illumina sequencing kit. fastp trims these off so only Tick DNA remains.

#Sequence Quality Control and Data Integrity Validation

Objective: Removal of Illumina adapters and low-quality bases ($Q < 20$) to ensure high-accuracy mapping to the reference genome.

Receptors (GRs, IRs, ORs) are long, complex proteins that weave in and out of the cell membrane seven times (7-Transmembrane domains).

Sequencing quality always drops at the end of a read. These "blurry" bases act like roadblocks for the assembler.

By trimming the low-quality ends, fastp allows the assembler to stitch reads together more smoothly, giving you full-length receptors instead of broken fragments.

#Run this Command for Replicate 1 (Rapp1)
```bash
fastp -i Rapp1_1.fastq.gz -I Rapp1_2.fastq.gz \
      -o clean_Rapp1_1.fq.gz -O clean_Rapp1_2.fq.gz \
      --html Rapp1_report.html --json Rapp1.json \
      --thread 4 --qualified_quality_phred 20
```
#Run this Command for Replicate 2 (Rapp2)
```bash
fastp -i Rapp2_1.fastq.gz -I Rapp2_2.fastq.gz \
      -o clean_Rapp2_1.fq.gz -O clean_Rapp2_2.fq.gz \
      --html Rapp2_report.html --json Rapp2.json \
      --thread 4 --qualified_quality_phred 20
```
#Run this Command for Replicate 3 (Rapp3)
```bash
fastp -i Rapp3_1.fastq.gz -I Rapp3_2.fastq.gz \
      -o clean_Rapp3_1.fq.gz -O clean_Rapp3_2.fq.gz \
      --html Rapp3_report.html --json Rapp3.json \
      --thread 4 --qualified_quality_phred 20
```

Following the pre-processing of each replicate, a formal data integrity check was performed. All compressed output files (.fq.gz) were verified using the gzip -t command to detect any potential truncation or bitstream errors.
```bash
gzip -t clean_Rapp*_*.fq.gz && echo "All chemosensory datasets verified healthy."
```
# Data Integrity & Crash Recovery Note

All processed files were verified for integrity using gzip -t. During the processing of Replicate 3 (Rapp3), a system interruption led to file truncation (error: unexpected end of file).

Resolution:

Corrupted outputs were removed, the process was re-executed with reduced computational overhead (--thread 1) to ensure stability.

All final datasets (Rapp1, Rapp2, and Rapp3) were re-verified and confirmed healthy before proceeding to genome-guided assembly.

