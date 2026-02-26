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

