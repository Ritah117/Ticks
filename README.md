 #TICKS
The first step is to  download the reference genome for Rhipicephalus appendiculatus (the brown ear tick), using  the NCBI Datasets tool, which is the modern interface for retrieving genomic data.

As of 2026, there is a high-quality assembly available under the accession GCA_030522465.2.
The download is through powershell 
--------------------------------------------------------------------------------
Follow these commands to create a workspace and download the Rhipicephalus appendiculatus data:
1. Move to a safe workspace
 # Create a folder on your C: drive for your genomic data
 mkdir C:\TickProject
 cd C:\TickProject
# Move the datasets.exe tool you just downloaded into this new folder
 move C:\WINDOWS\system32\datasets.exe C:\TickProject\
2. Run the Download Command
This command pulls the DNA sequence, both annotation formats (GFF3 and GTF), the proteins, and the coding sequences.
.\datasets.exe download genome accession GCA_030522465.2 --include genome,gff3,gtf,protein,cds
3. Unzip the Package
# Unzip the downloaded file into a folder named 'GenomeData'
Expand-Archive -Path ncbi_dataset.zip -DestinationPath .\GenomeData
-------------------------------------------------------------------------------
Verifying Your Files
Once unzipped, you can verify everything is there by navigating to the data folder:

PowerShell
cd .\GenomeData\ncbi_dataset\data\GCA_030522465.2\
dir
