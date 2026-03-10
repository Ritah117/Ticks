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
# 2. Transcriptome Assembly with Trinity
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



# 3.Reference Genome Acquisition
To perform a comparative analysis of chemosensory genes (ORs, IRs, GRs, CSPs, SNMPs, and Binding Proteins), a consolidated reference proteome was constructed. This allows for a standardized search across diverse arthropod lineages.

#3.1 Species Selection 

We selected 10 target arthropod species representing various ecological niches (vectors, agricultural pests, and model organisms).

## Dataset: Arthropods of Medical, Agricultural, and Research Importance

This dataset contains selected **arthropod species of medical, agricultural, and scientific importance**.  
The organisms are grouped into three categories:

- **Vectors of Disease** – Species that transmit pathogens to humans or animals  
- **Agricultural Pests** – Species that cause economic losses in crops, stored products, or livestock  
- **Model Organisms** – Species widely used in biological and genetic research  

### Species and Their Significance

| Category | Species | Common Name | Ecological/Economic Significance |
|---|---|---|---|
| Vectors of Disease | *Rhipicephalus appendiculatus* | Brown Ear Tick | Vector of East Coast Fever (Theileriosis). |
|  | *Aedes aegypti* | Yellow Fever Mosquito | Vector of Zika, Dengue, and Yellow Fever. |
|  | *Anopheles gambiae* | Malaria Mosquito | Primary vector of human malaria in Africa. |
|  | *Glossina morsitans* | Tsetse Fly | Vector of African Trypanosomiasis (Sleeping Sickness). |
|  | *Cimex lectularius* | Bed Bug | Significant urban nuisance pest and blood-feeder. |
| Agricultural Pests | *Stomoxys calcitrans* | Stable Fly | Livestock pest causing significant milk and meat production losses. |
|  | *Tribolium castaneum* | Red Flour Beetle | Major pest of stored grain products worldwide. |
|  | *Cryptotermes secundus* | Drywood Termite | Structural pest that damages timber and buildings. |
| Model Organisms | *Drosophila melanogaster* | Vinegar Fly | Widely used model organism for genetics and olfaction research. |
|  | *Bombyx mori* | Silk Moth | Model organism for Lepidopteran pheromone detection and silk production studies. |





# 3.2 Data Sourcing & Proteome Assembly 

Proteomes were retrieved in FASTA format from UniProt (Reference Proteomes) and NCBI RefSeq. To ensure a standardized comparative analysis, these individual proteomes were concatenated into a single master database.

Master File Name: Master_Arthropod_Reference.fasta

Total Sequence Count: 259,579 proteins

## Protein Sequence Dataset Summary

This table summarizes the number of protein sequences retrieved for each species and the primary biological database used as the source.

| Species | Sequence Count | Primary Source |
|---|---|---|
| *Rhipicephalus appendiculatus* | 56,508 | NCBI RefSeq |
| *Drosophila melanogaster* | 30,802 | UniProt |
| *Cryptotermes secundus* | 29,285 | NCBI RefSeq |
| *Stomoxys calcitrans* | 26,015 | NCBI RefSeq |
| *Cimex lectularius* | 24,194 | NCBI RefSeq |
| *Tribolium castaneum* | 22,610 | UniProt |
| *Bombyx mori* | 22,510 | UniProt |
| *Aedes aegypti* | 20,643 | UniProt |
| *Anopheles gambiae* | 14,102 | UniProt |
| *Glossina morsitans* | 12,910 | UniProt |



#Assembly Workflow (Bash)

The consolidated reference database was generated using standard command-line tools to merge individual species proteomes and verify the final sequence counts.
```bash
# 1. Navigating to the project workspace
cd /home/amukami/R.appendiculatus/insect_chemo_project/raw_proteomes

# 2. Consolidating 10 species-specific proteomes into a single Master Reference
# We used the wildcard *.fasta to capture all downloaded datasets
cat *.fasta > Master_Arthropod_Reference.fasta

# 3. Validation: Verifying sequence counts per species
# This script was used to generate the statistics for Table 1
for proteome in *.fasta; do
    echo "Species File: $proteome"
    grep -c ">" "$proteome"
done

# 4. Final verification of the Master Database
grep -c ">" Master_Arthropod_Reference.fasta
```

# 4. HMMER Analysis & Gene Discovery
To identify the chemosensory repertoire, we performed a profile-based search using HMMER 3.3.2. for detecting highly divergent homologs,  found in the tick genome.


#Targeted Gene Families

We searched for six distinct families essential for arthropod chemical sensing:

Odorant Receptors (OR)

Ionotropic Receptors (IR)

Gustatory Receptors (GR)

Odorant Binding Proteins (OBP)

Chemosensory Proteins (CSP)

Sensory Neuron Membrane Proteins (SNMP)

#The Search Command

We utilized the hmmsearch tool with a stringent E-value cutoff of 1e-5. The following loop was used to automate the search across all six families:

```bash
for family in OR IR GR OBP CSP SNMP; do
    echo "Processing family: $family"
    hmmsearch --tblout results/${family}_hits.tbl \
    --noali -E 1e-5 \
    hmm_profiles/${family}.hmm \
    raw_proteomes/Master_Arthropod_Reference.fasta
done
```

##  Comparative Genomics Results

The discovery phase revealed significant differences in the sensory evolution of ticks versus insects.

###  Summary Table: Identified Gene Counts

| Species | OR | IR | GR | OBP | CSP | SNMP |
|---|---|---|---|---|---|---|
| *Rhipicephalus appendiculatus* | 0 | 189 | 45 | 35 | 0 | 6 |
| *Aedes aegypti* | 85 | 93 | 53 | 17 | 48 | 25 |
| *Anopheles gambiae* | 80 | 46 | 92 | 17 | 7 | 17 |
| *Bombyx mori* | 109 | 52 | 26 | 29 | 43 | 35 |
| *Cimex lectularius* | 65 | 85 | 28 | 21 | 19 | 29 |
| *Cryptotermes secundus* | 104 | 166 | 59 | 52 | 10 | 30 |
| *Drosophila melanogaster* | 64 | 93 | 84 | 36 | 9 | 35 |
| *Glossina morsitans* | 46 | 28 | 15 | 13 | 5 | 15 |
| *Stomoxys calcitrans* | 75 | 116 | 69 | 27 | 13 | 34 |
| *Tribolium castaneum* | 213 | 91 | 176 | 24 | 25 | 36 |


# Data Extraction and Curated Databases
To facilitate downstream analysis, full protein sequences were extracted from the master file using the hit IDs from the HMMER results.

1. The Automated HMMER Search
This code runs the search for OR, IR, GR, OBP, CSP, and SNMP all at once. It looks for the .hmm profiles and searches them against your Master_Arthropod_Reference.fasta.

```bash
# Define the families we are looking for
for family in OR IR GR OBP CSP SNMP; do
    echo "Starting search for: $family"
    
    # Run hmmsearch with the 1e-5 threshold
    # --tblout saves the results in an easy-to-read table format
    hmmsearch --tblout results/${family}_hits.tbl \
    --noali -E 1e-5 \
    hmm_profiles/${family}.hmm \
    raw_proteomes/Master_Arthropod_Reference.fasta
done
```
2. Counting the Hits (To build your Table)
After the search finished, we used this "one-liner" to count how many sequences were found for each family, excluding the header lines (which start with #).

```bash
for file in results/*.tbl; do
    printf "$file: "
    grep -v "^#" "$file" | wc -l
done
```
3. Extracting the FASTA Sequences
Once we had the hit IDs in the .tbl files, we needed to pull the actual protein sequences out of the master file. s

This is the logic we used for all families:

```bash
for family in OR IR GR OBP CSP SNMP; do
    echo "Extracting sequences for $family..."
    
    # 1. Get the IDs from the first column of the HMMER table
    awk '!/^#/ {print $1}' results/${family}_hits.tbl > temp_ids.txt
    
    # 2. Use grep to pull those IDs (and the sequence line below them) 
    # from the master FASTA. We also filter out any contaminants.
    grep -Ff temp_ids.txt -A 1 raw_proteomes/Master_Arthropod_Reference.fasta | \
    grep -v "Solanum lycopersicum" | grep -v "^--" > family_databases/${family}_arthropod.fasta
done

# Clean up temporary file
rm temp_ids.txt
```


# Summary of Curated Databases
After executing the extraction pipeline above, the following family-specific FASTA databases were generated.
##  Sequence Files Summary

This table summarizes the **protein sequence files used for comparative genomics**, including the number of sequences and a brief description of each gene family.

| File Name | Sequence Count | Description |
|---|---|---|
| OR_arthropod.fasta | 841 | Odorant Receptors (Insect-specific) |
| IR_arthropod.fasta | 959 | Ionotropic Receptors (Ancestral system) |
| GR_arthropod.fasta | 647 | Gustatory Receptors (Taste/CO2) |
| OBP_arthropod.fasta | 271 | Odorant Binding Proteins |
| CSP_arthropod.fasta | 179 | Chemosensory Proteins |
| SNMP_arthropod.fasta | 262 | Sensory Neuron Membrane Proteins |




















