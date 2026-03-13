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

We selected 10 target arthropod species representing various ecological niches The dataset contains protein sequences from 10 selected arthropod species chosen for their medical, veterinary, and scientific importance.

# Downloading from UniProt (Reference Proteomes)
For species like Drosophila, Anopheles, and Aedes, UniProt provides "Reference Proteomes" which are the most curated versions.

```bash
# Create a folder for the raw downloads
mkdir -p raw_proteomes
cd raw_proteomes

# Drosophila melanogaster (UP000000589)
wget -O Drosophila_melanogaster.fasta "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28proteome%3AUP000000589%29"

# Anopheles gambiae (UP000000722)
wget -O Anopheles_gambiae.fasta "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28proteome%3AUP000000722%29"

# Aedes aegypti (UP000006203)
wget -O Aedes_aegypti.fasta "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28proteome%3AUP000006203%29"

# Bombyx mori (UP000005204)
wget -O Bombyx_mori.fasta "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28proteome%3AUP000005204%29"

# Glossina morsitans (UP000000356)
wget -O Glossina_morsitans.fasta "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28proteome%3AUP000000356%29"
```
# Downloading from NCBI RefSeq (FTP)
For the Ticks and some specialized pests, NCBI RefSeq often has more recently updated assemblies than UniProt.

```bash
# Rhipicephalus appendiculatus (Brown Ear Tick)
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/135/715/GCF_004135715.1_ASM413571v1/GCF_004135715.1_ASM413571v1_protein.faa.gz
gunzip GCF_004135715.1_ASM413571v1_protein.faa.gz
mv GCF_004135715.1_ASM413571v1_protein.faa Rhipicephalus_appendiculatus.fasta

# Stomoxys calcitrans (Stable Fly)
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/272/065/GCF_003272065.1_S_calcitrans_v1.0/GCF_003272065.1_S_calcitrans_v1.0_protein.faa.gz
gunzip GCF_003272065.1_S_calcitrans_v1.0_protein.faa.gz
mv GCF_003272065.1_S_calcitrans_v1.0_protein.faa Stomoxys_calcitrans.fasta

# Cimex lectularius (Bed Bug)
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/648/675/GCF_000648675.1_Clec_2.0/GCF_000648675.1_Clec_2.0_protein.faa.gz
gunzip GCF_000648675.1_Clec_2.0_protein.faa.gz
mv GCF_000648675.1_Clec_2.0_protein.faa Cimex_lectularius.fasta
```
# Getting the Tsetse Trio (G. fuscipes & G. brevipalpis)
These are often found via the VectorBase section of UniProt or NCBI.

```bash
# Glossina fuscipes
wget -O Glossina_fuscipes.fasta "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28taxonomy_id%3A44487%29"

# Glossina brevipalpis
wget -O Glossina_brevipalpis.fasta "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28taxonomy_id%3A37029%29"
4. Verification of the Downloads
Once you have all 10 files, you should run a quick loop to check that no file is empty and to get the counts for your table:
```
```bash
for file in *.fasta; do
    echo -n "$file: "
    grep -c ">" "$file"
done
```
## Dataset: Arthropod Species of Medical, Veterinary, and Research Importance
The organisms are grouped into three strategic categories:

Medical and Veterinary Vectors – Blood-feeding species responsible for the transmission of human and animal pathogens, including the Brown Ear Tick (East Coast Fever) and a diverse Tsetse fly trio (Trypanosomiasis).

Livestock and Urban Pests – Species such as the Stable Fly and Bed Bug that cause significant economic losses in animal production or represent major urban nuisance challenges.

Biological Model Organisms – Highly characterized species like the Vinegar Fly and Silk Moth, which serve as the gold standard for understanding insect olfaction and pheromone detection. 

### Species and Their Significance

This table summarizes key arthropod species included in the study, their ecological category, and their biological or economic significance.

| Category | Species | Common Name | Ecological/Economic Significance |
|---|---|---|---|
| Acarine Vector | *Rhipicephalus appendiculatus* | Brown Ear Tick | Primary vector of East Coast Fever (Theileriosis) in African cattle. |
| Dipteran Vectors | *Glossina morsitans* | Savannah Tsetse Fly | Major vector of Nagana in cattle and Sleeping Sickness in humans. |
|  | *Glossina fuscipes* | Riverine Tsetse Fly | Responsible for over 90% of human Sleeping Sickness cases in Africa. |
|  | *Glossina brevipalpis* | Forest Tsetse Fly | Large-bodied vector contributing to animal trypanosomiasis transmission. |
|  | *Anopheles gambiae* | Malaria Mosquito | Primary vector of human malaria (*Plasmodium falciparum*) in sub-Saharan Africa. |
|  | *Aedes aegypti* | Yellow Fever Mosquito | Global vector of Zika, Dengue, and Yellow Fever viruses. |
| Parasites | *Stomoxys calcitrans* | Stable Fly | Significant livestock pest causing severe milk and meat production losses. |
|  | *Cimex lectularius* | Bed Bug | Obligate blood-feeder and significant urban nuisance pest. |
| Model Organisms | *Drosophila melanogaster* | Vinegar Fly | Benchmark model for genetics and insect olfaction research. |
|  | *Bombyx mori* | Silk Moth | Model for Lepidopteran pheromone detection and silk production. |




# 3.2 Data Sourcing & Proteome Assembly 

Proteomes were retrieved in FASTA format from UniProt (Reference Proteomes) and NCBI RefSeq. To ensure a standardized comparative analysis, these individual proteomes were concatenated into a single master database.

Master File Name: Master_Arthropod_Reference_FINAL.fasta

Total Sequence Count: 234,439 proteins

## Protein Sequence Dataset Summary

This table summarizes the number of protein sequences retrieved for each species and the primary biological database used as the source.

## 8. Protein Sequence Sources by Species

This table summarizes the total number of protein sequences retrieved for each species and the primary biological database used.

| Species | Sequence Count | Primary Source |
|---|---|---|
| *Rhipicephalus appendiculatus* | 56,508 | NCBI RefSeq |
| *Drosophila melanogaster* | 30,802 | UniProt |
| *Stomoxys calcitrans* | 26,015 | NCBI RefSeq |
| *Cimex lectularius* | 24,194 | NCBI RefSeq |
| *Bombyx mori* | 22,510 | UniProt |
| *Aedes aegypti* | 20,643 | UniProt |
| *Anopheles gambiae* | 14,102 | UniProt |
| *Glossina morsitans* | 12,910 | UniProt |
| *Glossina fuscipes* | 13,500 | NCBI/UniProt |
| *Glossina brevipalpis* | 13,200 | NCBI/UniProt |

\*Approximate sequence counts based on combined database records.



#Assembly Workflow (Bash)

The consolidated reference database was generated using standard command-line tools to merge the final selection of 10 arthropod species, ensuring a balanced comparison between Chelicerates and Insects.

```bash
# 1. Navigating to the project workspace
cd /home/amukami/R.appendiculatus/insect_chemo_project

# 2. Consolidating the "Official 10" species-specific proteomes
# Note: We explicitly list the files to ensure Beetle/Termite are excluded
cat Rhipicephalus_appendiculatus.faa \
Glossina_morsitans.faa G_fuscipes.fasta G_brevipalpis.fasta \
Anopheles_gambiae.faa Aedes_aegypti.faa \
Cimex_lectularius.faa Stomoxys_calcitrans.faa \
Drosophila_melanogaster.faa Bombyx_mori.faa \
> Master_Arthropod_Reference_FINAL.fasta

# 3. Validation: Verifying sequence counts for the new roster
# Run this to get the exact counts for your Table above
for proteome in *.f*; do
    # Skip the Master file itself in the count
    if [[ "$proteome" != "Master"* ]]; then
        echo -n "Species File: $proteome | Count: "
        grep -c ">" "$proteome"
    fi
done

# 4. Final verification of the Master Database
echo -n "Total Master sequences: "
grep -c ">" Master_Arthropod_Referen
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
# 1. Create the results directory if it doesn't exist
mkdir -p results

# 2. Automated search loop across the Official 10-species Master Reference
for family in OR IR GR OBP CSP SNMP; do
    echo "Processing family: $family"
    
    # Using the full path to HMM profiles identified in our workspace
    hmmsearch --tblout results/${family}_hits_FINAL.tbl \
    --noali -E 1e-5 \
    insect_chemo_project/hmm_profiles/${family}.hmm \
    Master_Arthropod_Reference_FINAL.fasta
done
```

##  Comparative Genomics Results

The discovery phase revealed significant differences in the sensory evolution of ticks versus insects.

###  Summary Table: Identified Gene Counts

## 7. Comparative Chemosensory Gene Distribution

This table compares the distribution of major chemosensory gene families across selected arthropod species, including ticks, mosquitoes, flies, and model organisms.

| Species | OR | IR | GR | OBP | CSP | SNMP |
|---|---|---|---|---|---|---|
| *R. appendiculatus* (Tick) | 0 | 84 | 46 | 18 | 0 | 6 |
| *G. morsitans* (Savannah Tsetse) | 46 | 20 | 15 | 5 | 5 | 15 |
| *G. fuscipes* (Riverine Tsetse) | 42 | 28 | 15 | 6 | 5 | 15 |
| *G. brevipalpis* (Forest Tsetse) | 41 | 24 | 13 | 6 | 4 | 13 |
| *C. lectularius* (Bedbug) | 65 | 65 | 28 | 12 | 19 | 29 |
| *A. aegypti* (Mosquito) | 85 | 53 | 54 | 6 | 48 | 25 |
| *A. gambiae* (Mosquito) | 80 | 31 | 92 | 7 | 7 | 17 |
| *S. calcitrans* (Stable Fly) | 75 | 51 | 69 | 6 | 13 | 34 |
| *D. melanogaster* (Model Organism) | 64 | 53 | 84 | 9 | 9 | 35 |
| *B. mori* (Silkworm) | 109 | 43 | 26 | 5 | 43 | 35 |


# Data Extraction and Curated Databases
To facilitate downstream analysis, full protein sequences were extracted from the master file using the hit IDs from the HMMER results. This process creates specific FASTA files for each gene family, which are essential for building phylogenetic trees or performing BLAST searches.

#The Automated HMMER Search

This code runs the search for OR, IR, GR, OBP, CSP, and SNMP simultaneously. It uses the curated .hmm profiles against the 10-species consolidated proteome.

```bash
# Define the families we are looking for
for family in OR IR GR OBP CSP SNMP; do
    echo "Starting search for: $family"
    
    # Run hmmsearch with the 1e-5 threshold
    # --tblout saves the results in an easy-to-read table format
    # We point to the specific path identified in the workspace
    hmmsearch --tblout results/${family}_hits_FINAL.tbl \
    --noali -E 1e-5 \
    insect_chemo_project/hmm_profiles/${family}.hmm \
    Master_Arthropod_Reference_FINAL.fasta
done
```
#Counting the Hits (To Build the Comparative Table)

After the search, we use this command to count how many high-confidence sequences were found for each family.

```bash
for file in results/*_FINAL.tbl; do
    printf "$file: "
    grep -v "^#" "$file" | wc -l
done
```
#Extracting the FASTA Sequences

Once the hit IDs were identified in the .tbl files, we pulled the actual protein sequences out of the master file to create curated family-specific databases.

```bash
# Ensure the output directory exists
mkdir -p family_databases

for family in OR IR GR OBP CSP SNMP; do
    echo "Extracting sequences for $family..."
    
    # 1. Get the IDs from the first column of the HMMER table
    awk '!/^#/ {print $1}' results/${family}_hits_FINAL.tbl > temp_ids.txt
    
    # 2. Use the ID list to pull sequences from the Master Reference
    # We use -A 1 to grab the header and the sequence line immediately following it
    grep -Ff temp_ids.txt -A 1 Master_Arthropod_Reference_FINAL.fasta | \
    grep -v "^--" > family_databases/${family}_arthropod.fasta
done

# Clean up temporary file
rm temp_ids.txt
```
# Summary of Curated Databases
After executing the extraction pipeline above, the following family-specific FASTA databases were generated.
##  Sequence Files Summary

This table summarizes the **protein sequence files used for comparative genomics**, including the number of sequences and a brief description of each gene family.

## Chemosensory Sequence Files

The following FASTA files contain curated protein sequences used for comparative analysis of major arthropod chemosensory gene families.

## Chemosensory Gene Family Sequence Files

| File Name | Sequence Count | Description |
|---|---|---|
| OR_arthropod.fasta | 607 | Odorant Receptors: Insect-specific proteins; key for volatile sensing. |
| IR_arthropod.fasta | 452 | Ionotropic Receptors: Evolutionarily ancestral; found in all 10 species. |
| GR_arthropod.fasta | 442 | Gustatory Receptors: Taste (sugar/bitter) and CO₂ sensors. |
| OBP_arthropod.fasta | 80 | Odorant Binding Proteins: Small carriers for hydrophobic odors. |
| CSP_arthropod.fasta | 153 | Chemosensory Proteins: Broadly expressed soluble carriers. |
| SNMP_arthropod.fasta | 224 | Sensory Neuron Membrane Proteins: Essential pheromone co-factors. |

# Assembly Quality Control
The completeness of the R. appendiculatus transcriptome assembly was evaluated using BUSCO v6.0.0 against the Arachnida lineage.

#Execution Command
```bash
busco -i Assembly_Discovery/Trinity.fasta \
      -m transcriptome \
      -l arachnida_odb10 \
      -o busco_trinity_report \
      --cpu 12
```
# Quality Metrics
The analysis searched 2,934 orthologous groups specific to the Arachnida lineage (arachnida_odb10).

| Metric                | Percentage | Count |
|----------------------|-----------|-------|
| Complete (C)         | 76.7%     | 2250  |
| ├── Single-copy (S)  | 9.1%      | 268   |
| └── Duplicated (D)   | 67.6%     | 1982  |
| Fragmented (F)       | 4.4%      | 129   |
| Missing (M)          | 18.9%     | 555   |

# Protein Prediction (Open Reading Frames)
To identify potential protein-coding regions within the Trinity transcripts, we used TransDecoder v5.5.0. This process ensures that we are searching actual translated peptides rather than raw nucleotide sequences.

#Identification of Longest ORFs:

Extracted all possible Open Reading Frames (ORFs) with a minimum length of 100 amino acids.

```bash
TransDecoder.LongestOrfs -t Trinity.fasta
```
#Candidate Prediction
Predicted the most likely coding sequences based on coding score and homology.

```bash
TransDecoder.Predict -t Trinity.fasta
```
The slurm script 
```bash
#!/bin/bash
#SBATCH --job-name=busco_trans
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/home/amukami/nfs/amukami/R.appendiculatus/trans_%j.log
#SBATCH --error=/home/amukami/nfs/amukami/R.appendiculatus/trans_%j.err

# 1. The NEW path to the environment binaries on NFS
BIN="/home/amukami/nfs/amukami/tick_tools_env/bin"

# 2. Navigate to your data
cd /home/amukami/nfs/amukami/R.appendiculatus

# 3. Define the input
TRINITY_FILE="Assembly_Discovery/Trinity.fasta"

echo "--- Starting Pipeline with NFS-mounted tools ---"
date

# 4. Run BUSCO using the NFS binary
echo "Running BUSCO..."
$BIN/busco -i $TRINITY_FILE \
           -l arachnida_odb10 \
           -o busco_out \
           -m transcriptome \
           --cpu $SLURM_NTASKS \
           --force

# 5. Run TransDecoder using the NFS binaries
echo "Running TransDecoder..."
$BIN/TransDecoder.LongOrfs -t $TRINITY_FILE
$BIN/TransDecoder.Predict -t $TRINITY_FILE

date
echo "--- Pipeline Finished ---"
```

# Functional Annotation (Discovery Phase)
After generating the proteome, we performed a homology-based search to identify chemosensory and binding proteins.

1. Homology Search (BLASTp)

We used BLASTp (Protein-Protein BLAST) to compare our predicted proteome against a curated database of known arthropod chemoreceptors.

Database: A combined collection of sequences from Ixodes ricinus, Drosophila melanogaster, and other Chelicerata.

Threshold: E-value cutoff of $1 \times 10^{-5}$ to ensure high-confidence matches.

```bash
# 1. Define the targets (including the binding proteins)
FAMILIES=("OR" "GR" "IR" "OBP" "SNMP" "CSP")

# 2. Run the Master Search Loop
echo "--- Starting Master Chemosensory Search ---"

for FAM in "${FAMILIES[@]}"; do
     echo "Searching for ${FAM}s..."
     
     # Reciprocal BLASTp: searching known arthropod genes against your tick proteome
     blastp -query family_databases/${FAM}_arthropod.fasta \
            -db R_append_prot_db \
            -out ${FAM}_tick_matches.tab \
            -evalue 1e-10 \
            -outfmt "6 qseqid sseqid pident length evalue bitscore" \
            -num_threads 4
done

# 3. Generate the Final Results Summary
printf "\n===============================================\n"
printf "   R. appendiculatus DISCOVERY SCOREBOARD\n"
printf "===============================================\n"

for FAM in "${FAMILIES[@]}"; do
    # Count unique tick protein IDs (column 2 of the BLAST output)
    COUNT=$(cut -f2 ${FAM}_tick_matches.tab | sort | uniq | wc -l)
    printf "%-10s found: %d unique proteins\n" "${FAM}" "${COUNT}"
done
printf "===============================================\n"
```
| Protein Family                          | Unique Proteins Found |
| --------------------------------------- | --------------------- |
| IR (Ionotropic Receptors)               | 53                    |
| OBP (Odorant Binding Proteins)          | 28                    |
| SNMP (Sensory Neuron Membrane Proteins) | 20                    |
| GR (Gustatory Receptors)                | 15                    |
| OR (Odorant Receptors)                  | 0                     |
| CSP (Chemosensory Proteins)             | 0                     |



Sequence Extraction & File Generation
After identifying the candidate proteins via BLASTp, we extracted the full-length amino acid sequences for each chemosensory family. This ensured that our final datasets contained only the relevant predicted proteins for R. appendiculatus.

Extraction Workflow
We used a custom awk script within a bash loop to cross-reference the BLAST hits with the full proteome (Trinity.fasta.transdecoder.pep).

```bash
# 1. Create a directory for final output
mkdir -p final_chemo_sequences

# 2. Extract sequences for each family with identified hits
for FAM in IR GR OBP SNMP; do
     echo "Extracting ${FAM} sequences..."

     # Step A: Get unique protein IDs from the BLAST results
     cut -f2 ${FAM}_tick_matches.tab | sort | uniq > ${FAM}_ids.txt

     # Step B: Pull full sequences from the proteome using the ID list
     awk 'BEGIN{while((getline<"'${FAM}_ids.txt'")>0)l[">"$1]=1} /^>/{f=l[$1]} f' \
     Trinity.fasta.transdecoder.pep > final_chemo_sequences/R_append_${FAM}_candidates.faa

     echo "Done. Saved to final_chemo_sequences/R_append_${FAM}_candidates.faa"
done
```




