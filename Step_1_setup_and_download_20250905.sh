#!/bin/bash

#SBATCH --job-name=01_setup_and_download
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=240gb
#SBATCH --time=48:00:00
#SBATCH --output=./SLURM_OUT/%x_%A.out
#SBATCH --error=./SLURM_OUT/%x_%A.err

# This script downloads the raw sequencing data from the SRA database
# and prepares the FASTQ files for the analysis.

# Create a directory for the raw data if it doesn't exist
mkdir -p rawData

# List of SRA run accessions to download
sra_list=(
    "SRR7179504" "SRR7179505" "SRR7179506" "SRR7179507"
    "SRR7179508" "SRR7179509" "SRR7179510" "SRR7179511"
    "SRR7179520" "SRR7179521" "SRR7179522" "SRR7179523"
    "SRR7179524" "SRR7179525" "SRR7179526" "SRR7179527"
    "SRR7179536" "SRR7179537" "SRR7179540" "SRR7179541"
)

# Download the SRA files and convert them to gzipped FASTQ files

for sra_id in "${sra_list[@]}"; do
    echo "Downloading and processing: $sra_id"
    fastq-dump --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip --outdir rawData/ "$sra_id"
done

# Concatenate the FASTQ files for the LNCaP samples
echo "Concatenating LNCaP FASTQ files..."
cat rawData/SRR7179504_pass.fastq.gz rawData/SRR7179505_pass.fastq.gz rawData/SRR7179506_pass.fastq.gz rawData/SRR7179507_pass.fastq.gz > rawData/LNCAP_Normoxia_S1.fastq.gz
cat rawData/SRR7179508_pass.fastq.gz rawData/SRR7179509_pass.fastq.gz rawData/SRR7179510_pass.fastq.gz rawData/SRR7179511_pass.fastq.gz > rawData/LNCAP_Normoxia_S2.fastq.gz
cat rawData/SRR7179520_pass.fastq.gz rawData/SRR7179521_pass.fastq.gz rawData/SRR7179522_pass.fastq.gz rawData/SRR7179523_pass.fastq.gz > rawData/LNCAP_Hypoxia_S1.fastq.gz
cat rawData/SRR7179524_pass.fastq.gz rawData/SRR7179525_pass.fastq.gz rawData/SRR7179526_pass.fastq.gz rawData/SRR7179527_pass.fastq.gz > rawData/LNCAP_Hypoxia_S2.fastq.gz

# Rename the PC3 FASTQ files
echo "Renaming PC3 FASTQ files..."
mv rawData/SRR7179536_pass.fastq.gz rawData/PC3_Normoxia_S1.fastq.gz
mv rawData/SRR7179537_pass.fastq.gz rawData/PC3_Normoxia_S2.fastq.gz
mv rawData/SRR7179540_pass.fastq.gz rawData/PC3_Hypoxia_S1.fastq.gz
mv rawData/SRR7179541_pass.fastq.gz rawData/PC3_Hypoxia_S2.fastq.gz

# Clean up the individual SRA run files
echo "Cleaning up intermediate files..."
rm rawData/SRR*pass.fastq.gz

echo "Data download and preparation complete."