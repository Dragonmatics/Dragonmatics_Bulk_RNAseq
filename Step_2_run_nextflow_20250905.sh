#!/bin/bash

#SBATCH --job-name=02_run_nextflow
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=256gb
#SBATCH --time=24:00:00
#SBATCH --output=./SLURM_OUT/%x_%A.out
#SBATCH --error=./SLURM_OUT/%x_%A.err

# Print job info
echo "Starting job at: $(date)"
echo "Running on node: $(hostname)"
echo "Working directory: $(pwd)"
echo "Using ${SLURM_CPUS_ON_NODE} CPU cores"

# This script runs the nf-core/rnaseq pipeline using Nextflow.
# It assumes that the raw FASTQ files are located in the 'rawData' directory
# and that the reference genome files are in the 'references' directory.

# Create a samplesheet for the Nextflow pipeline
cat << EOF > rawData/sample_sheet.csv
sample,fastq_1,fastq_2,strandedness
LNCAP_Hypoxia_S1,`pwd`/rawData/LNCAP_Hypoxia_S1.fastq.gz,,auto
LNCAP_Hypoxia_S2,`pwd`/rawData/LNCAP_Hypoxia_S2.fastq.gz,,auto
LNCAP_Normoxia_S1,`pwd`/rawData/LNCAP_Normoxia_S1.fastq.gz,,auto
LNCAP_Normoxia_S2,`pwd`/rawData/LNCAP_Normoxia_S2.fastq.gz,,auto
PC3_Hypoxia_S1,`pwd`/rawData/PC3_Hypoxia_S1.fastq.gz,,auto
PC3_Hypoxia_S2,`pwd`/rawData/PC3_Hypoxia_S2.fastq.gz,,auto
PC3_Normoxia_S1,`pwd`/rawData/PC3_Normoxia_S1.fastq.gz,,auto
PC3_Normoxia_S2,`pwd`/rawData/PC3_Normoxia_S2.fastq.gz,,auto
EOF

# Get latest Ensembl release number
latest_release=$(curl -s 'http://rest.ensembl.org/info/software?content-type=application/json' | \
                grep -o '"release":[0-9]*' | cut -d: -f2)

# Download genome sequence
wget -L ftp://ftp.ensembl.org/pub/release-${latest_release}/fasta/homo_sapiens/dna/\
            Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz

# Download genome annotation
wget -L ftp://ftp.ensembl.org/pub/release-${latest_release}/gtf/homo_sapiens/\
            Homo_sapiens.GRCh38.${latest_release}.gtf.gz


spack load openjdk@21.0.0_35

# Run the Nextflow pipeline
nextflow run nf-core/rnaseq -r 3.16.0 \
-profile apptainer \
--input rawData/sample_sheet.csv \
--outdir results \
--fasta references/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--gtf references/Homo_sapiens.GRCh38.114.gtf \
--aligner star_salmon