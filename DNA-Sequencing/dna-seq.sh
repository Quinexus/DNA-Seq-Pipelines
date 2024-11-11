#!/bin/bash

<<<<<<< HEAD
cd ~/DNA-Sequencing

bowtie2-build reference_genome/GCF_000006765.1_ASM676v1_genomic.fna reference_genome/Bowtie2Idx/ASM676v1.108

mkdir QC
fastqc -o QC/ raw_fastq/*

conda create --name multiqc -c bioconda multiqc
conda activate multiqc
multiqc QC/ -o QC
conda deactivate

mkdir Alignment

time bowtie2 -p 1 \
    -x reference_genome/Bowtie2Idx/ASM676v1.108 \
    -1 raw_fastq/ERR3697100_1.fastq.gz \
    -2 raw_fastq/ERR3697100_2.fastq.gz |
    samtools sort -o Alignment
=======
# Sequencing pipeline
# 1. QC FASTQ files
# 2. Align to genome
# 3. Mark duplicates

# Ask user whick step to process
choice=''
while [[ "$choice" != 'q' ]]
do

echo "\nWhich step do you want to start?"
echo "Press 'a' to produce QC FASTQ files"
echo "Press 'b' to align to genome"
echo "Press 'c' to mark duplicates"
echo "Press 'q' to quit"
read choice

# Step 1: QC FASTQ files
if [[ $choice == 'a' ]]; then
    # Load required modules: fastqc
    echo "Loading fastqc and conda (make sure anaconda is installed)"
    if command -v module >/dev/null 2>&1; then
        module load fastqc
        module load anaconda3
        echo "module package found"
    elif command -v fastqc >/dev/null 2>&1; then
        echo "fastqc found"
    else
        echo "fastqc not found, installing fastqc"
        sudo apt install -y fastqc
    fi

    echo "Files will be stored in directory QC"
    mkdir QC

    echo "Provide location to raw files:"
    read raw_fasta_files

    fastqc -o QC/ ${raw_fasta_files}*

    echo "Loading conda environment to use multiqc"
    conda create --name multiqc_env -c bioconda multiqc -y
    source activate multiqc_env
    multiqc QC/ -o QC
    conda deactivate

# Step 2: Align to genome
elif [[ $choice == 'b' ]]; then
    # Load required modules: bowtie2, samtools
    echo "Loading required modules"
    if command -v module >/dev/null 2>&1; then
        module load bowtie2
        module load samtools
        echo "module package found"
    elif command -v bowtie >/dev/null 2>&1 && command -v samtools >/dev/null 2>&1; then
        echo "bowtie and samtools found"
    else
        echo "one or more packages not found"
        echo "installing all packages"
        conda create --name align_genome_env -c bioconda bowtie2 samtools -y
        source activate align_genome_env
    fi

    echo "Add read 1 Raw FASTA:"
    read fasta_raw1
    echo "Add read 2 Raw FASTA:"
    read fasta_raw2
    echo "Add reference genome fasta:"
    read reference_genome

    mkdir Bowtie2Idx
    bowtie2-build ${reference_genome} Bowtie2Idx/reference.108

    mkdir Alignment
    echo "Commencing Alignment (may take some time!)"q
    time bowtie2 -p 2 \
        -x Bowtie2Idx/reference.108 \
        -1 ${fasta_raw1} \
        -2 ${fasta_raw2} |
        samtools sort -o Alignment/aligned.bam -
    echo "Aligned bam stored as aligned.bam in Alignment folder" 
    conda deactivate

# Step 3: Mark duplicates
elif [[ $choice == 'c' ]]; then
    # Load required modules: gatk
    if command -v module >/dev/null 2>&1; then
        module load java
        module load gatk
        echo "module package found"
    elif command -v gatk >/dev/null 2>&1; then
        echo "gatk found"
    else
        echo "one or more packages not found"
        echo "installing all packages"
        conda create --name mark_duplicates_env -c bioconda gatk -y
        source activate mark_duplicates_env
    fi

    mkdir Alignment
    mkdir QC

    echo "Load bam file, if empty will load from Alignment/aligned.bam"
    read aligned_bam
    if [-z $aligned_bam]; then
        aligned_bam="Alignment/aligned.bam"
    fi

    echo "Marking Duplicates"
    gatk --java-options "-Xmx1G" MarkDuplicates \
        -I {$aligned_bam} \
        -M QC/markedQC.marked \
        -O Alignment/aligned.marked.bam

    echo "Marked duplicates, markedQC.marked stored in QC folder, aligned.marked.bam stored in BAM folder"

# Else
else
    if [[ $choice != 'q' ]]; then
        echo "Invalid command"
    fi
fi

done

echo "Finished"
>>>>>>> master
