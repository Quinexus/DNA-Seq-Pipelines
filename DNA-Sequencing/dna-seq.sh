#!/bin/bash

# Sequencing pipeline
# 1. QC FASTQ files
# 2. Align to genome
# 3. Mark duplicates

read1=$1
read1name=$"${read1%.*}"
read2=$2
read2name=$"${read2%.*}"
reference=$3
cores=$4
ramAllocation=$5
known_sites=$6

# Create folder for DNA Sequencing
if [ ! -d DNA_seq ]; then
    mkdir DNA_seq
fi
cd DNA_seq

# Create folder to store all logs
if [ ! -d logs ]; then
    mkdir logs
fi

# 1. QC FASTA files

# Create folder for QC files
if [ ! -d QC ]; then
    mkdir QC
fi

# Generating FASTQC analysis
if [ ! -f QC/multiqc_report.html ]; then
    module load fastqc
    module load anaconda3

    conda create --namae multiqc -c bioconda multiqc -y
    conda activate multiqc

    {
        fastqc -o QC/ read1
        fastqc -o QC/ read2
        multiqc QC/ -o QC
    } &> logs/fastqc_logs.txt

    conda deactivate multiqc
fi

# 2. Align to genome

# Create Bowtie2Idk
module load Bowtie2
if [ ! -d Bowtie2Idx ]; then
    mkdir Bowtie2Idx
fi
bowtie2-build reference Bowtie2IDx/reference.108

# Create folder for alignment
if [ ! -d Alignment ]; then
    mkdir Alignment
fi

if [ ! -f Aligment/${read1name}.bam ]; then
    {
    bowtie2 -p cores \
        --rg ID:$read1name \
        --rg SM:$read1name \
        --rg PL:ILLUMINA \
        --rg LB:$read1name \
        -x Bowtie2IDx/reference.108 \
        -1 read1 \
        -2 read2 |
        samtools sort -o Aligment/${read1name}.bam -
    } &> logs/read1name_alignment_logs.txt
fi

# 3. Mark Duplicates
if [ ! -f Alignment/${read1name}.marked.bam ]; then
    {
        gatk --java-options "-Xmx10G" MarkDuplicates \
                -I Alignment/${read1name}.bam \
                -M QC/${read1name}.marked \
                -O Alignment/${read1name}.marked.bam
    } &> logs/read1name_marked_logs.txt
fi

# 4. Run Base Score Quality Recalibration
if [ ! -f Alignment/${read1name}.recalib.bam ]; then 
        {
        gatk --java-options "-Xmx10G" BaseRecalibrator \
                -I Alignment/${read1name}.marked.bam \
                -R  $reference \
                --known-sites ${known_sites} \
                -O Alignment/${read1name}.table

        gatk --java-options "-Xmx10G" ApplyBQSR \
                -R $reference \
                -I Alignment/${read1name}.marked.bam \
                --bqsr-recal-file Alignment/${read1name}.table \
                -O Alignment/${read1name}.recalib.bam
        } &> logs/${read1name}_BSQR_logs.txt
fi

# Running Flagstat analysis
if [ ! -f logs/${read1name}_alignment_flagstat.txt ]; then 
        samtools flagstat Alignment/${read1name}.recalib.bam > logs/${read1name}_alignment_flagstat.txt
fi