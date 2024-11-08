#!/bin/bash

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
