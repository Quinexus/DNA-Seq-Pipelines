#!/bin/bash

mkdir -p rna-seq
cd rna-seq

# Parse Arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -r|--read1) FASTQ_R1=$2 shift 2 ;;
        -s|--read2) FASTQ_R2=$2 shift 2 ;;
        -t|--reference) REFERENCE=$2 shift 2 ;;
        -u|--gtf) GTF=$2 shift 2 ;;
        -v|--rsem) RSEM_REF=$2 shift 2 ;;
    esac
done

# Derive file names from the input files
BASENAME=$(basename "$FASTQ_R1" .fastq.gz))
PREFIX=${BASENAME%%_*}

# Output Files
TRIMMED_FASTQ_R1=trimmed/FASTQs/${PREFIX}_val_1.fq
TRIMMED_FASTQ_R2=trimmed/FASTQs/${PREFIX}_val_2.fq
STAR_BAM=${PREFIX}_Aligned.toTranscriptome.out.bam
READ_COUNTS_FOLDER=rsem_counts/

# Create folders
mkdir -p trimmed_FASTQs trimmed_fastQC rsem_counts

# Step 1: Quality Control

## Input: FASTQ_R1 FASTQ_R2
## Output: FASTQ.html

module load fastqc
fastqc FASTQ_R1
fastqc FASTQ_R2
# ? multiqc

# Step 2: Adapter Trimming

## Input: FASTQ_R1 FASTQ_R2
## Output: TRIMMED_FASTQ_R1 TRIMMED_FASTQ_R2

module load trimgalore

trim_galore --paired \
    --illumina \
    --dont_gzip \
    --fastqc_args "-o trimmed_fastQC" \
    -o trimmed_FASTQs \
    $FASTQ_R1 $FASTQ_R2

# Step 3: STAR Alignment

## Input: TRIMMED_FASTQ_R1 TRIMED_FASTQ_R2 REFERENCE GTF
## Output: STAR_BAM

module load star

STAR --runThreadN 2 \
     --genomeDir $REFERENCE \
     --sjdbGTFfile $GTF \
     --sjdbOverhang 149 \
     --outFilterType BySJout \
     --outFilterMultimapNmax 20 \
     --alignSJoverhangMin 8 \
     --alignSJDBoverhangMin 1 \
     --outFilterMismatchNmax 999 \
     --outFilterMismatchNoverReadLmax 0.04 \
     --alignIntronMin 20 \
     --alignIntronMax 1000000 \
     --alignMatesGapMax 1000000 \
     --readFilesIn $TRIMMED_FASTQ_R1 $TRIMMED_FASTQ_R2 \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode TranscriptomeSAM \
     --outFileNamePrefix ${PREFIX}_

# Step 4: Read Counting

## Input: STAR_BAM RSEM_REF
## Output: READ_COUNTS

module load rsem
rsem-calculate-expression --bam \
    --strandedness reverse \
    --no-bam-output \
    -p 2 \
    --paired-end \
    $STAR_BAM \
    $RSEM_REF \
    ${PREFIX}_Quant


# Step 5: R Analysis

# Input: READ_COUNTS
# Output: read counts folder

module load R
Rscript rna-seq.R $READ_COUNTS_FOLDER sample-groups.txt output