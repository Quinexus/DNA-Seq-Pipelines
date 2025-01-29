#!/bin/bash

set -e # exit if error

mkdir -p DNA_SEQ
cd DNA_SEQ

# Hardcoded reference paths 
cp -vR /data/teaching/bci_teaching/DNAseq/VarScan.v2.4.3.jar ./
HUMANDB_PATH="/data/teaching/bci_teaching/DNAseq/Reference/humandb/"
VARSCAN=VarScan.v2.4.3.jar

# Default message
usage() {
    echo $0 -r <FASTQ_R1> -s <FASTQ_R2> -t <REFERENCE_BOWTIE> -u <REFERENCE> -v <THREADS> -x <KNOWN_VARIANTS>
    exit 1
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -r|--read1)
            FASTQ_R1=$2
            shift 2
            ;;
        -s|--read2)
            FASTQ_R2=$2
            shift 2
            ;;
        -t|--reference_bowtie)
            REFERENCE_BOWTIE=$2
            shift 2
            ;;
        -u|--reference)
            REFERENCE=$2
            shift 2
            ;;
        -u|--threads)
            THREADS=$2
            shift 2
            ;;
        -v|--ram_usage)
            RAM=$2
            shift 2
            ;;
        -w|--known_variants)
            KNOWN_VARIANTS=$2
            shift 2
            ;;
        *)
            usage
            ;;
    esac
done

# Ensure required arguments are provided
if [[ -z $FASTQ_R1 || -z $FASTQ_R2 || -z $REFERENCE || -z $THREADS || -z $RAM || -z $KNOWN_VARIANTS || -z $REFERENCE_BOWTIE ]]; then
    usage
fi

# Derive file names from the input files
BASENAME_R1=$(basename "$FASTQ_R1" .fastq.gz)
BASENAME_R2=$(basename "$FASTQ_R2" .fastq.gz)
PREFIX=${BASENAME_R1%%_*}

# Output files
BAM_FILE=Alignment/${PREFIX}.bam
MARKED_FILE=QC/${PREFIX}.marked
MARKED_BAM=Alignment/${PREFIX}.marked.bam
RECALIB_BAM=Alignment/${PREFIX}.recalib.bam
RECALIB_TABLE=Alignment/${PREFIX}.table
RAW_VCF=VCF/${PREFIX}.vcf
PASS_VCF=VCF/${PREFIX}.pass.vcf
MULTIANNO=VCF/${PREFIX}.hg38_multianno.txt

# Quality Control - FASTQC analysis
mkdir -p QC

# Step 1: Quality Control

## Input: FASTQ_R1 FASTQ_R2
## Output: multiqc_report.html

echo "Step 1: Quality Control"

module load fastqc anaconda3
conda create --name multiqc -c bioconda multiqc -y # Not ideal to hardcode in conda environment

if [ ! -f QC/multiqc_report.html ]; then
    conda activate multiqc
    fastqc $FASTQ_R1 $FASTQ_R2 -o QC/ -t $THREADS
    multiqc QC/ -o QC
    conda deactivate
fi
module unload fastqc

# Step 2: Read Alignment

## Input: FASTQ_R1, FASTQ_R2, REFERENCE_BOWTIE
## Output: BAM_FILE

echo "Step 2: Read Alignment"
mkdir -p Alignment

module load bowtie2
if [ ! -f $BAM_FILE ]; then
    bowtie2 -p $THREADS \
        --rg ID:$PREFIX \
        --rg SM:$PREFIX \
        --rg PL:ILLUMINA \
        --rg LB:$PREFIX \
        -x $REFERENCE_BOWTIE \
        -1 $FASTQ_R1 \
        -2 $FASTQ_R2 |
        samtools sort -o $BAM_FILE -
fi
module unload bowtie2

# Step 3: Mark Duplicates

## Input: BAM_FILE
## Output: MARKED_BAM (also MARKED_FILE)

echo "Step 3: Mark Duplicates"

module load gatk
if [ ! -f $MARKED_BAM ]; then
        gatk --java-options "-Xmx${RAM}G" MarkDuplicates \
                -I $BAM_FILE \
                -M $MARKED_FILE \
                -O $MARKED_BAM
fi

# Step 4: Base Score Quality Recalibration

## Input: MARKED_BAM, REFERENCE, KNOWN_VARIANTS
## Output: $RECALIB_BAM (also $RECALIB_TABLE)

echo "Step 4: Base Score Quality Recalibration"

if [ ! -f $RECALIB_BAM ]; then 
        gatk --java-options "-Xmx${RAM}G" BaseRecalibrator \
                -I $MARKED_BAM \
                -R  $REFERENCE \
                --known-sites $KNOWN_VARIANTS \
                -O $RECALIB_TABLE
        gatk --java-options "-Xmx${RAM}G" ApplyBQSR \
                -R $REFERENCE \
                -I $MARKED_BAM \
                --bqsr-recal-file $RECALIB_TABLE \
                -O $RECALIB_BAM
fi
module unload gatk

# Step 5: Variant Calling

## Input: REFERENCE, RECALIB_BAM
## Output: RAW_VCF

echo "Step 5: Variant Calling"
mkdir -p VCF

module load java
if [ ! -f $RAW_VCF ]; then
    samtools mpileup \
        -q 20 \
        -f $REFERENCE \
        $RECALIB_BAM |
    java -jar VARSCAN mpileup2snp \
        --min-coverage 20 \
        --min-avg-qual 20 \
        --min-read2 4 \
        --p-value 0.2 \
        --min-var-freq 0.01 \
        --strand-filter 1 \
        --output-vcf 1 > $RAW_VCF
fi

# Step 6: Annotate Variants
echo "Step 6: Annotate Variants"

## Input: RAW_VCF, HUMANDB_PATH
## Output: MULTIANNO (+ other filtered vcfs)

module load annovar
if [ ! -f $PASS_VCF ]; then
    convert2annovar.pl --format vcf4 \
        $RAW_VCF \
        --includeinfo \
        --filter PASS \
        --outfile $PASS_VCF

    # Filter mutations with MAF > 0.01 from the 1000 Genomes Project
    annotate_variation.pl -filter \
        -dbtype 1000g2015aug_all \
        -buildver hg38 \
        -out VCF/${PREFIX} \
        $PASS_VCF \
        $HUMANDB_PATH \
        -maf 0.01

    # Filter mutations with MAF > 0.01 from the Exome Sequencing Project
    annotate_variation.pl -filter \
        -dbtype esp6500siv2_all \
        -buildver hg38 \
        -out VCF/${PREFIX} \
        VCF/${PREFIX}.hg38_ALL.sites.2015_08_filtered \
        "$HUMANDB_PATH" \
        -score_threshold 0.01

    # Annotate the variants with RefGene, dbSNP, COSMIC, and Cytoband information
    table_annovar.pl \
        -buildver hg38 \
        -out VCF/${PREFIX} \
        VCF/${PREFIX}.hg38_esp6500siv2_all_filtered \
        $HUMANDB_PATH \
        -remove \
        -otherinfo \
        -protocol refgene,avsnp150,cosmic92_coding,cytoband \
        -operation g,f,f,r \
        -nastring .

# Step 7: Filter out variants with variant allele frequency <10% using R script 

## Input: MULTIANNO
## Output: filtered.csv

module load R
Rscript dna-seq-filter.R $MULTIANNO