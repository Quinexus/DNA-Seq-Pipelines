#!/bin/bash

# Exit on error
set -e

# Default message
usage() {
    echo "$0 -r <FASTQ_R1> -s <FASTQ_R2> -t <REFERENCE>"
    exit 1
}

# Parse Arguments
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
        -t|--reference)
            REFERENCE=$2
            shift 2
            ;;
        -u|--known-sites)
            KNOWN_VARIANTS=$2
            shift 2
            ;;
        *)
            usage
            ;;
    esac
done

# Ensure required arguments are provided
if [[ -z $FASTQ_R1 || -z $FASTQ_R2 || -z $REFERENCE ]]; then
    usage
fi

# Derive file names from the input files
BASENAME_R1=$(basename "$FASTQ_R1" .fastq.gz)
BASENAME_R2=$(basename "$FASTQ_R2" .fastq.gz)
PREFIX=${BASENAME_R1%%_*}

# Output files
INDEX=bowtie-idx/${PREFIX}.reference
BAM_FILE=Alignment/${PREFIX}.bam
MARKED_FILE=QC/${PREFIX}.marked
MARKED_BAM=Alignment/${PREFIX}.marked.bam
TABLE_FILE=Alignment/${PREFIX}.table
RECALIB_BAM=Alignment/${PREFIX}.recalib.bam
RAW_VCF=Alignment/${PREFIX}.raw.vcf.gz
FILTERED_VCF=Alignment/${PREFIX}.filtered.vcf.gz

mkdir -p DNA-Seq-${PREFIX}
cd DNA-Seq-${PREFIX}

# Index reference genome
samtools faidx $REFERENCE
gatk CreateSequenceDictionary -R $REFERENCE -O ${REFERENCE%.fa}.dict
gatk IndexFeatureFile -I $KNOWN_VARIANTS

# Build Bowtie2 index
check_mkdir bowtie-idx
bowtie2-build $REFERENCE $INDEX

# Run Quality Control
mkdir -p QC
fastqc $READ_1 $READ_2 -o QC/
multiqc QC/

mkdir -p Alignment

# Align reads 
bowtie2 -p 4 \
    --rg ID:$PREFIX \
    --rg SM:$PREFIX \
    --rg PL:ILLUMINA \
    --rg LB:$PREFIX \
    -x bowtie-idx/${PREFIX}.reference \
    -1 $READ_1 \
    -2 $READ_2 | \
    samtools sort -o $BAM_FILE -

# Mark Duplicates
gatk MarkDuplicates \
    -I $BAM_FILE \
    -M $MARKED_FILE \
    -O $MARKED_BAM 

# Base Recalibration and Base Quality Score Recalibration
gatk BaseRecalibrator \
    -I $MARKED_BAM \
    -R $REFERENCE \
    --known-sites $KNOWN_VCF \
    -O $TABLE_FILE

gatk ApplyBQSR \
    -R $REFERENCE \
    -I $MARKED_BAM \
    --bqsr-recal-file $TABLE_FILE \
    -O $RECALIB_BAM

# Calling Variants
gatk HaplotypeCaller \
    -R $REFERENCE \
    -I $RECALIB_BAM \
    -O $RAW_VCF \
    --dont-use-soft-clipped-bases \
    --standard-min-confidence-threshold-for-calling 20

# Filter raw variants
gatk VariantFiltration \
    -R $REFERENCE \
    -V $RAW_VCF \
    -O $FILTERED_VCF \
    --filter-name "QDFilter" --filter-expression "QD < 2.0" \
    --filter-name "FSFilter" --filter-expression "FS > 30.0" \
    --filter-name "MQFilter" --filter-expression "MQ < 40.0"