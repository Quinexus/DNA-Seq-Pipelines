#!/bin/bash

# Exit on error
set -e

# Hardcoded reference paths
REFERENCE=/data/home/ha20830/Reference/Homo_sapiens.GRCh38.108.dna.chromosome.17.fa
KNOWN_VARIANTS=/data/teaching/bci_teaching/DNAseq/Reference/gatkResources/resources_broad_hg38_v0_1000G_omni2.5.hg38.noCHR.vcf

# Default message
usage() {
  echo $0 "-r <FASTQ_R1> -s <FASTQ_R2> -v <THREADS> -w <RAM_USAGE>"
  exit 1
}

# Parse Arguments
while [[ $# -gt 0 ]]; do
  case $1 in
  -r | --read1)
    FASTQ_R1=$2
    shift 2
    ;;
  -s | --read2)
    FASTQ_R2=$2
    shift 2
    ;;
  -v | --threads)
    THREADS=$2
    shift 2
    ;;
  -w | --ram_usage)
    RAM=$2
    shift 2
    ;;
  *)
    usage
    ;;
  esac
done

# Ensure required arguments are provided
if [[ -z $FASTQ_R1 || -z $FASTQ_R2 || -z $THREADS || -z $RAM ]]; then
  usage
fi

# Derive file names from the input files
BASENAME_R1=$(basename "$FASTQ_R1" .fastq.gz)
PREFIX=${BASENAME_R1%%_*}

# Output files
BAM_FILE=Alignment/${PREFIX}.bam
MARKED_FILE=QC/${PREFIX}.marked
MARKED_BAM=Alignment/${PREFIX}.marked.bam
RECALIB_BAM=Alignment/${PREFIX}.recalib.bam
RAW_VCF=VCF/${PREFIX}.raw.vcf.gz
FILTERED_VCF=VCF/${PREFIX}.filtered.vcf.gz
ANNOTATED_VCF=VCF/${PREFIX}.annotated.vcf
RECALIB_TABLE=Alignment/${PREFIX}.recalib.table

# Create necessary directories
mkdir -p DNA-Seq-${PREFIX}
cd DNA-Seq-${PREFIX}
mkdir -p QC Alignment VCF

# Step -1: Index reference genome
module load bwa
bwa index $REFERENCE
module unload bwa

# Step 0: Run Quality Control
module load fastqc
fastqc $FASTQ_R1 $FASTQ_R2 -o QC/
module unload fastqc

# Step 1: Align reads
if [ ! -f $BAM_FILE ]; then
  module load bwa samtools

  bwa mem \
    $REFERENCE \
    $FASTQ_R1 \
    $FASTQ_R2 |
    samtools sort -o $BAM_FILE -
  
  samtools index $BAM_FILE

  module unload bwa samtools
fi

# Step 2: Mark Duplicates
if [ ! -f $MARKED_BAM ]; then
  module load gatk samtools

  gatk --java-options "-Xmx${RAM}G" MarkDuplicates \
    -I $BAM_FILE \
    -M $MARKED_FILE \
    -O $MARKED_BAM
  
  samtools index $MARKED_BAM

  module unload gatk samtools
fi

# Step 3: Base Recalibration and Base Quality Score Recalibration
if [ ! -f $RECALIB_BAM ]; then
  module load gatk samtools

  gatk --java-options "-Xmx${RAM}G" BaseRecalibrator \
    -I $MARKED_BAM \
    -R $REFERENCE \
    --known-sites $KNOWN_VARIANTS \
    -O $RECALIB_TABLE

  gatk --java-options "-Xmx${RAM}G" ApplyBQSR \
    -R $REFERENCE \
    -I $MARKED_BAM \
    --bqsr-recal-file $RECALIB_TABLE \
    -O $RECALIB_BAM

  samtools index $RECALIB_BAM

  module unload gatk samtools
fi

# Step 4: Calling Variants
if [ ! -f $RAW_VCF ]; then
  module load gatk

  gatk --java-options "-Xmx${RAM}G" Mutect2 \
    -R $REFERENCE \
    -I $RECALIB_BAM \
    -O $RAW_VCF

  module unload gatk
fi

# Step 5: Filter variants
if [ ! -f $FILTERED_VCF ]; then
  module load gatk

  gatk --java-options "-Xmx${RAM}G" FilterMutectCalls \
    -R $REFERENCE \
    -V $RAW_VCF \
    -O $FILTERED_VCF

  module unload gatk
fi

# Step 6: Annotate Variants
if [ ! -f $ANNOTATED_VCF ]; then
  module load ensembl-vep

  vep --input_file $FILTERED_VCF \
    --output_file $ANNOTATED_VCF \
    --format vcf \
    --species homo_sapiens \
    --database \
    --fork 20 \
    --everything
  
  module unload ensembl-vep
fi