#!/bin/bash

set -euo pipefail  # Exit on error and handle unset variables

# Define constants
PREFIX="1M_SRR9336468"
READ_1="RNA-Seq_Sample_Files/1M_SRR9336468_1.fastq.gz"
READ_2="RNA-Seq_Sample_Files/1M_SRR9336468_2.fastq.gz"
REFERENCE="RNA-Seq_Sample_Files/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa"
KNOWN_VCF="RNA-Seq_Sample_Files/saccharomyces_cerevisiae.vcf.gz"

# Create required directories
mkdir -p RNA-Seq_Sample_Files Alignment QC bowtie_base

# Download sample RNA-Seq data and reference files
echo "Downloading RNA-Seq sample files..."
wget -q -O RNA-Seq_Sample_Files.zip https://bioinfogp.cnb.csic.es/files/samples/rnaseq/RNA-Seq_Sample_Files.zip
unzip -o RNA-Seq_Sample_Files.zip -d RNA-Seq_Sample_Files
extract ${REFERENCE}.gz

echo "Downloading reference VCF..."
wget -q -P RNA-Seq_Sample_Files/ https://ftp.ensembl.org/pub/release-113/variation/vcf/saccharomyces_cerevisiae/saccharomyces_cerevisiae.vcf.gz

# Index the reference genome
echo "Indexing reference genome..."
samtools faidx $REFERENCE
gatk CreateSequenceDictionary -R $REFERENCE -O ${REFERENCE%.fa}.dict
gatk IndexFeatureFile -I $KNOWN_VCF

# Build Bowtie2 index
echo "Building Bowtie2 index..."
bowtie2-build $REFERENCE bowtie_base/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel

# Perform quality control
echo "Running FastQC..."
fastqc $READ_1 $READ_2 -o QC/
cd QC
multiqc . 
cd ..

# Align reads with Bowtie2
echo "Aligning reads with Bowtie2..."
bowtie2 -p 4 \
    --rg ID:$PREFIX \
    --rg SM:$PREFIX \
    --rg PL:ILLUMINA \
    --rg LB:$PREFIX \
    -x bowtie_base/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel \
    -1 $READ_1 \
    -2 $READ_2 | \
    samtools sort -o Alignment/${PREFIX}.bam -

# Mark duplicates
echo "Marking duplicates..."
gatk MarkDuplicates \
    -I Alignment/${PREFIX}.bam \
    -M QC/${PREFIX}.marked \
    -O Alignment/${PREFIX}.marked.bam

# Base recalibration
echo "Running BaseRecalibrator..."
gatk BaseRecalibrator \
    -I Alignment/${PREFIX}.marked.bam \
    -R $REFERENCE \
    --known-sites $KNOWN_VCF \
    -O Alignment/${PREFIX}.table

# Apply base quality score recalibration
echo "Applying BQSR..."
gatk ApplyBQSR \
    -R $REFERENCE \
    -I Alignment/${PREFIX}.marked.bam \
    --bqsr-recal-file Alignment/${PREFIX}.table \
    -O Alignment/${PREFIX}.recalib.bam

# Call variants
echo "Calling variants with GATK HaplotypeCaller..."
gatk HaplotypeCaller \
    -R $REFERENCE \
    -I Alignment/${PREFIX}.recalib.bam \
    -O Alignment/${PREFIX}.raw.vcf.gz \
    --dont-use-soft-clipped-bases \
    --standard-min-confidence-threshold-for-calling 20

# Filter raw variants
echo "Filtering variants..."
gatk VariantFiltration \
    -R $REFERENCE \
    -V Alignment/${PREFIX}.raw.vcf.gz \
    -O Alignment/${PREFIX}.filtered.vcf.gz \
    --filter-name "QDFilter" --filter-expression "QD < 2.0" \
    --filter-name "FSFilter" --filter-expression "FS > 30.0" \
    --filter-name "MQFilter" --filter-expression "MQ < 40.0"

# Annotate variants with Ensembl VEP
echo "Annotating variants with Ensembl VEP..."
vep_install --SPECIES saccharomyces_cerevisiae --CACHE --ASSEMBLY R64-1-1 --VERSION 113
vep \
    --input_file Alignment/${PREFIX}.filtered.vcf.gz \
    --output_file Alignment/${PREFIX}.annotated.vcf \
    --format vcf \
    --vcf \
    --species saccharomyces_cerevisiae \
    --cache \
    --assembly R64-1-1 \
    --offline \
    --force_overwrite

# Generate statistics
echo "Generating variant statistics..."
bcftools stats Alignment/${PREFIX}.annotated.vcf > Alignment/${PREFIX}.variant_stats.txt

# Final message
echo "Pipeline completed successfully!"
echo "Annotated VCF file: Alignment/${PREFIX}.annotated.vcf"
echo "Variant statistics: Alignment/${PREFIX}.variant_stats.txt"
