#!/bin/bash

# Exit on error
set -e

# Hardcoded reference paths 
cp -vR /data/teaching/bci_teaching/DNAseq/VarScan.v2.4.3.jar ./
BOWTIE_INDEX=/data/teaching/bci_teaching/DNAseq/Reference/Bowtie2Idx/GRCh38.108.chr17
REFERENCE=/data/teaching/bci_teaching/DNAseq/Reference/Homo_sapiens.GRCh38.108.dna.chromosome.17.fa
KNOWN_VARIANTS=/data/teaching/bci_teaching/DNAseq/Reference/gatkResources/resources_broad_hg38_v0_1000G_omni2.5.hg38.noCHR.vcf
HUMANDB_PATH="/data/teaching/bci_teaching/DNAseq/Reference/humandb/"
VARSCAN="VarScan.v2.4.3.jar"

# Default message
usage() {
    echo $0 "-r <FASTQ_R1> -s <FASTQ_R2> -v <THREADS> -w <RAM_USAGE>"
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
        -v|--threads)
            THREADS=$2
            shift 2
            ;;
        -w|--ram_usage)
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
BASENAME_R2=$(basename "$FASTQ_R2" .fastq.gz)
PREFIX=${BASENAME_R1%%_*}

# Output files
INDEX=bowtie-idx/${PREFIX}.reference
BAM_FILE=Alignment/${PREFIX}.bam
MARKED_FILE=QC/${PREFIX}.marked
MARKED_BAM=Alignment/${PREFIX}.marked.bam
TABLE_FILE=Alignment/${PREFIX}.table
RECALIB_BAM=Alignment/${PREFIX}.recalib.bam
RAW_VCF=VCF/${PREFIX}.raw.vcf.gz
FILTERED_VCF=VCF/${PREFIX}.filtered.vcf.gz
ANNOTATED_VCF=VCF/${PREFIX}.annotated.vcf

mkdir -p DNA-Seq-${PREFIX}
cd DNA-Seq-${PREFIX}

# Step -2: Index reference genome
# samtools faidx $REFERENCE
#gatk CreateSequenceDictionary -R $REFERENCE -O ${REFERENCE%.fa}.dict
# gatk IndexFeatureFile -I $KNOWN_VARIANTS

# Step -1: Build Bowtie2 index
# check_mkdir bowtie-idx
# bowtie2-build $REFERENCE $INDEX

# Step 0: Run Quality Control 
# In CONDA environment
mkdir -p QC
fastqc $READ_1 $READ_2 -o QC/
multiqc QC/

# Step 1: Align reads 
mkdir -p Alignment

if [ ! -f $BAM_FILE ]; then
    bowtie2 -p $THREADS \
        --rg ID:$PREFIX \
        --rg SM:$PREFIX \
        --rg PL:ILLUMINA \
        --rg LB:$PREFIX \
        -x $BOWTIE_INDEX \
        -1 $FASTQ_R1 \
        -2 $FASTQ_R2 |
        samtools sort -o $BAM_FILE -
fi

# Step 2: Mark Duplicates
if [ ! -f $MARKED_BAM ]; then
    gatk --java-options "-Xmx${RAM}G" MarkDuplicates \
            -I $BAM_FILE \
            -M $MARKED_FILE \
            -O $MARKED_BAM
fi

# Step 3: Base Recalibration and Base Quality Score Recalibration
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

# Step 4: Calling Variants
samtools index Alignment/tumour.recalib.bam
gatk --java-options "-Xmx10G" HaplotypeCaller \
    -R $REFERENCE \
    -I $RECALIB_BAM \
    -O $RAW_VCF \
    --dont-use-soft-clipped-bases \
    --standard-min-confidence-threshold-for-calling 20

# Step 5: Filter variants
gatk --java-options "-Xmx10G" VariantFiltration \
    -R $REFERENCE \
    -V $RAW_VCF \
    -O $FILTERED_VCF \
    --filter-name "QDFilter" --filter-expression "QD < 2.0" \
    --filter-name "FSFilter" --filter-expression "FS > 30.0" \
    --filter-name "MQFilter" --filter-expression "MQ < 40.0"

# Step 6: Annotate Variants
vep --input_file $FILTERED_PASS_VCF \
    --output_file $ANNOTATED_VCF \
    --format vcf \
    --species homo_sapiens \
    --database \
    --fork 20 \
    --everything
