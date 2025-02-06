#!/bin/bash

set -e # exit if error

mkdir -p DNA_SEQ
cd DNA_SEQ

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

if [ ! -f QC/multiqc_report.html ]; then
    echo "Step 1: Quality Control"

    module load fastqc anaconda3
    conda create --name multiqc -c bioconda multiqc -y # Not ideal to hardcode in conda environment

    conda activate multiqc
    fastqc $FASTQ_R1 $FASTQ_R2 -o QC/ -t $THREADS
    multiqc QC/ -o QC

    conda deactivate
    module unload fastqc
fi

# Step 2: Read Alignment

## Input: FASTQ_R1, FASTQ_R2, BOWTIE_INDEX
## Output: BAM_FILE

mkdir -p Alignment

if [ ! -f $BAM_FILE ]; then
    echo "Step 2: Read Alignment"
    module load bowtie2 samtools

    bowtie2 -p $THREADS \
        --rg ID:$PREFIX \
        --rg SM:$PREFIX \
        --rg PL:ILLUMINA \
        --rg LB:$PREFIX \
        -x $BOWTIE_INDEX \
        -1 $FASTQ_R1 \
        -2 $FASTQ_R2 |
        samtools sort -o $BAM_FILE -
    
    module unload bowtie2 samtools
fi


# Step 3: Mark Duplicates

## Input: BAM_FILE
## Output: MARKED_BAM (also MARKED_FILE)

if [ ! -f $MARKED_BAM ]; then
    echo "Step 3: Mark Duplicates"

    module load gatk

        gatk --java-options "-Xmx${RAM}G" MarkDuplicates \
                -I $BAM_FILE \
                -M $MARKED_FILE \
                -O $MARKED_BAM
    
    module unload gatk
fi

# Step 4: Base Score Quality Recalibration

## Input: MARKED_BAM, REFERENCE, KNOWN_VARIANTS
## Output: $RECALIB_BAM (also $RECALIB_TABLE)

if [ ! -f $RECALIB_BAM ]; then 

    echo "Step 4: Base Score Quality Recalibration"

    module load gatk

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
    
    module unload gatk
fi


# Step 5: Variant Calling

## Input: REFERENCE, RECALIB_BAM
## Output: RAW_VCF

mkdir -p VCF
if [ ! -f $RAW_VCF ]; then
    echo "Step 5: Variant Calling"

    module load java samtools

    samtools mpileup \
        -q 20 \
        -f $REFERENCE \
        $RECALIB_BAM |
    java -jar $VARSCAN mpileup2snp \
        --min-coverage 20 \
        --min-avg-qual 20 \
        --min-read2 4 \
        --p-value 0.2 \
        --min-var-freq 0.01 \
        --strand-filter 1 \
        --output-vcf 1 > $RAW_VCF

    module unload java samtools
fi

# Step 6: Annotate Variants

## Input: RAW_VCF, HUMANDB_PATH
## Output: MULTIANNO (+ other filtered vcfs)

if [ ! -f $MULTIANNO ]; then
    echo "Step 6: Annotate Variants"

    module load annovar

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

fi

# Step 7: Filter out variants with variant allele frequency <10% using python script 
# NEEDS CONDA ENVIRONMENT!

## Input: MULTIANNO
## Output: filtered.csv

if [ ! -f "${PREFIX}_filtered.csv" ]; then
    python3 dna-seq-filter.py $MULTIANNO
fi