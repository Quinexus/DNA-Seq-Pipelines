#!/bin/bash

# Exit on error
set -e

# Organise folders
check_mkdir() { [ ! -d $1 ] && mkdir $1; }
check_mkdir Alignment
cd Alignment

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
BAM_FILE=${PREFIX}.bam
MARKED_FILE=QC/${PREFIX}.marked
MARKED_BAM=Alignment/${PREFIX}.marked.bam
RECALIBRATED_BAM=Alignment/${PREFIX}.recalib.bam