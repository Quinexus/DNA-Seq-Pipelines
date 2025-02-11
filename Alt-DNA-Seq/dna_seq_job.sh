#!/bin/bash
#$ -cwd           # Set the working directory for the job to the current directory
#$ -j y
#$ -pe smp 5     
#$ -l h_rt=8:0:0 
#$ -l h_vmem=8G

bash dna_seq.sh -r /data/home/ha20830/Practice/FASTQ_Raw/tumour_R1.fq.gz \
    -s /data/home/ha20830/Practice/FASTQ_Raw/tumour_R2.fq.gz \
    -v 5 \
    -w 10