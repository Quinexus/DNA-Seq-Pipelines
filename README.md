# Omics Pipelines

Omics pipelines for analysing raw sequencing files

## 1) DNA-Seq

- Inputs:
    - 2 paired reads (with extension .fastq.gz)
    - Number of threads
    - Allocated ram
- Hardcoded references/inputs
    - VarScan jar file
    - Bowtie Index
    - Known-variants
    - HumanDB (for annotation)
- Outputs
    - mutliqc analysis
    - recalibrated bam
    - multianno txt file
    - csv of filtered variants

Pipelines written in bash, nextflow and snakemake

- Steps
    1. Quality Control - modules needed fastqc, anaconda3
    2. Read Alignment - modules needed bowtie2, samtools
    3. Marked Duplicates - modules needed gatk
    4. Base Score Quality Recalibration - modules needed gatk
    5. Variant Calling - modules needed java samtools + VarScan
    6. Annotated Variants - modules needed annovar
    7. Filter Variants - modules needed anaconda3 (any conda environment with panda)

## 2) RNA-Seq

- Under Construction
- Untested
- Steps
    1. Quality Control
    2. Adapter Trimming
    3. STAR Alignment
    4. Read Counting
    5. R Analysis

## 3) Alternative DNA-Seq

- Under Construction
- Similar steps as DNA-Seq but using different modules
    - Variant Calling using gatk HaplotypeCaller
    - Variant Filtering using gatk VariantFiltration
    - Variant Annotation using Ensemble vep
