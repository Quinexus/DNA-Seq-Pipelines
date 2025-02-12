# Alternative DNA-Seq Pipeline

This folder contains scripts and workflows for running an alternative DNA-Seq analysis pipeline. 

## Pipeline Steps


*Step -1 - **Index reference genome:** Checks if genome has been indexed and if not will index using BWA, samtools and GATK*
1. **Quality Control**: Perform quality control on raw FASTQ files using FastQC.
2. **Read Alignment**: Align reads to the reference genome using BWA and sort the alignments using Samtools.
3. **Mark Duplicates**: Mark duplicate reads using GATK.
4. **Base Quality Score Recalibration**: Perform base quality score recalibration using GATK.
5. **Variant Calling**: Call variants using GATK Mutect2.
6. **Filter Variants**: Filter variants using GATK FilterMutectCalls.
7. **Annotate Variants**: Annotate variants using Ensembl VEP.

## Quick Start

### Bash Pipeline

To run the DNA-Seq pipeline using the Bash script, use the following command:

```sh
bash dna_seq.sh -r <FASTQ_R1> -s <FASTQ_R2> -v <THREADS> -w <RAM_USAGE>
```

## Configuration

The following configuration parameters should be adjusted before running the above pipelines:

- *Read Details:* Location of raw fastq.gz paired reads
- *Parameters*
    - threads
    - ram_usage
- *Hardcoded References*
    - reference (.fa file)
    - known_variants (.vcf file)
- *Tools*
    - ensembl-vep

- Pass in Read Details and Parameters as arguments
- Hardcoded References and tools hardcoded into shell script

## Software and Versions

- **FastQC**: 0.11.9
- **BWA**: 0.7.17
- **Samtools**: 1.10
- **GATK**: 4.5.0.0
- **Ensembl VEP**: 104


## Job Submission

To submit the job on a cluster, use the following command:

```sh
qsub dna_seq_job.sh
```

## TP53 Analysis

The `TP53_Analysis.Rmd` file contains R code for analyzing TP53 mutations from the annotated VCF file.

- Load the annotated VCF file locally
- Filter for relevant mutations in the TP53 gene