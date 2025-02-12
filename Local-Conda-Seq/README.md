# Local-Conda-Seq Pipeline

This folder contains scripts and workflows for running a DNA-Seq analysis pipeline using a local conda environment.

The pipeline can be executed using Snakemake.

## Pipeline Steps

1. **Quality Control**: Perform quality control on raw FASTQ files using FastQC and MultiQC.
2. **Read Alignment**: Align reads to the reference genome using BWA and sort the alignments using Samtools.
3. **Mark Duplicates**: Mark duplicate reads using GATK.
4. **Base Quality Score Recalibration**: Perform base quality score recalibration using GATK.
5. **Variant Calling**: Call variants using GATK Mutect2.
6. **Annotate Variants**: Annotate variants using snpEff.

*Optional: Index Reference Files: using GATK, Samtools and BWA*

## Quick Start

### Conda Environment

Create the conda environment using the provided `conda-env.yaml` file:

```sh
conda env create -f conda-env.yaml
conda activate dna-seq-pipeline
```

**See `conda-env.yaml` file for more details on module versions**

### Snakemake Pipeline

To run the DNA-Seq pipeline using Snakemake, use the following command:

```sh
snakemake 
```

If reference file has not been indexed previously, please run:

```sh
snakemake index_reference
```

## Configuration

The following configuration parameters should be adjusted before running the pipeline:

Pass in all of the above in `config.yaml` file

- *Read Details:* Location of raw fastq.gz paired reads
- *Parameters*
    - threads
    - ram_usage
- *Reference Files*
    - reference (.fa file)
    - known_variants (.vcf file)

## TP53 Analysis

The `TP53_Analysis.Rmd` file contains R code for analyzing TP53 mutations from the annotated VCF file and `snpEff_genes.txt` file

- Load `TP53_Analysis.Rmd` in the folder the snakemake was run.
- Filter for relevant mutations in the TP53 gene.