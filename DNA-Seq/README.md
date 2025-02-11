# DNA-Seq Pipeline

This folder contains scripts and workflows for running a DNA-Seq analysis pipeline. The pipeline can be executed using different workflow managers including Bash, Snakemake, and Nextflow.

## Quick Start

### Bash Pipeline

To run the DNA-Seq pipeline using the Bash script, use the following command:

```sh
bash dna-seq.sh -r <FASTQ_R1> -s <FASTQ_R2> -v <THREADS> -w <RAM_USAGE>
```

### Snakemake Pipeline

To run the DNA-Seq pipeline using Snakemake, use the following command in the snakemake folder.
Please make sure to activate a conda pipeline (needed for python functionality)

```sh
snakemake
```

### Config File Adjustments

Before running the Snakemake pipeline, you need to adjust the configuration file to match your input data and parameters. The configuration file is located at `/DNA-Seq/Snakemake/config.yaml`. Here are the key sections to modify:

- **samples**: Specify the paths to your FASTQ files.
- **threads**: Set the number of threads to use.
- **ram**: Set the amount of RAM to use.
- **reference_genome**: Provide the path to the reference genome.

Make sure to adjust these paths and parameters according to your specific setup before running the pipeline.



### Nextflow Pipeline

To run the DNA-Seq pipeline using Nextflow, use the following command:

```sh
nextflow run dna-seq.nf
```

Similarly adjust config params by viewing dna-seq.nf

## Software and Versions



- **FastQC**: 0.11.9
- **MultiQC**: 1.9
- **Bowtie2**: 2.4.2
- **Samtools**: 1.10
- **GATK**: 4.5.0.0
- **VarScan**: 2.4.3
- **Annovar**: 2019Oct24
- **Python**: 3.8.5 (for filtering script)
- **R**: 4.0.3 (for filtering script)

### Conda Environment

Create a conda environment (required for multiqc and running snakemake)

```sh
conda create --name dna-seq -c bioconda multiqc snakemake -y
conda activat dna-seq
```

## Pipeline Steps

1. **Quality Control**: Perform quality control on raw FASTQ files using FastQC and MultiQC.
2. **Read Alignment**: Align reads to the reference genome using Bowtie2 and sort the alignments using Samtools.
3. **Mark Duplicates**: Mark duplicate reads using GATK.
4. **Base Quality Score Recalibration**: Perform base quality score recalibration using GATK.
5. **Variant Calling**: Call variants using VarScan.
6. **Annotate Variants**: Annotate variants using Annovar.
7. **Filter Variants**: Filter variants based on specific criteria using a Python or R script.
   1. Currently R script is broken
