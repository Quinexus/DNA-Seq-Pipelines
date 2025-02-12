# DNA-Seq Pipeline

This folder contains scripts and workflows for running a DNA-Seq analysis pipeline. 

The pipeline can be executed using different workflow managers including Bash, Snakemake, and Nextflow.

This 

## Pipeline Steps

1. **Quality Control**: Perform quality control on raw FASTQ files using FastQC and MultiQC.
2. **Read Alignment**: Align reads to the reference genome using Bowtie2 and sort the alignments using Samtools.
3. **Mark Duplicates**: Mark duplicate reads using GATK.
4. **Base Quality Score Recalibration**: Perform base quality score recalibration using GATK.
5. **Variant Calling**: Call variants using VarScan.
6. **Annotate Variants**: Annotate variants using Annovar.
7. **Filter Variants**: Filter variants *automatically* based on specific criteria using a Python. *Note: Currently R script is broken*

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

### Nextflow Pipeline

To run the DNA-Seq pipeline using Nextflow, use the following command:

```sh
nextflow run dna-seq.nf
```

## Configuration

The following configuration parameters should be adjusted before running the above pipelines

- *Read Details:* Location of raw fastq.gz paired reads
- *Parameters*
  - threads
  - ram_usage
- *Hardcoded References*
  - bowtie_index (generated using bowtie index)
  - reference (.fa file)
  - known_variants (.vcf file)
  - humandb_path (for annotation)
- *Tools*
  - varscan_location
  - varscan (name of jar file)

### Location of configuration data

**Bash:**

- Pass in Read Details and Parameters as arguments
- Hardcoded References and tools hardcoded into shell script
  
**Snakemake:** Pass in all of the above in `config.yaml` file

**Nextflow:** Pass in all of the above at top of `dna-seq.nf`


## Software and Versions (based on QMUL HPC Apocrita)

- **FastQC**: 0.11.9
- **Bowtie2**: 2.4.5
- **Samtools**: 1.19.2
- **GATK**: 4.5.0.0
- **VarScan**: 2.4.3
- **Annovar**: 2018Apr16

### Conda Environment

Create a conda environment (required for multiqc and running snakemake)

```sh
conda create --name dna-seq -c bioconda multiqc snakemake -y
conda activat dna-seq
```