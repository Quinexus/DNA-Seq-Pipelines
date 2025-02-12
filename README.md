# Omics Pipelines

Omics pipelines for analysing raw sequencing files.

Please check individual folders for more details on how to use.

## Overview

|   | DNA-Seq | Alt-DNA-Seq | Local-Conda-Seq | RNA-Seq |
|---|---|---|---|---|
| Working? | Yes | Yes | Yes | No |
| Description | DNA seq analysis pipeline | Alternative DNA seq analysis pipeline | DNA Seq pipeline entirely on conda env able to run on apple silicon mac | RNA seq analysis pipeline
| Tested | On HPC | On HPC | On M3 Macbook air | Not tested |
| Conda env | For multiqc | For multiqc | Yes | No |
| Workflow | bash, nextflow, snakemake | bash | snakemake | bash |

## Modules Used

|   | DNA-Seq | Alt-DNA-Seq | Local-Conda-Seq | RNA-Seq |
|---|---|---|---|---|
| Alignment | Bowtie2 | Bwa | Bwa | Star |
| Variant Calling | VarScan | GATK Mutect2 | GATK Mutect2 | n/a |
| Variant Annotation | Annovar | Ensembl-VEP | snpEff | n/a |
| Analysis | Python | Rmd | Rmd | n/a |
