# Omics Pipelines

Omics pipelines for analysing raw sequencing files.

Please check individual folders for more details on how to use.

## Overview

|   | DNA-Seq | Alt-DNA-Seq | Local-Conda-Seq |
|---|---|---|---|
| Working? | Yes | Yes | Yes z
| Description | DNA seq analysis pipeline | Alternative DNA seq analysis pipeline | DNA Seq pipeline entirely on conda env able to run on apple silicon mac |
| Tested | On HPC | On HPC | On M3 Macbook air |
| Conda env | For multiqc | For multiqc | Yes |
| Workflow | bash, nextflow, snakemake | bash | snakemake |

## Modules Used

|   | DNA-Seq | Alt-DNA-Seq | Local-Conda-Seq |
|---|---|---|---|
| Alignment | Bowtie2 | Bwa | Bwa |
| Variant Calling | VarScan | GATK Mutect2 | GATK Mutect2 |
| Variant Annotation | Annovar | Ensembl-VEP | snpEff |
| Analysis | Python | Rmd | Rmd |
