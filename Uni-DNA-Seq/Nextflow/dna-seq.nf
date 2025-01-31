#!/usr/bin/env nextflow

/* 
 * CONFIG PARAMS
 */

// Reads
params.prefix = "tumour"
params.reads_r1 = "/data/home/ha20830/Practice/FASTQ_Raw/tumour_R1.fq.gz"
params.reads_r2 = "/data/home/ha20830/Practice/FASTQ_Raw/tumour_R2.fq.gz"

// Parameters
params.threads = 20
params.ram_usage = 15

// Reference Files
params.bowtie_index = "/data/teaching/bci_teaching/DNAseq/Reference/Bowtie2Idx/GRCh38.108.chr17"
params.reference = "/data/teaching/bci_teaching/DNAseq/Reference/Homo_sapiens.GRCh38.108.dna.chromosome.17.fa"
params.known_variants = "/data/teaching/bci_teaching/DNAseq/Reference/gatkResources/resources_broad_hg38_v0_1000G_omni2.5.hg38.noCHR.vcf"
params.humandb = "/data/teaching/bci_teaching/DNAseq/Reference/humandb/"

// Tools
params.varscan = "/data/home/ha20830/Practice/VarScan.v2.4.3.jar"

/*
 * WORKFLOW
 */

workflow {
    // quality_control(params.reads_r1, params.reads_r2) // Java conflict!

    read_alignment(params.reads_r1, params.reads_r2) | mark_duplicates | base_score_quality_recalibration | variant_calling | annotate_variants | filter_variants
}

process quality_control {
    input:
        path reads_r1
        path reads_r2 
    output:
        path "QC/multiqc_report.html"
    script:
        """
        fastqc $reads_r1 $reads_r2 -o QC/ -t $params.threads
        multiqc QC/ -o QC
        """ 
}

process read_alignment {
    input:
        path reads_r1
        path reads_r2
    output:
        path "Alignment/${params.prefix}.bam"
    script:
        """
        module load bowtie2 samtools

        mkdir -p Alignment

        bowtie2 -p ${params.threads} \
            --rg ID:${params.prefix} \
            --rg SM:${params.prefix} \
            --rg PL:ILLUMINA \
            --rg LB:${params.prefix} \
            -x ${params.bowtie_index} \
            -1 $reads_r1 \
            -2 $reads_r2 |
            samtools sort -o Alignment/${params.prefix}.bam -
        
        module unload bowtie2 samtools
        """
}

process mark_duplicates {
    input:
        path bam
    output:
        path "Alignment/${params.prefix}.marked.bam"
    script:
        """
        module load gatk

        mkdir -p Alignment QC

        gatk --java-options "-Xmx${params.ram_usage}G" MarkDuplicates \
            -I $bam \
            -O Alignment/${params.prefix}.marked.bam \
            -M QC/${params.prefix}.marked
        
        module unload gatk
        """
}

process base_score_quality_recalibration {
    input:
        path marked_bam
    output:
        path "Alignment/${params.prefix}.recalib.bam"
    script:
        """
        module load gatk

        mkdir -p Alignment

        gatk --java-options "-Xmx${params.ram_usage}G" BaseRecalibrator \
            -I $marked_bam \
            -R $params.reference \
            --known-sites $params.known_variants \
            -O Alignment/${params.prefix}.table
        
        gatk --java-options "-Xmx${params.ram_usage}G" ApplyBQSR \
            -R $params.reference \
            -I $marked_bam \
            --bqsr-recal-file Alignment/${params.prefix}.table \
            -O Alignment/${params.prefix}.recalib.bam
        
        module unload gatk
        """
}

process variant_calling {
    input:
        path recalib_bam
    output:
        path "VCF/${params.prefix}.vcf"
    script:
    """
    module load samtools

    mkdir -p VCF

    samtools mpileup \
        -q 20 \
        -f $params.reference \
        $recalib_bam |
    java -jar $params.varscan mpileup2snp \
        --min-coverage 20 \
        --min-avg-qual 20 \
        --min-read2 4 \
        --p-value 0.2 \
        --min-var-freq 0.01 \
        --strand-filter 1 \
        --output-vcf 1 > VCF/${params.prefix}.vcf
        
    module unload samtools
    """
}

process annotate_variants {
    input:
        path vcf
    output:
        path "VCF/${params.prefix}.hg38_multianno.txt"
    script:
        """
        module load annovar

        mkdir -p VCF

        convert2annovar.pl --format vcf4 \
            $vcf \
            --includeinfo \
            --filter PASS \
            --outfile VCF/${params.prefix}.pass.vcf

        annotate_variation.pl -filter \
            -dbtype 1000g2015aug_all \
            -buildver hg38 \
            -out VCF/${params.prefix} \
            VCF/${params.prefix}.pass.vcf \
            ${params.humandb} \
            -maf 0.01

        annotate_variation.pl -filter \
            -dbtype esp6500siv2_all \
            -buildver hg38 \
            -out VCF/${params.prefix} \
            VCF/${params.prefix}.hg38_ALL.sites.2015_08_filtered \
            ${params.humandb} \
            -score_threshold 0.01

        table_annovar.pl \
            -buildver hg38 \
            -out VCF/${params.prefix} \
            VCF/${params.prefix}.hg38_esp6500siv2_all_filtered \
            ${params.humandb} \
            -remove \
            -otherinfo \
            -protocol refgene,avsnp150,cosmic92_coding,cytoband \
            -operation g,f,f,r \
            -nastring .
        
        module unload annovar
        """
}

process filter_variants {
    input:
        path annotated_vcf
    output:
        path "${params.prefix}_filtered.csv"
    script:
        """
        #!/usr/bin/env python

        import pandas as pd
        import sys
        import os

        multianno_file = "${annotated_vcf}"
        # extracting the prefix
        prefix = os.path.basename(multianno_file).replace('.hg38_multianno.txt', '')
        # loading the data
        variants = pd.read_csv(multianno_file, skiprows=1, delimiter='\t', header=None, on_bad_lines='warn')

        # setting column names
        headings = ["Chr","Start","End","Ref","Alt","Func.refgene","Gene.refgene","GeneDetail.refgene","ExonicFunc.refgene","AAChange.refgene","avsnp150","cosmic92_coding","cytoband", 
                    "chr", "position", "id", "ref", "alt", "qual", "filter", "info", "format", "sample"]
        annotated_variants = variants.iloc[1:].copy()
        annotated_variants.columns = headings

        # splitting the format and sample columns
        format_headings = annotated_variants['format'].iloc[0].split(':')
        allele_counts = annotated_variants['sample'].str.split(':', expand=True)
        allele_counts.columns = format_headings
        allele_counts['FREQ'] = allele_counts['FREQ'].str.replace('%', '').astype(float)

        # combining the data
        annotated_variants = pd.concat([annotated_variants, allele_counts], axis=1)

        # filtering based on exonic variants + filtering out synonymous SNVs
        annotated_variants_exonic = annotated_variants[
            (annotated_variants['Func.refgene'] == 'exonic') & 
            (annotated_variants['ExonicFunc.refgene'] != 'synonymous SNV')
        ]

        # filtering based on VAF > 10%
        vaf_10 = annotated_variants_exonic[annotated_variants_exonic['FREQ'] > 10]
        print(prefix + '_filtered.csv')
        vaf_10.to_csv(prefix + '_filtered.csv', index=False)

        """
}
