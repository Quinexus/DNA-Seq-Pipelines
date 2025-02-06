#!/usr/bin/env python

import pandas as pd
import sys
import os

# parsing arguments
multianno_file = sys.argv[1]
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
vaf_10.to_csv(prefix + '_filtered.csv', index=False)
