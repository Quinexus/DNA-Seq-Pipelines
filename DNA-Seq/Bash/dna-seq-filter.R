#!/usr/bin/env Rscript

# loading required libraries
if (!require("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}

# parsing arguments
args <- commandArgs(trailingOnly = TRUE)
multianno_file <- args[1]

variants <- read.delim(multianno_file, header = FALSE)
headings <- c("chr", "position", "id", "ref", "alt", "qual", "filter", "info", "format", "sample") # nolint: line_length_linter.
annotated_variants <- setNames(variants[-1, ], c(variants[1, 1:13] %>% unlist(), headings)) # nolint: line_length_linter.

headings <- str_split(annotated_variants$format[1], ":") %>% unlist()
allele_counts <- str_split(annotated_variants$sample, ":") %>%
  do.call("rbind", .) %>%
  as.data.frame() %>%
  setNames(headings) %>%
  mutate(., FREQ = gsub("%", "", FREQ) %>% as.numeric())

annotated_variants <- cbind(annotated_variants, allele_counts)
separate(annotated_variants, sample, into = headings, sep = ":")

# filtering based on exonic variants + filtering out synonymous SNVs
annotated_variants_exonic <- subset(annotated_variants, Func.refgene == "exonic" & ExonicFunc.refgene != "synonymous SNV") # nolint: line_length_linter.

# filtering based on VAF > 10%
vaf_10 <- subset(annotated_variants_exonic, FREQ > 10)
write.csv(VAF_10, "filtered.csv")