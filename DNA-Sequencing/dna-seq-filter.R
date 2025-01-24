#!/usr/bin/env Rscript

# loading required libraries
library(tidyverse)

# parsing arguments
args <- commandArgs(trailingOnly = TRUE)
multianno_file <- args[1]

variants <- read.delim(multianno_file, header = FALSE)
headings <- c("chr", "position","id", "ref", "alt", "qual", "filter", "info", "format" ,"sample")
Annotated_variants <- setNames(variants[-1,], c(variants[1,1:13] %>% unlist(), headings))

headings <- str_split(Annotated_variants$format[1], ":") %>% unlist()
AlleleCounts <- str_split(Annotated_variants$sample, ":") %>% do.call("rbind", .) %>% as.data.frame() %>% setNames(headings) %>% mutate(., FREQ = gsub("%", "", FREQ) %>% as.numeric())

Annotated_variants <- cbind(Annotated_variants, AlleleCounts)
separate(Annotated_variants, sample, into = headings, sep = ":")

# filtering based on exonic variants + filtering out synonymous SNVs
Annotated_variants_exonic <- subset(Annotated_variants, Func.refgene == "exonic" & ExonicFunc.refgene != "synonymous SNV")

# filtering based on VAF > 10%
VAF_10 <- subset(Annotated_variants_exonic, FREQ > 10)
write.csv(VAF_10, “filtered.csv”)