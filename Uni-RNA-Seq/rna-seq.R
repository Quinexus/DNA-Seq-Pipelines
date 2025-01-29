# RNA-Seq Analysis Pipeline with Image Saving

# Parse Command-Line Arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript rna-seq.R <read_counts_folder> <sample_group_file> <output_folder>")
}

# Input arguments
input_dir <- args[1]      # Folder containing RSEM expected counts
sample_group_file <- args[2]  # Sample group file
output_dir <- args[3]     # Output directory

# Create output directory if not exists
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Load Required Libraries
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
required_packages <- c("DESeq2", "AnnotationDbi", "org.Hs.eg.db", "EnhancedVolcano", 
                       "clusterProfiler", "fgsea", "qusage", "msigdbr", 
                       "dplyr", "data.table", "pheatmap", "ggplot2")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, ask = FALSE)
}

library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(EnhancedVolcano)
library(clusterProfiler)
library(fgsea)
library(qusage)
library(msigdbr)
library(dplyr)
library(data.table)
library(pheatmap)
library(ggplot2)

# Load RNA-Seq Read Counts
read_counts_file <- file.path(input_dir, "rsem_counts_allsamples.txt")
raw_counts <- read.delim(read_counts_file, header = TRUE, row.names = 1)
raw_counts <- ceiling(raw_counts)  # Round expected_counts to integers for DESeq2

# Load Sample Groups
sample_data <- read.delim(sample_group_file, header = TRUE)

# Create DESeq2 Dataset
dds <- DESeqDataSetFromMatrix(countData = raw_counts, colData = sample_data, design = ~ group)

# Filter Low-Count Genes
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

# Run DESeq2 Analysis
dds <- DESeq(dds)

# Save Normalized Counts
norm_counts <- counts(dds, normalized = TRUE)
write.csv(norm_counts, file = file.path(output_dir, "normalized_counts.csv"), row.names = TRUE)

# Variance Stabilized Transformation
vsd <- vst(dds, blind = FALSE)

# Save PCA Plot
png(file = file.path(output_dir, "PCA_plot.png"))
plotPCA(vsd, intgroup = c("group"))
dev.off()

# Hierarchical Clustering
sample_dists <- dist(t(assay(vsd)))
png(file = file.path(output_dir, "Hierarchical_Clustering.png"))
plot(hclust(sample_dists, method = "complete"))
dev.off()

# Differential Expression Analysis
res <- results(dds, contrast = c("group", "chemo_sensitive", "chemo_resistant"))
write.csv(as.data.frame(res), file = file.path(output_dir, "DE_sensitiveVSresistant_results.csv"))

# Annotate Results with Gene Symbols
res_df <- as.data.frame(res)
res_annot <- mapIds(org.Hs.eg.db, keys = rownames(res_df), keytype = "ENSEMBL", column = "SYMBOL")
res_df_annot <- cbind(res_df, gene_symbol = res_annot)
write.csv(res_df_annot, file = file.path(output_dir, "DE_sensitiveVSresistant_results_wGeneSymbol.csv"))

# Save Volcano Plot
png(file = file.path(output_dir, "Volcano_plot.png"))
EnhancedVolcano(res_df_annot,
                lab = res_df_annot$gene_symbol,
                x = "log2FoldChange",
                y = "padj")
dev.off()

# Heatmap for Differentially Expressed Genes
de_sign <- subset(res_df_annot, padj < 0.05 & abs(log2FoldChange) > 1)
de_genes <- rownames(de_sign)
scaled_data <- t(scale(t(assay(vsd)[de_genes, ])))

png(file = file.path(output_dir, "Heatmap.png"))
pheatmap(scaled_data, cluster_rows = TRUE, show_rownames = FALSE, cluster_cols = TRUE)
dev.off()

# Over-Representation Analysis (ORA)
go_gene_sets <- msigdbr(species = "human", category = "C5")
msigdbr_t2g <- go_gene_sets %>% dplyr::distinct(gs_name, gene_symbol)

enrich_res <- enricher(gene = de_sign$gene_symbol, TERM2GENE = msigdbr_t2g, pvalueCutoff = 0.05)
write.csv(enrich_res@result, file = file.path(output_dir, "GO_sensitiveVSresistant_up.csv"))

# Pre-Ranked Gene Set Enrichment Analysis (GSEA)
gmt_file <- read.gmt("c5.go.bp.v7.4.symbols.gmt")
de_res_ranked <- res_df_annot[order(res_df_annot$log2FoldChange, decreasing = TRUE), ]
de_ranks <- setNames(de_res_ranked$log2FoldChange, de_res_ranked$gene_symbol)

fgsea_res <- fgsea(gmt_file, de_ranks, minSize = 15, maxSize = 500)
fwrite(fgsea_res, file = file.path(output_dir, "PrerankedGSEA_gobp_sensitiveVSresistant.txt"), sep = "\t")

# Save GSEA Enrichment Plots
png(file = file.path(output_dir, "GSEA_Epithelial_Development.png"))
plotEnrichment(gmt_file[["GOBP_EPITHELIUM_DEVELOPMENT"]], de_ranks) + 
  labs(title = "GOBP_EPITHELIUM_DEVELOPMENT") + 
  theme(text = element_text(size = 20))
dev.off()

png(file = file.path(output_dir, "GSEA_Oxidative_Phosphorylation.png"))
plotEnrichment(gmt_file[["GOBP_OXIDATIVE_PHOSPHORYLATION"]], de_ranks) + 
  labs(title = "GOBP_OXIDATIVE_PHOSPHORYLATION") + 
  theme(text = element_text(size = 20))
dev.off()

cat("RNA-Seq Analysis Completed Successfully! Results saved in:", output_dir, "\n")
