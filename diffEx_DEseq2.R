library(DESeq2)
library(ggplot2)
library(dplyr)

counts <- read.delim("GSE212787_raw_counts_GRCh38.p13_NCBI.tsv", row.names = 1)
counts <- as.matrix(counts)  # Ensure it's a matrix

#The metadata file should contain sample names and condition labels (e.g., Control vs. Treatment).
metadata <- read.csv("GSE212787_metadata.csv", row.names = 1, sep= ",", header=T)
rownames(metadata) <- gsub('^\"|\"$', '', rownames(metadata))
metadata[] <- lapply(metadata, function(x) if(is.character(x)) gsub('"""', '', x) else x)
all(colnames(counts) %in% rownames(metadata))  # Should return TRUE

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ DiseaseState  # Adjust based on your experimental design
)

dds <- DESeq(dds)
#metadata$Disease_State <- factor(metadata$Disease_State, levels = c("disease state: control", "disease state: endometriosis"))

res <- results(dds, contrast = c("DiseaseState", " endometriosis ", " control "))
res_sorted <- res[order(res$padj), ]

#5.3: Filter Significant Genes (Adjusted p-value < 0.05)
sig_genes <- subset(res_sorted, padj < 0.05)
#7. Save Results
write.csv(as.data.frame(res_sorted), "DESeq2_results.csv")
# Convert DESeq2 results to a data frame
res_df <- as.data.frame(res)

# Remove NA values (optional but recommended)
res_df <- res_df %>% filter(!is.na(padj))

# Define significance threshold
padj_cutoff <- 0.05
lfc_cutoff <- 1  # Adjust based on biological relevance

upregulated_genes <- res_df %>%
  filter(log2FoldChange > lfc_cutoff & padj < padj_cutoff)

# Extract downregulated genes (log2FoldChange < -1 & padj < 0.05)
downregulated_genes <- res_df %>%
  filter(log2FoldChange < -lfc_cutoff & padj < padj_cutoff)
write.table(upregulated_genes, "upregulated_genes.txt", sep="\t", quote=FALSE, row.names=TRUE)
write.table(downregulated_genes, "downregulated_genes.txt", sep="\t", quote=FALSE, row.names=TRUE)


# Add a new column for labeling significance
res_df <- res_df %>%
  mutate(Significance = case_when(
    log2FoldChange > lfc_cutoff & padj < padj_cutoff ~ "Upregulated",
    log2FoldChange < -lfc_cutoff & padj < padj_cutoff ~ "Downregulated",
    TRUE ~ "Not Significant"
  ))

# Add gene names (optional, if row names contain gene IDs)
res_df$Gene <- rownames(res_df)

top_genes <- res_df %>%
  filter(padj < padj_cutoff & abs(log2FoldChange) > 4)  # Adjust threshold as needed
# volcanoplot
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.8, size = 2) +
  geom_text(data = top_genes, aes(label = Gene), vjust=1.2, size=3, check_overlap = TRUE) +
  scale_color_manual(values = c("Upregulated" = "darkgreen", "Downregulated" = "red", "Not Significant" = "gray")) +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "black") +
  labs(title = "GSE212787: Control vs. Endometriosis",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +theme_bw()
