
library(maftools)

# Read the MAF file using read.maf()
maf_data <- read.maf("Annotated_with_Classificationtest.maf")

# Visualize the somatic mutations using a waterfall plot
#oncoplot(maf_data, 
#         top = 20,
#         draw_titv = TRUE)  # Adjust settings as necessary

# Filter out SNP variants (assuming 'Variant_Classification' contains SNP information)
# 
#snp_maf <- subset(maf_data, Variant_Classification %in% c("missense_variant", "synonymous_variant", "downstream_gene_variant", "intron_variant"))

# You can also select only the gene names (Hugo_Symbol)
#snp_genes <- snp_maf$Hugo_Symbol

# Show top 20 genes based on the number of SNPs
#top_snp_genes <- head(sort(table(snp_genes), decreasing = TRUE), 20)

# Create the Oncoplot with the top 10 genes showing SNP data
oncoplot(maf_data, 
         top = 10, 
         altered = TRUE,
         draw_titv = TRUE,
         showTumorSampleBarcodes = TRUE
)
