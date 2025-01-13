# Load maftools library
library(maftools)

# Read the MAF file using read.maf()
maf_data <- read.maf("Annotated_with_Classificationtest.maf")

# Visualize the somatic mutations using a waterfall plot
#oncoplot(maf_data, 
#         top = 20,
#         draw_titv = TRUE)  # Adjust settings as necessary


plotmafSummary(maf_data,  showBarcodes = TRUE, top=10)
