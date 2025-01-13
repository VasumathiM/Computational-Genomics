library(dplyr)
# Load the data
data <- read.csv("test", sep = ",", header = TRUE, stringsAsFactors = FALSE)

# Convert relevant columns to numeric
numeric_cols <- c("X.000", "X.001", "X.010", "X.011", "X.100", "X.101", "X.110", "X.111")  # Adjust column names if needed
data[numeric_cols] <- lapply(data[numeric_cols], as.numeric)

# Calculate Sequencing Depth
data$Depth <- rowSums(data[, numeric_cols])

# Optionally, filter data for Tissue #2 
tissue_data <- subset(data, Tissue == "Islet")

# Sort the tissue data by Depth and select the top 10
top_10_tissue_data <- tissue_data %>%
  arrange(desc(Depth)) %>%
  head(10)

# Display the top 10 entries for CpG coordinates and Depth
print(top_10_tissue_data[, c("CpG_Coordinates", "Depth")])

# Estimate the threshold for confident variant call
# Define parameters for variant call threshold (based on 1 million reads)
total_reads <- 1000000  # Sequencing depth (1 million reads)
confidence_level <- 0.95  # Desired confidence (95%)
vaf_threshold <- 0.05  # Assuming a variant allele frequency threshold of 5%

# Calculate the minimum number of variant reads required for confidence
min_variant_reads <- ceiling(total_reads * vaf_threshold)

# Binomial probability function to calculate p-value for detecting variants
p_value <- function(variant_reads, total_reads) {
  # Calculate p-value for binomial test (detecting variant_reads or more variants)
  pbinom(variant_reads - 1, total_reads, vaf_threshold, lower.tail = FALSE)
}

# Find the threshold reads where p-value is less than 0.05 (for 95% confidence)
threshold_reads <- min_variant_reads
while (p_value(threshold_reads, total_reads) > (1 - confidence_level)) {
  threshold_reads <- threshold_reads + 1
}

# Display the threshold number of variant reads required for 95% confidence
print(paste("Threshold reads required for confident variant call:", threshold_reads))

