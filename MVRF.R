# Load necessary libraries
library(dplyr)

# Read the CSV file
data <- read.csv("test.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)

# Convert relevant columns to numeric
numeric_cols <- c("X.000", "X.001", "X.010", "X.011", "X.100", "X.101", "X.110", "X.111")  # Adjust column names if necessary
data[numeric_cols] <- lapply(data[numeric_cols], as.numeric)

# Calculate Total Reads and Variant Counts
data <- data %>%
  mutate(
    Total_Reads = rowSums(select(., all_of(numeric_cols))),
    Variant_Counts = rowSums(select(., c("X.001", "X.010", "X.011", "X.100", "X.101", "X.110", "X.111"))),
    VRF = Variant_Counts / Total_Reads
  )

# Group by Tissue and CpG_Coordinates, and calculate mean VRF
mean_vrf <- data %>%
  group_by(Tissue, CpG_Coordinates) %>%
  summarise(Mean_VRF = mean(VRF, na.rm = TRUE), .groups = "drop")

# View results
print(mean_vrf)

# Optionally save the results to a CSV
write.csv(mean_vrf, "mean_vrf_results.csv", row.names = TRUE)

