# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)

print("Load the data")
data <- read.csv("test", sep = ",", header = TRUE, stringsAsFactors = FALSE)

# Filter data for cfDNA tissue
data_cfDNA <- data %>%
  filter(Tissue == "cfDNA")

# Split CpG_Coordinates into groups
data_cfDNA <- data_cfDNA %>%
  mutate(CpG_Group = sapply(strsplit(CpG_Coordinates, ":"), function(x) paste(x[1:2], collapse = ":")))

# Define methylation state columns
state_cols <- c("X.000", "X.001", "X.010", "X.011", "X.100", "X.101", "X.110", "X.111")

# Step 1: Calculate row sums
data_cfDNA <- data_cfDNA %>%
  mutate(row_sum = rowSums(select(., all_of(state_cols)), na.rm = TRUE)) %>%
  mutate(row_sum = ifelse(row_sum == 0, NA, row_sum))  # Replace 0 with NA to avoid division by zero

# Step 2: Calculate Variant Read Frequency (VRF) for each methylation state
data_cfDNA <- data_cfDNA %>%
  mutate(across(all_of(state_cols), 
                ~ . / row_sum, 
                .names = "{.col}_vrf"))

# Step 3: Apply VRF Threshold for High Specificity
# Adjusted Threshold
vrf_threshold <- 0.4  # Relax threshold if needed




data_cfDNA <- data_cfDNA %>%
  mutate(across(ends_with("_vrf"), 
                ~ ifelse(. >= vrf_threshold, "High", "Low"),
                .names = "{.col}_specificity"))

# Partial High-Specificity Filter
high_specificity_pmp <- data_cfDNA %>%
  filter(rowSums(select(., ends_with("_specificity")) == "High", na.rm = TRUE) >= 3)  # Adjust minimum states
 
print(nrow(high_specificity_pmp))
# Step 4: Identify High-Specificity PMPs
#high_specificity_pmp <- data_cfDNA %>%
#  filter(if_all(ends_with("_specificity"), ~ . == "High"))

# Step 5: Plot VRF for Top 10 PMPs
# Select Top 10 PMPs based on row_sum
top_10_pmps <- data_cfDNA %>%
  arrange(desc(row_sum)) %>%
  slice(1:10)
print(head(top_10_pmps)) 
# Gather data for plotting VRF
plot_data <- top_10_pmps %>%
  pivot_longer(cols = ends_with("_vrf"), 
               names_to = "Methylation_State", 
               values_to = "VRF") 

# Plotting
ggplot(plot_data, aes(x = CpG_Group, y = VRF, fill = Methylation_State)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "VRF for Top 10 PMPs",
       x = "CpG Group", y = "Variant Read Frequency (VRF)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme_classic() +
  scale_fill_brewer(palette = "Set3")

# Save the plot to a file
ggsave("vrf_top_10_pmps_plot30.png", width = 10, height = 8)

# Optional: Save data for high-specificity PMPs
write.csv(high_specificity_pmp, "high_specificity_pmp_data_vrf_30.csv", row.names = FALSE)

print("VRF-based specificity analysis completed successfully.")

