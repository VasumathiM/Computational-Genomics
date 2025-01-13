
library(dplyr)
library(tidyr)  # Load tidyr for separate_rows function
library(ggplot2)

# Step 1: Load Data
data <- read.csv("test", sep = ",", header = TRUE, stringsAsFactors = FALSE)

# Step 2: Filter Data for cfDNA Tissue
data_cfDNA <- data %>%
  filter(Tissue == "cfDNA")
state_cols <- c("X.000", "X.001", "X.010", "X.011", "X.100", "X.101", "X.110", "X.111")
# Step 3: Process Individual CpG Sites
data_individual <- data_cfDNA %>%
  separate_rows(CpG_Coordinates, sep = ":") %>%
  mutate(row_sum = rowSums(select(., all_of(state_cols)), na.rm = TRUE),  # Add row_sum
         Variant_Counts = X.001 + X.010 + X.011 + X.100 + X.101 + X.110 + X.111)  # Calculate Variant_Counts

# Step 4: Calculate VRF for Individual CpG Sites
data_individual_summary <- data_individual %>%
  group_by(CpG_Coordinates) %>%
  summarise(
    Total_Reads = sum(row_sum, na.rm = TRUE),
    Variant_Counts = sum(Variant_Counts, na.rm = TRUE),  # Sum Variant_Counts
    VRF = Variant_Counts / Total_Reads,
    .groups = "drop"
  )
# Step 5: Calculate VRF for Top 10 PMPs
data_pmp_summary <- data_cfDNA %>%
  mutate(row_sum = rowSums(select(., all_of(state_cols)), na.rm = TRUE),  # Calculate row_sum for cfDNA
         Variant_Counts = X.001 + X.010 + X.011 + X.100 + X.101 + X.110 + X.111) %>%
  mutate(VRF = Variant_Counts / row_sum) %>%
  arrange(desc(row_sum)) %>%
  slice(1:10) %>%
  select(CpG_Coordinates, VRF)

# Step 6: Statistical Test
# Combine PMP and CpG data for comparison
data_combined <- bind_rows(
  data_individual_summary %>% mutate(Group = "Individual CpG"),
  data_pmp_summary %>% mutate(Group = "Top 10 PMPs")
)

# Perform Wilcoxon Rank-Sum Test
wilcox_test <- wilcox.test(
  data_combined %>% filter(Group == "Top 10 PMPs") %>% pull(VRF),
  data_combined %>% filter(Group == "Individual CpG") %>% pull(VRF),
  alternative = "greater"
)
print(paste("Wilcoxon test p-value:", wilcox_test$p.value))

# Step 7: Visualization
comparison_plot <- ggplot(data_combined, aes(x = Group, y = VRF, fill = Group)) +
  geom_boxplot() +
  labs(
    title = "VRF Comparison: Top 10 PMPs vs Individual CpG Sites",
    x = "Group",
    y = "Variant Read Frequency (VRF)"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "none"
  )

# Save the Plot
ggsave("vrf_comparison_pmp_vs_cpg.png", plot = comparison_plot, width = 8, height = 6)

# Step 8: Interpret Results
if (wilcox_test$p.value < 0.05) {
  print("The hypothesis is validated: Top 10 PMPs exhibit higher specificity than individual CpG sites.")
} else {
  print("The hypothesis is not validated: No significant difference in specificity between top 10 PMPs and individual CpG sites.")
}

