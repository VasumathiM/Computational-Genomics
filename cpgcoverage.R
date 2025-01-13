# Load necessary libraries
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
# Load the CSV file
data <- read.csv("test.csv", sep=",", header=T)
colnames(data) <- trimws(colnames(data))
coverage_cols <- c("X.000", "X.001", "X.010", "X.011", "X.100", "X.101", "X.110", "X.111")

# Group by Tissue and CpG_Coordinates, then calculate median and CV
summary_stats <- data %>%
  group_by(Tissue, CpG_Coordinates) %>%
  summarise(
    Median = median(unlist(across(all_of(coverage_cols)))),
    Mean = mean(unlist(across(all_of(coverage_cols)))),
    SD = sd(unlist(across(all_of(coverage_cols)))),
    CV = ifelse(Mean > 0, SD / Mean, NA)  # Avoid division by zero
  )

# Save the summarized data to a new CSV
write.csv(summary_stats, "CpG_Summary_Statistics.csv", row.names = FALSE)

# Plot 1: Median coverage by Tissue
p1 <- ggplot(summary_stats, aes(x = Tissue, y = Median, fill = Tissue)) +
  geom_boxplot() +
  labs(title = "Median Coverage by Tissue", x = "Tissue", y = "Median Coverage") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot 2: CV by Tissue
p2 <- ggplot(summary_stats, aes(x = Tissue, y = CV, fill = Tissue)) +
  geom_boxplot() +
  labs(title = "Coefficient of Variation by Tissue", x = "Tissue", y = "CV") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Display plots
grid.arrange(p1, p2, nrow = 1)
