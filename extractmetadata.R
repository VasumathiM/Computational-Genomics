# Read the Series Matrix file
metadata_raw <- read.delim("GSE212787_series_matrix.txt", header = FALSE, stringsAsFactors = FALSE)

# Extract Sample IDs
sample_ids <- unlist(strsplit(metadata_raw[grep("!Sample_geo_accession", metadata_raw$V1), ], "\t"))[-1]

# Extract Sample Titles (Descriptive Labels)
sample_titles <- unlist(strsplit(metadata_raw[grep("!Sample_title", metadata_raw$V1), ], "\t"))[-1]

# Extract Disease State
disease_state <- unlist(strsplit(metadata_raw[grep("disease state", metadata_raw$V1), ], "\t"))[-1]

# Create a metadata table
metadata <- data.frame(
  Sample_ID = sample_ids,
  Title = sample_titles,
  Disease_State = disease_state,
  stringsAsFactors = FALSE
)

# Display first few rows
head(metadata)

write.csv(metadata, "GSE212787_metadata.csv", row.names = FALSE)
