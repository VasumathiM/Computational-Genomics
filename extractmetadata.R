# Read the Series Matrix file
# Read the file as a character vector (not a data frame)
metadata_raw <- readLines("GSE135485_series_matrix.txt")

# Extract Sample IDs
sample_ids <- unlist(strsplit(metadata_raw[grep("^!Sample_geo_accession", metadata_raw)], "\t"))[-1]

# Extract Sample Titles
sample_titles <- unlist(strsplit(metadata_raw[grep("^!Sample_title", metadata_raw)], "\t"))[-1]

# Extract Disease State
disease_state <- unlist(strsplit(metadata_raw[grep("tissue:", metadata_raw)], "\t"))[-1]

# Create a metadata table
metadata <- data.frame(
  Sample_ID = sample_ids,
  Title = sample_titles,
  Disease_State = disease_state,
  stringsAsFactors = FALSE
)

# Display first few rows
head(metadata)

write.csv(metadata, "GSE135485_metadata.csv", row.names = FALSE)

