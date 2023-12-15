# Integration of GWAS and scRNA-seq Data for Alzheimer’s Disease Research

# Load required libraries
library(dplyr)

# Set the file paths for input data
gwas_data_path <- "path/to/GWAS_data.csv"
scRNAseq_data_path <- "path/to/scRNAseq_data.csv"
output_path <- "path/to/final_merged_data.csv"

# Read in the GWAS data
gwas_data <- read.csv(gwas_data_path, stringsAsFactors = FALSE)

# Read in the Single-cell RNA-seq data
scRNAseq_data <- read.csv(scRNAseq_data_path, stringsAsFactors = FALSE)

# Merge the datasets on the gene identifier column
merged_data <- merge(gwas_data, scRNAseq_data, by = "gene_identifier")

# Filter for significant genes and annotate cell types
filtered_data <- merged_data %>%
  filter(p_value < 0.05) %>%
  mutate(cell_type = case_when(
    cell_type_condition == "AD" ~ "Alzheimer's",
    cell_type_condition == "HC" ~ "Healthy Control",
    TRUE ~ as.character(cell_type_condition)
  ))

# Write the final merged data to a CSV file
write.csv(filtered_data, output_path, row.names = FALSE)

# Print completion message
message("Data integration script has finished executing.")
