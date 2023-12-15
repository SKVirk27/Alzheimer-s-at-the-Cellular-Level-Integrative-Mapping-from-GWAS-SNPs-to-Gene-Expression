library(tidyverse)
library(LDlinkR)

# Read the CSV file into a dataframe
gene_rsids <- read.csv("/Users//Downloads/astrocytes_genes_rsids.csv", stringsAsFactors = FALSE)

# Add a new column to store the kept rsIDs
gene_rsids$Kept_rsIDs <- vector("list", length = nrow(gene_rsids))

# Set parameters
pop <- "EUR"
r2_threshold <- 0.1
maf_threshold <- 0.01
token <- 'token'  # Replace with your actual token
genome_build <- "grch37"

# Function to perform LD pruning on chunks of rsIDs
perform_ld_pruning <- function(rsids) {
  kept_variants <- vector("list", length = length(rsids))
  
  for (j in seq_along(rsids)) {
    snp_result <- SNPclip(snps = rsids[[j]], 
                          pop = pop, 
                          r2_threshold = r2_threshold, 
                          maf_threshold = maf_threshold, 
                          token = token, 
                          file = FALSE,
                          genome_build = genome_build)
    
    if ("Details" %in% names(snp_result)) {
      kept_variants[[j]] <- snp_result %>% 
        filter(grepl("Variant kept", Details)) %>% 
        pull(RS_Number)
    } else {
      warning(paste("Details column not found for chunk", j))
    }
  }
  
  return(kept_variants)
}

# Loop through each row of the dataframe
for (i in 1:nrow(gene_rsids)) {
  # Extract rsIDs for the current row
  rsids_string <- gene_rsids$rsIDs[i]
  
  # Split the rsIDs string into a vector and then into chunks
  rsids_vector <- unlist(str_split(rsids_string, pattern = ";\\s*"))
  rsids_chunks <- split(rsids_vector, ceiling(seq_along(rsids_vector)/5000))
  
  # Perform LD pruning on each chunk
  kept_rsids_list <- perform_ld_pruning(rsids_chunks)
  
  # Combine kept rsIDs from all chunks
  gene_rsids$Kept_rsIDs[i] <- paste(unlist(kept_rsids_list), collapse = "; ")
}

# Convert the list to a character vector for CSV output
gene_rsids$Kept_rsIDs <- sapply(gene_rsids$Kept_rsIDs, function(x) paste(x, collapse = "; "))

okl# Write the updated dataframe back to a CSV file
write.csv(gene_rsids, "/Users//Downloads/updated_rows_astrocytes_genes_rsids.csv", row.names = FALSE)

head(gene_rsids$Kept_rsIDs)