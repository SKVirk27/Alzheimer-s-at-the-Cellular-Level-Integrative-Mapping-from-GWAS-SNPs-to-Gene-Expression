# GWAS Data Formatting and Quality Control Script

# Load necessary libraries
library(EPIC)
library(dplyr)
library(data.table)
library(readr)
library(utils)

# Set file paths
bfile_path <- "path/to/bfile"
plink_path <- "path/to/plink"
input_file_path <- "path/to/AD_sumstats_file.txt"

# Read and prepare the GWAS summary statistics data
gwas_data <- fread(input_file_path, header = TRUE, sep = "\t")
colnames(gwas_data) <- c("chr", "pos", "A1", "A2", "rsid", "Zscore", "beta", "P", "EAF", "N", "se")

# Perform quality control and alignment
gwas_formatted <- format_gwas(gwas_data, bfile_path, plink_path)
gwas_formatted$MAC <- 2 * gwas_formatted$N * gwas_formatted$MAF

# Save the formatted GWAS data
saveRDS(gwas_formatted, "gwas_formatted.rds")

# Display the tail of the formatted GWAS data
tail(gwas_formatted)

# Optional: Read back the saved RDS file if needed
gwas_formatted <- readRDS("path/to/gwas_formatted.rds")

# Perform liftover from hg19 to hg38 using the ieugwasr package
# Make sure to install and load the ieugwasr package if not already done
library(ieugwasr)

# Liftover the GWAS data to the hg38 build
Gwas_lifting <- liftover_gwas(
  gwas_data,
  build = c("37", "38"),
  to = "38",
  chr_col = "chr",
  pos_col = "pos",
  snp_col = "rsid",
  ea_col = "A1",
  oa_col = "A2"
)

# Inspect the lifted GWAS data
head(Gwas_lifting)

# Citation for EPIC package and liftover method
# Add the citation for the EPIC package and the liftover package here as per the provided references

# Save the script as an R script file
writeLines(capture.output(sink()), "gwas_data_formatting_script.R")
