# Load necessary libraries
library(MAGMA.Celltyping)
library(data.table)

# Reading the GWAS formatted data
gwas_Alz <- readRDS("/path/to/gwas_formatted.rds")

# Renaming columns to match MAGMA.Celltyping requirements
setnames(gwas_Alz, old = c("rsid", "chr", "pos", "Zscore", "beta", "se", "P", "MAF", "N", "MAC", "A1", "A2"),
         new = c("SNP", "CHR", "BP", "Z", "beta", "se", "P", "MAF", "N", "MAC", "A1", "A2"))

# Saving the reformatted GWAS data as a TSV file
path_formatted <- "/path/to/gwas_Alz_formatted.tsv"
fwrite(gwas_Alz, file = path_formatted, sep = "\t")

# Mapping SNPs to genes using MAGMA.Celltyping
genesOutPath <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = path_formatted,
  genome_build = "GRCh37"
)

# Optional: Save the output of the gene mapping
fwrite(genesOutPath, file="/path/to/output_genes_mapped.tsv", sep="\t")


