# Load necessary libraries
library(Seurat)
library(tidyverse)
library(metap)

# Load your integrated Seurat object with cell type annotations
seurat_object <- readRDS("path/to/your/seurat_integrated_celltypes.rds")

# Assuming 'seurat_object' is your Seurat object, combine cell type and condition for DE analysis
seurat_object$celltype.cnd <- paste0(seurat_object$cell_types, '_', seurat_object$Type)

# Update cell identities in the Seurat object for the DE analysis
Idents(seurat_object) <- seurat_object$celltype.cnd

# Retrieve unique cell type and condition combinations for DE analysis
unique_celltype_conditions <- unique(seurat_object$celltype.cnd)

# Initialize a list to store DE analysis results
de_results_list <- list()
min_cells_per_group <- 4  # Minimum number of cells per group for DE analysis

# Perform DE analysis across all unique cell type and condition combinations
for (celltype_condition in unique_celltype_conditions) {
  # Perform DE analysis if both conditions have sufficient cell numbers
  ident1 <- paste0(cell_type, '_Alzheimer')
  ident2 <- paste0(cell_type, '_HealthyControl')
  
  if (sum(Idents(seurat_object) == ident1) >= min_cells_per_group & sum(Idents(seurat_object) == ident2) >= min_cells_per_group) {
    markers <- FindMarkers(seurat_object, ident.1 = ident1, ident.2 = ident2)
    de_results_list[[celltype_condition]] <- markers
    write.csv(markers, file = paste0("DE_markers_", celltype_condition, ".csv"))
  } else {
    cat("Skipping ", celltype_condition, " due to insufficient cell numbers.\n")
  }
}

# Combine and save DE results from all comparisons
combined_results <- bind_rows(lapply(de_results_list, function(x) x$markers))
write.csv(combined_results, file = "combined_DE_results.csv", row.names = FALSE)

# Print a message when the script completes
cat("Differential expression analysis completed and results saved.\n")
