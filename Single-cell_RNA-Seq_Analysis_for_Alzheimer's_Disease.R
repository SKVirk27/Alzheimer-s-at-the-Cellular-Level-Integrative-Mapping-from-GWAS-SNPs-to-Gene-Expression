# Single-cell RNA-Seq Analysis for Alzheimer's Disease Study

# Load necessary libraries
library(Seurat)
library(SingleR)
library(ggplot2)
library(tidyr)
library(dplyr)
library(DoubletFinder)
library(scran)
library(scRNAseq)

# Set the working directory
setwd("/path/to/your/data")

# Load Alzheimer's and control samples
# Define sample paths and prefixes
AD_samples <- c("AD04", "AD12", "AD16", "AD30")
HC_samples <- c("HC03", "HC07", "HC14", "HC19")
all_samples <- c(AD_samples, HC_samples)
file_prefixes <- c(...)  # Define the file prefixes for each sample

# Create Seurat objects for each sample
seurat_objects <- list()
for (sample in all_samples) {
  prefix <- file_prefixes[[sample]]
  seurat_objects[[sample]] <- Seurat::CreateSeuratObject(...)
}

# Merge and Quality Control
# Merge Seurat objects
merged_seurat <- Reduce(function(x, y) merge(x, y), seurat_objects)

# Quality control - mitochondrial, ribosomal, hemoglobin content
merged_seurat <- PercentageFeatureSet(merged_seurat, pattern = "^MT-", col.name = "percent_mito")
merged_seurat <- PercentageFeatureSet(merged_seurat, pattern = "^RP[SL]", col.name = "percent_ribo")
merged_seurat <- PercentageFeatureSet(merged_seurat, pattern = "^HB[^(P)]", col.name = "percent_hb")
# Data Filtration and Quality Control

# Normalize, identify variable features, scale data, and run PCA
merged_seurat <- NormalizeData(merged_seurat)
merged_seurat <- FindVariableFeatures(merged_seurat)
merged_seurat <- ScaleData(merged_seurat)
merged_seurat <- RunPCA(merged_seurat)
filtered_seurat <- subset(data.filt, subset = nCount_RNA > 1000 & nCount_RNA < 22850 & nFeature_RNA > 200 & nFeature_RNA < 6037)
data.filt <- filtered_seurat 
# Doublet Detection with DoubletFinder
# Doublet Detection with DoubletFinder
sweep.res.list <- DoubletFinder::paramSweep_v3(data.filt, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK <- bcmvn %>% filter(BCmetric == max(BCmetric)) %>% select(pK) %>% as.numeric(as.character(.))
homotypic.prop <- modelHomotypic(data.filt@meta.data$seurat_clusters)
nExp_poi.adj <- round(0.076 * nrow(data.filt@meta.data) * (1 - homotypic.prop))

data.filt <- doubletFinder_v3(data.filt, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
saveRDS(data.filt, "data_after_doublet_removal.rds")

# Data Integration
# Split object by patient or sample
object_list <- SplitObject(data.filt, split.by = 'orig.ident')

# Find integration anchors and integrate data
anchors <- FindIntegrationAnchors(object.list = object_list, ...)
data.filt <- IntegrateData(anchorset = anchors)

# UMAP visualization
data.filt <- RunUMAP(merged_seurat, dims = 1:20)



# Cell Type Annotation with SingleR
# Load reference datasets
# Datasets choosen is mostly related to brain 
hpca <- celldex::HumanPrimaryCellAtlasData()
sceBacher <- BacherTCellData()


sceDarmanis <- DarmanisBrainData()


sceZhong <- ZhongPrefrontalData()


original_counts <- GetAssayData(data.filt, assay = "RNA", slot = "counts")
# View cell metadata for sceZhong
cell_metadata_Zhong <- colData(sceZhong)
head(cell_metadata_Zhong)

# Run SingleR with the additional references
# Run SingleR with the additional references
resultlabel <- SingleR(
  test = original_counts, 
  ref = list(
    HPCA = hpca,
    Bacher = counts(sceBacher),
    Darmanis = counts(sceDarmanis),
    Zhong = counts(sceZhong)
  ),
  labels = list(
    hpca$label.main,
    sceBacher$new_cluster_names,
    sceDarmanis$cell.type,
    
    sceZhong$cell_types 
  ),
  de.method="wilcox"
)
resultlabel$labels <- tolower(resultlabel$labels)
# Example: Convert 'HSC_-G-CSF' to 'hsc-g-csf'
resultlabel$labels <- gsub("_-", "-", resultlabel$labels)

# Construct a table to have a look cell types annotated
table(resultlabel$labels)

# Correcting the labels
correct_labels <- function(labels) {
  # Replacements
  replace_map <- list(
    "astrocytes" = "astrocyte",
    "endothelial_cells" = "endothelial"
  )
  
  # Loop through and replace
  for (i in seq_along(labels)) {
    if (labels[i] %in% names(replace_map)) {
      labels[i] <- replace_map[[labels[i]]]
    }
  }
  return(labels)
}

# Applying the correction function
resultlabel$labels <- correct_labels(resultlabel$labels)

# Construct a table to have a look at the reorganized cell types annotated
table(resultlabel$labels)


# Visualization
DimPlot(merged_seurat, reduction = "umap")

# End of script
