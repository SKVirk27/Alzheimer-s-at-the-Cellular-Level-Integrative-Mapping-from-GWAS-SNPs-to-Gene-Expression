# Alzheimer's at the Cellular Level: Integrative Mapping from GWAS SNPs to Gene Expression

## Abstract
This project investigates the genetic underpinnings of Alzheimer's disease (AD) by integrating Genome-Wide Association Studies (GWAS) with single-cell RNA sequencing (scRNA-seq). The goal is to map genetic variants to specific brain cell types to understand their contributions to AD at a cellular level. Our analysis highlights significant expression changes in astrocytes and microglia linked to key genes such as APOE and SORL1, providing insights into their roles in AD pathology.

## Project Structure
The repository is structured as follows:

### Data Files
- **Combined_Plots.pdf**: Visual representation of combined results.
- **GWAS Data Formatting and Quality Control Script.R**: Script for formatting and quality control of GWAS data.
- **Genetic_Analysis_and_Marker_Integration_for_Alzheimer-1.ipynb**: Jupyter notebook for genetic analysis and marker integration.
- **Integration_of_GWAS_and_scRNA-seq.R**: R script for integrating GWAS data with scRNA-seq data.
- **LD_pruning_LDlink.R**: Script for linkage disequilibrium pruning using LDlinkR.
- **MAGMA_Files**: Various output files from MAGMA analysis.
- **Mapping SNPs to Genes for GWAS and Single-Cell Data Integration.R**: Script for mapping SNPs to genes.
- **Project_report_Methodology_results.pdf**: Detailed project report including methodology and results.
- **Single-cell_RNA-Seq_Analysis_for_Alzheimer's_Disease.R**: R script for analyzing single-cell RNA-Seq data.
- **astrocytes_genes_rsids.csv**: Gene and rsID data for astrocytes.
- **combined_markers.csv**: Combined marker data.
- **differential_markers.R**: Script for differential marker analysis.
- **final_merged_data_on_Gwas.csv**: Merged GWAS data.
- **gwas_celltype_filterted.csv**: Filtered GWAS data by cell type.
- **microglia_genes_rsids.csv**: Gene and rsID data for microglia.
- **oligodendrocytes_genes_rsids.csv**: Gene and rsID data for oligodendrocytes.
- **opc_genes_rsids.csv**: Gene and rsID data for oligodendrocyte progenitor cells (OPCs).
- **significant_genes_with_names.csv**: Significant genes with names.
- **updated_*_genes_rsids.csv**: Updated gene and rsID data for various cell types.

### Key Scripts
- **GWAS Data Formatting and Quality Control Script.R**: Formats and ensures quality control of GWAS data.
- **Genetic_Analysis_and_Marker_Integration_for_Alzheimer-1.ipynb**: Integrates GWAS with scRNA-seq data to identify genetic markers.
- **Integration_of_GWAS_and_scRNA-seq.R**: Combines GWAS data with scRNA-seq for comprehensive analysis.
- **LD_pruning_LDlink.R**: Performs linkage disequilibrium pruning to refine genetic data.
- **Single-cell_RNA-Seq_Analysis_for_Alzheimer's_Disease.R**: Analyzes single-cell RNA-Seq data to identify cell-specific gene expression changes.

## Methodology
### Data Preparation
1. **GWAS Data Formatting**: Ensuring GWAS datasets are formatted correctly and quality controlled.
2. **Single-Cell RNA-Seq Data**: Quality control and annotation of single-cell RNA-seq data using tools like DoubletFinder and SingleR.

### Integration and Analysis
1. **Mapping SNPs to Genes**: Using MAGMA to map SNPs from GWAS data to specific genes.
2. **Differential Expression Analysis**: Identifying differentially expressed genes between AD and control samples using Seurat.
3. **LD Pruning**: Reducing redundancy in genetic data through linkage disequilibrium pruning.
4. **Data Integration**: Combining GWAS and scRNA-seq data to identify cell-type-specific gene expression changes.

### Visualization
- **UMAP Clustering**: Visualizing the cellular composition of brain tissue from AD patients and controls.
- **Heatmaps**: Displaying cell type-specific gene expression changes and significant genes related to AD.

## Results
The integration of GWAS and scRNA-seq data revealed significant cell-type-specific gene expression changes in AD. Key findings include:
- **APOE and SORL1**: Significant expression changes in astrocytes and microglia, respectively.
- **Gene Expression Profiles**: Differential expression profiles highlighting the roles of various genes in AD pathology.

## Conclusion
This project underscores the importance of cell-type-specific genetic expressions in understanding Alzheimer's Disease. The integrative approach bridges the gap between genetic susceptibility and cellular dysfunction, paving the way for potential therapeutic targets.

For detailed methodology, results, and data, please refer to the [Project Report](./Project_report_Methodology_results.pdf).

## How to Use
1. **Clone the repository**: `git clone https://github.com/yourusername/Alzheimer-s-at-the-Cellular-Level.git`
2. **Navigate to the project directory**: `cd Alzheimer-s-at-the-Cellular-Level`
3. **Run scripts as needed** using R or Jupyter Notebook for analysis.

## References
Please refer to the project report for a comprehensive list of references and citations used in this study.

## Contact
For any questions or contributions, please contact Simranjit Kaur Virk at simivk1991@gmail.com

---

This repository provides all the necessary scripts and data files to reproduce the analysis and results presented in the study. Thank you for your interest in this project!
