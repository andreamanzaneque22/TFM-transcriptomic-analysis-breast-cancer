# Transcriptomic Analysis of Breast Cancer

## Overview

This repository contains the code and resources developed for a Master's Thesis in Bioinformatics focused on the transcriptomic analysis of RNA-seq data to investigate the association between gene expression profiles and response to neoadjuvant treatment in breast cancer.

The main objective is to characterize transcriptomic patterns associated with treatment response and explore potential molecular biomarkers and biological pathways related to response.

## Objectives

- Analyze RNA-seq data to characterize transcriptomic profiles associated with treatment response.
- Perform quality control, filtering, normalization, and exploratory data analysis.
- Investigate differential gene expression between response groups.
- Explore biological pathways associated with treatment response.
- Evaluate exploratory machine learning approaches for response classification.
- Identify candidate molecular biomarkers for further investigation.

## Repository Structure

```text
TFM-transcriptomic-analysis-breast-cancer/
├── datos/
│   ├── RNAseq-rawcounts.tsv.gz
│   └── Supplementary Tables 1-6.xlsx
├── resultados_TFM/
├── resultados_principales/
├── FINAL_CODIGO.R
├── README.md
└── .gitignore
```

- `datos/`: Input RNA-seq counts and supplementary clinical information.
- `FINAL_CODIGO.R`: Main R script containing the analysis pipeline.
- `resultados_TFM/`: Directory for intermediate and generated analysis outputs.
- `resultados_principales/`: Selected figures and tables intended to summarize the main findings.

## Data

The analysis uses publicly available RNA-seq data and supplementary information from the TransNEO study.

The input files are provided in the `datos/` directory to facilitate access to the data required for the analysis.

## Analysis Workflow

The analysis pipeline includes the following steps:

1. Import of RNA-seq raw counts and clinical metadata.
2. Quality control of samples and gene expression data.
3. Filtering of lowly expressed genes.
4. TMM normalization and log-CPM transformation.
5. Exploratory analysis, including principal component analysis (PCA).
6. Differential expression analysis using limma-voom and DESeq2.
7. Gene set and pathway analysis using GSVA.
8. Exploratory machine learning analyses.
9. Generation of figures, tables, and summary outputs.

## Software Requirements

The analysis was performed using R version 4.5.3.

The following package versions were recorded in the original analysis environment:

| Package | Version |
|---|---:|
| edgeR | 4.8.2 |
| limma | 3.66.0 |
| DESeq2 | 1.50.2 |
| GSVA | 2.4.9 |
| msigdbr | 26.1.1 |
| caret | 7.0.1 |
| pROC | 1.19.0.1 |
| glmnet | 5.0 |
| ranger | 0.18.0 |
| e1071 | 1.7.17 |
| pheatmap | 1.0.13 |
| corrplot | 0.95 |
| igraph | 2.3.3 |

Additional packages may be required by the script.

## Reproducibility

To run the analysis:

1. Install R 4.5.3 or a compatible version.
2. Download or clone this repository.
3. Ensure that the input files are available in the `datos/` directory.
4. Open `FINAL_CODIGO.R` in RStudio.
5. Set the working directory to the repository root.
6. Install the required packages and check their versions.
7. Run the script.

The script uses relative paths to locate input files and organize output directories. The exact results may depend on the software environment and package versions.

## Results

Selected figures and tables will be made available in the `resultados_principales/` directory.

The results are intended to characterize associations between transcriptomic profiles, biological pathways, and treatment response.

## Limitations

The findings should be interpreted in the context of the cohort size, available clinical information, and statistical limitations of the analyses.

Any candidate biomarkers identified in this study require further validation in independent datasets before potential clinical application.

## Author

**Andrea Manzaneque Pérez**  
Master's Degree in Bioinformatics

This repository contains the code and selected results from my Master's thesis, focused on transcriptomic analysis of RNA-seq data to investigate molecular signatures associated with response to neoadjuvant treatment in breast cancer.
