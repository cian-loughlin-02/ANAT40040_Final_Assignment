# ANAT40040_Final_Assignment
## ERBB2 Amplification and Transcrtiptomic Analysis in Breast Cancer

### Overview
This repository contains a R script which implements a bioinformatics pipeline to analyse breast cancer data from The Cancer Genome Atlas (TCGA)
### Data Sources
All datasets used for this assignment were obtained from the TCGA BRCA Pan-Cancer Atlas (2018):
- Clinical Patient Data
- RNA-seq (RSEM Normalised Counts)
- Copy Number Alteration (CNA) Data
Patient barcodes were matched across all three datasets to ensure accurate matching between clinical and molecular information.

### Analysis Workflow
The analysis workflow can be broken into the following steps:

#### 1. Exploration of Clinical Data
Patient age at diagnosis, tumour stage and overall survival distributions were obtained and visualised to gain insight into the clinical data.
#### 2. Data Integration & Preprocessing
- Clinical, RNA-seq and CNA datasets were matched using the TCGA identifiers.
- ERBB2 amplification statues is obtained from CNA data.
- This amplification status is the encoded as a binary metadata variable.
#### 3. Differential Gene Expression Analysis
- Differences in gene expression between the ERBB2 amplified and non-amplified tumours were evaluated using DESeq2.
  - This includes filtering low-count genes and Variance Stabilising Transformation (VST).
#### 4. Exploratory Data Analysis
Principle Component Analysis (PCA) and heatmaps were used to:
- Visualise expression patterns
- Create sample clusters based on ERBB2 amplification status.
#### 5. Pathway Enrichment Analysis
Differentially expressed genes were evaluated using:
- Gene Ontology Enrichment
- KEGG Pathway Enrichment
- Reactome Pathway Enrichment
All three enrichment methods were used to separately evaluate over- and under-expressed genes.
#### 6. Survival Modelling
A Lasso Regularised Cox Regression model using glmnet, was used to evaluate how the predictive relevance of differentially expressed genes with respect to overall survival

### R Packages Used
- DESeq2 - differential expression analysis
- clusterProfiler - enrichment analysis
- org.Hs.eg.db - gene annotation
- ReactomePA - Reactome Pathway Enrichment
- pathview - KEGG Pathway Enrichment
- glmnet - survival modelling
- survival - survival analysis
- pheatmap - creating heatmap

Each analysis was performed in R using Bioconductor packages. Package versions are specified where necessary to ensure reproducability.

### Output
The above analysis workflow generates the following:
- Differential Expression Analysis Results Table
- PCA Plots and Heatmaps
- GO, KEGG and Reactome Pathway Enrichment Plots
- Cross-validated Survival Model Performance Metrics
- List of Prognostically Revlevant Genes

### References
R package citations are provided at the end of the R script.
