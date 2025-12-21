# ANAT40040_Final_Assignment
## ERBB2 Amplification and Transcrtiptomic Analysis in Breast Cancer

### Overview

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
- This amplification status is the encoded as a binary metadata variable
#### 3. Differential Gene Expression Analysis
#### 4. Exploratory Data Analysis
#### 5. Pathway Enrichment Analysis
#### 6. Survival Modelling


### R Packages Used

### Output
