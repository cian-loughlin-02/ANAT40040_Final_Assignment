########## WILL TIDY THE CODE BEFORE CREATING GitHub REPOSITORY#########
########################################################################

#READING IN CLINICAL PATIENT DATA
data_patient  = read.delim("brca_tcga_pan_can_atlas_2018/data_clinical_patient.txt", skip = 4)
data_patient
#data_patient = data_patient[5:dim(data_patient)[1],]

#VISUALISING "AGE"
col_age = which(colnames(data_patient)=='AGE')
par(mfrow = c(1, 3))
hist(as.numeric(data_patient[,col_age]), main = "Histogram of Age at Diagnosis", breaks = 50, xlab = "Age at Diagnosis")

#VISUALISING 'STAGE'
col_stage = which(colnames(data_patient) == "AJCC_PATHOLOGIC_TUMOR_STAGE")
stage = data_patient[,col_stage]

stage = gsub("I[^VI]", "I", stage)
stage = gsub("V[^\\s]", "V", stage)
barplot(table(as.factor(stage)), las = 2, cex.names = 1, main = "Plot of Disease Stages" )


#VISUALISING OVERALL SURVIVAL TIME
col_overall_surv = which(colnames(data_patient) == "OS_STATUS")

uncensored = which(data_patient[,col_overall_surv] == "1:DECEASED")

col_overall_survival_time = which(colnames(data_patient) == "OS_MONTHS")
hist(as.numeric(data_patient[uncensored,col_overall_survival_time]), breaks = 50, main = "Histogram of Overall survival Months", xlab = "Survival [Months]")

#READING IN RNA SEQUENCE DATA
data_rna_seq  = read.delim("brca_tcga_pan_can_atlas_2018/data_mrna_seq_V2_rsem.txt")

#READING IN CNA DATA
data_cna  = read.delim("brca_tcga_pan_can_atlas_2018/data_cna.txt")
data_cna

#GETTING RNA ID's
assay = round(as.matrix(data_rna_seq[,-c(1,2)])) #columns 1 and 2 are the gene names
#assay
rownames(assay) = data_rna_seq[,1] # setting gene names as the row names
#rownames(assay)
n_samples = dim(assay)[2]
#n_samples
rna_id = rep("", n_samples)

for (i in 1:n_samples){
  patient_barcode = colnames(assay)[i]
  patient_barcode = substr(patient_barcode, 1, 12)
  patient_barcode = gsub("\\.", "-",patient_barcode)
  rna_id[i] = patient_barcode
}

#GETTING CNA ID's

cna_matrix = as.matrix(data_cna[,-1]) 
#cna_matrix
n_samples_cna = dim(cna_matrix)[2]
#n_samples_cna
cna_id = rep("", n_samples_cna)
#cna_id
for (i in 1:n_samples_cna){
  patient_barcode = colnames(cna_matrix)[i]
  patient_barcode = substr(patient_barcode, 1, 12)
  patient_barcode = gsub("\\.", "-",patient_barcode)
  cna_id[i] = patient_barcode
}
cna_id

# GETTING CLINICAL PATIENTS ID's

# Print the first column name
colnames(data_patient)[1]

pat_ids = data_patient[,1] #patient id is the first col on data_patient
n_samples_pat = length(pat_ids)
#n_samples_pat
clean_pat_ids = rep("", n_samples_pat)
#clean_pat_ids
for (i in 1:n_samples_pat){
  patient_barcode = pat_ids[i]
  patient_barcode = substr(patient_barcode, 1, 12)
  patient_barcode = gsub("\\.", "-",patient_barcode)
  clean_pat_ids[i] = patient_barcode
}
clean_pat_ids

# MATCHING 3 DATASETS
com_ids = intersect(intersect(rna_id, cna_id), clean_pat_ids)
length(com_ids)
length(rna_id)
length(cna_id)
length(clean_pat_ids) # this prints the amount of (number of) clean_pat_ids

#Creating subset of RNA-Seq assay for common ID's
assay = assay[, match(com_ids, rna_id)]
cna_matrix <- cna_matrix[, match(com_ids, cna_id)]

# CREATING METADATA FOR NORMALISATION

# Build metadata.

n_samples_com = length(com_ids) 
ERBB2_metadata = rep(0, n_samples_com) # creating a metadata matrix for ERBB2
idx_ERBB2 = which(data_cna[,1] == "ERBB2")
idx_ERBB2
for (i in 1:n_samples_com){
  pat_id = com_ids[i]
  idx_cna = which(cna_id == pat_id) # finding the matching cna cols
  cna_val = data_cna[idx_ERBB2, idx_cna +1] # obtaining cna value
  ERBB2_metadata[i] = 1*(as.numeric(cna_val) > 0) # assigning vals as amplified if cna val > 0
}

ERBB2_metadata[is.na(ERBB2_metadata)] = 0 #changes any na value to 0
sum(is.na(ERBB2_metadata))
metadata[is.na(ERBB2_metadata), ]
metadata = data.frame(ERBB2_amplified = ERBB2_metadata, row.names = com_ids) #Converting to metadata for DESeq2 use.
rownames(metadata) = com_ids

#metadata <- metadata[match(com_ids, rownames(metadata)), , drop=FALSE]
rownames(metadata) = colnames(assay)
metadata$ERBB2_amplified = as.factor(metadata$ERBB2_amplified)



# STARTING NORMALISATION
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.22") #had to specify version as I was having issues installing

BiocManager::install("DESeq2")
library(DESeq2)

assay = round(assay)
assay[is.na(assay)] = 0  # Impute 0 for the na values
assay[assay<0] = 0


smallestGroupSize = 3 
keep = rowSums(assay >= 10) >= smallestGroupSize # filter out genes with too many missing values.
assay = assay[keep,]

sum(is.na(metadata$ERBB2_amplified))##
metadata[is.na(metadata$ERBB2_amplified), ]

dds =  DESeqDataSetFromMatrix(countData = assay,
                              colData = metadata,
                              design = ~ ERBB2_amplified)
# Normalising
dds <- DESeq(dds)


resultsNames(dds) # this will list all of the coefficients

res_dds = results(dds)

# print Top 10 most differentially expressed

res_dds[order(res_dds$padj)[1:10],] # this prints the Top 10 most differentially expressed

#The top 10 genes by uing log2foldchange
topgenes <- res_dds[order(abs(res_dds$log2FoldChange), decreasing = TRUE), ] # ordering the top 10 genes
top10genes <- head(topgenes, 10)

#Ordering by log2FoldChange (largest positive at top)
topgenes_upreg <- res_dds[order(res_dds$log2FoldChange, decreasing = TRUE), ]
top10_upreg <- head(topgenes_upreg, 10)

#Ordering by log2FoldChange (largest negative at top)
topgenes_downreg <- res_dds[order(res_dds$log2FoldChange), ]
top10_downreg <- head(topgenes_downreg,10 )

print(top10_upreg)
print(top10_downreg)
print(top10genes)


head(res_dds, 10)

# Performing PCA
vst = vst(dds) #vst = variance stabilising transformation
vst_matrix <- assay(vst)

colnames(vst_matrix) <- gsub("\\.", "-", colnames(vst_matrix))
colnames(vst_matrix) <- substr(colnames(vst_matrix), 1, 12)

# Plotting PCA
plotPCA(vst, intgroup = "ERBB2_amplified")



#Instaling clusterProfile, org.Hs.eg.db & enrichplot

if (!requireNamespace("clusterProfiler", quietly = TRUE))
  BiocManager::install("clusterProfiler")

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
  BiocManager::install("org.Hs.eg.db")

if (!requireNamespace("enrichplot", quietly = TRUE))
  BiocManager::install("enrichplot")
  #install.packages("enrichplot")

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)


# Creating a subset of the differentially expressed genes.

res_sig = res_dds[res_dds$padj<0.05,]

# Separating into over and under expressed using log2foldchange

Diff_Exp_over = rownames(res_sig[res_sig$log2FoldChange>0,])
Diff_Exp_under = rownames(res_sig[res_sig$log2FoldChange<0,])

nrow(res_sig)
length(Diff_Exp_over)
length(Diff_Exp_under)

# Gene Ontology (GO) results for the over expressed
go_results_over = enrichGO(
  gene          = Diff_Exp_over,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",  
  ont           = "BP", 
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

# Printing and Plotting Results of Gene Ontology Enrichment OverExpressed

print(head(go_results_over)) 
dotplot(go_results_over, showCategory=10) + ggtitle("Gene Ontology Enrichment Over Expressed")

# Gene Ontology (GO) results for the under expressed
go_results_under = enrichGO(
  gene          = Diff_Exp_under,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",  
  ont           = "BP", 
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

# print and plot results
print(head(go_results_under)) 
dotplot(go_results_under, showCategory=10) + ggtitle("Gene Ontology Enrichment Under Expressed")


#Installing pathview and ReactomePA
if (!requireNamespace("pathview", quietly = TRUE))
  BiocManager::install("pathview")

if (!requireNamespace("ReactomePA", quietly = TRUE))
  BiocManager::install("ReactomePA")

library(ReactomePA)
library(pathview)

# Mapping into entrez for Reactome and Keggs
# Mapping for Over expressed
gene_entrez_over <- bitr(
  Diff_Exp_over,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

failed_over <- setdiff(Diff_Exp_over, gene_entrez_over$SYMBOL)
failed_over[1:20] # to see which genes failed to map

# Mapping for Under Expressed
gene_entrez_under <- bitr(
  Diff_Exp_under,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

failed_under <- setdiff(Diff_Exp_under, gene_entrez_under$SYMBOL)
failed_under[1:20]

#Kegg Results & Plots for both over and under expressed

kegg_res_over =  enrichKEGG(
  gene = gene_entrez_over$ENTREZID,
  organism = "human",   
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)


kegg_res_under =  enrichKEGG(
  gene          = gene_entrez_under[,2],
  organism      = "human",   
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

print(head(kegg_res_over))
dotplot(kegg_res_over, showCategory=10) + ggtitle("Kegg Pathway Enrichment Over Expressed")

print(head(kegg_res_under))
dotplot(kegg_res_under, showCategory=10) + ggtitle("Kegg Pathway Enrichment Under Expressed")

# Reactome Results & Plots for both over and under expressed
reactome_res_over =  enrichPathway(
  gene          = gene_entrez_over[,2],
  organism      = "human",   
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
)

reactome_res_under =  enrichPathway(
  gene          = gene_entrez_under[,2],
  organism      = "human",   
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
)


print(head(reactome_res_over))
dotplot(reactome_res_over, showCategory=10) + ggtitle("Reactome Pathway Enrichment Over Expressed")

print(head(reactome_res_under))
dotplot(reactome_res_under, showCategory=10) + ggtitle("Reactome Pathway Enrichment Under Expressed")


#Performing PCA
vst = vst(dds) #vst = variance stabilising transformation
vst_matrix <- assay(vst)

colnames(vst_matrix) <- gsub("\\.", "-", colnames(vst_matrix))
colnames(vst_matrix) <- substr(colnames(vst_matrix), 1, 12)
#Plotting PCA
plotPCA(vst, intgroup = "ERBB2_amplified")

#INSTALLING PHEATMAP
if (!requireNamespace("pheatmap", quietly = TRUE))
  install.packages("pheatmap")


library(pheatmap)

# Creating subset of data on DEG (Differentially Expressed Gene)
top_Diff_Exp = order(res_dds$padj)
vst_Diff_Exp = assay(vst)[top_Diff_Exp[1:20],]

# Creating the Heatmap
pheatmap(
  vst_Diff_Exp,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = 'row',
  show_colnames = FALSE,
  show_rownames = TRUE,
  annotation_col = as.data.frame(colData(vst)["ERBB2_amplified"])
)

# GENERATING AN OVERALL SURVIVAL MODEL
install.packages("glmnet")
library(glmnet)
library(survival)

clin_com <- data_patient[match(com_ids, clean_pat_ids), ] #finding OS_T for common patients

os_time <- as.numeric(clin_com$OS_MONTHS) #Finding Survival Time in months


os_stat <- ifelse(clin_com$OS_STATUS %in% c("1:DECEASED"), 1,
                  ifelse(clin_com$OS_STATUS %in% c("0:LIVING"), 0, NA))


valid_idx <- which(os_time > 0) #OS has to be positive 
os_time <- os_time[valid_idx]
os_stat <- os_stat[valid_idx]
 
y <- Surv(time = os_time, event =  os_stat) #Creating a survival object

vst_Diff_Exp <- vst_matrix[rownames(res_sig), match(com_ids, colnames(vst_matrix))] #vst vals for Diff Exp genes

x <- t(vst_Diff_Exp) #makes sure rows and genes are the cols (transpose)

x <- x[valid_idx, ] #this matches x to y


#Creating a Plot
set.seed(1)
cvfit <- cv.glmnet(x, y, family = "cox", type.measure = "C") #using "cox"
plot(cvfit)

sum(coef(cvfit, s = "lambda.min") != 0) #number non-zero coeff at lamda min
sum(coef(cvfit, s = "lambda.1se") != 0) #number non-zero coeff at lamda 1se (standard error)

peak_c_idx <- max(cvfit$cvm) #the peak C-index
peak_c_idx

peak_lambda <- cvfit$lambda[which.max(cvfit$cvm)] #finds lambda val at peak
peak_lambda

cvfit$lambda.min
cvfit$lambda.1se

-log(peak_lambda)

#Finding genes at lambda min and at 1 standard error

min_genes <- rownames(coef(cvfit, s = "lambda.min"))[as.numeric(coef(cvfit, s = "lambda.min")) != 0]
min_genes

genes_1se <- rownames(coef(cvfit, s = "lambda.1se"))[as.numeric(coef(cvfit, s = "lambda.1se")) != 0]
genes_1se

# References / Citations
citation("BiocManager")
citation("DESeq2")
citation("clusterProfiler")
print(toBibtex(citation("clusterProfiler")))
print(toBibtex(citation("org.Hs.eg.db")))
print(toBibtex(citation("enrichplot")))
print(toBibtex(citation("ReactomePA")))
print(toBibtex(citation("pathview")))
print(toBibtex(citation("pheatmap")))


