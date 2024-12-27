# Assignment 2
# Reference: Course content and practicals uploaded on brightspace.

### Step 0 - Load packages
library(stringr)

### Step 1 - Download the dataset from:
# https://www.cbioportal.org/study/summary?id=brca_tcga_pan_can_atlas_2018

### Step 2 - Untar the folder and extract the files
path= "C:\\Users\\Irum\\Documents\\MaliM\\UCD\\ANAT40040\\BIOa2\\assign2"

directory = "brca_tcga_pan_can_atlas_2018.tar.gz"
dir = paste(path, directory, sep = "/")
untar(dir)

# Set a new working directory
new_dir = paste(getwd(),"brca_tcga_pan_can_atlas_2018", sep = "/")
setwd(new_dir)

### Step 3 - Read RNA-seq file.
mrna_data = read.delim("data_mrna_seq_v2_rsem.txt")

### Step 4 - Read patient Data file.
patient_data = read.delim("data_clinical_patient.txt")

### Step 5 - Read CNA file.
cna_data = read.delim("data_cna.txt")

### Step 6 & 7 - Match patient IDs and build Metadata
assay = as.matrix(mrna_data[,-c(1,2)])
metadata = matrix(0, dim(assay)[2],1)
patientIDs = patient_data[,1]

# Used ctrl-f in the data_cna.txt file to locate which row is the ERBB2 row.
rowERBB2 = 19403

for (i in 1:dim(assay)[2]){
  pat_barcode = colnames(assay)[i]
  idx = which(pat_barcode == colnames(cna_data))
  if (length(idx) > 0) {
    metadata[i, 1] = 1*(as.numeric(cna_data[rowERBB2,idx])>0)
  }
 
}
metadata[is.na(metadata)] = 0
colnames(metadata) = "Early"

### Step 8 - Normalize with Dseq2
# First install BiocManager
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Next install DESeq2
if (!requireNamespace("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2")
library(DESeq2)

# Now we will normalize the data
assay[is.na(assay)] = 0
assay[assay<=0] = 0
dds <- DESeqDataSetFromMatrix(countData = round(assay),
                              colData = metadata,
                              design = ~ Early)
dds <- DESeq(dds)

### Step 9 - Differential gene expression analysis
resultsNames(dds)
res = results(dds)
res

res[order(res$padj)[1:10],]







