# HPAP scRNA-seq Analysis Pipeline
Here you can find all the code used to generate data and figures in our manuscript using publically available pancreatic islet scRNA-seq data from the Human Pancreas Analysis Program (HPAP).

Data and interactive web browsers can be found at www.isletgenomics.org

#### Elgamal RM et al. A reference map of cell type-specific gene expression in pancreatic islets in physiology and disease. 2023

<img src="https://github.com/Gaulton-Lab/HPAP-scRNA-seq/blob/56e78b18ccdc319b3d6767568526fe2e7bacf6f0/Images/CellxGene_UMAP.png" width="500" height="400" />

### Step 1: HPAP-Initial-Clustering.R 
Runs CellRanger on raw fastq files, basic filtering, batch correction and initial clustering

### Step 2: HPAP-SoupX.R 
Runs SoupX on individual samples to remove ambient RNA

### Step 3: HPAP-Scrublet.R 
Runs Scrublet to model and remove dooublets from individual samples

### Step 4: HPAP-Post-Processing.R 
Clean up the data and perform final clustering

### Step 5: HPAP-DESeq.R 
Generates pseudo-bulk matrices for each cluster to perform differential expression analyses

### Bonus: HPAP-Azimuth-Reference.R 
Instructions to format final Seurat object into an Azimuth reference
