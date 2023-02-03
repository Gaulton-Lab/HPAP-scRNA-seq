###This script provides instructions for converting your Seurat object to be compatible as an Azimuth reference
###See Azimuth reference format requirements here: https://github.com/satijalab/azimuth/wiki/Azimuth-Reference-Format

#Load necessary packages
suppressMessages(library(Azimuth))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(hdf5r))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(harmony))
suppressMessages(library(data.table))
suppressMessages(library(ggpubr))
suppressMessages(library(tidyverse))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(reticulate))

###############################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
#Load in final Seurat object
data <- readRDS('~/hpap/hpap_islet_scRNAseq.rds')

#Azimuth references are required to have sc transformed data (normalizes for library size)
data <- SCTransform(data, verbose = FALSE)
data <- RunPCA(data, verbose = FALSE)

#Run UMAP and find neighbors using SCT data. UMAP method must be  set to uwot and return.model = TRUE
data <- RunUMAP(data, dims = 1:20, assay = 'SCT', reduction = 'harmony', verbose = FALSE, umap.method = 'uwot', return.model = TRUE)
data <- FindNeighbors(data, dims = 1:20,assay = 'SCT', reduction = 'harmony', verbose = FALSE)

names(data)[names(data) == 'cell_type'] <- 'cell.type'

reference_data <- AzimuthReference(
 object = data,
 refUMAP = 'umap',
 refDR = 'pca',
 metadata = c('cell.type'), 
 dims = 1:50,
 k.param = 31)
 
saveRDS(reference_data, '~/hpap/ref.Rds')
SaveAnnoyIndex(reference_data[['refdr.annoy.neighbors']], '~/hpap/idx.annoy')

###############################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
sessionInfo()
R version 4.1.1 (2021-08-10)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.2 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
LAPACK: ~/.conda/envs/reticulate/lib/libmkl_rt.so.1

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] SoupX_1.6.1               reticulate_1.27          
 [3] patchwork_1.1.2           forcats_0.5.2            
 [5] stringr_1.5.0             purrr_1.0.1              
 [7] readr_2.1.2               tidyr_1.3.0              
 [9] tibble_3.1.8              tidyverse_1.3.2          
[11] future_1.30.0             ggpubr_0.5.0             
[13] data.table_1.14.6         harmony_0.1.1            
[15] Rcpp_1.0.10               Matrix_1.5-1             
[17] ggplot2_3.4.0             dplyr_1.0.10             
[19] EnsDb.Hsapiens.v86_2.99.0 ensembldb_2.18.4         
[21] AnnotationFilter_1.18.0   GenomicFeatures_1.46.5   
[23] AnnotationDbi_1.56.2      Biobase_2.54.0           
[25] GenomicRanges_1.46.1      GenomeInfoDb_1.30.1      
[27] IRanges_2.28.0            S4Vectors_0.32.4         
[29] BiocGenerics_0.40.0       Signac_1.7.0             
[31] SeuratObject_4.1.3        Seurat_4.3.0             
[33] hdf5r_1.3.8 
