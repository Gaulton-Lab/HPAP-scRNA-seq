#Note, perform this script after running SoupX to ensure use of counts corrected for ambient RNA
##Failure to account for ambient RNA can result in false positive identification of doublets

#Load necessary packages
suppressMessages(library(hdf5r))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(Matrix))
suppressMessages(library(harmony))
suppressMessages(library(data.table))
suppressMessages(library(ggpubr))
suppressMessages(library(future))
suppressMessages(library(tictoc))
suppressMessages(library(plyr))
suppressMessages(library(enrichR))
suppressMessages(library(stringr))
suppressMessages(library(ggrepel))
suppressMessages(library(reticulate)) #Must load reticulate to run Python scripts in R

########################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

samples <- c('HPAP-019','HPAP-020','HPAP-021','HPAP-022','HPAP-023','HPAP-024','HPAP-026','HPAP-028','HPAP-029','HPAP-032','HPAP-034','HPAP-035','HPAP-036','HPAP-037','HPAP-038','HPAP-039','HPAP-040','HPAP-042','HPAP-043','HPAP-044','HPAP-045','HPAP-047','HPAP-049','HPAP-050','HPAP-051','HPAP-052','HPAP-053','HPAP-054','HPAP-055','HPAP-056','HPAP-057','HPAP-058','HPAP-059','HPAP-061','HPAP-063','HPAP-064','HPAP-065','HPAP-070','HPAP-071','HPAP-072','HPAP-074','HPAP-075','HPAP-077','HPAP-079','HPAP-080','HPAP-081','HPAP-082','HPAP-083','HPAP-084','HPAP-085','HPAP-087','HPAP-088','HPAP-090','HPAP-091','HPAP-092','HPAP-099','HPAP-100','HPAP-101','HPAP-103','HPAP-104','HPAP-105','HPAP-106','HPAP-107','HPAP-108','HPAP-109')

#Loop through samples and output MatrixMarket files with count data
##Note that counts are now SoupX corrected
for (sample in samples){    
    #Read in post-SoupX Seurat objects
    data <- readRDS(file = paste0('~/hpap/SoupX/', sample, '_SoupX_roundToInt.rds'))
  
    DefaultAssay(data) <- 'RNA'
    rna_counts <- GetAssayData(adata,slot='counts')
    mm <- paste0('~/hpap/Scrublet/',sample,'_count_matrix_roundToInt.mtx')
    writeMM(rna_counts,mm)
    
    #Export gene list
    genes <- paste0('~/hpap/Scrublet/',sample,'_genes_roundToInt.tsv')
    write(row.names(rna_counts),genes,sep='\n')

    #Export barcode list
    barcodes <- paste0('~/hpap/Scrublet/',sample,'_barcodes_roundToInt.tsv')
    write(colnames(rna_counts),barcodes,sep='\n')
}

import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import time
from datetime import datetime

#########################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     ##
##              The following must be saved as a .py script before continuing (see example in repository). We will use reticulate to run the Python script in R.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                ##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     ##                                                      
#########################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

# samples = ['HPAP-019','HPAP-020','HPAP-021','HPAP-022','HPAP-023','HPAP-024','HPAP-026','HPAP-028','HPAP-029','HPAP-032','HPAP-034','HPAP-035','HPAP-036','HPAP-037','HPAP-038','HPAP-039','HPAP-040','HPAP-042','HPAP-043','HPAP-044','HPAP-045','HPAP-047','HPAP-049','HPAP-050','HPAP-051','HPAP-052','HPAP-053','HPAP-054','HPAP-055','HPAP-056','HPAP-057','HPAP-058','HPAP-059','HPAP-061','HPAP-063','HPAP-064','HPAP-065','HPAP-070','HPAP-071','HPAP-072','HPAP-074','HPAP-075','HPAP-077','HPAP-079','HPAP-080','HPAP-081','HPAP-082','HPAP-083','HPAP-084','HPAP-085','HPAP-087','HPAP-088','HPAP-090','HPAP-091','HPAP-092','HPAP-099','HPAP-100','HPAP-101','HPAP-103','HPAP-104','HPAP-105','HPAP-106','HPAP-107','HPAP-108','HPAP-109']
# for sample in samples:
#     #Read in sample counts and gene lists
#     input_dir = '~/hpap/Scrublet/' + sample
#     counts_matrix = scipy.io.mmread(input_dir + '_matrix_roundToInt.mtx').T.tocsc()
#     genes = np.array(scr.load_genes(input_dir + '_genes_roundToInt.tsv', delimiter='\t', column=0))
    
#     #Run Scrublet with default thresholds
#     scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)
#     doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
#                                                               min_cells=3, 
#                                                               min_gene_variability_pctl=85, 
#                                                               n_prin_comps=30)
    
#     #Pull out predicted doublets from Scrublet
#     predicted_doublets = scrub.call_doublets() 
    
#     #Save results
#     barcodes = input_dir + '_barcodes_roundToInt.tsv'
#     with open(barcodes, 'r') as f:
#         with open('~/hpap/Scrublet/scrublet_predicted_doublets_cutoff{}.txt'.format(scrub.threshold_), 'w') as o:
#             counter=0
#             for line in f:
#                 o.write('\t'.join((line.strip(), str(predicted_doublets[counter]), str(doublet_scores[counter]), '\n')))
#                 counter+=1
    
#     #Delete the Scrublet object to continue loop
#     del scrub

########################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

py_run_file('~/hpap/Scrublet/Scrublet_Python_Script.py')

#Next we want to combine the Scrublet doublet calls from all samples so we can remove them from our Seurat object
fin_scrub_df <- data.frame()

for (sample in samples){
    print(sample)
    #Extract sample Scrublet cutoff value
    outdir <- sprintf('~/hpap/Scrublet/')
    files <- list.files('~/hpap/Scrublet')
    file_path <- files[grep(sample,files)]
    cutoff <- str_split_fixed(file_path, 'cutoff', n=2)[2]
    cutoff <- gsub('.txt', '', cutoff)
    print(cutoff)
    
    #Read in Scrublet file
    scrub_fp <- file.path(outdir,sprintf('%sscrublet_predicted_doublets_cutoff%s.txt',sample,cutoff))
    scrub_df <- read.table(scrub_fp, sep='\t', header=FALSE)
    
    #Add on sample barcode prefix
    scrub_df$V1 <- paste(sample, '_', scrub_df$V1, sep='')
    
    fin_scrub_df <- rbind(fin_scrub_df,scrub_df)
}

write.table(fin_scrub_df[,c(1,2,3)],'~/hpap/Scrublet/hpap_combined_scrublet_predicted_doublets.txt', sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)

###It's time to remove doublets from our Seurat object!
#Read in post-SoupX merged Seurat object
soupx_merged_data <- readRDS('~/hpap/SoupX/hpap_SoupX.rds')

#Read in Scrublet doublets
scrub <- read.table('~/hpap/Scrublet/hpap_combined_scrublet_predicted_doublets.txt')
scrub_doublets <- scrub[scrub$V2=='True',1]

#Create vector of assignments with barcode names
cells <- Cells(soupx_merged_data)
scrublet_metadata <- rep('Singleton', length(cells))
names(scrublet_metadata) <- cells
scrublet_metadata[scrub_doublets] <- 'Doublet'
table(scrublet_metadata)
head(names(scrublet_metadata))

#Add doublet calls as metadata column
soupx_merged_data <- AddMetaData(soupx_merged_data, scrublet_metadata, col.name='Scrublet')

#Subset non-doublets (singletons)
soupx_merged_data <- subset(soupx_merged_data, subset=scrublet=='Singleton')

#Since we removed barcodes, it's important to re-run batch correction, dimensional reduction, and finding clusters
soupx_scrublet_merged <- soupx_merged_data
soupx_scrublet_merged <- NormalizeData(soupx_scrublet_merged, normalization.method = 'LogNormalize', scale.factor = 10000)
soupx_scrublet_merged <- FindVariableFeatures(soupx_scrublet_merged, selection.method = 'vst', nfeatures = 2000)
soupx_scrublet_merged <- ScaleData(soupx_scrublet_merged, verbose = FALSE) %>% 
    RunPCA(pc.genes = soupx_scrublet_merged@var.genes, npcs = 20, verbose = FALSE)

soupx_scrublet_merged <- RunHarmony(soupx_scrublet_merged,c('library', 'tissue_source', 'chemistry'), assay.use='RNA', plot_convergence = FALSE)

soupx_scrublet_merged <- soupx_scrublet_merged %>% 
    RunUMAP(reduction = 'harmony', dims = 1:20) %>% 
    FindNeighbors(reduction = 'harmony', dims = 1:20) %>% 
    FindClusters(algorithm=4,resolution = 0.5 ,method = "igraph") #Algorithm 4 is Leiden clustering

#Save post-Scrublet Seurat object
saveRDS(soupx_scrublet_merged, '~/hpap/hpap_SoupX_Scrublet.rds')

########################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
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
 [1] Azimuth_0.4.6             shinyBS_0.61.1           
 [3] SoupX_1.6.1               reticulate_1.27          
 [5] patchwork_1.1.2           forcats_0.5.2            
 [7] stringr_1.5.0             purrr_1.0.1              
 [9] readr_2.1.2               tidyr_1.3.0              
[11] tibble_3.1.8              tidyverse_1.3.2          
[13] future_1.30.0             ggpubr_0.5.0             
[15] data.table_1.14.6         harmony_0.1.1            
[17] Rcpp_1.0.10               Matrix_1.5-1             
[19] ggplot2_3.4.0             dplyr_1.0.10             
[21] EnsDb.Hsapiens.v86_2.99.0 ensembldb_2.18.4         
[23] AnnotationFilter_1.18.0   GenomicFeatures_1.46.5   
[25] AnnotationDbi_1.56.2      Biobase_2.54.0           
[27] GenomicRanges_1.46.1      GenomeInfoDb_1.30.1      
[29] IRanges_2.28.0            S4Vectors_0.32.4         
[31] BiocGenerics_0.40.0       Signac_1.7.0             
[33] SeuratObject_4.1.3        Seurat_4.3.0             
[35] hdf5r_1.3.8              
