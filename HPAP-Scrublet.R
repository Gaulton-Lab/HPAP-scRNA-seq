#Note, perform this script after running SoupX to ensure use of counts corrected for ambient RNA
##Failure to account for ambient RNA can result in false positive identification of doublets

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

###Set up a reticulate in conda environment to run Python packages in R
Sys.setenv(RETICULATE_PYTHON="~/.conda/envs/reticulate/bin/python")
reticulate::use_python("~/.conda/envs/reticulate/bin/python")
reticulate::use_condaenv("~/.conda/envs/reticulate")
reticulate::py_module_available(module='leidenalg')
reticulate::import('leidenalg') 

samples <- c('HPAP-019','HPAP-020','HPAP-021','HPAP-022','HPAP-023','HPAP-024','HPAP-026','HPAP-028','HPAP-029','HPAP-032','HPAP-034','HPAP-035','HPAP-036','HPAP-037','HPAP-038','HPAP-039','HPAP-040','HPAP-042','HPAP-043','HPAP-044','HPAP-045','HPAP-047','HPAP-049','HPAP-050','HPAP-051','HPAP-052','HPAP-053','HPAP-054','HPAP-055','HPAP-056','HPAP-057','HPAP-058','HPAP-059','HPAP-061','HPAP-063','HPAP-064','HPAP-065','HPAP-070','HPAP-071','HPAP-072','HPAP-074','HPAP-075','HPAP-077','HPAP-079','HPAP-080','HPAP-081','HPAP-082','HPAP-083','HPAP-084','HPAP-085','HPAP-087','HPAP-088','HPAP-090','HPAP-091','HPAP-092','HPAP-099','HPAP-100','HPAP-101','HPAP-103','HPAP-104','HPAP-105','HPAP-106','HPAP-107','HPAP-108','HPAP-109')

#Loop through samples and output MatrixMarket files with count data
##Note that counts are now SoupX corrected
for (sample in samples){    
    #Read in post-SoupX Seurat objects
    data <- readRDS(file = paste0("~/hpap/SoupX/", sample, "_SoupX_roundToInt.rds"))
  
    DefaultAssay(data) <- 'RNA'
    rna_counts <- GetAssayData(adata,slot='counts')
    mm <- paste0("~/hpap/Scrublet/",sample,'_count_matrix_roundToInt.mtx')
    writeMM(rna_counts,mm)
    
    #Export gene list
    genes <- paste0("~/hpap/Scrublet/",sample,'_genes_roundToInt.tsv')
    write(row.names(rna_counts),genes,sep='\n')

    #Export barcode list
    barcodes <- paste0("~/hpap/Scrublet/",sample,'_barcodes_roundToInt.tsv')
    write(colnames(rna_counts),barcodes,sep='\n')
}
