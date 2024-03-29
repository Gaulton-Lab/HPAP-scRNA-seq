###This script uses automated SoupX contamination estimates. SoupX has other methods for estimating contamination.

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
suppressMessages(library(tidyverse))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(reticulate))
suppressMessages(library(SoupX))

###############################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
samples <- c('HPAP-019','HPAP-020','HPAP-021','HPAP-022','HPAP-023','HPAP-024','HPAP-026','HPAP-027','HPAP-028','HPAP-029','HPAP-032','HPAP-034','HPAP-035','HPAP-036','HPAP-037','HPAP-038','HPAP-039','HPAP-040','HPAP-042','HPAP-043','HPAP-044','HPAP-045','HPAP-047','HPAP-049','HPAP-050','HPAP-051','HPAP-052','HPAP-053','HPAP-054','HPAP-055','HPAP-056','HPAP-057','HPAP-058','HPAP-059','HPAP-061','HPAP-063','HPAP-064','HPAP-065','HPAP-070','HPAP-071','HPAP-072','HPAP-074','HPAP-075','HPAP-077','HPAP-079','HPAP-080','HPAP-081','HPAP-082','HPAP-083','HPAP-084','HPAP-085','HPAP-087','HPAP-088','HPAP-090','HPAP-091','HPAP-092','HPAP-093','HPAP-099','HPAP-100','HPAP-101','HPAP-103','HPAP-104','HPAP-105','HPAP-106','HPAP-107','HPAP-108','HPAP-109')

for (sample in samples){
    wd <- sprintf('~/hpap/rna/%s/Upenn_scRNAseq/cellranger/%s/outs', samples, samples)
    }

contam_frac_results0 <- NULL
contam_frac_results <- data.frame()
for (x in wd){
    name <- str_split_fixed(x, "/", n=12)[9] #splits in the path to output just the sample name
    rna_counts <- Read10X_h5(file.path(x, 'raw_feature_bc_matrix.h5'))
    data <- CreateSeuratObject(counts=rna_counts)
    data[['percent.mt']] <- PercentageFeatureSet(data, pattern = '^MT-')
    data <- subset(x = data, subset = nFeature_RNA > 500)
    
    data <- NormalizeData(data, normalization.method = 'LogNormalize', scale.factor = 10000)
    data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
    all.genes <- rownames(data)
    data <- ScaleData(data, features = all.genes)
    
    data <- RunPCA(data, verbose = FALSE)
    data <- RunUMAP(data, dims = 1:30, verbose = FALSE)
    data <- FindNeighbors(data, dims = 1:30, verbose = FALSE)
    data <- FindClusters(data, algorithm=4, resolution = 1, verbose=FALSE)
    
    DefaultAssay(data) <- 'RNA'
    toc <- GetAssayData(object = data, slot = 'counts') #with nFeature >500 filter
    tod <- Seurat::Read10X_h5(file.path(x, 'raw_feature_bc_matrix.h5'))
    
    metadata <- (cbind(as.data.frame(data[["umap"]]@cell.embeddings),
                   as.data.frame(Idents(data)),
                   as.data.frame(Idents(data))))
    colnames(metadata) <- c('RD1','RD2','Cluster','Annotation')
    
    sc <- SoupChannel(tod,toc)
    sc <- setDR(sc,metadata[colnames(sc$toc),c('RD1','RD2')])
    sc <- setClusters(sc,setNames(metadata$Cluster,rownames(metadata)))
    sc <- autoEstCont(sc)
    save <- sc$fit$dd[,c('gene', 'soupExp')]
    colnames(save) <- c('gene', name)
    write.table(save, paste0('~/hpap/SoupX/', name, '_gene_level_soup_exp.tsv'), sep ='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
    
    contamination_fraction <- mean(sc$metaData$rho*100)
    contam_frac_results0$Sample <- name
    contam_frac_results0$ContaminationFraction <- contamination_fraction
    contam_frac_results <- rbind(contam_frac_results, contam_frac_results0) #This will output a table with calculated contamination fraction for each sample
    
    out <- adjustCounts(sc, roundToInt=TRUE) #Stochastically adjusts counts while rounding to maintain overall contamination fraction while outputting integer counts

    data2 <- CreateSeuratObject(out)
    data2[['percent.mt']] <- PercentageFeatureSet(data2, pattern = '^MT-')
    data2 <- NormalizeData(data2, normalization.method = 'LogNormalize', scale.factor = 10000)
    data2 <- FindVariableFeatures(data2, selection.method = "vst", nfeatures = 2000)
    all.genes <- rownames(data2)
    data2 <- ScaleData(data2, features = all.genes)
    data2 <- RunPCA(data2, verbose = FALSE)
    data2 <- RunUMAP(data2, dims = 1:30, verbose = FALSE)
    data2 <- FindNeighbors(data2, dims = 1:30, verbose = FALSE)
    data2 <- FindClusters(data2, algorithm=4, resolution = 1, verbose=FALSE)
    saveRDS(data2, file = sprintf('~/hpap/SoupX/%s_SoupX.rds',name))
}

#Create a merged Seurat object from the individual sample post-SoupX Seurat objects
setwd('~/hpap/SoupX/')
soupx_files <- list.files('/~/hpap/SoupX/', pattern='_SoupX.rds')
soupx_data <- list()

for (x in soupx_files){
    sample_name <- str_split_fixed(x, '_', n=4)[2]
    tmp <- readRDS(x)
    soupx_data[[sample_name]] <- tmp
}

soupx_merged_data <- merge(soupx_data[[samples[[1]]]], y=soupx_data[samples[2:length(samples)]], add.cell.ids=samples, project='HPAP')
soupx_merged_data$library <- substr(rownames(soupx_merged_data@meta.data),1,8)

#Add meta data to merged object
##Note that donor metadata files can be downloaded directly from PANC-DB with account registration
cond_t1d <- c('HPAP-020','HPAP-021','HPAP-023','HPAP-028','HPAP-032','HPAP-055','HPAP-064','HPAP-071','HPAP-084','HPAP-087')
cond_t2d <- c('HPAP-051','HPAP-057','HPAP-058','HPAP-061','HPAP-065','HPAP-070','HPAP-079','HPAP-081','HPAP-083','HPAP-085','HPAP-088','HPAP-090','HPAP-091','HPAP-100','HPAP-106','HPAP-108', 'HPAP-109')
aab <- c('HPAP-019', 'HPAP-024', 'HPAP-029', 'HPAP-038', 'HPAP-043', 'HPAP-045', 'HPAP-049', 'HPAP-050', 'HPAP-072', 'HPAP-092', 'HPAP-107')
sex_F <- c('HPAP-022','HPAP-027','HPAP-036', 'HPAP-037','HPAP-039','HPAP-044','HPAP-045','HPAP-050','HPAP-051','HPAP-053','HPAP-054','HPAP-057','HPAP-058','HPAP-061','HPAP-063', 'HPAP-074', 'HPAP-079','HPAP-081','HPAP-085','HPAP-090','HPAP-091','HPAP-093','HPAP-099','HPAP-101','HPAP-103','HPAP-105', 'HPAP-109')
penn <- c('HPAP-022','HPAP-027','HPAP-034','HPAP-035','HPAP-037','HPAP-040','HPAP-047','HPAP-051','HPAP-052','HPAP-053','HPAP-054','HPAP-055','HPAP-056','HPAP-057','HPAP-059','HPAP-061','HPAP-063','HPAP-064','HPAP-074','HPAP-075','HPAP-077','HPAP-083','HPAP-085','HPAP-099','HPAP-103','HPAP-104','HPAP-106')
v2 <- c('HPAP-019','HPAP-020','HPAP-021','HPAP-022','HPAP-023','HPAP-024','HPAP-026','HPAP-027','HPAP-028','HPAP-029','HPAP-032','HPAP-034','HPAP-035','HPAP-036','HPAP-037')

soupx_merged_data@meta.data$sex[soupx_merged_data@meta.data$library %in% sex_F] <- 'F'
soupx_merged_data@meta.data$sex[!soupx_merged_data@meta.data$library %in% sex_F] <- 'M'
soupx_merged_data@meta.data$condition <- 'ND'
soupx_merged_data@meta.data$condition[soupx_merged_data@meta.data$library %in% cond_t1d] <- 'T1D'
soupx_merged_data@meta.data$condition[soupx_merged_data@meta.data$library %in% cond_t2d] <- 'T2D'

soupx_merged_data@meta.data$condition2 <- soupx_merged_data@meta.data$condition
soupx_merged_data@meta.data$condition2[soupx_merged_data@meta.data$library %in% aab] <- 'AAB+'

soupx_merged_data@meta.data$tissue_source <- 'nPod'
soupx_merged_data@meta.data$tissue_source[soupx_merged_data@meta.data$library %in% penn] <- 'UPenn'

soupx_merged_data@meta.data$chemistry <- '10Xv3'
soupx_merged_data@meta.data$chemistry[soupx_merged_data@meta.data$library %in% v2] <- '10Xv2'

soupx_merged_data[['percent.mt']] <- PercentageFeatureSet(soupx_merged_data, pattern = '^MT-')

#Set percent mitochondrial read threshold, log normalize data and find variable features
soupx_merged_data <- subset(soupx_merged_data, subset = percent.mt < 15) #Keeping everything less than 15 percent mitochondrial reads
soupx_merged_data <- NormalizeData(soupx_merged_data, normalization.method = 'LogNormalize', scale.factor = 10000)
soupx_merged_data <- FindVariableFeatures(soupx_merged_data, selection.method = 'vst', nfeatures = 2000)
soupx_merged_data <- ScaleData(soupx_merged_data, verbose = FALSE) %>% 
    RunPCA(pc.genes = soupx_merged_data@var.genes, npcs = 20, verbose = FALSE)

#Run Harmony batch correction for donor ID, tissue source, and 10x kit chemistry
soupx_merged_data <- RunHarmony(soupx_merged_data,c('library', 'tissue_source', 'chemistry'), assay.use='RNA', plot_convergence = TRUE)

#Dimensional reduction and finding clusters
soupx_merged_data <- soupx_merged_data %>% 
    RunUMAP(reduction = 'harmony', dims = 1:20) %>% 
    FindNeighbors(reduction = 'harmony', dims = 1:20) %>% 
    FindClusters(algorithm=4,resolution = 0.5 ,method = 'igraph') #Algorithm 4 is Leiden

saveRDS(soupx_merged_data, '~/hpap/SoupX/hpap_SoupX.rds')

####################################################################################################################################################################
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
 [1] SoupX_1.6.1               knitr_1.41               
 [3] forcats_0.5.2             stringr_1.5.0            
 [5] purrr_1.0.1               readr_2.1.2              
 [7] tidyr_1.3.0               tibble_3.1.8             
 [9] tidyverse_1.3.2           ggpubr_0.5.0             
[11] data.table_1.14.6         harmony_0.1.1            
[13] Rcpp_1.0.10               Matrix_1.5-1             
[15] ggplot2_3.4.0             dplyr_1.0.10             
[17] EnsDb.Hsapiens.v86_2.99.0 ensembldb_2.18.4         
[19] AnnotationFilter_1.18.0   GenomicFeatures_1.46.5   
[21] AnnotationDbi_1.56.2      Biobase_2.54.0           
[23] Signac_1.7.0              SeuratObject_4.1.3       
[25] Seurat_4.3.0              hdf5r_1.3.8              
[27] GenomicRanges_1.46.1      GenomeInfoDb_1.30.1      
[29] IRanges_2.28.0            S4Vectors_0.32.4         
[31] BiocGenerics_0.40.0      
