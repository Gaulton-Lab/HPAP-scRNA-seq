###This script will go over how to process and cluster your final Seurat object after performing SoupX and Scrublet

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

###Set up a reticulate in conda environment to run Python packages in R
Sys.setenv(RETICULATE_PYTHON='~/.conda/envs/reticulate/bin/python')
reticulate::use_python('~/.conda/envs/reticulate/bin/python')
reticulate::use_condaenv('~/.conda/envs/reticulate')
reticulate::py_module_available(module='leidenalg')
reticulate::import('leidenalg')

###############################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
###Perform manual removal of remaining doublets, and subcluster if a single cluster has marker genes from two distinct cell types

#Read in post-SoupX, post-Scrublet Seurat object
data <- readRDS('~/hpap/hpap_SoupX_Scrublet.rds')

#Use the following function for subclustering
FindSubCluster1 <- function(
  object,
  cluster,
  graph.name,
  subcluster.name = 'sub.cluster',
  resolution = 0.5,
  algorithm = 1
) {
  sub.cell <- WhichCells(object = object, idents = cluster)
  sub.graph <- as.Graph(x = object[[graph.name]][sub.cell, sub.cell])
  sub.clusters <- FindClusters(
    object = sub.graph,
    resolution = resolution,
    algorithm = algorithm,
    method = 'igraph'
  )
  sub.clusters[, 1] <- paste(cluster,  sub.clusters[, 1], sep = '_')
  object[[subcluster.name]] <- as.character(x = Idents(object = object))
  object[[subcluster.name]][sub.cell, ] <- sub.clusters[, 1]
  return(object)
}

Idents(data) <- data@meta.data$seurat_clusters

#Pick a cluster to sub-cluster
sc.num <- 9

data_subcluster <- FindSubCluster1(data, cluster=sc.num, algorithm=4, resolution=0.25, graph.name = 'RNA_snn')  
data_subset <- subset(data_subcluster, subset=seurat_clusters==sc.num)

options(repr.plot.width=10, repr.plot.height=8)
DimPlot(data_subset, reduction='umap', group.by='sub.cluster', label=TRUE, label.size=6, repel=TRUE)


#Pick marker genes that you are trying to separate in your subclustering
marker.genes <- c('INS', 'GCG')

options(repr.plot.width=25, repr.plot.height=12)
DotPlot(data_subset, assay='RNA', group.by='sub.cluster', features=marker.genes, cluster.idents=TRUE) + theme(axis.text.x=element_text(angle=45, hjust=1)) + xlab('') + ylab('')

#If you are happy with subclustering, add to object metadata
###For first subcluster:
data@meta.data$seurat_clusters <- data_subcluster@meta.data$sub.cluster

###For remaining subclusters:
for (x in 1:nrow(data_subcluster@meta.data)){
    if (data_subcluster@meta.data[x,]$seurat_clusters ==  9){
        data@meta.data[x,]$seurat_clusters <-  data_subcluster@meta.data[x,]$sub.cluster
    }    
}

#Manually remove low quality clusters or doublets
subclusters_to_remove <- c('3_4','4_2','7_1','7_4', '9_3', '9_4', '10_3','13','14','16_1','16_2','19','20')
barcodes_to_remove <- data[[]][data$seurat_clusters %in% subclusters_to_remove, c()]
barcodes_to_remove <- rownames(barcodes_to_remove)

data$remove_cells <- (Cells(data) %in% barcodes_to_remove)
data <- subset(data, subset=remove_cells==FALSE)

#Since we removed barcodes, it's important to re-run batch correction, dimensional reduction, and finding clusters
data <- RunHarmony(data,c('library', 'chemistry', 'tissue_source'), assay.use='RNA', plot_convergence = FALSE)

data <- data %>% 
    RunUMAP(reduction = 'harmony', dims = 1:20) %>% 
    FindNeighbors(reduction = 'harmony', dims = 1:20) %>% 
    FindClusters(algorithm=4,resolution = 0.30 ,method = 'igraph') #Algorithm 4 is Leiden clustering

#########At this point we can perform subclustering again but only for annotation, not to remove any barcodes#########
#Pancreatic islet cell types have well established marker genes. Add metadata column for cell type annotations based on marker gene expression

data@meta.data$cell_type <-  NA

data@meta.data[which(data@meta.data$seurat_clusters == 1),]$cell_type <- 'Alpha_1'
data@meta.data[which(data@meta.data$seurat_clusters == 2),]$cell_type <- 'Acinar'
data@meta.data[which(data@meta.data$seurat_clusters == 3),]$cell_type <- 'Beta_1'
data@meta.data[which(data@meta.data$seurat_clusters == 4),]$cell_type <- 'Ductal'
data@meta.data[which(data@meta.data$seurat_clusters == 5),]$cell_type <- 'Active_Stellate'
data@meta.data[which(data@meta.data$seurat_clusters == 6),]$cell_type <- 'Endothelial'
data@meta.data[which(data@meta.data$seurat_clusters == '7_1'),]$cell_type <- 'Delta'
data@meta.data[which(data@meta.data$seurat_clusters == '7_2'),]$cell_type <- 'Gamma+Epsilon'
data@meta.data[which(data@meta.data$seurat_clusters == '8_1'),]$cell_type <- 'Beta_2'
data@meta.data[which(data@meta.data$seurat_clusters == '8_2'),]$cell_type <- 'Alpha_2' 
data@meta.data[which(data@meta.data$seurat_clusters == 9),]$cell_type <- 'Quiescent_Stellate' 
data@meta.data[which(data@meta.data$seurat_clusters == 10),]$cell_type <- 'Macrophage' 
data@meta.data[which(data@meta.data$seurat_clusters == 11),]$cell_type <- 'MUC5B+_Ductal'
data@meta.data[which(data@meta.data$seurat_clusters == 12),]$cell_type <- 'Mast' 
data@meta.data[which(data@meta.data$seurat_clusters == 13),]$cell_type <- 'Cycling_Alpha'

###You made it! Plot your final UMAP and favorite marker genes.
Idents(data) <- data@meta.data$cell_type
options(repr.plot.height = 7, repr.plot.width = 9)
DimPlot(data, reduction = 'umap', label = FALSE, pt.size = 0.01, raster=FALSE)

marker.genes <- c('INS','IAPP', #Beta
                  'GCG', #Alpha
                  'SST', #Delta
                  'PPY', #Gamma
                  'GHRL', #Epsilon
                  'CFTR','MUC5B', #Ductal
                  'REG1A','CTRB2','PRSS1', 'PRSS2', #Acinar
                  'PDGFRB','SPARC','COL6A1','RGS5', #Stellate
                  'PLVAP','ESAM','VWF', #Endothelial
                  'KIT', #Mast
                  'CD69','C1QB', 'C1QC', 'C1QA', #Immune (Macrophages),
                  'MKI67','CDK1') #Dividing cells
options(repr.plot.width=12, repr.plot.height=6)
DotPlot(data, assay='RNA', features=marker.genes, cluster.idents=TRUE) + theme(axis.text.x=element_text(angle=45, hjust=1)) + xlab('') + ylab('')

#Be sure to save your final Seurat object!
saveRDS(data, '~/hpap/hpap_islet_scRNAseq.rds')

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
