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
suppressMessages(library(ggplot2))
suppressMessages(library(reticulate))
suppressMessages(library(SoupX))
suppressMessages(library(Azimuth))

###Set up a reticulate in conda environment to run Python packages in R
Sys.setenv(RETICULATE_PYTHON="~/.conda/envs/reticulate/bin/python")
reticulate::use_python("~/.conda/envs/reticulate/bin/python")
reticulate::use_condaenv("~/.conda/envs/reticulate")
reticulate::py_module_available(module='leidenalg')
reticulate::import('leidenalg')

###############################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

#Run CellRanger on all individual samples, aligning to GRCh38 genome 
#!/bin/bash
~/cellranger-6.0.1/cellranger count --id ${SAMPLE} --fastqs ~/hpap/rna/${SAMPLE}/Upenn_scRNAseq/fastq/ --sample ${SAMPLE} --transcriptome ~/refdata-gex-GRCh38-2020-A/ --localcores 24 --localmem 150 --disable-ui; done

#Create list of sample directories where CellRanger outputs are stored
samples <- c('HPAP-019','HPAP-020','HPAP-021','HPAP-022','HPAP-023','HPAP-024','HPAP-026','HPAP-027','HPAP-028','HPAP-029','HPAP-032','HPAP-034','HPAP-035','HPAP-036','HPAP-037','HPAP-038','HPAP-039','HPAP-040','HPAP-042','HPAP-043','HPAP-044','HPAP-045','HPAP-047','HPAP-049','HPAP-050','HPAP-051','HPAP-052','HPAP-053','HPAP-054','HPAP-055','HPAP-056','HPAP-057','HPAP-058','HPAP-059','HPAP-061','HPAP-063','HPAP-064','HPAP-065','HPAP-070','HPAP-071','HPAP-072','HPAP-074','HPAP-075','HPAP-077','HPAP-079','HPAP-080','HPAP-081','HPAP-082','HPAP-083','HPAP-084','HPAP-085','HPAP-087','HPAP-088','HPAP-090','HPAP-091','HPAP-092','HPAP-093','HPAP-099','HPAP-100','HPAP-101','HPAP-103','HPAP-104','HPAP-105','HPAP-106','HPAP-107','HPAP-108','HPAP-109')
samples2 <- c('HPAP-019_10xscRNA_71783','HPAP-020_10xscRNA_73876','HPAP-021_10xscRNA_75162','HPAP-022_10xscRNA_75751','HPAP-023_10xscRNA_75798','HPAP-024_10xscRNA_75799','HPAP-026_10xscRNA_77407','HPAP-027_10xscRNA_78239','HPAP-028_10xscRNA_78241','HPAP-029_10xscRNA_78242','HPAP-032_10xscRNA_79863','HPAP-034_10xscRNA_81423','HPAP-035_10xscRNA_81832','HPAP-036_10xscRNA_81833','HPAP-037_10xscRNA_82543','HPAP-038_10xscRNA_84764','HPAP-039_10xscRNA_86672','HPAP-040_10xscRNA_86673','HPAP-042_10xscRNA_86674','HPAP-043_10xscRNA_86675','HPAP-044_10xscRNA_86676','HPAP-045_10xscRNA_86677','HPAP-047_10xscRNA_97169','HPAP-049_10xscRNA_97170','HPAP-050_10xscRNA_97171','HPAP-051_10xscRNA_97173','HPAP-052_10xscRNA_98611','HPAP-053_10xscRNA_98612','HPAP-054_10xscRNA_98613','HPAP-055_10xscRNA_98614','HPAP-056_10xscRNA_98615','HPAP-057_10xscRNA_98616','HPAP-058_10xscRNA_98617','HPAP-059_10xscRNA_98618','HPAP-061_99868','HPAP063_99870','HPAP064_99871','HPAP065_99872','HPAP070_99875','HPAP071_103407','HPAP072_103408','HPAP074_103409','HPAP075_103410','HPAP077_103411','HPAP079_103412','HPAP-080_103413','HPAP081_103414','HPAP082_103415','HPAP083_103416','HPAP-084_105149','HPAP085_103417','HPAP-087_103418','HPAP088_105150','HPAP-090_105151','HPAP091_105152','HPAP-092_105153','HPAP-093_105154','HPAP-099_106838','HPAP-100_106839','HPAP-101_106840','HPAP-103_106841','HPAP-104_107485','HPAP-105_107486','HPAP-106_107487','HPAP-107_109359','HPAP-108_109360','HPAP-109_109361')

for (sample in samples){
    wd <- sprintf('~/hpap/rna/%s/Upenn_scRNAseq/cellranger/%s/outs/raw_feature_bc_matrix/', samples, samples2)
    }
#Read in the raw data for each sample
data <- list()
for (x in wd){
    sample_name <- str_split_fixed(x, "/", n=13)[9]
    hpap.data <- Read10X(data.dir = file.path(x))
    #Initialize the Seurat object with the raw (non-normalized data)
    hpap <- CreateSeuratObject(counts = hpap.data, project = "HPAP", min.cells = 0, min.features = 500)
    data[[sample_name]] <- hpap   
}
data

#Merge Seurat objects for individual samples
merged_data <- merge(data[[samples[[1]]]], y=data[samples[2:length(samples)]], add.cell.ids=samples, project='HPAP')
merged_data$library <- substr(rownames(merged_data@meta.data),1,8)
merged_data

#Add meta data to merged object
##Note that donor metadata files can be downloaded directly from PANC-DB with account registration
cond_t1d <- c('HPAP-020','HPAP-021','HPAP-023','HPAP-028','HPAP-032','HPAP-055','HPAP-064','HPAP-071','HPAP-084','HPAP-087')
cond_t2d <- c('HPAP-051','HPAP-057','HPAP-058','HPAP-061','HPAP-065','HPAP-070','HPAP-079','HPAP-081','HPAP-083','HPAP-085','HPAP-088','HPAP-090','HPAP-091','HPAP-100','HPAP-106','HPAP-108', 'HPAP-109')
aab <- c('HPAP-019', 'HPAP-024', 'HPAP-029', 'HPAP-038', 'HPAP-043', 'HPAP-045', 'HPAP-049', 'HPAP-050', 'HPAP-072', 'HPAP-092', 'HPAP-107')
sex_F <- c('HPAP-022','HPAP-027','HPAP-036', 'HPAP-037','HPAP-039','HPAP-044','HPAP-045','HPAP-050','HPAP-051','HPAP-053','HPAP-054', 'HPAP-056','HPAP-057','HPAP-058','HPAP-061','HPAP-063', 'HPAP-074', 'HPAP-079','HPAP-081','HPAP-085','HPAP-090','HPAP-091','HPAP-093','HPAP-099','HPAP-101','HPAP-103','HPAP-105', 'HPAP-109')
penn <- c('HPAP-022','HPAP-027','HPAP-034','HPAP-035','HPAP-037','HPAP-040','HPAP-047','HPAP-051','HPAP-052','HPAP-053','HPAP-054','HPAP-055','HPAP-056','HPAP-057','HPAP-059','HPAP-061','HPAP-063','HPAP-064','HPAP-074','HPAP-075','HPAP-077','HPAP-083','HPAP-085','HPAP-099','HPAP-103','HPAP-104','HPAP-106')
v2 <- c('HPAP-019','HPAP-020','HPAP-021','HPAP-022','HPAP-023','HPAP-024','HPAP-026','HPAP-027','HPAP-028','HPAP-029','HPAP-032','HPAP-034','HPAP-035','HPAP-036','HPAP-037')

merged_data@meta.data$sex[merged_data@meta.data$library %in% sex_F] <- 'F'
merged_data@meta.data$sex[!merged_data@meta.data$library %in% sex_F] <- 'M'

merged_data@meta.data$condition <- 'ND'
merged_data@meta.data$condition[merged_data@meta.data$library %in% cond_t1d] <- 'T1D'
merged_data@meta.data$condition[merged_data@meta.data$library %in% cond_t2d] <- 'T2D'
merged_data@meta.data$condition[merged_data@meta.data$library %in% aab] <- 'AAB+'

merged_data@meta.data$tissue_source <- 'nPod'
merged_data@meta.data$tissue_source[merged_data@meta.data$library %in% penn] <- 'UPenn'

merged_data@meta.data$chemistry <- '10Xv3'
merged_data@meta.data$chemistry[merged_data@meta.data$library %in% v2] <- '10Xv2'

merged_data[['percent.mt']] <- PercentageFeatureSet(merged_data, pattern = '^MT-')

#Density plot can be generated to help determine a percent mitochondrial read threshold to set
quants <- quantile(merged_data@meta.data[,'percent.mt'])
quants
p <- ggplot(merged_data@meta.data, aes(x=percent.mt)) + 
  geom_density() + theme_classic()
p + geom_vline(xintercept=c(quants[2],quants[3],quants[4]), colour=c('blue', 'red', 'black'),linetype = 'longdash')

#Visualize QC metrics *Note that HPAP-027 and HPAP-093 were removed from downstream analyses since mean nFeatures < 1000
options(repr.plot.height = 20, repr.plot.width = 30)
VlnPlot(merged_data, features = c('nFeature_RNA', 'nCount_RNA','percent.mt'), ncol = 1, group.by = 'library',
  log = TRUE, pt.size = 0) + NoLegend()

#Perform basic threshold filtering and log normalization
merged_data <- subset(merged_data, subset = percent.mt < 15) #Keeping every cell with less than 15 percent mitochondrial reads
merged_data <- NormalizeData(merged_data, normalization.method = "LogNormalize", scale.factor = 10000)
merged_data <- FindVariableFeatures(merged_data, selection.method = "vst", nfeatures = 2000)
merged_data <- ScaleData(merged_data, verbose = FALSE) %>% 
    RunPCA(pc.genes = merged_data@var.genes, npcs = 20, verbose = FALSE)

#Run Harmony batch correction with library, tissue source, and 10X kit chemistry as covariates
merged_data <- RunHarmony(merged_data,c('library','tissue_source','chemistry'), assay.use='RNA', plot_convergence = TRUE)

#Can now use harmony as a reduction in running UMAP and finding neighbors
merged_data <- merged_data %>% 
    RunUMAP(reduction = "harmony", dims = 1:20) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
    FindClusters(algorithm=4,resolution = 0.5 ,method = "igraph") #Algorithm 4 is Leiden clustering

#Plot UMAP
options(repr.plot.height = 7, repr.plot.width = 8)
DimPlot(merged_data, reduction = "umap", label = TRUE, pt.size = .1)

#Plot your favorite marker genes
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
p1 <- DotPlot(merged_data, assay='RNA', features=marker.genes, cluster.idents=TRUE) 
p1 <- p1 + theme(axis.text.x=element_text(angle=45, hjust=1)) + xlab('') + ylab('')
p1

#Congratulations- you've completed basic clustering of the HPAP dataset!

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
