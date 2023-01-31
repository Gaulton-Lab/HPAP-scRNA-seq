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

#This block uses automated contamination estimates
samples <- c('HPAP-019','HPAP-020','HPAP-021','HPAP-022','HPAP-023','HPAP-024','HPAP-026','HPAP-027','HPAP-028','HPAP-029','HPAP-032','HPAP-034','HPAP-035','HPAP-036','HPAP-037','HPAP-038','HPAP-039','HPAP-040','HPAP-042','HPAP-043','HPAP-044','HPAP-045','HPAP-047','HPAP-049','HPAP-050','HPAP-051','HPAP-052','HPAP-053','HPAP-054','HPAP-055','HPAP-056','HPAP-057','HPAP-058','HPAP-059','HPAP-061','HPAP-063','HPAP-064','HPAP-065','HPAP-070','HPAP-071','HPAP-072','HPAP-074','HPAP-075','HPAP-077','HPAP-079','HPAP-080','HPAP-081','HPAP-082','HPAP-083','HPAP-084','HPAP-085','HPAP-087','HPAP-088','HPAP-090','HPAP-091','HPAP-092','HPAP-093','HPAP-099','HPAP-100','HPAP-101','HPAP-103','HPAP-104','HPAP-105','HPAP-106','HPAP-107','HPAP-108','HPAP-109')
samples2 <- c('HPAP-019_10xscRNA_71783','HPAP-020_10xscRNA_73876','HPAP-021_10xscRNA_75162','HPAP-022_10xscRNA_75751','HPAP-023_10xscRNA_75798','HPAP-024_10xscRNA_75799','HPAP-026_10xscRNA_77407','HPAP-027_10xscRNA_78239','HPAP-028_10xscRNA_78241','HPAP-029_10xscRNA_78242','HPAP-032_10xscRNA_79863','HPAP-034_10xscRNA_81423','HPAP-035_10xscRNA_81832','HPAP-036_10xscRNA_81833','HPAP-037_10xscRNA_82543','HPAP-038_10xscRNA_84764','HPAP-039_10xscRNA_86672','HPAP-040_10xscRNA_86673','HPAP-042_10xscRNA_86674','HPAP-043_10xscRNA_86675','HPAP-044_10xscRNA_86676','HPAP-045_10xscRNA_86677','HPAP-047_10xscRNA_97169','HPAP-049_10xscRNA_97170','HPAP-050_10xscRNA_97171','HPAP-051_10xscRNA_97173','HPAP-052_10xscRNA_98611','HPAP-053_10xscRNA_98612','HPAP-054_10xscRNA_98613','HPAP-055_10xscRNA_98614','HPAP-056_10xscRNA_98615','HPAP-057_10xscRNA_98616','HPAP-058_10xscRNA_98617','HPAP-059_10xscRNA_98618','HPAP-061_99868','HPAP063_99870','HPAP064_99871','HPAP065_99872','HPAP070_99875','HPAP071_103407','HPAP072_103408','HPAP074_103409','HPAP075_103410','HPAP077_103411','HPAP079_103412','HPAP-080_103413','HPAP081_103414','HPAP082_103415','HPAP083_103416','HPAP-084_105149','HPAP085_103417','HPAP-087_103418','HPAP088_105150','HPAP-090_105151','HPAP091_105152','HPAP-092_105153','HPAP-093_105154','HPAP-099_106838','HPAP-100_106839','HPAP-101_106840','HPAP-103_106841','HPAP-104_107485','HPAP-105_107486','HPAP-106_107487','HPAP-107_109359','HPAP-108_109360','HPAP-109_109361')

for (sample in samples){
    wd <- sprintf('~/hpap/rna/%s/Upenn_scRNAseq/cellranger/%s/outs', samples, samples2)
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
    
    out <- adjustCounts(sc, roundToInt=TRUE) #Stochastically adjusts counts while rounding to maintain overally contamination fraction while outputting integer counts

    data2 <- CreateSeuratObject(out)
    data2[['percent.mt']] <- PercentageFeatureSet(data2, pattern = '^MT-')
    data2 <- NormalizeData(data2, normalization.method = "LogNormalize", scale.factor = 10000)
    data2 <- FindVariableFeatures(data2, selection.method = "vst", nfeatures = 2000)
    all.genes <- rownames(data2)
    data2 <- ScaleData(data2, features = all.genes)
    data2 <- RunPCA(data2, verbose = FALSE)
    data2 <- RunUMAP(data2, dims = 1:30, verbose = FALSE)
    data2 <- FindNeighbors(data2, dims = 1:30, verbose = FALSE)
    data2 <- FindClusters(data2, algorithm=4, resolution = 1, verbose=FALSE)
    saveRDS(data2, file = sprintf("~/hpap/SoupX/soupX_%s_automated_roundToInt.rds",name))
    }

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

loaded via a namespace (and not attached):
  [1] utf8_1.2.2                  spatstat.explore_3.0-5     
  [3] reticulate_1.27             tidyselect_1.2.0           
  [5] RSQLite_2.2.19              htmlwidgets_1.6.1          
  [7] grid_4.1.1                  BiocParallel_1.28.3        
  [9] Rtsne_0.16                  munsell_0.5.0              
 [11] codetools_0.2-18            ica_1.0-3                  
 [13] pbdZMQ_0.3-8                future_1.30.0              
 [15] miniUI_0.1.1.1              withr_2.5.0                
 [17] spatstat.random_3.1-2       colorspace_2.1-0           
 [19] progressr_0.13.0            filelock_1.0.2             
 [21] uuid_1.1-0                  ROCR_1.0-11                
 [23] ggsignif_0.6.4              tensor_1.5                 
 [25] listenv_0.9.0               MatrixGenerics_1.6.0       
 [27] repr_1.1.4                  GenomeInfoDbData_1.2.7     
 [29] polyclip_1.10-4             bit64_4.0.5                
 [31] parallelly_1.34.0           vctrs_0.5.2                
 [33] generics_0.1.3              xfun_0.36                  
 [35] timechange_0.1.1            BiocFileCache_2.2.1        
 [37] R6_2.5.1                    bitops_1.0-7               
 [39] spatstat.utils_3.0-1        cachem_1.0.6               
 [41] DelayedArray_0.20.0         assertthat_0.2.1           
 [43] promises_1.2.0.1            BiocIO_1.4.0               
 [45] scales_1.2.1                googlesheets4_1.0.1        
 [47] rgeos_0.5-9                 gtable_0.3.1               
 [49] globals_0.16.2              goftest_1.2-3              
 [51] rlang_1.0.6                 RcppRoll_0.3.0             
 [53] splines_4.1.1               rtracklayer_1.54.0         
 [55] rstatix_0.7.1               lazyeval_0.2.2             
 [57] gargle_1.2.1                broom_1.0.1                
 [59] spatstat.geom_3.0-5         modelr_0.1.10              
 [61] yaml_2.3.7                  reshape2_1.4.4             
 [63] abind_1.4-5                 backports_1.4.1            
 [65] httpuv_1.6.8                tools_4.1.1                
 [67] ellipsis_0.3.2              RColorBrewer_1.1-3         
 [69] ggridges_0.5.4              plyr_1.8.8                 
 [71] base64enc_0.1-3             progress_1.2.2             
 [73] zlibbioc_1.40.0             RCurl_1.98-1.9             
 [75] prettyunits_1.1.1           deldir_1.0-6               
 [77] pbapply_1.7-0               cowplot_1.1.1              
 [79] zoo_1.8-11                  haven_2.5.1                
 [81] SummarizedExperiment_1.24.0 ggrepel_0.9.2              
 [83] cluster_2.1.4               fs_1.6.0                   
 [85] magrittr_2.0.3              scattermore_0.8            
 [87] reprex_2.0.2                lmtest_0.9-40              
 [89] RANN_2.6.1                  googledrive_2.0.0          
 [91] ProtGenerics_1.26.0         fitdistrplus_1.1-8         
 [93] matrixStats_0.63.0          hms_1.1.2                  
 [95] patchwork_1.1.2             mime_0.12                  
 [97] evaluate_0.20               xtable_1.8-4               
 [99] XML_3.99-0.13               readxl_1.4.1               
[101] gridExtra_2.3               compiler_4.1.1             
[103] biomaRt_2.50.3              KernSmooth_2.23-20         
[105] crayon_1.5.2                htmltools_0.5.4            
[107] tzdb_0.3.0                  later_1.3.0                
[109] lubridate_1.9.0             DBI_1.1.3                  
[111] dbplyr_2.2.1                MASS_7.3-58.1              
[113] rappdirs_0.3.3              car_3.1-1                  
[115] cli_3.6.0                   parallel_4.1.1             
[117] igraph_1.3.5                pkgconfig_2.0.3            
[119] GenomicAlignments_1.30.0    sp_1.6-0                   
[121] IRdisplay_1.1               plotly_4.10.1              
[123] spatstat.sparse_3.0-0       xml2_1.3.3                 
[125] XVector_0.34.0              rvest_1.0.3                
[127] digest_0.6.31               sctransform_0.3.5          
[129] RcppAnnoy_0.0.20            spatstat.data_3.0-0        
[131] Biostrings_2.62.0           cellranger_1.1.0           
[133] leiden_0.4.3                fastmatch_1.1-3            
[135] uwot_0.1.14                 restfulr_0.0.15            
[137] curl_5.0.0                  shiny_1.7.4                
[139] Rsamtools_2.10.0            rjson_0.2.21               
[141] lifecycle_1.0.3             nlme_3.1-160               
[143] jsonlite_1.8.4              carData_3.0-5              
[145] viridisLite_0.4.1           fansi_1.0.4                
[147] pillar_1.8.1                lattice_0.20-45            
[149] KEGGREST_1.34.0             fastmap_1.1.0              
[151] httr_1.4.4                  survival_3.4-0             
[153] glue_1.6.2                  png_0.1-8                  
[155] bit_4.0.5                   stringi_1.7.12             
[157] blob_1.2.3                  memoise_2.0.1              
[159] IRkernel_1.3.1              irlba_2.3.5.1              
[161] future.apply_1.10.0 
