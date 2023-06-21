#Load necessary packages
suppressMessages(library(Seurat))
suppressMessages(library(stringr))
suppressMessages(library(parallel))
suppressMessages(library(readr))
suppressMessages(library(DESeq2))
suppressMessages(library(beeswarm))
suppressMessages(library(limma))
suppressMessages(library(edgeR))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(data.table))

########################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
###Step 1: Make Pseudobulk Matrices###
#Read in final Seurat object
data <- readRDS('~/hpap/hpap_islet_scRNAseq.rds')
Idents(data) <- data@meta.data$cell_type

#Pull out list of all cell types
unique_cell_types <- unique(data$cell_type)

DefaultAssay(data) <- 'RNA'
#Pull out gene expression counts
gex.counts <- GetAssayData(data,slot='counts')

##Pull out sample barcodes
sample_bcs <- list()
samples <- c('HPAP-019','HPAP-020','HPAP-021','HPAP-022','HPAP-023','HPAP-024','HPAP-026','HPAP-028','HPAP-029','HPAP-032','HPAP-034','HPAP-035','HPAP-036','HPAP-037','HPAP-038','HPAP-039','HPAP-040','HPAP-042','HPAP-043','HPAP-044','HPAP-045','HPAP-047','HPAP-049','HPAP-050','HPAP-051','HPAP-052','HPAP-053','HPAP-054','HPAP-055','HPAP-056','HPAP-057','HPAP-058','HPAP-059','HPAP-061','HPAP-063','HPAP-064','HPAP-065','HPAP-070','HPAP-071','HPAP-072','HPAP-074','HPAP-075','HPAP-077','HPAP-079','HPAP-080','HPAP-081','HPAP-082','HPAP-083','HPAP-084','HPAP-085','HPAP-087','HPAP-088','HPAP-090','HPAP-091','HPAP-092','HPAP-099','HPAP-100','HPAP-101','HPAP-103','HPAP-104','HPAP-105','HPAP-106','HPAP-107','HPAP-108','HPAP-109')
for (sample in samples){
    sample_bcs[[sample]] <- row.names(data[[]][data[[]]$library == sample,])
}

#Write function that will save sample gene expression count matrices for each unique cell type
##Note: Since we ran SoupX with the roundToInt flag, counts are SoupX corrected but still integers for DESeq's negative binomial model
data_matrices <- data

get_per_sample_gex_SUMS <- function(cell.type, mtx.fp){
    print(cell.type)
    #Pull out rows of gex.counts where barcode Ident matches cell.type
    bcs <- names(Idents(data_matrices)[Idents(data_matrices) == cell.type])
    counts <- gex.counts[,colnames(gex.counts) %in% bcs]
    print(dim(counts))

    #Initialize the sample gex matrix
    counts.df <- as.data.frame(rep(0,length(row.names(gex.counts))))
    row.names(counts.df) <- row.names(gex.counts)
    colnames(counts.df) <- c('temp')

    #Loop through samples and calculate sum of gex values
    for (sample in samples){
        sample_cols <- colnames(counts) %in% sample_bcs[[sample]]
        counts.cut <- counts[,sample_cols]
        
        #If only one barcode, this becomes a vector which is an issue
        if (typeof(counts.cut) == 'double'){
            sum.counts <- counts.cut
        #If there are no barcodes, this will return NA (just return 0 for everything)
        } else if(length(colnames(counts.cut)) == 0){
            sum.counts <- rep(0,length(row.names(counts)))
        } else {
            sum.counts <- rowSums(counts.cut)
        }
        counts.df <- cbind(counts.df,as.data.frame(sum.counts))
     }
    fin.counts.df <- counts.df[,-c(1)]
    colnames(fin.counts.df) <- samples
    head(fin.counts.df)

    #Export cell type specific gene by sample matrices
    mtx.fp <- sprintf('~/hpap/deseq/%s_sample_gex_total_counts.txt',cell.type)
    write.table(fin.counts.df,mtx.fp,sep='\t',quote=FALSE)
}

#Run function to make matrices
for (cell.type in unique_cell_types){
    fp = sprintf('~/hpap/deseq/%s_pseudobulk.txt',cell.type)
    get_per_sample_gex_SUMS(cell.type, fp)
}

###Step 2: Make TPM Matrices###
#We want to normalize for gene size. Need to extract gene sizes from Gencode v38 publicly available annotations
#Pull out gene exon info and calculate effective length
suppressMessages(txdb <- makeTxDbFromGFF('~/publicdata/gencode_v38/gencode.v38.annotation_comprehensive_CHR.gtf',format='gtf')) #This gtffile can be downloaded from Gencode
exons.list.per.gene <- exonsBy(txdb,by='gene') #Collect the exons per gene_id

#Reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then
exonic.gene.sizes <- sum(width(reduce(exons.list.per.gene)))

#This creates a separate gene info table so we can add gene names to exonic sizes later
gene.info <- rtracklayer::import('~/publicdata/gencode_v38/gencode.v38.annotation_comprehensive_CHR.gtf')
gene.info <- as.data.frame(gene.info)
gene.info <- gene.info[,c('gene_id', 'gene_name', 'gene_type', 'seqnames', 'start', 'end', 'strand', 'source', 'level')]
colnames(gene.info) <- c('gene_id', 'gene_name', 'gene_type', 'chrom', 'start', 'end', 'strand', 'source', 'level')

#Add the effective lengths to the original gene.info dataframe
temp_df <- gene.info
rownames(temp_df) <- gene.info$gene_id
temp_df2 <- as.data.frame(exonic.gene.sizes)
temp_df2$gene_id <- rownames(temp_df2)
new_df <- merge(temp_df,temp_df2, by='gene_id', all=TRUE)

#Remove duplicate rows from gene info df
fin.gene.info <- new_df[!duplicated(new_df$gene_name),]
write.table(fin.gene.info, '~/publicdata/gencode_v38/gene_info_withExonicGeneSizes.tsv', quote=FALSE, col.names=TRUE, row.names=FALSE, sep='\t')

#Now that we have extracted  gene sizes, can make our TPM matrices
#Read in psedobulk gex matrices from above 
dir <- '~/hpap/deseq/'
files <- list.files(dir, pattern ='sample_gex_total_counts.txt')
cells <- gsub('_sample_gex_total_counts.txt','', files)

make_tpm <- function(raw_counts, gene_sizes){
    rpk <- raw_counts / gene_sizes
    tpm <- rpk
    for (i in 1:ncol(rpk)){
        tpm[,i] <- rpk[,i]/(sum(rpk[,i])/1e6)
        }
    return(tpm)
    }

for (x in files){
    cell <- cells[which(files == x)]
    raw_counts <- read.table(paste0(dir, x), row.names=1)

    raw_counts <- subset(raw_counts ,rownames(raw_counts) %in% fin.gene.info$gene_name)
    gene_sizes <- fin.gene.info$exonic.gene.sizes[match(rownames(raw_counts), fin.gene.info$gene_name )]

    tpm_mat <- make_tpm(raw_counts, gene_sizes)
    write.table(tpm_mat, paste0('~/hpap/deseq/', cell , '_TPM_per_sample.txt'), sep='\t', quote=FALSE)
    }

###Step 3: Run DESeq###
meta <- read.csv('~/hpap/Donor_Summary_127.csv') #This file can be downloaded from PANC-DB
meta <- meta[which(meta$donor_ID %in% data@meta.data$library),]
meta$library <- gsub('-', '.', meta$donor_ID)

#Add some metadata columns from Seurat object to HPAP metadata (mostly for how it is formatted)
diabetes_status  <- data@meta.data[,c('library', 'sex', 'condition', 'chemistry', 'tissue_source')]
rownames(diabetes_status) <- NULL
meta2 <- merge(meta2, diabetes_status, by.x = 'donor_ID', by.y='library')

#Scale age and BMI since they are numeric covariates- helps with DESeq model fitting
meta2$scaled_bmi <- scale(meta2$bmi)
meta2$scaled_age_years <- scale(meta2$age_years)

#Get list of pseudobulk files
files <- list.files('~/hpap/deseq/',pattern='_sample_gex_total_counts.txt')
#Create outdir for results
outdir <- '~/hpap/deseq/results/'
dir.create(outdir)

#Remove file suffixes to get cell type names
cells <- gsub('_sample_gex_total_counts.txt','', files) 

#Need sample diabetes status for filtering counts
nd  <- unique(data@meta.data[which(data@meta.data$condition == 'ND'),]$library)
aab  <- unique(data@meta.data[which(data@meta.data$condition == 'AAB+'),]$library)
t1d  <- unique(data@meta.data[which(data@meta.data$condition == 'T1D'),]$library)
t2d  <- unique(data@meta.data[which(data@meta.data$condition == 'T2D'),]$library)

t1d <- gsub('-', '.', cond_t1d)
t2d <- gsub('-', '.', cond_t2d)
aab <- gsub('-', '.', aab)
nd <- gsub('-', '.', nd)

###Example code for running T1D vs. non-diabetic. Can perform pairwise tests for all other conditions as well.
for (x in files) {
    cell <- gsub('_sample_gex_total_counts.txt', '', x)
    raw_counts <- read.table(paste0(dir, x), row.names=1)
    raw_counts <- raw_counts[,(colSums(raw_counts != 0) > 0)]
    meta_cell <- subset(meta2, library %in% colnames(raw_counts))

    raw_counts <- raw_counts[,which(colnames(raw_counts) %in% meta2$library)]

    t1d_raw_counts <- raw_counts[which(colnames(raw_counts)%in% t1d)]
    nd_raw_counts <- raw_counts[which(colnames(raw_counts)%in% nd)]

    n_nd <- floor(ncol(nd_raw_counts)/2)  #Here I am calculating half the sample size for each condition, rounded down. I will use this as a minimum number of samples to meet gene count thresholds
    n_t1d <- floor(ncol(t1d_raw_counts)/2)
    
    #The  next few lines of code subset the count matrix to only include genes where at least half the samples per condition have greater than 5 counts
    ##Genes must meet this criteria in both conditions (see intersect below). This helps avoid significance due to one condition largely having 0 counts
    nd_genes_to_keep <- c()
        for (i in 1:nrow(nd_raw_counts)) {
          if (sum(nd_raw_counts[i, ] >= 5) >= n_nd) {
            nd_genes_to_keep <- c(nd_genes_to_keep, rownames(nd_raw_counts[i, ]))
          }
        }

    t1d_genes_to_keep <- c()
        for (i in 1:nrow(t1d_raw_counts)) {
          if (sum(t1d_raw_counts[i, ] >= 5) >= n_t1d) {
            t1d_genes_to_keep <- c(t1d_genes_to_keep, rownames(t1d_raw_counts[i, ]))
          }
        }

    nd_t1d <- intersect(nd_genes_to_keep, t1d_genes_to_keep)

        if ('T1D' %in% data@meta.data$condition){
            print(cell)
            print('All conditions found!')
            counts <- raw_counts[which(rownames(raw_counts) %in% nd_t1d),] #Subsets full count matrix to only genes that met our filtering criteria

            my_design <- as.formula ('~ condition + scaled_bmi + sex + scaled_age_years + chemistry + tissue_source')
            dds <- DESeqDataSetFromMatrix(round(counts), colData = meta_cell, design = my_design) #colData is where design columns are found
            dds <- estimateSizeFactors(dds)
            dds <- estimateDispersions(dds)

            ### Pairwise Wald test: conditon vs control  
            dds <- nbinomWaldTest(dds)
            tests <- c('T1D') #Conditions to test compared to control
                for (y in tests){
                    res <- results(dds, contrast=c('condition',y,'ND'))
                    res <- as.data.frame(res)
                    res <- res[order(res$pvalue),]
                    outfile <- paste0( cell, '.deseq.WaldTest.', y , '.tsv')
                    write.table(res,paste0('~/hpap/deseq/results/', outfile) , sep='\t', quote=FALSE)           
                }     
        }
} 

########################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
sessionInfo()
R version 4.1.1 (2021-08-10)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.2 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] data.table_1.14.6           SeuratObject_4.1.3         
 [3] Seurat_4.3.0                GenomicFeatures_1.46.5     
 [5] AnnotationDbi_1.56.2        edgeR_3.36.0               
 [7] limma_3.50.3                beeswarm_0.4.0             
 [9] DESeq2_1.34.0               SummarizedExperiment_1.24.0
[11] Biobase_2.54.0              MatrixGenerics_1.6.0       
[13] matrixStats_0.63.0          GenomicRanges_1.46.1       
[15] GenomeInfoDb_1.30.1         IRanges_2.28.0             
[17] S4Vectors_0.32.4            BiocGenerics_0.40.0        
[19] readr_2.1.2                 stringr_1.5.0 
