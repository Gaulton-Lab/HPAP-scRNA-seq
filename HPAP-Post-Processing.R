#This script will go over how to process and cluster your final Seurat object after performing SoupX and Scrublet

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
