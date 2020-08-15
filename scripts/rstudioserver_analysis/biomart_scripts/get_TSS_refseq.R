# Jake Yeung
# Date of Creation: 2020-08-06
# File: ~/projects/scchic/scripts/rstudioserver_analysis/biomart_scripts/get_TSS_refseq.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(biomaRt)



library(GenomicRanges)

# Constants ---------------------------------------------------------------
# 
#   
# winsize <- 50000L
# 
# 
# # Biomart init ------------------------------------------------------------
# 
# 
# mart.obj <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'drerio_gene_ensembl')
# 
# gos <- getBM(
#   attributes=c("external_gene_name", "chromosome_name", "transcription_start_site", "start_position", "end_position
#   mart=mart.obj)
# 
# #
# print(unique(gos$chromosome_name))
# 
# # chromos <- c(seq(19), c("X", "Y"))
# chromos <- seq(25)
# chromos.withprefix <- paste("chr", chromos, sep = "")
# 
# 
# 