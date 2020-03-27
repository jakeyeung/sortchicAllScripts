# Jake Yeung
# Date of Creation: 2020-03-27
# File: ~/projects/scchic/scripts/rstudioserver_analysis/stemcell_analysis/create_bed_annotated_regions.R
# Create a list of TSS, enhancer, etc and make count table out of it.

rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


inf.annot <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/ENCODE/genomic_regions_10kbpromoters.txt"
dat.annot <- fread(inf.annot)


