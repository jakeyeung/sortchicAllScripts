# Jake Yeung
# Date of Creation: 2019-03-26
# File: ~/projects/scchic/scripts/scripts_analysis/make_tables/create_bam_list_for_peak_analysis_build95.R
# Peak analysis for build 95


rm(list=ls())

library(ggplot2)
library(ggrepel)

library(dplyr)
library(hash)

library(umap)
library(igraph)

source("scripts/Rfunctions/MaraDownstream.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/PlotFunctions.R")
source("scripts/Rfunctions/GetMetaCellHash.R")

# outdir <- "~/data/scchic/tables/bamlist_for_merging_build95"
outdir <- "~/data/scchic/tables/bamlist_for_peak_analysis_build95"
dir.create(outdir)

# load("~/data/scchic/robjs/TFactivity_genelevels_objects.RData", v=T)
load("~/data/scchic/robjs/gene_levels_build95.Rdata", v=T)

jmarks.all <- list("H3K4me1" = "H3K4me1", "H3K4me3" = "H3K4me3", "H3K27me3" = "H3K27me3", "H3K9me3" = "H3K9me3")

# need new experihash

# change name to cell name
# write table summary for all 
cellhash <- hash(rownames(barcodes), unlist(barcodes))
cellhash.bc <- hash(unlist(barcodes), paste("cell", rownames(barcodes), sep = ""))

dat.merge <- bind_rows(dat.umap.long.new.lst) %>% dplyr::select(cell, louvain, mark)
dat.merge$cellnew <- sapply(dat.merge$cell, MakeNewCellName.rev, experihash, cellhash)
dat.merge <- dat.merge %>% arrange(mark, louvain, cell)


# Write bam files for each mark -------------------------------------------

for (jmark in jmarks.all){
  dat.tmp <- subset(dat.merge, mark == jmark, select = c(cellnew))
  fwrite(dat.tmp, file = file.path(outdir, paste0("JY_", jmark, "_bamnames.out")), sep = "\t", col.names = FALSE)
}

