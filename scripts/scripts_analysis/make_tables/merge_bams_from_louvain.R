# Jake Yeung
# Date of Creation: 2019-03-19
# File: ~/projects/scchic/scripts/scripts_analysis/make_tables/merge_bams_from_louvain.R
# Do Louvain and then creaet list for bams to merge

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

outdir <- "~/data/scchic/tables/bamlist_for_merging"
dir.create(outdir)

load("~/data/scchic/robjs/TFactivity_genelevels_objects.RData", v=T)

jmarks.all <- list("H3K4me1" = "H3K4me1", "H3K4me3" = "H3K4me3", "H3K27me3" = "H3K27me3", "H3K9me3" = "H3K9me3")


barcodes <- read.table("data/barcode_summaries/barcodes/maya_384NLA.bc", col.names = "Barcodes", stringsAsFactors = FALSE)
experis <- read.table("data/outputs/experinames.txt", stringsAsFactors = FALSE)
experihash <- hash(sapply(experis$V1, function(x) strsplit(x, "_")[[1]][[1]]), sapply(experis$V1, function(x) paste(strsplit(x, "_")[[1]][2:3], collapse="_")))
cellhash <- hash(rownames(barcodes), unlist(barcodes))
cellhash.bc <- hash(unlist(barcodes), paste("cell", rownames(barcodes), sep = ""))


# Do louvains for clustering  ------------------------------------------------------------

nn.louv <- c(27, 27, 60, 60)
jmetric.louv='euclidean' 
jmindist.louv=0.4
jseed.louv=123

custom.settings.louv.lst <- lapply(nn.louv, function(x) GetUmapSettings(x, jmetric.louv, jmindist.louv, jseed.louv))

# topics.mat.lst <- lapply(out.objs, function(x) x$tm.result$topics)
topics.mat.lst <- lapply(tm.result.lst, function(tm.result) tm.result$topics)

dat.umap.lst <- mapply(function(custom.settings, topics.mat){
  dat.umap <- umap(topics.mat, config = custom.settings) 
  return(dat.umap)
}, custom.settings.lst, topics.mat.lst, SIMPLIFY = FALSE)
names(dat.umap.lst) <- jmarks.all

topics.mat.lst <- lapply(tm.result.lst, function(tm.result) tm.result$topics)

clstr.hash.lst <- mapply(function(topics.mat, custom.settings.louv) DoLouvain(topics.mat, custom.settings.louv, dat.umap.long = NULL), topics.mat.lst, custom.settings.louv.lst, SIMPLIFY = FALSE)

# assign cluster to dat umap

dat.umap.long.lst <- lapply(jmarks.all, function(jmark){
  dat.umap <- dat.umap.lst[[jmark]]
  dat.umap.long <- data.frame(umap1 = dat.umap$layout[, 1], umap2 = dat.umap$layout[, 2], cell = rownames(dat.umap$layout), stringsAsFactors = FALSE)
  dat.umap.long$louvain <- as.character(sapply(dat.umap.long$cell, function(x) clstr.hash.lst[[jmark]][[x]]))
  dat.umap.long$mark <- jmark
  return(dat.umap.long)
})

# get bam to louvain list
jmark <- "H3K4me1"
jmark <- "H3K4me3"
jmark <- "H3K27me3"
jmark <- "H3K9me3"

print(head(dat.umap.long.lst[[jmark]]))

# print to file

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
pdf("~/Dropbox/scCHiC_figs/FIG4_BM/primetime_plots/louvain_colors.pdf")
for (jmark in jmarks.all){
  m1 <- ggplot(dat.umap.long.lst[[jmark]], aes(x = umap1, y = umap2, color = louvain)) + geom_point(size=2.5) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jmark) + scale_color_manual(values = cbPalette)
  print(m1)
}
dev.off()

# change name to cell name
# write table summary for all 
dat.merge <- bind_rows(dat.umap.long.lst) %>% dplyr::select(cell, louvain, mark)
dat.merge$cellnew <- sapply(dat.merge$cell, MakeNewCellName.rev, experihash, cellhash)
dat.merge <- dat.merge %>% arrange(mark, louvain, cell)
data.table::fwrite(dat.merge, file = file.path(outdir, paste0("cell_to_bam_summary.txt")), sep = "\t", col.names = FALSE)

# write UMAP coordinates 
data.table::fwrite(bind_rows(dat.umap.long.lst), file = file.path(outdir, paste0("cell_to_bam_summary_withUmapCoordinates.txt")), sep = "\t", col.names = TRUE)

# write 


# write to output
for (jmark in jmarks.all){
  dat.sub <- subset(dat.umap.long.lst[[jmark]], select = c(cell, louvain))
  dat.sub$cellnew <- sapply(dat.sub$cell, function(x) MakeNewCellName.rev(x, experihash, cellhash))
  # outf <- paste0("~/data/scchic/tables/", jmark, "_cell_to_louvain.txt")
  # data.table::fwrite(dat.sub, file = outf, sep = "\t", col.names = FALSE)
  for (jclst in unique(base::sort(dat.sub$louvain))){
    dat.subsub <- subset(dat.sub, louvain == jclst)
    outf.sub <- file.path(outdir, paste0(paste("bamlist", jclst, jmark, sep = "-"), ".txt"))
    data.table::fwrite(dat.subsub %>% dplyr::select(cellnew), file = outf.sub, sep = "\t", col.names = FALSE)
  }
}

# 
# for (i in seq(nrow(dat.sub))){
#   x <- dat.sub$cell[[i]]
#   x <- "BM_H3K9me3_m1_rep2_cell364"
#   indx <- gsub("cell", "", strsplit(x, "_")[[1]][[5]])
#   bc <- AssignHash(indx, cellhash)
#   jsci <- "PZ"
#   jtiss <- strsplit(x, "_")[[1]][[1]]
#   jmark <- strsplit(x, "_")[[1]][[2]]
#   jmouse <- strsplit(x, "_")[[1]][[3]]
#   jrepl <- gsub("rep", "", strsplit(x, "_")[[1]][[4]])
#   jkey <- paste(jsci, jtiss, jmouse, jmark, jrepl, sep = "-")
#   jval <- AssignHash(jkey, experihash)
#   if (is.na(jval)){
#     print(x)
#   }
#   # print(jval)
# }
