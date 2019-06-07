# Jake Yeung
# Date of Creation: 2019-06-06
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/filter_stringent_H3K4me3/compare_pseudobulk_stringent.R
# Compare pseudobulk with new H3K4me3 dataset



rm(list=ls())

library(ggplot2)
library(dplyr)
library(reshape2)
library(hash)
library(JFuncs)
library(data.table)

source("scripts/Rfunctions/PseudobulkComparisons.R")

# Load umap ---------------------------------------------------------------

inf.umap <- "/Users/yeung/data/scchic/robjs/B6_objs_stringent/traj_objs_H3K4me3.stringent.Rdata"
load(inf.umap, v=T)

# Constants ---------------------------------------------------------------

jylims <- c(0, 0.95)

dirmain <- "/Users/yeung/data/scchic/from_cluster/pseudobulk_comparisons/merged_softlinks_textfile_build95_B6_stringent_mergedDir"
assertthat::assert_that(dir.exists(dirmain))
dirmainCompare <- dirmain
assertthat::assert_that(dir.exists(dirmainCompare))

# identify celltypes from Lara-Astiaso
la.ctypes <- data.table::fread("/Users/yeung/projects/scchic/data/lara-astiaso_celltypes_macbook.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE, col.names = c("mark_ctype"))
la.ctypes$ctype <- sapply(la.ctypes$mark_ctype, function(x) strsplit(x, "_")[[1]][[2]])
la.ctypes$mark <- sapply(la.ctypes$mark_ctype, function(x) strsplit(x, "_")[[1]][[1]])

themesize <- 24
xvar <- "umap1"; yvar <- "umap2"; jsize <- 2; jcol <- "louvain"

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#C0C0C0", "#32CD32", "D3D3D3")



# Do H3K4me3 on public data -----------------------------------------------



# jthres <- 0.998
jthres <- 0.998

jmark <- "H3K4me3"
infs.h3k4me3 <- list.files(path = dirmain, pattern = paste0(jmark, ".*comparison.txt"))
# get celltype and replicate out of it to form the names of list
infs.h3k4me3.ctypes <- sapply(infs.h3k4me3, function(x) strsplit(x, "_")[[1]][[2]])
infs.h3k4me3.reps <- sapply(infs.h3k4me3, function(x){
  jrep <- strsplit(x, "_")[[1]][[3]]  # sometimes no rep, it will be called comparison.txt, switch this to rep1
  return(ifelse(jrep == "comparison.txt", "rep1", jrep))
})
infs.h3k4me3.names <- paste(infs.h3k4me3.ctypes, infs.h3k4me3.reps, sep = "_")
# attach dirmain to filename
infs.h3k4me3 <- paste(dirmain, infs.h3k4me3, sep = "/")
names(infs.h3k4me3) <- infs.h3k4me3.names
out.lst.h3k4me3 <- lapply(names(infs.h3k4me3), function(jname) ComparePublicLinear(infs.h3k4me3[[jname]], thres=jthres, lab = jname, filter.min = FALSE))

# reame clusters
keys <- as.character(seq(8))
# vals <- c("HSC?", "Eryth", "Bcell", "Neut", "Mega", "NK", "NeutProg", "BcellProg", "NeutProgProg")
vals <- rep("", length(keys))
# vals <- paste(vals, seq(length(vals)), sep = "")
vals <- paste(vals, seq(length(vals)), sep = "")
cols <- cbPalette[1:length(vals)]
kval.lst <- c(vals)
names(kval.lst) <- keys
clstr.hash.h3k4me3 <- hash(kval.lst)
# col.hash.h3k4me3 <- hash(vals, cols)
col.hash.h3k4me3 <- hash(kval.lst[order(names(kval.lst))], cols)


cors.merge <- purrr::reduce(.x = lapply(out.lst.h3k4me3, function(x) x$dat.cors), .f = bind_rows) %>%
  rowwise %>%
  mutate(Sample = strsplit(strsplit(as.character(Sample), "_")[[1]][[3]], "\\.")[[1]][[1]])

cors.merge.h3k4me3 <- cors.merge %>%
  rowwise() %>%
  mutate(Sample = clstr.hash.h3k4me3[[as.character(Sample)]],
         ctype = strsplit(compare, "_")[[1]][[1]],
         in.ido.amit = ifelse(ctype %in% la.ctypes$ctype, TRUE, FALSE)) %>%
  as.data.table() 


cors.merge.h3k4me3$color <- sapply(cors.merge.h3k4me3$Sample, function(x) col.hash.h3k4me3[[x]])

cors.merge.h3k4me3[, ord := sprintf("%02i", frank(cors.merge.h3k4me3, compare, -corr.pears, ties.method = "first"))]

print(subset(cors.merge.h3k4me3, ctype == "Mono") %>% arrange(desc(corr.pears)))


pdf(paste0("~/data/scchic/pdfs/B6_figures/compare_with_bulk_stringent/", jmark, "_ido_amit_comparison_stringent.pdf"), useDingbats = FALSE)
m <- ggplot(dat.umap.long.trajs[[jmark]], aes_string(x = xvar, y = yvar, color = jcol)) + 
  # ggrastr::geom_point_rast(size = jsize) + 
  geom_point(size = jsize) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.border = element_blank(), legend.position = "bottom") + 
  xlab("") + ylab("") + scale_color_manual(values=cbPalette)
print(m)

m.h3k4me3.pears.filt <- ggplot(cors.merge.h3k4me3 %>% filter(in.ido.amit), aes(x = ord, y = corr.pears, fill = color)) + facet_wrap(~ctype, scales = "free_x") + 
  geom_bar(stat = "identity") + theme_bw(themesize/3) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") + 
  scale_x_discrete(labels = cors.merge.h3k4me3[, setNames(as.character(Sample), ord)]) + scale_fill_identity() + 
  xlab("") + ylab("Pearson Correlation") 
print(m.h3k4me3.pears.filt)

# individual celltypes
for (jcompare in unique(cors.merge.h3k4me3$ctype)){
  print(jcompare)
  jtmp <- subset(cors.merge.h3k4me3, ctype == jcompare)
  jtmp <- OrderDecreasing(jtmp, "ord", "corr.pears")
  m.h3k4me1.filt.pears <- ggplot(jtmp, aes(x = ord, y = corr.pears, fill = color)) + 
    geom_bar(stat = "identity") + theme_bw(themesize) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") + 
    scale_x_discrete(labels = cors.merge.h3k4me3[, setNames(as.character(Sample), ord)]) + 
    ggtitle(jcompare) + xlab("") + ylab("Pearson Correlation") +
    ylim(jylims) + scale_fill_identity()
  print(m.h3k4me1.filt.pears)
}
dev.off()