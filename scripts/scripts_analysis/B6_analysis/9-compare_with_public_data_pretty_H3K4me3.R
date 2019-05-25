# Jake Yeung
# Date of Creation: 2019-04-30
# File: ~/projects/scchic/scripts/scripts_analysis/sorted_bulk_analysis/pseudobulk_vs_sorted_chipseq_pretty.R
# Summarize it

# Some are log some are not ??

rm(list=ls())

library(ggplot2)
library(dplyr)
library(reshape2)
library(hash)
library(JFuncs)
library(data.table)

source("scripts/Rfunctions/PseudobulkComparisons.R")

# Load umap ---------------------------------------------------------------

inf.umap <- "/Users/yeung/data/scchic/robjs/B6_objs/traj_objs_all_marks.Rdata"
load(inf.umap, v=T)

# Constants ---------------------------------------------------------------

jylims <- c(0, 0.95)

dirmain.log <- paste0("/Users/yeung/data/scchic/from_cluster/pseudobulk_comparisons/merged_softlinks_textfile_log1p_build95_B6")
assertthat::assert_that(dir.exists(dirmain.log))
dirmain <- "/Users/yeung/data/scchic/from_cluster/pseudobulk_comparisons/merged_softlinks_textfile_build95_B6"
assertthat::assert_that(dir.exists(dirmain))
dirmainCompare <- dirmain
assertthat::assert_that(dir.exists(dirmainCompare))

# identify celltypes from Lara-Astiaso
la.ctypes <- data.table::fread("/Users/yeung/projects/scchic/data/lara-astiaso_celltypes_macbook.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE, col.names = c("mark_ctype"))
la.ctypes$ctype <- sapply(la.ctypes$mark_ctype, function(x) strsplit(x, "_")[[1]][[2]])
la.ctypes$mark <- sapply(la.ctypes$mark_ctype, function(x) strsplit(x, "_")[[1]][[1]])

themesize <- 18
xvar <- "umap1"; yvar <- "umap2"; jsize <- 2; jcol <- "louvain"

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#C0C0C0", "#32CD32", "D3D3D3")

# First do the activators -------------------------------------------------

# comparisons with Amit data does not need log


# Summarize H3K4me1 and H3K4me3 with Bcell and Neutrohil public da --------


# Summarize H3K4me1, H3K4me3 in linear comparisons
# for Bcell and Neutrophil


jmark <- "H3K4me3"
jthres <- 0.998

infs.h3k4me1 <- list.files(path = dirmain, pattern = paste0(jmark, ".*comparison.txt"))
# get celltype and replicate out of it to form the names of list
infs.h3k4me1.ctypes <- sapply(infs.h3k4me1, function(x) strsplit(x, "_")[[1]][[2]])
infs.h3k4me1.reps <- sapply(infs.h3k4me1, function(x){
  jrep <- strsplit(x, "_")[[1]][[3]]  # sometimes no rep, it will be called comparison.txt, switch this to rep1
  return(ifelse(jrep == "comparison.txt", "rep1", jrep))
})
infs.h3k4me1.names <- paste(infs.h3k4me1.ctypes, infs.h3k4me1.reps, sep = "_")
# attach dirmain to filename
infs.h3k4me1 <- paste(dirmain, infs.h3k4me1, sep = "/")
names(infs.h3k4me1) <- infs.h3k4me1.names
# out.lst.h3k4me1 <- lapply(names(infs.h3k4me1), function(jname) ComparePublicLinear(infs.h3k4me1[[jname]], thres=jthres, lab = jname))
out.lst.h3k4me1 <- lapply(names(infs.h3k4me1), function(jname) ComparePublicLinear(infs.h3k4me1[[jname]], thres=0.99, lab = jname))


# rename clusters
keys <- as.character(seq(11))
# keys <- sort(as.character(seq(11)))
vals <- c("NeutProgProg", "Eryth", "BcellProg", "NeutProg", "Neut", "NeutIsland", "Mega", "Tcell", "NK", "HSC?", "Bcell")
reorder.i <- sort(as.character(seq(length(vals))))
vals <- paste(vals, keys, sep = "_")
cols <- cbPalette[1:length(vals)]
kval.lst <- c(vals)
names(kval.lst) <- keys
clstr.hash.h3k4me1 <- hash(kval.lst)
col.hash.h3k4me1 <- hash(kval.lst[order(names(kval.lst))],  cols)

cors.merge <- purrr::reduce(.x = lapply(out.lst.h3k4me1, function(x) x$dat.cors), .f = bind_rows) %>%
  rowwise %>%
  mutate(Sample = strsplit(strsplit(as.character(Sample), "_")[[1]][[3]], "\\.")[[1]][[1]])

cors.merge.h3k4me1 <- cors.merge %>%
  rowwise() %>%
  mutate(Sample = clstr.hash.h3k4me1[[as.character(Sample)]],
         ctype = strsplit(compare, "_")[[1]][[1]],
         in.ido.amit = ifelse(ctype %in% la.ctypes$ctype, TRUE, FALSE)) %>%
  as.data.table() 

cors.merge.h3k4me1$ctype <- gsub("_rep1", "", cors.merge.h3k4me1$compare)
cors.merge.h3k4me1$color <- sapply(cors.merge.h3k4me1$Sample, function(x) col.hash.h3k4me1[[x]])


pdf(paste0("~/data/scchic/pdfs/B6_figures/compare_with_bulk/", jmark, "_ido_amit_comparison.pdf"), useDingbats = FALSE)
  
  m <- ggplot(dat.umap.long.trajs[[jmark]], aes_string(x = xvar, y = yvar, color = jcol)) + 
    # ggrastr::geom_point_rast(size = jsize) + 
    geom_point(size = jsize) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.border = element_blank(), legend.position = "bottom") + 
    xlab("") + ylab("") + scale_color_manual(values=cbPalette)
  print(m)
  
  cors.merge.h3k4me1[, ord := sprintf("%02i", frank(cors.merge.h3k4me1, compare, -corr.pears, ties.method = "first"))]
  m.h3k4me1.filt.pears <- ggplot(cors.merge.h3k4me1 %>% filter(in.ido.amit), aes(x = ord, y = corr.pears, fill = color)) + facet_wrap(~ctype, scales = "free_x") + 
    geom_bar(stat = "identity") + theme_bw(themesize/3) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") + 
    scale_x_discrete(labels = cors.merge.h3k4me1[, setNames(as.character(Sample), ord)]) + scale_fill_identity() + 
    xlab("") + ylab("Pearson Correlation") 
  print(m.h3k4me1.filt.pears)
  
  # plot individual
  compares <- unique(subset(cors.merge.h3k4me1 %>% filter(in.ido.amit), in.ido.amit)$ctype)
  
  for (jcompare in compares){
    print(jcompare)
    jtmp <- subset(cors.merge.h3k4me1, ctype == jcompare)
    jtmp <- OrderDecreasing(jtmp, "ord", "corr.pears")
    m.h3k4me1.filt.pears <- ggplot(jtmp, aes(x = ord, y = corr.pears, fill = color)) + 
      geom_bar(stat = "identity") + theme_bw(themesize) + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") + 
      scale_x_discrete(labels = cors.merge.h3k4me1[, setNames(as.character(Sample), ord)]) + 
      ggtitle(jcompare) + xlab("") + ylab("Pearson Correlation") +
      ylim(jylims) + scale_fill_identity()
    print(m.h3k4me1.filt.pears)
  }

dev.off()


# Do H3K4me3 on public data -----------------------------------------------



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
out.lst.h3k4me3 <- lapply(names(infs.h3k4me3), function(jname) ComparePublicLinear(infs.h3k4me3[[jname]], thres=jthres, lab = jname))

# reame clusters
keys <- as.character(seq(9))
vals <- c("HSC?", "Eryth", "Bcell", "Neut", "Mega", "NK", "NeutProg", "BcellProg", "NeutProgProg")
vals <- paste(vals, seq(length(vals)), sep = "_")
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


pdf(paste0("~/data/scchic/pdfs/B6_figures/compare_with_bulk/", jmark, "_ido_amit_comparison.pdf"), useDingbats = FALSE)

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

# Do H3k27me3 on erythryoblastsS? -----------------------------------------



jthres <- 0.99
jmark <- "H3K27me3"

m <- ggplot(dat.umap.long.trajs[[jmark]], aes_string(x = xvar, y = yvar, color = jcol)) + 
  # ggrastr::geom_point_rast(size = jsize) + 
  geom_point(size = jsize) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.border = element_blank()) + 
  xlab("") + ylab("") + scale_color_manual(values=cbPalette)
print(m)

# label clusters 1 to 7 to meaningful names
keys <- as.character(seq(8))
vals <- c("Mega", "Neut", "Eryth", "Tcell", "Bcell", "HSC?", "BcellProg", "NeutProg")
vals <- paste(vals, seq(length(vals)), sep = "_")
kval.lst <- c(vals)

cols <- cbPalette[1:length(keys)]

names(kval.lst) <- keys
clstr.hash.h3k27me3 <- hash(kval.lst)
# col.hash.h3k27me3 <- hash(vals, cols)
col.hash.h3k27me3 <- hash(kval.lst[order(names(kval.lst))], cols)

# LOG SPACE
jmark <- "H3K27me3"
compares_keep <- c("MatBcell", "ProB", "HSC")
# compares_keep <- c("Erythrobl", "Megakar")
compares_grepstr <- paste(compares_keep, collapse = "|")

infs.h3k27me3 <- list.files(path = dirmain.log, pattern = paste0(jmark, ".*comparison.txt"))
# get celltype and replicate out of it to form the names of list
infs.h3k27me3.ctypes <- sapply(infs.h3k27me3, function(x) strsplit(x, "_")[[1]][[2]])
infs.h3k27me3.reps <- sapply(infs.h3k27me3, function(x){
  jrep <- strsplit(x, "_")[[1]][[3]]  # sometimes no rep, it will be called comparison.txt, switch this to rep1
  return(ifelse(jrep == "comparison.txt", "rep1", jrep))
})
infs.h3k27me3.names <- paste(infs.h3k27me3.ctypes, infs.h3k27me3.reps, sep = "_")
# attach dirmain to filename
infs.h3k27me3 <- paste(dirmain.log, infs.h3k27me3, sep = "/")  # MAKE SURE THE PATH IS CORRECT 2019-04-30
names(infs.h3k27me3) <- infs.h3k27me3.names

out.lst <- lapply(names(infs.h3k27me3), function(jname) ComparePublicLog(infs.h3k27me3[[jname]], thres=jthres, lab = jname, clip.chic.only = TRUE))
# summarize for each cell type and rep
cors.merge <- purrr::reduce(.x = lapply(out.lst, function(x) x$dat.cors), .f = bind_rows) %>%
  rowwise %>%
  mutate(Sample = strsplit(strsplit(as.character(Sample), "_")[[1]][[3]], "\\.")[[1]][[1]])
# rename cors.merge
cors.merge.h3k27me3 <- cors.merge %>%
  rowwise() %>%
  mutate(Sample = kval.lst[[as.character(Sample)]]) %>%
  as.data.table()

cors.merge.h3k27me3$color <- sapply(cors.merge.h3k27me3$Sample, function(x) col.hash.h3k27me3[[x]])

cors.merge.h3k27me3[, ord := sprintf("%02i", frank(cors.merge.h3k27me3, "compare", -"corr.pears", ties.method = "first"))]
m.h3k27me3 <- ggplot(cors.merge.h3k27me3 %>% filter(grepl(compares_grepstr, compare)), aes_string(x = "ord", y = "corr.pears", fill = "color")) + facet_wrap(~compare, scales = "free_x") + 
  geom_bar(stat = "identity") + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_x_discrete(labels = cors.merge.h3k27me3[, setNames(as.character(Sample), ord)]) + 
  ggtitle(jmark) + xlab("") + ylab("Pearson Correlation") + 
  scale_fill_identity() 
print(m.h3k27me3)



# LINEAR SPACE
jthres <- 0.999
jmark <- "H3K27me3"
compares_keeplin <- c("Erythrobl", "Megakar")
# compares_keeplin <- c("MatBcell", "ProB", "HSC")
compares_grepstrlin <- paste(compares_keeplin, collapse = "|")

dirmain.linear <- "/Users/yeung/data/scchic/from_cluster/pseudobulk_comparisons/merged_softlinks_textfile"
assertthat::assert_that(dir.exists(dirmain.linear))
infs.h3k27me3.linear <- list.files(path = dirmain.linear, pattern = paste0(jmark, ".*comparison.txt"))
# get celltype and replicate out of it to form the names of list
infs.h3k27me3.ctypes.linear <- sapply(infs.h3k27me3.linear, function(x) strsplit(x, "_")[[1]][[2]])
infs.h3k27me3.reps.linear <- sapply(infs.h3k27me3.linear, function(x){
  jrep <- strsplit(x, "_")[[1]][[3]]  # sometimes no rep, it will be called comparison.txt, switch this to rep1
  return(ifelse(jrep == "comparison.txt", "rep1", jrep))
})
infs.h3k27me3.names.linear <- paste(infs.h3k27me3.ctypes.linear, infs.h3k27me3.reps.linear, sep = "_")
# attach dirmain to filename
infs.h3k27me3.linear <- paste(dirmain.linear, infs.h3k27me3.linear, sep = "/")
names(infs.h3k27me3.linear) <- infs.h3k27me3.names.linear

out.lst.linear <- lapply(names(infs.h3k27me3.linear), function(jname) ComparePublicLinear(infs.h3k27me3.linear[[jname]], thres=jthres, lab = jname))
# summarize for each cell type and rep
cors.merge.linear <- purrr::reduce(.x = lapply(out.lst.linear, function(x) x$dat.cors), .f = bind_rows) %>%
  rowwise %>%
  mutate(Sample = strsplit(strsplit(as.character(Sample), "_")[[1]][[3]], "\\.")[[1]][[1]])
# rename cors.merge
cors.merge.h3k27me3.linear <- cors.merge.linear %>%
  rowwise() %>%
  mutate(Sample = kval.lst[[as.character(Sample)]]) %>%
  as.data.table()
cors.merge.h3k27me3.linear$color <- sapply(cors.merge.h3k27me3.linear$Sample, function(x) col.hash.h3k27me3[[x]])

# compare linear
cors.merge.h3k27me3.linear[, ord := sprintf("%02i", frank(cors.merge.h3k27me3.linear, "compare", -"corr.pears", ties.method = "first"))]
m.h3k27me3.linear <- ggplot(cors.merge.h3k27me3.linear %>% filter(grepl(compares_grepstrlin, compare)), aes_string(x = "ord", y = "corr.pears", fill = "color")) + 
  facet_wrap(~compare, scales = "free_x") + 
  geom_bar(stat = "identity") + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_x_discrete(labels = cors.merge.h3k27me3.linear[, setNames(as.character(Sample), ord)]) + 
  ggtitle(jmark) + xlab("") + ylab("Pearson Correlation") + 
  scale_fill_identity()
print(m.h3k27me3.linear)

m.umap.h3k27me3 <- ggplot(dat.umap.long.trajs[[jmark]], aes_string(x = xvar, y = yvar, color = jcol)) + 
  # ggrastr::geom_point_rast(size = jsize) + 
  geom_point(size = jsize) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.border = element_blank(), legend.position = "bottom") + 
  xlab("") + ylab("") + scale_color_manual(values=cbPalette)
print(m.umap.h3k27me3)


# plot outputs along with the UMAP

pdf(paste0("~/data/scchic/pdfs/B6_figures/compare_with_bulk/", jmark, "_ido_amit_comparison.pdf"), useDingbats = FALSE)
print(m.umap.h3k27me3)
print(m.h3k27me3)
print(m.h3k27me3.linear)
dev.off()


# H3K9me3 -----------------------------------------------------------------

# no public chip-seq, but we try to compare


# Summarize H3K9me3 using pseudobulk correlations 
# for Bcell, Neutrophil, and Erythryoblasts 
# compare H3K9me3 

# jthres <- 0.999
jthres <- 0.999
jmark <- "H3K9me3"
bsize <- "100000"
refmark <- "H3K4me1"

inf.h3k9me3 <- list("11" = file.path(dirmain, paste0(jmark, "_Multimark_H3K4me1_vs_", jmark, "_comparison_binsize-", bsize, ".txt")),
                    "5" = file.path(dirmain, paste0(jmark, "_Multimark_H3K4me1_vs_", jmark, "_comparison_binsize-", bsize, ".txt")),
                    "2" = file.path(dirmain, paste0(jmark, "_Multimark_H3K4me1_vs_", jmark, "_comparison_binsize-", bsize, ".txt")),
                    "6" = file.path(dirmain, paste0(jmark, "_Multimark_H3K4me1_vs_", jmark, "_comparison_binsize-", bsize, ".txt")))

lapply(inf.h3k9me3, function(x) assertthat::assert_that(file.exists(x)))

# clster to name
labs <- paste("H3K4me1", sapply(names(inf.h3k9me3), function(x) clstr.hash.h3k4me1[[x]]), sep = " ")
# out.lst.multi <- lapply(names(inf.h3k9me3), function(jname) CompareAcrossMarks(inf.h3k9me3[[jname]], jmark, clstr = as.numeric(jname), 
#                                                                                thres=0.95, lab = jname))



# out.lst.multi.h3k9me3 <- mapply(function(jname, jlab) CompareAcrossMarks(inf.h3k9me3[[jname]], jmark, refmark, clstr = as.numeric(jname), thres = jthres, lab = jlab, clip.top.only = TRUE),
                                # names(inf.h3k9me3), labs, SIMPLIFY = FALSE)

out.lst.multi.h3k9me3 <- mapply(function(jname, jlab) CompareAcrossMarks(inf.h3k9me3[[jname]], jmark, refmark, clstr = as.numeric(jname), thres = 0.9, lab = jlab, clip.top.only = FALSE),
                                names(inf.h3k9me3), labs, SIMPLIFY = FALSE)

# out.lst.multi.h3k9me3 <- mapply(function(jname, jlab) CompareAcrossMarks(inf.h3k9me3[[jname]], jmark, refmark, clstr = as.numeric(jname), thres = 1, lab = jlab, clip.top.only = FALSE),
#                                 names(inf.h3k9me3), labs, SIMPLIFY = FALSE)

nclst.h3k9me3 <- 7
keys <- as.character(seq(nclst.h3k9me3))
vals <- c("NeutroProg", "Center", "Bcell", "Neutro", "BcellProg", "NeutroProgProg", "Eryth")
vals <- paste(vals, seq(length(vals)), sep = "_")
kval.lst <- c(vals)
names(kval.lst) <- keys
cols <- cbPalette[1:length(keys)]

clstr.hash.multi.h3k9me3 <- hash(kval.lst)
col.hash.multi.h3k9me3 <- hash(kval.lst[order(names(kval.lst))], cols)

cors.merge.multi.h3k9me3 <- bind_rows(out.lst.multi.h3k9me3[[1]]$dat.cors, out.lst.multi.h3k9me3[[2]]$dat.cors, out.lst.multi.h3k9me3[[3]]$dat.cors, out.lst.multi.h3k9me3[[4]]$dat.cors) %>%
  rowwise() %>%
  mutate(Sample = strsplit(strsplit(as.character(Sample), "_")[[1]][[3]], "\\.")[[1]][[1]]) %>%
  mutate(Sample = clstr.hash.multi.h3k9me3[[as.character(Sample)]])

cors.merge.multi.h3k9me3 <- cors.merge.multi.h3k9me3 %>%
  group_by(compare) %>%
  do(OrderDecreasing(., jfactor = "Sample", jval = "corr.spear.log2"))
cors.merge.multi.h3k9me3 <- as.data.table(cors.merge.multi.h3k9me3)
# print(cors.merge.multi.h3k9me3)

cors.merge.multi.h3k9me3$color <- sapply(as.character(cors.merge.multi.h3k9me3$Sample), function(x) col.hash.multi.h3k9me3[[x]])
cors.merge.multi.h3k9me3[, ord := sprintf("%02i", frank(cors.merge.multi.h3k9me3, "compare", "corr.pears", ties.method = "first"))]

m.multi.h3k9me3 <- ggplot(cors.merge.multi.h3k9me3, aes(x = ord, y = corr.pears, fill = color)) + facet_wrap(~compare, scales = "free_x") + geom_bar(stat = "identity") + 
  theme_bw(themesize/3) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_x_discrete(labels = cors.merge.multi.h3k9me3[, setNames(as.character(Sample), ord)]) + 
  ggtitle(jmark) + scale_fill_identity() + xlab("") + ylab("Pearson Correlation")
print(m.multi.h3k9me3)

cors.merge.multi.h3k9me3[, ord := sprintf("%02i", frank(cors.merge.multi.h3k9me3, "compare", "corr.pears.log2", ties.method = "first"))]
m.multi.h3k9me3.log2 <- ggplot(cors.merge.multi.h3k9me3, aes(x = ord, y = corr.pears.log2, fill = color)) + facet_wrap(~compare, scales = "free_x") + geom_bar(stat = "identity") +
  theme_bw(themesize/3) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = cors.merge.multi.h3k9me3[, setNames(as.character(Sample), ord)]) +
  ggtitle(jmark) + scale_fill_identity() + xlab("") + ylab("Pearson Correlation")
print(m.multi.h3k9me3.log2)


pdf(paste0("~/data/scchic/pdfs/B6_figures/compare_with_bulk/", jmark, "_multimark_comparison.pdf"), useDingbats = FALSE)
  m <- ggplot(dat.umap.long.trajs[[jmark]], aes_string(x = xvar, y = yvar, color = jcol)) + 
    # ggrastr::geom_point_rast(size = jsize) + 
    geom_point(size = jsize) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.border = element_blank(), legend.position = "bottom") + 
    xlab("") + ylab("") + scale_color_manual(values=cbPalette)
  print(m)
  print(m.multi.h3k9me3)
  for (ctype in unique(cors.merge.multi.h3k9me3$compare)){
    print(ctype)
    jsub <- as.data.table(cors.merge.multi.h3k9me3 %>% filter(compare == ctype))
    m.multi.h3k9me3.indiv <- ggplot(jsub, aes(x = ord, y = corr.pears.log2, fill = color)) + facet_wrap(~compare, scales = "free_x") + geom_bar(stat = "identity") + 
      theme_bw(themesize) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
      scale_x_discrete(labels = jsub[, setNames(as.character(Sample), ord)]) + 
      ggtitle(jmark) + scale_fill_identity() + xlab("") + ylab("Pearson Correlation")
    print(m.multi.h3k9me3.indiv)
  }
dev.off()


# 
# 
# # Debug -------------------------------------------------------------------
# 
# 
# # debug
# ggplot(out.lst.multi.h3k9me3[[2]]$dat.multi.long.filt, aes(x = chic_ref, y = chic)) + geom_point(alpha = 0.25) +
#   scale_x_log10() + scale_y_log10() +
#   facet_wrap(~Sample) + theme_bw() +
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# ggplot(out.lst.multi.h3k9me3[[3]]$dat.multi.long.filt, aes(x = chic_ref, y = chic)) + geom_point(alpha = 0.25) +
#   scale_x_log10() + scale_y_log10() +
#   facet_wrap(~Sample) + theme_bw() +
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# 
# ggplot(out.lst.multi.h3k9me3[[4]]$dat.multi.long.filt, aes(x = chic_ref, y = chic)) + geom_point(alpha = 0.25) +
#   scale_x_log10() + scale_y_log10() +
#   facet_wrap(~Sample) + theme_bw() +
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# 
# # clip differently on both sides
# jsub <- subset(out.lst.multi.h3k9me3[[3]]$dat.multi.long.filt, chic > 0.01 & chic < 2 & chic_ref > 0.1 & chic_ref < 2)
# 
# thres <- 0.9
# jsub <- subset(out.lst.multi.h3k9me3[[2]]$dat.multi.long.filt) %>%
#   group_by(Sample) %>%
#   filter(!is.nan(chic) & !is.nan(chic_ref)) %>%  # otherwise fails
#   filter(chic <= quantile(chic, probs = thres) & chic >= quantile(chic, probs = 1 - thres)) %>%
#   filter(chic_ref <= quantile(chic_ref, probs = thres) & chic_ref >= quantile(chic_ref, probs = 1 - thres))
# 
# ggplot(jsub, aes(x = chic_ref, y = chic)) + geom_point(alpha = 0.25) +
#   scale_x_log10() + scale_y_log10() +
#   facet_wrap(~Sample) + theme_bw() +
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   geom_smooth(method = "lm", se = FALSE)
# 
# jcors <- jsub %>%
#   group_by(Sample) %>%
#   summarise(cor.pears = cor(chic, chic_ref),
#             cor.pears.log2 = cor(log2(chic), log2(chic_ref)))
# print(jcors)
# 
# barplot(jcors$cor.pears.log2)
