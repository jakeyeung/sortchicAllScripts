# Jake Yeung
# Date of Creation: 2019-05-15
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/9-compare_with_public_data.R
# Compare with public data


rm(list=ls())

tstart <- Sys.time() 

library(ggplot2)
library(dplyr)
library(reshape2)
library(hash)
library(JFuncs)
library(data.table)

# Functions ---------------------------------------------------------------

suffix <- ""
refmark <- "H3K4me1"
source("scripts/Rfunctions/PseudobulkComparisons.R")

plot.to.file <- TRUE

if (plot.to.file){
  pdf(paste0("~/data/scchic/pdfs/B6_figures/compare_pseudo_build95_with_amit", suffix, ".pdf"), useDingbats = FALSE)
}

# Set constants -----------------------------------------------------------

jthres <- 0.995

# identify celltypes from Lara-Astiaso
la.ctypes <- data.table::fread("/Users/yeung/projects/scchic/data/lara-astiaso_celltypes_macbook.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE, col.names = c("mark_ctype"))
la.ctypes$ctype <- sapply(la.ctypes$mark_ctype, function(x) strsplit(x, "_")[[1]][[2]])
la.ctypes$mark <- sapply(la.ctypes$mark_ctype, function(x) strsplit(x, "_")[[1]][[1]])

dirmain <- paste0("/Users/yeung/data/scchic/from_cluster/pseudobulk_comparisons/merged_softlinks_textfile", suffix, "_build95_B6")
assertthat::assert_that(dir.exists(dirmain))

# Compare H3K27me3 with public data ---------------------------------------

# compare H3K27me3 data for all public data

jmark <- "H3K27me3"

infs.h3k27me3 <- list.files(path = dirmain, pattern = paste0(jmark, ".*comparison.txt"))
# get celltype and replicate out of it to form the names of list
infs.h3k27me3.ctypes <- sapply(infs.h3k27me3, function(x) strsplit(x, "_")[[1]][[2]])
infs.h3k27me3.reps <- sapply(infs.h3k27me3, function(x){
  jrep <- strsplit(x, "_")[[1]][[3]]  # sometimes no rep, it will be called comparison.txt, switch this to rep1
  return(ifelse(jrep == "comparison.txt", "rep1", jrep))
})
infs.h3k27me3.names <- paste(infs.h3k27me3.ctypes, infs.h3k27me3.reps, sep = "_")
# attach dirmain to filename
infs.h3k27me3 <- paste(dirmain, infs.h3k27me3, sep = "/")
names(infs.h3k27me3) <- infs.h3k27me3.names

out.lst <- lapply(names(infs.h3k27me3), function(jname) ComparePublicLinear(infs.h3k27me3[[jname]], thres=jthres, lab = jname))

# summarize for each cell type and rep

cors.merge <- purrr::reduce(.x = lapply(out.lst, function(x) x$dat.cors), .f = bind_rows) %>%
  rowwise %>%
  mutate(Sample = strsplit(strsplit(as.character(Sample), "_")[[1]][[3]], "\\.")[[1]][[1]])

print(cors.merge %>% filter(compare == "MatBcell_rep1"))

# label clusters 1 to 8 to meaningful names
keys <- as.character(seq(8))
vals <- c("Mega", "Neut", "Eryth", "Tcell", "Bcell", "HSC?", "BcellProg", "NeutProg")
vals <- paste(vals, seq(length(vals)), sep = "_")
# vals <- c("Eryth", "BetweenNeuAndBcells?", "HSC?", "Bcell", "Neutrophil", "NeutroProgs?", "BcellProgs?")
kval.lst <- c(vals)
names(kval.lst) <- keys
clstr.hash.h3k27me3 <- hash(kval.lst)

# rename cors.merge
cors.merge.h3k27me3 <- cors.merge %>%
  rowwise() %>%
  mutate(Sample = kval.lst[[as.character(Sample)]]) %>%
  as.data.table()

print(subset(cors.merge.h3k27me3, compare == "MatBcell_rep1"))

ggplot(cors.merge.h3k27me3, aes(x = Sample, y = corr.spear)) + geom_col() + facet_wrap(~compare, scales = "free_x") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

cors.merge.h3k27me3[, ord := sprintf("%02i", frank(cors.merge.h3k27me3, "compare", -"corr.pears", ties.method = "first"))]
m.h3k27me3 <- ggplot(cors.merge.h3k27me3, aes_string(x = "ord", y = "corr.pears")) + facet_wrap(~compare, scales = "free_x") + 
  geom_bar(stat = "identity") + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_x_discrete(labels = cors.merge.h3k27me3[, setNames(as.character(Sample), ord)]) + 
  ggtitle(jmark) + xlab("Cluster")
print(m.h3k27me3)

cors.merge.h3k27me3[, ord := sprintf("%02i", frank(cors.merge.h3k27me3, "compare", -"corr.pears.log", ties.method = "first"))]
m.h3k27me3 <- ggplot(cors.merge.h3k27me3, aes_string(x = "ord", y = "corr.pears.log")) + facet_wrap(~compare, scales = "free_x") + 
  geom_bar(stat = "identity") + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_x_discrete(labels = cors.merge.h3k27me3[, setNames(as.character(Sample), ord)]) + 
  ggtitle(jmark) + xlab("Cluster")
print(m.h3k27me3)

cors.merge.h3k27me3[, ord := sprintf("%02i", frank(cors.merge.h3k27me3, "compare", -"corr.spear", ties.method = "first"))]
m.h3k27me3 <- ggplot(cors.merge.h3k27me3, aes_string(x = "ord", y = "corr.spear")) + facet_wrap(~compare, scales = "free_x") + 
  geom_bar(stat = "identity") + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_x_discrete(labels = cors.merge.h3k27me3[, setNames(as.character(Sample), ord)]) + 
  ggtitle(jmark) + xlab("Cluster")
print(m.h3k27me3)


cors.merge.h3k27me3[, ord := sprintf("%02i", frank(cors.merge.h3k27me3, "compare", -"corr.spear.log", ties.method = "first"))]
m.h3k27me3 <- ggplot(cors.merge.h3k27me3, aes_string(x = "ord", y = "corr.spear.log")) + facet_wrap(~compare, scales = "free_x") + 
  geom_bar(stat = "identity") + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_x_discrete(labels = cors.merge.h3k27me3[, setNames(as.character(Sample), ord)]) + 
  ggtitle(jmark) + xlab("Cluster")
print(m.h3k27me3) 



jmark <- "H3K4me1"
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
out.lst.h3k4me1 <- lapply(names(infs.h3k4me1), function(jname) ComparePublicLinear(infs.h3k4me1[[jname]], thres=jthres, lab = jname))

# reame clusters
keys <- as.character(seq(11))
vals <- c("NeutProgProg", "Eryth", "BcellProg", "NeutProg", "Neut", "NeutIsland", "Mega", "Tcell", "NK", "HSC?", "Bcell")
vals <- paste(vals, seq(length(vals)), sep = "_")
kval.lst <- c(vals)
names(kval.lst) <- keys
clstr.hash.h3k4me1 <- hash(kval.lst)

cors.merge <- purrr::reduce(.x = lapply(out.lst.h3k4me1, function(x) x$dat.cors), .f = bind_rows) %>%
  rowwise %>%
  mutate(Sample = strsplit(strsplit(as.character(Sample), "_")[[1]][[3]], "\\.")[[1]][[1]])

cors.merge.h3k4me1 <- cors.merge %>%
  rowwise() %>%
  mutate(Sample = clstr.hash.h3k4me1[[as.character(Sample)]],
         ctype = strsplit(compare, "_")[[1]][[1]],
         in.ido.amit = ifelse(ctype %in% la.ctypes$ctype, TRUE, FALSE)) %>%
  as.data.table() 


cors.merge.h3k4me1[, ord := sprintf("%02i", frank(cors.merge.h3k4me1, compare, -corr.spear, ties.method = "first"))]
m.h3k4me1 <- ggplot(cors.merge.h3k4me1 %>% filter(grepl("rep1", compare)), aes(x = ord, y = corr.spear)) + facet_wrap(~compare, scales = "free_x") + 
  geom_bar(stat = "identity") + theme_bw(8) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_x_discrete(labels = cors.merge.h3k4me1[, setNames(as.character(Sample), ord)]) + 
  ggtitle(jmark) + xlab("Cluster")
print(m.h3k4me1)


cors.merge.h3k4me1[, ord := sprintf("%02i", frank(cors.merge.h3k4me1, compare, -corr.pears, ties.method = "first"))]
m.h3k4me1 <- ggplot(cors.merge.h3k4me1 %>% filter(grepl("rep1", compare)), aes(x = ord, y = corr.pears)) + facet_wrap(~compare, scales = "free_x") + 
  geom_bar(stat = "identity") + theme_bw(8) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_x_discrete(labels = cors.merge.h3k4me1[, setNames(as.character(Sample), ord)]) + 
  ggtitle(jmark) + xlab("Cluster")
print(m.h3k4me1)

# filter for ido amit study
cors.merge.h3k4me1[, ord := sprintf("%02i", frank(cors.merge.h3k4me1, compare, -corr.spear, ties.method = "first"))]
m.h3k4me1.filt <- ggplot(cors.merge.h3k4me1 %>% filter(in.ido.amit), aes(x = ord, y = corr.spear)) + facet_wrap(~compare, scales = "free_x") + 
  geom_bar(stat = "identity") + theme_bw(8) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_x_discrete(labels = cors.merge.h3k4me1[, setNames(as.character(Sample), ord)]) + 
  ggtitle(jmark) + xlab("Cluster")
print(m.h3k4me1.filt)

cors.merge.h3k4me1[, ord := sprintf("%02i", frank(cors.merge.h3k4me1, compare, -corr.pears, ties.method = "first"))]
m.h3k4me1.filt.pears <- ggplot(cors.merge.h3k4me1 %>% filter(in.ido.amit), aes(x = ord, y = corr.pears)) + facet_wrap(~compare, scales = "free_x") + 
  geom_bar(stat = "identity") + theme_bw(8) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_x_discrete(labels = cors.merge.h3k4me1[, setNames(as.character(Sample), ord)]) + 
  ggtitle(jmark) + xlab("Cluster")
print(m.h3k4me1.filt.pears)

# Compare H3K4me3  --------------------------------------------------------


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
# vals <- c("Bcell", "NeutroProg", "BcellProg", "Eryth", "Neutro", "HSC")
vals <- c("HSC?", "Eryth", "Bcell", "Neut", "Mega", "NK", "NeutProg", "BcellProg", "NeutProgProg")
vals <- paste(vals, seq(length(vals)), sep = "_")
kval.lst <- c(vals)
names(kval.lst) <- keys
clstr.hash.h3k4me3 <- hash(kval.lst)

cors.merge <- purrr::reduce(.x = lapply(out.lst.h3k4me3, function(x) x$dat.cors), .f = bind_rows) %>%
  rowwise %>%
  mutate(Sample = strsplit(strsplit(as.character(Sample), "_")[[1]][[3]], "\\.")[[1]][[1]])

cors.merge.h3k4me3 <- cors.merge %>%
  rowwise() %>%
  mutate(Sample = clstr.hash.h3k4me3[[as.character(Sample)]],
         ctype = strsplit(compare, "_")[[1]][[1]],
         in.ido.amit = ifelse(ctype %in% la.ctypes$ctype, TRUE, FALSE)) %>%
  as.data.table() 


# pdf("/tmp/h3k4me3.pdf", useDingbats = FALSE)

cors.merge.h3k4me3[, ord := sprintf("%02i", frank(cors.merge.h3k4me3, compare, -corr.spear, ties.method = "first"))]
m.h3k4me3.spear <- ggplot(cors.merge.h3k4me3 %>% filter(grepl("rep1", compare)), aes(x = ord, y = corr.spear)) + facet_wrap(~compare, scales = "free_x") + 
  geom_bar(stat = "identity") + theme_bw(8) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_x_discrete(labels = cors.merge.h3k4me3[, setNames(as.character(Sample), ord)]) + 
  ggtitle(jmark) + xlab("Cluster")
print(m.h3k4me3.spear)

cors.merge.h3k4me3[, ord := sprintf("%02i", frank(cors.merge.h3k4me3, compare, -corr.pears, ties.method = "first"))]
m.h3k4me3.pears <- ggplot(cors.merge.h3k4me3 %>% filter(grepl("rep1", compare)), aes(x = ord, y = corr.pears)) + facet_wrap(~compare, scales = "free_x") + 
  geom_bar(stat = "identity") + theme_bw(8) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_x_discrete(labels = cors.merge.h3k4me3[, setNames(as.character(Sample), ord)]) + 
  ggtitle(jmark) + xlab("Cluster")
print(m.h3k4me3.pears)

cors.merge.h3k4me3[, ord := sprintf("%02i", frank(cors.merge.h3k4me3, compare, -corr.spear, ties.method = "first"))]
m.h3k4me3.spear.filt <- ggplot(cors.merge.h3k4me3 %>% filter(in.ido.amit), aes(x = ord, y = corr.spear)) + facet_wrap(~compare, scales = "free_x") + 
  geom_bar(stat = "identity") + theme_bw(8) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_x_discrete(labels = cors.merge.h3k4me3[, setNames(as.character(Sample), ord)]) + 
  ggtitle(jmark) + xlab("Cluster")
print(m.h3k4me3.spear.filt)


cors.merge.h3k4me3[, ord := sprintf("%02i", frank(cors.merge.h3k4me3, compare, -corr.pears, ties.method = "first"))]
m.h3k4me3.pears.filt <- ggplot(cors.merge.h3k4me3 %>% filter(in.ido.amit), aes(x = ord, y = corr.pears)) + facet_wrap(~compare, scales = "free_x") + 
  geom_bar(stat = "identity") + theme_bw(8) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_x_discrete(labels = cors.merge.h3k4me3[, setNames(as.character(Sample), ord)]) + 
  ggtitle(jmark) + xlab("Cluster")
print(m.h3k4me3.pears.filt)


# SummarizE H3K9me3 in pseudobulk by comapring across marks ---------------

# Summarize H3K9me3 using pseudobulk correlations 
# for Bcell, Neutrophil, and Erythryoblasts 
# compare H3K9me3 

jmark <- "H3K9me3"
bsize <- "100000"
# inf2 <- "/Users/yeung/data/scchic/public_data/bigwig_compares_output_tables/H3K9me3_Multimark_H3K4me1_vs_H3K9me3_comparison.tsv"
inf.h3k9me3 <- list("11" = file.path(dirmain, paste0(jmark, "_Multimark_H3K4me1_vs_", jmark, "_comparison_binsize-", bsize, ".txt")),
                    "5" = file.path(dirmain, paste0(jmark, "_Multimark_H3K4me1_vs_", jmark, "_comparison_binsize-", bsize, ".txt")),
                    "2" = file.path(dirmain, paste0(jmark, "_Multimark_H3K4me1_vs_", jmark, "_comparison_binsize-", bsize, ".txt")))

lapply(inf.h3k9me3, function(x) assertthat::assert_that(file.exists(x)))

# clster to name
labs <- paste("H3K4me1", sapply(names(inf.h3k9me3), function(x) clstr.hash.h3k4me1[[x]]), sep = " ")
# out.lst.multi <- lapply(names(inf.h3k9me3), function(jname) CompareAcrossMarks(inf.h3k9me3[[jname]], jmark, clstr = as.numeric(jname), 
#                                                                                thres=0.95, lab = jname))

out.lst.multi.h3k9me3 <- mapply(function(jname, jlab) CompareAcrossMarks(inf.h3k9me3[[jname]], jmark, refmark, clstr = as.numeric(jname), thres = jthres, lab = jlab), 
                                names(inf.h3k9me3), labs, SIMPLIFY = FALSE)


# rename clusters in Sample to match H3K9em3
nclst.h3k9me3 <- 7
keys <- as.character(seq(nclst.h3k9me3))
# vals <- c("Bcell", "BcellProg", "Neutro", "Erythro", "NeutroProg", "HSC?")
vals <- c("NeutroProg", "Center", "Bcell", "Neutro", "BcellProg", "NeutroProgProg", "Eryth")
kval.lst <- c(vals)
names(kval.lst) <- keys
clstr.hash.multi.h3k9me3 <- hash(kval.lst)


cors.merge.multi.h3k9me3 <- bind_rows(out.lst.multi.h3k9me3[[1]]$dat.cors, out.lst.multi.h3k9me3[[2]]$dat.cors, out.lst.multi.h3k9me3[[3]]$dat.cors) %>%
  rowwise() %>%
  mutate(Sample = strsplit(strsplit(as.character(Sample), "_")[[1]][[3]], "\\.")[[1]][[1]]) %>%
  mutate(Sample = clstr.hash.multi.h3k9me3[[as.character(Sample)]])

cors.merge.multi.h3k9me3 <- cors.merge.multi.h3k9me3 %>%
  group_by(compare) %>%
  do(OrderDecreasing(., jfactor = "Sample", jval = "corr.spear.log2"))
# 
# cors.merge.multi.h3k9me3.bcell <- cors.merge.multi.h3k9me3 %>%
#   filter(compare == "H3K4me1 Bcell_3") %>%
#   do(OrderDecreasing(., jfactor = "Sample", jval = "corr.spear.log2"))
# 
# cors.merge.multi.h3k9me3.eryth <- cors.merge.multi.h3k9me3 %>%
#   filter(compare == "H3K4me1 Eryth_4") %>%
#   do(OrderDecreasing(., jfactor = "Sample", jval = "corr.spear.log2"))
# 
# cors.merge.multi.h3k9me3.neut <- cors.merge.multi.h3k9me3 %>%
#   filter(compare == "H3K4me1 Neutro_8") %>%
#   do(OrderDecreasing(., jfactor = "Sample", jval = "corr.spear.log2"))

m.multi.h3k9me3 <- ggplot(cors.merge.multi.h3k9me3, aes(x = Sample, y = corr.spear.log2)) + facet_wrap(~compare) + geom_bar(stat = "identity") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle(jmark)
print(m.multi.h3k9me3)
# 
# 
# m.multi.h3k9me3.bcell <- ggplot(cors.merge.multi.h3k9me3.bcell, aes(x = Sample, y = corr.spear.log2)) + facet_wrap(~compare) + geom_bar(stat = "identity") + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
#   ggtitle(jmark)
# 
# m.multi.h3k9me3.eryth <- ggplot(cors.merge.multi.h3k9me3.eryth, aes(x = Sample, y = corr.spear.log2)) + facet_wrap(~compare) + geom_bar(stat = "identity") + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
#   ggtitle(jmark)
# 
# m.multi.h3k9me3.neut <- ggplot(cors.merge.multi.h3k9me3.neut, aes(x = Sample, y = corr.spear.log2)) + facet_wrap(~compare) + geom_bar(stat = "identity") + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
#   ggtitle(jmark)
# 
# multiplot(m.multi.h3k9me3.bcell, m.multi.h3k9me3.eryth, m.multi.h3k9me3.neut, cols = 3)


# Cmprae H3k27me3 vs H3K4me1 ----------------------------------------------


# compare H3K27me3 muulti 

jmark <- "H3K27me3"

# for Bcell, Neutrophil, and Erythryoblasts 

bsize <- "100000"
# inf2 <- "/Users/yeung/data/scchic/public_data/bigwig_compares_output_tables/H3K9me3_Multimark_H3K4me1_vs_H3K9me3_comparison.tsv"
inf.multi.h3k27me3 <- list("11" = file.path(dirmain, paste0(jmark, "_Multimark_H3K4me1_vs_", jmark, "_comparison_binsize-", bsize, ".txt")),
                           "5" = file.path(dirmain, paste0(jmark, "_Multimark_H3K4me1_vs_", jmark, "_comparison_binsize-", bsize, ".txt")),
                           "2" = file.path(dirmain, paste0(jmark, "_Multimark_H3K4me1_vs_", jmark, "_comparison_binsize-", bsize, ".txt")))


# out.lst.multi.h3k27me3 <- lapply(names(inf.multi.h3k27me3), function(jname) CompareAcrossMarks(inf.multi.h3k27me3[[jname]], jmark, clstr = as.numeric(jname), 
# thres=0.95, lab = jname))
labs <- paste("H3K4me1", sapply(names(inf.h3k9me3), function(x) clstr.hash.h3k4me1[[x]]))
out.lst.multi.h3k27me3 <- mapply(function(jname, jlab) CompareAcrossMarks(inf.multi.h3k27me3[[jname]], jmark, refmark, clstr = as.numeric(jname), thres = jthres, lab = jlab), names(inf.multi.h3k27me3), labs, SIMPLIFY = FALSE)


cors.merge.multi.h3k27me3 <- bind_rows(out.lst.multi.h3k27me3[[1]]$dat.cors, out.lst.multi.h3k27me3[[2]]$dat.cors, out.lst.multi.h3k27me3[[3]]$dat.cors) %>%
  rowwise() %>%
  mutate(Sample = strsplit(strsplit(as.character(Sample), "_")[[1]][[3]], "\\.")[[1]][[1]]) %>%
  mutate(Sample = clstr.hash.h3k27me3[[as.character(Sample)]])

cors.merge.multi.h3k27me3 <- cors.merge.multi.h3k27me3 %>%
  group_by(compare) %>%
  do(OrderDecreasing(., jfactor = "Sample", jval = "corr.spear.log2"))


m.multi.h3k27me3 <- ggplot(cors.merge.multi.h3k27me3, aes(x = Sample, y = corr.spear.log2)) + facet_wrap(~compare) + geom_bar(stat = "identity") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle(jmark)
print(m.multi.h3k27me3)



if (plot.to.file){
  dev.off()
}
