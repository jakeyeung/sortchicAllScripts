# Jake Yeung
# Date of Creation: 2019-04-02
# File: ~/projects/scchic/scripts/scripts_analysis/pseudobulk_analysis/plot_correlations_summarize_everything_build95_log1p.R
# log1p is better for H3K27me3 probably


rm(list=ls())

tstart <- Sys.time() 

library(ggplot2)
library(dplyr)
library(reshape2)
library(hash)
library(JFuncs)
library(data.table)

# Functions ---------------------------------------------------------------

source("scripts/Rfunctions/PseudobulkComparisons.R")

plot.to.file <- TRUE

if (plot.to.file){
  pdf("~/data/scchic/pdfs/compare_pseudo_build95_log1p.pdf", useDingbats = FALSE)
}

# Set constants -----------------------------------------------------------

jthres <- 0.995
# jthres <- 1

# identify celltypes from Lara-Astiaso
la.ctypes <- data.table::fread("/Users/yeung/projects/scchic/data/lara-astiaso_celltypes_macbook.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE, col.names = c("mark_ctype"))
la.ctypes$ctype <- sapply(la.ctypes$mark_ctype, function(x) strsplit(x, "_")[[1]][[2]])
la.ctypes$mark <- sapply(la.ctypes$mark_ctype, function(x) strsplit(x, "_")[[1]][[1]])

# Compare H3K27me3 with public data ---------------------------------------

# compare H3K27me3 data for all public data


# label clusters 1 to 7 to meaningful names
keys <- as.character(seq(7))
vals <- c("Eryth", "Bcell", "HSC?", "NeutroProgProg", "Neutro", "NeutroProg", "BcellProg")
vals <- paste(vals, seq(length(vals)), sep = "_")
# vals <- c("Eryth", "BetweenNeuAndBcells?", "HSC?", "Bcell", "Neutrophil", "NeutroProgs?", "BcellProgs?")
kval.lst <- c(vals)
names(kval.lst) <- keys
clstr.hash.h3k27me3 <- hash(kval.lst)


# LOG SPACE
jmark <- "H3K27me3"
compares_keep <- c("MatBcell", "ProB", "HSC")
compares_grepstr <- paste(compares_keep, collapse = "|")

dirmain <- "/Users/yeung/data/scchic/from_cluster/pseudobulk_comparisons/merged_softlinks_textfile_log1p"
assertthat::assert_that(dir.exists(dirmain))
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

out.lst <- lapply(names(infs.h3k27me3), function(jname) ComparePublicLog(infs.h3k27me3[[jname]], thres=jthres, lab = jname))
# summarize for each cell type and rep
cors.merge <- purrr::reduce(.x = lapply(out.lst, function(x) x$dat.cors), .f = bind_rows) %>%
  rowwise %>%
  mutate(Sample = strsplit(strsplit(as.character(Sample), "_")[[1]][[3]], "\\.")[[1]][[1]])
# rename cors.merge
cors.merge.h3k27me3 <- cors.merge %>%
  rowwise() %>%
  mutate(Sample = kval.lst[[as.character(Sample)]]) %>%
  as.data.table()

cors.merge.h3k27me3[, ord := sprintf("%02i", frank(cors.merge.h3k27me3, "compare", -"corr.spear", ties.method = "first"))]
m.h3k27me3 <- ggplot(cors.merge.h3k27me3 %>% filter(grepl(compares_grepstr, compare)), aes_string(x = "ord", y = "corr.spear")) + facet_wrap(~compare, scales = "free_x") + 
  geom_bar(stat = "identity") + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_x_discrete(labels = cors.merge.h3k27me3[, setNames(as.character(Sample), ord)]) + 
  ggtitle(jmark) + xlab("Cluster")
print(m.h3k27me3)


cors.merge.h3k27me3[, ord := sprintf("%02i", frank(cors.merge.h3k27me3, "compare", -"corr.pears", ties.method = "first"))]
m.h3k27me3 <- ggplot(cors.merge.h3k27me3 %>% filter(grepl(compares_grepstr, compare)), aes_string(x = "ord", y = "corr.pears")) + facet_wrap(~compare, scales = "free_x") + 
  geom_bar(stat = "identity") + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_x_discrete(labels = cors.merge.h3k27me3[, setNames(as.character(Sample), ord)]) + 
  ggtitle(jmark) + xlab("Cluster")
print(m.h3k27me3)



# LINEAR SPACE
jmark <- "H3K27me3"
compares_keeplin <- c("Erythrobl", "Megakar")
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

# compare linear
cors.merge.h3k27me3.linear[, ord := sprintf("%02i", frank(cors.merge.h3k27me3.linear, "compare", -"corr.pears", ties.method = "first"))]
m.h3k27me3.linear <- ggplot(cors.merge.h3k27me3.linear %>% filter(grepl(compares_grepstrlin, compare)), aes_string(x = "ord", y = "corr.pears")) + 
  facet_wrap(~compare, scales = "free_x") + 
  geom_bar(stat = "identity") + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_x_discrete(labels = cors.merge.h3k27me3.linear[, setNames(as.character(Sample), ord)]) + 
  ggtitle(jmark) + xlab("Cluster")
print(m.h3k27me3.linear)


cors.merge.h3k27me3.linear[, ord := sprintf("%02i", frank(cors.merge.h3k27me3.linear, "compare", -"corr.spear", ties.method = "first"))]
m.h3k27me3.linear.spear <- ggplot(cors.merge.h3k27me3.linear %>% filter(grepl(compares_grepstrlin, compare)), aes_string(x = "ord", y = "corr.spear")) + 
  facet_wrap(~compare, scales = "free_x") + 
  geom_bar(stat = "identity") + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_x_discrete(labels = cors.merge.h3k27me3.linear[, setNames(as.character(Sample), ord)]) + 
  ggtitle(jmark) + xlab("Cluster")
print(m.h3k27me3.linear.spear)


# Why we get pos and neg corr for Erythro? --------------------------------

ggplot(out.lst[[1]]$dat.long.filt %>% filter(chic > 0.1), aes(x = chip, y = chic)) + geom_point() + facet_wrap(~Sample) + theme_bw()  + geom_smooth(method = 'lm', se = FALSE)
ggplot(out.lst[[2]]$dat.long.filt %>% filter(chic > 0.1), aes(x = chip, y = chic)) + geom_point() + facet_wrap(~Sample) + theme_bw()  + geom_smooth(method = 'lm', se = FALSE)

# B cells?
ggplot(out.lst[[3]]$dat.long.filt %>% filter(chic < Inf & chic > 0), aes(x = chip, y = chic)) + geom_point() + facet_wrap(~Sample) + theme_bw()  + geom_smooth(method = 'lm', se = FALSE)
ggplot(out.lst[[4]]$dat.long.filt %>% filter(chic < Inf & chic > 0), aes(x = chip, y = chic)) + geom_point() + facet_wrap(~Sample) + theme_bw()  + geom_smooth(method = 'lm', se = FALSE)

ggplot(out.lst[[3]]$dat.long.filt %>% filter(chic <= quantile(chic, probs = jthres) & chip <= quantile(chip, probs = jthres) & chic >= quantile(chic, probs = 1 - jthres) & chip >= quantile(chip, probs = 1-jthres)), 
       aes(x = chip, y = chic)) + geom_point() + facet_wrap(~Sample) + theme_bw()  + geom_smooth(method = 'lm', se = FALSE)
ggplot(out.lst[[4]]$dat.long.filt, aes(x = chip, y = chic)) + geom_point() + facet_wrap(~Sample) + theme_bw()  + geom_smooth(method = 'lm', se = FALSE)

# Megakaryos???
ggplot(out.lst[[5]]$dat.long.filt %>% filter(chic > 0.1), aes(x = chip, y = chic)) + geom_point() + facet_wrap(~Sample) + theme_bw()  + geom_smooth(method = 'lm', se = FALSE)
ggplot(out.lst[[6]]$dat.long.filt %>% filter(chic > 0.1), aes(x = chip, y = chic)) + geom_point() + facet_wrap(~Sample) + theme_bw()  + geom_smooth(method = 'lm', se = FALSE)


# Run it off --------------------------------------------------------------



if (plot.to.file){
  dev.off()
}

print(Sys.time() - tstart)