# Jake Yeung
# Date of Creation: 2019-03-25
# File: ~/projects/scchic/scripts/scripts_analysis/pseudobulk_analysis/plot_correlations_summarize_everything.R
# After exploring, summarize everything 

rm(list=ls())

library(ggplot2)
library(dplyr)
library(reshape2)
library(hash)
library(JFuncs)

# Functions ---------------------------------------------------------------


ComparePublicLog <- function(inf, thres = 0.995, lab = "MyLabel"){
  dat <- data.table::fread(inf, sep = "\t")
  
  datref <- dat[, 1]
  colnames(datref) <- c("chip")
  datcompare <- dat[, 2:ncol(dat)]
  
  datref$peakid <- seq(nrow(datref))
  datcompare$peakid <- seq(nrow(datcompare))
  
  # plot 1 vs all
  datcompare.long <- melt(datcompare, variable.name = "Sample", value.name = "chic", id.vars = "peakid")
  
  dat.long <- dplyr::left_join(datcompare.long, datref)
  
  dat.long.filt <- dat.long %>%
    filter(!is.nan(chip)) %>%
    group_by(Sample) %>%
    filter(chic < quantile(chic, probs = thres, na.rm = TRUE)) %>%
    mutate(compare = lab)
  
  
  dat.cors <- dat.long.filt %>%
    group_by(Sample) %>%
    summarise(corr.pears = cor(chip, chic, method = "pearson"),
              corr.spear = cor(chip, chic, method = "spearman")) %>%
    mutate(compare = lab)
  # remoev outliers?
  
  # m <- ggplot(dat, aes(x = H3K4me1_MatBcell_rep1_reltoinput.bw, y = H3K4me1_cluster_1.log1p.bw)) + geom_point() + 
  #   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  m <- ggplot(dat.long.filt, aes(x = chip, y = chic)) + 
    geom_hex(bins = 25) +
    # geom_point(alpha = 0.01) + 
    # stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~Sample)
  return(list(dat.long.filt = dat.long.filt, dat.cors = dat.cors, m = m))
}

ComparePublicLinear <- function(inf, thres = 0.995, lab = "MyLabel"){
  dat <- data.table::fread(inf, sep = "\t")
  
  datref <- dat[, 1]
  colnames(datref) <- c("chip")
  datcompare <- dat[, 2:ncol(dat)]
  
  datref$peakid <- seq(nrow(datref))
  datcompare$peakid <- seq(nrow(datcompare))
  
  # plot 1 vs all
  datcompare.long <- melt(datcompare, variable.name = "Sample", value.name = "chic", id.vars = "peakid")
  
  dat.long <- dplyr::left_join(datcompare.long, datref)
  
  # remove NaNs
  
  dat.long.filt <- dat.long %>%
    filter(!is.nan(chip)) %>%
    group_by(Sample) %>%
    filter(log(chic) < quantile(log(chic), probs = thres, na.rm = TRUE) & log(chic) > quantile(log(chic), probs = 1 - thres)) %>%
    filter(log(chip) < quantile(log(chip), probs = thres, na.rm = TRUE) & log(chip) > quantile(log(chip), probs = 1 - thres)) %>%
    mutate(compare = lab)
  
  # dat.long.filt <- subset(dat.long.filt, !is.nan(chip))
  
  dat.cors <- dat.long.filt %>%
    group_by(Sample) %>%
    summarise(corr.pears = cor(chip, chic, method = "pearson"),
              corr.spear = cor(chip, chic, method = "spearman"),
              corr.pears.log = cor(log(chip), log(chic), method = "pearson"),
              corr.spear.log = cor(log(chip), log(chic), method = "spearman")) %>%
    mutate(compare = lab)
  # remoev outliers?
  
  # m <- ggplot(dat, aes(x = H3K4me1_MatBcell_rep1_reltoinput.bw, y = H3K4me1_cluster_1.log1p.bw)) + geom_point() + 
  #   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  m <- ggplot(dat.long.filt, aes(x = log10(chip), y = log10(chic))) + 
    geom_hex(bins = 25) +
    # geom_point(alpha = 0.01) + 
    # stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~Sample)
  return(list(dat.long.filt = dat.long.filt, dat.cors = dat.cors, m = m))
}

CompareAcrossMarks <- function(inf2, jmark, clstr, thres = 0.995, lab = "MyLabel"){
  dat.multi <- data.table::fread(inf2, sep = "\t")
  # compare erythryroblasts
  datref.multi <- dat.multi[, ..clstr]
  colnames(datref.multi) <- c("chic_ref")
  cols.compare <- grepl(jmark, colnames(dat.multi))
  datcompare.multi <- dat.multi[, ..cols.compare]
  
  datref.multi$peakid <- seq(nrow(datref.multi))
  datcompare.multi$peakid <- seq(nrow(datcompare.multi))
  
  datcompare.multi.long <- melt(datcompare.multi, id.vars = "peakid", variable.name = "Sample", value.name = "chic")
  
  dat.multi.long <- dplyr::left_join(datcompare.multi.long, datref.multi)
  
  # remove outliers
  dat.multi.long.filt <- dat.multi.long %>%
    group_by(Sample) %>%
    filter(chic < quantile(chic, probs = thres) & chic > quantile(chic, probs = 1 - thres)) %>%
    filter(chic_ref < quantile(chic_ref, probs = thres) & chic_ref > quantile(chic_ref, probs = 1 - thres)) %>%
    mutate(compare = lab)
  
  m <- ggplot(dat.multi.long.filt, aes(x = chic_ref, y = chic)) + geom_point(alpha = 0.25) + 
    scale_x_log10() + scale_y_log10() + facet_wrap(~Sample) + theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  dat.cors <- dat.multi.long.filt %>%
    group_by(Sample) %>%
    summarise(corr.pears = cor(chic, chic_ref, method = "pearson"),
              corr.spear = cor(chic, chic_ref, method = "spearman"),
              corr.pears.log2 = cor(log2(chic), log2(chic_ref), method = "pearson"),
              corr.spear.log2 = cor(log2(chic), log2(chic_ref), method = "spearman")) %>%
    
    mutate(compare = lab)
  # print(dat.cors)
  return(list(dat.multi.long.filt = dat.multi.long.filt, dat.cors = dat.cors, m = m))
}


# Summarize H3K27me3 with Bcell and Neutrophil data -----------------------


# Summarize H3K27me3 in log comparisons
# for Bcell and Neutrophil

jmark <- "H3K27me3"
infs.h3k27me3 <- list("MatBcell" = paste0("/Users/yeung/data/scchic/public_data/bigwig_compares_output_tables/", jmark, "_MatBcell_rep1_comparison_InputNorm.tsv"),
                      "Neu" = paste0("/Users/yeung/data/scchic/public_data/bigwig_compares_output_tables/", jmark, "_Neu_comparison_InputNorm.tsv"))
out.lst <- lapply(names(infs.h3k27me3), function(jname) ComparePublicLog(infs.h3k27me3[[jname]], thres=0.995, lab = jname))
# summarize correlations
cors.merge <- bind_rows(out.lst[[1]]$dat.cors, out.lst[[2]]$dat.cors) %>%
  rowwise() %>%
  mutate(Sample = strsplit(strsplit(as.character(Sample), "_")[[1]][[3]], "\\.")[[1]][[1]])

# label clusters 1 to 7 to meaningful names
keys <- as.character(seq(7))
vals <- c("Eryth", "BetweenNeuAndBcells?", "HSC?", "Bcell", "Neutrophil", "NeutroProgs?", "BcellProgs?")
kval.lst <- c(vals)
names(kval.lst) <- keys
clstr.hash.h3k27me3 <- hash(kval.lst)

# rename cors.merge
cors.merge.h3k27me3 <- cors.merge %>%
  rowwise() %>%
  mutate(Sample = kval.lst[[as.character(Sample)]])

cors.merge.h3k27me3 <- cors.merge.h3k27me3 %>%
  group_by(compare) %>%
  do(OrderDecreasing(., jfactor = "Sample", jval = "corr.pears"))

m.h3k27me3 <- ggplot(cors.merge.h3k27me3, aes(x = Sample, y = corr.pears)) + facet_wrap(~compare) + 
  geom_bar(stat = "identity") + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle(jmark) + xlab("Cluster")

print(m.h3k27me3)

# split into two then plot again
cors.merge.h3k27me3.bcell <- cors.merge.h3k27me3 %>%
  filter(compare == "MatBcell") %>%
  do(OrderDecreasing(., jfactor = "Sample", jval = "corr.pears"))

cors.merge.h3k27me3.neut <- cors.merge.h3k27me3 %>%
  filter(compare == "Neu") %>%
  do(OrderDecreasing(., jfactor = "Sample", jval = "corr.pears"))

m.h3k27me3.bcell <- ggplot(cors.merge.h3k27me3.bcell, aes(x = Sample, y = corr.pears)) + facet_wrap(~compare) + 
  geom_bar(stat = "identity") + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle(jmark) + xlab("Cluster")

m.h3k27me3.neut <- ggplot(cors.merge.h3k27me3.neut, aes(x = Sample, y = corr.pears)) + facet_wrap(~compare) + 
  geom_bar(stat = "identity") + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle(jmark) + xlab("Cluster")

multiplot(m.h3k27me3.bcell, m.h3k27me3.neut, cols = 2)



# Summarize H3K4me1 and H3K4me3 with Bcell and Neutrohil public da --------


# Summarize H3K4me1, H3K4me3 in linear comparisons
# for Bcell and Neutrophil

jmark <- "H3K4me1"
infs.h3k4me1 <- list("MatBcell" = paste0("/Users/yeung/data/scchic/public_data/bigwig_compares_output_tables/", jmark, "_MatBcell_MatBcell_comparison.tsv"),
                      "Neu" = paste0("/Users/yeung/data/scchic/public_data/bigwig_compares_output_tables/", jmark, "_Neu_comparison.tsv"))
out.lst <- lapply(names(infs.h3k4me1), function(jname) ComparePublicLinear(infs.h3k4me1[[jname]], thres=0.995, lab = jname))

# reame clusters
keys <- as.character(seq(9))
vals <- c("Bcell", "BcellProg", "Neutro", "NK", "Tcells?", "Eryth", "HSCs?", "BTCellProg?", "NeutroProg")
kval.lst <- c(vals)
names(kval.lst) <- keys
clstr.hash.h3k4me1 <- hash(kval.lst)

# summarize correlations
cors.merge <- bind_rows(out.lst[[1]]$dat.cors, out.lst[[2]]$dat.cors) %>%
  rowwise() %>%
  mutate(Sample = strsplit(strsplit(as.character(Sample), "_")[[1]][[3]], "\\.")[[1]][[1]]) 

cors.merge.h3k4me1 <- cors.merge %>%
  rowwise() %>%
  mutate(Sample = clstr.hash.h3k4me1[[as.character(Sample)]])

cors.merge.h3k4me1 <- cors.merge.h3k4me1 %>%
  group_by(compare) %>%
  do(OrderDecreasing(., jfactor = "Sample", jval = "corr.pears"))

cors.merge.h3k4me1.bcell <- cors.merge.h3k4me1 %>% 
  filter(compare == "MatBcell") %>%
  do(OrderDecreasing(., jfactor = "Sample", jval = "corr.pears"))

cors.merge.h3k4me1.neut <- cors.merge.h3k4me1 %>% 
  filter(compare == "Neu") %>%
  do(OrderDecreasing(., jfactor = "Sample", jval = "corr.pears"))



m.h3k4me1 <- ggplot(cors.merge.h3k4me1, aes(x = Sample, y = corr.pears)) + facet_wrap(~compare) + 
  geom_bar(stat = "identity") + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle(jmark) + xlab("Cluster")
print(m.h3k4me1)

m.h3k4me1.bcell <- ggplot(cors.merge.h3k4me1.bcell, aes(x = Sample, y = corr.pears)) + facet_wrap(~compare) + 
  geom_bar(stat = "identity") + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle(jmark) + xlab("Cluster")
m.h3k4me1.neut <- ggplot(cors.merge.h3k4me1.neut, aes(x = Sample, y = corr.pears)) + facet_wrap(~compare) + 
  geom_bar(stat = "identity") + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle(jmark) + xlab("Cluster")

multiplot(m.h3k4me1.bcell, m.h3k4me1.neut, cols = 2)

jmark <- "H3K4me3"
infs.h3k4me3 <- list("MatBcell" = paste0("/Users/yeung/data/scchic/public_data/bigwig_compares_output_tables/", jmark, "_MatBcell_MatBcell_comparison.tsv"),
                     "Neu" = paste0("/Users/yeung/data/scchic/public_data/bigwig_compares_output_tables/", jmark, "_Neu_comparison.tsv"))
out.lst <- lapply(names(infs.h3k4me3), function(jname) ComparePublicLinear(infs.h3k4me3[[jname]], thres=0.995, lab = jname))
# summarize correlations
keys <- as.character(seq(6))
vals <- c("Bcell", "NeutroProg", "BcellProg", "BcellProgProg", "Neutro", "Eryth")
kval.lst <- c(vals)
names(kval.lst) <- keys
clstr.hash.h3k4me3 <- hash(kval.lst)

cors.merge <- bind_rows(out.lst[[1]]$dat.cors, out.lst[[2]]$dat.cors) %>%
  rowwise() %>%
  mutate(Sample = strsplit(strsplit(as.character(Sample), "_")[[1]][[3]], "\\.")[[1]][[1]])

cors.merge.h3k4me3 <- cors.merge %>%
  mutate(Sample = clstr.hash.h3k4me3[[as.character(Sample)]])

cors.merge.h3k4me3 <- cors.merge.h3k4me3 %>%
  group_by(compare) %>%
  do(OrderDecreasing(., jfactor = "Sample", jval = "corr.pears"))

m.h3k4me3 <- ggplot(cors.merge.h3k4me3, aes(x = Sample, y = corr.pears)) + facet_wrap(~compare) + 
  geom_bar(stat = "identity") + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle(jmark) + xlab("Cluster")
print(m.h3k4me3)

cors.merge.h3k4me3.bcell <- cors.merge.h3k4me3 %>%
  filter(compare == "MatBcell") %>%
  do(OrderDecreasing(., jfactor = "Sample", jval = "corr.pears"))

cors.merge.h3k4me3.neut <- cors.merge.h3k4me3 %>%
  filter(compare == "Neu") %>%
  do(OrderDecreasing(., jfactor = "Sample", jval = "corr.pears"))

m.h3k4me3.bcell <- ggplot(cors.merge.h3k4me3.bcell, aes(x = Sample, y = corr.pears)) + facet_wrap(~compare) + 
  geom_bar(stat = "identity") + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle(jmark) + xlab("Cluster")
m.h3k4me3.neut <- ggplot(cors.merge.h3k4me3.neut, aes(x = Sample, y = corr.pears)) + facet_wrap(~compare) + 
  geom_bar(stat = "identity") + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle(jmark) + xlab("Cluster")

multiplot(m.h3k4me3.bcell, m.h3k4me3.neut, cols = 2)

# SummarizE H3K9me3 in pseudobulk by comapring across marks ---------------


# Summarize H3K9me3 using pseudobulk correlations 
# for Bcell, Neutrophil, and Erythryoblasts 
# compare H3K9me3 

jmark <- "H3K9me3"
# inf2 <- "/Users/yeung/data/scchic/public_data/bigwig_compares_output_tables/H3K9me3_Multimark_H3K4me1_vs_H3K9me3_comparison.tsv"
inf.h3k9me3 <- list("3" = paste0("/Users/yeung/data/scchic/public_data/bigwig_compares_output_tables/", jmark, "_Multimark_H3K4me1_vs_", jmark, "_comparison.tsv"),
                    "1" = paste0("/Users/yeung/data/scchic/public_data/bigwig_compares_output_tables/", jmark, "_Multimark_H3K4me1_vs_", jmark, "_comparison.tsv"),
                    "6" = paste0("/Users/yeung/data/scchic/public_data/bigwig_compares_output_tables/", jmark, "_Multimark_H3K4me1_vs_", jmark, "_comparison.tsv"))

# clster to name
labs <- paste("H3K4me1", sapply(names(inf.h3k9me3), function(x) clstr.hash.h3k4me1[[x]]))
# out.lst.multi <- lapply(names(inf.h3k9me3), function(jname) CompareAcrossMarks(inf.h3k9me3[[jname]], jmark, clstr = as.numeric(jname), 
#                                                                                thres=0.95, lab = jname))

out.lst.multi.h3k9me3 <- mapply(function(jname, jlab) CompareAcrossMarks(inf.h3k9me3[[jname]], jmark, clstr = as.numeric(jname), thres = 0.95, lab = jlab), names(inf.h3k9me3), labs, SIMPLIFY = FALSE)

# rename clusters in Sample to match H3K9em3
keys <- as.character(seq(7))
vals <- c("NeutroProg", "Neutro", "Bcell", "BcellProg", "HSCs?", "Eryth", "NeutroProgProg")
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

cors.merge.multi.h3k9me3.bcell <- cors.merge.multi.h3k9me3 %>%
  filter(compare == "H3K4me1 Bcell") %>%
  do(OrderDecreasing(., jfactor = "Sample", jval = "corr.spear.log2"))

cors.merge.multi.h3k9me3.eryth <- cors.merge.multi.h3k9me3 %>%
  filter(compare == "H3K4me1 Eryth") %>%
  do(OrderDecreasing(., jfactor = "Sample", jval = "corr.spear.log2"))

cors.merge.multi.h3k9me3.neut <- cors.merge.multi.h3k9me3 %>%
  filter(compare == "H3K4me1 Neutro") %>%
  do(OrderDecreasing(., jfactor = "Sample", jval = "corr.spear.log2"))


m.multi.h3k9me3 <- ggplot(cors.merge.multi.h3k9me3, aes(x = Sample, y = corr.spear.log2)) + facet_wrap(~compare) + geom_bar(stat = "identity") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle(jmark)
print(m.multi.h3k9me3)


m.multi.h3k9me3.bcell <- ggplot(cors.merge.multi.h3k9me3.bcell, aes(x = Sample, y = corr.spear.log2)) + facet_wrap(~compare) + geom_bar(stat = "identity") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle(jmark)

m.multi.h3k9me3.eryth <- ggplot(cors.merge.multi.h3k9me3.eryth, aes(x = Sample, y = corr.spear.log2)) + facet_wrap(~compare) + geom_bar(stat = "identity") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle(jmark)

m.multi.h3k9me3.neut <- ggplot(cors.merge.multi.h3k9me3.neut, aes(x = Sample, y = corr.spear.log2)) + facet_wrap(~compare) + geom_bar(stat = "identity") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle(jmark)

multiplot(m.multi.h3k9me3.bcell, m.multi.h3k9me3.eryth, m.multi.h3k9me3.neut, cols = 3)


# Cmprae H3k27me3 vs H3K4me1 ----------------------------------------------


# compare H3K27me3 muulti 

jmark <- "H3K27me3"
# inf2 <- "/Users/yeung/data/scchic/public_data/bigwig_compares_output_tables/H3K9me3_Multimark_H3K4me1_vs_H3K9me3_comparison.tsv"
inf.multi.h3k27me3 <- list("3" = paste0("/Users/yeung/data/scchic/public_data/bigwig_compares_output_tables/", jmark, "_Multimark_H3K4me1_vs_", jmark, "_comparison.tsv"),
                           "1" = paste0("/Users/yeung/data/scchic/public_data/bigwig_compares_output_tables/", jmark, "_Multimark_H3K4me1_vs_", jmark, "_comparison.tsv"),
                           "6" = paste0("/Users/yeung/data/scchic/public_data/bigwig_compares_output_tables/", jmark, "_Multimark_H3K4me1_vs_", jmark, "_comparison.tsv"))
# out.lst.multi.h3k27me3 <- lapply(names(inf.multi.h3k27me3), function(jname) CompareAcrossMarks(inf.multi.h3k27me3[[jname]], jmark, clstr = as.numeric(jname), 
                                                                                               # thres=0.95, lab = jname))
labs <- paste("H3K4me1", sapply(names(inf.h3k9me3), function(x) clstr.hash.h3k4me1[[x]]))
out.lst.multi.h3k27me3 <- mapply(function(jname, jlab) CompareAcrossMarks(inf.multi.h3k27me3[[jname]], jmark, clstr = as.numeric(jname), thres = 0.95, lab = jlab), names(inf.multi.h3k27me3), labs, SIMPLIFY = FALSE)



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


# Plot outputs ------------------------------------------------------------

pdf(paste0("~/Dropbox/scCHiC_figs/FIG4_BM/analyses/", Sys.Date(), ".pseudobulk_public_data_comparisons.pdf"), useDingbats = FALSE)

print(m.h3k27me3)
multiplot(m.h3k27me3.bcell, m.h3k27me3.neut, cols = 2)
print(m.h3k4me1)
multiplot(m.h3k4me1.bcell, m.h3k4me1.neut, cols = 2)
print(m.h3k4me3)
multiplot(m.h3k4me3.bcell, m.h3k4me1.neut, cols = 2)
print(m.multi.h3k9me3)
multiplot(m.multi.h3k9me3.bcell, m.multi.h3k9me3.eryth, m.multi.h3k9me3.neut, cols = 3)
print(m.multi.h3k27me3)
multiplot(m.multi.h3k9me3.bcell, m.multi.h3k9me3.eryth, m.multi.h3k9me3.neut, cols = 3)


dev.off()


