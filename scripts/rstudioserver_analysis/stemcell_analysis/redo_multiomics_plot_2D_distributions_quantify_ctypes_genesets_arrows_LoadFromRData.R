# Jake Yeung
# Date of Creation: 2020-05-05
# File: ~/projects/scchic/scripts/rstudioserver_analysis/stemcell_analysis/redo_multiomics_plot_2D_distributions_quantify_ctypes.R
# 
# https://stackoverflow.com/questions/4126326/how-to-quickly-form-groups-quartiles-deciles-etc-by-ordering-columns-in-a



rm(list=ls())

library(ggrastr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)

library(preprocessCore)

library(mixtools)

library(scchicFuncs)
library(JFuncs)

CalculateGeometricMedian <- function(x1, x2, x3, cname.out, cnames = c("H3K4me1", "H3K4me3", "H3K27me3")){
  # create matrix
  X <- as.matrix(data.frame(x1 = x1, x2 = x2, x3 = x3, stringsAsFactors = FALSE))
  colnames(X) <- cnames
  gmed.out <- Gmedian::Gmedian(X)
  colnames(gmed.out) <- cnames
  gmed.out <- as.data.frame(gmed.out)
  return(gmed.out[[cname.out]])
}

# Constants ---------------------------------------------------------------

make.plots <- TRUE

# jmeth <- "GeometricMedian"
jmeth <- "MarginalMedian"

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks
cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")


outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/stemcell_analysis/plot_2D_distributions_fit_threshold_signal_to_noise"
outname.rdata <- paste0("BM_2D_active_vs_repressed_fit_thresholds.", Sys.Date(), ".RData")
outname.pdf <- paste0("BM_2D_active_vs_repressed_fit_thresholds.WithGeneSets.ArrowMeth.", jmeth, ".", Sys.Date(), ".pdf")
outname.txt <- paste0("BM_2D_active_vs_repressed_fit_thresholds.WithGeneSets.ArrowMeth.", jmeth, ".", Sys.Date(), ".txt")
outf.rdata <- file.path(outdir, outname.rdata)
outf.txt <- file.path(outdir, outname.txt)

if (make.plots){
  outf <- file.path(outdir, outname.pdf)
  pdf(outf, useDingbats = FALSE)
}

# Load objects from previous analysis  ------------------------------------


if (!file.exists(outf.rdata)){
  print(paste(" rdata does not exist exists:", outf.rdata))
} else {
  load(outf.rdata, v=T)
}

jmerged.s2n.fc <- reshape2::dcast(as.data.frame(jlong.thres.lst %>% bind_rows()), formula = "bin + cluster ~ mark", value.var = "s2n.fc")
jmerged.s2n.sqrt <- reshape2::dcast(as.data.frame(jlong.thres.lst %>% bind_rows()), formula = "bin + cluster ~ mark", value.var = "s2n.diffsqrt")
jmerged.s2n.lin <- reshape2::dcast(as.data.frame(jlong.thres.lst %>% bind_rows()), formula = "bin + cluster ~ mark", value.var = "s2n.diff")

jmerged.counts <- reshape2::dcast(as.data.frame(jlong.thres.lst %>% bind_rows()), formula = "bin + cluster ~ mark", value.var = "counts")


# Analyze "flow" in different gene sets -------------------------------------

# load gene sets
# mouse BM
inf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/rdata_robjs/integrated_analysis_3_marks.stringentDE/integrated_analysis_3_marks.stringentDE.WithHighLowExprs.pseudobulk.fewerk27me3_TRUE.DSreads_TRUEDScells.FALSE.2020-05-04.RData"
load(inf, v=T)
jlong.diff.genesets.BM <- jlong.diff.genesets
print(unique(jlong.diff.genesets.BM$cluster))
jlong.diff.genesets.BM$cluster <- factor(jlong.diff.genesets.BM$cluster, c("HSPCs", "Granulocytes", "Bcells", "Erythroblasts"))
jlong.diff.genesets.BM$mark <- factor(jlong.diff.genesets.BM$mark, jmarks)
jlong.diff.genesets.BM$geneset <- factor(jlong.diff.genesets.BM$geneset, c("HSCs", "Neutrophil", "Bcell", "Erythroblast", "HighExprs", "LowExprs", "zOther"))
gsetfilt.BM <- c("HSCs", "Neutrophil", "Bcell", "Erythroblast")

head(jlong.diff.genesets)

b2e.dat <- subset(jlong.diff.genesets.BM, select = c(bin, ens, geneset))
b2e.dat <- b2e.dat[!duplicated(b2e.dat), ]

# first show the two end states

jlong.thres.annot <- as.data.frame(jlong.thres.lst %>% bind_rows()) %>%
  left_join(., b2e.dat) %>%
  filter(!is.na(ens))

jlong.thres.annot.diff <- jlong.thres.annot %>%
  group_by(bin, mark) %>%
  mutate(log2p1counts = log2(counts + 1),
         log2fc = log2p1counts - mean(log2p1counts),
         zscore = log2fc / (sd(log2p1counts)))
jlong.thres.annot.diff$cluster <- factor(jlong.thres.annot.diff$cluster, levels = c("HSPCs", "Granulocytes", "Bcells", "Erythroblasts"))
jlong.thres.annot.diff$mark <- factor(jlong.thres.annot.diff$mark, levels = jmarks)

# Saniy check with boxplots -----------------------------------------------

# matches
ggplot(jlong.thres.annot.diff, aes(x = mark, y = log2p1counts, fill = cluster)) + 
  facet_wrap(~geneset, ncol = 4) +
  geom_boxplot() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  scale_fill_manual(values = cbPalette)

ggplot(jlong.thres.annot.diff, aes(x = mark, y = counts, fill = cluster)) + 
  facet_wrap(~geneset, ncol = 4) +
  geom_boxplot() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  scale_fill_manual(values = cbPalette)

ggplot(jlong.thres.annot.diff, aes(x = counts, fill = cluster)) + 
  facet_grid(mark~geneset) +
  geom_density(alpha = 0.4) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  scale_fill_manual(values = cbPalette)


# Now look at flows -------------------------------------------------------

print(colnames(jlong.thres.annot))
# jmerged.wide <- reshape2::dcast(jlong.thres.annot, formula = "bin + ens + geneset ~ cluster + mark", value.var = "s2n.fc")

# pint(colnames(jmerge))
jmerged.s2n.annot <- left_join(jmerged.s2n.fc, b2e.dat)  %>%
  filter(!is.na(ens))

jmerged.s2n.annot$cluster <- factor(jmerged.s2n.annot$cluster, levels = c("HSPCs", "Granulocytes", "Bcells", "Erythroblasts"))
# jmerged.s2n.annot$geneset <- c("HSPCs", "Granulocytes", "Bcells", "Erythroblasts")

head(jmerged.s2n.annot)

# should be identical 
m1 <- ggplot(jmerged.s2n.annot, aes(x = H3K4me1, fill = cluster)) + facet_wrap(~geneset, nrow = 2) +  
  geom_density(alpha = 0.25) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_fill_manual(values = cbPalette) + ggtitle("log2 fold change over background (log2 signal to noise)")

m2 <- ggplot(jlong.thres.annot.diff %>% filter(mark == "H3K4me1"), aes(x = log2(counts + 1), fill = cluster)) + 
  facet_wrap(~geneset, nrow = 2) +
  geom_density(alpha = 0.4) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  scale_fill_manual(values = cbPalette) + ggtitle("log2 counts (background may vary across celltypes)")

multiplot(m2, m1, cols = 1)
 
# Sanity checks  ----------------------------------------------------------

# plot dot plots with contours
ggplot(jmerged.s2n.annot, aes(x = H3K4me3, y = H3K27me3, color = geneset)) + 
  facet_grid(geneset~cluster) + geom_point(alpha = 0.25) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_vline(xintercept = 0, linetype = "dotted") + geom_hline(yintercept = 0, linetype = "dotted") + 
  geom_density_2d(color = 'grey85') + ggtitle("log2 signal to noise. Four pseudobulks (x) across genesets (y)")

ggplot(jmerged.s2n.annot, aes(x = H3K4me1, y = H3K27me3, color = geneset)) + 
  facet_grid(geneset~cluster) + geom_point(alpha = 0.25) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_vline(xintercept = 0, linetype = "dotted") + geom_hline(yintercept = 0, linetype = "dotted") + 
  geom_density_2d(color = 'grey85') + ggtitle("log2 signal to noise 2D distribution")

# plot marginals 
ggplot(jmerged.s2n.annot, aes(x = H3K4me3, fill = geneset)) + 
  facet_grid(geneset~cluster) + geom_density() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_vline(xintercept = 0, linetype = "dotted")  + ggtitle("log2 signal to noise marginal distributions")

ggplot(jmerged.s2n.annot, aes(x = H3K4me1, fill = geneset)) + 
  facet_grid(geneset~cluster) + geom_density() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_vline(xintercept = 0, linetype = "dotted")  + ggtitle("log2 signal to noise marginal distributions")

ggplot(jmerged.s2n.annot, aes(x = H3K27me3, fill = geneset)) + 
  facet_grid(geneset~cluster) + geom_density() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_vline(xintercept = 0, linetype = "dotted")  + ggtitle("log2 signal to noised marginal distributions")


# marginals as boxplots
ggplot(jmerged.s2n.annot, aes(y = H3K4me3, fill = cluster, x = geneset)) + 
  geom_boxplot() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_vline(xintercept = 0, linetype = "dotted")  + ggtitle("log2 signal to noise marginal distributions as boxplot")

ggplot(jmerged.s2n.annot, aes(y = H3K4me1, fill = cluster, x = geneset)) + 
  geom_boxplot() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_vline(xintercept = 0, linetype = "dotted") + ggtitle("log2 signal to noise marginal distributions as boxplot")

ggplot(jmerged.s2n.annot, aes(y = H3K27me3, fill = cluster, x = geneset)) + 
  geom_boxplot() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_vline(xintercept = 0, linetype = "dotted") + ggtitle("log2 signal to noise marginal distributions as boxplot")



# Proceed -----------------------------------------------------------------



if (jmeth == "GeometricMedian"){
  # calculate geometric median
  # jsub <- subset(jmerged.s2n.annot, select = c(H3K4me1, H3K4me3, H3K27me3, cluster, geneset))
  # jsub.split <- split(x = jsub, f = jsub$cluster)  # cluster specific median
  # jsub.split.split <- lapply(jsub.split, function(jsub.tmp) split(x = jsub.tmp, f = jsub.tmp$geneset))
  # 
  # jclsts.lst <- names(jsub.split); names(jclsts.lst) <- jclsts.lst
  # jgenesets.lst <- as.character(unique(jmerged.s2n.annot$geneset)); names(jgenesets.lst) <- jgenesets.lst
  # 
  # gmed.out.lst.lst <- lapply(jclsts.lst, function(jclst){
  #   print(jclst)
  #   gmed.out.tmp.lst <- lapply(jgenesets.lst, function(jgset){
  #     print(jgset)
  #     jmat.tmp <- as.matrix(subset(jsub.split.split[[jclst]][[jgset]], select = c(H3K4me1, H3K4me3, H3K27me3)))
  #     gmed.out.tmp <- Gmedian::Gmedian(jmat.tmp)
  #     colnames(gmed.out.tmp) <- c("H3K4me1", "H3K4me3", "H3K27me3")
  #     gmed.out.tmp <- as.data.frame(gmed.out.tmp)
  #     gmed.out.tmp$cluster <- jclst
  #     return(gmed.out.tmp)
  #   })
  # })
  # 
  # jmerged.s2n.annot.sum <- jmerged.s2n.annot %>%
  #   group_by(geneset, cluster) %>%
  #   summarise(H3K4me3 = gmed.out.lst.lst[[unique(cluster)]][[unique(geneset)]]$H3K4me3,
  #             H3K27me3 = gmed.out.lst.lst[[unique(cluster)]][[unique(geneset)]]$H3K27me3,
  #             H3K4me1 = gmed.out.lst.lst[[unique(cluster)]][[unique(geneset)]]$H3K4me1,
  #             H3K4me1.start = 0,
  #             H3K4me3.start = 0,
  #             H3K27me3.start = 0)
  
  jmerged.s2n.annot.sum <- jmerged.s2n.annot %>%
    group_by(geneset, cluster) %>%
    summarise(H3K4me3 = CalculateGeometricMedian(x1 = H3K4me1, x2 = H3K4me3, x3 = H3K27me3, cname.out = "H3K4me3", cnames = c("H3K4me1", "H3K4me3", "H3K27me3")),
              H3K27me3 = CalculateGeometricMedian(x1 = H3K4me1, x2 = H3K4me3, x3 = H3K27me3, cname.out = "H3K27me3", cnames = c("H3K4me1", "H3K4me3", "H3K27me3")),
              H3K4me1 = CalculateGeometricMedian(x1 = H3K4me1, x2 = H3K4me3, x3 = H3K27me3, cname.out = "H3K4me1", cnames = c("H3K4me1", "H3K4me3", "H3K27me3")),
              H3K4me1.start = 0,
              H3K4me3.start = 0,
              H3K27me3.start = 0)
  
  jmerged.s2n.annot.sum.check2 <- jmerged.s2n.annot %>%
    group_by(geneset, cluster) %>%
    summarise(H3K4me3 = median(H3K4me3),
              H3K27me3 = median(H3K27me3),
              H3K4me1 = median(H3K4me1),
              H3K4me1.start = 0,
              H3K4me3.start = 0,
              H3K27me3.start = 0)
  
} else if (jmeth == "MarginalMedian"){
  # do average for now, simple 
  jmerged.s2n.annot.sum <- jmerged.s2n.annot %>%
    group_by(geneset, cluster) %>%
    summarise(H3K4me3 = median(H3K4me3),
              H3K27me3 = median(H3K27me3),
              H3K4me1 = median(H3K4me1),
              H3K4me1.start = 0,
              H3K4me3.start = 0,
              H3K27me3.start = 0)
} else {
  stop("jmeth must be GeometricMean or MarginalMedian, found, ", jmeth)
}


ggplot(jmerged.s2n.annot.sum, aes(x = H3K4me3, y = H3K27me3, color = geneset)) + 
  geom_point(aes(x = H3K4me3, y = H3K27me3), data = jmerged.s2n.annot, alpha = 0.1) + 
  geom_segment(aes(xend = H3K4me3.start, yend = H3K27me3.start), arrow = arrow(length=unit(0.10,"cm"), ends = "first"), color = "black") +
  geom_vline(xintercept = 0, linetype = "dotted") + geom_hline(yintercept = 0, linetype = "dotted") +  
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_grid(geneset ~ cluster) + ggtitle(paste("log2 signal to noise, arrow pointing to:", jmeth))


ggplot(jmerged.s2n.annot.sum, aes(x = H3K4me1, y = H3K27me3, color = geneset)) + 
  geom_point(aes(x = H3K4me1, y = H3K27me3), data = jmerged.s2n.annot, alpha = 0.1) + 
  geom_segment(aes(xend = H3K4me1.start, yend = H3K27me3.start), arrow = arrow(length=unit(0.10,"cm"), ends = "first"), color = "black") +
  geom_vline(xintercept = 0, linetype = "dotted") + geom_hline(yintercept = 0, linetype = "dotted") +  
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_grid(geneset ~ cluster) + ggtitle(paste("log2 signal to noise, arrow pointing to:", jmeth))

# 
# print(colnames(jmerged.wide))
# # do subtraction arrows to show "flow" from HSPCs to granulocytes for example 
# ggplot(jmerged.wide, aes(color = geneset)) + 
#   geom_segment(aes(x = HSPCs_H3K4me3, y = HSPCs_H3K27me3, xend = Granulocytes_H3K4me3, yend = Granulocytes_H3K27me3), arrow = arrow(length=unit(0.10,"cm"), ends = "first"), alpha = 0.25) + 
#   # geom_segment(aes(x = 0, y = 0, xend = HSPCs_H3K4me3, yend = HSPCs_H3K27me3), arrow = arrow(length=unit(0.10,"cm"), ends = "first"), alpha = 0.25) + 
#   facet_wrap(~geneset) + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   geom_vline(xintercept = 0, linetype = "dotted") + geom_hline(yintercept = 0, linetype = "dotted")
# 
# ggplot(jmerged.wide, aes(color = geneset)) + 
#   geom_segment(aes(x = HSPCs_H3K4me3, y = HSPCs_H3K27me3, xend = Erythroblasts_H3K4me3, yend = Erythroblasts_H3K27me3), arrow = arrow(length=unit(0.10,"cm"), ends = "first"), alpha = 0.25) + 
#   facet_wrap(~geneset) + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   geom_vline(xintercept = 0, linetype = "dotted") + geom_hline(yintercept = 0, linetype = "dotted")
#   
# ggplot(jmerged.wide %>% filter(geneset %in% c("Neutrophil")), aes(color = geneset)) + 
#   geom_segment(aes(x = HSPCs_H3K4me3, y = HSPCs_H3K27me3, xend = Granulocytes_H3K4me3, yend = Granulocytes_H3K27me3), arrow = arrow(length=unit(0.10,"cm"), ends = "first"), alpha = 0.25) + 
#   # geom_segment(aes(x = 0, y = 0, xend = HSPCs_H3K4me3, yend = HSPCs_H3K27me3), arrow = arrow(length=unit(0.10,"cm"), ends = "first"), alpha = 0.25) + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   geom_vline(xintercept = 0, linetype = "dotted") + geom_hline(yintercept = 0, linetype = "dotted")


# too messy, try just HSPCs to three celltypes for each geneset with an average arrow 

arrow.ref <- subset(jmerged.s2n.annot.sum, cluster == "HSPCs")
k4me1.ref <- hash::hash(arrow.ref$geneset, arrow.ref$H3K4me1)
k4me3.ref <- hash::hash(arrow.ref$geneset, arrow.ref$H3K4me3)
k27me3.ref <- hash::hash(arrow.ref$geneset, arrow.ref$H3K27me3)

jmerged.s2n.annot.sum2 <- jmerged.s2n.annot.sum %>%
  rowwise() %>%
  mutate(geneset = as.character(geneset)) %>%
  mutate(H3K4me1.start = ifelse(cluster != "HSPCs", k4me1.ref[[geneset]], 0),
         H3K4me3.start = ifelse(cluster != "HSPCs", k4me3.ref[[geneset]], 0),
         H3K27me3.start = ifelse(cluster != "HSPCs", k27me3.ref[[geneset]], 0))
jmerged.s2n.annot.sum2$geneset <- factor(jmerged.s2n.annot.sum2$geneset, levels = levels(jmerged.s2n.annot.sum$geneset))

jclst.end <- "Granulocytes"
jgsets <- c("HSCs", "Bcell", "Erythroblast", "Neutrophil")
jclsts.end <- c("Granulocytes", "Bcells", "Erythroblasts")

for (jclst.end in jclsts.end){
  print(jclst.end)
  jclsts <- c("HSPCs", jclst.end)
  
  m.k4me3.arrows <- ggplot(mapping = aes(x = H3K4me3, y = H3K27me3, color = cluster)) + 
    geom_point(mapping = aes(x = H3K4me3, y = H3K27me3), data = subset(jmerged.s2n.annot, geneset %in% jgsets & cluster %in% jclsts), alpha = 0.1) + 
    geom_segment(mapping = aes(xend = H3K4me3.start, yend = H3K27me3.start), 
                 arrow = arrow(length=unit(0.10,"cm"), ends = "first"), 
                 color = "black",
                 data = subset(jmerged.s2n.annot.sum2, geneset %in% jgsets & cluster %in% jclsts)) +
    geom_vline(xintercept = 0, linetype = "dotted") + geom_hline(yintercept = 0, linetype = "dotted") +  
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    facet_wrap(~geneset, nrow = 1) + 
    ggtitle(paste("From HSPCs to", jclst.end, "across 4 gene sets", jmeth))
  
  m.k4me1.arrows <- ggplot(mapping = aes(x = H3K4me1, y = H3K27me3, color = cluster)) + 
    geom_point(mapping = aes(x = H3K4me1, y = H3K27me3), data = subset(jmerged.s2n.annot, geneset %in% jgsets & cluster %in% jclsts), alpha = 0.1) + 
    geom_segment(mapping = aes(xend = H3K4me1.start, yend = H3K27me3.start), 
                 arrow = arrow(length=unit(0.10,"cm"), ends = "first"), 
                 color = "black",
                 data = subset(jmerged.s2n.annot.sum2, geneset %in% jgsets & cluster %in% jclsts)) +
    geom_vline(xintercept = 0, linetype = "dotted") + geom_hline(yintercept = 0, linetype = "dotted") +  
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    facet_wrap(~geneset, nrow = 1) + 
    ggtitle(paste("From HSPCs to", jclst.end, "across 4 gene sets", jmeth))
  
  print(m.k4me3.arrows)
  print(m.k4me1.arrows)
  
}

# write genes to output let Peter check
jclsts.vec <- names(de.ens.lst)
names(jclsts.vec) <- jclsts.vec
de.ens.dat <- lapply(jclsts.vec, function(jclst){
  de.ens.tab <- data.frame(geneset = jclst, ens = de.ens.lst[[jclst]], stringsAsFactors = FALSE) %>%
    left_join(., b2e.dat)
}) %>%
  bind_rows()

fwrite(de.ens.dat, file = outf.txt, sep = "\t")

if (make.plots){
  dev.off()
}
