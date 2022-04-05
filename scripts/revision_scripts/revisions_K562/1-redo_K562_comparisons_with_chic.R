# Jake Yeung
# Date of Creation: 2022-03-22
# File: ~/projects/scchic/scripts/revision_scripts/revisions_K562/1-redo_K562_comparisons.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(GGally)

library(heatmap3)
library(gplots)

library(ggrastr)
my_facet <- function(data, mapping, ...) {
  ggplot(data, mapping) +
    geom_point_rast(alpha = 0.1, size = 0.1, color = 'grey85') +
    geom_density_2d(color = 'blue')
}

bsizes <- c("50000")
# jmarks <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarks) <- jmarks

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks


outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_cluster/K562_new/plots"

bsize <- bsizes[[1]]

outpdf1 <- file.path(outdir, paste0("K562_correlations_with_chic_chip.", Sys.Date(), ".check.again.norm_nonorm8.heatmap.r1r2.bsize_", bsize, ".pdf"))
outtxt <- file.path(outdir, paste0("K562_correlations_with_chic_chip.", Sys.Date(), ".check.again.norm_nonorm8.heatmap.r1r2.bsize_", bsize, ".txt"))

indir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_cluster/K562_new/compare_with_chipseq_K562"
inf <- file.path(indir, paste0("K562_chipseq_vs_chic_comparison.bsize_", bsize, ".txt"))
assertthat::assert_that(file.exists(inf))

dat.out <- fread(inf)
print(colnames(dat.out))
rnames <- paste(unlist(dat.out[, 1]), paste(unlist(dat.out[, 2]), unlist(dat.out[, 3]), sep = "-"), sep = ":")

# should log
x1 <- unlist(dat.out[, 8])
names(x1) <- rnames
xfilt1 <- x1[!is.nan(x1)]
plot(density(log2(xfilt1)))

# no log
x2 <- unlist(dat.out[, 7])
names(x2) <- rnames
xfilt2 <- x2[!is.nan(x2)]
plot(density(xfilt2))

common.rows <- intersect(names(xfilt1), names(xfilt2))
plot(xfilt1[common.rows], xfilt2[common.rows])



# Set up a matrix for ggally  ---------------------------------------------


cnames.keep <- grepl(".bw", colnames(dat.out))
xmat <- as.matrix(dat.out[, ..cnames.keep])

cnames.new <- gsub("'|.bw|K562_AllMerged_|.sorted.tagged.G1filt.sorted.bsize_1000|_dedup_index.bsize_1000.mapq_60|.1000|.tagged.GoodCells.sorted.bsize_1000.mapq_60.r1only.bw|K562_", "", colnames(xmat))
cnames.new <- gsub("_plates_merged", "_chic", cnames.new)
cnames.new <- gsub("k4me1", "H3K4me1", cnames.new)
cnames.new <- gsub("k4me3", "H3K4me3", cnames.new)
cnames.new <- gsub("k27me3", "H3K27me3", cnames.new)
cnames.new <- gsub("k9me3", "H3K9me3", cnames.new)

colnames(xmat) <- cnames.new
rownames(xmat) <- rnames

#mes.new get rows with no NaNs
xmat.filt <- xmat[complete.cases(xmat), ]

# log2 scale the chic
cnames.chic <- grepl("chic|ChIP", colnames(xmat.filt))

xmat.filt[, cnames.chic] <- log2(xmat.filt[, cnames.chic])

# plot by chromos
jchromos <- paste("chr", c(seq(22), "X", "Y"), sep = "")

for (jchromo in jchromos[1:10]){
  rkeep <- grepl(jchromo, rownames(xmat.filt))
  plot(density(xmat.filt[rkeep, "H3K9me3_chic"]), main = jchromo)
}

jcheck <- data.frame(bin = rownames(xmat.filt), log2signal = xmat.filt[, "H3K9me3_chic"], stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(chromo = scchicFuncs::GetChromo(bin))



m.check <- ggplot(jcheck, aes(x = log2signal)) +
  geom_density(alpha = 0.25) +
  geom_vline(xintercept = c(-1, -0.5)) +
  theme_bw() +
  ggtitle("K562 H3K9me3 signal all chromosomes") +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")

m.check2 <- ggplot(jcheck, aes(x = log2signal, fill = chromo)) +
  geom_density(alpha = 0.25) +
  theme_bw() +
  geom_vline(xintercept = c(-1, -0.5), linetype = "dotted") +
  facet_wrap(~chromo, nrow = 3) +
  ggtitle("K562 H3K9me3 signal by chromosome") +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")

print(m.check)
print(m.check2)



# Show correlations norm vs nonorm ----------------------------------------

jmark <- "H3K9me3"
print(colnames(xmat.filt))

cnamex.nonorm <- paste0(jmark, "_ChIP")
cnamex.norm <- paste0(jmark, "_input_normalized")
cnamey <- paste0(jmark, "_chic")

m1 <- ggplot(as.data.frame(xmat.filt), aes_string(x = cnamex.nonorm)) +
  geom_density() +
  theme_bw() + ggtitle(cnamex.nonorm) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

m2 <- ggplot(as.data.frame(xmat.filt), aes_string(x = cnamey)) +
  geom_density() +
  theme_bw() + ggtitle(cnamey) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

m <- ggplot(as.data.frame(xmat.filt), aes_string(x = cnamex.nonorm, y = cnamey)) +
  geom_point_rast(alpha = 0.2) +
  geom_density_2d() +
  theme_bw() + ggtitle(paste(cnamex.nonorm, "vs", cnamey)) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

print(m1)
print(m2)
print(m)


m1.nonorm <- ggplot(as.data.frame(xmat.filt), aes_string(x = cnamex.norm)) +
  geom_density() +
  theme_bw() + ggtitle(cnamex.norm) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m1.nonorm)

m.nonorm <- ggplot(as.data.frame(xmat.filt), aes_string(x = cnamex.norm, y = cnamey)) +
  geom_point_rast(alpha = 0.2) +
  geom_density_2d() +
  theme_bw() + ggtitle(paste(cnamex.norm, "vs", cnamey)) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m1.nonorm)
print(m.nonorm)

print(colnames(xmat.filt))
cnames.final <- c("H3K4me1_ChIP", "H3K4me3_ChIP", "H3K27me3_ChIP", "H3K9me3_input_normalized", "H3K4me1_chic", "H3K4me3_chic", "H3K27me3_chic", "H3K9me3_chic")
xmat.final <- xmat.filt[, cnames.final]


pdf(outpdf1, useDingbats = FALSE)

print(colnames(xmat.final))
cnames.reorder <- c("H3K9me3_input_normalized", "H3K9me3_chic", "H3K4me1_ChIP", "H3K4me1_chic", "H3K4me3_ChIP", "H3K4me3_chic", "H3K27me3_ChIP", "H3K27me3_chic")

x <- cor(xmat.final[, cnames.reorder])
jcol <- colByValue(seq(-1, 1, length.out = 1024), col = colorRampPalette(c("navy", "white", "firebrick3"))(21), range = c(-1, 1), las = 0)
heatmap3::heatmap3(x, Rowv = NA, Colv = NA, scale = "none", col = jcol, balanceColor = TRUE)

gplots::heatmap.2(x, Rowv = NULL, Colv = NULL, dendrogram = "none", scale = "none", col = jcol, cellnote = signif(x, digits = 2), trace = "none", density.info = "none", notecol = "#000000")
gplots::heatmap.2(x, Rowv = NULL, Colv = NULL, dendrogram = "none", scale = "none", col = jcol, cellnote = signif(x, digits = 2), trace = "none", density.info = "none", notecol = "#000000", key = FALSE)


for (jmark in jmarks){
  
  cnamex.nonorm <- paste0(jmark, "_ChIP")
  cnamex.norm <- paste0(jmark, "_input_normalized")
  cnamey <- paste0(jmark, "_chic")
  
  m1 <- ggplot(as.data.frame(xmat.filt), aes_string(x = cnamex.norm)) +
    geom_density() +
    theme_bw() + ggtitle(jmark, cnamex.norm) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  m2 <- ggplot(as.data.frame(xmat.filt), aes_string(x = cnamey)) +
    geom_density() +
    theme_bw() + ggtitle(jmark, cnamey) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  m <- ggplot(as.data.frame(xmat.filt), aes_string(x = cnamex.norm, y = cnamey)) +
    geom_point_rast(size = 1, color = "#8D918E") +
    geom_density_2d(color = "#D8292D", bins = 30) +
    theme_bw() + ggtitle(jmark, paste(cnamex.norm, "vs", cnamey)) +
    theme(aspect.ratio=1)
  
  print(m1)
  print(m2)
  print(m)
  
  m1.nonorm <- ggplot(as.data.frame(xmat.filt), aes_string(x = cnamex.nonorm)) +
    geom_density() +
    theme_bw() + ggtitle(jmark, cnamex.nonorm) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m1.nonorm)
  
  m.nonorm <- ggplot(as.data.frame(xmat.filt), aes_string(x = cnamex.nonorm, y = cnamey)) +
    geom_point_rast(size = 1, color = "#8D918E") +
    geom_density_2d(color = "#D8292D", bins = 30) +
    theme_bw() + ggtitle(jmark, paste(cnamex.nonorm, "vs", cnamey)) +
    theme(aspect.ratio=1)
  print(m1.nonorm)
  print(m.nonorm)
}
dev.off()
write.csv(xmat.final, file = outtxt)

