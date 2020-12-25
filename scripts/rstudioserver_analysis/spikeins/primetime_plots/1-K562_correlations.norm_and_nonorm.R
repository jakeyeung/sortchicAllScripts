# Jake Yeung
# Date of Creation: 2020-11-27
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/primetime_plots/1-K562_correlations.norm_and_nonorm.R
# Combine norm and no norm 



rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(GGally)


library(ggrastr)
my_facet <- function(data, mapping, ...) {
  ggplot(data, mapping) +
    geom_point_rast(alpha = 0.1, size = 0.1, color = 'grey85') +
    geom_density_2d(color = 'blue')
}

hubprefix <- "/home/jyeung/hub_oudenaarden"

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2"
outpdf1 <- file.path(indir, paste0("K562_correlations_with_chic_chip.", Sys.Date(), ".check.again.norm_nonorm.pdf"))
outpdf2 <- file.path(indir, paste0("K562_correlations_with_chic_chip.", Sys.Date(), ".again.norm_nonorm.pdf"))


# Set up norm -------------------------------------------------------------

# inf <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/compare_with_chipseq_K562/multibigwigsummary.log2inputs.1kb/K562_chipseq_vs_chic_comparison.txt")
inf <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/compare_with_chipseq_K562/multibigwigsummary.newbams.log2input_nonorm_r1only.binsizes/K562_chipseq_vs_chic_comparison.bsize_25000.txt")
dat.out <- fread(inf)
print(colnames(dat.out))
rnames <- paste(unlist(dat.out[, 1]), paste(unlist(dat.out[, 2]), unlist(dat.out[, 3]), sep = "-"), sep = ":")

# should log
x1 <- unlist(dat.out[, 8])
names(x1) <- rnames
xfilt1 <- x1[!is.nan(x1)]
plot(density(log2(xfilt1)))

# no log
x2 <- unlist(dat.out[, 4])
names(x2) <- rnames
xfilt2 <- x2[!is.nan(x2)]
plot(density(xfilt2))

common.rows <- intersect(names(xfilt1), names(xfilt2))
plot(xfilt1[common.rows], xfilt2[common.rows])



# Set up a matrix for ggally  ---------------------------------------------


cnames.keep <- grepl(".bw", colnames(dat.out))
xmat <- as.matrix(dat.out[, ..cnames.keep])

cnames.new <- gsub("'|.bw|K562_AllMerged_|.sorted.tagged.G1filt.sorted.bsize_1000", "", colnames(xmat))
cnames.new <- gsub(".merged", "_chic", cnames.new)

colnames(xmat) <- cnames.new
rownames(xmat) <- rnames

#mes.new get rows with no NaNs
xmat.filt <- xmat[complete.cases(xmat), ]

# log2 scale the chic 
cnames.chic <- grepl("chic|ChIP", colnames(xmat.filt))

xmat.filt[, cnames.chic] <- log2(xmat.filt[, cnames.chic])


# i <- 1
# j <- i + 4
# plot(xmat.filt[, i], xmat.filt[, j], main = paste(colnames(xmat.filt)[c(i, j)], collapse = ","), pch = 20)
# 
# plot(xmat.filt[, i], xmat.filt[, j], pch = 20)
# plot(xmat.filt[, i], xmat.filt[, j], pch = 20)


# Set up no norm ----------------------------------------------------------


inf.nonorm <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/compare_with_chipseq_K562/multibigwigsummary/K562_chipseq_vs_chic_comparison.txt")
dat.out.nonorm <- fread(inf.nonorm)

rnames.nonorm <- paste(unlist(dat.out.nonorm[, 1]), paste(unlist(dat.out.nonorm[, 2]), unlist(dat.out.nonorm[, 3]), sep = "-"), sep = ":")


cnames.keep.nonorm <- grepl(".bw", colnames(dat.out.nonorm))
xmat.nonorm <- as.matrix(dat.out.nonorm[, ..cnames.keep.nonorm])

cnames.new.nonorm <- gsub("'|.bw|K562_AllMerged_|.sorted.tagged.G1filt.sorted.bsize_1000|_dedup_index.bsize_1000", "", colnames(xmat.nonorm))
cnames.new.nonorm <- gsub(".merged", "_chic", cnames.new.nonorm)

colnames(xmat.nonorm) <- cnames.new.nonorm
rownames(xmat.nonorm) <- rnames.nonorm

#mes.new get rows with no NaNs
xmat.filt.nonorm <- xmat.nonorm[complete.cases(xmat.nonorm), ]

# log2 scale the chic 
cnames.tolog.nonorm <- grepl("chic|ChIP", colnames(xmat.filt.nonorm))
xmat.filt.nonorm[, cnames.tolog.nonorm] <- log2(xmat.filt.nonorm[, cnames.tolog.nonorm])


# Show correlations norm vs nonorm ----------------------------------------

jmark <- "H3K4me3"

m1 <- ggplot(as.data.frame(xmat.filt), aes_string(x = paste0(jmark, "log2Input_1kb"))) + 
  geom_density() + 
  theme_bw() + ggtitle(paste0(jmark, "log2Input_1kb")) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

m2 <- ggplot(as.data.frame(xmat.filt), aes_string(x = paste0(jmark, "_chic"))) + 
  geom_density() + 
  theme_bw() + ggtitle(paste0(jmark, "_chic")) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

m <- ggplot(as.data.frame(xmat.filt), aes_string(x = paste0(jmark, "log2Input_1kb"), y = paste0(jmark, "_chic"))) + 
  geom_point_rast(alpha = 0.2) + 
  geom_density_2d() + 
  theme_bw() + ggtitle(jmark) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m1)
print(m2)
print(m)


m1.nonorm <- ggplot(as.data.frame(xmat.filt.nonorm), aes_string(x = paste0(jmark, "_ChIP"))) + 
  geom_density() + 
  theme_bw() + ggtitle(paste0(jmark, "_ChIP")) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m1.nonorm)

m2.nonorm <- ggplot(as.data.frame(xmat.filt.nonorm), aes_string(x = paste0(jmark, "_chic"))) + 
  geom_density() + 
  theme_bw() + ggtitle(paste0(jmark, "_chic")) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

m.nonorm <- ggplot(as.data.frame(xmat.filt.nonorm), aes_string(x = paste0(jmark, "_ChIP"), y = paste0(jmark, "_chic"))) + 
  geom_point_rast(alpha = 0.2) + 
  geom_density_2d() + 
  theme_bw() + ggtitle(jmark) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m1.nonorm)
print(m2.nonorm)
print(m.nonorm)

print(m.nonorm)


# Write pdf ---------------------------------------------------------------



pdf(outpdf1, useDingbats = FALSE)
for (jmark in jmarks){
  m1 <- ggplot(as.data.frame(xmat.filt), aes_string(x = paste0(jmark, "log2Input_1kb"))) + 
    geom_density() + 
    theme_bw() + ggtitle(paste0(jmark, "log2Input_1kb")) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  m2 <- ggplot(as.data.frame(xmat.filt), aes_string(x = paste0(jmark, "_chic"))) + 
    geom_density() + 
    theme_bw() + ggtitle(paste0(jmark, "_chic")) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  m <- ggplot(as.data.frame(xmat.filt), aes_string(x = paste0(jmark, "log2Input_1kb"), y = paste0(jmark, "_chic"))) + 
    geom_point_rast(alpha = 0.2) + 
    geom_density_2d() + 
    theme_bw() + ggtitle(jmark) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m1)
  print(m2)
  print(m)
  
  m1.nonorm <- ggplot(as.data.frame(xmat.filt.nonorm), aes_string(x = paste0(jmark, "_ChIP"))) + 
    geom_density() + 
    theme_bw() + ggtitle(paste0(jmark, "_ChIP")) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m1.nonorm)
  
  m2.nonorm <- ggplot(as.data.frame(xmat.filt.nonorm), aes_string(x = paste0(jmark, "_chic"))) + 
    geom_density() + 
    theme_bw() + ggtitle(paste0(jmark, "_chic")) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  m.nonorm <- ggplot(as.data.frame(xmat.filt.nonorm), aes_string(x = paste0(jmark, "_ChIP"), y = paste0(jmark, "_chic"))) + 
    geom_point_rast(alpha = 0.2) + 
    geom_density_2d() + 
    theme_bw() + ggtitle(jmark) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m1.nonorm)
  print(m2.nonorm)
  print(m.nonorm)
}
dev.off()

# for (i in seq(ncol(xmat.filt))){
#   print(i)
#   plot(density(xmat.filt[, i]), main = paste0(colnames(xmat.filt)[[i]]), xlab = "log signal")
# }


# plot(density(xmat.filt[, 3]))
# plot(density(xmat.filt[, 4]))
# plot(density(xmat.filt[, 5]))
# plot(density(xmat.filt[, 6]))
# plot(density(xmat.filt[, 7]))
# plot(density(xmat.filt[, 8]))

# 
# # Check without input normalization ---------------------------------------
# 



print("Writing ggpairs")

pdf(outpdf2, useDingbats=FALSE)
  ggpairs(as.data.frame(xmat.filt), lower = list(continuous = my_facet)) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  ggpairs(as.data.frame(xmat.filt.nonorm), lower = list(continuous = my_facet)) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

