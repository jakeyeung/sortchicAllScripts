# Jake Yeung
# Date of Creation: 2020-11-22
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/primetime_plots/1-K562_correlations.R
# Correlate bigwigs

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
outpdf1 <- file.path(indir, paste0("K562_correlations_with_chic_chip.", Sys.Date(), ".check.again.first_submission.pdf"))
outpdf2 <- file.path(indir, paste0("K562_correlations_with_chic_chip.", Sys.Date(), ".again.first_submission.pdf"))


# Set up norm -------------------------------------------------------------


inf <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/compare_with_chipseq_K562/multibigwigsummary.log2inputs.first_submission/K562_chipseq_vs_chic_comparison.first_submission.txt")
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
x2 <- unlist(dat.out[, 4])
names(x2) <- rnames
xfilt2 <- x2[!is.nan(x2)]
plot(density(xfilt2))

common.rows <- intersect(names(xfilt1), names(xfilt2))
plot(xfilt1[common.rows], xfilt2[common.rows])


# Set up a matrix for ggally  ---------------------------------------------


cnames.keep <- grepl(".bw", colnames(dat.out))
xmat <- as.matrix(dat.out[, ..cnames.keep])

cnames.new <- gsub("'|.bw|K562_AllMerged_|.sorted.tagged.G1filt.sorted.bsize_1000|.1000", "", colnames(xmat))
cnames.new <- gsub(".merged", "_chic", cnames.new)
cnames.new <- gsub("input_normalized", "chip", cnames.new)

colnames(xmat) <- cnames.new
rownames(xmat) <- rnames

#mes.new get rows with no NaNs
xmat.filt <- xmat[complete.cases(xmat), ]

# log2 scale the chic 
cnames.chic <- grepl("chic", colnames(xmat.filt))

xmat.filt[, cnames.chic] <- log2(xmat.filt[, cnames.chic])

print(colnames(xmat.filt))


pdf(outpdf1, useDingbats = FALSE)
for (jmark in jmarks){
  m1 <- ggplot(as.data.frame(xmat.filt), aes_string(x = paste0(jmark, "_chip"))) + 
    geom_density() + 
    theme_bw() + ggtitle(paste0(jmark, "log2Input_1kb")) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  m2 <- ggplot(as.data.frame(xmat.filt), aes_string(x = paste0(jmark, "_chic"))) + 
    geom_density() + 
    theme_bw() + ggtitle(paste0(jmark, "_chic")) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  m <- ggplot(as.data.frame(xmat.filt), aes_string(x = paste0(jmark, "_chip"), y = paste0(jmark, "_chic"))) + 
    geom_point_rast(alpha = 0.2) + 
    geom_density_2d() + 
    theme_bw() + ggtitle(jmark) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m1)
  print(m2)
  print(m)
  
}
dev.off()



print("Writing ggpairs")

pdf(outpdf2, useDingbats=FALSE)
  ggpairs(as.data.frame(xmat.filt), lower = list(continuous = my_facet)) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

