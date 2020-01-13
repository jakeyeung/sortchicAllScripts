# Jake Yeung
# Date of Creation: 2020-01-13
# File: ~/projects/scchic/scripts/rstudioserver_analysis/debug/show_var_on_plate.R
# Show var on plate: BM K4me3 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(topicmodels)
library(hash)
library(igraph)
library(umap)

library(scchicFuncs)


library(platetools)
library(ggplot2)
library(viridis)

inf.well <- "~/data/from_rstudioserver/scchic/well_names.txt"
assertthat::assert_that(file.exists(inf.well))

well <- read.table(file=inf.well, header = F, stringsAsFactors = FALSE)
well$V1 <- sapply(well$V1, function(x) as.numeric(gsub("^X", "", x)))
# shift by one
well$V1 <- as.character(well$V1 - 1)
well.hash <- hash(well$V1, well$V2)

# Load data ---------------------------------------------------------------

inf <- "/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2019-12-16/lda_outputs.B6BM_AllMerged_H3K4me3.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-16.K-50.binarize.FALSE/ldaOut.B6BM_AllMerged_H3K4me3.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-16.K-50.Robj"

load(inf, v=T)

topics.mat <- posterior(out.lda)$topics
terms.mat <- posterior(out.lda)$terms

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
umap.out <- umap(topics.mat, config = jsettings)
dat.umap.long <- data.frame(cell = rownames(umap.out[["layout"]]), umap1 = umap.out[["layout"]][, 1], umap2 = umap.out[["layout"]][, 2], stringsAsFactors = FALSE)
dat.umap.long <- DoLouvain(topics.mat, jsettings, dat.umap.long)
cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette)


dat.impute.log <- log2(t(topics.mat %*% terms.mat))
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")

dat.var <- CalculateVarAll(dat.impute.log, jchromos)

dat.merge <- left_join(dat.umap.long, dat.var)

dat.merge$indx <- sapply(dat.merge$cell, function(x) strsplit(x, "_")[[1]][[2]])
dat.merge$well <- sapply(dat.merge$indx, function(x) well.hash[[x]])

# show one plate

dat.merge <- dat.merge %>%
  rowwise() %>%
  mutate(experi = ClipLast(cell, jsep = "_"))

experis <- unique(dat.merge$experi)
jprefix <- "B6-13W1-BM-H3K4me3"
experis <- experis[grepl(jprefix, experis)]

pdf(paste0("/home/jyeung/hpc/scChiC/from_rstudioserver/pdfs_all/variance_along_plates/scchicseq_", jprefix, ".", Sys.Date(), ".plate_effect_variance", ".pdf"), useDingbats = FALSE)

# check that there is plate effect?
jsub <- subset(dat.merge, grepl(jprefix, experi))
m.dens <- ggplot(jsub, aes(x = cell.var.within.sum.norm, fill = experi)) + geom_density(alpha = 0.25) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_fill_manual(values = cbPalette) + facet_wrap(~experi, ncol = 1)
print(m.dens)
jrange <- range(jsub$cell.var.within.sum.norm)
for (jexperi in experis){
  dat.sub <- subset(dat.merge, experi == jexperi)
  dat.sub$well <- sapply(dat.sub$indx, function(x) well.hash[[x]])
  print(paste(jexperi, nrow(dat.sub)))
  m <- raw_map(data = dat.sub$cell.var.within.sum.norm, well = dat.sub$well, plate = 384) +
    ggtitle("Single cell 384-well plate", jexperi) +
    theme_dark() +
    # scale_fill_viridis(direction = -1, begin = 0, end = 1) + 
    scale_fill_gradientn(colors = viridis_pal()(9), limits=jrange, 
                         na.value = "#FDE725FF")
  print(m)
}
dev.off()


