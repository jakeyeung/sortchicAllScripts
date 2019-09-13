# Jake Yeung
# Date of Creation: 2019-07-24
# File: ~/projects/scchic/scripts/scripts_analysis/Linneg_analysis/3-LDA_downstream_from_projection_other_marks.R
# Other baits are in different directory from H3K4me3 because of srtringent filtering used in H3K4me3
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

library(umap)

source("scripts/Rfunctions/VariabilityFunctions.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/AuxB6.R")
source("scripts/Rfunctions/PlotFunctions.R")


# Constants for umap ------------------------------------------------------


mindist <- 0.4
nneigh <- 34
randomstate <- 123
jmetric <- "euclidean"

# umap.settings <- GetUmapSettings(nn = nneigh, jmindist = mindist, jmetric = jmetric, seed = randomstate)
umap.settings <- umap.defaults
umap.settings$min_dist <- mindist; umap.settings$n_neighbors <- nneigh


# Load original data  -----------------------------------------------------


jmark <- "H3K27me3"
# jmark <- "H3K4me1"
jbin <- "FALSE"
outpdf <- paste0("/Users/yeung/data/scchic/pdfs/B6_figures/linneg/", jmark, "_enrichment.pdf")

assertthat::assert_that(!file.exists(outpdf))


# inf <- "/Users/yeung/data/scchic/from_cluster/lda_outputs_linneg_from_projection/lda_out_meanfilt.B6_H3K4me3_pcutoff_0.CountThres0.K-25_30_35_50_mindist_0.4_mindist_processed_lda.Rdata"
# inf <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6/dat_umap_long_with_louvain.H3K4me3.RData"
# inf <- "/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/bamlist_for_merging_build95_B6_stringent/dat_umap_long_with_louvain.H3K4me3.stringent_filter.RData"
# inf <- paste0("/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6_stringent/dat_umap_long_with_louvain.", ".stringent_filter.RData")
inf <- paste0("/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6/dat_umap_long_with_louvain.", jmark, ".RData")
assertthat::assert_that(file.exists(inf))
load(inf, v=T)

# umap.orig <- umap(out.objs$tm.result$topics, config = umap.settings)
umap.orig <- umap(out.objs$tm.result$topics, config = custom.settings)

umap.orig.long <- data.table(umap1 = umap.orig$layout[, 1], umap2 = -1 * umap.orig$layout[, 2], samp = rownames(umap.orig$layout), batch = "Ctrl")

ggplot(umap.orig.long, aes(x = umap1, y = umap2)) + geom_point()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

# load new data
# inf.new <- "/Users/yeung/data/scchic/from_cluster/lda_outputs_linneg_from_projection/bin_PZ-Bl6-BM-Linneg_H3K4me3_binfilt_cellfilt.2019-06-17.RData"
inf.new <- paste0("/Users/yeung/data/scchic/from_cluster/lda_outputs_linneg_from_projection/bin_", jbin, "_PZ-Bl6-BM-Linneg_", jmark, "_binfilt_cellfilt.from_louvain.RData")
load(inf.new, v=T)


umap.new.predict <- predict(umap.orig, out.lda.predict$topics)

umap.new.predict.long <- data.table(umap1 = umap.new.predict[, 1], umap2 = umap.new.predict[, 2], samp = rownames(umap.new.predict), batch = "Linneg") %>%
  mutate(umap2 = umap2 * -1)

ggplot(umap.new.predict.long, aes(x = umap1, y = umap2)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

umap.merged <- bind_rows(umap.orig.long, umap.new.predict.long) %>%
  dplyr::rename(cell = samp)


ggplot(umap.merged, aes(x = umap1, y = umap2, color = batch)) + geom_point(alpha = 0.7, size = 1.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle(jmark)



# Calculate intra chromosomal variance  -----------------------------------

# Get datmat --------------------------------------------------------------


# calculate intrachromosomal original 
jfac <- 10^6
jpseudo <- 0
dat.mat.predict <- t(out.lda.predict$terms) %*% t(out.lda.predict$topics)
# log2 transform
dat.mat.predict <- log2(dat.mat.predict * jfac + jpseudo)
cells.sd.predict <- GetCellSd(dat.mat.predict, "", log2.scale = FALSE, fn = var) %>%
  mutate(batch = "Linneg") %>%
  dplyr::select(-label)

dat.mat.orig <- t(out.objs$tm.result$terms) %*% t(out.objs$tm.result$topics)
dat.mat.orig <- log2(dat.mat.orig * jfac + jpseudo)
cells.sd.orig <- GetCellSd(dat.mat.orig, "", log2.scale = FALSE, fn = var)  %>%
  mutate(batch = "Ctrl") %>%
  dplyr::select(-label)

cells.sd.merged <- bind_rows(cells.sd.predict, cells.sd.orig)

# replot umap 
umap.merged.var <- left_join(umap.merged, cells.sd.merged)

m <- ggplot(umap.merged.var, aes(x = umap1, y = umap2, color = cell.sd)) + geom_point(size = 1.5) + facet_wrap(~batch) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle(jmark)
print(m)


# Fraction of high low variance  ------------------------------------------

m.hist <- ggplot(umap.merged.var, aes(x = cell.sd, fill = batch)) + geom_density(alpha = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab("Genome-wide Variance")
print(m.hist)


# Within chromosomes  -----------------------------------------------------

jchromos <- c(paste("chr", seq(19), sep = ""), "chrX", "chrY")
jchromos.grep <- paste(jchromos, ":", sep = "")

common.rows <- intersect(rownames(dat.mat.orig), rownames(dat.mat.predict))
# dat.merged <- list("H3K4me3" = cbind(dat.mat.orig[common.rows, ], dat.mat.predict[common.rows, ]))
dat.merged <- list("Ctrl" = dat.mat.orig,
                   "Linneg" = dat.mat.predict)
# jmarks <- list("H3K4me3" = "H3K4me3")
batches <- list("Ctrl" = "Ctrl", "Linneg" = "Linneg")

# calculate variance within chromosomes
cells.var.chromo.within <- lapply(batches, function(batch){
  cells.var.chromo <- lapply(jchromos.grep, function(jstr){
    cells.sd <- GetCellSd(dat.merged[[batch]], grep.str = jstr, log2.scale = FALSE, fn = SumSqrDev) %>%
      mutate(batch = batch)
    return(cells.sd)
  }) %>%
    bind_rows()
}) %>%
  bind_rows() %>%
  dplyr::rename(cell.var.within = cell.sd)

# summarize across chromosomes
cells.var.chromo.within.sum <- cells.var.chromo.within %>%
  group_by(cell, batch) %>%
  summarise(cell.var.within.sum = sum(cell.var.within))



# Across chromosomes  -----------------------------------------------------

# calculate variance across chromosomes
cells.chromo.mean <- lapply(batches, function(batch){
  cells.chromo.mean.tmp <- lapply(jchromos.grep, function(jstr){
    dat.mat.all <- dat.merged[[batch]]
    rows.keep <- grepl(jstr, rownames(dat.mat.all))
    dat.mat <- dat.mat.all[rows.keep, ]
    chromo.avg <- colMeans(dat.mat)
    # chromo.avg <- colSums(dat.mat)
    chromo.avg.dat <- data.frame(cell = names(chromo.avg), chromo.mean = chromo.avg, label = jstr, nbins = nrow(dat.mat))
    return(chromo.avg.dat)
  }) %>%
    bind_rows()
}) %>%
  bind_rows()

# calculate global mean
dat.mat.global.mean <- lapply(batches, function(batch){
  dat.mat.global <- dat.merged[[batch]]
  out <- colMeans(dat.mat.global)
  return(out)
}) %>%
  unlist()

# make sure cell names align with cells.chromo.mean.wide.mat
names(dat.mat.global.mean) <- sapply(names(dat.mat.global.mean), function(x) strsplit(x, "\\.")[[1]][[2]])
dat.mat.global.mean <- dat.mat.global.mean[order(names(dat.mat.global.mean))]

chromo.constant <- cells.chromo.mean %>%
  group_by(label) %>%
  summarise(nbins = unique(nbins))
# normalization factor
chromo.constant.sum <- chromo.constant %>%
  summarise(nbins = sum(nbins))

cells.var.chromo.within.sum <- cells.var.chromo.within.sum %>%
  group_by(intravar = cell.var.within.sum / chromo.constant.sum$nbins)


# Plot intra chromosomal variance  ----------------------------------------

umap.merged.var.intra <- left_join(umap.merged.var, cells.var.chromo.within.sum)

ggplot(umap.merged.var.intra, aes(x = umap1, y = umap2, color = intravar)) + geom_point(size = 5) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~batch)

ggplot(umap.merged.var.intra, aes(x = intravar, fill = batch)) + geom_density(alpha = 0.5, size = 2) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 


# normalize 
# normalization factor
chromo.constant <- cells.chromo.mean %>%
  group_by(label) %>%
  summarise(nbins = unique(nbins))

chromo.constant.sum <- chromo.constant %>%
  summarise(nbins = sum(nbins))

# Check enrichment and depletion on louvain cluster -----------------------

head(dat.umap.long)

# cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#C0C0C0", "#32CD32", "D3D3D3")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#C0C0C0", "#32CD32", "#D3D3D3")
ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point()+ 
  theme_minimal() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)

# Do umap with 2 batches --------------------------------------------------

umap.merged.louvain <- left_join(umap.merged, dat.umap.long %>% dplyr::select(cell, louvain)) %>%
  rowwise() %>%
  mutate(louvain = ifelse(batch == "Linneg", "Linneg", louvain))

ggplot(umap.merged.louvain, aes(x = umap1, y = umap2, color = louvain)) + geom_point()+ 
  theme_minimal() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette) + facet_wrap(~batch)


# Statistical side: find change in Gaussian  ------------------------------



# PDFs  -------------------------------------------------------------------


m.scatter <- ggplot(umap.merged, aes(x = umap1, y = umap2, color = batch)) + geom_point(alpha = 0.7, size = 1.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle(jmark)
m.louvain <- ggplot(umap.merged.louvain, aes(x = umap1, y = umap2, color = louvain)) + geom_point()+ 
  theme_minimal() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette) + facet_wrap(~batch)

m.var <- ggplot(umap.merged.var.intra, aes(x = umap1, y = umap2, color = intravar)) + geom_point(size = 2.5) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~batch)


ksout <- ks.test(subset(umap.merged.var.intra, batch == "Ctrl")$intravar, subset(umap.merged.var.intra, batch == "Linneg")$intravar)
ksout$p.value

m.var.density <- ggplot(umap.merged.var.intra, aes(x = intravar, fill = batch)) + geom_density(alpha = 0.5, size = 1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  ggtitle("KS test Pval < 2.2e-16")

pdf(outpdf, useDingbats = FALSE)
  print(m.scatter)
  print(m.louvain)
  print(m.var)
  print(m.var.density)
dev.off()





# plot 
# umap.merged.var <- left_join(umap.merged.var, cells.var.chromo.within.sum)




# merge the two datasets and then go

# 
# 
# # Load data ---------------------------------------------------------------
# 
# inf <- "/Users/yeung/data/scchic/from_cluster/lda_outputs_linneg/lda_out_meanfilt.PZ-Bl6-BM-Linneg_H3K4me3_binfilt_cellfilt.CountThres0.K-30_35_50_mindist_0.4_mindist_processed_lda.Rdata"
# assertthat::assert_that(file.exists(inf))
# 
# load(inf, v=T)
# 
# 
# # Plot output  ------------------------------------------------------------
# 
# dim(out.objs$tm.result$topics)
# 
# settings <- umap.defaults
# settings$n_neighbors <- 58
# settings$min_dist <- 0.4
# settings$random_state <- 123
# umap.out <- umap(out.objs$tm.result$topics, config = settings)
# 
# umap.long <- data.frame(umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2], batch = sapply(rownames(umap.out$layout), function(x) strsplit(x, split = "-")[[1]][[1]]))
# 
# ggplot(umap.long, aes(x = umap1, y = umap2, color = batch)) + geom_point() + 
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# 
# # Plot intra chromosomal variance  ----------------------------------------
# 
# 
