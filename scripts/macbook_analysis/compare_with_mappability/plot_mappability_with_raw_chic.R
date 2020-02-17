# Jake Yeung
# Date of Creation: 2020-01-30
# File: ~/projects/scchic/scripts/macbook_analysis/compare_with_mappability/plot_mappability_with_raw_chic.R
# Compare count mat to 

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(topicmodels)
library(hash)
library(igraph)
library(umap)

library(glmpca)


# Functions ---------------------------------------------------------------

CalculateInformationContent <- function(xvec, xvec.bg, binsize=50000, kernel.bandwidth = 10^6, jalpha = 0.5){
  xpos <- seq(length(xvec)) * binsize
  x.smooth <- ksmooth(x = xpos, y = xvec, kernel = "normal", bandwidth = kernel.bandwidth)
  Xinput <- cbind(x.smooth$y, xvec.bg)
  bg.weights <- calc_scale_factor(Xinput, alpha = jalpha)
  bg.ratio <- bg.weights[1] / bg.weights[2]
  fgbg.ratio <- sum(Xinput[, 1]) / sum(Xinput[, 2])
  ic <- fgbg.ratio / bg.ratio
  return(ic)
}

CalculateVarSmooth <- function(xvec, binsize=50000, ktype = "box", kernel.bandwidth = 10^6){
  xpos <- seq(length(xvec)) * binsize
  x.smooth <- ksmooth(x = xpos, y = xvec, kernel = ktype, bandwidth = kernel.bandwidth)
  var.raw.smooth <- var(x = x.smooth$y)
  return(var.raw.smooth)
}

CalculateVarBoxSmooth <- function(count.mat, merge.size = 1000, jpseudocount = 1, jscale = 10^6, calculate.ncuts = TRUE){
  if (calculate.ncuts){
    ncuts.dat <- data.frame(cell = colnames(count.mat), ncuts = colSums(count.mat), stringsAsFactors = FALSE)
  }
  xpos <- seq(nrow(count.mat))
  rnames <- rownames(count.mat)
  jcells <- colnames(count.mat)
  
  dat.smooth <- lapply(jcells, function(jcell){
    xsmooth <- ksmooth(x = xpos, y = count.mat[, jcell], kernel = "box", bandwidth = merge.size)
    out.dat <- data.frame(coord = rnames, cell = jcell, ncuts = xsmooth$y)
  }) %>%
    bind_rows() 
  
  dat.sum.bigbins.var <- dat.smooth %>%
    group_by(cell) %>%
    summarise(ncuts.var = var(log2(ncuts * jscale + jpseudocount))) %>%
    left_join(., ncuts.dat)
  return(dat.sum.bigbins.var)
}

# 
# function(count.mat, merge.size = 1000, chromo.exclude.grep = "^chrX|^chrY", jpseudocount = 1, jscale = 10^6, calculate.ncuts = TRUE){
#   if (calculate.ncuts){
#     ncuts.dat <- data.frame(cell = colnames(count.mat), ncuts = colSums(count.mat), stringsAsFactors = FALSE)
#   }
#   count.mat.filt <- sweep(count.mat, MARGIN = 2, STATS = colSums(count.mat), FUN = "/")
#   count.mat.filt.autosome <- count.mat.filt[!grepl(chromo.exclude.grep, rownames(count.mat.filt)), ]
#   bins.all <- rownames(count.mat.filt.autosome)
#   bin.name <- paste("chr", ceiling(seq(length(bins.all)) / merge.size), sep = "_")
#   jchromos.bin <- unique(bin.name)
#   bin.newname <- paste(bin.name, seq(length(bin.name)), sep = ":")
#   rownames(count.mat.filt.autosome) <- bin.newname
#   dat.sum.bigbins <- SumAcrossChromos(count.mat = count.mat.filt.autosome, jchromos.bin, mean)
#   dat.sum.bigbins.var <- dat.sum.bigbins %>%
#     group_by(cell) %>%
#     summarise(ncuts.var = var(log2(ncuts * jscale + jpseudocount))) %>%
#     left_join(., ncuts.dat)
#   return(dat.sum.bigbins.var)
# }
# 



# Settings ----------------------------------------------------------------



jchromos.auto <- paste("chr", c(seq(19)), sep = "")

inf.mappability <- "/Users/yeung/data/mappability_from_count_mat/ensembl97_gem_65_mappability.txt"

dat.mappability.filt <- fread(inf.mappability, col.names = c("chromo", "start", "end", "signal")) %>%
  rowwise() %>%
  mutate(chromo = paste("chr", chromo, sep = ""))  %>%
  rowwise() %>%
  mutate(bname = paste(chromo, paste(start, end, sep = "-"), sep = ":")) %>%
  ungroup() %>%
  filter(chromo %in% jchromos.auto) 

dat.mappability.filt <- dat.mappability.filt[gtools::mixedorder(dat.mappability.filt$bname), ]

jchromo <- "chr19"
x <- subset(dat.mappability.filt, chromo == jchromo)$signal
plot(x, type = "l", main = jchromo)


# Download count mat  -----------------------------------------------------

inf <- "/Users/yeung/data/scchic/from_cluster/LDA_outputs/ldaAnalysisBins_B6BM_All_allmarks.2020-01-12.bsizestepsize_50000_50000.NoSliding/lda_outputs.count_mat.H3K4me1.countcutoff_1000-500-1000-1000.TAcutoff_0.5.K-30.binarize.FALSE/ldaOut.count_mat.H3K4me1.countcutoff_1000-500-1000-1000.TAcutoff_0.5.K-30.Robj"
load(inf, v=T)


# merge the two count matrices --------------------------------------------


# sort rows by chromosomes
jrows.sorted <- gtools::mixedsort(rownames(count.mat))
(unique(sapply(jrows.sorted, function(x) strsplit(x, ":")[[1]][[1]])))
# filter autosomes only
chromo.names.vec <- sapply(jrows.sorted, function(x) strsplit(x, ":")[[1]][[1]])
jrows.sorted.filt <- jrows.sorted[which(chromo.names.vec %in% jchromos.auto)]

count.mat.filt <- count.mat[jrows.sorted.filt, ]


# Plot UMAP with variance  ------------------------------------------------

tm.result <- posterior(out.lda)
topics.mat <- tm.result$topics
jsettings <- umap.defaults
jsettings[["n_neighbors"]] <- 30
jsettings[["min_dist"]] <- 0.1
jsettings[["random_state"]] <- 123

dat.umap.long <- DoUmapAndLouvain(topics.mat, jsettings)

dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms))

dat.var <- CalculateVarAll(dat.impute.log, jchromos.auto)


# Calculate information signal --------------------------------------------


binsize <- 50000
bad.cell <- as.character(subset(dat.var, cell.var.within.sum.norm == min(cell.var.within.sum.norm))$cell)
# good.cell <- as.character(subset(dat.var, cell.var.within.sum.norm == max(cell.var.within.sum.norm))$cell)
random.cell <- as.character(sample(dat.var$cell, size = 1))
# random.cell <- bad.cell
cell.var <- subset(dat.var, cell == random.cell)$cell.var.within.sum.norm

x.bad <- count.mat.filt[, random.cell]
x.pos <- seq(length(x.bad)) * binsize

par(mfrow=c(2,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
plot(x.bad[2000:3000], type = "l")
plot(dat.mappability.filt$signal[2000:3000], type = "l")
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
  
# smooth signal by 1 MB (20 windows)
x.smooth <- ksmooth(x = x.pos, y = x.bad, kernel = "normal", bandwidth = 1e6)
x.smooth.bg <- ksmooth(x = x.pos, y = dat.mappability.filt$signal, kernel = "normal", bandwidth = 1e6)

plot(x.smooth.bg, type = "l")

par(mfrow=c(4,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
plot(x.pos[2500:3000], x.bad[2500:3000], type = "l")
plot(x.smooth$x[2500:3000], x.smooth$y[2500:3000], type = "l")
plot(x.smooth$x[2500:3000], dat.mappability.filt$signal[2500:3000], type = "l")
plot(x.smooth$x[2500:3000], x.smooth.bg$y[2500:3000], type = "l")
# plot(dat.mappability.filt$signal[2000:3000], type = "l")
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)


# Calculate variance ------------------------------------------------------

par(mfrow=c(2,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
plot(log2(x.smooth$y * 10^3 + 1)[500:1500], type = "l")
plot(log2(x.smooth.bg$y * 10^3 + 1)[500:1500], type = "l")
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)


system.time(
  var.smooth.raw <- apply(count.mat.filt, MARGIN = 2, FUN = function(jcol){
    CalculateVarSmooth(xvec = jcol, binsize = 50000, ktype = "box", kernel.bandwidth = 5e7)
  })
)

cell.sizes <- colSums(count.mat.filt)
names(cell.sizes) <- colnames(count.mat.filt)
cell.sizes.dat <- data.frame(cell = names(cell.sizes), cf = cell.sizes)

var.smooth.raw.dat <- data.frame(cell = names(var.smooth.raw), var.smooth.raw = var.smooth.raw, stringsAsFactors = FALSE) %>%
  left_join(., dat.var) %>%
  left_join(., cell.sizes.dat)

ggplot(var.smooth.raw.dat, aes(x = var.smooth.raw, y = cell.var.within.sum.norm)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_log10() + scale_y_log10() 

ggplot(var.smooth.raw.dat, aes(x = var.smooth.raw, y = cf)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_log10() + scale_y_log10() 

var.smooth.raw <- CalculateVarRaw(as.matrix(count.mat.filt), merge.size = 100, chromo.exclude.grep = "^chrX|^chrY", jpseudocount = 1, jscale = 10^3, calculate.ncuts = TRUE)

var.smooth.raw.dat <- var.smooth.raw %>%
  left_join(., dat.var) %>%
  left_join(., cell.sizes.dat)

ggplot(var.smooth.raw.dat, aes(y = cell.var.within.sum.norm, x = ncuts.var, color = log10(cf))) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_log10() + scale_y_log10() + scale_color_viridis_c()
ggplot(var.smooth.raw.dat, aes(y = ncuts.var, x = log10(cf), color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_log10() + scale_y_log10() + scale_color_viridis_c(direction = -1)

# check with fake data
set.seed(0)
p.bg <- dat.mappability.filt$signal / sum(dat.mappability.filt$signal)
cell.sizes <- colSums(count.mat.filt)
names(cell.sizes) <- colnames(count.mat.filt)
fake.data <- lapply(cell.sizes, function(sf) rmultinom(n = 1, size = sf, prob = p.bg))
fake.mat <- fake.data %>% bind_cols() %>% as.data.frame()
rownames(fake.mat) <- rownames(count.mat.filt)

var.smooth.raw.bg <- CalculateVarRaw(as.matrix(fake.mat), merge.size = 100, chromo.exclude.grep = "^chrX|^chrY", jpseudocount = 1, jscale = 10^3, calculate.ncuts = TRUE)

var.smooth.raw.bg.dat <- var.smooth.raw.bg %>%
  left_join(., cell.sizes.dat)

ggplot(var.smooth.raw.bg.dat, aes(y = ncuts.var, x = log10(cf))) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_log10() + scale_y_log10() + scale_color_viridis_c(direction = -1)

var.smooth.raw.merge <- rbind(var.smooth.raw.dat %>% dplyr::select(cell, ncuts.var, cf) %>% mutate(experi = "real"), 
                              var.smooth.raw.bg.dat %>% dplyr::mutate(cell = paste(cell, "fake", sep = "")) %>% dplyr::select(cell, ncuts.var, cf) %>% mutate(experi = "mappability"))


ggplot(var.smooth.raw.merge, aes(y = ncuts.var, x = cf, color = experi)) + geom_point(alpha = 0.25) + theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_log10() + scale_y_log10()  + xlab("log10(ncuts") + ylab("Raw Intrachrom Var")

ggplot(var.smooth.raw.merge, aes(y = ncuts.var, x = cf, color = experi)) + geom_point(alpha = 0.25) + theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
   scale_x_log10()  + xlab("ncuts") + ylab("Raw Intrachrom Var")


var.smooth.raw.box <- CalculateVarBoxSmooth(as.matrix(count.mat.filt), merge.size = 1000, jpseudocount = 1, jscale = 10^3, calculate.ncuts = TRUE)
# saveRDS(var.smooth.raw.box, file = paste0("/Users/yeung/data/scchic/robjs/mappability_objs/smoothing_window_1000_bins.", Sys.Date(), ".rds"))
ggplot(var.smooth.raw.box %>% left_join(., dat.var), aes(x = cell.var.within.sum.norm, y = ncuts.var)) + geom_point()  + scale_x_log10() + scale_y_log10()



# Calculate scaling factors -----------------------------------------------

library(reticulate)
source_python("~/projects/scchic/scripts/macbook_analysis/compare_with_mappability/calc_scale_factor.py")

# Xinput <- cbind(x.smooth$y, dat.mappability.filt$signal)
Xinput <- cbind(x.smooth$y, x.smooth.bg$y)
# Xinput <- sweep(Xinput, MARGIN = 2, STATS = colSums(Xinput))
# par(mfrow=c(2,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
plot(Xinput[, 1][1500:3500], type = "l", col = 'blue')
lines(Xinput[, 2][1500:3500], type = "l", col = 'red')
# par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)

(bg.weights <- calc_scale_factor(Xinput, alpha = 0.4, return_nbins = FALSE))
bg.ratio <- bg.weights[1] / bg.weights[2]
fgbg.ratio <- sum(Xinput[, 1]) / sum(Xinput[, 2])


plot(bg.weights[[1]] * Xinput[, 1][2500:3500], type = "l", col = 'blue', ylim = c(0, ceiling(max(bg.weights[[1]] * Xinput))))
lines(bg.weights[[2]] * Xinput[, 2][2500:3500], type = "l", col = 'red')

ic <- fgbg.ratio / bg.ratio
print(random.cell)
print(cell.var)
print(paste("IC:", signif(ic, digits = 2)))


total.mapab <- sum(Xinput[, 2])

sf <- total.mapab * bg.weights[1] / bg.weights[2]

par(mfrow=c(2,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
plot(Xinput[, 1][1:5000] / sf, type = "l", col = 'red')
lines(Xinput[, 2][1:5000] / total.mapab, type = "l", col = 'blue')
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)


p.bg <- dat.mappability.filt$signal / sum(dat.mappability.filt$signal)
mdev <- apply(count.mat.filt, MARGIN = 2, function(xvec) multinomial_deviance(xvec, p = p.bg))

mdev.dat <- data.frame(cell = names(mdev), dev = mdev, stringsAsFactors = FALSE) %>%
  left_join(., dat.var)

ggplot(mdev.dat, aes(x = cell.var.within.sum.norm, y = dev)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_log10() + scale_y_log10()

# Simulate data from background data --------------------------------------

set.seed(0)
cell.sizes <- colSums(count.mat.filt)
names(cell.sizes) <- colnames(count.mat.filt)
fake.data <- lapply(cell.sizes, function(sf) rmultinom(n = 1, size = sf, prob = p.bg))

fake.mat <- fake.data %>% bind_cols() %>% as.data.frame()
rownames(fake.mat) <- rownames(count.mat.filt)

# fake.mdev <- apply(fake.mat, MARGIN = 2, function(xvec) multinomial_deviance(xvec, p = p.bg))
# fake.mdev.dat <- data.frame(cell = names(fake.mdev), dev = mdev, stringsAsFactors = FALSE) 

# smooth data and compare variance???



# Calculate ic for every cell and plot  -----------------------------------



xvec.bg <- dat.mappability.filt$signal
assertthat::assert_that(all(dat.mappability.filt$bname == rownames(count.mat.filt)))

ic.vec <- apply(count.mat.filt, MARGIN = 2, FUN = function(xvec){
  ic <- CalculateInformationContent(xvec, xvec.bg, binsize = 50000, kernel.bandwidth = 3e6, jalpha = 0.5)
  return(ic)
})

#save to output takes long time
# saveRDS(ic.vec, file = paste0("/Users/yeung/data/scchic/robjs/mappability_objs/information_content_output.", Sys.Date(), ".rds"))
# plot IC vs log counts

ic.vec.dat <- data.frame(cell = names(ic.vec), ic = ic.vec, stringsAsFactors = FALSE) %>%
  left_join(., dat.var) %>%
  left_join(., cell.sizes.dat)

ggplot(ic.vec.dat, aes(x = ic, y = cell.var.within.sum.norm, color = log10(cf))) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_log10() + scale_y_log10() + scale_color_viridis_c()

ggplot(ic.vec.dat, aes(y = ic, x = log10(cf), color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_log10() + scale_y_log10() + scale_color_viridis_c(direction = -1)

ic.vec.dat <- ic.vec.dat %>%
  rowwise() %>%
  mutate(plate = ClipLast(cell, jsep = "_"))

ggplot(ic.vec.dat, aes(y = ic, x = log10(cf), color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_log10() + scale_y_log10() + scale_color_viridis_c(direction = -1) + facet_wrap(~plate)


