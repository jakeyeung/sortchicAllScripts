# Jake Yeung
# Date of Creation: 2020-06-23
# File: ~/projects/scchic/scripts/rstudioserver_analysis/WKM_and_BM_together/LDA_TSS_downstream.R
# 


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
library(JFuncs)

RenameClusterBM <- function(clstr.orig, bm.rename){
  # clstr.orig <- "Bcells-Cd83_topic10"
  clstr.new <- paste0("z", clstr.orig)
  for (cname in names(bm.rename)){
    if (startsWith(clstr.orig, prefix = cname)){
      clstr.new <- bm.rename[[cname]]
    } else{
    }
  }
  return(clstr.new)
}



FitClstVar <- function(jrow, jmeta){
  fit.input <- data.frame(exprs = jrow, xvar = jmeta$xvar, clst = jmeta$clst, stringsAsFactors = TRUE)
  jfit <- lm(exprs ~ xvar:clst + clst, fit.input)
  return(jfit)
}


jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")


bm.rename <- as.list(hash(c("Bcells", "Eryth", "HSCs", "Neutrophils"), c("lymph", "eryth", "HSPCs", "granu")))
# Load data ---------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"

jdate <- "2020-06-23"
jspecies <- "MouseBM"
jdist <- 10000
jmark <- "H3K4me3"

out.objs <- lapply(jmarks, function(jmark){
  indir <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisTSS_BM_WKM_dists/lda_TSS.", jspecies, ".", jmark, ".TSSdist_", jdist, ".", jdate, ".K-30.binarize.FALSE"))
  assertthat::assert_that(dir.exists(indir))
  fname <- paste0("ldaOut.", jspecies, ".", jmark, ".TSSdist_", jdist, ".", jdate, ".K-30.Robj")
  inf <- file.path(indir, fname)
  load(inf, v=T)
  return(list(out.lda = out.lda, count.mat = count.mat))
})

out.ldas <- lapply(out.objs, function(x) x$out.lda)
count.mat <- lapply(out.objs, function(x) x$count.mat)



# Plot output -------------------------------------------------------------

tm.result.lst <- lapply(out.ldas, function(out.lda){
  tm.result <- posterior(out.lda)
  tm.result <- AddTopicToTmResult(tm.result, jsep = "")
})

dat.umaps.lst <- lapply(tm.result.lst, function(tm.result){
  dat.umap <- DoUmapAndLouvain(tm.result$topics, jsettings) %>%
    rowwise() %>%
    mutate(cond = ifelse(grepl("cd41", cell, ignore.case = TRUE), "zCd41Enriched", "Unenriched"))
  return(dat.umap)
})

m.umaps.lst <- lapply(jmarks, function(jmark){
  dat.umap <- dat.umaps.lst[[jmark]]
  m <- ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain)) + geom_point(alpha = 0.3) + theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    ggtitle(paste(jmark, "dist:", jdist))
})

JFuncs::multiplot(m.umaps.lst$H3K4me1, m.umaps.lst$H3K4me3, m.umaps.lst$H3K27me3, cols = 3)


dat.imputed.lst <- lapply(tm.result.lst, function(tm.result){
  dat.imput <- t(log(tm.result$topics %*% tm.result$terms))  # natural log links with Poisson regression 
})

dat.annot.lst.BM <- lapply(jmarks, function(jmark){
  inf.annot <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2.celltyping/GLMPCA_celltyping.", jmark, ".AllMerged.mergesize_1000.nbins_1000.penalty_1.covar_ncuts.var.log2.CenteredAndScaled.RData")
  assertthat::assert_that(file.exists(inf.annot))
  load(inf.annot, v=T)
  dat.umap.glm.fillNAs <- subset(dat.umap.glm.fillNAs, !is.na(cluster))
  dat.umap.glm.fillNAs$cluster <- sapply(dat.umap.glm.fillNAs$cluster, RenameClusterBM, bm.rename)
  return(dat.umap.glm.fillNAs)
})

dat.umap.merged <- lapply(jmarks, function(jmark){
  left_join(dat.umaps.lst[[jmark]], subset(dat.annot.lst.BM[[jmark]], select = c(cell, cluster, cluster.orig)))
})

dat.vars.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  dat.impute <- dat.imputed.lst[[jmark]]
  dat.tmp <- CalculateVarAll(dat.impute, jchromos = jchromos) %>%
    left_join(., dat.umap.merged[[jmark]])
})

ggplot(dat.annot.lst.BM$H3K27me3, aes(x = umap1, y = umap2)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Load fits ---------------------------------------------------------------

jmark <- "H3K27me3"
jctype <- "HSPCs"
hubpath <- "/home/jyeung/hub_oudenaarden"
inf.fits <- file.path(hubprefix, paste0("jyeung/data/scChiC/cut_distances/", jmark, "-BM_AllMerged.", jctype, ".sorted/fit.csv"))
jfits <- fread(inf.fits)

inf.counts <- file.path(hubprefix, paste0("jyeung/data/scChiC/cut_distances/", jmark, "-BM_AllMerged.", jctype, ".sorted/counts.csv"))
jcounts <- fread(inf.counts)
jcounts$V1 <- NULL


plot(density(unlist(jcounts)))



maxcounts <- 3
jcounts.filt <- as.matrix(jcounts)
# remove some cells
print(dim(jcounts.filt))
jcounts.filt <- jcounts.filt[, which(colSums(jcounts.filt) > 100)]
print(dim(jcounts.filt))


cmeans <- colMeans(jcounts.filt[20:nrow(jcounts.filt), ])
# cmeans <- colMeans(jcounts.filt[1:nrow(jcounts.filt), ])
jcounts.filt <- sweep(jcounts.filt, MARGIN = 2, STATS = cmeans, FUN = "/")

jcounts.filt[which(jcounts.filt > maxcounts, arr.ind = TRUE)] <- maxcounts

jcounts.filt <- jcounts.filt[, which(!is.nan(colSums(jcounts.filt)))]

# plot(density(unlist(jcounts.filt)))
# heatmap3::heatmap3(t(jcounts.filt), Rowv = NULL, Colv = NA, scale = "none")

# order by variance
jcells <- colnames(jcounts.filt)
dat.vars.filt <- subset(dat.vars.lst[[jmark]], cell %in% jcells) %>%
  arrange(cell.var.within.sum.norm)

dat.ncuts <- data.frame(ncuts = colSums(count.mat[[jmark]]), cell = colnames(count.mat[[jmark]]), stringsAsFactors = FALSE) %>%
  arrange(ncuts) %>%
  filter(cell %in% jcells)

# jcells.ordered <- dat.vars.filt$cell
jcells.ordered <- dat.ncuts$cell


colvec <- viridis::viridis(100, direction = 1)[as.numeric(cut(dat.vars.filt$cell.var.within.sum.norm, breaks = 100))]

# jcounts.filt <- log2(jcounts.filt + 1)

heatmap3::heatmap3(t(jcounts.filt[, jcells.ordered]), Rowv = NA, Colv = NA, scale = "none", RowSideColors = colvec, RowSideLabs = "TotalCuts", main = paste(jmark, jctype), margins = c(5, 20), col = viridis::viridis(50, end = 1))

heatmap3::heatmap3(t(jcounts.filt[, jcells.ordered]), Rowv = NULL, Colv = NA, scale = "none", RowSideColors = colvec, RowSideLabs = "TotalCuts", main = paste(jmark, jctype), margin = c(5, 20), col = viridis::viridis(50, end = 1), method = 'ward.D')

# 
# # write cell ordering
# # outf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/cut_distances/H3K27me3-BM_AllMerged.HSPCs.sorted/ordered_count.csv"
# # fwrite(jcounts[, ..jcells.ordered], file = outf, quote = FALSE, sep = ",")
# 
# 
# # Load cut distances ------------------------------------------------------
# 
# # refit?
# plot(unlist(jcounts[20:nrow(jcounts), 1]), pch = 20)
# plot(unlist(jcounts[20:nrow(jcounts), 8]), pch = 20)
# 
# plot(unlist(jcounts.filt[, 1]), pch = 20)
# 
# library(zoo)
# 
# jcell <- "PZ-ChIC-B6BMSC-H3K27me3-G1-2_299"
# x <- unlist(jcounts[20:nrow(jcounts), ..jcell])
# xsmooth <- zoo::rollapply(x, width = 50, FUN = mean, align = "left")
# plot(xsmooth, type = "l")
# # points(x, pch = 20)
# 
# # plot(log(x + 1), pch=20)
# 
# # Fit GLM -----------------------------------------------------------------
# 
# input.dat <- data.frame(counts = x, x = 20:nrow(jcounts), stringsAsFactors = FALSE)
# jfit <- glm(formula = counts ~ cos(w * x) + sin(x), data = input.dat)




