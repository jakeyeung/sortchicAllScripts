# Jake Yeung
# Date of Creation: 2019-03-11
# File: ~/projects/scchic/scripts/scripts_analysis/integrate_datasets/multi_canonical_correlation_rewrite_redo.R
# Redo with own implementation 

rm(list=ls())

setwd("~/projects/scchic")

library(GGally)
library(purrr)

library(ggplot2)
library(ggrepel)
library(tidyr)
library(umap)
library(data.table)
library(dplyr)
library(hash)
library(JFuncs)
library(topicmodels)
library(scales)

library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

library(cowplot)



library(GGally)

# use Seurat v3
# devtools::install_github(repo = "satijalab/seurat", ref = "release/3.0")
library(Seurat)

source("scripts/Rfunctions/MaraDownstream.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")



# MultiCCA helper function - calculates critical value (when to stop iterating
# in the while loop)
#
# Modified from PMA package
# @references Witten, Tibshirani, and Hastie, Biostatistics 2009
# @references \url{https://github.com/cran/PMA/blob/master/R/MultiCCA.R}
#
# @param mat.list list of matrices
# @param ws vector of projection vectors
# @param num.sets number of datasets
#
# @return returns updated critical value
#
GetCrit <- function(mat.list, ws, num.sets){
  crit <- 0
  for(i in 2:num.sets){
    for(j in 1:(i-1)){
      crit <- crit + t(ws[[i]])%*%t(mat.list[[i]])%*%mat.list[[j]]%*%ws[[j]]
    }
  }
  return(crit)
}

# MultiCCA helper function - updates W
#
# Modified from PMA package
# @references Witten, Tibshirani, and Hastie, Biostatistics 2009
# @references \url{https://github.com/cran/PMA/blob/master/R/MultiCCA.R}
#
# @param mat.list list of matrices
# @param i index of current matrix
# @param num.sets number of datasets
# @param ws initial vector of projection vectors
# @param ws.final final vector of projection vectors
#
# @return returns updated w value
#
UpdateW <- function(mat.list, i, num.sets, ws, ws.final){
  tots <- 0
  for(j in (1:num.sets)[-i]){
    diagmat <- (t(ws.final[[i]])%*%t(mat.list[[i]]))%*%(mat.list[[j]]%*%ws.final[[j]])
    diagmat[row(diagmat)!=col(diagmat)] <- 0
    tots <- tots + t(mat.list[[i]])%*%(mat.list[[j]]%*%ws[[j]]) - ws.final[[i]]%*%(diagmat%*%(t(ws.final[[j]])%*%ws[[j]]))
  }
  w <- tots/l2n(tots)
  return(w)
}

# Calculates the l2-norm of a vector
#
# Modified from PMA package
# @references Witten, Tibshirani, and Hastie, Biostatistics 2009
# @references \url{https://github.com/cran/PMA/blob/master/R/PMD.R}
#
# @param vec numeric vector
#
# @return returns the l2-norm.
#
l2n <- function(vec){
  a <- sqrt(sum(vec^2))
  if(a==0){
    a <- .05
  }
  return(a)
}

# MultiCCA helper function - calculates correlation
#
# Modified from PMA package
# @references Witten, Tibshirani, and Hastie, Biostatistics 2009
# @references \url{https://github.com/cran/PMA/blob/master/R/MultiCCA.R}
#
# @param mat.list list of matrices to calculate correlation
# @param ws vector of projection vectors
# @param num.sets number of datasets
#
# @return total correlation
#
GetCors <- function(mat.list, ws, num.sets){
  cors <- 0
  for(i in 2:num.sets){
    for(j in 1:(i-1)){
      thiscor  <-  cor(mat.list[[i]]%*%ws[[i]], mat.list[[j]]%*%ws[[j]])
      if(is.na(thiscor)) thiscor <- 0
      cors <- cors + thiscor
    }
  }
  return(cors)
}

jMultiCCA <- function(mat.list, num.ccs = 10, niter = 25){
  num.sets <- length(mat.list)
  ws <- list()
  for (i in 1:num.sets){
    ws[[i]] <- irlba(mat.list[[i]], nv = num.ccs)$v[, 1:num.ccs, drop = F]
  }
  ws.init <- ws
  ws.final <- list()
  cors <- NULL
  for(i in 1:length(ws)){
    ws.final[[i]] <- matrix(0, nrow=ncol(mat.list[[i]]), ncol=num.ccs)
  }
  for (cc in 1:num.ccs){
    print(paste0("Computing CC ", cc))
    ws <- list()
    for (i in 1:length(ws.init)){
      ws[[i]] <- ws.init[[i]][, cc]
    }
    cur.iter <- 1
    crit.old <- -10
    crit <- -20
    storecrits <- NULL
    while(cur.iter <= niter && abs(crit.old - crit)/abs(crit.old) > 0.001 && crit.old !=0){
      crit.old <- crit
      crit <- GetCrit(mat.list, ws, num.sets)
      storecrits <- c(storecrits, crit)
      cur.iter <- cur.iter + 1
      for(i in 1:num.sets){
        ws[[i]] <- UpdateW(mat.list, i, num.sets, ws, ws.final)
      }
    }
    for(i in 1:length(ws)){
      ws.final[[i]][, cc] <- ws[[i]]
    }
    cors <- c(cors, GetCors(mat.list, ws, num.sets))
  }
  results <- list(ws=ws.final, ws.init=ws.init, num.sets = num.sets, cors=cors)
}


# get minimum absolute number, return with actual sign
Vectorize(SelectAbsMin <- function(x1, x2){
  return(c(x1, x2)[[which.min(c(x1, x2))]])
}, vectorize.args = c("x1", "x2"), SIMPLIFY = FALSE, USE.NAMES = TRUE)

Vectorize(SelectAbsMin2 <- function(...){
  xlst <- list(...)
  return(xlst[[which.min(unlist(xlst))]])
}, SIMPLIFY = FALSE, USE.NAMES = TRUE)


# Get dirs ----------------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3")
mdirs <- lapply(jmarks, function(jmark){
  mdir <- paste0("/Users/yeung/data/scchic/from_cluster/mara_analysis/hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_", 
                 jmark, 
                 ".filt_0.99.center_TRUE-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0-",
                 "/",
                 "hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_", jmark, ".filt_0.99.center_TRUE")
  assertthat::assert_that(dir.exists(mdir))
  return(mdir)
})
jmarks.repress <- c("H3K27me3", "H3K9me3")
mdirs.repress <- lapply(jmarks.repress, function(jmark){
  mdir <- paste0("/Users/yeung/data/scchic/from_cluster/mara_analysis/hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_", 
                 jmark, 
                 ".filt_0.99.center_TRUE-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.bugfix-",
                 "/",
                 "hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_", jmark, ".filt_0.99.center_TRUE")
  print(mdir)
  assertthat::assert_that(dir.exists(mdir))
  return(mdir)
})
jmarks.all <- c(jmarks, jmarks.repress)
names(jmarks.all) <- jmarks.all
mdirs <- c(mdirs, mdirs.repress)
names(mdirs) <- c(jmarks, jmarks.repress)
mara.outs <- lapply(mdirs, LoadMARA)
head(mara.outs[[1]]$act.long)


# CCA analysis ------------------------------------------------------------

jmark1 <- "H3K4me1"
jmark2 <- "H3K4me3"
jmark3 <- "H3K27me3"
jmark4 <- "H3K9me3"

# do on activation marks
X <- as.matrix(mara.outs[[jmark1]]$act.mat %>% dplyr::select(-motif)); rownames(X) <- mara.outs[[jmark1]]$act.mat$motif
Y <- as.matrix(mara.outs[[jmark2]]$act.mat %>% dplyr::select(-motif)); rownames(Y) <- mara.outs[[jmark2]]$act.mat$motif
Z <- as.matrix(mara.outs[[jmark3]]$act.mat %>% dplyr::select(-motif)); rownames(Z) <- mara.outs[[jmark3]]$act.mat$motif
A <- as.matrix(mara.outs[[jmark4]]$act.mat %>% dplyr::select(-motif)); rownames(A) <- mara.outs[[jmark4]]$act.mat$motif

dat <- list(X = X, Y = Y, Z = Z, A = A)

out <- jMultiCCA(dat, num.ccs = 10, niter = 25)

# analyze output???
embeds.lst <- out$ws
loadings.lst <- mapply(function(W, X) return(X %*% W), out$ws, dat, SIMPLIFY = FALSE)

# take average?? or min??
loadings.min <- matrix(mapply(SelectAbsMin2, loadings.lst[[1]], loadings.lst[[2]], loadings.lst[[3]], loadings.lst[[4]], SIMPLIFY = FALSE, USE.NAMES = TRUE), 
                       nrow = nrow(loadings.lst[[1]]), ncol = ncol(loadings.lst[[1]]), 
                       dimnames = list(rownames(loadings.lst[[1]]), colnames(loadings.lst[[2]])))

loadings.3d <- array(data = unlist(loadings.lst), c(677, 10, 4), dimnames = list((rownames(loadings.lst[[1]]))))
loadings.mean <- apply( loadings.3d , 1:2 , mean ); rownames(loadings.mean) <- rownames(loadings.3d)
loadings.median <- apply( loadings.3d , 1:2 , median ); rownames(loadings.median) <- rownames(loadings.3d)
# loadings.mean <- matrix(mapply(mean, loadings.lst[[1]], loadings.lst[[2]], loadings.lst[[3]], loadings.lst[[4]], SIMPLIFY = FALSE, USE.NAMES = TRUE), 
#                        nrow = nrow(loadings.lst[[1]]), ncol = ncol(loadings.lst[[1]]), 
#                        dimnames = list(rownames(loadings.lst[[1]]), colnames(loadings.lst[[2]])))
# 
# loadings.median <- matrix(mapply(median, loadings.lst[[1]], loadings.lst[[2]], loadings.lst[[3]], loadings.lst[[4]], SIMPLIFY = FALSE, USE.NAMES = TRUE), 
#                        nrow = nrow(loadings.lst[[1]]), ncol = ncol(loadings.lst[[1]]), 
#                        dimnames = list(rownames(loadings.lst[[1]]), colnames(loadings.lst[[2]])))

cc1 <- 1
cc2 <- 2

cc1 <- 1
cc2 <- 2
plot(loadings.min[, cc1], loadings.min[, cc2], pch = 20, main = paste(cc1, "vs", cc2))
text(loadings.min[, cc1], loadings.min[, cc2], labels = rownames(loadings.min))

cc1 <- 5
cc2 <- 6
plot(loadings.mean[, cc1], loadings.mean[, cc2], pch = 20, main = paste(cc1, "vs", cc2))
text(loadings.mean[, cc1], loadings.mean[, cc2], labels = rownames(loadings.mean))

plot(loadings.median[, cc1], loadings.median[, cc2], pch = 20, main = paste(cc1, "vs", cc2))
text(loadings.median[, cc1], loadings.median[, cc2], labels = rownames(loadings.median))


# Plot examples  ----------------------------------------------------------

nn.vec <- c(40, 35, 40, 40)
jmindist.vec <- c(0.2, 0.1, 0.2, 0.1)
jmetric <- "euclidean"
jseed=123
# custom.settings.lst <- lapply(nn.vec, function(nn) GetUmapSettings(nn, jmetric, jmindist, jseed))
custom.settings.lst <- mapply(function(nn, jmindist) GetUmapSettings(nn, jmetric, jmindist, jseed), nn.vec, jmindist.vec, SIMPLIFY = FALSE)

meanfilt <- 10

Kstr.bin <- "15_20_25_30_35"
Kstr.nobin <- "15_20_25_30"

infs.nobin <- lapply(jmarks.all, function(jmark){
  inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_", meanfilt, 
                ".cellmin_100.cellmax_500000.binarize.", "FALSE", ".no_filt",
                "/lda_out_meanfilt.BM-", jmark, 
                ".CountThres0.K-", Kstr.nobin, ".Robj")
  assertthat::assert_that(file.exists(inf))
  return(inf)
})
infs.bin <- lapply(jmarks.all, function(jmark){
  inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_", meanfilt, 
                ".cellmin_100.cellmax_500000.binarize.", "TRUE", ".no_filt",
                "/lda_out_meanfilt.BM-", jmark, 
                ".CountThres0.K-", Kstr.bin, ".Robj")
  assertthat::assert_that(file.exists(inf))
  return(inf)
})

infs <- c(infs.bin[c("H3K4me1", "H3K4me3")], infs.nobin[c("H3K27me3", "H3K9me3")])
# out.objs <- mapply(function(jmark, inf) LoadLDABins(jmark, inf=inf), jmarks.all, infs, SIMPLIFY = FALSE)

out.lda.lst <- lapply(infs, LoadLDA)

# out.lda.lst <- out.objs

umap.lda.lst <- mapply(DoUmapFromLDA, out.lda.lst, custom.settings.lst, SIMPLIFY = FALSE)

names(umap.lda.lst) <- jmarks.all
names(mara.outs) <- jmarks.all

dat.merged.lst <- lapply(jmarks.all, function(jmark) left_join(umap.lda.lst[[jmark]], mara.outs[[jmark]]$act.long))
names(dat.merged.lst) <- jmarks.all

jmotif <- "Hoxc6"

jmotif <- "Cebpb"

jmotif <- "Gata3"
jmotif <- "Gfi1"
jmotif <- "Hoxa2"
jmotif <- "Bcl3"

jmotif <- "Tal1"
jmotif <- "Zbtb16"
jmotif <- "Hoxa2"
jmotif <- "Hoxa1"
jmotif <- "Zfp110"
jmotif <- "Gfi1"
jmotif <- "Foxc1"


jmotif <- "Hoxb6"

jmotif <- "Ebf1"
jmotif <- "Ezh2"

jmotif <- "Ezh2"

jmotif <- "Ebf1"

jmotif <- "Nfatc1"
jmotif <- "Hoxa1"
jmotif <- "Zbtb16"

jmotif <- "Stat6"
m1 <- PlotMotifInUmap(jmotif, dat.merged.lst[[jmark1]], mara.outs[[jmark1]]$zscores, jmark1, jsize = 0.75)
m2 <- PlotMotifInUmap(jmotif, dat.merged.lst[[jmark2]], mara.outs[[jmark2]]$zscores, jmark2, jsize = 0.75)
m3 <- PlotMotifInUmap(jmotif, dat.merged.lst[[jmark3]], mara.outs[[jmark3]]$zscores, jmark3, jsize = 0.75)
m4 <- PlotMotifInUmap(jmotif, dat.merged.lst[[jmark4]], mara.outs[[jmark4]]$zscores, jmark4, jsize = 0.75)
multiplot(m1, m2, m3, m4, cols = 2)

cc1 <- 1
cc2 <- 2
plot(loadings.min[, cc1], loadings.min[, cc2], pch = 20, main = paste(cc1, "vs", cc2))
text(loadings.min[, cc1], loadings.min[, cc2], labels = rownames(loadings.min))

plot(loadings.mean[, cc1], loadings.min[, cc2], pch = 20, main = paste(cc1, "vs", cc2))
text(loadings.mean[, cc1], loadings.min[, cc2], labels = rownames(loadings.mean))

plot(loadings.median[, cc1], loadings.median[, cc2], pch = 20, main = paste(cc1, "vs", cc2))
text(loadings.median[, cc1], loadings.median[, cc2], labels = rownames(loadings.median))


# loads.umap <- umap(loads.min[, 1:10])
# plot(loads.umap$layout[, 1], loads.umap$layout[, 2], pch = 20)
# text(loads.umap$layout[, 1], loads.umap$layout[, 2], labels = rownames(loads.umap$layout))

# Now let’s do multiCCA ---------------------------------------------------



