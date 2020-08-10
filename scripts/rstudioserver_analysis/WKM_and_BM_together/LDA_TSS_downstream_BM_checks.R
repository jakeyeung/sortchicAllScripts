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

FitClstVar <- function(jrow, jmeta){
  fit.input <- data.frame(exprs = jrow, xvar = jmeta$xvar, clst = jmeta$clst, stringsAsFactors = TRUE)
  jfit <- lm(exprs ~ xvar:clst + clst, fit.input)
  return(jfit)
}

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/var_slope_estimates"
fname <- paste0("MouseBM_log_lm_fits.", Sys.Date(), ".CleanUpEryth.RData")
fnamepdf <- paste0("MouseBM_log_lm_fits.", Sys.Date(), ".CleanUpEryth.pdf")
outf <- file.path(outdir, fname)

# pdf(file = file.path(outdir, fnamepdf), useDingbats = FALSE)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")


# Load annots -------------------------------------------------------------

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

bm.rename <- as.list(hash(c("Bcells", "Eryth", "HSCs", "Neutrophils"), c("lymph", "eryth", "HSPCs", "granu")))

dat.annot.lst.BM <- lapply(jmarks, function(jmark){
  inf.annot <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2.celltyping/GLMPCA_celltyping.", jmark, ".AllMerged.mergesize_1000.nbins_1000.penalty_1.covar_ncuts.var.log2.CenteredAndScaled.RData")
  assertthat::assert_that(file.exists(inf.annot))
  load(inf.annot, v=T)
  dat.umap.glm.fillNAs <- subset(dat.umap.glm.fillNAs, !is.na(cluster))
  dat.umap.glm.fillNAs$cluster <- sapply(dat.umap.glm.fillNAs$cluster, RenameClusterBM, bm.rename)
  return(dat.umap.glm.fillNAs)
})




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

# make umaps, connect to annots 
dat.umaps.lst <- lapply(tm.result.lst, function(tm.result){
  dat.umap <- DoUmapAndLouvain(tm.result$topics, jsettings) %>%
    rowwise() %>%
    mutate(cond = ifelse(grepl("cd41", cell, ignore.case = TRUE), "zCd41Enriched", "Unenriched"))
  return(dat.umap)
})

dat.umap.merged <- lapply(jmarks, function(jmark){
  left_join(dat.umaps.lst[[jmark]], subset(dat.annot.lst.BM[[jmark]], select = c(cell, cluster, cluster.orig)))
})

m.umaps.lst <- lapply(jmarks, function(jmark){
  dat.umap <- dat.umaps.lst[[jmark]]
  m <- ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain)) + geom_point(alpha = 0.3) + theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    ggtitle(paste(jmark, "dist:", jdist))
})

JFuncs::multiplot(m.umaps.lst$H3K4me1, m.umaps.lst$H3K4me3, m.umaps.lst$H3K27me3, cols = 3)





m.umaps.lst <- lapply(jmarks, function(jmark){
  dat.umap <- dat.umaps.lst[[jmark]]
  m <- ggplot(dat.umap, aes(x = umap1, y = umap2, color = cond)) + geom_point(alpha = 0.3) + theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    ggtitle(paste(jmark, "dist:", jdist))
})

JFuncs::multiplot(m.umaps.lst$H3K4me1, m.umaps.lst$H3K4me3, m.umaps.lst$H3K27me3, cols = 3)


# Plot imputed ------------------------------------------------------------

dat.imputed.lst <- lapply(tm.result.lst, function(tm.result){
  dat.imput <- t(log(tm.result$topics %*% tm.result$terms))  # natural log links with Poisson regression 
})


# Load RNAseq genes  ------------------------------------------------------

inf.de <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/rdata_robjs/de_genes_stringent_objects/de_genes_sorted_and_giladi.WithHouseKeepAndNotExpressed.FixExprsOther.RData"
load(inf.de, v=T)



# Plot gene sets  ---------------------------------------------------------

rnames.all <- Reduce(intersect, lapply(dat.imputed.lst, rownames))
jgenes.all <- sapply(rnames.all, function(x) strsplit(x, ";")[[1]][[2]])
jens.all <- sapply(jgenes.all, AssignHash, g2e.hash)

e2g.hash <- hash::invert(g2e.hash)

# filter by gene sets? 
rnames.bygset <- lapply(de.ens.sorted.stringent, function(ens.vec){
  jgenes <- sapply(ens.vec, AssignHash, e2g.hash)
  jgenes.filt <- jgenes.all %in% jgenes
  rnames.filt <- rnames.all[jgenes.filt]
})


# plot gene sets onto UMAP 
gset.dat <- lapply(rnames.bygset, function(rnames){
  lapply(jmarks, function(jmark){
    jmat <- dat.imputed.lst[[jmark]]
    jmat.filt <- jmat[rnames, ]
    jmat.means <- colMeans(jmat.filt)
    exprs.dat <- data.frame(exprs = jmat.means, cell = names(jmat.means), mark= jmark, stringsAsFactors = FALSE) 
    left_join(dat.umap.merged[[jmark]], exprs.dat)
  })
})
dat.vars.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  dat.impute <- dat.imputed.lst[[jmark]]
  dat.tmp <- CalculateVarAll(dat.impute, jchromos = jchromos) %>%
    left_join(., dat.umap.merged[[jmark]])
    # left_join(., subset(gset.dat[[jgset]][[jmark]] %>% dplyr::select(c(cell, exprs))))
})







jgene <- "Ednra"

jgene <- "F3"
jgene <- "Aldh3b2"
jgene <- "Mogat2"
jgene <- "Ddx60"
jgene <- "Ifnlr1"

jgene <- "Rnd1"

jgene <- "Celsr3"
jgene <- "Sox6"
jgene <- "Pax5"

jgene <- "Sgms2"

jgene <- "Lgals4"

jgene <- "Cebpb"

jgene <- "S100a7a"

jgene <- "Hbb-y"
jgene <- "Sox6"

jgene <- "Irf4"
# jgene <- "Retnlg"
# jgene <- ""

exprs.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  jmat <- dat.imputed.lst[[jmark]]
  (rnames.keep <- grep(jgene, rownames(jmat), value = TRUE))
  if (length(rnames.keep) == 1){
    exprs.vec <- jmat[rnames.keep, ]
  } else {
    exprs.vec <- colMeans(jmat[rnames.keep, ])
  }
  exprs.dat <- data.frame(exprs = exprs.vec, cell = names(exprs.vec), rname = rnames.keep, stringsAsFactors = FALSE) 
  left_join(dat.umaps.lst[[jmark]], exprs.dat)
})

# exprs.lst <- lapply(jmarks, function(jmark){
#   print(jmark)
#   jmat <- dat.imputed.lst[[jmark]]
#   exprs.vec <- colMeans(jmat)
#   exprs.dat <- data.frame(exprs = exprs.vec, cell = names(exprs.vec), stringsAsFactors = FALSE) 
#   left_join(dat.umaps.lst[[jmark]], exprs.dat)
# })

m.exprs.lst <- lapply(jmarks, function(jmark){
  dat.umap <- exprs.lst[[jmark]]
  m <- ggplot(dat.umap, aes(x = umap1, y = umap2, color = exprs)) + geom_point(alpha = 0.3) + theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    ggtitle(paste(jmark, jgene))
  if (jmark == "H3K27me3"){
    m <- m + scale_color_viridis_c(direction = 1)
  } else {
    m <- m + scale_color_viridis_c(direction = 1)
  }
})

JFuncs::multiplot(m.exprs.lst$H3K4me1, m.exprs.lst$H3K4me3, m.exprs.lst$H3K27me3, cols = 3)



# Downstream --------------------------------------------------------------



# dat.an

dat.umap.merged2 <- lapply(jmarks, function(jmark){
  jtmp <- left_join(dat.umap.merged[[jmark]], exprs.lst[[jmark]]) %>%
    left_join(., dat.vars.lst[[jmark]]) %>%
    filter(cluster %in% c("eryth", "granu", "HSPCs", "lymph")) %>%
    ungroup() %>%
    mutate(xvar = max(cell.var.within.sum.norm) - cell.var.within.sum.norm)
})

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
ggplot(dat.umap.merged2$H3K27me3, aes(x = cluster, y = exprs)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap.merged2$H3K27me3, aes(x = cell.var.within.sum.norm, y = exprs, color = cluster)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap.merged2$H3K27me3, aes(x = xvar, y = exprs, color = cluster)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap.merged2$H3K27me3, aes(x = xvar, y = exprs, color = cluster)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + facet_wrap(~cluster)

ggplot(dat.umap.merged2$H3K27me3, aes(x = cluster, y = exprs, fill = cluster)) + geom_boxplot() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap.merged2$H3K4me1, aes(x = xvar, y = exprs, color = cluster)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap.merged2$H3K27me3 %>% filter(cluster %in% c("eryth", "granu", "HSPCs", "lymph")), aes(x = cell.var.within.sum.norm, y = exprs, color = cluster)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap.merged2$H3K27me3 %>% filter(cluster %in% c("eryth", "granu", "HSPCs", "lymph")), aes(x = cell.var.within.sum.norm, y = exprs, color = cluster)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

