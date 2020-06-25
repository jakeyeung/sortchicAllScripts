# Jake Yeung
# Date of Creation: 2020-06-24
# File: ~/projects/scchic/scripts/rstudioserver_analysis/WKM_and_BM_together/fit_poisson_model_gene_by_gene.ZF_WithVarianceSlope.R
# Add variance slope? 



rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(hash)
library(igraph)
library(umap)

library(scchicFuncs)
library(topicmodels)


# Load pois fits ----------------------------------------------------------


inf.poisfits <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/poisson_fits.OtherNorm.DiffDists/fit_poisson_model_on_TSS.ZebrafishWKM.NormMeth_ncuts.alltss.dist_10000.RData"
load(inf.poisfits, v=T)

# Load set up  ------------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

# Load data ---------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"

jdate <- "2020-06-22"
jspecies <- "ZebrafishWKM"
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
count.mat.lst <- lapply(out.objs, function(x) x$count.mat)


# Load annots -------------------------------------------------------------


wkm.rename <- as.list(hash(c("Bcells", "Eryth", "HSCs", "Neutrophils"), c("lymph", "eryth", "HSPCs", "granu")))

dat.annot.lst.WKM <- lapply(jmarks, function(jmark){
  print(jmark)
  # filter by previously defined celltypes? 
  inf.annot.louv <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/LDA_downstream/LDA_downstream_ZF.2020-04-23.imputevarfilt.lessstringent/ZF_LDA_output.", jmark, ".keepn_150.final.ClusterTables.txt")
  assertthat::assert_that(file.exists(inf.annot.louv))
  
  
  inf.annot.glmpca <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/celltyping/from_louvain/cell_to_cluster_table.", jmark, ".2020-04-13.txt")
  assertthat::assert_that(file.exists(inf.annot.glmpca))
  
  annot.louv <- fread(inf.annot.louv)
  annot.louv$clusterplate <- paste(annot.louv$cluster, annot.louv$plate, "_")
  
  annot.glmpca <- fread(inf.annot.glmpca)
  annot.glmpca.filt <- subset(annot.glmpca, cell %in% annot.louv$cell) %>%
    rowwise() %>%
    mutate(clusterplate = paste(cluster, plate, sep = "_")) %>%
    mutate(cluster.orig = cluster) %>%
    mutate(cluster = ifelse(cluster %in% c("lymph1", "lymph2"), "lymph", cluster),
           cluster = ifelse(cluster %in% c("eryth1", "eryth2"), "eryth", cluster),
           cluster = ifelse(cluster %in% c("HSC1", "HSC2"), "HSC", cluster)) %>%   # rename lymph1 and lymph2 into lymph
    ungroup() %>%
    filter(cluster != "Unknown") %>%
    left_join(., subset(annot.louv, select = c(cell, var.imputed)))
  # rename clusters
  annot.glmpca.filt$cluster <- sapply(annot.glmpca.filt$cluster, function(jclst) AssignHash(jclst, wkm.rename, null.fill = jclst))
  print("annot glmpca filt")
  print(annot.glmpca.filt)
  return(annot.glmpca.filt)
})


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

dat.umap.merged <- lapply(jmarks, function(jmark){
  left_join(dat.umaps.lst[[jmark]], subset(dat.annot.lst.WKM[[jmark]], select = c(cell, cluster, cluster.orig)))
})



dat.imputed.lst <- lapply(tm.result.lst, function(tm.result){
  dat.imput <- t(log(tm.result$topics %*% tm.result$terms))  # natural log links with Poisson regression 
})



# Set up meta -------------------------------------------------------------


jchromos <- paste("chr", seq(25), sep = "")

clsts.keep <- c("HSC", "lymph", "monocyte", "eryth")
# fit genome-wide 
metadat.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  dat.impute <- dat.imputed.lst[[jmark]]
  dat.tmp <- CalculateVarAll(dat.impute, jchromos = jchromos) %>%
    left_join(., dat.umap.merged[[jmark]]) %>%
    filter(cluster %in% clsts.keep)
})

# Fit genes  --------------------------------------------------------------






# pax5 as a test?
# jmark <- "H3K27me3"
jmark <- "H3K27me3"
jmark <- "H3K4me3"
jmark <- "H3K4me1"

jrow.i <- grep("pax5", rownames(count.mat.lst[[jmark]]), value = TRUE)

jmat <- dat.imputed.lst[[jmark]]
jmat.raw <- count.mat.lst[[jmark]]

cnames.keep <- metadat.lst[[jmark]]$cell
jmat.filt <- jmat[, cnames.keep]
jmat.raw.filt <- jmat.raw[, cnames.keep]

jmeta <- metadat.lst[[jmark]] %>%
  dplyr::rename(xvar = cell.var.within.sum.norm, clst = cluster) %>%
  dplyr::select(cell, xvar, clst) %>%
  mutate(xvar.orig = xvar,
         xvar = max(xvar) - xvar)

jmeta.ordered <- jmeta[match(colnames(jmat.filt), jmeta$cell), ]
assertthat::assert_that(identical(jmeta.ordered$cell, colnames(jmat.filt)))
assertthat::assert_that(identical(jmeta.ordered$cell, colnames(jmat.raw.filt)))
print(unique(jmeta.ordered$clst))
jmeta.ordered$clst <- factor(jmeta.ordered$clst, levels = clsts.keep)

jrow <- jmat.filt[jrow.i, ]
jrow.raw <- jmat.raw.filt[jrow.i, ]

fit.input <- data.frame(exprs = jrow.raw, exprs.imputed = exp(jrow), xvar = jmeta$xvar, clst = jmeta$clst, stringsAsFactors = TRUE)

fit.input.wide <- reshape2::melt(fit.input, id.vars = c("xvar", "clst"), measure.vars = c("exprs", "exprs.imputed"), variable.name = "jtype", value.name = "jexprs")

# ggplot(fit.input.wide, aes(x = xvar, color = clst, y = jexprs)) + geom_point() 
# 
# ggplot(fit.input, aes(x = xvar, y = exprs, color = clst)) + geom_point() + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# ggplot(fit.input, aes(x = xvar, y = log(exprs.imputed), color = clst)) + geom_point() + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# ggplot(fit.input, aes(x = clst, y = log(exprs.imputed), color = clst)) + geom_boxplot() + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# ggplot(fit.input, aes(x = clst, y = exprs, fill = clst)) + geom_boxplot() + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

m.var <- ggplot(fit.input, aes(x = clst, y = xvar, fill = clst)) + geom_boxplot() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +  
  ylab("MNaseActivity") + xlab("") + 
  ggtitle(jmark)


# Load previous poisson fits to compare xvar with Poisson estimates ---------------------------------------------



jfits.pois <- jfits.lst.bymark[[jmark]] %>%
  bind_rows()
jfits.pois$rname <- names(jfits.lst.bymark[[jmark]])

jfits.pois.long <- reshape2::melt(jfits.pois, id.vars = c("rname"), measure.vars = c("Clustereryth", "Clustergranu", "Clusterlymph"), variable.name = "clst", value.name = "logfc")

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

# ggplot(jfits.pois.long %>% filter(abs(logfc) < 5), aes(x = logfc, fill = clst)) + geom_density(alpha = 0.25) + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   geom_vline(xintercept = 0) + scale_fill_manual(values = cbPalette)
# 

m.fits <- ggplot(jfits.pois.long %>% filter(abs(logfc) < 5), aes(x = clst, y = logfc)) + geom_boxplot() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_fill_manual(values = cbPalette) + geom_hline(yintercept = 0, linetype = "dotted") + ggtitle(jmark)


# m.fits <- ggplot(jfits.pois.long %>% filter(abs(logfc) < 5), aes(x = logfc, fill = clst)) + geom_histogram(alpha = 0.33, position = "dodge") + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   scale_fill_manual(values = cbPalette) + ggtitle(jmark)

# print(m.fits)

JFuncs::multiplot(m.var, m.fits, cols = 2)




# Check mouse BM ----------------------------------------------------------


inf.bm <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/poisson_fits.OtherNorm.DiffDists/fit_poisson_model_on_TSS.MouseBM.NormMeth_ncuts.inbins.bsize_10000.RData"
inf.bm <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/poisson_fits.OtherNorm.DiffDists/fit_poisson_model_on_TSS.ZebrafishWKM.NormMeth_ByHetero.dist_10000.RData"


inf.bm <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/poisson_fits.OtherNorm.DiffDists/fit_poisson_model_on_TSS.MouseBM.NormMeth_ByTotalFromBins.bsize_10000.RData"
inf.bm <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/poisson_fits.OtherNorm.DiffDists/fit_poisson_model_on_TSS.ZebrafishWKM.NormMeth_ByTotalFromBins.dist_10000.RData"
inf.bm <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/poisson_fits.OtherNorm.DiffDists/fit_poisson_model_on_TSS.MouseBM.NormMeth_ncuts.alltss.bsize_10000.RData"
load(inf.bm, v=T)

jfits.pois.bm.lst <- lapply(jfits.lst.bymark, function(jdat){
  jnames <- names(jdat)
  jdat <- jdat %>% bind_rows()
  jdat$rname <- jnames
  return(jdat)
})

if (grepl("Zebrafish", inf.bm)){
  mvars <- c("Clustereryth", "Clustergranu", "Clusterlymph")
} else if (grepl("MouseBM", inf.bm)){
  mvars <- c("ClusterBcells", "ClusterErythroblasts", "ClusterGranulocytes")
}

jfits.pois.long.bm.lst <- lapply(jfits.pois.bm.lst, function(jdat){
  jfits.pois.long <- reshape2::melt(jdat, id.vars = c("rname"), measure.vars = mvars, variable.name = "clst", value.name = "logfc")
})

m.lst <- lapply(jmarks, function(jmark){
  jfits.pois.long.bm <- jfits.pois.long.bm.lst[[jmark]]
  m <- ggplot(jfits.pois.long.bm %>% filter(abs(logfc) < 5), aes(x = clst, y = logfc)) + geom_boxplot()  + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
    geom_hline(yintercept = 0) + 
    ggtitle(jmark)
})

m.lst
