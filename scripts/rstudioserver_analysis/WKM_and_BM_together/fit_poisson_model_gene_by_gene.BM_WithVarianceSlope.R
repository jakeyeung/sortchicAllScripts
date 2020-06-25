# Jake Yeung
# Date of Creation: 2020-06-25
# File: ~/projects/scchic/scripts/rstudioserver_analysis/WKM_and_BM_together/fit_poisson_model_gene_by_gene.BM_WithVarianceSlope.R
# 



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

# compare with lm? 
inf.lm <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/var_slope_estimates/MouseBM_log_lm_fits.2020-06-25.RData"
load(inf.lm, v=T)

jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
clsts.keep <- c("HSPCs", "lymph", "granu", "eryth")

# Load pois fits ----------------------------------------------------------


inf.poisfits <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/poisson_fits.OtherNorm.DiffDists/fit_poisson_model_on_TSS.MouseBM.NormMeth_ByHetero.bsize_10000.RData"
load(inf.poisfits, v=T)

# Load set up  ------------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

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
count.mat.lst <- lapply(out.objs, function(x) x$count.mat)


# Get ncuts total (offset) ------------------------------------------------

ncuts.dat <- lapply(count.mat.lst, function(jmat){
  data.frame(cell = colnames(jmat), ncuts.intss = colSums(jmat), stringsAsFactors = FALSE)
})


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


# Set up  -----------------------------------------------------------------


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
  left_join(dat.umaps.lst[[jmark]], subset(dat.annot.lst.BM[[jmark]], select = c(cell, cluster, cluster.orig)))
})



dat.imputed.lst <- lapply(tm.result.lst, function(tm.result){
  dat.imput <- t(log(tm.result$topics %*% tm.result$terms))  # natural log links with Poisson regression 
})



# Set up meta -------------------------------------------------------------


# fit genome-wide 
metadat.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  dat.impute <- dat.imputed.lst[[jmark]]
  dat.tmp <- CalculateVarAll(dat.impute, jchromos = jchromos) %>%
    left_join(., dat.umap.merged[[jmark]]) %>%
    filter(cluster %in% clsts.keep) %>%
    left_join(., ncuts.dat[[jmark]])
})

# Fit genes  --------------------------------------------------------------


# pax5 as a test?
# jmark <- "H3K27me3"
jmark <- "H3K4me3"
jmark <- "H3K27me3"


jmark <- "H3K4me3"

# jrow.i <- grep("F3", rownames(count.mat.lst[[jmark]]), value = TRUE)
jrow.i <- grep("Cebpd", rownames(count.mat.lst[[jmark]]), value = TRUE)
jrow.i <- grep("Hoxa4", rownames(count.mat.lst[[jmark]]), value = TRUE)

jrow.i <- grep("Hlf", rownames(count.mat.lst[[jmark]]), value = TRUE)


jmat <- dat.imputed.lst[[jmark]]
jmat.raw <- count.mat.lst[[jmark]]

cnames.keep <- metadat.lst[[jmark]]$cell
jmat.filt <- jmat[, cnames.keep]
jmat.raw.filt <- jmat.raw[, cnames.keep]

jmeta <- metadat.lst[[jmark]] %>%
  dplyr::rename(xvar = cell.var.within.sum.norm, clst = cluster) %>%
  dplyr::select(cell, xvar, clst, ncuts.intss) %>%
  group_by(clst) %>%
  mutate(xvar.orig = xvar,
         xvar = max(xvar) - xvar)

jmeta.ordered <- jmeta[match(colnames(jmat.filt), jmeta$cell), ]
assertthat::assert_that(identical(jmeta.ordered$cell, colnames(jmat.filt)))
assertthat::assert_that(identical(jmeta.ordered$cell, colnames(jmat.raw.filt)))
print(unique(jmeta.ordered$clst))
jmeta.ordered$clst <- factor(jmeta.ordered$clst, levels = clsts.keep)

jrow <- jmat.filt[jrow.i, ]
jrow.raw <- jmat.raw.filt[jrow.i, ]

fit.input <- data.frame(cell = jmeta.ordered$cell, ncuts = jrow.raw, exprs.imputed = exp(jrow), xvar = jmeta.ordered$xvar, clst = jmeta.ordered$clst, totalcuts = jmeta.ordered$ncuts.intss, stringsAsFactors = TRUE)



fit.input.wide <- reshape2::melt(fit.input, id.vars = c("xvar", "clst"), measure.vars = c("ncuts", "exprs.imputed"), variable.name = "jtype", value.name = "jexprs")

m.var <- ggplot(fit.input, aes(x = clst, y = xvar, fill = clst)) + geom_boxplot() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +  
  ylab("MNaseActivity") + xlab("") + 
  ggtitle(jmark)

print(m.var)

ggplot(fit.input, aes(x = xvar, y = exprs.imputed, color = clst)) + geom_point()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(fit.input, aes(x = xvar, y = ncuts / totalcuts, color = clst)) + geom_point()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_y_log10() 

# Aadd total counts -------------------------------------------------------

# Fit poiisson model and compare with linear fits -----------------------

head(fit.input)

jfit.pois.raw <- glm(formula = ncuts ~ 1 + clst + xvar:clst + offset(log(totalcuts)), data = fit.input)

jparams <- subset(dat.params.all.lst[[jmark]], rname == jrow.i) %>%
  rowwise() %>%
  mutate(clst = strsplit(params.mean, "clst")[[1]][[2]]) %>%
  mutate(clst = ifelse(clst == "HSC", "HSPCs", clst)) %>%
  ungroup() %>%
  mutate(clst = factor(clst, levels = clsts.keep))

# add the slope estimates directly into the offset and rereun 

head(fit.input)

fit.input.wlm <- left_join(fit.input, jparams) %>%
  rowwise() %>%
  mutate(joffset = log(totalcuts) + xvar * jslope)

jfit.pois.raw.clstonly.orig <- glm(formula = ncuts ~ 1 + clst + offset(log(totalcuts)), data = fit.input.wlm)
jfit.pois.raw.clstonly <- glm(formula = ncuts ~ 1 + clst + offset(joffset), data = fit.input.wlm)





