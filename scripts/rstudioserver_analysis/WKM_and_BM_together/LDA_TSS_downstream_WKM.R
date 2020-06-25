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

library(hash)
library(igraph)
library(umap)

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


jgene <- "tal1"

jgene <- "pax5"

jgene <- "hbb"

exprs.lst <- lapply(jmarks, function(jmark){
  jmat <- dat.imputed.lst[[jmark]]
  rnames.keep <- grepl(jgene, rownames(jmat))
  if (length(rnames.keep) == 1){
    exprs.vec <- jmat[rnames.keep, ]
  } else {
    exprs.vec <- colMeans(jmat[rnames.keep, ])
  }
  exprs.dat <- data.frame(exprs = exprs.vec, cell = names(exprs.vec), stringsAsFactors = FALSE) 
  left_join(dat.umaps.lst[[jmark]], exprs.dat)
})

m.exprs.lst <- lapply(jmarks, function(jmark){
  dat.umap <- exprs.lst[[jmark]]
  m <- ggplot(dat.umap, aes(x = umap1, y = umap2, color = exprs)) + geom_point(alpha = 0.3) + theme_bw() + 
    scale_color_viridis_c() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    ggtitle(paste(jmark, jgene))
})

JFuncs::multiplot(m.exprs.lst$H3K4me1, m.exprs.lst$H3K4me3, m.exprs.lst$H3K27me3, cols = 3)


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

dat.umap.merged <- lapply(jmarks, function(jmark){
  left_join(dat.umaps.lst[[jmark]], subset(dat.annot.lst.WKM[[jmark]], select = c(cell, cluster, cluster.orig)))
})


# Overlay old labels  -----------------------------------------------------

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
m.overlay <- lapply(jmarks, function(jmark){
  ggplot(dat.umap.merged[[jmark]], aes(x = umap1, y = umap2, color = cluster)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    scale_color_manual(values = cbPalette, na.value = "grey85") + ggtitle(jmark)
})
JFuncs::multiplot(m.overlay$H3K4me1, m.overlay$H3K4me3, m.overlay$H3K27me3, cols = 3)

m.overlay.orig <- lapply(jmarks, function(jmark){
  ggplot(dat.umap.merged[[jmark]], aes(x = umap1, y = umap2, color = cluster.orig)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    scale_color_manual(values = cbPalette, na.value = "grey85") + ggtitle(jmark)
})
JFuncs::multiplot(m.overlay.orig$H3K4me1, m.overlay.orig$H3K4me3, m.overlay.orig$H3K27me3, cols = 3)



# Load RNAseq genes  ------------------------------------------------------

inf.scrnaseq <- "/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/pseudobulk_analysis_all/pseudobulk_analysis.get_DE_genes_from_pbulk_scrnaseq/de_genes_sorted_and_giladi.WithHouseKeepAndNotExpressed.FixExprsOther.exprsmax_5.logfcmin_2.logfcmax_1.RData"
load(inf.scrnaseq, v=T)


# Take top bins in H3K4me3 and see in H3K27me3 ----------------------------

topics.ordered.lst <- lapply(tm.result.lst, function(tm.result){
  OrderTopicsByEntropy(tm.result = tm.result, jquantile = 0.99)
})

# plot a topic? 
jmark <- "H3K4me1"
jtop <- "topic5"
jtop <- "topic1"
jtop <- "topic16"
jtop <- "topic25"
jtop <- "topic19"
jtop <- "topic3"

jtop <- "topic5"
jmark <- "H3K4me3"
jtop <- "topic21"
jtop <- "topic10"
jtop <- "topic15"

jbins <- sort(tm.result.lst[[jmark]]$terms[jtop, ], decreasing = TRUE)[1:150]
print(head(jbins))

# plot topic weights


# plot top bins

# integrate with RNAseq dataset



# Plot gene sets  ---------------------------------------------------------

rnames.all <- Reduce(intersect, lapply(dat.imputed.lst, rownames))
jgenes.all <- sapply(rnames.all, function(x) strsplit(x, ";")[[1]][[2]])
jens.all <- sapply(jgenes.all, AssignHash, g2e.zf)

e2g.zf <- hash::invert(g2e.zf)

# filter by gene sets? 
rnames.bygset <- lapply(de.ens.zf.stringent, function(ens.vec){
  jgenes <- sapply(ens.vec, AssignHash, e2g.zf)
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

gsetnames <- names(gset.dat); names(gsetnames) <- gsetnames
m.gsets <- lapply(gsetnames, function(gsetname){
  m.lst <- lapply(jmarks, function(jmark){
    m.tmp <- ggplot(gset.dat[[gsetname]][[jmark]], aes(x = umap1, y = umap2, color = exprs)) + 
      geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
      ggtitle(paste(jmark, gsetname))
    if (jmark == "H3K27me3"){
      m.tmp <- m.tmp + scale_color_viridis_c(direction = 1) 
    } else {
      m.tmp <- m.tmp + scale_color_viridis_c(direction = 1) 
    }
  })
})

for (gsetname in gsetnames){
  JFuncs::multiplot(m.gsets[[gsetname]]$H3K4me1, m.gsets[[gsetname]]$H3K4me3, m.gsets[[gsetname]]$H3K27me3, cols = 3)
}

jgene <- "tal1"
jgene <- "pax5"
jgene <- "hbb"

exprs.lst <- lapply(jmarks, function(jmark){
  jmat <- dat.imputed.lst[[jmark]]
  rnames.keep <- grepl(jgene, rownames(jmat))
  if (length(rnames.keep) == 1){
    exprs.vec <- jmat[rnames.keep, ]
  } else {
    exprs.vec <- colMeans(jmat[rnames.keep, ])
  }
  exprs.dat <- data.frame(exprs = exprs.vec, cell = names(exprs.vec), stringsAsFactors = FALSE) 
  left_join(dat.umaps.lst[[jmark]], exprs.dat)
})



# Compare exprs with variance?  -------------------------------------------

jchromos <- paste("chr", c(seq(25)), sep = "")

print(names(gset.dat))
jgset <- "lymphocytes"
jgset <- "lymphocytes"
jgset <- "HighExprs"


jgset <- "HSPCs"
dat.vars.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  dat.impute <- dat.imputed.lst[[jmark]]
  dat.tmp <- CalculateVarAll(dat.impute, jchromos = jchromos) %>%
    left_join(., dat.umap.merged[[jmark]]) %>%
    left_join(., subset(gset.dat[[jgset]][[jmark]] %>% dplyr::select(c(cell, exprs))))
})

# ggplot(dat.vars.lst$H3K27me3, aes(x = umap1, y = umap2, color = exprs)) + 
#   scale_color_viridis_c(direction = 1) + 
#   geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.vars.lst$H3K27me3, aes(x = exprs, y = cell.var.within.sum.norm, color = cluster)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette) + 
  ggtitle(jgset)
  # facet_wrap(~cluster)

ggplot(dat.vars.lst$H3K4me1, aes(x = exprs, y = cell.var.within.sum.norm, color = cluster.orig)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette) + 
  ggtitle(jgset)

ggplot(dat.vars.lst$H3K4me3, aes(x = exprs, y = cell.var.within.sum.norm, color = cluster.orig)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette) + 
  ggtitle(jgset)


ggplot(dat.vars.lst$H3K27me3, aes(y = exprs, x = cell.var.within.sum.norm, color = cluster)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette) + 
  ggtitle(jgset)
  # facet_wrap(~cluster)

ggplot(dat.vars.lst$H3K4me1, aes(y = exprs, x = cell.var.within.sum.norm, color = cluster.orig)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette) + 
  ggtitle(jgset)

ggplot(dat.vars.lst$H3K4me3, aes(y = exprs, x = cell.var.within.sum.norm, color = cluster.orig)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette) + 
  ggtitle(jgset)

# Try to correct by fitting linear model  ---------------------------------

jsub.lst <- lapply(dat.vars.lst, function(dat.vars){
  subset(dat.vars, !is.na(cluster.orig))
})

jmark <- "H3K27me3"
y <- jsub.lst[[jmark]]$exprs
xvar <- jsub.lst[[jmark]]$cell.var.within.sum.norm
xclstr <- as.factor(jsub.lst[[jmark]]$cluster)

fit.input <- data.frame(exprs = y, xvar = xvar, clst = xclstr, stringsAsFactors = FALSE)

jfit <- lm(exprs ~ xvar + clst, fit.input)

fit.pred <- predict(jfit)

fit.input$pred <- fit.pred

ggplot(fit.input, aes(x = xvar, y = exprs, color = clst)) + geom_point()  +  geom_line(mapping = aes(x = xvar, y = pred, color = clst)) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~clst) + ggtitle(jgset)

ggplot(fit.input, aes(x = xvar, y = exprs, color = clst)) + geom_point()  +  
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~clst) + ggtitle(jgset)

ggplot(fit.input, aes(x = max(xvar) - xvar, y = exprs, color = clst)) + geom_point()  +  
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~clst) + ggtitle(jgset)

# Can I predict when slopes go up or down?  -------------------------------

clsts.keep <- c("HSC", "lymph", "monocyte", "eryth")
# fit genome-wide 
metadat.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  dat.impute <- dat.imputed.lst[[jmark]]
  dat.tmp <- CalculateVarAll(dat.impute, jchromos = jchromos) %>%
    left_join(., dat.umap.merged[[jmark]]) %>%
    filter(cluster %in% clsts.keep)
})

# for H3K4me1 first as sanity check 
jmark <- "H3K4me1"
jmat <- dat.imputed.lst[[jmark]]
cnames.keep <- metadat.lst[[jmark]]$cell
jmat.filt <- jmat[, cnames.keep]

jmeta <- metadat.lst[[jmark]] %>%
  dplyr::rename(xvar = cell.var.within.sum.norm, clst = cluster) %>%
  dplyr::select(cell, xvar, clst) %>%
  # group_by(clst) %>%
  mutate(xvar.orig = xvar,
         xvar = max(xvar) - xvar)
jmeta.ordered <- jmeta[match(colnames(jmat.filt), jmeta$cell), ]
assertthat::assert_that(identical(jmeta.ordered$cell, colnames(jmat.filt)))
print(unique(jmeta.ordered$clst))
jmeta.ordered$clst <- factor(jmeta.ordered$clst, levels = clsts.keep)

system.time(
  jfits.rows <- apply(jmat.filt, 1, function(jrow, jmeta){
    FitClstVar(jrow, jmeta)
  }, jmeta = jmeta.ordered)
)


# Plot params ------------------------------------------------------------

# jfits.rows$`chr1:-15081-34919;rpl24;2`
jrowname <- names(jfits.rows)[[2]]
jrowname <- grep("tal1", names(jfits.rows), value = TRUE)
jrowname <- grep("pax5", names(jfits.rows), value = TRUE)
jrowname <- grep("sox6", names(jfits.rows), value = TRUE)
# plot inputs and fits
y <- jmat.filt[jrowname, ]
xvar <- jmeta.ordered$xvar
xclstr <- as.factor(jsub.lst[[jmark]]$cluster)
jfit <- jfits.rows[[jrowname]]
# jfit <- FitClstVar(y, jmeta)
jpred <- predict(jfit)
fit.input <- data.frame(exprs = y, xvar = jmeta.ordered$xvar, xvar.orig = jmeta.ordered$xvar.orig, clst = jmeta.ordered$clst, pred = jpred, stringsAsFactors = FALSE)

ggplot(fit.input, aes(x = xvar, y = exprs, color = clst)) + geom_point()  +  geom_line(aes(x = xvar, y = pred, color = clst)) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~clst) + ggtitle(jrowname)

ggplot(fit.input, aes(x = xvar.orig, y = exprs, color = clst)) + geom_point()  +  
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~clst) + ggtitle(jrowname)

# collect slope and mean exprs for every gene
jint <- coefficients(jfit)[["(Intercept)"]]  # HSPC intercept
jmeans <- jint + coefficients(jfit)[which(startsWith(x = names(coefficients(jfit)), prefix = "clst"))]
jslopes <- coefficients(jfit)[which(startsWith(x = names(coefficients(jfit)), prefix = "xvar"))]
dat.params <- data.frame(params.mean = c("HSC", names(jmeans)), jmean = c(jint, jmeans), params.slope = names(jslopes), jslope = jslopes, stringsAsFactors = FALSE)



# Collect all parameters --------------------------------------------------


jnames <- names(jfits.rows); names(jnames) <- jnames
dat.params.all <- lapply(jnames, function(jname){
  jfit <- jfits.rows[[jname]]
  jint <- coefficients(jfit)[["(Intercept)"]]  # HSPC intercept
  jmeans <- jint + coefficients(jfit)[which(startsWith(x = names(coefficients(jfit)), prefix = "clst"))]
  jslopes <- coefficients(jfit)[which(startsWith(x = names(coefficients(jfit)), prefix = "xvar"))]
  dat.params <- data.frame(params.mean = c("clstHSC", names(jmeans)), jmean = c(jint, jmeans), params.slope = names(jslopes), jslope = jslopes, stringsAsFactors = FALSE)
  dat.params$rname <- jname
  return(dat.params)
}) %>%
  bind_rows()

# Plot when slopes turn positive or negative  -----------------------------

ggplot(dat.params.all, aes(x = jmean, y = jslope, color = params.mean)) + geom_point(alpha = 0.1)  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~params.mean)

ggplot(dat.params.all, aes(x = jmean, y = jslope, color = params.mean)) + geom_point(alpha = 0.1)  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Fit all 3 marks and save ------------------------------------------------


# for H3K4me1 first as sanity check 
system.time(
  jfits.rows.lst <- lapply(jmarks, function(jmark){
    print(jmark)
    jmat <- dat.imputed.lst[[jmark]]
    cnames.keep <- metadat.lst[[jmark]]$cell
    jmat.filt <- jmat[, cnames.keep]
    
    jmeta <- metadat.lst[[jmark]] %>%
      dplyr::rename(xvar = cell.var.within.sum.norm, clst = cluster) %>%
      dplyr::select(cell, xvar, clst) %>%
      # group_by(clst) %>%
      mutate(xvar.orig = xvar,
             xvar = max(xvar) - xvar)
    jmeta.ordered <- jmeta[match(colnames(jmat.filt), jmeta$cell), ]
    assertthat::assert_that(identical(jmeta.ordered$cell, colnames(jmat.filt)))
    print(unique(jmeta.ordered$clst))
    jmeta.ordered$clst <- factor(jmeta.ordered$clst, levels = clsts.keep)
    
    jfits.rows <- apply(jmat.filt, 1, function(jrow, jmeta){
      FitClstVar(jrow, jmeta)
    }, jmeta = jmeta.ordered)
  })
)



dat.params.all.lst <- lapply(jfits.rows.lst, function(jfits.rows){
  jnames <- names(jfits.rows); names(jnames) <- jnames
  dat.params.all <- lapply(jnames, function(jname){
    jfit <- jfits.rows[[jname]]
    jint <- coefficients(jfit)[["(Intercept)"]]  # HSPC intercept
    jmeans <- jint + coefficients(jfit)[which(startsWith(x = names(coefficients(jfit)), prefix = "clst"))]
    jslopes <- coefficients(jfit)[which(startsWith(x = names(coefficients(jfit)), prefix = "xvar"))]
    dat.params <- data.frame(params.mean = c("clstHSC", names(jmeans)), jmean = c(jint, jmeans), params.slope = names(jslopes), jslope = jslopes, stringsAsFactors = FALSE)
    dat.params$rname <- jname
    return(dat.params)
  }) %>%
    bind_rows()
  return(dat.params.all)
})

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/var_slope_estimates"
fname <- paste0("ZebrafishWKM_log_lm_fits.", Sys.Date(), ".RData")
outf <- file.path(outdir, fname)
save(dat.params.all.lst, jmat.filt, metadat.lst, file = outf)

# saveRDS(jfits.rows.lst, outf)








