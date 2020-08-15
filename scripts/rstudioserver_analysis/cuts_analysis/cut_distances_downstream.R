# Jake Yeung
# Date of Creation: 2020-07-29
# File: ~/projects/scchic/scripts/rstudioserver_analysis/cuts_analysis/cut_distances_downstream.R
# Downstream

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)

library(zoo)

library(hash)
library(igraph)
library(umap)

library(topicmodels)


# Functions ---------------------------------------------------------------



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


# Load data ---------------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")


bm.rename <- as.list(hash(c("Bcells", "Eryth", "HSCs", "Neutrophils"), c("lymph", "eryth", "HSPCs", "granu")))

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

ggplot(dat.annot.lst.BM$H3K27me3, aes(x = umap1, y = umap2, color = cluster)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Load outputs ------------------------------------------------------------

jmark <- "H3K27me3"
infmain <-  file.path(hubprefix, "jyeung/data/scChiC/cut_distances_all_clusters2.rerun.withraw")
assertthat::assert_that(dir.exists(infmain))

pdf("/home/jyeung/data/from_rstudioserver/distances/distances_summary_BM.pdf", useDingbats = FALSE)

jout2 <- lapply(jmarks, function(jmark){
  print(jmark)
  
  # inf.counts.dirs <- list.dirs(path = )
  inf.counts.dirs <- list.files(infmain, pattern = paste0(jmark), full.names = TRUE)
  
  jout <- lapply(inf.counts.dirs, function(jdir){
    basedir <- basename(jdir)
    jctype <- strsplit(basedir, split = "\\.")[[1]][[2]]
    if (jctype == ""){
      jctype <- "NA"
    }
    
    jtitle <- paste(jmark, jctype)
    inf.counts <- file.path(jdir, "strand_unspecific_counts_raw.csv")
    print(jdir)
    assertthat::assert_that(file.exists(inf.counts))
    
    
    jcounts <- fread(inf.counts)
    jcounts$V1 <- NULL
    jcounts <- as.matrix(jcounts)
    
    
    jcounts.filt <- apply(as.matrix(jcounts), MARGIN = 2, FUN = function(jcol){
      zoo::rollapply(jcol, width = 35, FUN = sum)
    })
    
    
    jcounts.filt <- jcounts.filt[, which(colSums(jcounts.filt) > 0)]
    
    jcounts.filt <- sweep(jcounts.filt, MARGIN = 2, STATS = colSums(jcounts.filt), FUN = "/")
    
    # remove first 20 bases?
    
    jcounts.filt2 <- DescTools::Winsorize(jcounts.filt, probs = c(0, 0.95))
    
    plot(rowSums(jcounts.filt)[1:nrow(jcounts.filt)], main = jtitle)
    plot(rowSums(jcounts.filt2)[1:nrow(jcounts.filt2)], main = jtitle)
    
    jcounts.long <- data.frame(distance = seq(nrow(jcounts.filt)), jcounts.filt, stringsAsFactors = FALSE) %>%
      melt(id.vars = "distance", variable.name = "cell", value.name = "counts.norm") %>%
      ungroup() %>%
      mutate(cell = gsub("\\.", "-", cell),
             counts.norm.wins = DescTools::Winsorize(counts.norm, probs = c(0, 0.95)))
    
    jcounts.sum <- jcounts.long %>%
      group_by(distance) %>%
      summarise(counts.norm.mean = mean(counts.norm),
                counts.norm.se = sd(counts.norm),
                counts.norm.wins.mean = mean(counts.norm.wins),
                counts.norm.wins.se = sd(counts.norm.wins))
    
    m <- ggplot(jcounts.sum %>% filter(distance > 20), aes(x = distance, y = counts.norm.wins.mean, 
                                                           ymin = counts.norm.wins.mean - 1 * counts.norm.wins.se, ymax = counts.norm.wins.mean + 1 * counts.norm.wins.se)) + 
      geom_line() + 
      geom_ribbon(alpha = 0.25, fill = 'green') + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      xlab("Distance to nearest cut") + 
      ylab("Normalized counts (winsorized)") + 
      ggtitle(jtitle)
    print(m)
    
    
    # order by variance
    jcells <- colnames(jcounts.filt)
    dat.vars.filt <- subset(dat.vars.lst[[jmark]], cell %in% jcells) %>%
      arrange(cell.var.within.sum.norm)
    
    dat.ncuts <- data.frame(ncuts = colSums(count.mat[[jmark]]), cell = colnames(count.mat[[jmark]]), stringsAsFactors = FALSE) %>%
      arrange(ncuts) %>%
      filter(cell %in% jcells)
    
    # jcells.ordered <- dat.vars.filt$cell
    # jcells.ordered <- dat.ncuts$cell
    
    dat.var.sub <- as.data.frame(subset(dat.var.lst[[jmark]], cell %in% jcells))
    rownames(dat.var.sub) <- dat.var.sub$cell
    dat.var.sub <- dat.var.sub[jcells, ]
    
    colvec <- viridis::viridis(100, direction = 1)[as.numeric(cut(dat.var.sub$cell.var.within.sum.norm, breaks = 100))]
    colvec.totalcuts <- viridis::viridis(100, direction = 1)[as.numeric(cut(log(dat.var.sub$totalcounts), breaks = 100))]
    heatmap3::heatmap3(t(jcounts.filt2), Rowv = NULL, Colv = NA, scale = "none", RowSideColors = colvec, RowSideLabs = "IntrachromVar", margin = c(5, 20), col = viridis::viridis(50), method = 'ward.D', main = jtitle)
    heatmap3::heatmap3(t(jcounts.filt2), Rowv = NULL, Colv = NA, scale = "none", RowSideColors = colvec.totalcuts, RowSideLabs = "TotalCuts", margin = c(5, 20), col = viridis::viridis(50), method = 'ward.D', main = jtitle)
    
  })
  
})
dev.off()

