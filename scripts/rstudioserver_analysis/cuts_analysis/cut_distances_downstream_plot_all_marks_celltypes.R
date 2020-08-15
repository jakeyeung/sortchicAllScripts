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



# Load UMAPs containing intrachromvar and totalcuts -----------------------



# Constants ---------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"


cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123


# set conditions ----------------------------------------------------------

mergesize <- "1000"
nbins <- "1000"
jcovar.cname <- "ncuts.var.log2.CenteredAndScaled"
jpenalty <- 1

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
jconds <- c("AllMerged")


# jcondsmarks <- levels(interaction(jconds, jmarks, sep = "_"))
jcondsmarks <- jmarks
names(jcondsmarks) <- jcondsmarks

jchromos <- paste("chr", seq(19), sep = "")

jexperi <- "AllMerged"
dat.var.lst <- lapply(jcondsmarks, function(jcondmark){
  print(jcondmark)
  
  print(jcondmark)
  # jexperi <- strsplit(jcondmark, "_")[[1]][[1]] 
  # jmark <- strsplit(jcondmark, "_")[[1]][[2]]
  jmark <- jcondmark
  
  # Load LDA output ---------------------------------------------------------
  
  
  inf.lda <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-02-11.var_filt.UnenrichedAndAllMerged.KeepBestPlates2/lda_outputs.BM_", jmark, "_varfilt_countmat.2020-02-11.", jexperi, ".K-30.binarize.FALSE/ldaOut.BM_", jmark, "_varfilt_countmat.2020-02-11.", jexperi, ".K-30.Robj"))
  assertthat::assert_that(file.exists(inf.lda))
  load(inf.lda, v=T)
  
  tm.result <- posterior(out.lda)
  colnames(tm.result$topics) <- paste("topic", colnames(tm.result$topics), sep = "")
  rownames(tm.result$terms) <- paste("topic", rownames(tm.result$terms), sep = "")
  
  dat.impute.log <- t(log2(tm.result$topics %*% tm.result$terms))
  
  dat.var <- CalculateVarAll(dat.impute.log, jchromos)
  
  # add total counts
  dat.counts <- data.frame(cell = colnames(count.mat), totalcounts = colSums(count.mat), stringsAsFactors = FALSE) 
  
  dat.var.merged <- left_join(dat.var, dat.counts, by = "cell")
  
  return(dat.var.merged)
})


# Load nucleosome positioning ---------------------------------------------


# Load outputs ------------------------------------------------------------

# jmark <- "H3K27me3"
infmain <-  file.path(hubprefix, "jyeung/data/scChiC/cut_distances_all_clusters2.rerun.withraw")
assertthat::assert_that(dir.exists(infmain))

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/cut_distances"


jout2 <- lapply(jmarks, function(jmark){
  
  pdf(file.path(outdir, paste0("distances_summary_BM.", jmark, ".pdf")), useDingbats = FALSE)
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
    dat.var.filt <- subset(dat.var.lst[[jmark]], cell %in% jcells) %>%
      arrange(cell.var.within.sum.norm)
    
    
    dat.var.sub <- as.data.frame(subset(dat.var.lst[[jmark]], cell %in% jcells))
    rownames(dat.var.sub) <- dat.var.sub$cell
    dat.var.sub <- dat.var.sub[jcells, ]
    
    colvec <- viridis::viridis(100, direction = 1)[as.numeric(cut(dat.var.sub$cell.var.within.sum.norm, breaks = 100))]
    colvec.totalcuts <- viridis::viridis(100, direction = 1)[as.numeric(cut(log(dat.var.sub$totalcounts), breaks = 100))]
    heatmap3::heatmap3(t(jcounts.filt2), Rowv = NULL, Colv = NA, scale = "none", RowSideColors = colvec, RowSideLabs = "IntrachromVar", margin = c(5, 20), col = viridis::viridis(50), method = 'ward.D', main = jtitle)
    heatmap3::heatmap3(t(jcounts.filt2), Rowv = NULL, Colv = NA, scale = "none", RowSideColors = colvec.totalcuts, RowSideLabs = "TotalCuts", margin = c(5, 20), col = viridis::viridis(50), method = 'ward.D', main = jtitle)
    
  })
  dev.off()
  
})


# 
# 
# # Load outputs ------------------------------------------------------------
# 
# 
# inmain <- file.path(hubprefix, paste0("jyeung/data/scChiC/cut_distances_all_clusters2.rerun.withraw"))
# 
# jmark <- "H3K27me3"
# 
# inf.counts <- file.path(hubprefix, paste0("jyeung/data/scChiC/cut_distances_all_clusters2.rerun.withraw/", jmark, "-BM_AllMerged.HSCs-Tead1-_topic9.sorted/strand_unspecific_counts_raw.csv"))
# 
# jcounts <- fread(inf.counts)
# jcounts$V1 <- NULL
# jcounts <- as.matrix(jcounts)
# 
# jcounts.filt <- apply(as.matrix(jcounts), MARGIN = 2, FUN = function(jcol){
#   zoo::rollapply(jcol, width = 35, FUN = sum)
# })
# 
# 
# jcounts.filt <- jcounts.filt[, which(colSums(jcounts.filt) > 0)]
# 
# jcounts.filt <- sweep(jcounts.filt, MARGIN = 2, STATS = colSums(jcounts.filt), FUN = "/")
# 
# plot(rowSums(jcounts.filt)[20:nrow(jcounts.filt)])
# 
# 
# 
# # plot(density(unlist(jcounts.filt)))
# # heatmap3::heatmap3(t(jcounts.filt), Rowv = NULL, Colv = NA, scale = "none")
# 
# # order by variance
# jcells <- colnames(jcounts.filt)
# dat.var.filt <- subset(dat.var.lst[[jmark]], cell %in% jcells) %>%
#   arrange(cell.var.within.sum.norm)
# 
# dat.ncuts <- data.frame(ncuts = colSums(count.mat[[jmark]]), cell = colnames(count.mat[[jmark]]), stringsAsFactors = FALSE) %>%
#   arrange(ncuts) %>%
#   filter(cell %in% jcells)
# 
# # jcells.ordered <- dat.var.filt$cell
# jcells.ordered <- dat.ncuts$cell
# 
# colvec <- viridis::viridis(100, direction = 1)[as.numeric(cut(dat.var.filt$cell.var.within.sum.norm, breaks = 100))]
# # jcounts.filt <- log2(jcounts.filt + 1)
# jctype <- "HSCs"
# 
# # heatmap3::heatmap3(t(jcounts.filt[, jcells.ordered]), Rowv = NA, Colv = NA, scale = "none", RowSideColors = colvec, RowSideLabs = "TotalCuts", main = paste(jmark, jctype), margins = c(5, 20), col = viridis::viridis(50, end = 1))
# # 
# # heatmap3::heatmap3(t(jcounts.filt[, jcells.ordered]), Rowv = NULL, Colv = NA, scale = "none", RowSideColors = colvec, RowSideLabs = "TotalCuts", main = paste(jmark, jctype), margin = c(5, 20), col = viridis::viridis(50, end = 1), method = 'ward.D')
# # 
# # jcounts.filt2 <- log2(t(jcounts.filt[, jcells.ordered]) * 1e3 + 1)
# 
# jcounts.filt2 <- DescTools::Winsorize(jcounts.filt2, probs = c(0, 0.95))
# plot(density(jcounts.filt2))
# 
# heatmap3::heatmap3(jcounts.filt2, Rowv = NULL, Colv = NA, scale = "none", RowSideColors = colvec, RowSideLabs = "TotalCuts", main = paste(jmark, jctype), margin = c(5, 20), col = viridis::viridis(50), method = 'ward.D')


