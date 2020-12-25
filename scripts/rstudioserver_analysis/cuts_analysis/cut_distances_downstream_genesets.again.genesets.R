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

library(GGally)

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
# infmain <-  file.path(hubprefix, "jyeung/data/scChiC/cut_distances_all_clusters2.rerun.withraw")

gset.dummy <- "geneset_HSCs"
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

# "cut_distances_all_clusters2.rerun.bedfile_Topics_winsize_50000.topn_5000.withraw.fourCtypes/"

jwinsize <- 50000



# jtopn <- 5000
jtopn <- "topn_2000"
# jtopn <- "heterochromatin"


for (jmark in jmarks){
  
  infmain <- file.path(hubprefix, paste0("jyeung/data/scChiC/cut_distances_all_clusters2.rerun.bedfile_GeneSets_winsize_", jwinsize, ".withraw.fourCtypes"))
  assertthat::assert_that(dir.exists(infmain))
  outdir <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/pdfs_all/cut_distances.GeneSets.winsize_", jwinsize, ".withraw.fourCtypes"))
  dir.create(outdir)
  
  # ctypes <- c("Bcells_topic16", "Eryth-Gf1-_topic17", "Eryth-Slc7a6-_topic1", "Eryth-Sox6-_topic6", "HSCs-Tead1-_topic9", "InnateLymph_topic27", "Neutrophils_topic22")
  
  indirs.base <- list.files(infmain, pattern = paste0(gset.dummy, "_", jmark))
  ctypes <- sapply(indirs.base, function(indir) strsplit(indir, "\\.")[[1]][[2]], USE.NAMES = FALSE)
  names(ctypes) <- ctypes
  
  
  
  for (ctype in ctypes){
    if (ctype == ""){
      next
    }
    
    # jdir <- file.path(infmain, paste0(gset, "_", jmark, "-BM_AllMerged.", ctype, ".sorted"))
    # assertthat::assert_that(dir.exists(jdir))
    
    outpdf <- file.path(outdir, paste0("distances_summary_BM.", jmark, ".", ctype, ".", Sys.Date(), ".pdf"))
    if (file.exists(outpdf)){
      print(paste("Skipping...", outpdf))
      next
    }
    pdf(outpdf, useDingbats = FALSE)
    
    # pdf(file.path(outdir, paste0("distances_summary_BM.", jmark, ".", ctype, ".", gset, ".pdf")), useDingbats = FALSE)
    print(jmark)
    
    # inf.counts.dirs <- list.dirs(path = )
    inf.counts.dirs <- list.files(infmain, pattern = paste0(jmark, "-BM_AllMerged.", ctype, ".sorted.cleaned"), full.names = TRUE)
    names(inf.counts.dirs) <- sapply(inf.counts.dirs, function(x) strsplit(basename(x), split = "_")[[1]][[2]])
    
    jout <- lapply(inf.counts.dirs, function(jdir){
      basedir <- basename(jdir)
      ctype <- strsplit(basedir, split = "\\.")[[1]][[2]]
      if (ctype == ""){
        ctype <- "NA"
      }
      
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
      return(jcounts.filt)
    })
    
    jtitle <- paste(jmark, ctype)
    jout.csums <- lapply(jout, rowMeans)
    
    jmat <- as.data.frame(do.call(cbind, jout.csums))
    rownames(jmat) <- seq(nrow(jmat))
    
    jmat.filt <- jmat[20:nrow(jmat), ]
    
    x <- t(as.matrix(BinarizeMatrix(as.matrix(jout$Bcell))))
    gplots::heatmap.2(log(t(jout$Bcell) + 1), Rowv = NA, Colv = NA, trace = "none")
    gplots::heatmap.2(as.matrix(BinarizeMatrix(as.matrix(jout$Bcell))), Rowv = NA, Colv = NA, trace = "none")
    gplots::heatmap.2(x), Rowv = NA, Colv = NA, trace = "none")
    plot(density(as.matrix(BinarizeMatrix(as.matrix(jout$Bcell)))))
    
    
    counts.long <- as.matrix(jmat.filt) %>%
      melt()
    
    cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
    
    m1 <- ggplot(counts.long, aes(x = Var1, y = value, color = Var2)) + 
      geom_point(alpha = 0.25, size = 0.1) + geom_line(alpha = 0.25) + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
      scale_color_manual(values = cbPalette) + 
      ggtitle(jtitle)
    
    m2 <- ggplot(counts.long %>% filter(Var1 < 500), aes(x = Var1, y = value, color = Var2)) + 
      geom_point(alpha = 0.25, size = 0.1) + geom_line(alpha = 0.25) +
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
      scale_color_manual(values = cbPalette) + 
      xlab("Distance from Cut") + 
      ylab("Normalized counts") + 
      ggtitle(jtitle)
    
    m3 <- ggplot(counts.long %>% filter(Var1 < 1000), aes(x = Var1, y = value, color = Var2)) + 
      geom_point(alpha = 0.25, size = 0.1) + geom_line(alpha = 0.25) +
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
      scale_color_manual(values = cbPalette) + 
      xlab("Distance from Cut") + 
      ylab("Normalized counts") + 
      ggtitle(jtitle)
    
    m4 <- GGally::ggpairs(as.data.frame(jmat.filt)) + ggtitle(jtitle, paste("Distance:", nrow(jmat.filt)))
    m5 <- GGally::ggpairs(as.data.frame(jmat.filt[1:250, ])) + ggtitle(jtitle, paste("Distance:", nrow(jmat.filt[1:250, ])))
    
    print(m1)
    print(m2)
    print(m3)
    print(m4)
    print(m5)
    
    
    dev.off()
  }
  
  
  
}








