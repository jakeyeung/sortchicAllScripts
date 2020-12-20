# Jake Yeung
# Date of Creation: 2020-10-31
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged/1-LDA_downstream.WithRatios.again.K9me3_clean_with_celltypes.ForPeter.R
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

library(ggrepel)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)


SplitGetLast <- function(x, jsplit = "-"){
  # get last element after splitting
  xsplit <- strsplit(x, jsplit)[[1]]
  xnew <- xsplit[[length(xsplit)]]
  return(xnew)
}



DoFitsDownstream <- function(jfit, jfit.null, jgrep = "^experi|Intercept"){
  jcompare <- anova(jfit.null, jfit)
  jfit.effects <- data.frame(param = rownames(summary(jfit)$coefficients), summary(jfit)$coefficients, stringsAsFactors = FALSE) %>%
    filter(!grepl(jgrep, param))
  jfit.int <- data.frame(param = rownames(summary(jfit)$coefficients), summary(jfit)$coefficients, stringsAsFactors = FALSE) %>%
    filter(param == "(Intercept)") %>%
    rowwise()
  jint <- jfit.int$Estimate
  jint.se <- jfit.int$Std..Error
  
  jfit.int <- data.frame(param = rownames(summary(jfit)$coefficients), summary(jfit)$coefficients, stringsAsFactors = FALSE) %>%
    filter(param == "(Intercept)") %>%
    rowwise() %>%
    mutate(est = jint,
           est.se = jint.se)
  
  jfit.dat <- jfit.effects %>%
    rowwise() %>%
    mutate(est = Estimate + jint,
           est.se = sqrt(Std..Error ^ 2 + jint.se ^ 2))
  
  jfit.merge <- rbind(jfit.dat, jfit.int)
  return(list(jfit = jfit, jfit.null = jfit.null, jcompare = jcompare, jfit.merge = jfit.merge))
}

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

# outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/spikeins_mouse.BMround2_umaps_and_ratios"
indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2"
jdate <- "2020-11-01"
jmark <- "H3K27me3"
inf <- file.path(indir, paste0("spikeins_mouse.BMround1and2_umaps_and_ratios.colfix.celltyping.", jdate, ".WithRelLevels.mark_", jmark, ".cell_cluster_tables.txt"))
dat <- fread(inf)

dats.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  inf <- file.path(indir, paste0("spikeins_mouse.BMround1and2_umaps_and_ratios.colfix.celltyping.", jdate, ".WithRelLevels.mark_", jmark, ".cell_cluster_tables.txt"))
  dat <- fread(inf)
  if (jmark == "H3K9me3"){
    dat <- dat %>%
      rowwise() %>%
      mutate(plate = ClipLast(cell, jsep = "_"),
             experi = ClipLast(cell, jsep = "-"),
             batch = ClipLast(cell, jsep = "-"))
  }
  dat <- dat %>%
    rowwise() %>%
    mutate(batch = ClipLast(plate, jsep = "-"))
  return(dat)
})

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", '#abdf76', "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

m.lst <- lapply(jmarks, function(jmark){
  m <- ggplot(dats.lst[[jmark]], aes(x = umap1, y = umap2, color = cluster)) + 
    geom_point() + 
    scale_color_manual(values = cbPalette, na.value = "grey95") + 
    facet_wrap(~batch) + 
    ggtitle(jmark) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  return(m)
})

print(m.lst)
