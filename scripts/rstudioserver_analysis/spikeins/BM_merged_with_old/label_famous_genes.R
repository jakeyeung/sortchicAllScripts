# Jake Yeung
# Date of Creation: 2020-12-09
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged_with_old/16-add_geneset_markers_to_UMAP.use_genes_from_topics.heatmap_with_genenames.MakeGoodGenes.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)

library(hash)
library(igraph)
library(umap)

library(scchicFuncs)
library(topicmodels)

library(heatmap3)

library(DescTools)

stypecols <- c("grey", "red", "blue")

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
jdist <- 10000
hubprefix <- "/home/jyeung/hub_oudenaarden"
niter <- 500
binskeep <- 0

# winsorize constants
jprobmin <- 0.04
jprobmax <- 0.96

# keepn <- 150
keepn <- 400

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/genesets_on_umap.genes_refmarks"
dir.create(outdir)
# outpdf <- file.path(outdir, paste0("geneset_on_umap.binskeep_", binskeep, ".niter_", niter, ".", Sys.Date(), ".from_LDA_topics.condensed.heatmap.famousgenes.keepn_", keepn, ".refmark_", refmark, ".", Sys.Date(), ".pdf"))
# outbase <- file.path(outdir, paste0("geneset_on_umap.binskeep_", binskeep, ".niter_", niter, ".", Sys.Date(), ".from_LDA_topics.metadata.condensed.heatmap"))

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks


# Add spikein annots  -----------------------------------------------------

indir.meta <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2.redo_after_varfilt_with_spikeins"
jdate <- "2020-11-18"
dat.metas.init <- lapply(jmarks, function(jmark){
  fname <- paste0("cell_cluster_table_with_spikeins.", jmark, ".", jdate, ".dupfilt.txt")
  inf.meta <- file.path(indir.meta, fname)
  dat.tmp <- fread(inf.meta)
  dat.tmp <- subset(dat.tmp, select = -c(umap1, umap2))
})

# fix stype for round1 and round2 
dat.round1.lst <- lapply(dat.metas.init, function(x) subset(x, batch != "Round2"))
dat.round2.lst <- lapply(dat.metas.init, function(x) subset(x, batch == "Round2"))

dat.round2.reannot.lst <- lapply(jmarks, function(jmark){
  jreannot <- dat.round2.lst[[jmark]] %>%
    rowwise() %>%
    mutate(experi = ClipLast(cell, jsep = "_"),
           plate = as.numeric(SplitGetLast(experi, jsplit = "-")),
           rowcoord = AddPlateCoordinates(cell)$rowcoord,
           colcoord = AddPlateCoordinates(cell)$colcoord,
           jrep = GetRepBM(experiname = experi), 
           batch = AnnotateSortFromLayoutBMall(plate, rowcoord, colcoord, jrep, jmark))
})
dat.round1.reannot.lst <- lapply(jmarks, function(jmark){
  jreannot <- dat.round1.lst[[jmark]] %>%
    rowwise() %>%
    mutate(experi = ClipLast(cell, jsep = "_"),
           plate = as.numeric(SplitGetLast(experi, jsplit = "-")),
           rowcoord = AddPlateCoordinates(cell)$rowcoord,
           colcoord = AddPlateCoordinates(cell)$colcoord,
           jrep = "rep1old")
})

dat.metas <- lapply(jmarks, function(jmark){
  jreannot1 <- dat.round1.reannot.lst[[jmark]]
  jreannot2 <- dat.round2.reannot.lst[[jmark]]
  jreannot <- rbind(jreannot1, jreannot2) %>%
    ungroup() %>%
    mutate(batch = gsub("LinNeg", "Linneg", batch),
           batch = gsub("LSK", "StemCell", batch))
  jreannot$stype <- jreannot$batch
  return(jreannot)
})

# set up colors
clstrs <- sort(unique(dat.metas$H3K4me1$cluster))  # use reference mark with all the celltypes
clstrs.col <- cbPalette[1:length(clstrs)]
clstrs.hash <- hash::hash(clstrs, clstrs.col)
head(unique(dat.round1.lst$H3K4me1$stype))
head(unique(dat.round2.lst$H3K4me1$stype))

# Load GLMPCA -------------------------------------------------------------

binskeep <- 0
iter <- 500

indir.glmpca <- file.path(hubprefix, "jyeung/data/scChiC/glmpca_outputs/same_annot_file_rerun")
jsuffix <- paste0("bincutoff_0.binskeep_", binskeep, ".byplate.szname_none.niter_", niter, ".reorder_rownames.dupfilt")

# jmark <- "H3K27me3"

dat.glmpca.annot.lst <- lapply(jmarks, function(jmark){
  fname <- paste0("glmpca.", jmark, ".", jsuffix, ".RData")
  inf.glmpca <- file.path(indir.glmpca, fname)
  print(jmark)
  print(inf.glmpca)
  load(inf.glmpca, v=T)
  dat.glmpca <- DoUmapAndLouvain(glm.out$factors, jsettings)
  dat.glmpca.annot <- left_join(dat.glmpca, subset(dat.metas[[jmark]], select = c(cell, cluster, plate, batch, cuts_in_peak, cuts_total, spikein_cuts, rowcoord, colcoord, jrep)))
  # cuts_in_peak,cuts_total,spikein_cuts,ctype,plate.orig,Cluster,batch,experi,rowcoord,colcoord,jrep
  dat.glmpca.annot$mark <- jmark
  return(dat.glmpca.annot)
})


# Load TSS counts ---------------------------------------------------------

lda.out.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  # indir <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.from_TSS.dist_", jdist))
  indir <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.from_TSS.dist_10000.RemoveDupRows/lda_outputs.count_mat_from_TSS.", jmark, ".dist_10000.K-30.binarize.FALSE"))
  assertthat::assert_that(dir.exists(indir))
  # fname <- paste0("lda_outputs.count_mat_from_TSS.", jmark, ".dist_", jdist, ".K-30.binarize.FALSE/ldaOut.count_mat_from_TSS.", jmark, ".dist_", jdist, ".K-30.Robj")
  # fname <- paste0("ldaOut.PZ-BM-rep3-", jmark, "-rep2rep3reseq.TSS.K-30.Robj")
  fname <- paste0("ldaOut.count_mat_from_TSS.", jmark, ".dist_10000.K-30.Robj")
  
  load(file.path(indir, fname), v=T)
  tm.result <- posterior(out.lda)
  tm.result <- AddTopicToTmResult(tm.result)
  return(list(count.mat = count.mat, tm.result = tm.result))
})


# Load LDA outputs from K4me3 bins to define topics -----------------------


refmark <- "H3K4me3"

for (refmark in jmarks){
  pdftmp <- file.path(outdir, paste0("TSS_", jdist, "_refmark_", refmark, ".pdf"))
  pdf(pdftmp, useDingbats = FALSE)
  
  # order topics
  topics.lst <- lapply(lda.out.lst, function(x){
    OrderTopicsByEntropy(x$tm.result)
  })
  
  topics.keep <- topics.lst[[refmark]]
  topics.keep$cluster <- NA
  genes.from.lda <- list()
  
  # check loadings
  # jtopic <- "topic7"
  topic.loadings.long <- lda.out.lst[[refmark]]$tm.result$topics %>%
    melt()
  colnames(topic.loadings.long) <- c("cell", "topic", "loadings")
  
  # topic.loadings <- data.frame(topicloadings = lda.out.lst[[refmark]]$tm.result$topics[, jtopic], stringsAsFactors = FALSE)
  jmerge <- left_join(dat.glmpca.annot.lst[[refmark]], topic.loadings.long)
  
  for (i in 1:nrow(topics.keep)){
    print(i)
    (jtopic <- topics.keep$topic[[i]])
    outlist.topic <- file.path(outdir, paste0("geneset_on_umap.binskeep_", binskeep, ".niter_", niter, ".", Sys.Date(), ".from_LDA_topics.condensed.heatmap.famousgenes.keepn_", keepn, ".refmark_", refmark, ".topic_", jtopic, ".", Sys.Date(), ".txt"))
    m <- ggplot(subset(jmerge, topic == jtopic), aes(x = umap1, y = umap2, color = loadings)) + 
      geom_point() + 
      ggtitle(refmark, jtopic) +
      theme_bw() + 
      scale_color_viridis_c() + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    print(m)
    gvec <- names(sort(lda.out.lst[[refmark]]$tm.result$terms[jtopic, ], decreasing = TRUE)[1:keepn])
    datout <- data.frame(topic = jtopic, gene = gvec, stringsAsFactors = FALSE)
    fwrite(datout, file = outlist.topic)
  }
  dev.off()
  
}



