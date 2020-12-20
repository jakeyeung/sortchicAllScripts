# Jake Yeung
# Date of Creation: 2020-12-04
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged_with_old/16-add_geneset_markers_to_UMAP.use_genes_from_topics.R
# Use topics from LDA to define gene sets



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
jprobmin <- 0.03
jprobmax <- 0.97

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/genesets_on_umap"
outpdf <- file.path(outdir, paste0("geneset_on_umap.binskeep_", binskeep, ".niter_", niter, ".", Sys.Date(), ".from_LDA_topics.condensed.heatmap.pdf"))
outbase <- file.path(outdir, paste0("geneset_on_umap.binskeep_", binskeep, ".niter_", niter, ".", Sys.Date(), ".from_LDA_topics.metadata.condensed.heatmap"))

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
  indir <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.from_TSS.dist_", jdist))
  fname <- paste0("lda_outputs.count_mat_from_TSS.", jmark, ".dist_", jdist, ".K-30.binarize.FALSE/ldaOut.count_mat_from_TSS.", jmark, ".dist_", jdist, ".K-30.Robj")
  load(file.path(indir, fname), v=T)
  tm.result <- posterior(out.lda)
  tm.result <- AddTopicToTmResult(tm.result)
  return(list(count.mat = count.mat, tm.result = tm.result))
})


# Load LDA outputs from K4me3 bins to define topics -----------------------

keepn <- 150
refmark <- "H3K4me3"
# order topics
topics.lst <- lapply(lda.out.lst, function(x){
  OrderTopicsByEntropy(x$tm.result)
})

topics.keep <- topics.lst[[refmark]]
topics.keep$cluster <- NA
genes.from.lda <- list()

# check loadings
jtopic <- "topic7"

topic.loadings.long <- lda.out.lst[[refmark]]$tm.result$topics %>%
  melt()
colnames(topic.loadings.long) <- c("cell", "topic", "loadings")

# topic.loadings <- data.frame(topicloadings = lda.out.lst[[refmark]]$tm.result$topics[, jtopic], stringsAsFactors = FALSE)
jmerge <- left_join(dat.glmpca.annot.lst[[refmark]], topic.loadings.long)


pdftmp <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/pdfs_all/genesets_on_umap/TSS_", jdist, "_refmark_", refmark, ".pdf"))
pdf(pdftmp, useDingbats = FALSE)
for (i in 1:nrow(topics.keep)){
  print(i)
  (jtopic <- topics.keep$topic[[i]])
  m <- ggplot(subset(jmerge, topic == jtopic), aes(x = umap1, y = umap2, color = loadings)) + 
    geom_point() + 
    ggtitle(refmark, jtopic) +
    theme_bw() + 
    scale_color_viridis_c() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
}
dev.off()


jannot <- list()
jannot[["topic7"]] <- "pDCs"
jannot[["topic9"]] <- "Eryths"
jannot[["topic17"]] <- "NKs"
jannot[["topic4"]] <- "DCs"
jannot[["topic20"]] <- "Bcells"
jannot[["topic28"]] <- "HSPCs"
jannot[["topic19"]] <- "Granulocytes"
jannot[["topic12"]] <- "Basophils"
jannot[["topic8"]] <- "ErythsProg2"
jannot[["topic30"]] <- "ErythsProg"
jannot[["topic15"]] <- "Eryths2"

jannot.dat <- data.frame(topic = names(jannot), gsetname = unlist(jannot), stringsAsFactors = FALSE)

topics.filt <- names(jannot)

# get gvec
for (jtopic in topics.filt){
  gvec <- names(sort(lda.out.lst[[refmark]]$tm.result$terms[jtopic, ], decreasing = TRUE)[1:keepn])
  jset <- subset(jannot.dat, topic == jtopic)$gsetname[[1]]
  genes.from.lda[[jset]] <- sapply(gvec, function(g) strsplit(g, "\\.")[[1]][[4]])
}



# UMAP with celltypes  ----------------------------------------------------

pdf(file = outpdf, useDingbats = FALSE)

mcheck <- lapply(jmarks, function(jmark){
  m <- ggplot(dat.glmpca.annot.lst[[jmark]], 
              aes(x = umap1, y = umap2, color = batch)) + 
    geom_point(alpha = 0.2) + 
    ggtitle(jmark) + 
    theme_bw() + 
    scale_color_manual(values = stypecols) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          legend.position = "bottom") 
  return(m)
})
print(mcheck)

mctypes <- ggplot(dat.glmpca.annot.lst %>% bind_rows(), aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point(alpha = 0.2, size = 0.5) + 
  theme_bw() + 
  facet_wrap(~mark, scales = "free") + 
  scale_color_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
print(mctypes)

# mctypes <- lapply(jmarks, function(jmark){
#   m <- ggplot(dat.glmpca.annot.lst[[jmark]], 
#               aes(x = umap1, y = umap2, color = cluster)) + 
#     geom_point(alpha = 0.2) + 
#     ggtitle(jmark) + 
#     theme_bw() + 
#     scale_color_manual(values = cbPalette) + 
#     theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#           legend.position = "bottom") 
#   return(m)
# })


# Make UMAPs with raw TSS counts  -----------------------------------------


# get gene set: HSPCs
# print(names(de.ens.sorted.stringent))

# jsets <- c("Erythroblast", "Bcell", "HSCs", "Neutrophil", "NKcell")
jsets <- jannot.dat$gsetname
names(jsets) <- jsets

count.mat.lst.all <- lapply(lda.out.lst, function(x) x$count.mat)
rnames.common <- Reduce(intersect, lapply(count.mat.lst.all, function(x) rownames(x)))
count.mat.lst <- lapply(lda.out.lst, function(x) x$count.mat[rnames.common, ])
rnames.all.genes <- sapply(rnames.common, function(x) strsplit(x, "\\.")[[1]][[4]])

pseudogene.lst.lst <- lapply(jsets, function(jset){
  
  # jtopic <- subset(jannot.dat, gsetname == jset)$topic[[1]]
  eset <- genes.from.lda[[jset]]
  # eset <- as.character(de.ens.sorted.stringent[[jset]])
  
  rnames.filt <- names(rnames.all.genes[rnames.all.genes %in% eset])
  
  # rnames.filt <- names(rnames.all.ens[rnames.all.ens %in% eset])
  
  # get pseudogene across cells 
  pseudogene.lst <- lapply(count.mat.lst, function(jmat){
    totalsum <- colSums(jmat)
    genesum <- colSums(jmat[rnames.filt, ])
    genefrac <- genesum / totalsum
    jdat <- data.frame(cell = names(genefrac), genefrac = genefrac, stringsAsFactors = FALSE)
  })
  
  # add to umap and plot
  
  m.lst <- lapply(jmarks, function(jmark){
    print(jmark)
    jdat.merge <- left_join(dat.glmpca.annot.lst[[jmark]], pseudogene.lst[[jmark]]) %>%
      ungroup() %>%
      mutate(log2.genefrac.win = log2(Winsorize(genefrac, probs = c(jprobmin, jprobmax))),
             log2.genefrac = log2(genefrac))
    m.tmp <- ggplot(jdat.merge, aes(x = umap1, y = umap2, color = log2.genefrac.win)) + 
      geom_point() + 
      theme_bw() + 
      scale_color_viridis_c(direction = 1) + 
      ggtitle(jmark, jset) + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
    print(m.tmp)
    m.tmp.box <- ggplot(jdat.merge, aes(x = cluster, y = log2.genefrac.win)) + 
      geom_boxplot(outlier.shape = NA) + 
      geom_point() + 
      theme_bw() + 
      ggtitle(jmark, jset) + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
    print(m.tmp)
    m.tmp2 <- ggplot(jdat.merge, aes(x = umap1, y = umap2, color = log2.genefrac)) + 
      geom_point() + 
      theme_bw() + 
      scale_color_viridis_c(direction = 1) + 
      ggtitle(jmark, jset) + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
    print(m.tmp2)
    m.tmp2.box <- ggplot(jdat.merge, aes(x = cluster, y = log2.genefrac)) + 
      geom_boxplot(outlier.shape = NA) + 
      geom_point() + 
      theme_bw() + 
      ggtitle(jmark, jset) + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
    print(m.tmp2.box)
    return(m.tmp)
  })
  return(pseudogene.lst)
})




# make heatmap
jsets.keep <- c("pDCs", "Eryths", "NKs", "DCs", "Bcells", "HSPCs", "Granulocytes", "Basophils")
names(jsets.keep) <- jsets.keep

rnames.vec <- lapply(jsets.keep, function(jset){
  x <- names(genes.from.lda[[jset]])
  # sample(x, 50)
}) %>%
  unlist() %>%
  unique()

# cluster heatmap
# jmark.test <- "H3K4me3"
jmark.test <- "H3K4me1"
jmark.test <- "H3K27me3"

for (jmark.test in jmarks){
  cells.keep <- dat.glmpca.annot.lst[[jmark.test]]$cell
  cols.keep <- colnames(count.mat.lst[[jmark.test]]) %in% cells.keep
  rnames.keep <- rownames(count.mat.lst[[jmark.test]]) %in% rnames.vec
  count.keep <- count.mat.lst[[jmark.test]][rnames.keep, cols.keep]
  
  cellsums <- colSums(count.mat.lst[[jmark.test]])
  
  count.keep.norm <- log2(sweep(count.keep, MARGIN = 2, STATS = cellsums, FUN = "/") + 1)
  # heatmap3(as.matrix(count.keep.norm))
  
  # do pseudocuts across genes
  jmat <- lapply(pseudogene.lst.lst, function(x) x[[jmark.test]]$genefrac) %>%
    bind_cols()
  jmat <- as.matrix(jmat[, jsets.keep])
  rownames(jmat) <- colnames(count.mat.lst[[jmark.test]])
  
  jannot <- dat.metas[[jmark.test]] %>%
    arrange(cluster)
  jannot$col <- sapply(jannot$cluster, function(x) AssignHash(x, clstrs.hash, null.fill = NA))
  
  clstrs.keep <- unique(jannot$cluster)
  jmat.sorted <- jmat[jannot$cell, clstrs.keep]
  
  hm.out <- heatmap3(jmat.sorted, margins = c(5, 8), cexCol = 2, Colv = NA, Rowv = NA, 
                     RowSideColors = jannot$col, 
                     RowSideLabs = "cells", 
                     main = jmark.test,
                     labRow = FALSE, scale = "row", revC = TRUE)
  
  # hm.out <- heatmap3(jmat.sorted, margins = c(5, 8), cexCol = 2, Colv = NA, Rowv = TRUE, 
  #                    RowSideColors = jannot$col, 
  #                    RowSideLabs = "cells", 
  #                    main = jmark.test,
  #                    labRow = FALSE, scale = "row", revC = TRUE)
  
  hm.out <- heatmap3(t(jmat.sorted), margins = c(5, 8), cexRow = 2, Colv = NA, Rowv = NA, 
                     ColSideColors = jannot$col, 
                     ColSideLabs = "cells", 
                     main = jmark.test,
                     labCol = FALSE, scale = "column", revC = TRUE)
  
  
}


dev.off()







# # write metadata
# print("Writing metas")
# for (jmark in jmarks){
#   print(jmark)
#   outtxt <- paste0(outbase, ".", jmark, ".txt")
#   fwrite(dat.glmpca.annot.lst[[jmark]], file = outtxt, sep = "\t")
# }
# 
