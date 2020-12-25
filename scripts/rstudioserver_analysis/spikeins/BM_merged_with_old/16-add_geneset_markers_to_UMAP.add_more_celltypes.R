# Jake Yeung
# Date of Creation: 2020-12-04
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged_with_old/16-add_geneset_markers_to_UMAP.add_more_celltypes.R
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

library(DescTools)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
jdist <- 10000
hubprefix <- "/home/jyeung/hub_oudenaarden"
niter <- 500
binskeep <- 0

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/genesets_on_umap"
outpdf <- file.path(outdir, paste0("geneset_on_umap.binskeep_", binskeep, ".niter_", niter, ".", Sys.Date(), ".more_celltypes.pdf"))
outbase <- file.path(outdir, paste0("geneset_on_umap.binskeep_", binskeep, ".niter_", niter, ".", Sys.Date(), ".more_celltypes.metadta"))

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




# Load celltype specific gnes  --------------------------------------------

# infrdata2 <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/rdata_robjs/de_genes_stringent_objects/de_genes_sorted_and_giladi.WithHouseKeepAndNotExpressed.FixExprsOther.RData"
infrdata2 <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/rdata_robjs/de_genes_stringent_objects/de_genes_sorted_and_giladi.WithHouseKeepAndNotExpressed.FixExprsOther.add_NK_other_cells.2020-12-04.RData"
load(infrdata2, v=T)



# UMAP with celltypes  ----------------------------------------------------

pdf(file = outpdf, useDingbats = FALSE)

stypecols <- c("grey", "red", "blue")
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
  geom_point(alpha = 0.2) + 
  theme_bw() + 
  facet_wrap(~mark) + 
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
print(names(de.ens.sorted.stringent))

jsets <- c("Erythroblast", "Bcell", "HSCs", "Neutrophil", "NKcell")
names(jsets) <- jsets

count.mat.lst.all <- lapply(lda.out.lst, function(x) x$count.mat)
rnames.common <- Reduce(intersect, lapply(count.mat.lst.all, function(x) rownames(x)))
count.mat.lst <- lapply(lda.out.lst, function(x) x$count.mat[rnames.common, ])
rnames.all.genes <- sapply(rnames.common, function(x) strsplit(x, "\\.")[[1]][[4]])
rnames.all.ens <- sapply(rnames.all.genes, function(x) AssignHash(x = x, jhash = g2e.hash, null.fill = x))

m.lst <- lapply(jsets, function(jset){
  
  eset <- as.character(de.ens.sorted.stringent[[jset]])
  
  
  rnames.filt <- names(rnames.all.ens[rnames.all.ens %in% eset])
  
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
      mutate(log2.genefrac.win = log2(Winsorize(genefrac, probs = c(0, 0.99))),
             log2.genefrac = log2(genefrac))
    m.tmp <- ggplot(jdat.merge, aes(x = umap1, y = umap2, color = log2.genefrac.win)) + 
      geom_point() + 
      theme_bw() + 
      scale_color_viridis_c(direction = 1) + 
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
    return(m.tmp)
  })
  return(m.lst)
})

dev.off()

# write metadata
print("Writing metas")
for (jmark in jmarks){
  print(jmark)
  outtxt <- paste0(outbase, ".", jmark, ".txt")
  fwrite(dat.glmpca.annot.lst[[jmark]], file = outtxt, sep = "\t")
}

