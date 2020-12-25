# Jake Yeung
# Date of Creation: 2020-12-16
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/primetime_plots/3b-make_umaps_bone_marrow_K27me3_cleaned.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)

# add also K9me3 might as well 

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

stypecols <- c("grey", "red", "blue")
# stypecols.hex <- c("#808080", "#FF0000", "#0000FF")
cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
# Load meta data  ---------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"

outdir <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2/umaps_final")

outpdf <- file.path(outdir, paste0("umaps_final_get_marks.", Sys.Date(), ".pdf"))

pdf(outpdf, useDingbats = FALSE)

dat.meta.lst <- lapply(jmarks, function(jmark){
  print(jmark)
 if (jmark == "H3K9me3"){
   # take the old one
   # inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2/BM_celltypes.bincutoff_0.binskeep_0.byplate.szname_none.niter_1000.reorder_rownames.dupfilt.2020-11-23.", jmark, ".txt"))
   inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2/BM_celltypes.bincutoff_0.binskeep_1000.byplate.szname_none.niter_500.reorder_rownames.dupfilt.2020-11-23.", jmark, ".txt"))
   dat.meta <- fread(inf)
   
   dat.round1 <- subset(dat.meta, batch != "Round2")
   dat.round2 <- subset(dat.meta, batch == "Round2")
   
   dat.round2.reannot <- dat.round2 %>%
     rowwise() %>%
     mutate(experi = ClipLast(cell, jsep = "_"),
            plate = as.numeric(SplitGetLast(experi, jsplit = "-")),
            rowcoord = AddPlateCoordinates(cell)$rowcoord,
            colcoord = AddPlateCoordinates(cell)$colcoord,
            jrep = GetRepBM(experiname = experi),
            batch = AnnotateSortFromLayoutBMall(plate, rowcoord, colcoord, jrep, jmark))
   
   dat.round1.reannot <- dat.round1 %>%
     rowwise() %>%
     mutate(experi = ClipLast(cell, jsep = "_"),
            plate = as.numeric(SplitGetLast(experi, jsplit = "-")),
            rowcoord = AddPlateCoordinates(cell)$rowcoord,
            colcoord = AddPlateCoordinates(cell)$colcoord,
            jrep = "rep1old")
   
   dat.meta.reannot <- rbind(dat.round1.reannot, dat.round2.reannot) %>%
     ungroup() %>%
     mutate(batch = gsub("LinNeg", "Linneg", batch),
            batch = gsub("LSK", "StemCell", batch))
   dat.meta.reannot$stype <- dat.meta.reannot$batch
   # rename clusters to match active marks
   dat.meta.reannot <- dat.meta.reannot %>%
     rowwise() %>%
     mutate(cluster = gsub(pattern = "Lymphoid", replacement = "Bcells", x = cluster),
            cluster = gsub(pattern = "Eryth", replacement = "Eryths", x = cluster))
   
   return(dat.meta.reannot)
   
 } else {
   inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/metadata_umap_celltype_cuts.", jmark, ".txt"))
   dat.meta <- fread(inf)
 }
  if (jmark == "H3K4me3"){
    # flip UMAP1 axis
    dat.meta$umap1 <- -1 * dat.meta$umap1
  }
  return(dat.meta)
})


# Get color code for cluster and stype  -----------------------------------

jmark.ref <- "H3K4me1"
jdat <- dat.meta.lst[[jmark.ref]]

clstr2col <- hash::hash(jdat$cluster, jdat$colorcode)
stype2col <- hash::hash(sort(unique(jdat$batch)), stypecols)

dat.meta.lst <- lapply(dat.meta.lst, function(jdat){
  jdat$clustercol <- sapply(jdat$cluster, function(x) AssignHash(x = x, jhash = clstr2col, null.fill = x))
  jdat$stypecol <- sapply(jdat$batch, function(x) AssignHash(x = x, jhash = stype2col, null.fill = x))
  return(jdat)
})

# Plot UMAPs and stype ----------------------------------------------------


for (jmark in jmarks){
  print(jmark)
  jsub <- dat.meta.lst[[jmark]] %>% arrange(batch)
  jsub$batchfac <- factor(jsub$batch, levels = c("Unenriched", "Linneg", "StemCell"))
  jsub <- jsub %>% arrange(batchfac)
  m <- ggplot(jsub, aes(x = umap1, y = umap2, color = stypecol)) + 
    geom_point(alpha = 0.2) + 
    theme_bw() + 
    ggtitle(jmark) + 
    scale_color_identity(guide = "legend", 
                         labels = unique(jsub$batch),
                         breaks = unique(jsub$stypecol)) + 
    # scale_color_manual(values = stypecols) + 
    # facet_wrap(~jrep) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m)
}


for (jmark in jmarks){
  jsub <- dat.meta.lst[[jmark]] %>% arrange(batch)
  jsub$batchfac <- factor(jsub$batch, levels = c("Unenriched", "Linneg", "StemCell"))
  jsub <- jsub %>% arrange(batchfac)
  m <- ggplot(jsub, aes(x = umap1, y = umap2, color = stypecol)) + 
    geom_point(alpha = 1) + 
    theme_bw() + 
    facet_wrap(~jrep) + 
    ggtitle(jmark) + 
    scale_color_identity(guide = "legend", 
                         labels = unique(jsub$batch),
                         breaks = unique(jsub$stypecol)) + 
    # scale_color_manual(values = stypecols) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m)
}

for (jmark in jmarks){
  m <- ggplot(dat.meta.lst[[jmark]], aes(x = umap1, y = umap2, color = clustercol)) + 
    geom_point(alpha = 0.2) + 
    theme_bw() + 
    ggtitle(jmark) + 
    scale_color_identity(guide = "legend", 
                         labels = unique(dat.meta.lst[[jmark]]$cluster)) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m)
}

for (jmark in jmarks){
  m <- ggplot(dat.meta.lst[[jmark]], aes(x = umap1, y = umap2, color = clustercol)) + 
    geom_point(alpha = 1) + 
    theme_bw() + 
    facet_wrap(~jrep) + 
    ggtitle(jmark) + 
    scale_color_identity(guide = "legend", 
                         labels = unique(dat.meta.lst[[jmark]]$cluster)) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m)
}


# Calculation fraction of stype for each cluster --------------------------

for (jmark in jmarks){
  jsum <- dat.meta.lst[[jmark]] %>%
    group_by(cluster, batch) %>%
    summarise(ncell = length(cell)) %>%
    group_by(cluster) %>%
    mutate(nfrac = ncell / sum(ncell)) %>%
    rowwise() %>% 
    mutate(stypecol = AssignHash(x = batch, jhash = stype2col, null.fill = batch))
  # order b y decreasing StemCell fraction
  jorder <- (subset(jsum, batch == "StemCell") %>% arrange(desc(nfrac)))$cluster
  jsum$cluster <- factor(jsum$cluster, levels = jorder)
  
  m <- ggplot(jsum, aes(x = cluster, y = nfrac, fill = stypecol)) + 
    # geom_col(position = "dodge") + 
    geom_col() + 
    theme_bw(12) + 
    xlab("") + ylab("Fraction") +
    ggtitle(jmark) + 
    scale_fill_identity(guide = "legend", 
                         labels = unique(jsub$batch),
                         breaks = sapply(unique(jsub$batch), function(x) AssignHash(x, stype2col, null.fill = x))) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "bottom")
  print(m)
}



dev.off()

# Write to outputs --------------------------------------------------------

for (jmark in jmarks){
  fname <- paste0("umaps_final_get_marks.", jmark, ".metadata.", Sys.Date(), ".txt")
  outf <- file.path(outdir, fname)
  fwrite(dat.meta.lst[[jmark]], outf)
}


