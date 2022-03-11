# Jake Yeung
# Date of Creation: 2022-02-12
# File: ~/projects/scchic/scripts/scripts_analysis/revisions_from_istbea/13-make_UMAPs_pretty.R
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

library(JFuncs)
library(scchicFuncs)

jmarks <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarks) <- jmarks

outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/varfilt_final"
dir.create(outdir)

# Load metadata -----------------------------------------------------------

indir.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/multinom_celltyping_varfilt_final"

infs.meta <- lapply(jmarks, function(jmark){
  file.path(indir.meta, paste0("metadata_reannotate_from_LLmat_dynamicbins.", jmark, ".2022-02-15.txt"))
})

dat.meta.lst <- lapply(infs.meta, function(jinf){
  fread(jinf) 
    # filter(ctype.from.LL != "Tcells") %>%
    # filter(!is.na(ctype.from.LL)) 
})


# Load objs ---------------------------------------------------------------

# inf.rds <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/umaps/nn_50.md_0.1.2022-02-11.rds"
# inf.rds <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/umaps/nn_75.md_0.9.2022-02-12.rds"
inf.rds <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/umaps_noNA_varfilt_final/nn_75.md_0.9.2022-02-15.rds"
umap.out.lst <- readRDS(inf.rds)


# Make UMAP  --------------------------------------------------------------

dat.umap.lst <- lapply(jmarks, function(jmark){
  umap.out <- umap.out.lst[[jmark]]
  dat.umap.louv.long <- data.frame(cell = rownames(umap.out$layout), 
                              umap1 = umap.out$layout[, 1], 
                              umap2 = umap.out$layout[, 2], stringsAsFactors = FALSE) %>%
    left_join(., dat.meta.lst[[jmark]] %>% dplyr::select(-umap1, -umap2)) %>%
    filter(!is.na(ctype.from.LL))
  
  dat.umap.louv <- umap.out
  dat.umap.long <- dat.umap.louv.long
  clstr.cname <- "louvain"
  
  cell.indx <- hash(rownames(dat.umap.louv$knn$indexes), dat.umap.louv$knn$indexes[, 1])
  cell.indx.rev <- hash(dat.umap.louv$knn$indexes[, 1], rownames(dat.umap.louv$knn$indexes))
  nr <- nrow(dat.umap.louv$knn$indexes)
  nc <- ncol(dat.umap.louv$knn$indexes)
  edgelist <- matrix(NA, nrow = nr * nc, ncol = 2)
  colnames(edgelist) <- c("from", "to")
  for (vertex.i in seq(nr)) {
    istart <- nc * (vertex.i - 1) + 1
    iend <- nc * vertex.i
    edgelist[istart:iend, 1] <- cell.indx.rev[[as.character(vertex.i)]]
    edgelist[istart:iend, 2] <- sapply(dat.umap.louv$knn$indexes[vertex.i, 
                                                                 1:nc], function(x) cell.indx.rev[[as.character(x)]])
  }
  g <- graph_from_data_frame(edgelist, directed = FALSE)
  g.out <- cluster_louvain(g, weights = NULL)
  V(g)$color <- g.out$membership
  clstr <- hash(g.out$names, g.out$membership)
  if (is.data.frame(dat.umap.long)) {
    if (clstr.cname == "louvain") {
      dat.umap.long[[clstr.cname]] <- as.character(sapply(dat.umap.long$cell, 
                                                          function(x) clstr[[as.character(x)]]))
    }
    else {
      dat.umap.long[[clstr.cname]] <- as.character(sapply(dat.umap.long$cell, 
                                                          function(x) paste0("louvain_", clstr[[as.character(x)]])))
    }
  }
  else {
    dat.umap.long <- clstr
  }
  return(dat.umap.long)
})


cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
# m.lst <- lapply(jmarks, function(jmark){
#   jdat <- dat.umap.lst[[jmark]]
#   ggplot(jdat, aes(x = umap1, y = umap2, color = ctype.from.LL)) + 
#     geom_point() + 
#     ggtitle(jmark) + 
#     # facet_wrap(~ctype.from.LL) + 
#     scale_color_manual(values = cbPalette, na.value = "grey85") + 
#     theme_bw() + 
#     theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
# })
# print(m.lst)
# 
# 
# m.batch.lst <- lapply(jmarks, function(jmark){
#   ggplot(dat.umap.lst[[jmark]], aes(x = umap1, y = umap2, color = ctype.from.LL)) + 
#     geom_point() + 
#     theme_bw() + 
#     facet_wrap(~batch) + 
#     ggtitle(jmark) + 
#     scale_color_manual(values = cbPalette, na.value = "grey85") + 
#     theme_bw() + 
#     theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
# })
# print(m.batch.lst)
# 
# 
# # cb <- palette.colors(palette = "SteppedSequential5Steps", n = 25)
# m.batch.lst <- lapply(jmarks, function(jmark){
#   ggplot(dat.umap.lst[[jmark]] %>% filter(batch == "Old"), 
#          aes(x = umap1, y = umap2, color = ctype.from.LL)) + 
#     geom_point() + 
#     theme_bw() + 
#     facet_wrap(~ctype) + 
#     ggtitle(jmark) + 
#     scale_color_manual(values = cbPalette, na.value = "grey85") + 
#     theme_bw() + 
#     theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
# })
# print(m.batch.lst$k27me3)
# 
# m.ctype.lst <- lapply(jmarks, function(jmark){
#   jdat <- dat.umap.lst[[jmark]]
#   ggplot(jdat, aes(x = umap1, y = umap2, color = ctype.from.LL)) + 
#     geom_point() + 
#     ggtitle(jmark) + 
#     facet_wrap(~ctype.from.LL) + 
#     scale_color_manual(values = cbPalette, na.value = "grey85") + 
#     theme_bw() + 
#     theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
# })
# print(m.ctype.lst)
# 
# m.ctype.lst$k27me3

# Fix colors  -------------------------------------------------------------


inf.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/first_submission_data/metadata/metadata_batch_corrected.arranged_by_lineage.H3K4me1.2021-01-02.txt"
dat.meta.colors <- fread(inf.meta)

clsts.old <- unique(dat.meta.colors$cluster)
colors.old <- unique(dat.meta.colors$clustercol)

clst2color.old <- hash::hash(clsts.old, colors.old)

# add new ctypes
print(unique(sort(dat.meta.lst$k4me1$ctype.from.LL)))

clst2color.old <- hash::hash(clsts.old, colors.old)
clst2color.old["GMP"] <- "#d88543"  # same as Granus
clst2color.old["Baso/Eo"] <- "#696969"
clst2color.old["CMP"] <- "#daac89"
clst2color.old["HSCs"] <- "#CC79A7"
clst2color.old["LT"] <- "#d093b5"
clst2color.old["ST"] <- "#fd87c9"
clst2color.old["MEP"] <- "#4c90b7"
clst2color.old["Monocytes"] <- "#d8ce40"
clst2color.old["MPPs"] <- "#d02b87"

dat.annot.lst <- lapply(dat.umap.lst, function(jdat){
  dat.annot <- jdat %>%
    rowwise() %>%
    mutate(colcode = AssignHash(ctype.from.LL, clst2color.old, NA))
})

m.batch.lst <- lapply(jmarks, function(jmark){
  m <- ggplot(dat.annot.lst[[jmark]], aes(x = umap1, y = umap2, color = batch)) + 
    geom_point(alpha = 0.5) + 
    theme_bw() + 
    ggtitle(jmark) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  return(m)
})

m.pretty.lst <- lapply(jmarks, function(jmark){
  m <- ggplot(dat.annot.lst[[jmark]], aes(x = umap1, y = umap2, color = colcode)) + 
    geom_point(alpha = 0.5) + 
    theme_bw() + 
    ggtitle(jmark) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_identity() 
  return(m)
})


m.pretty.annot.lst <- lapply(jmarks, function(jmark){
  jsub <- dat.annot.lst[[jmark]] %>%
    dplyr::select(ctype.from.LL, colcode) 
  jsub <- jsub[!duplicated(jsub), ]
  m <- ggplot(dat.annot.lst[[jmark]], aes(x = umap1, y = umap2, color = colcode)) + 
    geom_point(alpha = 1) + 
    theme_bw() + 
    ggtitle(jmark) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_identity( labels = jsub$ctype.from.LL, breaks = jsub$colcode,
                          guide = "legend") 
  return(m)
})


# By celltype, highlight the new  ------------------------------------------


for (jmark in jmarks){

  outpdf <- file.path(outdir, paste0("celltyping_umaps_primetime.", jmark, ".", Sys.Date(), ".pdf"))

  pdf(outpdf, useDingbats = FALSE)
  print(m.batch.lst[[jmark]])
  print(m.pretty.lst[[jmark]])
  print(m.pretty.annot.lst[[jmark]])

  fwrite(dat.annot.lst[[jmark]], file = file.path(outdir, paste0("metadata_with_colors.", jmark, ".", Sys.Date(), ".txt")))
  jctypes <- sort(unique(dat.annot.lst[[jmark]]$ctype.from.LL))
  for (jctype in jctypes){
    # color ctype, otherwise gray 
    jdat <- dat.annot.lst[[jmark]] %>% 
      mutate(is.ctype.LL = ctype.from.LL == jctype & batch == "New",
             is.ctype = ctype == jctype & batch == "New") %>%
      arrange(is.ctype.LL)
    m <- ggplot(jdat, aes(x = umap1, y = umap2, color = is.ctype.LL)) + 
      geom_point() + 
      theme_bw() + 
      ggtitle(paste(jmark, jctype)) + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

    print(m)
    
    m <- ggplot(jdat, aes(x = umap1, y = umap2, color = is.ctype)) + 
      geom_point() + 
      theme_bw() + 
      ggtitle(paste(jmark, jctype)) + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
    print(m)
  }
	dev.off()
  
  
}
# save output

