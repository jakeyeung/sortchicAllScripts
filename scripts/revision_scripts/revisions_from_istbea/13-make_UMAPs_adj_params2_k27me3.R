# Jake Yeung
# Date of Creation: 2022-04-14
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/13-make_UMAPs_adj_params2_k27me3.R
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

jsettings <- umap.defaults
jsettings[["n_neighbors"]] <- 50
jsettings[["min_dist"]] <- 0.1
jsettings[["random_state"]] <- 123
jsettings[["spread"]] <- 8

jmarks <- c("k27me3"); names(jmarks) <- jmarks
jmark <- "k27me3"


# Load metadata -----------------------------------------------------------

# indir.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/multinom_celltyping/rename"
# infs.meta <- lapply(jmarks, function(jmark){
#   # file.path(indir.meta, paste0("metadata_", jmark, ".txt"))
# })

dat.meta <- fread(paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/celltyping/repressive_cleaned_clean_eryths/metadata_celltyping_", jmark, ".dynamicbins.2022-04-14.txt"))



# Load LDA objs -----------------------------------------------------------


inf.lda <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_k27me3_clean_eryths/lda_outputs.count_mat_merged_with_old_dynbins.k27me3.2022-04-12/ldaOut.count_mat_merged_with_old_dynbins.", jmark, ".2022-04-12.Robj")
load(inf.lda, v=T)


nns <- c(18, 22, 28, 30, 33)
mds <- c(0.01, 0.9)

combos <- levels(interaction(nns, mds, sep = "_")); names(combos) <- combos

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597", "#1283fe", "#943e1b", "#01e976", "#726a5e", "#0e849a", "#d49b41", "#b6fc93", "#7b80ad", "#1a58b5", "#cb2849", "#def887")
outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/umaps_adj_params"
dir.create(outdir)

pdf(file.path(outdir, paste0("nns_mds_outputs.", Sys.Date(), ".pdf")))

m.lst.lst <- parallel::mclapply(combos, function(combo){
  print(combo)
  
  nn <- as.integer(strsplit(combo, split = "_")[[1]][[1]])
  md <- as.numeric(strsplit(combo, split = "_")[[1]][[2]])
  
  jsettings.tmp <- jsettings
  
  jsettings.tmp[["n_neighbors"]] <- nn
  jsettings.tmp[["min_dist"]] <- md
  
  
  outf.tmp <- file.path(outdir, paste0("nn_", nn, ".md_", md, ".", Sys.Date(), ".rds"))
  
  tm.result <- topicmodels::posterior(out.lda)
  topics.mat <- tm.result$topics
  umap.out <- umap(topics.mat, config = jsettings.tmp)
  
  dat.umap.long <- data.frame(cell = rownames(umap.out$layout), 
                              umap1 = umap.out$layout[, 1], 
                              umap2 = umap.out$layout[, 2], stringsAsFactors = FALSE) %>%
    left_join(., dat.meta %>% dplyr::select(-umap1, -umap2))
  
  jtitle <- paste(jmark, "NN", nn, "MD", md)
  jdat <- dat.umap.long %>%
    rowwise() %>%
    mutate(ctype.eryths = ifelse(ctype %in% c("Eryths", "MEP"), ctype, "zNot"))
  m1 <- ggplot(jdat, aes(x = umap1, y = umap2, color = ctype)) + 
    geom_point() + 
    ggtitle(jtitle) + 
    scale_color_manual(values = cbPalette, na.value = "grey85") + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  m2 <- ggplot(jdat, aes(x = umap1, y = umap2, color = ctype)) + 
    geom_point() + 
    ggtitle(jtitle) + 
    facet_wrap(~batch) + 
    scale_color_manual(values = cbPalette, na.value = "grey85") + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  
  m3 <- ggplot(jdat, aes(x = umap1, y = umap2, color = ctype.eryths)) + 
    geom_point() + 
    ggtitle(jtitle) + 
    facet_wrap(~batch) + 
    scale_color_manual(values = cbPalette, na.value = "grey85") + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  
    
  m4 <- ggplot(jdat, aes(x = umap1, y = umap2, color = ctype.eryths)) + 
    geom_point() + 
    ggtitle(jtitle) + 
    scale_color_manual(values = cbPalette, na.value = "grey85") + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  
  # save output
  saveRDS(umap.out, file = outf.tmp)
 return(list(m1, m2, m3, m4)) 
  
}, mc.cores = length(combos))
# for (nn in nns){
#   for (md in mds){

print(m.lst.lst)

#   }
# }
dev.off()



# 
# 
# # Clean up with nearest neighbor ------------------------------------------
# 
# 
# # jmark <- "k4me3"
# jmark <- "k4me1"
# cell2act <- hash::hash(dat.umap.lst[[jmark]]$cell, dat.umap.lst[[jmark]]$ctype.from.LL)
# 
# umap.out <- umap.out.lst[[jmark]]
# knn.out.act <- umap.out$knn$indexes
# rnames.orig.act <- rownames(knn.out.act); names(rnames.orig.act) <- rnames.orig.act
# rownames(knn.out.act) <- sapply(rownames(knn.out.act), function(x) AssignHash(x, cell2act, null.fill = x))
# 
# ctype.neighbors <- sapply(rnames.orig.act, function(jcell){
#   indx <- which(rnames.orig.act == jcell)
#   xsummary <- sort(table(rownames(knn.out.act)[knn.out.act[indx, ]]), decreasing = TRUE)
#   return(names(xsummary)[[1]])
# })
# 
# 
# coords.dbl.impute <- dat.umap.lst[[jmark]]
# coords.dbl.impute$ctype.impute <- sapply(coords.dbl.impute$cell, function(x){
#   ctype <- ctype.neighbors[[x]]
#   if (is.null(ctype)){
#     return(NA)
#   } else{
#     return(ctype)
#   }
# })
# 
# m1 <- ggplot(coords.dbl.impute, aes(x = umap1, y = umap2, color = ctype.impute)) + 
#   geom_point() + 
#   # facet_wrap(~ctype.impute) + 
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
# 
# m2 <- ggplot(coords.dbl.impute, aes(x = umap1, y = umap2, color = ctype.from.LL)) + 
#   geom_point() + 
#   # facet_wrap(~ctype) + 
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
# 
# m3 <- ggplot(coords.dbl.impute, aes(x = umap1, y = umap2, color = ctype)) + 
#   geom_point() + 
#   # facet_wrap(~ctype) + 
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
# 
# JFuncs::multiplot(m1, m2, m3, cols = 3)
# 
# 
