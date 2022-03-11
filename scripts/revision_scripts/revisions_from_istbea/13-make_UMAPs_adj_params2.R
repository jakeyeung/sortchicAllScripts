# Jake Yeung
# Date of Creation: 2022-02-11
# File: ~/projects/scchic/scripts/scripts_analysis/revisions_from_istbea/13-make_UMAPs_adj_params.R
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

jmarks <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarks) <- jmarks




# Load metadata -----------------------------------------------------------

indir.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/multinom_celltyping/rename"
infs.meta <- lapply(jmarks, function(jmark){
  file.path(indir.meta, paste0("metadata_", jmark, ".txt"))
})

dat.meta.lst <- lapply(infs.meta, function(jinf){
  fread(jinf) %>%
    filter(ctype.from.LL != "Tcells") %>%
    filter(!is.na(ctype.from.LL))
})



# Load LDA objs -----------------------------------------------------------

infs.lda <- lapply(jmarks, function(jmark){
  inf <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_varfilt/ldaAnalysis_fripfilt_varfilt/lda_outputs.count_mat_var_filt_dynamicbins.out_dynamic_bins_new_only.varcutoff.", jmark, ".2022-01-28/ldaOut.count_mat_var_filt_dynamicbins.out_dynamic_bins_new_only.varcutoff.", jmark, ".2022-01-28.Robj")
  assertthat::assert_that(file.exists(inf))
  return(inf)
})

ldas.out <- lapply(infs.lda, function(jinf){
  load(jinf, v=T)
  # remove Tcells 
  
  return(list(count.mat = count.mat, out.lda = out.lda))
})



nns <- c(25, 50, 75)
mds <- c(0.01, 0.1, 0.9)

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/umaps_noNA"
dir.create(outdir)

pdf(file.path(outdir, paste0("nns_mds_outputs.", Sys.Date(), ".pdf")))
for (nn in nns){
  for (md in mds){
    
    
    jsettings.tmp <- jsettings
    
    jsettings.tmp[["n_neighbors"]] <- nn
    jsettings.tmp[["min_dist"]] <- md
    
    
    outf.tmp <- file.path(outdir, paste0("nn_", nn, ".md_", md, ".", Sys.Date(), ".rds"))
    
    umap.out.lst <- lapply(ldas.out, function(jout){
      out.lda <- jout$out.lda
      tm.result <- topicmodels::posterior(out.lda)
      topics.mat <- tm.result$topics
      dat.umap <- umap(topics.mat, config = jsettings.tmp)
    })
    
    
    dat.umap.lst <- lapply(jmarks, function(jmark){
      umap.out <- umap.out.lst[[jmark]]
      dat.umap.long <- data.frame(cell = rownames(umap.out$layout), 
                                  umap1 = umap.out$layout[, 1], 
                                  umap2 = umap.out$layout[, 2], stringsAsFactors = FALSE) %>%
        left_join(., dat.meta.lst[[jmark]] %>% dplyr::select(-umap1, -umap2))
      return(dat.umap.long)
    })
    
    m.lst <- lapply(jmarks, function(jmark){
      jtitle <- paste(jmark, "NN", nn, "MD", md)
      jdat <- dat.umap.lst[[jmark]]
      ggplot(jdat, aes(x = umap1, y = umap2, color = ctype.from.LL)) + 
        geom_point() + 
        ggtitle(jtitle) + 
        # facet_wrap(~ctype.from.LL) + 
        scale_color_manual(values = cbPalette, na.value = "grey85") + 
        theme_bw() + 
        theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
    })
    print(m.lst)
    
    # save output
    saveRDS(umap.out.lst, file = outf.tmp)
  }
}
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
