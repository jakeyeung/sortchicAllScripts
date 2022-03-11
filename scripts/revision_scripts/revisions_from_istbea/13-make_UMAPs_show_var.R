# Jake Yeung
# Date of Creation: 2022-02-12
# File: ~/projects/scchic/scripts/scripts_analysis/revisions_from_istbea/13-make_UMAPs_show_var.R
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

library(topicmodels)

jsettings <- umap.defaults
jsettings[["n_neighbors"]] <- 75
jsettings[["min_dist"]] <- 0.9
jsettings[["random_state"]] <- 123
jsettings[["spread"]] <- 8

# jmarks <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarks) <- jmarks
jmarks <- c("k27me3", "k9me3"); names(jmarks) <- jmarks


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
  inf <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_varfilt/ldaAnalysis_fripfilt_varfilt/lda_outputs.count_mat_var_filt_allbins.out_dynamic_bins_new_only.varcutoff.", jmark, ".2022-01-28/ldaOut.count_mat_var_filt_allbins.out_dynamic_bins_new_only.varcutoff.", jmark, ".2022-01-28.Robj")
  assertthat::assert_that(file.exists(inf))
  return(inf)
})

ldas.out <- lapply(infs.lda, function(jinf){
  load(jinf, v=T)
  # remove Tcells 
  
  return(list(count.mat = count.mat, out.lda = out.lda))
})


cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/umaps_noNA_var"
dir.create(outdir)

jsettings.tmp <- jsettings
pdf(file.path(outdir, paste0("UMAP_with_var.", Sys.Date(), ".pdf")))

for (jmark in jmarks){
  print(jmark)
  
  outtxt.tmp <- file.path(outdir, paste0("metadata_with_var.", jmark, ".", Sys.Date(), ".txt"))
  
  # tm.lst <- lapply(ldas.out, function(jout){
  #   out.lda <- jout$out.lda
  #   tm.result <- topicmodels::posterior(out.lda)
  # })
  
  out.lda <- ldas.out[[jmark]]$out.lda
  tm.result <- topicmodels::posterior(out.lda)
  topics.mat <- tm.result$topics
  umap.out <- umap(topics.mat, config = jsettings.tmp)
  dat.umap.long <- data.frame(cell = rownames(umap.out$layout), 
                              umap1 = umap.out$layout[, 1], 
                              umap2 = umap.out$layout[, 2], stringsAsFactors = FALSE) %>%
    left_join(., dat.meta.lst[[jmark]] %>% dplyr::select(-umap1, -umap2))
  
    tm <- tm.result
    dat.impute.log <- t(log2(tm$topics %*% tm$terms))
    # jchromos <- sort(unique(sapply(rownames(dat.impute.log), function(x) strsplit(x, ":")[[1]][[1]])))
    # jchromos <- jchromos[!grepl("chrX|chrY", jchromos)]
    
    jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
    jchromos <- jchromos[!grepl("chrX|chrY", jchromos)]
    
	  print(jchromos)
    print(dat.impute.log[1:5, 1:5])
    dat.var <- CalculateVarAll(dat.impute.log, jchromos)
    dat.meta.new.var <- left_join(dat.umap.long, dat.var)
    
    jdat <- dat.meta.new.var

	jtitle <- jmark
    m.ctype <- ggplot(jdat, aes(x = umap1, y = umap2, color = ctype.from.LL)) + 
      geom_point() + 
      ggtitle(jtitle) + 
      # facet_wrap(~ctype.from.LL) + 
      scale_color_manual(values = cbPalette, na.value = "grey85") + 
      theme_bw() + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
    print(m.ctype)
    
    m <- ggplot(dat.meta.new.var, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
      geom_point() + 
      theme_bw() + 
      scale_color_viridis_c(direction = -1) + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    print(m)
    
    fwrite(dat.meta.new.var, file = outtxt.tmp, sep = "\t")
    
  
  # umap.out.lst <- lapply(ldas.out, function(jout){
  #   out.lda <- jout$out.lda
  #   tm.result <- topicmodels::posterior(out.lda)
  #   topics.mat <- tm.result$topics
  #   dat.umap <- umap(topics.mat, config = jsettings.tmp)
  # })
  
  # dat.umap.lst <- lapply(jmarks, function(jmark){
  #   return(dat.umap.long)
  # })
  
  
  # dat.var <- lapply(jmarks, function(jmark){
  #   tm <- tm.lst[[jmark]]
  #   dat.impute.log <- t(log2(tm$topics %*% tm$terms))
  #   jchromos <- sort(unique(sapply(rownames(dat.impute.log), function(x) strsplit(x, ":")[[1]][[1]])))
  #   jchromos <- jchromos[!grepl("chrX|chrY", jchromos)]
  #   print(dat.impute.log[1:5, 1:5])
  #   dat.var <- CalculateVarAll(dat.impute.log, jchromos)
  #   return(dat.var)
  #   dat.meta.new.var <- left_join(dat.umap.lst[[jmark]], dat.var)
  # })
  # 
  # 
  # 
  # m.lst <- lapply(jmarks, function(jmark){
  #   jtitle <- paste(jmark, "NN", nn, "MD", md)
  #   jdat <- dat.umap.lst[[jmark]]
  #   ggplot(jdat, aes(x = umap1, y = umap2, color = ctype.from.LL)) + 
  #     geom_point() + 
  #     ggtitle(jtitle) + 
  #     # facet_wrap(~ctype.from.LL) + 
  #     scale_color_manual(values = cbPalette, na.value = "grey85") + 
  #     theme_bw() + 
  #     theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  # })
  # print(m.lst)
  # 
  # m.var <- lapply(jmarks, function(jmark){
  #   dat.meta.new.var <- dat.var[[jmark]]
  #   m <- ggplot(dat.meta.new.var, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
  #     geom_point() + 
  #     theme_bw() + 
  #     scale_color_viridis_c(direction = -1) + 
  #     theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  #   print(m)
  # })
  # # # save output
  # # saveRDS(umap.out.lst, file = outf.tmp)
  
  
  
  
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
