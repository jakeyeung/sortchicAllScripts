# Jake Yeung
# Date of Creation: 2022-04-19
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/23-correct_batch_effect_k27me3.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(topicmodels)
library(scchicFuncs)

jstart <- Sys.time() 

outrds <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/batch_effect_corrections/mat_wide_k27me3_batch_corrected.", Sys.Date(), ".rds")
outrdsraw <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/batch_effect_corrections/mat_wide_k27me3_raw_counts.", Sys.Date(), ".rds")

AssignEffect2 <- function(batch, cluster, jeffects, jtype = c("plate", "cluster"), plateprefix = "jrep2", clusterprefix = "cluster"){
  if (jtype == "plate"){
    jname <- paste0(plateprefix, batch)
    if (jname %in% names(jeffects)){
      # print(paste("Adjusting plate nonzero because", batch, cluster))
      indx <- names(jeffects) %in% jname
      adj <- jeffects[indx]
      # adj <-  jeffects[[jeffects[jname]]]
    }  else {
      # print(paste("No adjusting plate because", batch, cluster))
      adj <- 0
    }
  } else if (jtype == "cluster"){
    jname <- paste0(plateprefix, batch, ":", clusterprefix,  cluster)
    if (jname %in% names(jeffects)){
      # print(paste("Adjusting plate nonzero because", batch, cluster))
      # print(jname)
      # print(names(jeffects))
      indx <- names(jeffects) %in% jname
      adj <- jeffects[indx]
    } else {
      # if HSCs this is zero because it is the reference cell type, therefore no interaction effect
      # print(paste("No adjusting plate because", batch, cluster))
      # print("Adjusting cluster zero")
      adj <- 0
    }
  } else {
    print(paste("jtype must be plate or cluster", jtype))
  }
  return(adj)
}



AdjustBatchEffect2 <- function(jdat, plateprefix = "jrep2", clusterprefix = "cluster", platesuffix = "old"){
  jout <- lm(data = jdat, formula = log2exprs ~ 1 + jrep2 + cluster + jrep2:cluster)
  plateeffect.i <- grep(paste0(plateprefix, platesuffix, "$"), names(jout$coefficients), value = FALSE)
  plateeffect <- jout$coefficients[plateeffect.i]
  assertthat::assert_that(length(plateeffect) > 0)
  interactioneffect.i <- grep(paste0(plateprefix, ".*.:", clusterprefix), names(jout$coefficients), value = FALSE)
  interactioneffect <- jout$coefficients[interactioneffect.i]
  assertthat::assert_that(length(interactioneffect) > 0)
  
  jdat <- jdat %>%
    rowwise() %>%
    mutate(plateadj2 = AssignEffect2(jrep2, cluster, plateeffect, jtype = "plate"),
           clstradj2 = AssignEffect2(jrep2, cluster, interactioneffect, jtype = "cluster"),
           log2exprsadj = log2exprs - plateadj2 - clstradj2)
  return(jdat)
}


jmark <- "k27me3"
jmarkold <- "H3K27me3"

# Load LDA ----------------------------------------------------------------

jsuffix2 <- "merged_with_old_dynbins"
indir <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_", jmark, "_clean_eryths")
dname <- paste0("count_mat_", jsuffix2, ".", jmark, ".2022-04-15")
inf <- file.path(indir, paste0("lda_outputs.", dname), paste0("ldaOut.", dname, ".Robj"))
load(inf, v=T)

tm <- posterior(out.lda)
dat.impute.log2 <- t(log2(tm$topics %*% tm$terms))

# Load meta ---------------------------------------------------------------

indir.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/multinom_celltyping_update_ctypes_from_LDA_k4me3_cleaned_k27me3_eryths2"
inf.meta <- file.path(indir.meta, paste0("metadata_reannotate_from_LLmat_fix_ctypes_by_batch_dynamicbins.k27me3.txt"))
dat.meta <- fread(inf.meta) %>%
  dplyr::select(cell, ctype.from.LL, batch) %>%
  rowwise() %>%
  mutate(cluster = ctype.from.LL,
         jrep2 = ifelse(batch == "New", "anew", "old"),
         cluster = ifelse(cluster == "HSCs", "aHSCs", cluster))

dat.meta$ctype.from.LL <- NULL
dat.meta$batch <- NULL
  
# 
# # Get a dynamic bins ------------------------------------------------------
# 
# inf.bins <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/first_submission_data/count_tables.BM.dynamic_bins_TSS_TES_regions/DE_bins_all_marks_padjcutoff_dists_to_TSS.annot_table.jmidbug_fixed.", jmarkold, ".2021-02-19.txt")
# dat.bins <- fread(inf.bins)
# 
# 
# # Show a region  ----------------------------------------------------------
# 
# set.seed(0) 
# jregion <- sample(dat.bins$CoordOriginal, size = 1)
# print(jregion)
# 
# # Plot hits ---------------------------------------------------------------
# 
# dat.exprs <- data.frame(cell = colnames(dat.impute.log2), log2exprs = dat.impute.log2[jregion, ], stringsAsFactors = FALSE) %>%
#   left_join(., dat.meta, by = "cell") %>%
#   rowwise() %>%
#   mutate(cluster = ctype.from.LL,
#          jrep2 = ifelse(batch == "New", "anew", "old"),
#          cluster = ifelse(cluster == "HSCs", "aHSCs", cluster))
# 
# m1 <- ggplot(dat.exprs, aes(x = cluster, y = log2exprs, fill = batch)) + 
#   geom_point() +
#   geom_boxplot() + 
#   ggtitle(jregion) + 
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

# 
# # Correct bath effects ----------------------------------------------------
# 
# 
# dat.exprs.adj <- AdjustBatchEffect2(dat.exprs)
# 
# m2 <- ggplot(dat.exprs.adj, aes(x = cluster, y = log2exprsadj, fill = batch)) + 
#   geom_point() +
#   geom_boxplot() + 
#   ggtitle(jregion) + 
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
# 
# JFuncs::multiplot(m1, m2, cols = 2)


# Go genome-wide ----------------------------------------------------------

rnames.keep <- rownames(dat.impute.log2)  # keep all
jmat.long <- dat.impute.log2[rnames.keep, ] %>%
  data.table::melt()
colnames(jmat.long) <- c("rname", "cell", "log2exprs")
jmat.long <- jmat.long %>%
  left_join(., dat.meta)

print("Correcting batch effects... multicore")

# multicore
jmat.long.lst <- split(jmat.long, jmat.long$rname)

# clean up 
rm(jmat.long, dat.impute.log2, out.lda)

jmat.wide.adj.lst <- parallel::mclapply(jmat.long.lst, function(jmat.long){
  jmat.long.adj <- jmat.long %>%
    group_by(rname) %>%
    do(AdjustBatchEffect2(.))
  
  # # check
  # jregion <- jmat.long.adj$rname[[1]] 
  # m.before <- ggplot(jmat.long.adj, aes(x = cluster, y = log2exprs, fill = jrep2)) + 
  #    geom_point() +
  #    geom_boxplot() + 
  #    ggtitle(jregion) + 
  #    theme_bw() + 
  #    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  # m.after <- ggplot(jmat.long.adj, aes(x = cluster, y = log2exprsadj, fill = jrep2)) + 
  #   geom_point() +
  #   geom_boxplot() + 
  #   ggtitle(jregion) + 
  #   theme_bw() + 
  #   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  # JFuncs::multiplot(m.before, m.after, cols = 2)
  
  # make wide vector
  jmat <- jmat.long.adj %>%
    reshape2::dcast(data = ., formula = "rname ~ cell", value.var = "log2exprsadj")
  return(jmat)
  # create 
}, mc.cores = 16) %>%
  bind_rows()

# jmat.wide.adj.lst <- readRDS("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/batch_effect_corrections/mat_long_k27me3_batch_corrected.2022-04-19.rds")

jmat.wide <- jmat.wide.adj.lst %>%
  bind_rows()
rnames.keep <- jmat.wide$rname
jmat.wide$rname <- NULL
jmat.wide <- as.matrix(jmat.wide)
rownames(jmat.wide) <- rnames.keep

# arrange colnames to match
cnames <- colnames(count.mat)
rnames <- rownames(count.mat)
jmat.wide <- jmat.wide[rnames, cnames]

assertthat::assert_that(all(rownames(count.mat) == rownames(jmat.wide)))
assertthat::assert_that(all(colnames(count.mat) == colnames(jmat.wide)))

saveRDS(jmat.wide, file = outrds)

# save raw counts
saveRDS(count.mat, file = outrdsraw)

print("Correcting batch effects... done")

# jmat.long.adj <- jmat.long %>%
#   group_by(rname) %>%
#   do(AdjustBatchEffect(.))
# print("Correcting batch effects... done")
# saveRDS(jmat.long.adj, outrds)

print(Sys.time() - jstart)
