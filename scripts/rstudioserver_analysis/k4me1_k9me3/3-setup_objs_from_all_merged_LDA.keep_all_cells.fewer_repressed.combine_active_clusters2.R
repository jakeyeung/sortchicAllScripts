# Jake Yeung
# Date of Creation: 2020-10-24
# File: ~/projects/scchic/scripts/rstudioserver_analysis/k4me1_k9me3/3-setup_objs_from_all_merged_LDA.keep_all_cells.fewer_repressed.combine_active_clusters.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(topicmodels)

library(scchicFuncs)

library(hash)
library(igraph)
library(umap)

# Load all merged LDA  ----------------------------------------------------

outsuffix <- "UnionRows_KeepAllCells_FewerRepressedClusters_MergeActiveClusters2"

bad.louvs.repressed <- c("louvain4", "louvain5", "louvain2")
bad.louvs.active <- c("louvain11")
merge.louvs.active.hash <- list("louvain15" = "louvain15x2",  # eryths
                              "louvain2" = "louvain15x2",
                              "louvain13" = "louvain13x6",  # dends
                              "louvain6" = "louvain13x6",
                              "louvain7" = "louvain7x4x12",  # granus
                              "louvain4" = "louvain7x4x12",
                              "louvain12" = "louvain7x4x12",
                              "louvain5" = "louvain5x1",  # HSPCs
                              "louvain1" = "louvain5x1")



remove.na <- TRUE
jmark1 <- "H3K4me1"
jmark2 <- "H3K9me3"
jmarkdbl <- "H3K4me1xH3K9me3"

jmarks <- c(jmark1, jmark2); names(jmarks) <- jmarks

# indir.lda <- ""
# /home/jyeung/hub_oudenaarden

inf.lda1 <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld/lda_outputs.count_mat_old_merged_with_new.", jmark1, ".K-30.binarize.FALSE/ldaOut.count_mat_old_merged_with_new.", jmark1, ".K-30.Robj")
assertthat::assert_that(file.exists(inf.lda1))

inf.lda2 <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld/lda_outputs.count_mat_old_merged_with_new.", jmark2, ".K-30.binarize.FALSE/ldaOut.count_mat_old_merged_with_new.", jmark2, ".K-30.Robj")
assertthat::assert_that(file.exists(inf.lda2))

# inf.dbl.clean <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_BM_round2_all.single_match_dbl/cellfilt_binfilt/count_mat.H3K4me1xH3K9me3.match_dbl.cellfilt.binfilt.rds"
inf.dbl.clean <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_BM_round2_all.single_match_dbl/count_mat.H3K4me1xH3K9me3.match_dbl.rds"

outdir <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/H3K4me1_H3K9me3_analysis/double_staining_input/SetupObjs_AllMerged_", outsuffix)
dir.create(outdir)
outfobjs <- file.path(outdir, paste0("SetupObjs_AllMerged_", outsuffix, ".clstr_by_louvain_", jmarkdbl, ".removeNA_", remove.na, ".RData"))
assertthat::assert_that(!file.exists(outfobjs))
outpdf <- file.path(outdir, paste0("SetupObjs_AllMerged_", outsuffix, ".clstr_by_louvain_", jmarkdbl, ".removeNA_", remove.na, ".pdf"))
assertthat::assert_that(!file.exists(outpdf))

# 
# outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/H3K4me1_H3K9me3_analysis/double_staining_input/SetupObjs_AllMerged_UnionRows_KeepAllCells_FewerRepressedClusters_MergeActiveClusters"
# dir.create(outdir)
# outfobjs <- file.path(outdir, paste0("SetupObjs_AllMerged_UnionRows_KeepAllCells.clstr_by_louvain_", jmarkdbl, ".removeNA_", remove.na, ".RData"))
# assertthat::assert_that(!file.exists(outfobjs))
# outpdf <- file.path(outdir, paste0("SetupObjs_AllMerged_UnionRows_KeepAllCells.clstr_by_louvain_", jmarkdbl, ".removeNA_", remove.na, ".pdf"))
# assertthat::assert_that(!file.exists(outpdf))

pdf(outpdf, useDingbats = FALSE)

inf.lda.lst <- list(inf.lda1, inf.lda2); names(inf.lda.lst) <- jmarks

outldas <- lapply(jmarks, function(jmark){
  load(inf.lda.lst[[jmark]], v=T)  # out.lda, count.mat, count.mat.orig
  tm.result <- posterior(out.lda)
  tm.result <- AddTopicToTmResult(tm.result)
  return(list(tm.result = tm.result, count.mat = count.mat))
})

count.mat.lst <- lapply(jmarks, function(jmark){
  return(outldas[[jmark]]$count.mat)
})

tm.result.lst <- lapply(jmarks, function(jmark){
  tm.result <- outldas[[jmark]]$tm.result
  return(tm.result)
})

# tm.result.lst <- lapply(jmarks, function(jmark){
#   load(inf.lda.lst[[jmark]], v=T)  # out.lda, count.mat, count.mat.orig
#   tm.result <- posterior(out.lda)
#   tm.result <- AddTopicToTmResult(tm.result)
#   return(tm.result)
# })

# Do louvains -------------------------------------------------------------

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

dat.umap.louvs <- lapply(tm.result.lst, function(tm.result){
  jdat <- DoUmapAndLouvain(topics.mat = tm.result$topics, jsettings = jsettings)
  jdat$louvain <- paste("louvain", jdat$louvain, sep = "")
  return(jdat)
})


# Plot output -------------------------------------------------------------

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf", "#28f9ff", "#88497e", "#bcf5c3", "#86f115", "#c3c89d", "#ff010b", "#664754", "#2af022", "#3afde0", "#b9b2a8", "#f6af7c", "#c3f582", "#3b3a9e", "#71a1ee", "#df5ba4", "#3a592e", "#010233", "#686cc2", "#9b114d", "#e6e6ba", "#b9f6c5")

m.louvs <- lapply(dat.umap.louvs, function(jdat){
  m <- ggplot(jdat, mapping = aes(x = umap1, y = umap2, color = louvain)) + 
    geom_point() + 
    scale_color_manual(values = cbPalette) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
})


# Check gene regions are common or not ------------------------------------

rnames.lst <- lapply(jmarks, function(jmark){
  return(colnames(tm.result.lst[[jmark]]$terms))
})

lapply(rnames.lst, length)

length(Reduce(f = intersect, x = rnames.lst))

rnames.keep <- unique(unlist(rnames.lst))

# Create pseudobulks from imputed ------------------------------------------------------

dat.imputes <- lapply(jmarks, function(jmark){
  dat.impute <- tm.result.lst[[jmark]]$topics %*% tm.result.lst[[jmark]]$terms  # cell by genes
  # handle missing gene regions by imputing with smallest dat.impute
  pseudocount <- min(dat.impute)
  rnames.to.add <- rnames.keep[which(!rnames.keep %in% colnames(dat.impute))]
  # create mat to append
  mat.to.add <- matrix(data = pseudocount, nrow = nrow(dat.impute), ncol = length(rnames.to.add), dimnames = list(rownames(dat.impute), rnames.to.add))
  # append and return
  dat.impute.append <- cbind(dat.impute, mat.to.add)
  # renormalize
  dat.impute.append.renorm <- sweep(dat.impute.append, MARGIN = 1, STATS = rowSums(dat.impute.append), FUN = "/")
  # order rows to be the ssame order as rnames.keep
  dat.impute.append.renorm <- dat.impute.append.renorm[, rnames.keep]
  return(dat.impute.append.renorm)
})

# Create cluster table ----------------------------------------------------


dat.louv <- lapply(jmarks, function(jmark){
  dat.louv.tmp <- dat.umap.louvs[[jmark]]
  dat.louv.tmp$cluster <- dat.louv.tmp$louvain
  dat.louv.tmp$mark <- jmark
  dat.louv.tmp$cluster <- gsub(pattern = "_", replacement = "", dat.louv.tmp$cluster)  # alow only one underscore  # cal lit CLUSTER
  return(dat.louv.tmp)
})

# clean up louvs
dat.louv$H3K9me3 <- dat.louv$H3K9me3 %>%
  rowwise() %>%
  mutate(cluster = ifelse(cluster %in% bad.louvs.repressed, NA, cluster))

dat.louv$H3K4me1 <- dat.louv$H3K4me1 %>%
  rowwise() %>%
  mutate(cluster = ifelse(cluster %in% bad.louvs.active, NA, cluster))

# subset(dat.louv$H3K9me3, !cluster %in% bad.louvs.repressed)
# dat.louv$H3K4me1 <- subset(dat.louv$H3K4me1, !cluster %in% bad.louvs.active)

dat.louv$H3K4me1 <- dat.louv$H3K4me1 %>%
  rowwise() %>%
  mutate(cluster = AssignHash(cluster, merge.louvs.active.hash, null.fill = cluster))

# remove NA
if (remove.na){
  dat.louv <- lapply(dat.louv, function(jdat){
    return(subset(jdat, !is.na(cluster)))
  })
}

# plot output
m.louvs.after <- lapply(dat.louv, function(jdat){
  m <- ggplot(jdat, mapping = aes(x = umap1, y = umap2, color = cluster)) + 
    geom_point() + 
    scale_color_manual(values = cbPalette) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
})


# average across cells
cell.avgs <- lapply(jmarks, function(jmark){
  jsplit <- split(dat.louv[[jmark]], dat.louv[[jmark]]$cluster)
  return(jsplit)
})



# Get cell-averaged proportions for each louvain  -------------------------

dat.impute.repress.lst <- lapply(cell.avgs[[jmark2]], function(jsub.louvain){
  jcells <- as.character(jsub.louvain$cell)
  jrow <- colMeans(dat.imputes[[jmark2]][jcells, ])
  return(jrow)
})

# make dat.active.mat: columns are genomic regions, rows are louvains?
dat.impute.active <- lapply(cell.avgs[[jmark1]], function(jsub.louvain){
  jcells <- as.character(jsub.louvain$cell)
  jrow <- colMeans(dat.imputes[[jmark1]][jcells, ])
}) %>%
  bind_rows() %>%
  as.data.frame() %>%
  t()



# Load raw counts add missing rows ----------------------------------------

indir.raw.mats <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5234_VAN5235_VAN5236_BM/mouse/tagged_bams/dbl_stains/countTablesAndRZr1only_ByChromo.NewFilters.blfix"
infs.raw.mats <- list.files(indir.raw.mats, pattern = "*.csv", full.names = TRUE)

count.mats.raw <- lapply(infs.raw.mats, function(inf){
  mat.tmp <- ReadMatSlideWinFormat(inf)
})

lapply(count.mats.raw, dim)

# filter good cells
mat.dbl.clean <- readRDS(inf.dbl.clean)

cells.keep <- colnames(mat.dbl.clean)
# keep all cells


count.mats.raw.cellfilt <- lapply(count.mats.raw, function(jmat){
  cols.keep <- colnames(jmat) %in% cells.keep
  assertthat::assert_that(length(which(cols.keep)) > 0)
  jmat[, cols.keep]
})

# fill in empty rows

count.mats.raw.cellfilt.binfilt <- lapply(count.mats.raw.cellfilt, function(jmat){
  rows.keep <- rownames(jmat) %in% rnames.keep
  assertthat::assert_that(length(which(rows.keep)) > 0)
  jmat.rowsfilt <- jmat[rows.keep, ]
  # fill in missing rows
  rows.to.add <- rnames.keep[rnames.keep %in% rownames(jmat)]
  mat.to.add <- matrix(data = 0, nrow = length(rows.to.add), ncol = ncol(jmat.rowsfilt), dimnames = list(rows.to.add, colnames(jmat.rowsfilt)))
  jmat.rowsfilt2 <- rbind(jmat.rowsfilt, mat.to.add)
  jmat.rowsfilt2 <- jmat.rowsfilt2[rnames.keep, ]
  return(jmat.rowsfilt2)
}) 

# bind everyting together to get one countmat
lapply(count.mats.raw.cellfilt.binfilt, dim)

count.mat <- do.call(cbind, count.mats.raw.cellfilt.binfilt)
dim(count.mat)
count.dat <- list()
count.dat$counts <- count.mat


# outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/H3K4me1_H3K9me3_analyses/setup_objs_from_all_merged_union_rows"

print(length(dat.impute.repress.lst))
print(dim(dat.impute.active))
save(count.dat, dat.impute.repress.lst, dat.impute.active, dat.louv, file = outfobjs)


dev.off()