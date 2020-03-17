# Jake Yeung
# Date of Creation: 2020-03-06
# File: ~/projects/scchic/scripts/rstudioserver_analysis/Blood_with_intestines/1-filter_low_var_cells.R
# use same filter as before


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(topicmodels)
library(hash)
library(igraph)
library(umap)

binsize <- 50000
mergesize <- 1000
bigbinsize <- 50000 * mergesize
jsystem <- "BoneMarrow"
jcutoff.ncuts.var <- 0.3
jdate <- "2020-03-06"


outdir <- "/home/jyeung/hpc/scChiC/from_rstudioserver/quality_control_postLDA_var_blood_varfilt"
dir.create(outdir)


# jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
jmarks <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarks) <- jmarks


# Load LDA outputs --------------------------------------------------------

inmain <- "/home/jyeung/hpc/intestinal_scchic/LDA_outputs/topicmodels/ldaAnalysisBins_intestinesWithBlood.2020-02-29"

out.lst <- lapply(jmarks, function(jmark){
  indir <- file.path(inmain, paste0("lda_outputs.count_mat.blood_", jmark, ".countcutoff_1000-1000-1000-1000-1000.TAcutoff_0.5.K-30.binarize.FALSE"))
  fname <- paste0("ldaOut.count_mat.blood_", jmark, ".countcutoff_1000-1000-1000-1000-1000.TAcutoff_0.5.K-30.Robj")
  inf <- file.path(indir, fname)
  assertthat::assert_that(file.exists(inf))
  
  load(inf, v=T)
  return(list(out.lda = out.lda, count.mat = count.mat))
})


# Calculate impute and raw var --------------------------------------------

# jmark <- "k4me1"

jchromos <- paste("chr", c(seq(19)), sep = "")

dat.var.impute <- lapply(jmarks, function(jmark){
  tm.result <- posterior(out.lst[[jmark]]$out.lda)
  dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms))
  dat.var <- CalculateVarAll(dat.impute.log, jchromos)
  return(dat.var)
})

dat.var.raw <- lapply(jmarks, function(jmark){
  dat.var.raw <- CalculateVarRaw(out.lst[[jmark]]$count.mat, merge.size = mergesize, 
                  chromo.exclude.grep = "^chrX|^chrY", 
                  jpseudocount = 1, jscale = 10^6, calculate.ncuts = TRUE)
  return(dat.var.raw)
})


jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

  
dat.umaps <- lapply(jmarks, function(jmark){
  dat.umap <- DoUmapAndLouvain(topics.mat = posterior(out.lst[[jmark]]$out.lda)$topics, jsettings = jsettings)
  # merge with LDA and raw var
  dat.umap <- left_join(dat.umap, dat.var.impute[[jmark]])
  dat.umap <- left_join(dat.umap, dat.var.raw[[jmark]])
  return(dat.umap)
})



# Define good cells -------------------------------------------------------

cells.keep <- lapply(jmarks, function(jm){
  subset(dat.umaps[[jm]], ncuts.var > jcutoff.ncuts.var)$cell
})
names(cells.keep) <- jmarks


# Select good cells -------------------------------------------------------

jmarksorig <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarksorig) <- jmarks

out.lda.orig <- lapply(jmarksorig, function(jm){
  inf.lda <- paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-02-11.var_filt.UnenrichedAndAllMerged.KeepBestPlates2/lda_outputs.BM_", jm, "_varfilt_countmat.2020-02-11.AllMerged.K-30.binarize.FALSE/ldaOut.BM_", jm, "_varfilt_countmat.2020-02-11.AllMerged.K-30.Robj")
  load(inf.lda)
  return(out.lda)
})

count.mat.orig <- lapply(jmarksorig, function(jm){
  inf.lda <- paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-02-11.var_filt.UnenrichedAndAllMerged.KeepBestPlates2/lda_outputs.BM_", jm, "_varfilt_countmat.2020-02-11.AllMerged.K-30.binarize.FALSE/ldaOut.BM_", jm, "_varfilt_countmat.2020-02-11.AllMerged.K-30.Robj")
  load(inf.lda)
  return(count.mat)
})


terms.orig <- lapply(out.lda.orig, function(x) x@terms)
terms.new <- lapply(out.lst, function(x) x$out.lda@terms)
terms.keep <- lapply(jmarks, function(jm){
  terms.overlap <- intersect(terms.orig[[jm]], terms.new[[jm]])
})

lapply(terms.keep, length)

# Filter good cells and good bins -----------------------------------------

count.mat.merge <- lapply(jmarks, function(jm){
  count.mat.new <- out.lst[[jm]]$count.mat[terms.keep[[jm]], cells.keep[[jm]]]
  count.mat.orig <- count.mat.orig[[jm]][terms.keep[[jm]], ]
  coount.mat.merge <- cbind(count.mat.orig, count.mat.new)
})

# Combine with BM and write to file  --------------------------------------


for (jmark in jmarks){
  outname <- paste0("countmat.Blood.", jmark, ".varfilt_", jcutoff.ncuts.var, ".", jdate, ".MergeWithBM.rds")
  count.mat.final <- count.mat.merge[[jmark]]
  print(count.mat.final[1:5, 1:5])
  print(dim(count.mat.final))
  saveRDS(count.mat.final, file = file.path(outdir, outname))
}


# Write pdfs --------------------------------------------------------------


pdfname <- paste0("plots.Blood.AllMarks.varfilt_", jcutoff.ncuts.var, ".", jdate, ".MergeWithBM.pdf")
outpdf <- file.path(outdir, pdfname)
pdf(outpdf, useDingbats=FALSE)
lapply(jmarks, function(jmark){
  m <- ggplot(dat.umaps[[jmark]], aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c(direction = -1) + ggtitle(jmark)
  print(m)
})
lapply(jmarks, function(jmark){
  m <- ggplot(dat.umaps[[jmark]], aes(x = umap1, y = umap2, color = log10(ncuts))) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c(direction = 1) + ggtitle(jmark)
  print(m)
})

lapply(jmarks, function(jmark){
  m <- ggplot(dat.umaps[[jmark]], aes(x = ncuts.var, y = cell.var.within.sum.norm)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c(direction = -1) + 
    scale_x_log10() + scale_y_log10() + 
    geom_vline(xintercept = jcutoff.ncuts.var) + ggtitle(jmark)
  print(m)
})
dev.off()