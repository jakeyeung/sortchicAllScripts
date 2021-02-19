# Jake Yeung
# Date of Creation: 2020-12-28
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_projection_analysis/2-projections_downstream.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)
library(topicmodels)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
jsettings$spread <- 8
cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")


jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

hubprefix <- "/home/jyeung/hub_oudenaarden"

# jmark <- "H3K9me3"; prefix1 <- "old"; prefix2 <- "new"
# jmark <- "H3K27me3"; prefix1 <- "new"; prefix2 <- "old"
# jmark <- "H3K4me3"; prefix1 <- "new"; prefix2 <- "old"
jmark <- "H3K4me1"; prefix1 <- "new"; prefix2 <- "old"


outdir <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/tm_results_from_projs")
fname <- paste0("tm_result_", prefix2, "_to_", prefix1, ".", jmark, ".", Sys.Date(), ".RData")
fnamepdf <- paste0("tm_result_", prefix2, "_to_", prefix1, ".", jmark, ".", Sys.Date(), ".pdf")
outf <- file.path(outdir, fname)
outpdf <- file.path(outdir, fnamepdf)


pdf(outpdf, useDingbats = FALSE)

# Loadmetas  --------------------------------------------------------------

dat.metas <- lapply(jmarks, function(jmark){
  inf.meta <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/count_tables.BMround2.from_peaks.sitecount_mat.split_old_and_new/count_mat_from_sitecount_mat.", jmark, ".filtNAcells_allbins.from_same_annot_file.metadata.2020-12-27.txt"))
  # if (jmark != "H3K27me3"){
  #   inf.meta <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2.redo_after_varfilt_with_spikeins.2020-12-22.umap_spread.H3K27me3_cleaned/cell_cluster_table_with_spikeins.", jmark, ".2020-12-27.umap_spread.final.order_by_cuts_to_spikeins.txt"))
  # } else {
  # }
  fread(inf.meta)
})

m <- ggplot(dat.metas[[jmark]], aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  facet_wrap(~jrep) + 
  theme_bw() + 
  scale_color_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
print(m)


# Check outputs -----------------------------------------------------------

inflda <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_old_new_separate_2020-12-27.projections/BM_project_", prefix2, "_onto_", prefix1, ".", jmark, ".RData"))
load(inflda, v=T)


# Do UMAP reference -------------------------------------------------------


tm.result.orig <- posterior(out.objs$out.lda)
umap.out <- umap(tm.result.orig$topics, config = jsettings)

dat.umap.long <- data.frame(cell = rownames(umap.out[["layout"]]), umap1 = umap.out[["layout"]][, 1], umap2 = umap.out[["layout"]][, 2], stringsAsFactors = FALSE) %>%
  left_join(., subset(dat.metas[[jmark]], select = c(cell, cluster, batch, jrep)), by = c("cell"))


ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette) + 
  facet_wrap(~jrep)
  


# Add new later -----------------------------------------------------------

tm.result.new <- out.lda.predict

umap.out.pred <- predict(umap.out, data = tm.result.new$topics)

dat.umap.long.pred <- data.frame(cell = rownames(umap.out.pred), umap1 = umap.out.pred[, 1], umap2 = umap.out.pred[, 2], stringsAsFactors = FALSE) %>%
  left_join(., subset(dat.metas[[jmark]], select = c(cell, cluster, batch, jrep)), by = c("cell"))

ggplot(dat.umap.long.pred, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette) + 
  facet_wrap(~jrep)

# combine it all

dat.umap.long.merge <- rbind(dat.umap.long, dat.umap.long.pred)

ggplot(dat.umap.long.merge, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette) + 
  facet_wrap(~jrep)

ggplot(dat.umap.long.merge, aes(x = umap1, y = umap2, color = batch)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette) + 
  facet_wrap(~jrep)

ggplot(dat.umap.long.merge, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette, na.value = "grey85") 



# Save RData with terms and topics mat  ----------------------------------------

tm.result <- list()

tm.result$topics <- rbind(tm.result.orig$topics, tm.result.new$topics)
tm.result$terms <- tm.result.orig$terms


# Get count mat -----------------------------------------------------------

indir.mats <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/count_tables.BMround2.from_peaks.sitecount_mat.split_old_and_new")
assertthat::assert_that(dir.exists(indir.mats))
infs.mats <- list.files(indir.mats, pattern = paste0(".*.", jmark, ".*.rds"), full.names = TRUE)

count.mat.lst <- lapply(infs.mats, function(inf){
  count.mat.tmp <- readRDS(inf)
})

lapply(count.mat.lst, dim)

count.mat <- do.call(cbind, count.mat.lst)

# save output
save(tm.result, count.mat, file = outf)
# saveRDS(tm.result, file = outf)

dev.off()