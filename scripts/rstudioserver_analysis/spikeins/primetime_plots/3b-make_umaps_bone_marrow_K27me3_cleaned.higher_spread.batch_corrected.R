# Jake Yeung
# Date of Creation: 2020-12-28
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/primetime_plots/3b-make_umaps_bone_marrow_K27me3_cleaned.higher_spread.batch_correct.R
# description

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

library(parallel)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$spread <- 8
jsettings$random_state <- 123

ctypes <- c("Eryths", "Bcells", "NKs", "Granulocytes", "Basophils", "pDCs", "DCs", "HSPCs")

# outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround3_umaps_and_ratios.Round1Round2.redo_after_varfilt_with_spikeins.2020-12-22.umap_spread.H3K27me3_cleaned.debug"
outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2/batch_corrected"
dir.create(outdir)
outpdf <- file.path(outdir, paste0("umap_spread_H3K27me3_cleaned.rep2rep3reseq.", Sys.Date(), ".batch_corrected.pdf"))

hubprefix <- "/home/jyeung/hub_oudenaarden"

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

stypecols <- c("grey", "red", "blue")
cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")


# Load metas  -------------------------------------------------------------

indir.meta <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28")

dat.merged.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  inf.meta <- file.path(indir.meta, paste0("metadata_batch_corrected.", jmark, ".2020-12-28.txt"))
  print(inf.meta)
  dat.meta <- fread(inf.meta)
})


# # Flip x axes for K4me3 K9me3  --------------------------------------------
# 
# 
# dat.merged.lst <- lapply(jmarks, function(jmark){
#   jtmp <- dat.merged.lst[[jmark]]
#   if (jmark %in% c("H3K4me3", "H3K9me3")){
#     jtmp$umap1 <- -1 * jtmp$umap1
#   }
#   return(jtmp)
# })
# 


# Make outputs ------------------------------------------------------------



pdf(outpdf, useDingbats = FALSE)

mlst <- lapply(jmarks, function(jmark){
  jmerge <- dat.merged.lst[[jmark]]
  m <- ggplot(jmerge, aes(x = umap1, y = umap2, color = clustercol)) + 
    geom_point(size = 0.5, alpha = 0.5) + 
    theme_minimal(2) + 
    scale_color_identity() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})

JFuncs::multiplot(mlst[[1]], mlst[[2]], mlst[[3]], mlst[[4]], cols = 4)

mlst <- lapply(jmarks, function(jmark){
  jmerge <- dat.merged.lst[[jmark]]
  jmerge$batch <- factor(jmerge$batch, levels = c("Unenriched", "Linneg", "StemCell"))
  jmerge <- jmerge %>% 
    arrange(batch)
  m <- ggplot(jmerge, aes(x = umap1, y = umap2, color = batch)) + 
    geom_point(size = 0.5, alpha = 0.25) + 
    theme_minimal(2) + 
    scale_color_manual(values = c("blue", "grey", "red")) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
})
JFuncs::multiplot(mlst[[1]], mlst[[2]], mlst[[3]], mlst[[4]], cols = 4)


# plot one with cluster color legend
m <- ggplot(dat.merged.lst$H3K4me1, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point(size = 0.5, alpha = 1) + 
  theme_minimal(12) + 
  scale_color_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m)

m <- ggplot(dat.merged.lst$H3K4me1, aes(x = umap1, y = umap2, color = batch)) + 
  geom_point(size = 0.5, alpha = 1) + 
  theme_minimal(12) + 
  scale_color_manual(values = stypecols) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m)

dev.off()

# Write new UMAP information to output ------------------------------------

m <- ggplot(dat.merged.lst$H3K4me1 %>% mutate(cluster = cluster == "Basophils"), aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point(size = 0.5, alpha = 0.5) + 
  theme_minimal() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

print(m)


for (jmark in jmarks){
  print(jmark)
  outftxt <- file.path(outdir, paste0("cell_cluster_table_with_spikeins.", jmark, ".", Sys.Date(), ".umap_spread.final.order_by_cuts_to_spikeins.batch_corrected.txt"))
  fwrite(dat.merged.lst[[jmark]], file = outftxt, sep = "\t", quote = FALSE)
}
