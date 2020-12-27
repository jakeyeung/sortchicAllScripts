# Jake Yeung
# Date of Creation: 2020-12-22
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/primetime_plots/3b-make_umaps_bone_marrow_H3K27me3_cleaned.higher_spread.R
# Higher spread also add spikein information


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

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround3_umaps_and_ratios.Round1Round2.redo_after_varfilt_with_spikeins.2020-12-22.umap_spread.H3K27me3_cleaned.debug"
dir.create(outdir)
outpdf <- file.path(outdir, paste0("umap_spread_H3K27me3_cleaned.rep2rep3reseq.", Sys.Date(), ".pdf"))

hubprefix <- "/home/jyeung/hub_oudenaarden"

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

stypecols <- c("grey", "red", "blue")
cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")

# Load old metas for cluster information ----------------------------------

indir <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2/umaps_final")

dat.metas.before <- lapply(jmarks, function(jmark){
  print(jmark)
  inf.meta <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2/umaps_final/umaps_final_get_marks.", jmark, ".metadata.2020-12-17.txt"))
  print(inf.meta)
  dat.meta <- fread(inf.meta)
})

# k9me3 only need to call peaks (Ive done this already probably)
inf.meta.k9.old <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2.redo_after_varfilt_with_spikeins/cell_cluster_table_with_spikeins.H3K9me3.2020-11-18.dupfilt.txt"
dat.meta.k9.old <- subset(fread(inf.meta.k9.old), select = c(cell, cuts_in_peak, cuts_total, spikein_cuts))
dat.metas.before$H3K9me3 <- left_join(dat.metas.before$H3K9me3, dat.meta.k9.old)

# K27me3 needs to be re-updated 
indir.k27.cuts <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/BM_H3K27me3_tech_rep_merged/raw_demultiplexed/bams_split_by_cluster/bams_merged_by_cluster/counts_in_peaks_vs_nonpeaks_vs_blacklist")
infs.k27.peaks <- list.files(indir.k27.cuts, pattern = "*cuts_in_peaks.csv", full.names = FALSE)
infs.k27.genome  <- list.files(indir.k27.cuts, pattern = "*cuts_in_genome.csv", full.names = FALSE)

ctypes <- sapply(infs.k27.peaks, function(x) strsplit(x, split = "-")[[1]][[4]])
names(ctypes) <- ctypes

names(infs.k27.peaks) <- sapply(infs.k27.peaks, function(x) strsplit(x, split = "-")[[1]][[4]])
names(infs.k27.genome) <- sapply(infs.k27.genome, function(x) strsplit(x, split = "-")[[1]][[4]])

jchromos <- paste("", c(seq(19), "X", "Y"), sep = "")
dat.metas.k27.cuts <- lapply(ctypes, function(jctype){
  print(jctype)
  
  # paste0("BM_round1_round2_merged_H3K9me3_HSPCs.HiddenDomains.Lymphoid_2500_H3K9me3.cuts_in_peaks.csv")
  inf1 <- file.path(indir.k27.cuts, infs.k27.peaks[[jctype]])
  inf2 <- file.path(indir.k27.cuts, infs.k27.genome[[jctype]])
  
  # jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
  dat1 <- fread(inf1) %>%
    group_by(samp) %>%
    summarise(cuts_in_peak = sum(cuts_in_peak))
  
  dat2 <- fread(inf2) %>%
    filter(chromosome %in% jchromos) %>%
    group_by(samp) %>%
    summarise(cuts_total = sum(cuts_in_chromo))
  
  spikeinchromo <- "J02459.1"
  spikeins <- fread(inf2) %>%
    filter(chromosome == spikeinchromo) %>%
    group_by(samp) %>%
    summarise(spikein_cuts = sum(cuts_in_blacklist))
  
  # NA spikeincuts if round1, nice!
  dat.merge <- dat1 %>%
    left_join(., dat2, by = "samp")  %>%
    left_join(., spikeins, by = "samp") %>%
    ungroup() %>%
    mutate(ctype = jctype)
  
  # dat.merge.annot <- left_join(dat.merge, dat.metas.before$H3K27me3, by = c("samp" = "cell")) %>%
  #   rowwise() %>%
  #   mutate(mark = "H3K27me3",
  #          plate = ClipLast(samp, jsep = "_"))
  return(dat.merge)
}) %>%
  bind_rows()


# compare old and new quantifications 
# old l2r
head(dat.metas.before$H3K27me3)

dat.l2r.old <- dat.metas.before$H3K27me3 %>% 
  filter(!is.na(spikein_cuts)) %>%
  mutate(l2r.old = log2(cuts_in_peak / spikein_cuts)) %>%
  dplyr::select(c(cell, l2r.old, cuts_in_peak, spikein_cuts, cuts_total))

dat.l2r.new <- dat.metas.k27.cuts %>% 
  filter(!is.na(spikein_cuts)) %>%
  mutate(l2r.new = log2(cuts_in_peak / spikein_cuts)) %>%
  dplyr::select(c(samp, l2r.new, cuts_in_peak, spikein_cuts, cuts_total))

dat.l2r.merge <- left_join(dat.l2r.old, dat.l2r.new, by = c("cell" = "samp")) %>%
  rowwise() %>%
  mutate(experi = strsplit(cell, split = "-")[[1]][[3]])

m.check <- ggplot(dat.l2r.merge, aes(x = cuts_total.x, y = cuts_total.y, color = experi)) + 
  geom_point(alpha = 0.15) + 
  theme_bw() + 
  ggtitle("H3K27me3 check reseq") + 
  geom_abline(slope = 1, intercept = 0, linetype = "dotted") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# update new k27me3 
dat.meta.k27.before <- subset(dat.metas.before$H3K27me3, select = -c(cuts_in_peak, cuts_total, spikein_cuts))
dat.meta.k27.after <- dat.metas.k27.cuts %>% dplyr::rename(cell = samp)
dat.meta.k27.remerged <- left_join(dat.meta.k27.before, dat.meta.k27.after, by = c("cell")) %>%
  mutate(mark = "H3K27me3")

dat.metas.before$H3K27me3 <- dat.meta.k27.remerged


# Load glmpcas ------------------------------------------------------------

# load precalculate glmpca factors
indir.factors <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2/umaps_final_glmpca_tweak_umap.H3K27me3_rep2rep3reseq"

dat.factors <- lapply(jmarks, function(jmark){
  fname <- paste0("glmpca_factors.cleaned.", jmark, ".2020-12-23.txt")
  inf.factors <- file.path(indir.factors, fname)
  dat.factors <- fread(inf.factors) %>%
    as.data.frame()
  rnames <- dat.factors$V1
  dat.factors$V1 <- NULL
  rownames(dat.factors) <- rnames
  dat.factors <- as.matrix(dat.factors)
  return(dat.factors)
})


# Make UMAP with new settings ---------------------------------------------

dat.umaps.new.lst <- parallel::mclapply(jmarks, function(jmark){
  dat.umap <- DoUmapAndLouvain(topics.mat = dat.factors[[jmark]], jsettings = jsettings)
}, mc.cores = 4)



# Check umaps  ------------------------------------------------------------


dat.merged.lst <- lapply(jmarks, function(jmark){
  jmerge <- left_join(dat.umaps.new.lst[[jmark]], subset(dat.metas.before[[jmark]], select = c(-umap1, -umap2, -louvain)))
  return(jmerge)
})


# Flip x axes for K4me3 K9me3  --------------------------------------------


dat.merged.lst <- lapply(jmarks, function(jmark){
  jtmp <- dat.merged.lst[[jmark]]
  if (jmark %in% c("H3K4me3", "H3K9me3")){
    jtmp$umap1 <- -1 * jtmp$umap1
  }
  return(jtmp)
})


# Arrange by cluster, then by cuts to spikein ratio -----------------------

ctypes <- c("Eryths", "Bcells", "NKs", "Granulocytes", "Basophils", "pDCs", "DCs", "HSPCs")
dat.merged.lst <- lapply(jmarks, function(jmark){
  jtmp <- dat.merged.lst[[jmark]]
  jtmp$cluster <- factor(jtmp$cluster, levels = ctypes)
  jtmp <- jtmp %>%
    arrange(cluster, log2(cuts_in_peak / spikein_cuts))
  print("Peaking...")
  print(head(jtmp))
  return(jtmp)
})




# Make outputs ------------------------------------------------------------



pdf(outpdf, useDingbats = FALSE)
m.check

mlst <- lapply(jmarks, function(jmark){
  jmerge <- dat.merged.lst[[jmark]]
  m <- ggplot(jmerge, aes(x = umap1, y = umap2, color = clustercol)) + 
    geom_point(size = 0.5, alpha = 0.5) + 
    # theme_bw(size = 3) + 
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
    # theme_bw(size = 3) + 
    theme_minimal(2) + 
    scale_color_manual(values = c("blue", "grey", "red")) + 
    # scale_color_identity() + 
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
  outftxt <- file.path(outdir, paste0("cell_cluster_table_with_spikeins.", jmark, ".", Sys.Date(), ".umap_spread.final.order_by_cuts_to_spikeins.txt"))
  fwrite(dat.merged.lst[[jmark]], file = outftxt, sep = "\t", quote = FALSE)
}
