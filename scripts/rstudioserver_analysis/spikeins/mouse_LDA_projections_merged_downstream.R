# Jake Yeung
# Date of Creation: 2020-08-23
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/mouse_LDA_projections_merged_downstream.R
# 





rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

hubprefix <- "/home/jyeung/hub_oudenaarden"


jsuffix <- ".chromo2spikeinfilt"

outdir <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.merged_old_and_new", jsuffix)
dir.create(outdir)
outf <- file.path(outdir, paste0("cell_cluster_merged_with_spikein_info", jsuffix, ".txt"))
outpdf <- file.path(outdir, paste0("cell_cluster_merged_with_spikein_info", jsuffix, ".pdf"))


inf.annot <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables/BM_AllMerged.H3K4me3.cell_cluster_table.txt")
# inf.glmpca <- file.path(hubprefix, "jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/glmpcaPois_mouse_spikein/H3K4me3_BM.dev_filt.glmpcaout.penalty_5.RData")
inf.glmpca <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/glmpcaPois_mouse_spikein", jsuffix, "/H3K4me3_BM.dev_filt.glmpcaout.penalty_5.RData"))

assertthat::assert_that(file.exists(inf.annot))
assertthat::assert_that(file.exists(inf.glmpca))



hubprefix <- "/home/jyeung/hub_oudenaarden"

# inf.proj <- file.path(hubprefix, "jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_projection_onto_old/ldaOut.BM_H3K4me3_varfilt_countmat.2020-02-11.AllMerged.K-30.x.H3K4me3_padded_zeros_for_projections.RData")
inf.proj <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_projection_onto_old", jsuffix, "/ldaOut.BM_H3K4me3_varfilt_countmat.2020-02-11.AllMerged.K-30.x.H3K4me3_BM.match_rownames_with_old.RData"))
assertthat::assert_that(file.exists(inf.proj))

load(inf.proj, v=T)


pdf(outpdf, width = 1020/72, height = 815/72, useDingbats = FALSE)


jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

tm.result <- posterior(out.objs$out.lda)
topics.mat <- tm.result$topics

umap.out <- umap(topics.mat, config = jsettings)

dat.umap <- DoUmapAndLouvain(topics.mat, jsettings)

dat.umap <- dat.umap %>%
  rowwise() %>%
  mutate(mark = GetMarkFromStr(cell),
         cond = GetCondFromSamp(samp = cell, mark = mark))

ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point(alpha = 0.25) + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
ggplot(dat.umap, aes(x = umap1, y = umap2, color = cond)) + 
  geom_point(alpha = 0.25) + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)


topics.mat.proj <- out.lda.predict$topics

umap.out.proj <- predict(umap.out, data = topics.mat.proj)


dat.umap.proj <- data.frame(cell = rownames(umap.out.proj), umap1 = umap.out.proj[, 1], umap2 = umap.out.proj[, 2], stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(mark = GetMarkFromStr(cell), 
         cond = "spikein",
         louvain = "99")

dat.umap.merge <- rbind(dat.umap, dat.umap.proj)

ggplot(dat.umap.proj, aes(x = umap1, y = umap2)) + 
  geom_point()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap.merge, aes(x = umap1, y = umap2, color = cond)) + 
  geom_point(alpha = 0.5)  +  
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

ggplot(dat.umap.merge, aes(x = umap1, y = umap2, color = cond)) + 
  geom_point(alpha = 0.5)  +  
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  facet_wrap(~cond)



# Check intrachrom var ----------------------------------------------------




# add celltype annotations  ---------------------------------------------------

dat.annot <- fread(inf.annot)

# ggplot(dat.annot, aes(x = umap1, y = umap2)) + 
#   geom_point() 



# Check GLMPCA ------------------------------------------------------------

load(inf.glmpca, v=T)

head(glmpcaout$factors)

dat.umap.glmpca <- DoUmapAndLouvain(glmpcaout$factors, jsettings)

ggplot(dat.umap.glmpca, aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

head(spikeincounts)

inf.spikeins <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/quality_control_count_tables_mouse_for_lda/H3K4me3_BM.spikeins.txt")
dat.spikeins <- fread(inf.spikeins)

dat.umap.glmpca <- dat.umap.glmpca %>%
  left_join(., dat.spikeins, by = c("cell" = "samp"))

ggplot(dat.umap.glmpca, aes(x = umap1, y = umap2, color = log2(totalcounts / spikeincounts))) + 
  geom_point() +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c()

dat.dim.glmpca <- data.frame(cell = rownames(glmpcaout$factors), dim1 = glmpcaout$factors$dim1, dim2 = glmpcaout$factors$dim2, stringsAsFactors = FALSE) %>%
  left_join(., dat.spikeins, by = c("cell" = "samp"))

ggplot(dat.dim.glmpca, aes(x = dim1, y = dim2, color = log2(totalcounts / spikeincounts))) +
  geom_point()  + 
  scale_color_viridis_c() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Add celltype annotation  ------------------------------------------------

dat.umap.merge.annot <- left_join(dat.umap.merge, subset(dat.annot, select = c(-umap1, -umap2), by = "cell"))




ggplot(dat.umap.merge.annot, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("Old and new together")

ggplot(dat.umap.merge.annot, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette) + 
  facet_wrap(~cond) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  ggtitle("Old and new together")

ggplot(dat.umap.merge.annot %>% filter(cond != "spikein"), aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("Old only") 

ggplot(dat.umap.merge.annot %>% filter(cond != "spikein"), aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  ggtitle("Old only") + 
  facet_wrap(~cond)



# redo louvain and assign new cluster? 

topics.mat.merge <- rbind(topics.mat, topics.mat.proj)

dat.louvain <- DoLouvain(topics.mat = topics.mat.merge, custom.settings.louv = jsettings, dat.umap.long = dat.umap.merge.annot)

ggplot(dat.louvain, aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())




# Check neutrophil maturation?  -------------------------------------------

SelectHighestCounts <- function(jdat){
  xvec <- jdat$cluster
  jtab <- table(xvec)
  indx <- which.max(jtab)
  jcount <- jtab[indx] 
  jfrac <- jcount / sum(jtab)
  jname <- names(jtab)[indx]
  return(data.frame(cnt = jcount, frac = jfrac, cluster = jname, stringsAsFactors = FALSE))
}

dat.louvain.sum <- dat.louvain %>%
  group_by(louvain) %>%
  do(SelectHighestCounts(.))

louv2cluster <- hash::hash(as.character(dat.louvain.sum$louvain), dat.louvain.sum$cluster)

# infer cluster 
dat.louvain$cluster.infer <- sapply(dat.louvain$louvain, function(x) AssignHash(x = x, jhash = louv2cluster, null.fill = x))

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
ggplot(dat.louvain, aes(x = umap1, y = umap2, color = cluster.infer)) + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~cond) +  
  scale_color_manual(values = cbPalette)

# write new cell tablles? 
dat.louvain.final <- subset(dat.louvain, select = c(cell, umap1, umap2, louvain, cluster.infer, cond, plate)) %>%
  dplyr::rename(cluster = cluster.infer) %>%
  left_join(., dat.spikeins, by = c("cell" = "samp"))

ggplot(dat.louvain, aes(x = umap1, y = umap2, color = cluster.infer)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("Old and new together")

ggplot(dat.louvain, aes(x = umap1, y = umap2, color = cluster.infer)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette) + 
  facet_wrap(~cond) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  ggtitle("Old and new together")

ggplot(dat.louvain %>% filter(cond != "spikein"), aes(x = umap1, y = umap2, color = cluster.infer)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("Old only") 

ggplot(dat.louvain %>% filter(cond != "spikein"), aes(x = umap1, y = umap2, color = cluster.infer)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  ggtitle("Old only") + 
  facet_wrap(~cond)

dev.off()


if (!file.exists(outf)){
  fwrite(dat.louvain.final, file = outf)
} else {
  print(paste(outf, "exists, skipping"))
}



# add metadat? 



# 
# 
# table(subset(dat.louvain, grepl("Neutrophils", cluster))$louvain) / sum(table(subset(dat.louvain, grepl("Neutrophils", cluster))$louvain))
# 
# 
# table(subset(dat.louvain, louvain == "1")$cluster)

