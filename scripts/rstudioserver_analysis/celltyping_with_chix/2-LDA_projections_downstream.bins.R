# Jake Yeung
# Date of Creation: 2020-08-30
# File: ~/projects/scchic/scripts/rstudioserver_analysis/celltyping_with_chix/2-LDA_projections_downstream.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(topicmodels)

library(hash)
library(igraph)
library(umap)

library(scchicFuncs)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

hubprefix <- "/home/jyeung/hub_oudenaarden"

jmark <- "H3K27me3"

outpdf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.K27me3.fromChIX/H3K27me3_cluster_renamed_from_chix_projections.pdf"

pdf(outpdf, useDingbats = FALSE)

# Load annots for double chic ---------------------------------------------



# Load annots -------------------------------------------------------------

inf.annot <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables/BM_AllMerged.", jmark, ".cell_cluster_table.txt"))
dat.annot <- fread(inf.annot)

# Load LDA projections ----------------------------------------------------

# inf.lda <- file.path(hubprefix, "jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-02-11.var_filt.UnenrichedAndAllMerged.KeepBestPlates2/lda_outputs.BM_H3K27me3_varfilt_countmat.2020-02-11.AllMerged.K-30.binarize.FALSE/ldaOut.BM_H3K27me3_varfilt_countmat.2020-02-11.AllMerged.K-30.Robj")
# load(inf.lda, v=T)
# mat.AllMerged <- count.mat
# rnames.AllMerged <- rownames(mat.AllMerged)
# load(inf.lda, v=T)


inf.proj <- file.path(hubprefix, "jyeung/data/scChiC/count_mat_B6_from_chix/LDA_outputs_projections/H3K27me3_project_on_BM_AllMerged.RData")
load(inf.proj, v=T)

topics.mat <- posterior(out.objs$out.lda)$topics
topics.mat.proj <- out.lda.predict$topics

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

umap.orig <- umap(topics.mat, jsettings)
umap.proj <- predict(umap.orig, data = topics.mat.proj)

dat.umap.proj <- data.frame(cell = rownames(umap.proj), umap1 = umap.proj[, 1], umap2 = umap.proj[, 2], type = "proj", stringsAsFactors = FALSE)

dat.umap.orig <- DoUmapAndLouvain(topics.mat, jsettings = jsettings)

dat.umap.orig$type <- "orig"


ggplot(dat.umap.orig, aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dat.umap.orig.annot <- dat.umap.orig %>%
  rowwise() %>%
  mutate(mark = GetMarkFromStr(cell),
         cond = GetCondFromSamp(samp = cell, mark = mark))

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

ggplot(dat.umap.orig.annot, aes(x = umap1, y = umap2, color = cond)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dat.umap.merged <- rbind(subset(dat.umap.orig, select = c(cell, umap1, umap2, type)), subset(dat.umap.proj, select = c(cell, umap1, umap2, type)))

ggplot(dat.umap.merged, aes(x = umap1, y = umap2, color = type)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Add annotations  --------------------------------------------------------

dat.umap.merged.annot <- dat.umap.merged %>%
  left_join(., dat.annot, by = "cell")

ggplot(dat.umap.merged.annot, aes(x = umap1.x, y = umap2.x, color = cluster)) + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~type) + 
  scale_color_manual(values = cbPalette, na.value = 'grey85')


# Load dbl annots ---------------------------------------------------------

inf.dblannot <- "/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/from_rstudio/pdfs/primetime.BM_unfixed_NAimputed.ByTopics/BM_UnfixedTopics.FinalCellClusterTable.NAimputed.single.K27m3.2020-04-28.txt"
dat.dblannot <- fread(inf.dblannot)

dat.annot.merged <- rbind(subset(dat.annot, select = c(cell, cluster)), subset(dat.dblannot, select = c(cell, cluster)))

dat.umap.merged.annot2 <- dat.umap.merged %>%
  left_join(., dat.annot.merged, by = "cell")

ggplot(dat.umap.merged.annot2 %>%
         rowwise() %>% 
         mutate(cluster = ifelse(type == "orig", NA, cluster)), 
       aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~type) + 
  scale_color_manual(values = cbPalette, na.value = 'grey85')


# Rename lineage neg with scChIX -----------------------------------------------------

# redo louvain
topics.mat.merged <- rbind(topics.mat, topics.mat.proj)

dat.umap.long.merged <- DoLouvain(topics.mat.merged, jsettings, dat.umap.merged.annot2)

ggplot(dat.umap.long.merged, aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)

ggplot(dat.umap.long.merged %>% rowwise() %>% mutate(louvain = ifelse(louvain != "9", NA, louvain)), aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette, na.value = "grey85")

jsum <- dat.umap.long.merged %>% 
  group_by(type) %>%
  mutate(ncells.total = length(cell)) %>%
  group_by(type, cluster, louvain) %>%
  summarise(ncells = length(cell)) %>%
  group_by(type, cluster) %>%
  mutate(ncells.frac = ncells / sum(ncells))

jclst <- "LinnegIsland2_topic29"
jclst <- "LinnegIsland2_topic29"
jclst <- "LinnegIsland_topic8"
ggplot(dat.umap.merged.annot %>% mutate(cluster = ifelse(cluster == jclst, cluster, NA)), aes(x = umap1.x, y = umap2.x, color = cluster)) + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~type) + 
  scale_color_manual(values = cbPalette, na.value = 'grey85')
dev.off()


# link using louvain  -----------------------------------------------------

library(reshape2)
jsum2 <- jsum %>%
  group_by(type, louvain) %>%
  filter(ncells.frac == max(ncells.frac)) %>%
  rowwise() %>%
  mutate(cluster = ifelse(ncells.frac < 0.1, paste0(cluster, "_TooFewCells"), cluster)) %>%
  dcast(data = ., formula = louvain ~ type, value.var = "cluster")


# Write new tables to output ----------------------------------------------

outf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.K27me3.fromChIX/H3K27me3_cluster_renamed_from_chix_projections.txt"
fwrite(jsum2, file = outf, sep = "\t")
