# Jake Yeung
# Date of Creation: 2020-10-04
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/mouse/H3K27me3_merged/project_new_K27me3_onto_old.R
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

# Load data  ---------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"
inf <- file.path(hubprefix, "jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_projection_onto_old.VAN5046_VAN5230_BM/ldaOut.BM_H3K27me3_varfilt_countmat.2020-02-11.AllMerged.K-30.x.count_mat_H3K27me3_l2r_filt.2020-10-03.minl2r_0.varfilt_1.RData")
load(inf, v=T)



# # Plot old --------------------------------------------------------------

tm.result.old <- posterior(out.objs$out.lda)
tm.result.old <- AddTopicToTmResult(tm.result.old)

topics.mat.old <- tm.result.old$topics

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

umap.out <- umap(topics.mat.old, config = jsettings)

dat.umap.long <- data.frame(cell = rownames(umap.out[["layout"]]), umap1 = umap.out[["layout"]][, 1], umap2 = umap.out[["layout"]][, 2], stringsAsFactors = FALSE)
dat.umap.long <- DoLouvain(topics.mat.old, jsettings, dat.umap.long)
cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette)

dat.umap.long$batch <- "old"
dat.umap.long$stype <- sapply(dat.umap.long$cell, function(x) GetCondFromSamp(x, mark = "H3K27me3"))

dat.umap.long <- dat.umap.long %>%
  mutate(stype = gsub("Linneg", "LinNeg", stype),
         stype = gsub("StemCell", "LSK", stype))


# Plot new ----------------------------------------------------------------

topics.mat.proj <- out.lda.predict$topics


umap.pred <- predict(umap.out, data = topics.mat.proj)

dat.umap.long.pred <- data.frame(cell = rownames(umap.pred), umap1 = umap.pred[, 1], umap2 = umap.pred[, 2], stringsAsFactors = FALSE)
dat.umap.long.pred$batch <- "new"

dat.umap.long.pred <- scchicFuncs::AnnotateSortFromLayout.dat(dat.umap.long.pred)

ggplot(dat.umap.long.pred, aes(x = umap1, y = umap2)) + 
  geom_point()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


dat.umap.merge <- bind_rows(dat.umap.long %>% dplyr::select(c(cell, umap1, umap2, batch, stype)), dat.umap.long.pred %>% dplyr::select(c(cell, umap1, umap2, batch, stype)))

ggplot(dat.umap.merge, aes(x = umap1, y = umap2, color = batch)) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Check var ---------------------------------------------------------------





# Load annots -------------------------------------------------------------

inf.annot <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/GLMPCA_peaks_primetime/H3K4me1_H3K4me3_H3K27me3_glmpca_peaks_primetime.2020-09-29.H3K27me3.txt"
dat.annot <- fread(inf.annot)


dat.umap.merge.annot <- left_join(dat.umap.merge, subset(dat.annot, select = c(cell, cluster.renamed)))

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
ggplot(dat.umap.merge.annot, aes(x = umap1, y = umap2, color = cluster.renamed)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values = cbPalette, na.value = "grey85") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  facet_wrap(~batch)



# Get sorted celltype -----------------------------------------------------


dat.umap.long.pred.annot <- dat.umap.long.pred

dat.umap.merge.annot2 <- left_join(dat.umap.merge.annot, subset(dat.umap.long.pred.annot, select = -c(umap1, umap2, batch)))

ggplot(dat.umap.merge.annot2 %>% filter(batch == "old"), aes(x = umap1, y = umap2, color = stype)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values = cbPalette, na.value = "grey85") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")  + 
  facet_wrap(~batch)

ggplot(dat.umap.merge.annot2 %>% filter(batch == "old"), aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values = cbPalette, na.value = "grey85") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")  + 
  facet_wrap(~batch)

ggplot(dat.umap.merge.annot2, aes(x = umap1, y = umap2, color = stype)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values = cbPalette, na.value = "grey85") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")  + 
  facet_wrap(~batch)


# Label new cells by nearest neighbors of old  ----------------------------

UpdateNAs <- function(jsub, jknn){
  # jsub from jsub <- subset(umap.out.merge.final.annots, mark == "K27m3")
  cells.nas <- subset(jsub, is.na(cluster))$cell
  names(cells.nas) <- cells.nas
  if (length(cells.nas) == 0){
    print("Done")
    return(data.frame(NULL))
  }
  clst.new <- lapply(cells.nas, function(jcell){
    cells.keep.i <- jknn[jcell, ]
    cells.keep <- rownames(jknn)[cells.keep.i]
    jtmp <- subset(jsub, cell %in% cells.keep) %>%
      group_by(cluster) %>%
      summarise(counts = length(cell)) %>%
      ungroup() %>%
      filter(!is.na(cluster)) %>%  # guarantees convergence
      filter(counts == max(counts))
    if (nrow(jtmp) == 0){
      print("No non-NA clusters nearby... returning NA")
      return(NA)
    } else {
      # break ties randomly
      return(sample(jtmp$cluster, size = 1))
    }
  })
  clst.new.hash <- hash::hash(names(clst.new), clst.new)
  # update jsub
  jsub.new <- jsub %>%
    rowwise() %>%
    mutate(cluster = AssignHash(cell, clst.new.hash, null.fill = cluster))
  return(jsub.new)
}

dat.umap.merge.annot2$cluster <- dat.umap.merge.annot2$cluster.renamed

jsettingsForImpute <- umap.defaults
jsettingsForImpute$n_neighbors <- 50
jsettingsForImpute$min_dist <- 0.1
jsettingsForImpute$random_state <- 123

topics.mat.merge <- rbind(topics.mat.old, topics.mat.proj)

umap.all <- umap(topics.mat.merge, config = jsettingsForImpute)

dat.umap.merge.annot2.imputed <- UpdateNAs(dat.umap.merge.annot2, jknn = umap.all$knn$indexes)

ggplot(dat.umap.merge.annot2.imputed, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette, na.value = "grey85") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  # facet_wrap(~batch)
  # facet_grid(stype~batch)
  facet_grid(batch~stype)

ggplot(dat.umap.merge.annot2.imputed, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette, na.value = "grey85") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Add spikein  ------------------------------------------------------------

inf.spikeins <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_count_tables_mouse_for_lda.chromo2spikeinfilt.merged_across_runs/spikeins_dat_H3K27me3_merged.txt"
dat.spikeins <- fread(inf.spikeins)

dat.umap.merge.annot2.imputed.spikeins <- left_join(subset(dat.umap.merge.annot2.imputed, batch == "new"), subset(dat.spikeins, select = c(samp, spikeincounts, chromocounts, plate)), by = c("cell" = "samp")) %>%
  ungroup() %>%
  mutate(l2r = log2(chromocounts / spikeincounts), 
         l2r.wins = DescTools::Winsorize(l2r, probs = c(0.001, 0.999)))

ggplot(dat.umap.merge.annot2.imputed.spikeins, aes(x = umap1, y = umap2, color = log2(chromocounts / spikeincounts))) + 
  geom_point() + 
  theme_bw() + 
  scale_color_viridis_c(direction = 1) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap.merge.annot2.imputed.spikeins, aes(x = umap1, y = umap2, color = l2r.wins)) + 
  geom_point() + 
  theme_bw(24) + 
  scale_color_viridis_c(direction = 1) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

library(forcats)


ggplot(dat.umap.merge.annot2.imputed.spikeins %>% filter(!is.na(cluster)), 
       aes(x = forcats::fct_reorder(.f = cluster, .x = log2(chromocounts / spikeincounts), .fun = median, .desc = TRUE), y = log2(chromocounts/spikeincounts))) + 
  geom_boxplot() + 
  geom_point() + 
  theme_bw() + 
  xlab("") + 
  # ylab("log2(chromo to spikein)") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) 

ggplot(dat.umap.merge.annot2.imputed.spikeins %>% filter(!is.na(cluster)), 
       aes(x = forcats::fct_reorder(.f = cluster, .x = log2(chromocounts), .fun = median, .desc = TRUE), y = log2(chromocounts))) + 
  geom_boxplot() + 
  geom_point() + 
  theme_bw() + 
  xlab("") + 
  # ylab("log2(chromo to spikein)") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) 


ggplot(dat.umap.merge.annot2.imputed.spikeins %>% filter(!is.na(cluster)), aes(x = forcats::fct_reorder(.f = cluster, .x = l2r, .fun = median, .desc = TRUE), y = l2r)) + 
  geom_boxplot() + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_wrap(~experi) 


# Check glmpca ------------------------------------------------------------

# inf.glmpca <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/glmpcaPois_mouse_spikein_VAN5046_VAN5230_BM_varfilt/count_mat_H3K27me3_l2r_filt.2020-10-03.minl2r_0.varfilt_1.5.glmpcaout.penalty_1.maxiter_1000.stochastic.avagrad.tol_1e-6.devfilt_5000.varfilt_1.5.RData"
inf.glmpca <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/glmpcaPois_mouse_spikein_VAN5046_VAN5230_BM_varfilt/count_mat_H3K27me3_l2r_filt.2020-10-03.minl2r_0.varfilt_1.glmpcaout.penalty_1.maxiter_1000.stochastic.avagrad.tol_1e-6.devfilt_5000.varfilt_1.RData"
assertthat::assert_that(file.exists(inf.glmpca))
load(inf.glmpca, v=T)

dat.umap.glmpca <- DoUmapAndLouvain(glmpcaout$factors, jsettings)

dat.umap.glmpca.annot <- left_join(dat.umap.glmpca, subset(dat.umap.merge.annot2.imputed.spikeins, select = -c(umap1, umap2))) %>%
  rowwise() %>%
  mutate(l2r = log2(chromocounts / spikeincounts)) 

ggplot(dat.umap.glmpca.annot, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values = cbPalette, na.value = "grey85") +  
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~stype) 

dat.umap.glmpca.annot$l2r.wins <- DescTools::Winsorize(dat.umap.glmpca.annot$l2r, probs = c(0.01, 0.99))

ggplot(dat.umap.glmpca.annot, aes(x = umap1, y = umap2, color = l2r.wins)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  facet_wrap(~stype) 

ggplot(dat.umap.glmpca.annot %>% filter(grepl("rep3", experi)), aes(x = umap1, y = umap2, color = l2r.wins)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  facet_wrap(~stype) 

ggplot(dat.umap.glmpca.annot, aes(y = rowcoord, x = colcoord, color = l2r.wins, shape = stype)) + 
  geom_point(size = 3) + 
  theme_bw() + 
  theme(aspect.ratio=2/3, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_y_reverse() + 
  scale_color_viridis_c() + 
  facet_wrap(~experi) 



# Check lda ---------------------------------------------------------------


inf.lda <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_VAN5046_VAN5230_BM_varfilt/lda_outputs.count_mat_H3K27me3_l2r_filt.2020-10-03.minl2r_0.varfilt_1.K-30.binarize.FALSE/ldaOut.count_mat_H3K27me3_l2r_filt.2020-10-03.minl2r_0.varfilt_1.K-30.Robj"
load(inf.lda, v=T)

tm.result <- posterior(out.lda)
tm.result <- AddTopicToTmResult(tm.result, jsep = "")

dat.umap.lda <- DoUmapAndLouvain(tm.result$topics, jsettings)

dat.impute.log <- t(log2(tm.result$topics %*% tm.result$terms))

jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
dat.var <- CalculateVarAll(dat.impute.log, jchromos)

dat.umap.lda.annot <- left_join(dat.umap.lda, subset(dat.umap.merge.annot2.imputed.spikeins, select = -c(umap1, umap2))) %>%
  ungroup() %>%
  mutate(l2r = log2(chromocounts / spikeincounts),
         l2r.wins = DescTools::Winsorize(l2r, probs = c(0.01, 0.99))) %>%
  left_join(., dat.var)

dat.umap.lda.annot$stype <- factor(dat.umap.lda.annot$stype, levels = c("LSK", "LinNeg", "Unenriched"))

ggplot(dat.umap.lda.annot, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap.lda.annot, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap.lda.annot, aes(x = umap1, y = umap2, color = stype)) + 
  geom_point(size = 2.5) + 
  scale_color_manual(values = cbPalette) + 
  theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

ggplot(dat.umap.lda.annot, aes(x = umap1, y = umap2, color = stype)) + 
  geom_point() + 
  facet_wrap(~experi, nrow = 2) + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  
ggplot(dat.umap.lda.annot, aes(x = umap1, y = umap2, color = l2r.wins)) + 
  geom_point() + 
  scale_color_viridis_c() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap.lda.annot, aes(x = umap1, y = umap2, color = l2r.wins)) + 
  geom_point() + 
  facet_wrap(~stype) + 
  scale_color_viridis_c() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap.lda.annot, aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + 
  facet_wrap(~stype) + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap.lda.annot, aes(x = louvain, y = l2r.wins)) + 
  geom_boxplot() + 
  geom_point() + 
  facet_wrap(~stype) + 
  scale_color_viridis_c() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap.lda.annot, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
  geom_point() + 
  scale_color_viridis_c(direction = -1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap.lda.annot, aes(x = l2r, y = cell.var.within.sum.norm, color = stype)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ylab("Intrachrom Var") + 
  xlab("log2(chromocounts / spikeincounts)")

m0 <- ggplot(dat.umap.lda.annot, aes(x = log2(chromocounts), fill = stype)) + 
  geom_density(alpha = 0.25) + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") 
print(m0)

m1 <- ggplot(dat.umap.lda.annot, aes(x = l2r, fill = stype)) + 
  geom_density(alpha = 0.25) + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  xlab("log2(chromocounts / spikeincounts)")

m2 <- ggplot(dat.umap.lda.annot, aes(x = cell.var.within.sum.norm, fill = stype)) + 
  geom_density(alpha = 0.25) + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  xlab("Intrachrom Var")

JFuncs::multiplot(m1, m2, cols = 2)
JFuncs::multiplot(m1, m0, m2, cols = 3)

jmerge <- left_join(dat.umap.merge.annot2.imputed, dat.umap.lda.annot, by = "cell")

ggplot(jmerge, aes(x = umap1.x, y = umap2.x, color = umap1.y)) + 
  geom_point() + 
  scale_color_viridis_c() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

