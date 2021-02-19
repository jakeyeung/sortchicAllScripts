# Jake Yeung
# Date of Creation: 2021-02-14
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K9me3_deeper/check_Florescu_data.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(topicmodels)
library(scchicFuncs)
library(JFuncs)
library(hash)
library(igraph)
library(umap)


hubprefix <- "/home/jyeung/hub_oudenaarden"

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

# Get metadata ------------------------------------------------------------

indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage"

dat.metas <- lapply(jmarks, function(jmark){
  print(jmark)
  fread(file.path(indir, paste0("metadata_batch_corrected.arranged_by_lineage.", jmark, ".2021-01-02.txt")))
})



# Load Florescu 50kb k9me3 data -------------------------------------------

inf.etoh <- file.path(hubprefix, "jyeung/data/dblchic/from_cluster/count_mats_for_LDA.2020-04-03.EtOH_NoTcells.varfilt/count_mat.K9m3.countcutoffmin_10000-3162-10000.TAcutoff_0.5.rds")

mat.etoh <- readRDS(inf.etoh)


# Load H3K9me3 dynamic bins  ----------------------------------------------

inf.lda <- file.path(hubprefix, "jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BM_dbl_reseq.varfilt.k4_k9_dynamic_bins/lda_outputs.count_name.H3K9me3.k4_k9_dynamic_bins.2021-01-30.K-30.binarize.FALSE/ldaOut.count_name.H3K9me3.k4_k9_dynamic_bins.2021-01-30.K-30.Robj")
load(inf.lda, v=T)



# Common rows -------------------------------------------------------------


rnames.k9 <- sapply(rownames(count.mat), function(x) paste("chr", strsplit(x, ";")[[1]][[1]], sep = ""))
rownames(count.mat) <- rnames.k9
common.rows <- intersect(rownames(mat.etoh), rnames.k9)


# Merge count mats --------------------------------------------------------

mat.merge <- cbind(mat.etoh[common.rows, ], count.mat[common.rows, ])


outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.H3K9me3_checks"
outfile <- file.path(outdir, paste0("merged_mat_H3K9me3.", Sys.Date(), ".rds"))
saveRDS(mat.merge, outfile)


# Check LL fits -----------------------------------------------------------


tm.result <- posterior(out.lda)
dat.imputed <- t(tm.result$topics %*% tm.result$terms)

rownames(dat.imputed) <- sapply(rownames(dat.imputed), function(x) paste("chr", strsplit(x, ";")[[1]][[1]], sep = ""))

dat.imputed <- dat.imputed[common.rows, ]

# Do fits -----------------------------------------------------------------

jmark <- "H3K9me3"

# avg across pseudbulk?
cnames.keep.lst <- split(x = dat.metas[[jmark]]$cell, f = dat.metas[[jmark]]$cluster)

dat.imputed.sum <- SumAcrossClusters(dat.imputed, cnames.keep.lst)
dat.imputed.sum <- do.call(cbind, dat.imputed.sum)
dat.imputed.mean <- sweep(dat.imputed.sum, MARGIN = 2, STATS = colSums(dat.imputed.sum), FUN = "/")

probs.lst.filt <- as.list(as.data.frame(dat.imputed.mean))


# # take raw
# jsub <- (dat.metas[[jmark]] %>% arrange(cuts_total))
# jcell <- jsub$cell[[1]]
# count.vec <- count.mat[, jcell]

# all.cells <- rownames(tm.result$topics); names(all.cells) <- all.cells
all.cells <- colnames(mat.etoh); names(all.cells) <- all.cells

mfits <- FitMultinoms(count.filt = as.matrix(mat.etoh[common.rows, ]), all.cells = all.cells, probs.lst.filt = probs.lst.filt, exppower = 1)

LL.ctype.lst <- mfits

LL.dat <- SummarizeMultinomFits(LL.ctype.lst, mat.etoh[common.rows, ], all.cells)



# Load objects for plotting?  ---------------------------------------------


annots.dat <- LoadCellAnnotsEtOH.MorePlates()

LL.dat.annot <- left_join(LL.dat, annots.dat, by = c("cell" = "cell.orig"))



# Load UMAP of double marks -----------------------------------------------

inf.lda.dbl <- "/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/from_cluster/LDA_outputs/ldaAnalysisBins_TopBins_autosomesOnly/lda_outputs.count_mat.K27m3xK9m3.KeepTopBins_500.KeepAllPlates.K-30.binarize.FALSE/ldaOut.count_mat.K27m3xK9m3.KeepTopBins_500.KeepAllPlates.K-30.Robj"
load(inf.lda.dbl, v=T)

out.lda.dbl <- out.lda
tm.result.dbl <- posterior(out.lda.dbl)



jsuffix <- ".FewerTopics2"
indir.clsts <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/from_cluster/LDA_downstream/clusterbytopics_TopBins_autosomesOnly", jsuffix)

inf.clsts <- file.path(indir.clsts, paste0("clusterbytopics_K9m3_TopBins_autosomesOnly_KeepAllPlates.csv"))
dat.clsts <- fread(inf.clsts) %>%
  left_join(., LL.dat.annot)

ggplot(dat.clsts, aes(x = umap1, y = umap2, color = celltype)) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.clsts, aes(x = umap1, y = umap2, color = ctype.pred)) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# dat.umap.dbl <- DoUmapAndLouvain(tm.result.dbl$topics, jsettings)
# dat.umap.dbl <- dat.clsts.lst[[jmarks.dbl]]





# Check LDA of both combined ----------------------------------------------


inf.lda.merged <- file.path(hubprefix, "jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_H3K9me3_dynamic_bins_MF_PZ/lda_outputs.merged_mat_H3K9me3.2021-02-14.K-30.binarize.FALSE/ldaOut.merged_mat_H3K9me3.2021-02-14.K-30.Robj")
load(inf.lda.merged, v=T)


tm.result.merged <- posterior(out.lda)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 1
jsettings$spread <- 8 

dat.umap.MFPZ <- DoUmapAndLouvain(tm.result.merged$topics, jsettings = jsettings)

dat.umap.MFPZ.annot <- dat.umap.MFPZ %>%
  rowwise() %>%
  mutate(batch = ifelse(grepl("BM-EtOH", cell), "MF", "PZ")) %>%
  left_join(., subset(dat.clsts, select = -c(umap1, umap2, louvain)))

ggplot(dat.umap.MFPZ, aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap.MFPZ.annot, aes(x = umap1, y = umap2, color = batch)) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap.MFPZ.annot, aes(x = umap1, y = umap2, color = celltype)) + 
  geom_point() + 
  facet_wrap(~batch) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
ggplot(dat.umap.MFPZ.annot %>% arrange(desc(batch)), aes(x = umap1, y = umap2, color = celltype)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values = cbPalette, na.value = "grey85") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap.MFPZ.annot %>% arrange(desc(batch)), aes(x = umap1, y = umap2, color = ctype.pred)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette, na.value = "grey85") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
