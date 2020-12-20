# Jake Yeung
# Date of Creation: 2020-12-09
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K27me3_techrepmerged/1-analyze_LDA_downstream.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(topicmodels)
library(JFuncs)
library(scchicFuncs)

library(hash)
library(igraph)
library(umap)

# Load --------------------------------------------------------------------

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

# create new metadata with filtered cells
outrds <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.BMrep2rep3reseq/varfilt_2/BM_rep2rep3reseq_H3K27me3.cleaned.varfilt_2.rds"
outtxt <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.BMrep2rep3reseq/varfilt_2/BM_rep2rep3reseq_H3K27me3.cleaned.varfilt_2.metadata.txt"
outpdf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.BMrep2rep3reseq/varfilt_2/BM_rep2rep3reseq_H3K27me3.cleaned.varfilt_2.pdf"

pdf(outpdf, useDingbats = FALSE)

# inf.lda <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_rep2rep3reseq_TSS/lda_outputs.PZ-BM-rep3-H3K27me3-rep2rep3reseq.TSS.K-30.binarize.FALSE/ldaOut.PZ-BM-rep3-H3K27me3-rep2rep3reseq.TSS.K-30.Robj"
inf.lda <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_rep2rep3reseq/lda_outputs.BM_rep2rep3reseq_H3K27me3.cleaned.K-30.binarize.FALSE/ldaOut.BM_rep2rep3reseq_H3K27me3.cleaned.K-30.Robj"
load(inf.lda, v=T)

tm.result <- posterior(out.lda)
dat.umap <- DoUmapAndLouvain(tm.result$topics, jsettings)

ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Load gene sets to celltype  ---------------------------------------------

inf.annot <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/genesets_on_umap/geneset_on_umap.binskeep_0.niter_500.2020-12-04.more_celltypes.metadta.H3K27me3.txt"
dat.annot <- fread(inf.annot)

dat.merge <- left_join(dat.umap, subset(dat.annot, select = c("cell", "cluster", "jrep"))) %>%
  rowwise() %>%
  mutate(experi = ClipLast(cell, jsep = "_"))

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

ggplot(dat.merge, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  facet_wrap(~louvain) + 
  theme_bw() + 
  scale_color_manual(values = cbPalette, na.value = "grey85") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


dat.imputed.log <- t(log2(tm.result$topics %*% tm.result$terms))

jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
dat.var <- CalculateVarAll(dat.impute.log = dat.imputed.log, jchromos = jchromos)

dat.merge2 <- left_join(dat.merge, dat.var)

ggplot(dat.merge2, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
  geom_point() + 
  theme_bw() + 
  facet_wrap(~louvain) + 
  scale_color_viridis_c(direction = -1) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.merge2, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
  geom_point() + 
  theme_bw() + 
  facet_wrap(~experi) + 
  scale_color_viridis_c(direction = -1) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


hubprefix <- "/home/jyeung/hub_oudenaarden"
inf.annot2 <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/count_tables.BMrep2rep3reseq/BM_rep2rep3reseq_H3K27me3.cleaned.l2rmin_-4tamin_0.5cutsmin_1000fraczeromax_0.25.metadata_good_cells.txt")
dat.annot2 <- fread(inf.annot2)


dat.compare <- left_join(dat.merge2, subset(dat.annot2, select = c(cell, l2r, frac.nzeros, totalcounts, spikeincounts, chromocounts)))

ggplot(dat.compare, aes(x = frac.nzeros, y = cell.var.within.sum.norm)) + 
  geom_point() 

ggplot(dat.compare, aes(x = cell.var.within.sum.norm, y = log2(chromocounts))) + 
  geom_point() 
                         
ggplot(dat.compare, aes(x = cell.var.within.sum.norm, y = l2r)) + 
  geom_point() 


ggplot(dat.merge, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  theme_bw() + 
  # scale_color_manual(values = cbPalette, na.value = "grey55") + 
  facet_wrap(~louvain) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Keep old cells but keep new cluster ------------------------------------

jlouv <- "10"
jvar.min <- 2
ggplot(dat.merge2, aes(x = cell.var.within.sum.norm, fill = cluster)) + 
  geom_density(alpha = 0.25)  + 
  theme_bw() + 
  facet_wrap(~cluster)  + 
  geom_vline(xintercept = jvar.min) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

cells.to.add <- subset(dat.merge2, louvain == jlouv & cell.var.within.sum.norm > jvar.min & is.na(cluster))$cell
print(length(cells.to.add))

# check variance of cells we kept
cells.firstbatch <- subset(dat.merge2, is.na(cluster))$cell

cells.clean <- c(cells.firstbatch, cells.to.add)


dat.merge2.showfilt <- dat.merge2 %>% mutate(highvar = cell.var.within.sum.norm > jvar.min)

ggplot(dat.merge2.showfilt, aes(x = umap1, y = umap2, color = highvar)) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Write new count tables  -------------------------------------------------

dat.filt <- subset(dat.merge2, cell.var.within.sum.norm > jvar.min)

cells.keep <- dat.filt$cell
assertthat::assert_that(length(cells.keep) == length(unique(cells.keep)))

# filter mat
mat.filt <- count.mat[, cells.keep]
dim(count.mat)
dim(mat.filt)

inf.old <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/count_tables.BMrep2rep3reseq/BM_rep2rep3reseq_H3K27me3.cleaned.l2rmin_-4tamin_0.5cutsmin_1000fraczeromax_0.25.metadata_good_cells.txt")
dat.annot.old <- fread(inf.old)

dat.annot.new <- subset(dat.annot.old, cell %in% cells.keep)

saveRDS(mat.filt, file = outrds)
fwrite(dat.annot.new, file = outtxt, sep = "\t", quote = FALSE)

dev.off()
