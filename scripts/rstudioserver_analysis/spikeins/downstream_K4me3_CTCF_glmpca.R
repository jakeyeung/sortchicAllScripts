# Jake Yeung
# Date of Creation: 2020-08-18
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/downstream_K4me3_CTCF_glmpca.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(hash)
library(igraph)
library(umap)


# Get chromosome converte -------------------------------------------------

inf.annot <- "/home/jyeung/hub_oudenaarden/jyeung/data/databases/refseq/MmRefseqTss.chromorenamed.50000.again.bed"
dat.annot <- fread(inf.annot, col.names = c("chromo", "start", "end", "origname", "strand"))
dat.annot$gene <- sapply(dat.annot$origname, function(x) strsplit(x, "\\.\\.")[[1]][[2]])
dat.annot$coord <- paste(paste(dat.annot$chromo, paste(dat.annot$start, dat.annot$end, sep = "-"), sep = ":"), dat.annot$gene, sep = ";")

rname.hash <- hash::hash(dat.annot$origname, dat.annot$coord)

# Load data  --------------------------------------------------------------

inf.lda <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/Mouse_DblChIC_CTCF_K4me3_dChIC_run-1/mats_for_LDA/LDA/lda_outputs.count_mat_K36me3_CTCF_dbl_TSS.devmin_400.K-30.binarize.FALSE/ldaOut.count_mat_K36me3_CTCF_dbl_TSS.devmin_400.K-30.Robj"
# inf.lda <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/Mouse_DblChIC_CTCF_K4me3_dChIC_run-1/mats_for_LDA/LDA/lda_outputs.count_mat_K36me3_CTCF_dbl_TSS.K-30.binarize.FALSE/ldaOut.count_mat_K36me3_CTCF_dbl_TSS.K-30.Robj"
inf.glmpca <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/Mouse_DblChIC_CTCF_K4me3_dChIC_run-1/mats_for_LDA/glmpcaout/count_mat_K36me3_CTCF_dbl_TSS.devmin_400.glmpcaout.penalty_1.by_plate.devmin_400.RData"

assertthat::assert_that(file.exists(inf.lda))
assertthat::assert_that(file.exists(inf.glmpca))

load(inf.lda, v=T)
load(inf.glmpca, v=T)

dat.spikeins <- data.frame(cuts.spikeins = spikeincounts, cell = names(spikeincounts), stringsAsFactors = FALSE)
dat.totalcuts <- data.frame(cuts.total = colSums(count.mat), cell = colnames(count.mat), stringsAsFactors = FALSE)

# rename rownames
out.lda@terms <- sapply(out.lda@terms, function(x) AssignHash(x, rname.hash, null.fill = x))

rownames(glmpcaout$loadings) <- sapply(rownames(glmpcaout$loadings), function(x) AssignHash(x, rname.hash, null.fill = x))

# Load LDA  ---------------------------------------------------------------


jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

tm.result <- posterior(out.lda)

dat.umap <- DoUmapAndLouvain(tm.result$topics, jsettings)

ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# do cell var

dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms))
# colnames(dat.impute.log) <- make.names(colnames(dat.impute.log))

jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")


# 
# cells.var.chromo.within.sum <- CalculateVarWithinChromo(dat.impute.log = dat.impute.log, jchromos = jchromos)
# cells.var.chromo.across <- CalculateVarAcrossChromo(dat.mat.log = dat.impute.log, jchromos = jchromos)
# cells.var.total <- CalculateVarTotal(dat.impute.log)

dat.var <- CalculateVarAll(dat.impute.log, jchromos = jchromos)


ggplot(left_join(dat.umap, dat.var), aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
  geom_point() + 
  scale_color_viridis_c(direction = -1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Load GLMPCA -------------------------------------------------------------

dat.umap.glmpca <- DoUmapAndLouvain(glmpcaout$factors, jsettings)

dat.glmpca <- data.frame(glmpcaout$factor, cell = rownames(glmpcaout$factor), stringsAsFactors = FALSE)

loadings.dat <- data.frame(dim1 = glmpcaout$loadings$dim1, gene = rownames(glmpcaout$loadings), stringsAsFactors = FALSE) %>%
  arrange(desc(dim1))

dat.umap.glmpca <- left_join(dat.umap.glmpca, dat.glmpca)

ggplot(left_join(dat.umap.glmpca, dat.var), aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(left_join(dat.umap.glmpca, dat.var), aes(x = dim1, y = dim2, color = louvain)) + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(left_join(dat.umap.glmpca, dat.var), aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
  geom_point() + 
  scale_color_viridis_c(direction = -1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Add grun data  ----------------------------------------------------------


inf.grun <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/bone_marrow_grun/pseudobulk_downsampled_neutrophil_trajectory.2020-05-22.WithTtest.RData"
load(inf.grun, v=T)

print(ttests.neutroprogs.out)

jgene <- "S100a8"
jgene <- "Ngp"
jgene <- "Mpo"
jgene <- "Prtn3"
jgene <- "S100a9"
jgene <- "Txn1"
jgene <- "Tal1"
jgene <- "Sox6"
jgene <- "S100a8"
jgene <- "Cebpd"
jgene <- "Bach2"

jgene <- "S100a9"
jgene <- "Elane"
jgene <- "S100a4"

jgene <- "Retnlg"
jgene <- "Sox6"
jgene <- "Krt81"
jgene <- "S100a8"
jlong.sub <- subset(jlong.grun.pbulk, grepl(jgene, gene))

ggplot(jlong.sub %>% filter(pbulk != "stem"), aes(x = pbulk, y = log2counts)) + 
  geom_point()  + ggtitle(jgene)


jrow <- subset(dat.annot, gene == jgene)$origname[[1]]
print(jrow)
assertthat::assert_that(length(jrow) != 0)

jvec <- count.mat[jrow, ]

jdat <- data.frame(cuts = jvec, cell = names(jvec), stringsAsFactors = FALSE) %>%
  left_join(., dat.spikeins) %>%
  left_join(., dat.totalcuts)

jdat.merge <- left_join(left_join(dat.umap.glmpca, dat.var), jdat)
jdat.merge.lda <- left_join(left_join(dat.umap, dat.var), jdat)

ggplot(jdat.merge, aes(x = umap1, y = umap2, color = log2(cuts / cuts.spikeins))) + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c() + ggtitle(jrow)

ggplot(jdat.merge, aes(x = dim1, y = dim2, color = log2(cuts / cuts.spikeins))) + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c() + ggtitle(jrow)

ggplot(jdat.merge, aes(x = umap1, y = umap2, color = log2(cuts / cuts.total))) + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c() + ggtitle(jrow)

ggplot(jdat.merge.lda, aes(x = umap1, y = umap2, color = log2(cuts / cuts.spikeins))) + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c() + ggtitle(jrow)

ggplot(jdat.merge.lda, aes(x = umap1, y = umap2, color = log2(cuts / cuts.total))) + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c() + ggtitle(jrow)

