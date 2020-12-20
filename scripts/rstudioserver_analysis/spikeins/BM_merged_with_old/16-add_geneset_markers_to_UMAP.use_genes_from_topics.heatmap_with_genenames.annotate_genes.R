# Jake Yeung
# Date of Creation: 2020-12-04
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged_with_old/16-add_geneset_markers_to_UMAP.use_genes_from_topics.R
# Use topics from LDA to define gene sets



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

library(scchicFuncs)
library(topicmodels)

library(heatmap3)

library(DescTools)
stypecols <- c("grey", "red", "blue")

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
# cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")

jdist <- 10000
hubprefix <- "/home/jyeung/hub_oudenaarden"
niter <- 500
binskeep <- 0

# winsorize constants
jprobmin <- 0.04
jprobmax <- 0.96

# keepn <- 150
keepn <- 400
refmark <- "H3K4me3"

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/genesets_on_umap"
outpdf <- file.path(outdir, paste0("geneset_on_umap.binskeep_", binskeep, ".niter_", niter, ".", Sys.Date(), ".from_LDA_topics.condensed.heatmap.famousgenes.keepn_", keepn, ".refmark_", refmark, ".", Sys.Date(), ".pdf"))
outlist <- file.path(outdir, paste0("geneset_on_umap.binskeep_", binskeep, ".niter_", niter, ".", Sys.Date(), ".from_LDA_topics.condensed.heatmap.famousgenes.keepn_", keepn, ".refmark_", refmark, ".", Sys.Date(), ".txt"))
# outbase <- file.path(outdir, paste0("geneset_on_umap.binskeep_", binskeep, ".niter_", niter, ".", Sys.Date(), ".from_LDA_topics.metadata.condensed.heatmap"))

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks


# Add spikein annots  -----------------------------------------------------

indir.meta <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2.redo_after_varfilt_with_spikeins"
jdate <- "2020-11-18"
dat.metas.init <- lapply(jmarks, function(jmark){
  fname <- paste0("cell_cluster_table_with_spikeins.", jmark, ".", jdate, ".dupfilt.txt")
  inf.meta <- file.path(indir.meta, fname)
  dat.tmp <- fread(inf.meta)
  dat.tmp <- subset(dat.tmp, select = -c(umap1, umap2))
})

# fix stype for round1 and round2 
dat.round1.lst <- lapply(dat.metas.init, function(x) subset(x, batch != "Round2"))
dat.round2.lst <- lapply(dat.metas.init, function(x) subset(x, batch == "Round2"))

dat.round2.reannot.lst <- lapply(jmarks, function(jmark){
  jreannot <- dat.round2.lst[[jmark]] %>%
    rowwise() %>%
    mutate(experi = ClipLast(cell, jsep = "_"),
           plate = as.numeric(SplitGetLast(experi, jsplit = "-")),
           rowcoord = AddPlateCoordinates(cell)$rowcoord,
           colcoord = AddPlateCoordinates(cell)$colcoord,
           jrep = GetRepBM(experiname = experi), 
           batch = AnnotateSortFromLayoutBMall(plate, rowcoord, colcoord, jrep, jmark))
})
dat.round1.reannot.lst <- lapply(jmarks, function(jmark){
  jreannot <- dat.round1.lst[[jmark]] %>%
    rowwise() %>%
    mutate(experi = ClipLast(cell, jsep = "_"),
           plate = as.numeric(SplitGetLast(experi, jsplit = "-")),
           rowcoord = AddPlateCoordinates(cell)$rowcoord,
           colcoord = AddPlateCoordinates(cell)$colcoord,
           jrep = "rep1old")
})

dat.metas <- lapply(jmarks, function(jmark){
  jreannot1 <- dat.round1.reannot.lst[[jmark]]
  jreannot2 <- dat.round2.reannot.lst[[jmark]]
  jreannot <- rbind(jreannot1, jreannot2) %>%
    ungroup() %>%
    mutate(batch = gsub("LinNeg", "Linneg", batch),
           batch = gsub("LSK", "StemCell", batch))
  jreannot$stype <- jreannot$batch
  return(jreannot)
})

# set up colors
clstrs <- sort(unique(dat.metas$H3K4me1$cluster))  # use reference mark with all the celltypes
clstrs.col <- cbPalette[1:length(clstrs)]
clstrs.hash <- hash::hash(clstrs, clstrs.col)
head(unique(dat.round1.lst$H3K4me1$stype))
head(unique(dat.round2.lst$H3K4me1$stype))

# Load GLMPCA -------------------------------------------------------------

binskeep <- 0
iter <- 500

indir.glmpca <- file.path(hubprefix, "jyeung/data/scChiC/glmpca_outputs/same_annot_file_rerun")
jsuffix <- paste0("bincutoff_0.binskeep_", binskeep, ".byplate.szname_none.niter_", niter, ".reorder_rownames.dupfilt")

# jmark <- "H3K27me3"

dat.glmpca.annot.lst <- lapply(jmarks, function(jmark){
  fname <- paste0("glmpca.", jmark, ".", jsuffix, ".RData")
  inf.glmpca <- file.path(indir.glmpca, fname)
  print(jmark)
  print(inf.glmpca)
  load(inf.glmpca, v=T)
  dat.glmpca <- DoUmapAndLouvain(glm.out$factors, jsettings)
  dat.glmpca.annot <- left_join(dat.glmpca, subset(dat.metas[[jmark]], select = c(cell, cluster, plate, batch, cuts_in_peak, cuts_total, spikein_cuts, rowcoord, colcoord, jrep)))
  # cuts_in_peak,cuts_total,spikein_cuts,ctype,plate.orig,Cluster,batch,experi,rowcoord,colcoord,jrep
  dat.glmpca.annot$mark <- jmark
  return(dat.glmpca.annot)
})


# Load TSS counts ---------------------------------------------------------

lda.out.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  indir <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.from_TSS.dist_", jdist))
  fname <- paste0("lda_outputs.count_mat_from_TSS.", jmark, ".dist_", jdist, ".K-30.binarize.FALSE/ldaOut.count_mat_from_TSS.", jmark, ".dist_", jdist, ".K-30.Robj")
  load(file.path(indir, fname), v=T)
  tm.result <- posterior(out.lda)
  tm.result <- AddTopicToTmResult(tm.result)
  return(list(count.mat = count.mat, tm.result = tm.result))
})


# Load genes --------------------------------------------------------------

inf.genes <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/pdfs_all/genesets_on_umap/geneset_on_umap.binskeep_0.niter_500.2020-12-09.from_LDA_topics.condensed.heatmap.famousgenes.keepn_400.refmark_H3K4me3.2020-12-09.txt")
dat.genes <- fread(inf.genes)



# Compare H3K4me1 imputed  ------------------------------------------------



jmark.test <- "H3K27me3"
# dat.imputed.lda <- t(log2(lda.out.lst[[jmark.test]]$tm.result$topics %*% lda.out.lst[[jmark.test]]$tm.result$terms)) dat.imputed.lda <- t(log2(lda.out.lst[[jmark.test]]$tm.result$topics %*% lda.out.lst[[jmark.test]]$tm.result$terms))
dat.imputed.lda <- t(lda.out.lst[[jmark.test]]$tm.result$topics %*% lda.out.lst[[jmark.test]]$tm.result$terms)

# indir.glmpca.tss <- file.path(hubprefix, "jyeung/data/scChiC/glmpca_outputs/TSS_filtgenes")
# assertthat::assert_that(dir.exists(indir.glmpca.tss))
# 
# # inf.glmpca.tss <- file.path(indir.glmpca.tss, paste0("glmpca.", jmark.test, ".bincutoff_0.binskeep_0.platename_jrep.szname_none.niter_1000.reorder_rownames.dupfilt.RData"))
# inf.glmpca.tss <- file.path(indir.glmpca.tss, paste0("glmpca.", jmark.test, ".bincutoff_0.binskeep_0.platename_plate.szname_none.niter_1000.reorder_rownames.dupfilt.RData"))
# print(inf.glmpca.tss)
# assertthat::assert_that(file.exists(inf.glmpca.tss))
# 
# load(inf.glmpca.tss, v=T)
# dat.imputed.glmpca <- as.matrix(glm.out$loadings) %*% as.matrix(t(glm.out$factors))


# Show heatmap  -----------------------------------------------------------

# order cells 


jsub.annot <- dat.glmpca.annot.lst[[jmark.test]] %>% arrange(cluster, jrep)
cells.keep <- jsub.annot$cell

clstrs.uniq <- sort(unique(dat.genes$jset))
batch.uniq <- unique(jsub.annot$jrep)
colors.uniq <- cbPalette[1:length(clstrs.uniq)]
colors.uniq.batch <- cbPalette[10:(9 + length(unique(jsub.annot$jrep)))]
cols.hash <- hash::hash(clstrs.uniq, colors.uniq)
cols.hash.batch <- hash::hash(batch.uniq, colors.uniq.batch)
colsvec <- sapply(jsub.annot$cluster, function(x) AssignHash(x, cols.hash, null.fill = NA))
colsvec.batch <- sapply(jsub.annot$jrep, function(x) AssignHash(x, cols.hash.batch, null.fill = NA))
colsvec.row <- sapply(dat.genes$jset, function(x) AssignHash(x, cols.hash, null.fill = NA))


mat.filt.lda <- dat.imputed.lda[dat.genes$gene, cells.keep]
# mat.filt.glmpca <- dat.imputed.glmpca[dat.genes$gene, cells.keep]


heatmap3::heatmap3(as.matrix(mat.filt.lda), Rowv = NA, Colv = NA, ColSideColors = colsvec, scale = "row", RowSideColors = colsvec.row, revC = TRUE, main = jmark.test)

heatmap3::heatmap3(Winsorize(as.matrix(mat.filt.lda), probs = c(0, 0.96)), 
                   Rowv = NA, Colv = NA, ColSideColors = colsvec, scale = "row", RowSideColors = colsvec.row, revC = TRUE, main = jmark.test)

# heatmap3::heatmap3(as.matrix(log2(mat.filt.lda)), Rowv = NA, Colv = NA, ColSideColors = colsvec, scale = "row", RowSideColors = colsvec.row, revC = TRUE, main = jmark.test)
# heatmap3::heatmap3(as.matrix(mat.filt.glmpca), Rowv = NA, Colv = NA, ColSideColors = colsvec, scale = "row", RowSideColors = colsvec.row, revC = TRUE, main = jmark.test)




# Check new imputed -------------------------------------------------------


jsuffix <- "TES"
# inf.lda.new <- paste0(hubprefix, "jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_rep2rep3reseq_TES/lda_outputs.PZ-BM-rep3-H3K27me3-rep2rep3reseq.TES.varfilt.K-30.binarize.FALSE/ldaOut.PZ-BM-rep3-H3K27me3-rep2rep3reseq.TES.varfilt.K-30.Robj")
# inf.lda.new <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_rep2rep3reseq_TSS/lda_outputs.PZ-BM-rep3-H3K27me3-rep2rep3reseq.TSS.varfilt.K-30.binarize.FALSE/ldaOut.PZ-BM-rep3-H3K27me3-rep2rep3reseq.TSS.varfilt.K-30.Robj"
# inf.lda.new <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_rep2rep3reseq_", jsuffix, "/lda_outputs.PZ-BM-rep3-", jsuffix, "-rep2rep3reseq.", jsuffix, ".varfilt.K-30.binarize.FALSE/ldaOut.PZ-BM-rep3-H3K27me3-rep2rep3reseq.", jsuffix, ".varfilt.K-30.Robj")
inf.lda.new <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_rep2rep3reseq_TES/lda_outputs.PZ-BM-rep3-H3K27me3-rep2rep3reseq.TES.varfilt.K-30.binarize.FALSE/ldaOut.PZ-BM-rep3-H3K27me3-rep2rep3reseq.TES.varfilt.K-30.Robj"
assertthat::assert_that(file.exists(inf.lda.new))


load(inf.lda.new, v=T)

tm.result.new <- posterior(out.lda)
dat.umap.new <- DoUmapAndLouvain(topics.mat = tm.result.new$topics, jsettings = jsettings)

ggplot(dat.umap.new, aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dat.imputed.k27me3 <- t(tm.result.new$topics %*% tm.result.new$terms)

cells.keep.k27me3 <- (dat.glmpca.annot.lst[[jmark.test]] %>% arrange(cluster, jrep))$cell


inf.annot.new <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.BMrep2rep3reseq/varfilt_2/BM_rep2rep3reseq_H3K27me3.cleaned.varfilt_2.metadata.txt"
dat.annot.new <- fread(inf.annot.new)

dat.umap.new.annot <- left_join(dat.umap.new, subset(dat.metas$H3K27me3, select = c("cell", "cluster"))) %>%
  rowwise() %>%
  mutate(experi = ClipLast(cell, "_"),
         jrep = GetRepBM(experi))

ggplot(dat.umap.new.annot, aes(x = umap1 , y = umap2, color = cluster)) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# make annotatiosn from new louvains
dat.sum <- dat.umap.new.annot %>%
  group_by(louvain, cluster) %>%
  summarise(ncells = length(cluster)) 

dat.sum.best <- dat.sum %>%
  group_by(louvain) %>%
  filter(ncells == max(ncells)) %>%
  dplyr::rename(cluster2 = cluster) %>%
  mutate(cluster2 = ifelse(is.na(cluster2), "Basophils", cluster2))

dat.umap.new.annot2 <- left_join(dat.umap.new.annot, dat.sum.best)

ggplot(dat.umap.new.annot2, aes(x = umap1, y = umap2, color = cluster2)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# check heatmap


jsub.annot.new <- dat.umap.new.annot2 %>% arrange(cluster2, jrep)
cells.keep.new <- jsub.annot.new$cell

clstrs.uniq.new <- sort(unique(dat.genes$jset))
batch.uniq.new <- unique(jsub.annot.new$jrep)
colors.uniq.new <- cbPalette[1:length(clstrs.uniq.new)]
colors.uniq.batch.new <- cbPalette[10:(9 + length(unique(jsub.annot.new$jrep)))]
cols.hash.new <- hash::hash(clstrs.uniq.new, colors.uniq.new)
cols.hash.batch.new <- hash::hash(batch.uniq.new, colors.uniq.batch.new)
colsvec.new <- sapply(jsub.annot.new$cluster2, function(x) AssignHash(x, cols.hash, null.fill = NA))
colsvec.batch.new <- sapply(jsub.annot.new$jrep, function(x) AssignHash(x, cols.hash.batch, null.fill = NA))
colsvec.row.new <- sapply(dat.genes$jset, function(x) AssignHash(x, cols.hash, null.fill = NA))

cells.keep.k27me3 <- jsub.annot.new$cell
genes.keep.k27me3 <- paste("chr", sapply(dat.genes$gene, function(g) strsplit(g, ";")[[1]][[1]]), sep = "")

jfilt <- dat.imputed.k27me3[genes.keep.k27me3, cells.keep.k27me3]

jfilt.raw <- count.mat[genes.keep.k27me3, cells.keep.k27me3]

dat.genes.new <- dat.genes
dat.genes.new$coord <- genes.keep.k27me3

# mat.filt.lda <- dat.imputed.lda[dat.genes$gene, cells.keep]
heatmap3::heatmap3(as.matrix(BinarizeMatrix(jfilt.raw)), Rowv = NA, Colv = NA, ColSideColors = colsvec.new, scale = "row", RowSideColors = colsvec.row.new, revC = TRUE, main = "H3K27me3 new")



heatmap3::heatmap3(as.matrix(jfilt), Rowv = NA, Colv = NA, ColSideColors = colsvec.new, scale = "row", RowSideColors = colsvec.row.new, revC = TRUE, main = "H3K27me3 new")
heatmap3::heatmap3(Winsorize(as.matrix(jfilt), probs = c(0.00, 1)), 
                   Rowv = NA, Colv = NA, ColSideColors = colsvec.new, scale = "row", RowSideColors = colsvec.row.new, revC = TRUE, main = "H3K27me3 new")
heatmap3::heatmap3(log2(as.matrix(jfilt)), Rowv = NA, Colv = NA, ColSideColors = colsvec.new, scale = "row", RowSideColors = colsvec.row.new, revC = TRUE, main = "H3K27me3 new")

# in basophil genes, which cluster has lowest expression? 

jsetcheck <- "Granulocytes"
jfilt.long <- data.frame(cell = colnames(jfilt), expression = colMeans(jfilt[subset(dat.genes.new, jset == jsetcheck)$coord, ]), stringsAsFactors = FALSE) %>%
  left_join(., jsub.annot.new, by = c("cell"))

ggplot(jfilt.long, aes(x = cluster2, y = expression)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter() + 
  theme_bw() + 
  ggtitle(jsetcheck) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())




# Try pseudobulk again  ---------------------------------------------------

bins.lst <- split(dat.genes.new$coord, dat.genes.new$jset)

jfilt.raw.norm <- sweep(jfilt.raw, MARGIN = 2, STATS = colSums(jfilt.raw), FUN = "/")
pbulk <- SumAcrossClusters(t(jfilt.raw.norm), bins.lst)

pbulk.mat <- do.call(rbind, pbulk)

heatmap3::heatmap3(log2(as.matrix(pbulk.mat) + 1), Rowv = NA, Colv = NA, ColSideColors = colsvec.new, scale = "row", RowSideColors = colors.uniq.new, revC = TRUE, main = "H3K27me3 new")


# Pseudobulk by cluster ---------------------------------------------------

clstrs.lst <- split(dat.umap.new.annot2$cell, dat.umap.new.annot2$cluster2)

pbulk.byclst <- SumAcrossClusters(jfilt.raw.norm, clstrs.lst)

pbulk.mat.byclst <- do.call(cbind, pbulk.byclst)

heatmap3::heatmap3(log2(as.matrix(pbulk.mat.byclst) + 1), Rowv = NA, Colv = NA, ColSideColors = colors.uniq.new, scale = "row", RowSideColors = colsvec.row.new, revC = TRUE, main = "H3K27me3 new")



# Check new UMAP? ---------------------------------------------------------

inf.bins.new <- file.path(hubprefix, "jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_rep2rep3reseq_varfilt/lda_outputs.BM_rep2rep3reseq_H3K27me3.cleaned.varfilt_2.K-30.binarize.FALSE/ldaOut.BM_rep2rep3reseq_H3K27me3.cleaned.varfilt_2.K-30.Robj")
load(inf.bins.new, v=T)

tm.result.bins.new <- posterior(out.lda)
dat.umap.bins.new <- DoUmapAndLouvain(tm.result.bins.new$topics, jsettings)

dat.umap.bins.new.annot <- left_join(dat.umap.bins.new, subset(dat.umap.new.annot2, select = c("cell", "cluster2")))

ggplot(dat.umap.bins.new.annot, aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette) +
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap.bins.new.annot, aes(x = umap1, y = umap2, color = cluster2)) + 
  geom_point() + 
  facet_wrap(~cluster2) + 
  scale_color_manual(values = cbPalette) +
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
dat.impute.log <- log2(t(tm.result.bins.new$topics %*% tm.result.bins.new$terms))
dat.var <- CalculateVarAll(dat.impute.log, jchromos)

dat.umap.bins.new.annot2  <- left_join(dat.umap.bins.new.annot, dat.var)



ggplot(dat.umap.bins.new.annot2, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
  geom_point() + 
  scale_color_viridis_c(direction = -1) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())







