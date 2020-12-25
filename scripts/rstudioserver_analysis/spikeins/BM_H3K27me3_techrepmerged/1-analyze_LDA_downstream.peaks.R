# Jake Yeung
# Date of Creation: 2020-12-10
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K27me3_techrepmerged/1-analyze_LDA_downstream.peaks.R
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

hubprefix <- "/home/jyeung/hub_oudenaarden"



# Load --------------------------------------------------------------------

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

# create new metadata with filtered cells
# outrds <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.BMrep2rep3reseq/varfilt_2/BM_rep2rep3reseq_H3K27me3.cleaned.varfilt_2.peaks.rds"
# outtxt <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.BMrep2rep3reseq/varfilt_2/BM_rep2rep3reseq_H3K27me3.cleaned.varfilt_2.peaks.metadata.txt"
outtxt <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.H3K27me3_techrepmerged/BM_rep2_rep3reseq_H3K27me3.", Sys.Date(), ".txt")
outpdf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.H3K27me3_techrepmerged/BM_rep2_rep3reseq_H3K27me3.", Sys.Date(), ".pdf")
# outpdf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.BMrep2rep3reseq/varfilt_2/BM_rep2rep3reseq_H3K27me3.cleaned.varfilt_2.peaks.pdf")

pdf(outpdf, useDingbats = FALSE)

# inf.lda <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_rep2rep3reseq_TSS/lda_outputs.PZ-BM-rep3-H3K27me3-rep2rep3reseq.TSS.K-30.binarize.FALSE/ldaOut.PZ-BM-rep3-H3K27me3-rep2rep3reseq.TSS.K-30.Robj"
# inf.lda <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_rep2rep3reseq/lda_outputs.BM_rep2rep3reseq_H3K27me3.cleaned.K-30.binarize.FALSE/ldaOut.BM_rep2rep3reseq_H3K27me3.cleaned.K-30.Robj"
inf.lda <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_rep2rep3reseq_varfilt_peaks/lda_outputs.PZ-BM-rep2rep3reseq-H3K27me3.peaks.varfilt.K-30.binarize.FALSE/ldaOut.PZ-BM-rep2rep3reseq-H3K27me3.peaks.varfilt.K-30.Robj"
load(inf.lda, v=T)

tm.result <- posterior(out.lda)

tm.result <- AddTopicToTmResult(tm.result)
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
  # facet_wrap(~louvain) + 
  theme_bw() + 
  scale_color_manual(values = cbPalette, na.value = "grey25") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.merge, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  facet_wrap(~experi) + 
  theme_bw() + 
  scale_color_manual(values = cbPalette, na.value = "grey25") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Plot factor weights -----------------------------------------------------

topics.ordered <- OrderTopicsByEntropy(tm.result = tm.result)

outpdf <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/pdfs_all/topics_loadings_H3K27me3_techrepmerged/topics_H3K27me3_peaks_topics.pdf")
pdf(outpdf, useDingbats = FALSE)

for (jtopic in topics.ordered$topic){
  cell.factors <- data.frame(cell = rownames(tm.result$topics), loadings = tm.result$topics[, jtopic], stringsAsFactors = FALSE)
  m <- ggplot(left_join(dat.merge, cell.factors), aes(x = umap1, y = umap2, color = loadings)) + 
    geom_point() + 
    theme_bw() + 
    scale_color_viridis_c() + 
    ggtitle(jtopic) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
}


# Get var -----------------------------------------------------------------


dat.imputed.log <- t(log2(tm.result$topics %*% tm.result$terms))

jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
dat.var <- CalculateVarAll(dat.impute.log = dat.imputed.log, jchromos = jchromos)


# Load gene lists ----------------------------------------------------------

inf.genes <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/genesets_on_umap/geneset_on_umap.binskeep_0.niter_500.2020-12-09.from_LDA_topics.condensed.heatmap.famousgenes.keepn_400.refmark_H3K4me3.2020-12-09.txt"
genes.dat <- fread(inf.genes)
genes.dat$coord <- paste("chr", sapply(genes.dat$gene, function(g) strsplit(g, ";")[[1]][[1]]), sep = "")


# Check famous genes ------------------------------------------------------

inf.tss <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_H3K27me3_tech_rep_merged/raw_demultiplexed/count_tables_from_TSS/dist_10000/for_LDA/PZ-BM-rep3-H3K27me3-rep2rep3reseq.TSS.varfilt.rds"
dat.tss <- readRDS(inf.tss)

dat.tss.norm <- sweep(dat.tss, MARGIN = 2, STATS = colSums(dat.tss), FUN = "/")

jgene <- "Plec"
jcoord <- subset(genes.dat, grepl(jgene, gene))$coord[[1]]

print(unique(genes.dat$jset))

jjset <- "pDCs"
jjset <- "NKs"
jjset <- "Basophils"
jjset <- "pDCs"
jjset <- "Bcells"
jjset <- "Eryths"
jjset <- "Granulocytes"
jjset <- "Bcells"
jjset <- "HSPCs"

for (jjset in unique(genes.dat$jset)){
  
  coords <- subset(genes.dat, jset == jjset)$coord
  dat.tss.norm.filt <- data.frame(cell = colnames(dat.tss.norm), cutsnorm = colSums(dat.tss.norm[coords, ]), stringsAsFactors = FALSE)
  m <- ggplot(left_join(dat.merge, dat.tss.norm.filt), mapping = aes(x = umap1, y = umap2, color = log2(cutsnorm))) + 
    geom_point() + 
    theme_bw() + 
    scale_color_viridis_c(direction = -1) +
    ggtitle(paste("geneset:", jjset, "expressed genes")) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
}

dev.off()




# Check glmpca correction -------------------------------------------------

inf.glmpca <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/glmpca_outputs/H3K27me3_rep2rep3reseq.peaks.varfilt/glmpca.H3K27me3.bincutoff_0.binskeep_0.platename_plate.szname_none.niter_1000.reorder_rownames.dupfilt.suffix_peaks.RData"
load(inf.glmpca, v=T)

dat.glmpca <- DoUmapAndLouvain(glm.out$factors, jsettings = jsettings)


dat.glmpca.merged <- left_join(dat.glmpca, subset(dat.annot, select = c(-umap1, -umap2, -louvain))) %>%
  rowwise() %>%
  mutate(experi = ClipLast(cell, jsep = "_"),
         plate = as.numeric(SplitGetLast(experi, jsplit = "-")),
         rowcoord = AddPlateCoordinates(cell)$rowcoord,
         colcoord = AddPlateCoordinates(cell)$colcoord,
         jrep = GetRepBM(experiname = experi), 
         batch = AnnotateSortFromLayoutBMall(plate, rowcoord, colcoord, jrep, jmark = "H3K27me3"),
         batch = gsub("LinNeg", "Linneg", batch),
         batch = gsub("LSK", "StemCell", batch))


cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
ggplot(dat.glmpca.merged, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values = cbPalette, na.value = "grey85") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# See if plate effects are significant ------------------------------------


ggplot(dat.glmpca.merged, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  facet_wrap(~experi) + 
  theme_bw() + 
  scale_color_manual(values = cbPalette, na.value = "grey85") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.glmpca.merged, aes(x = umap1, y = umap2, color = batch)) + 
  geom_point(alpha = 0.25) + 
  facet_wrap(~cluster) + 
  theme_bw() + 
  scale_color_manual(values = cbPalette, na.value = "grey85") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Assign louvains  --------------------------------------------------------


dat.sum <- dat.glmpca.merged %>%
  group_by(louvain, cluster) %>% 
  summarise(ncells = length(cluster)) %>%
  group_by(louvain) %>%
  mutate(frac.cells = ncells / sum(ncells))

# impute NAs
dat.sum <- dat.sum %>%
  filter(ncells == max(ncells))

# set all NAs to basophils
na.cluster  <- "Basophils"

dat.sum.imputed <- dat.sum %>%
  rowwise() %>%
  mutate(cluster = ifelse(is.na(cluster), na.cluster, cluster))

clstr.new.hash <- hash(dat.sum.imputed$louvain, dat.sum.imputed$cluster)


# Make new annotations ----------------------------------------------------


dat.glmpca.merged2 <- dat.glmpca.merged %>%
  rowwise() %>%
  mutate(cluster2 = AssignHash(x = as.character(louvain), jhash = clstr.new.hash, null.fill = NULL))

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")

ggplot(dat.glmpca.merged2, aes(x = umap1, y = umap2, color = cluster2)) + 
  geom_point() +
  theme_bw() + 
  scale_color_manual(values = cbPalette, na.value = "grey85") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggplot(dat.glmpca.merged2, aes(x = umap1, y = umap2, color = batch)) + 
  geom_point(alpha = 0.25) +
  facet_wrap(~cluster2) + 
  theme_bw() + 
  # scale_color_manual(values = cbPalette, na.value = "grey85") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

jsum <- subset(dat.glmpca.merged2, select = c(louvain, cluster2)) %>%
  group_by(louvain, cluster2) %>%
  summarise(cluster3 = length(cluster2))

ggplot(dat.glmpca.merged2, aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() +
  theme_bw() + 
  scale_color_manual(values = cbPalette, na.value = "grey85") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# set louvain9 -> Basophils1
# set louvain7 -> Basophils2

dat.glmpca.merged2 <- dat.glmpca.merged2 %>%
  rowwise() %>%
  mutate(cluster3 = ifelse(louvain == "9", "Basophils1", cluster2),
         cluster3 = ifelse(louvain == "7", "Basophils2", cluster2))


ggplot(dat.glmpca.merged2, aes(x = umap1, y = umap2, color = cluster3)) + 
  geom_point() +
  theme_bw() + 
  scale_color_manual(values = cbPalette, na.value = "grey85") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Rename colnames ---------------------------------------------------------

dat.glmpca.merged.pretty <- dat.glmpca.merged2 %>%
  dplyr::rename(cluster.old = cluster, 
                cluster = cluster3,
                cluster.fewer = cluster2)

ggplot(dat.glmpca.merged.pretty, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() +
  theme_bw() + 
  scale_color_manual(values = cbPalette, na.value = "grey85") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

stypecols <- c("grey", "red", "blue")
ggplot(dat.glmpca.merged.pretty, aes(x = umap1, y = umap2, color = batch)) + 
  geom_point() +
  theme_bw() + 
  scale_color_manual(values = stypecols, na.value = "grey85") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Write outputs -----------------------------------------------------------

fwrite(dat.glmpca.merged.pretty, file = outtxt, sep = "\t")


dev.off()




