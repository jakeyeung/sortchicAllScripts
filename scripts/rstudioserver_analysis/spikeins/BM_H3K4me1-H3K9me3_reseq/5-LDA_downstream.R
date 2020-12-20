# Jake Yeung
# Date of Creation: 2020-12-12
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K4me1-H3K9me3_reseq/5-LDA_downstream.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(topicmodels)

library(hash)
library(igraph)
library(umap)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)

library(DescTools)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

hubprefix <- "/home/jyeung/hub_oudenaarden"

inf <- file.path(hubprefix, "jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BM_dbl_reseq/lda_outputs.count_mat_cleaned_reseq.H3K4me1_H3K9me3.K-30.binarize.FALSE/ldaOut.count_mat_cleaned_reseq.H3K4me1_H3K9me3.K-30.Robj")
# inf <- file.path(hubprefix, "jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BM_dbl_reseq.varfilt/lda_outputs.count_mat_cleaned_reseq.H3K4me1_H3K9me3.varfilt.K-30.binarize.FALSE/ldaOut.count_mat_cleaned_reseq.H3K4me1_H3K9me3.varfilt.K-30.Robj")
load(inf, v=T)

tm.result <- posterior(out.lda)
dat.umap <- DoUmapAndLouvain(tm.result$topics, jsettings)

ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Check marker genes?  ----------------------------------------------------

jinf.tss <- file.path(hubprefix, "jyeung/data/databases/gene_tss/gene_tss_winsize.50000.bed")

jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")

coords.vec <- rownames(count.mat)

annot.bins <- AnnotateCoordsFromList.GeneWise(coords.vec = coords.vec, inf.tss = jinf.tss, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db", chromos.keep = jchromos)


# Plot raw data  ----------------------------------------------------------

# indir.raw <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/BM_H3K4me1-H3K9me3_tech_rep_merged/raw_demultiplexed/count_tables_from_H3K4me1peaks")
indir.raw <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/BM_H3K4me1-H3K9me3_tech_rep_merged/raw_demultiplexed/count_tables_from_TSS/dist_10000")
infs.raw <- list.files(indir.raw, pattern = "*.txt", full.names = TRUE)

count.mat.raw <- lapply(infs.raw, function(inf){
  ReadMatTSSFormat(inf)
}) %>%
  cbind.fill.lst(., unique(unlist(lapply(., rownames))))

count.mat.norm <- sweep(count.mat.raw, MARGIN = 2, STATS = colSums(count.mat.raw), FUN = "/")

# count.mat.norm <- sweep(count.mat, MARGIN = 2, STATS = colSums(count.mat), FUN = "/")

jgene <- "S100a8"
jgene <- "Ebf1"
jgene <- "Meis1"
# jgene <- "Hlf"
jgene <- "Pax5"
jgene <- "Sox6"
# jrow <- subset(annot.bins$out2.df, gene == jgene)$region_coord[[1]]
jrow <- grep(jgene, rownames(count.mat.norm), value = TRUE)[[1]]

raw.cuts <- data.frame(cuts = unlist(count.mat.norm[jrow, ]), cell = colnames(count.mat.norm), stringsAsFactors = FALSE)

dat.merge <- left_join(dat.umap, raw.cuts)

ggplot(dat.merge, aes(x = umap1 , y = umap2, color = ifelse(cuts > 0, 1, 0), size = log2(cuts + 1))) + 
  geom_point() + 
  theme_bw() + 
  scale_color_viridis_c() + 
  ggtitle(jgene, jrow) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Load gene sets  ---------------------------------------------------------

# celltype 
inf.genesets <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/pdfs_all/genesets_on_umap/geneset_on_umap.binskeep_0.niter_500.2020-12-09.from_LDA_topics.condensed.heatmap.famousgenes.keepn_400.refmark_H3K4me3.2020-12-09.txt")

dat.genesets <- fread(inf.genesets) %>%
  rowwise() %>%
  mutate(symbol = strsplit(gene, "\\.")[[1]][[4]],
         tssname = strsplit(gene, ";")[[1]][[2]])



# find NK cells
print(unique(dat.genesets$jset))
jjset <- "pDCs"
jjset <- "Eryths"
jjset <- "Eryths"
jjset <- "NKs"
jjset <- "DCs"
jjset <- "Bcells"
jjset <- "Granulocytes"

jjset <- "pDCs"
# jjset <- "Basophils"
jgenes <- subset(dat.genesets, jset == jjset)$tssname

jsum <- data.frame(cell = colnames(count.mat.norm), pseudogene.cuts = colMeans(count.mat.norm[jgenes, ]), stringsAsFactors = FALSE)

dat.merge2 <- left_join(dat.umap, jsum) %>%
  filter(!is.na(pseudogene.cuts))

ggplot(dat.merge2, aes(x = umap1, y = umap2, color = Winsorize(log2(pseudogene.cuts + 1), probs = c(0, 0.95)))) + 
  geom_point() + 
  theme_bw() + 
  ggtitle(jjset) + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Calculate Variance  -----------------------------------------------------

dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms))
dat.var <- CalculateVarAll(dat.impute.log = dat.impute.log, jchromos = jchromos)

dat.merge.var <- left_join(dat.umap, dat.var)

ggplot(dat.merge.var, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) +
  geom_point() + 
  theme_bw() + 
  scale_color_viridis_c(direction = -1) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Remove low variance -----------------------------------------------------

dat.merge.var <- dat.merge.var %>%
  rowwise() %>%
  mutate(low.var = cell.var.within.sum.norm < 1)

ggplot(dat.merge.var, aes(x = umap1, y = umap2, color = low.var)) +
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


cells.highvar <- subset(dat.merge.var, !low.var)$cell

count.mat.filt <- count.mat[, cells.highvar]

# outrds.filt <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.H3K4me1_H3K9me3.reseq/count_mat_cleaned_reseq.H3K4me1_H3K9me3.varfilt.rds"
# saveRDS(count.mat.filt, file = outrds.filt)



# Plot UMAP without cells -------------------------------------------------

unique(dat.genesets$jset)

jjset <- "pDCs"
jjset <- "DCs"
jjset <- "Basophils"
jjset <- "Granulocytes"
jjset <- "Eryths"
jjset <- "HSPCs"
jjset <- "Basophils"
jjset <- "Bcells"
jjset <- "NKs"
jjset <- "DCs"
jjset <- "pDCs"

jgenes <- subset(dat.genesets, jset == jjset)$tssname

jsum.filt <- data.frame(cell = colnames(count.mat.norm[, cells.highvar]), pseudogene.cuts = colMeans(count.mat.norm[jgenes, cells.highvar]), stringsAsFactors = FALSE)

dat.merge.filt <- left_join(dat.umap, jsum.filt) %>%
  filter(!is.na(pseudogene.cuts))

ggplot(dat.merge.filt, aes(x = umap1, y = umap2, color = Winsorize(log2(pseudogene.cuts + 1), probs = c(0, 0.99)))) + 
  geom_point() + 
  theme_bw() + 
  ggtitle(jjset) + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Load metadata -----------------------------------------------------------

inf.meta <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.H3K4me1_H3K9me3.reseq/metadata.H3K4me1_H3K9me3.txt"
dat.meta <- fread(inf.meta)

dat.merge.filt2 <- left_join(dat.merge.filt, dat.meta, by = c("cell" = "samp"))


ggplot(dat.merge.filt2, aes(x = umap1, y = umap2, color = stype)) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Set up bins  count ------------------------------------------------------


count.mat.bin.norm <- sweep(count.mat, MARGIN = 2, STATS = colSums(count.mat), FUN = "/")

# Get celltypes from bins  ------------------------------------------------

inf.bins1 <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.k4me1k9me3/poisson_fit_bins.H3K4me1.2020-12-12.newannot2.witherrors.MoreBins.RData")
load(inf.bins1, v=T)
jfits.lst1 <- jfits.lst

inf.bins2 <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.k4me1k9me3/poisson_fit_bins.H3K9me3.2020-12-12.newannot2.witherrors.MoreBins.RData")
load(inf.bins2, v=T)
jfits.lst2 <- jfits.lst

params.lst1 <- lapply(jfits.lst1, function(x){
  xkeep <- grepl("^Cluster.*.Estimate$", x = names(x))
  x[xkeep]
})


params.lst2 <- lapply(jfits.lst2, function(x){
  xkeep <- grepl("^Cluster.*.Estimate$", x = names(x))
  x[xkeep]
})

pvals.lst2 <- lapply(jfits.lst2, function(x) x$pval)


# find diff genes in k9me3
k9.bins <- which(pvals.lst2 < 1e-10)

bins.keep <- unique(names(pvals.lst2)[k9.bins])

jnames <- bins.keep
names(jnames) <- jnames

params.dat1 <- lapply(jnames, function(jname){
  jparams <- params.lst1[[jname]]
  if (is.null(jparams)){
    return(data.frame(NULL))
  }
  # assertthat::assert_that(!is.null(jparams))
  data.frame(bin = jname, param = names(jparams), estimate1 = unlist(jparams), stringsAsFactors = FALSE)
}) %>%
  bind_rows()

params.dat2 <- lapply(jnames, function(jname){
  jparams <- params.lst2[[jname]]
  data.frame(bin = jname, param = names(jparams), estimate2 = unlist(jparams), stringsAsFactors = FALSE)
}) %>%
  bind_rows()


params.dat2 <- params.dat2 %>%
  mutate(param = gsub("Eryth", "Eryths", param),
         param = gsub("Lymphoid", "Bcells", param))


# find eryth specific bins
bset.dat <- params.dat2 %>%
  group_by(bin) %>%
  filter(estimate2[[1]] > 0 & median(estimate2) < 0)

bset.dat <- params.dat2 %>%
  group_by(bin) %>%
  filter(estimate2[[1]] < 0 & median(estimate2) > 0)





# Assign by louvains  -----------------------------------------------------

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette) 

ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain == "12")) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette) 

ctype2louv <- list("pDCs" = "1", 
                   "Bcells1" = "2",
                   "Granulocytes" = "3", 
                   "DCs" = "4", 
                   "BadCells" = "5", 
                   "Granulocytes" = "6",
                   "Bcells2" = "7",
                   "NKs" = "8",
                   "Eryths" = "9",
                   "Bcells2" = "10",
                   "Basophils" = "11",
                   "Eryths" = "12")

# ctype2louv.hash <- hash::hash(ctype2louv)
louv2ctype.hash <- hash::hash(unlist(ctype2louv), names(ctype2louv))

dat.umap$ctype <- sapply(as.character(dat.umap$louvain), function(x) AssignHash(x = x, jhash = louv2ctype.hash, null.fill = x))

ggplot(dat.umap, aes(x = umap1, y = umap2, color = ctype)) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette) 

dat.umap.filt <- subset(dat.umap, ctype != "BadCells")



# Inf annots --------------------------------------------------------------



inf.k4me1 <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/pdfs_all/genesets_on_umap/geneset_on_umap.binskeep_0.niter_500.2020-12-04.more_celltypes.metadta.H3K4me1.txt")
inf.k9me3 <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2/BM_celltypes.bincutoff_0.binskeep_1000.byplate.szname_none.niter_500.reorder_rownames.dupfilt.2020-11-23.H3K9me3.txt")
assertthat::assert_that(file.exists(inf.k4me1))
assertthat::assert_that(file.exists(inf.k9me3))

dat.annot.k4me1 <- fread(inf.k4me1)
dat.annot.k9me3 <- fread(inf.k9me3)

# Load heterochromatins only  ---------------------------------------------

# inf.objs <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.H3K4me1_H3K9me3.reseq/K4me1_K9me3_objs.RData"  
inf.objs <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.H3K4me1_H3K9me3.reseq/K4me1_K9me3_objs.more_bins.twofilt.RData"
load(inf.objs, v=T)

inf.k4me1.morebins <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.H3K4me1_H3K9me3.reseq/H3K4me1_count_mat_all_bins.rds"
inf.k9me3.morebins <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.H3K4me1_H3K9me3.reseq/H3K9me3_count_mat_all_bins.rds"

count.k4me1.morebins <- readRDS(inf.k4me1.morebins)
count.k4me1.morebins <- sweep(count.k4me1.morebins, MARGIN = 2, STATS = colSums(count.k4me1.morebins), FUN = "/")

count.k9me3.morebins <- readRDS(inf.k9me3.morebins)
count.k9me3.morebins <- sweep(count.k9me3.morebins, MARGIN = 2, STATS = colSums(count.k9me3.morebins), FUN = "/")

params.dat.wide <- reshape2::dcast(params.dat2, formula = bin ~ param, value.var = "estimate2") %>%
  rowwise() %>%
  mutate(bcell.effect = ClusterBcells.Estimate - mean(c(ClusterEryths.Estimate, ClusterGranulocytes.Estimate)),
         eryth.effect = ClusterEryths.Estimate - mean(c(ClusterBcells.Estimate, ClusterGranulocytes.Estimate)),
         granu.effect = ClusterGranulocytes.Estimate - mean(c(ClusterBcells.Estimate, ClusterEryths.Estimate)))


bset.sum <- params.dat2 %>%
  group_by(bin) %>%
  # summarise(diff.est = estimate2[[1]] - mean(c(estimate2[[2]], estimate2[[3]]))) %>%
  # summarise(diff.est = estimate2[[2]] - mean(c(estimate2[[1]], estimate2[[3]]))) %>%
  # summarise(diff.est = estimate2[[3]] - mean(c(estimate2[[1]], estimate2[[2]]))) %>%
  # summarise(diff.est = estimate2[[2]] - estimate2[[3]])  %>%
  # summarise(diff.est = estimate2[[3]] - estimate2[[2]])  %>%
  # summarise(diff.est = estimate2[[1]] - estimate2[[2]])  %>%
  # summarise(diff.est = estimate2[[3]] - mean(c(estimate2[[2]], estimate2[[1]]))) %>%
  # arrange(desc(diff.est)) %>%
  arrange(diff.est) %>%
  filter(abs(diff.est) > 0.8) %>%
  filter(bin %in% bins.k9$rnames) %>%
  filter(bin %in% bins.keep)


# bset.sum <- params.dat2 %>%
#   group_by(bin) %>%
#   # filter(estimate2[[1]] < 0 & median(estimate2) > 0) %>%
#   # filter(estimate2[[2]] < 0 & median(estimate2) > 0) %>%
#   # filter(estimate2[[3]] < 0 & median(estimate2) > 0) %>%
#   # filter(estimate2[[1]] > 0 & median(estimate2) < 0) %>%
#   filter(bin %in% bins.keep) %>%
#   filter(bin %in% bins.k9$rnames)
# 

bset.sum <- params.dat2 %>%
  group_by(bin) %>%
  # filter(estimate2[[1]] < 0 & median(estimate2) > 0) %>%
  filter(estimate2[[2]] < 0 & median(estimate2) > 0) 
  # filter(estimate2[[3]] < 0 & median(estimate2) > 0) %>%
  # filter(estimate2[[1]] > 0 & median(estimate2) < 0) 


print(dim(bset.sum))

# 
# bset.sum <- params.dat.wide %>%
#   rowwise() %>%
#   mutate(bcell.effect = ClusterBcells.Estimate - mean(c(ClusterEryths.Estimate, ClusterGranulocytes.Estimate)),
#          eryth.effect = ClusterEryths.Estimate - mean(c(ClusterBcells.Estimate, ClusterGranulocytes.Estimate)),
#          granu.effect = ClusterGranulocytes.Estimate - mean(c(ClusterBcells.Estimate, ClusterEryths.Estimate)))
# 


bset.sum <- params.dat.wide %>%
  arrange(eryth.effect) %>%
  filter(bin %in% bins.keep) %>%
  filter(bin %in% bins.k9$rnames)

# bset.sum <- params.dat.wide %>%
#   arrange(granu.effect) %>%
#   filter(bin %in% bins.keep) %>%
#   filter(bin %in% bins.k9$rnames)

bset.sum <- params.dat.wide %>%
  arrange(desc(bcell.effect)) %>%
  filter(bin %in% bins.keep) %>%
  filter(bin %in% bins.k9$rnames)

bset.sum <- bset.sum[1:100, ]

# Plot raw cuts on UMAP of dbl --------------------------------------------

# cuts.norm.bin.filt <- colMeans(count.mat.bin.norm[unique(bset.dat$bin), cells.highvar])
jbins <- unique(bset.sum$bin)
cuts.norm.bin.filt <- colMeans(count.mat.bin.norm[jbins, cells.highvar])
# count.mat.binfilt <- colSums(count.mat[unique(bset.dat$bin), cells.highvar]

dat.merge.filt3 <- data.frame(cuts = cuts.norm.bin.filt, cell = names(cuts.norm.bin.filt), stringsAsFactors = FALSE) %>%
  left_join(., dat.merge.filt2) %>%
  left_join(., subset(dat.umap.filt, select = c(cell, ctype)))

ggplot(dat.merge.filt3, aes(x = umap1, y = umap2, color = Winsorize(log2(cuts + 1), probs = c(0, 0.99)))) + 
  geom_point() + 
  theme_bw() +
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
ggplot(dat.merge.filt3, aes(x = ctype, y = Winsorize(log2(cuts + 1), probs = c(0, 0.99)))) + 
  geom_boxplot() + 
  geom_jitter(width = 0.1) + 
  theme_bw() +
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# ggplot(dat.merge.filt3, aes(x = ctype, y = log2(cuts + 1))) + 
#   geom_boxplot() + 
#   geom_jitter(width = 0.1) + 
#   theme_bw() +
#   scale_color_viridis_c() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Check that these cuts are low in k4me1  ---------------------------------

k4me1.check <- data.frame(cuts = colMeans(count.k4me1.morebins[jbins, ]), cell = colnames(count.k4me1.morebins), stringsAsFactors = FALSE)
k9me3.check <- data.frame(cuts = colMeans(count.k9me3.morebins[jbins, ]), cell = colnames(count.k9me3.morebins), stringsAsFactors = FALSE)
# k9me3.check <- colMeans(count.k9me3.morebins[jbins, ])

# pseudobulk.filt <- colMeans(pseudobulk[bset.sum$bin, ])

plot(density(log2(k4me1.check$cuts)))
plot(density(log2(k9me3.check$cuts)))

jmerge1 <- left_join(k4me1.check, dat.annot.k4me1)
jmerge2 <- left_join(k9me3.check, dat.annot.k9me3)

# ggplot(jmerge1, aes(x = umap1, y = umap2, color = ifelse(cuts > 0, 1, 0), size = log2(cuts + 1))) +
ggplot(jmerge1, aes(x = umap1, y = umap2, color = log2(cuts))) +
  geom_point()  + 
  theme_bw() + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# ggplot(jmerge2, aes(x = umap1, y = umap2, color = ifelse(cuts > 0, 1, 0), size = log2(cuts + 1))) +
ggplot(jmerge2, aes(x = umap1, y = umap2, color = log2(cuts))) +
  geom_point()  + 
  theme_bw() + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Find celltype-specific bins  --------------------------------------------


