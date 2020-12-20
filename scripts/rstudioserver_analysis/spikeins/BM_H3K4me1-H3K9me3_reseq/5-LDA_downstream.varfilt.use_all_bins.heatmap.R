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

library(JFuncs)

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

pdfout <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/H3K4me1_H3K9me3_differential_expression_outputs/differential_expression_output_H3K4me1_H3K9me3.FilterByK4me1.", Sys.Date(), ".pdf")

pdf(pdfout, useDingbats = FALSE)

# inf <- file.path(hubprefix, "jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BM_dbl_reseq/lda_outputs.count_mat_cleaned_reseq.H3K4me1_H3K9me3.K-30.binarize.FALSE/ldaOut.count_mat_cleaned_reseq.H3K4me1_H3K9me3.K-30.Robj")
inf <- file.path(hubprefix, "jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BM_dbl_reseq.varfilt/lda_outputs.count_mat_cleaned_reseq.H3K4me1_H3K9me3.varfilt.K-30.binarize.FALSE/ldaOut.count_mat_cleaned_reseq.H3K4me1_H3K9me3.varfilt.K-30.Robj")
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

# jgene <- "Hlf"
jgene <- "Sox6"
jgene <- "Ebf1"
jgene <- "S100a8"
jgene <- "Pax5"
jgene <- "Meis1"
jgene <- "Hlf"
jgene <- "Tal1"
jgene <- "Hbb-y"
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
jjset <- "Bcells"
jjset <- "Granulocytes"
jjset <- "pDCs"
jjset <- "NKs"
jjset <- "DCs"
jjset <- "Eryths"
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
jjset <- "Granulocytes"
jjset <- "Eryths"
jjset <- "HSPCs"
jjset <- "Basophils"
jjset <- "Bcells"
jjset <- "NKs"
jjset <- "pDCs"
jjset <- "Basophils"

jjset <- "DCs"
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

pvals.lst1 <- lapply(jfits.lst1, function(x) x$pval)
pvals.lst2 <- lapply(jfits.lst2, function(x) x$pval)


# find diff genes in k9me3
k9.bins <- which(pvals.lst2 < 1e-10)
k4.bins <- which(pvals.lst1 < 1e-100)

k9.bins.names <- names(k9.bins)
k4.bins.names <- names(k4.bins)

bins.keep <- unique(names(pvals.lst2)[k9.bins])

bins.all <- unique(names(pvals.lst2))

jnames.all <- bins.all
names(jnames.all) <- jnames.all

jnames <- bins.keep
names(jnames) <- jnames

params.dat1.all <- lapply(jnames.all, function(jname){
  jparams <- params.lst1[[jname]]
  if (is.null(jparams)){
    return(data.frame(NULL))
  }
  # assertthat::assert_that(!is.null(jparams))
  data.frame(bin = jname, param = names(jparams), estimate1 = unlist(jparams), stringsAsFactors = FALSE)
}) %>%
  bind_rows()

params.dat1 <- lapply(jnames, function(jname){
  jparams <- params.lst1[[jname]]
  if (is.null(jparams)){
    return(data.frame(NULL))
  }
  # assertthat::assert_that(!is.null(jparams))
  data.frame(bin = jname, param = names(jparams), estimate1 = unlist(jparams), stringsAsFactors = FALSE)
}) %>%
  bind_rows()


params.dat2.all <- lapply(jnames.all, function(jname){
  jparams <- params.lst2[[jname]]
  data.frame(bin = jname, param = names(jparams), estimate2 = unlist(jparams), stringsAsFactors = FALSE)
}) %>%
  bind_rows() %>%
  mutate(param = gsub("Eryth", "Eryths", param),
         param = gsub("Lymphoid", "Bcells", param))

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



ctype2louv <- list("Bcells" = "1", 
                   "Bcells" = "2",
                   "Eryths" = "3", 
                   "Granulocytes" = "4", 
                   "Eryths" = "5", 
                   "NKs" = "6",
                   "DCs" = "7",
                   "Granulocytes" = "8",
                   "pDCs" = "9")


ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette) 


# ctype2louv.hash <- hash::hash(ctype2louv)
louv2ctype.hash <- hash::hash(unlist(ctype2louv), names(ctype2louv))

dat.umap$ctype <- sapply(as.character(dat.umap$louvain), function(x) AssignHash(x = x, jhash = louv2ctype.hash, null.fill = x))

ggplot(dat.umap, aes(x = umap1, y = umap2, color = ctype)) + 
  geom_point() + 
  theme_bw() + 
  ggtitle("H3K4me1-H3K9me3 dbl mark UMAP") + 
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
  summarise(diff.est = estimate2[[1]] - mean(c(estimate2[[2]], estimate2[[3]]))) %>%
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

# bset.sum <- params.dat.wide %>%
#   arrange(desc(eryth.effect)) %>%
#   filter(bin %in% bins.keep) %>%
#   filter(bin %in% bins.k9$rnames)

bset.sum <- params.dat.wide %>%
  arrange(granu.effect) %>%
  filter(bin %in% bins.keep) %>%
  filter(bin %in% bins.k9$rnames)

# bset.sum <- params.dat.wide %>%
#   arrange(desc(bcell.effect)) %>%
#   filter(bin %in% bins.keep) %>%
#   filter(bin %in% bins.k9$rnames)
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



# Check genome-wide H3K4me1 and H3K9me3  ----------------------------------

ests.keep <- c("ClusterBcells.Estimate", "ClusterEryths.Estimate", "ClusterGranulocytes.Estimate")
params.dat.merge <- full_join(params.dat1 %>% filter(param %in% ests.keep), params.dat2) %>%
  mutate(estimate.k4me1 = ifelse(is.na(estimate1), 0, estimate1),
         estimate.k9me3 = ifelse(is.na(estimate2), 0, estimate2))

ggplot(params.dat.merge %>% filter(estimate2 > -5), aes(x = estimate1, y = estimate2)) + 
  geom_point(alpha = 0.25) + 
  facet_wrap(~param) + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  theme_bw() + 
  xlab("H3K4me1 log2FC relative to HSPC") + 
  ylab("H3K9me3 log2FC relative to HSPC") + 
  geom_density_2d() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# do PCA to summarize
params.dat1.wide <- reshape2::dcast(params.dat1 %>% filter(param %in% ests.keep), formula = "bin ~ param", value.var = "estimate1")
params.dat2.wide <- reshape2::dcast(params.dat2, formula = "bin ~ param", value.var = "estimate2")

params.dat1.wide.all <- reshape2::dcast(params.dat1 %>% filter(param %in% ests.keep & bin %in% unique(c(names(k9.bins), names(k4.bins)))), formula = "bin ~ param", value.var = "estimate1")
params.dat2.wide.all <- reshape2::dcast(params.dat2 %>% filter(bin %in% unique(c(names(k9.bins), names(k4.bins)))), formula = "bin ~ param", value.var = "estimate2")

params.dat.wide.merge <- as.data.frame(left_join(params.dat1.wide, params.dat2.wide, by = c("bin")))
rownames(params.dat.wide.merge) <- params.dat.wide.merge$bin
params.dat.wide.merge$bin <- NULL

pca.out <- prcomp(t(params.dat.wide.merge), center = TRUE, scale. = TRUE)
dat.pca <- data.frame(cell = rownames(pca.out$x), as.data.frame(pca.out$x), stringsAsFactors = FALSE)

ggplot(dat.pca, aes(x = PC1, y = PC2, label = cell)) + 
  geom_point() + 
  geom_text() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(params.dat.wide.merge, aes(x = ClusterBcells.Estimate.y, y = ClusterGranulocytes.Estimate.y)) + 
  geom_point(alpha = 0.25) + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  theme_bw() + 
  xlab("H3K9me3 log2FC HSPC->Bcell") + 
  ylab("H3K9me3 log2FC HSPC->Granu") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(params.dat.wide.merge %>% filter(ClusterEryths.Estimate.y > -5), aes(x = ClusterBcells.Estimate.y, y = ClusterEryths.Estimate.y)) + 
  geom_point(alpha = 0.25) + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  theme_bw() + 
  xlab("H3K9me3 log2FC HSPC->Bcell") + 
  ylab("H3K9me3 log2FC HSPC->Eryth") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(params.dat.wide.merge %>% filter(ClusterEryths.Estimate.y > -5), aes(x = ClusterGranulocytes.Estimate.y, y = ClusterEryths.Estimate.y)) + 
  geom_point(alpha = 0.25) + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  theme_bw() + 
  ggtitle("K9me3 celltype-specific bins (50kb) pval < 10^-10") + 
  xlab("H3K9me3 log2FC HSPC->Bcell") + 
  ylab("H3K9me3 log2FC HSPC->Eryth") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

params.dat.merge.all <- full_join(params.dat1.all %>% filter(param %in% ests.keep), params.dat2.all) %>%
  mutate(estimate.k4me1 = ifelse(is.na(estimate1), 0, estimate1),
         estimate.k9me3 = ifelse(is.na(estimate2), 0, estimate2))

ggplot(params.dat.merge.all %>% 
         filter(estimate1 > -5 & estimate2 > -5) %>%
         filter(bin %in% c(k9.bins.names)), 
       aes(x = estimate1, y = estimate2)) + 
  geom_point(alpha = 0.25) + 
  facet_wrap(~param) +
  geom_density_2d() + 
  geom_hline(yintercept = 0, linetype = 'dotted') + 
  geom_vline(xintercept = 0, linetype = 'dotted') + 
  xlab("H3K4me1 log2FC relative to HSPC") + 
  ylab("H3K9me3 log2FC relative to HSPC") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



dev.off()



# Save objects so we can add imputed values later and make heatmap  ------------------------

outrdata <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_H3K4me1_H3K9me3_analysis/de_analysis_H3K4me1_H3K9me3.RData"
save(params.dat1.all, params.dat2.all, pvals.lst1, pvals.lst2, params.dat1.wide, params.dat2.wide, file = outrdata)



# Get bins high in each celltype  -----------------------------------------


# Get imputed  ------------------------------------------------------------


# k9.bins2 <- which(pvals.lst2 < 1e-5)

ggplot(params.dat2 %>% filter(abs(estimate2) < 5), aes(x = estimate2, fill = param)) + 
  geom_density(alpha = 0.25) + 
  facet_wrap(~param, ncol = 1) + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(params.dat.wide, aes(x = ClusterBcells.Estimate, y = ClusterGranulocytes.Estimate)) + 
  geom_point(alpha = 0.1)  + 
  theme_bw() + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(params.dat.wide, aes(x = ClusterGranulocytes.Estimate, y = ClusterBcells.Estimate)) + 
  geom_point(alpha = 0.1)  + 
  theme_bw() + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(params.dat.wide %>% filter(abs(ClusterEryths.Estimate) < 5), aes(x = ClusterEryths.Estimate, y = ClusterGranulocytes.Estimate)) + 
  geom_point(alpha = 0.1)  + 
  theme_bw() + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(params.dat.wide %>% filter(abs(ClusterEryths.Estimate) < 5), aes(x = ClusterEryths.Estimate, y = ClusterBcells.Estimate)) + 
  geom_point(alpha = 0.1)  + 
  theme_bw() + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


params.dat2.filt <- subset(params.dat2.all, bin %in% names(k9.bins2))

head(print(params.dat2.filt))

eryth.bins <- params.dat2 %>% 
  group_by(bin) %>%
  filter(estimate2[[1]] > 0 & median(estimate2) < 0)

granu.bins <- params.dat2 %>% 
  group_by(bin) %>%
  filter(estimate2[[2]] > 0 & median(estimate2) < 0)

bcells.bins <- params.dat2 %>% 
  group_by(bin) %>%
  filter(estimate2[[3]] > 0 & median(estimate2) < 0)

jsort <- params.dat.wide %>%
  group_by(bin) %>%
  # arrange(desc(bcell.effect))
  arrange(bcell.effect)

jsort <- params.dat.wide %>%
  group_by(bin) %>%
  # arrange(desc(granu.effect))
  arrange(granu.effect)

jsort <- params.dat.wide %>%
  group_by(bin) %>%
  # arrange(desc(granu.effect))
  arrange(eryth.effect)

ggplot(jsort[1:500, ], aes(x = ClusterGranulocytes.Estimate, y = ClusterBcells.Estimate)) + 
  geom_point(alpha = 0.25) + 
  theme_bw() +
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

jbins <- granu.bins$bin
jbins <- bcells.bins$bin

jbins <- jsort$bin[1:500]
# jbins <- jbins[jbins %in% names(k9.bins2)]

print(length(jbins))

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


# PCA of jsut K9me3  ------------------------------------------------------

mat.k9me3.effects <- as.matrix(params.dat.wide[, grepl("^Cluster", colnames(params.dat.wide))])
rownames(mat.k9me3.effects) <- params.dat.wide$bin

pca.out.k9me3 <- prcomp(t(mat.k9me3.effects), center = TRUE, scale. = TRUE)
dat.pca.k9me3 <- data.frame(cell = rownames(pca.out.k9me3$x), as.data.frame(pca.out.k9me3$x), stringsAsFactors = FALSE)

ggplot(dat.pca.k9me3, aes(x = PC1, y = PC2, label = cell)) +
  geom_point() + 
  geom_text() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

